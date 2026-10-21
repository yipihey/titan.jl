# S6a (production): Pop III protostar accreting at constant Mdot on a self-similar prescribed mass grid with SPS/OP03-type surface boundary conditions
# (see s6b_setup.jl). Phases: (1) relaxation of the isentropic initial star with accretion and burning off (dt grows 1.5x per step up to 3 t_dyn, until max|u|/c_s < 2e-3),
# (2) accretion with composition advection, burning (StarGas) and convective mixing; step control on max|d ln T, rho, r| of the zones with q > 1e-5. Output lib/popiii/s6b_<OUT>.log / .csv.
#   TITAN_NOLIB=1 julia -t 4 --project=env s6b_run.jl        options (environment, defaults):
#   star/BC:   MDOT [1.1e-3] Msun/yr, RHOC [5e-4] and TC [1e5] (central density/temperature of the embryo isentrope, ends where P = P_ram), NIN/NENV [70/50] zones, MFLOOR [3e-7], QSPLIT [0.05],
#              QC [1.12], ALPHA [1.5] MLT, OPAC [stellar|z0], KLSTAR [15] (L* from this many interfaces below the surface), ESW [0.5] / ESLIM [1.15] (damping of E_s), BCHOOK [1] (E_s(R), P_ram(R) in the
#              Jacobian), BCNR/BCSPAN, TFLOOROPAC [3000], ABSTRAD [1], WOLD [1], HSAD [1], NLIM [0.3], CONV [1e-7], XMEXT [1e-5]
#   steps:     TOL [0.03] (use 0.05-0.15), GROW [1.6], FMAX [3e-3] (cap of Mdot dt/M per step), NSKIP [20], MEND [Msun], TEND [yr], MAXSTEPS, RELAX [1], RELAXSTEPS, RELAXT [20], RELAXU [2e-3]
#   run:       EVERY, OUT, SAVE/RESTART/CKEVERY (checkpoints in /tmp/s6), PROFM [list of Msun: structure dumps], CSVEVERY, ACC [1], PSI [1] (hybrid mesh: Lagrangian interior, mass grid following the accreted mass from q = QHI [0.9] smoothly to q = QLO [1e-4]; PSI=0 = all interfaces self-similar: spurious heating, see notes), DTWHY, TRACE, SCAN, LEDGER, MASSDBG
include(joinpath(@__DIR__, "s4_common.jl"))
include(joinpath(@__DIR__, "s6_star.jl"))
include(joinpath(@__DIR__, "s6_setup.jl"))
include(joinpath(@__DIR__, "s6b_setup.jl"))
envf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
envi(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
Mdot = envf("MDOT", 1.1e-3) * Msun / YR
tag = get(ENV, "OUT", "run"); logio = open(joinpath(@__DIR__, "..", "lib", "popiii", "s6b_$(tag).log"), "w")
pr(f, args...) = (s = Printf.format(Printf.Format(f), args...); print(s); print(logio, s); flush(logio); flush(stdout))
setup = build_ale_star(; rhoc = envf("RHOC", 5e-4), Mdot, Tc = envf("TC", 1e5), nin = envi("NIN", 70), nenv = envi("NENV", 50), mfloor = envf("MFLOOR", 3e-7), qsplit = envf("QSPLIT", 0.05),
                       qc = envf("QC", 1.12), Lstar = envf("LSTAR", 3.0) * Lsun, alpha = envf("ALPHA", 1.5), nuc = get(ENV, "NUC", "1") == "1", efrac = envf("EFRAC", 0.5), dXmax = envf("DXMAX", 4e-3), xmext_fac = envf("XMEXT", 1e-5), opac = Symbol(get(ENV, "OPAC", "stellar")))
mat = setup.mat; Y0 = setup.Y0
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
TF.setscalar!(:thet, 1.0); TF.setscalar!(:stol, 0.3); TF.setscalar!(:tmax, 1e30); TF.setscalar!(:niter, envi("NITER", 25))
TF.setscalar!(:conv, envf("CONV", 1e-7))
dn, tn, rn, to, do_, ro = TF.var(:dn), TF.var(:tn), TF.var(:rn), TF.var(:to), TF.var(:do), TF.var(:ro)
tol = envf("TOL", 0.03); const G = Z0Physics.G
acc_on = get(ENV, "ACC", "1") == "1"; acc_final = acc_on; phase = Ref(:evol)

# ----------------------------------------------------------------------------------- helpers
intmass(k) = (TF.getscalar(:xmtot) - TF.var(:xmen)[k]) * 4π            # interior mass of interface k [g]
zone_masses() = [4π * dn[k] * TF.var(:dvoln)[k] for k in a:b]
BC = Ref{Any}(nothing)
"Surface boundary conditions for the next step from the current state (lagged): P_ext, inflow, incident radiation; ALE rate."
const LASTSD = Ref{Any}(nothing); const ES_OLD = Ref(0.0)
function update_bc!()
    R = rn[b+1]; M = intmass(b + 1); frn = TF.var(:frn)
    kL = b + 1 - envi("KLSTAR", 15)         # L* from a few interfaces below the surface (the surface flux itself is the stiff response to the boundary condition)
    Lstar = 4π * rn[kL]^2 * (frn[kL] + (LASTSD[] === nothing ? 0.0 : LASTSD[].fcv[kL]))
    vff = sqrt(2G * M / R); Pram = Mdot * vff / (4π * R^2)
    Lacc = G * M * Mdot / R
    pc = precursor(R, M, Mdot, max(Lstar, 0.0) + Lacc; X = SM.hmass(mat, b))
    TX.PEXT[] = Pram
    TF.setscalar!(:phir0, acc_on ? -Mdot / (4π * R^2) : 0.0)
    Es = ES_OLD[] > 0 ? exp((1 - envf("ESW", 0.5)) * log(ES_OLD[]) + envf("ESW", 0.5) * log(pc.Es)) : pc.Es
    ES_OLD[] > 0 && (Es = clamp(Es, ES_OLD[] / envf("ESLIM", 1.15), ES_OLD[] * envf("ESLIM", 1.15)))       # the precursor flips to an optically thick, ionised state (T_s ×4): follow it in steps of at most 15 % in E_s
    ES_OLD[] = Es; Tsh = (Es / Z0Physics.arad)^0.25
    # the accreted gas enters in the state of the post-shock gas: pressure P_ram, temperature T_s of the precursor radiation (cold accretion: no shock heat is kept; the
    # shock energy G M Ṁ/R is radiated). (Entering with the top zone's own e_in made the top zone cool without bound: its gas is not coupled to the radiation, κ_P(T_gas) → 0.)
    X_ = SM.hmass(mat, b); TX.ZONE[] = b
    ρin = SE.rho_PT(Pram, Tsh; X = X_, Y = 1 - X_, Z = 0.0)
    ein = TX.material_gas(mat, ρin, Tsh)[2]
    TF.setscalar!(:phir2, acc_on ? -Mdot / (4π * R^2) * ein : 0.0)
    if get(ENV, "TEFFBC", "1") == "1"
        g = TF.getscalar(:geddr); cc = TF.getscalar(:cc); cpi = TF.getscalar(:cpi)
        fac_ = cc * g / ((2g + 1) * cpi)
        TF.setscalar!(:ximr, fac_ * Es)
        if get(ENV, "BCHOOK", "1") == "1" && acc_on
            # E_s(R) and P_ram(R) as functions of the radius of the iterate: precursor integrals on a grid of radii around R, interpolated in ln R (linear, differentiable);
            # the iterate then sees the dependence on the radius (steep: opacity cliff) instead of the value of the previous step
            nR = envi("BCNR", 13); span = envf("BCSPAN", 0.06)
            lrs = collect(range(log(R) - span, log(R) + span; length = nR)); lEs = zeros(nR)
            for (i, lr) in enumerate(lrs)
                Ri = exp(lr)
                pci = precursor(Ri, M, Mdot, max(Lstar, 0.0) + G * M * Mdot / Ri; X = X_, n = envi("BCNPTS", 120))
                lEs[i] = log(pci.Es)
            end
            # same damping as the value: shift the curve so that it passes through the damped E_s at R
            lEs .+= log(Es) - lEs[(nR + 1) ÷ 2]
            Mi = M
            TX.XIMR_HOOK[] = r -> begin
                x = (log(r) - lrs[1]) / (lrs[end] - lrs[1]) * (nR - 1) + 1
                i = clamp(floor(Int, ForwardDiff.value(x)), 1, nR - 1); w = x - i
                fac_ * exp((1 - w) * lEs[i] + w * lEs[i+1])
            end
            TX.PEXT_HOOK[] = r -> Mdot * sqrt(2G * Mi / r) / (4π * r^2)
        end
    end
    TD.ALE_MDOT[] = acc_on ? Mdot : 0.0
    BC[] = (; R, M, Lstar, Lacc, pc, Pram, Tsh, ρin, ein)
end
"Shift the exterior-mass offset so that x_me of the surface stays at x_mext0 (high precision of the thin surface zones)."
function shift_xmext!(xmext0)
    xmen = TF.var(:xmen); Δ = xmen[b+1] - xmext0
    for k in a-1:b+2; xmen[k] -= Δ; end
    for k in a-1:b+2; TF.var(:xmeo)[k] -= Δ; end
    TF.setscalar!(:xmtot, TF.getscalar(:xmtot) - Δ); TF.setscalar!(:xmext, TF.getscalar(:xmext) - Δ)
    nothing
end
xmext0 = TF.var(:xmen)[b+1]
"Before a step: raise the exterior-mass offset by 1.5 Mdot dt/4π so that x_me of the surface stays positive while the accreted mass leaves it (shifted back after the step)."
prep_xmext!(dt) = acc_on ? shift_xmext!(xmext0 + 1.5 * Mdot * dt / (4π)) : nothing

function diag(sd)
    frn, ern = TF.var(:frn), TF.var(:ern)
    M = intmass(b + 1); R = rn[b+1]
    kL = b + 1 - envi("KLSTAR", 15)
    Lstar = 4π * rn[kL]^2 * (frn[kL] + sd.fcv[kL]); Lacc = G * M * Mdot / R
    pc = BC[].pc
    Ltot = Lstar + Lacc
    zm = zone_masses(); mcum = cumsum(zm)
    kT = a - 1 + argmax(tn[a:b]); Tmax = tn[kT]; mTmax = (mcum[kT-a+1] - 0.5zm[kT-a+1]) / M
    Lpp = 0.0; Ld = 0.0; Lcno = 0.0; Lhe = 0.0; YD = 0.0
    for k in a:b
        Yk = view(mat.Y, :, k)
        if tn[k] > 3e5
            eg = NN.energy(Yk, NN.rates_k(tn[k], dn[k], NN.Plasma(dn[k], Yk), NN.ye_of(Yk)))[3]
            Lpp += zm[k-a+1] * eg.pp; Ld += zm[k-a+1] * eg.d; Lcno += zm[k-a+1] * eg.cno; Lhe += zm[k-a+1] * eg.he
        else
            eg = NN.energy(Yk, NN.rates_k(tn[k], dn[k], NN.Plasma(dn[k], Yk), NN.ye_of(Yk)))[3]
            Lpp += zm[k-a+1] * eg.pp; Ld += zm[k-a+1] * eg.d; Lcno += zm[k-a+1] * eg.cno; Lhe += zm[k-a+1] * eg.he
        end
        YD += zm[k-a+1] * Yk[SM.ID]
    end
    fD = YD / (sum(zm) * Y0[SM.ID])
    conv = String[]; inreg = false; m0 = 0.0
    for k in a+1:b
        c = sd.fcv[k] > 0; mk_ = mcum[k-a] / M
        if c && !inreg; inreg = true; m0 = mk_; end
        if !c && inreg; inreg = false; push!(conv, @sprintf("%.3f-%.3f", m0, mk_)); end
    end
    inreg && push!(conv, @sprintf("%.3f-1", m0))
    tKH = G * M^2 / (R * max(Lstar, 1e-30)); tacc = M / Mdot
    Teff = pc.Tph; rph = pc.rph
    x0 = 13.6 * 1.602177e-12 / (Z0Physics.kB * Teff)
    I = sum(exp(-n * x0) * (x0^2 / n + 2x0 / n^2 + 2 / n^3) for n in 1:200)
    Sion = 4π * rph^2 * 2π * (Z0Physics.kB * Teff)^3 / (6.62607e-27^3 * Z0Physics.c^2) * I
    (; t = TF.getscalar(:timen), M, R, Lstar, Lacc, Ltot, Teff, rph, Ts = pc.Ts, tau = pc.tau, Tc = tn[a], ρc = dn[a], Tmax, mTmax, Lpp, Ld, Lcno, Lhe, fD, conv = join(conv, ","), tKH, tacc, Sion)
end

# ----------------------------------------------------------------------------------- energy ledger
"Total energy components of the current state [erg]: internal + radiation, kinetic, gravitational (as S5)."
function energies()
    egrn, un, dvn = TF.var(:egrn), TF.var(:un), TF.var(:dvoln)
    Eint = 0.0; Ekin = 0.0; Eg = 0.0; mcum = 0.0
    dmz = [4π * dn[k] * dvn[k] for k in a:b]
    fine = get(ENV, "LEDGERFINE", "1") == "1"
    xg = (-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363)
    wg = (0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763)
    for k in a:b
        j = k - a + 1
        Eint += dmz[j] * egrn[k]; Ekin += dmz[j] * 0.25 * (un[k]^2 + un[k+1]^2)
        if fine      # potential energy W = −G ∫ m dm/r(m): 8-point Gauss–Legendre in every zone, r³ linear in m inside the zone
            r13, r23 = rn[k]^3, rn[k+1]^3
            for (x, w) in zip(xg, wg)
                f = 0.5 * (1 + x); dm_ = f * dmz[j]; rr = cbrt(r13 + f * (r23 - r13))
                Eg -= 0.5 * w * dmz[j] * G * (mcum + dm_) / rr
            end
        end
        mcum += dmz[j]
        if !fine
            Mk = k < b ? 0.5 * (dmz[j] + dmz[j+1]) : 0.5 * dmz[j]
            Eg -= G * mcum * Mk / rn[k+1]
        end
    end
    (; Eint, Ekin, Eg, Etot = Eint + Ekin + Eg)
end
mutable struct Ledger; Eint_prev::Float64; Ekin_prev::Float64; Eg_prev::Float64; E_prev::Float64; E0::Float64; Ein::Float64; Eout::Float64; Ebind::Float64; Einj::Float64; Ebur::Float64; err::Float64; errabs::Float64; Eabs::Float64; end
led = Ledger(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
# ----------------------------------------------------------------------------------- step control
const COOL = Ref(0); const DTW = Ref((0.0, 0, ""))
Tdyn() = sqrt(rn[b+1]^3 / (G * intmass(b + 1)))
function dt_rule(dt0)
    if phase[] === :relax
        gm = COOL[] > 0 ? 1.15 : envf("RGROW", 1.5); COOL[] = max(COOL[] - 1, 0)
        return min(gm * dt0, envf("RDTMAX", 3.0) * Tdyn())
    end
    d = 0.0; dk = 0; dv = ""
    nskip = envi("NSKIP", 20)       # the optically thin top zones (q < 1e-5) follow the boundary conditions at once: their change per step is not a measure of the accuracy of the evolution
    @inbounds for k in a:b-nskip
        x = abs(log(tn[k] / to[k])); x > d && (d = x; dk = k - a + 1; dv = "T")
        x = abs(log(dn[k] / do_[k])); x > d && (d = x; dk = k - a + 1; dv = "rho")
    end
    @inbounds for k in a+1:b+1-nskip; x = abs(log(rn[k] / ro[k])); x > d && (d = x; dk = k - a + 1; dv = "r"); end
    DTW[] = (d, dk, dv)
    gmax = COOL[] > 0 ? 1.15 : envf("GROW", 1.6); COOL[] = max(COOL[] - 1, 0)
    f = clamp(0.9 * tol / max(d, 1e-12), 0.5, gmax)
    TF.getscalar(:iter) > envi("ITSLOW", 20) && (f = min(f, 0.8))
    dtn = f * dt0
    acc_on && (dtn = min(dtn, envf("FMAX", 3e-3) * intmass(b + 1) / Mdot))
    haskey(ENV, "DTMAX") && (dtn = min(dtn, parse(Float64, ENV["DTMAX"])))        # accreted mass per step ≤ FMAX × M
    dtn
end
TC.DT_RULE[] = get(ENV, "NORULE", "0") == "1" ? nothing : dt_rule
TD.FAIL_TFAC[] = envf("FAILTFAC", 0.5)
TX.MLT_RCAP[] = true; TX.MLT_HSDISC[] = true; TX.MLT_WOLD[] = get(ENV, "WOLD", "1") == "1"; TX.MLT_HS[] = get(ENV, "MLTHS", "1") == "1"; TX.TAU_MLT[] = envf("TAUMLT", 1.0); TX.MLT_HSAD[] = get(ENV, "HSAD", "1") == "1"
TC.EXTRAP_FRAC[] = envf("EXTRAP", 0.3)
TC.PRED_LAGRHO[] = false; SM.TFLOOR[] = envf("TFLOOR", 0.0); TX.ABSORB_TRAD[] = get(ENV, "ABSTRAD", "1") == "1"; SM.TFLOOR_OPAC[] = envf("TFLOOROPAC", 3000.0); TC.NEWTON_LIM[] = envf("NLIM", 0.3)
SM.REGION_CAP[] = true; SM.install_stage_hook!(mat); SM.install_dt_hook!(mat); get(ENV, "PREDICT", "1") == "1" && SM.install_state_predictor!()
TF.setscalar!(:dtime, envf("DT0", 1.0))
get(ENV, "NOVIS", "0") == "1" && TF.setscalar!(:cqvis, 0.0)
get(ENV, "NOCONV", "0") == "1" && (TX.CONVECTION[] = false)
TD.STEPFN[] = TD.step!

if get(ENV, "TRACE", "0") == "1"
    oldmd = TD.ROUTINES[:maxdel]
    ixn = Dict(Int(TF.getscalar(s)) => string(s) for s in (:ir, :im, :id, :iu, :it, :ie, :if))
    TD.ROUTINES[:maxdel] = () -> begin
        oldmd(); kx = TF.var(:kx); dx = TF.var(:dx)
        if get(ENV, "SCAN", "0") == "1" && Int(TF.getscalar(:iter)) in (2, 3)
            c = TX.context(); neqn = c.neqn; nb = c.nb; k0 = c.ngrs - 2
            rhs = TF.var(:rhs); δ = zeros(neqn * nb)
            for kb in 1:nb, i in 1:neqn; δ[(kb-1)*neqn+i] = rhs[i, k0+kb-1]; end
            Rm = zeros(neqn, nb); ws = TX.getws(Float64, c.mgr)
            ixs = [Int(TF.getscalar(s)) for s in (:ir, :im, :id, :iu, :it, :ie, :if)]
            pr("      scan along the Newton direction (iter %d): rows of zones b-2..b  [id iu it ie if], interface rows\n", Int(TF.getscalar(:iter)))
            for λ in (-0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5)
                TX.residual!(Rm, ws, c, λ .* δ)
                row(k) = join([@sprintf("%+.3e", Rm[ixs[j], k - k0 + 1]) for j in (3, 4, 5, 6, 7)], " ")
                pr("      lam=%5.2f  z%d: %s | z%d: %s | z%d: %s\n", λ, b - 2 - a + 1, row(b - 2), b - 1 - a + 1, row(b - 1), b - a + 1, row(b))
            end
        end
        pr("   it %2d dt=%.3e dmax=%.3e  %s\n", Int(TF.getscalar(:iter)), TF.getscalar(:dtime), TF.getscalar(:dmax),
           join([@sprintf("%s:%+.2e@%d", ixn[i], dx[i], kx[i]) for i in 1:Int(TF.getscalar(:neqn))], " "))
    end
end
nstep = Ref(0); t0wall = time(); series = Tuple[]; const MAXCONS = Ref(0.0); const FAILS = Ref(0)

function robust_step!()
    snap = TD.snapshot(); dt0 = TF.getscalar(:dtime)
    for attempt in 0:7
        TD.restore!(snap); TF.setscalar!(:dtime, dt0 * 0.3^attempt)
        res = nothing
        try
            res = SM.solve_step!(mat, TD.STEPFN[]; tol = envf("SPLITTOL", 3e-3), maxit = envi("SPLITIT", 4), ztol = envf("ZSPLITTOL", 3e-3))
        catch e
            e isa InterruptException && rethrow()
            FAILS[] += 1; pr("   step failed (%s): retry with dt x 0.3\n", sprint(showerror, e)[1:min(end, 100)])
            continue
        end
        if res[1] > envf("SPLITREJ", 0.02) || res[3] > envf("ZSPLITREJ", 1e-2)
            pr("   step rejected: source/network mismatch %.3f (zone %.1e), retry with dt x 0.3\n", res[1], res[3]); continue
        end
        # consistency of the zone masses: Δx_me (prescribed grid + mass row) against ρΔV, top zones
        xm_ = TF.var(:xmen); mmis = 0.0
        for kk in b-4:b; mmis = max(mmis, abs((xm_[kk] - xm_[kk+1]) / (dn[kk] * TF.var(:dvoln)[kk]) - 1)); end
        get(ENV, "MASSDBG", "0") == "1" && pr("      mass mismatch per zone (b-6..b): %s; xmen top: %s\n", join([@sprintf("%.2e", (xm_[kk] - xm_[kk+1]) / (dn[kk] * TF.var(:dvoln)[kk]) - 1) for kk in b-6:b], " "), join([@sprintf("%.10e", xm_[kk]) for kk in b-2:b+1], " "))
        if !(mmis < envf("MASSMIS", 5e-2)) || any(!isfinite, view(tn, a:b)) || any(<=(0), view(dn, a:b))
            pr("   step rejected: zone-mass mismatch %.2e (or non-finite state) at dt = %.3e s, retry with dt x 0.3\n", mmis, TF.getscalar(:timen) - TF.getscalar(:timeo)); continue
        end
        dmax = 0.0
        for kk in a:b-envi("NSKIP", 20); dmax = max(dmax, abs(log(tn[kk] / to[kk])), abs(log(dn[kk] / do_[kk]))); end
        if dmax > envf("REJECT", 0.3) && attempt < 3
            pr("   step rejected: max relative change %.3f at dt = %.3e s, retry with dt x 0.3\n", dmax, TF.getscalar(:timen) - TF.getscalar(:timeo)); continue
        end
        attempt > 0 && (COOL[] = 12)
        return res
    end
    error("robust_step!: 8 attempts failed")
end

XME_START = Ref{Any}(nothing)
function after_step(k, xme_start, Mz_start)
    XME_START[] = xme_start
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    sd = TX.step_diag(dt); LASTSD[] = sd
    Eh = energies(); acc = mat.acc
    Einj_h = sum(acc.Einj[a:b]); Fout_h = acc.Fout
    Ebind = acc_on ? G * (intmass(b + 1) - 0.5 * Mdot * dt) * Mdot * dt / (0.5 * (rn[b+1] + ro[b+1])) : 0.0       # gravitational energy released by the added mass (radiated at the shock): G M δm/R
    err = (Eh.Etot - led.E_prev) - (Einj_h - Fout_h - Ebind)
    if get(ENV, "LEDGER", "0") == "1"       # split of the balance: (1) ΔE_int + ΔE_kin against boundary fluxes + gravity work on the fluid; (2) ΔW against the gravity work and the binding of the new mass
        un_ = TF.var(:un); gw = 0.0; mm = 0.0
        dmz_ = [4π * dn[kk] * TF.var(:dvoln)[kk] for kk in a:b]
        for kk in a+1:b+1
            mm += dmz_[kk-a]; Mk = kk <= b ? 0.5 * (dmz_[kk-a] + dmz_[kk-a+1]) : 0.5 * dmz_[kk-a]
            gw += -G * (intmass(kk)) * Mk * un_[kk] / rn[kk]^2 * dt       # work of gravity on the fluid at interface kk (new-time values), erg
        end
        # Lagrangian displacement of the matter: old r at the NEW interface masses (r³ linear in m between the old interfaces)
        xmtot_ = TF.getscalar(:xmtot); mold = [(xmtot_ - XME_START[][j]) * 4π for j in 1:b-a+2]; rold = [ro[kk] for kk in a:b+1]
        r3old = rold .^ 3; gwL = 0.0
        for kk in a+1:b+1
            mnew = intmass(kk); jj = clamp(searchsortedlast(mold, mnew), 1, length(mold) - 1)
            w_ = (mnew - mold[jj]) / (mold[jj+1] - mold[jj]); rpar = cbrt((1 - w_) * r3old[jj] + w_ * r3old[jj+1])
            Mk = kk <= b ? 0.5 * (dmz_[kk-a] + dmz_[kk-a+1]) : 0.5 * dmz_[kk-a]
            gwL += -G * mnew * Mk * (rn[kk] - rpar) / rn[kk]^2
        end
        pr("   LEDGER4: gravity work from Lagrangian displacement %.3f Lsun vs sum F u = %.3f Lsun\n", gwL / dt / Lsun, gw / dt / Lsun)
        S1 = 0.0; S2 = 0.0; Mtot = intmass(b + 1)
        for kk in a+1:b
            Mk = 0.5 * (dmz_[kk-a] + dmz_[kk-a+1]); mk = intmass(kk); Fk = Mdot * mk / Mtot
            S1 += G * mk / rn[kk] * Fk                  # potential energy carried by the mass that crosses interface kk (per second)
            S2 += G * mk * Mdot * Mk / Mtot / rn[kk]    # same with the half-zone masses growing as M(t)
        end
        pr("   LEDGER3: sum_k G m_k/r_k F_k = %.3e Lsun, sum_k G m_k M_k/(r_k) Mdot/M = %.3e Lsun\n", S1 / Lsun, S2 / Lsun)
        Ei = Eh.Eint + Eh.Ekin - (led.Eint_prev + led.Ekin_prev)
        pr("   LEDGER2 dt=%.3e s M=%.3f R=%.1f E_int=%.3e E_kin=%.3e W=%.3e erg; dE_int=%.3e dW=%.3e\n", dt, intmass(b + 1) / Msun, rn[b+1] / Rsun, Eh.Eint, Eh.Ekin, Eh.Eg, Eh.Eint - led.Eint_prev, Eh.Eg - led.Eg_prev)
        pr("   LEDGER2 %d: dE_int+dE_kin=%.3f  Einj-Fout+Wgrav=%.3f (Wgrav=%.3f) -> resid %.3f ; dW=%.3f vs -Wgrav-Ebind=%.3f -> resid %.3f [Lsun]\n", k, Ei / dt / Lsun, (Einj_h - Fout_h + gw) / dt / Lsun, gw / dt / Lsun, (Ei - (Einj_h - Fout_h + gw)) / dt / Lsun,
           (Eh.Eg - led.Eg_prev) / dt / Lsun, (-gw - Ebind) / dt / Lsun, ((Eh.Eg - led.Eg_prev) - (-gw - Ebind)) / dt / Lsun)
    end
    get(ENV, "LEDGER", "0") == "1" && pr("   LEDGER %d: dE/dt=%.3f Einj/dt=%.3f Fout/dt=%.3f Ebind/dt=%.3f err/dt=%.3f [Lsun]; L*=%.1f Mdot*ein=%.2f Mdot*P/rho=%.2f\n", k, (Eh.Etot - led.E_prev) / dt / Lsun, Einj_h / dt / Lsun, Fout_h / dt / Lsun, Ebind / dt / Lsun, err / dt / Lsun, 4π * rn[b+1]^2 * TF.var(:frn)[b+1] / Lsun, Mdot * BC[].ein / Lsun, Mdot * BC[].Pram / BC[].ρin / Lsun)
    r = SM.burn!(mat; sd)
    SM.mix!(mat; sd, capmix = false, bridge = envi("BRIDGE", 4), tmax = envf("MIXTMAX", 0.02))
    led.Ein += 0.0; led.Eout += Fout_h; led.Ebind += Ebind; led.Einj += r.Einj; led.Ebur += r.Ebur; led.err += err; led.errabs += abs(err); led.Eabs += abs(Einj_h) + abs(Fout_h) + abs(Ebind)
    adv = acc_on ? SM.advect!(mat, xme_start, Mz_start; xin = Y0, order = envi("ADVORDER", 2), fac = Mdot * dt / ((TF.getscalar(:xmtot) - xme_start[end]) * 4π)) : (; cons = 0.0, cfl = 0.0, dTmax = 0.0, nsub = 0)
    MAXCONS[] = max(MAXCONS[], adv.cons)
    shift_xmext!(xmext0)
    let Eaft = energies(); led.E_prev = Eaft.Etot; led.Eint_prev = Eaft.Eint; led.Ekin_prev = Eaft.Ekin; led.Eg_prev = Eaft.Eg; end
    d = diag(sd)
    get(ENV, "DTWHY", "0") == "1" && pr("   DTWHY %d dt=%.3e s -> next %.3e s: max rel change %.4f (%s zone %d), iter %d, caps: dtcap %.3e predict %.3e; adv cfl %.2g; BC: Ts_used=%.1f Lstar=%.2f Pram=%.4e ein=%.3e rho_in=%.3e Ttop=%.1f rho_top=%.3e u_top=%.3e\n", k, dt, TF.getscalar(:dtime), DTW[][1], DTW[][3], DTW[][2], Int(TF.getscalar(:iter)), mat.dtcap, SM.predict_cap(mat), adv.cfl, BC[].Tsh, BC[].Lstar / Lsun, BC[].Pram, BC[].ein, BC[].ρin, tn[b], dn[b], TF.var(:un)[b+1])
    push!(series, (d.t / YR, dt / YR, d.M / Msun, d.R / Rsun, d.Lstar / Lsun, d.Lacc / Lsun, d.Teff, d.Tc, d.ρc, d.Tmax, d.mTmax, d.Ld / Lsun, d.Lpp / Lsun, d.Lcno / Lsun, d.fD, d.tKH / YR, d.tacc / YR, d.Sion, d.rph / Rsun, d.Ts, d.tau, d.Lhe / Lsun))
    if k % envi("EVERY", 10) == 0 || k <= 5
        pr("%5d t=%.4e yr dt=%.3e s M=%.4f R=%.2f Rsun L*=%.3f Lacc=%.2f Teff=%.0f Ts=%.0f tau=%.2f Tc=%.4e rho_c=%.3e Tmax=%.3e@%.3f LD=%.2e fD=%.4f it=%d adv(cons %.1e cfl %.1e dT %.1e) ledger(err/|terms| step %.1e cum %.1e) lim=%.3f %s@%d conv=%s\n", k, d.t / YR, dt, d.M / Msun, d.R / Rsun, d.Lstar / Lsun, d.Lacc / Lsun,
           d.Teff, d.Ts, d.tau, d.Tc, d.ρc, d.Tmax, d.mTmax, d.Ld / Lsun, d.fD, Int(TF.getscalar(:iter)), adv.cons, adv.cfl, adv.dTmax, err / max(abs(Einj_h) + abs(Fout_h) + abs(Ebind), 1e-300), led.errabs / max(led.Eabs, 1e-300), DTW[][1], DTW[][3], DTW[][2], d.conv)
    end
    d
end

function write_csv()
    io = open(joinpath(@__DIR__, "..", "lib", "popiii", "s6b_$(tag).csv"), "w")
    println(io, "t_yr,dt_yr,M,R_Rsun,Lstar_Lsun,Lacc_Lsun,Teff,Tc,rho_c,Tmax,mTmax,LD,Lpp,Lcno,fD,tKH_yr,tacc_yr,Sion,rph_Rsun,Ts,tau_prec,Lhe")
    for row in series; println(io, join(row, ",")); end
    close(io)
end
profm = haskey(ENV, "PROFM") ? sort(parse.(Float64, split(ENV["PROFM"], ","))) : Float64[]
"Structure dump (interior mass fraction, r, rho, T, P, L, convective flag, X_H, X_D) when M passes the next value of PROFM."
function dump_profile(d)
    sd = LASTSD[]; frn = TF.var(:frn); pgn = TF.var(:pgn)
    io = open(joinpath(@__DIR__, "..", "lib", "popiii", "s6b_$(tag)_prof_$(@sprintf("%.1f", d.M / Msun)).csv"), "w")
    println(io, "# t = $(d.t / YR) yr, M = $(d.M / Msun) Msun, R = $(d.R / Rsun) Rsun, Teff_ph = $(d.Teff), L* = $(d.Lstar / Lsun) Lsun")
    println(io, "m_over_M,r_Rsun,rho,T,P_gas,L_rad_Lsun,L_conv_Lsun,conv,X_H,X_D_over_X_D0")
    mc = 0.0
    for k in a:b
        mc += 4π * dn[k] * TF.var(:dvoln)[k]
        println(io, join([mc / d.M, rn[k+1] / Rsun, dn[k], tn[k], pgn[k], 4π * rn[k+1]^2 * frn[k+1] / Lsun, 4π * rn[k+1]^2 * sd.fcv[k+1] / Lsun, sd.fcv[k+1] > 0 ? 1 : 0, SM.hmass(mat, k), mat.Y[SM.ID, k] / Y0[SM.ID]], ","))
    end
    close(io)
end
# ----------------------------------------------------------------------------------- checkpoints
using Serialization
ckdir = get(ENV, "CKDIR", "/tmp/s6")
function checkpoint(name)
    mkpath(ckdir)
    TF.save_state(joinpath(ckdir, "s6b_$(name).tstate.gz"))
    serialize(joinpath(ckdir, "s6b_$(name).aux"), (; Y = mat.Y, Yprev = mat.Yprev, fscale = mat.fscale, fslope = mat.fslope, Dmix = mat.Dmix, ledger = mat.ledger, dtcap = mat.dtcap, capE = mat.capE, series, nstep = nstep[], pext = TX.PEXT[], xmext0))
end
function restart!(name)
    TF.load_state!(joinpath(ckdir, "s6b_$(name).tstate.gz"))
    aux = deserialize(joinpath(ckdir, "s6b_$(name).aux"))
    mat.Y .= aux.Y; mat.Yprev .= aux.Yprev; mat.fscale .= aux.fscale; mat.fslope .= aux.fslope; mat.Dmix = aux.Dmix; mat.ledger = aux.ledger; mat.dtcap = aux.dtcap; mat.capE = aux.capE
    append!(series, aux.series); nstep[] = aux.nstep; TX.PEXT[] = aux.pext
    global dn, tn, rn, to, do_, ro = TF.var(:dn), TF.var(:tn), TF.var(:rn), TF.var(:to), TF.var(:do), TF.var(:ro)
    TF.setscalar!(:thet, 1.0); TF.setscalar!(:stol, 0.3); TF.setscalar!(:tmax, 1e30); TF.setscalar!(:niter, envi("NITER", 25))
    TX.REGULAR_CENTRE[] = true; TX.XMINT[] = Float64[]; TX.ALPHA_MLT[] = envf("ALPHA", 1.5); TC.FREE_OUTER[] = true
    TD.EOS[](); TD.OPAC[](); SM.refresh!()
end
haskey(ENV, "RESTART") && (restart!(ENV["RESTART"]); pr("restarted from %s at step %d\n", ENV["RESTART"], nstep[]))

# ----------------------------------------------------------------------------------- run
# hybrid mesh: interfaces deep in the star are Lagrangian (psi = 0), those in the envelope follow the accreted mass (psi = 1), smooth ramp (smoothstep in ln q) between QHI and QLO.
# With psi = 1 everywhere (the first version) the energy ledger had a spurious source proportional to Mdot G M/R_core (8 % of L_acc at 0.2 Msun, 70 % of L* at 50 Msun); this mesh closes it to 2e-3.
if get(ENV, "PSI", "1") == "1"
    qhi = envf("QHI", 0.9); qlo = envf("QLO", 1e-4); Mi0 = intmass(b + 1)
    psi = ones(length(TF.var(:xmen)))
    for k in a:b+1
        q = max((Mi0 - intmass(k)) / Mi0, 1e-30)
        x = clamp(log(qhi / q) / log(qhi / qlo), 0.0, 1.0); psi[k] = x * x * (3 - 2x)
    end
    TD.ALE_PSI[] = psi
    pr("hybrid mesh: psi = 0 for q > %.2f (%d interfaces), 1 for q < %.3f (%d interfaces)\n", qhi, count(==(0.0), psi[a:b+1]), qlo, count(==(1.0), psi[a:b+1]))
end
update_bc!()
let E0_ = energies(); led.E_prev = E0_.Etot; led.Eint_prev = E0_.Eint; led.Ekin_prev = E0_.Ekin; led.Eg_prev = E0_.Eg; end; led.E0 = led.E_prev; const tstart0 = Ref(TF.getscalar(:timen)); const Mstart0 = Ref(intmass(b + 1))
let d = diag(TX.step_diag(1.0))
    pr("S6b start: M* = %.4f Msun, R* = %.2f Rsun, Mdot = %.3e Msun/yr, %d zones, T_s(precursor) = %.0f K, tau_prec = %.2f, P_ram = %.3e, accretion %s\n", d.M / Msun, d.R / Rsun, Mdot * YR / Msun, setup.n, d.Ts, d.tau, BC[].Pram, acc_on)
end
function profile_top(nz = 30)
    chif, pgn, ern, frn = TF.var(:chifn), TF.var(:pgn), TF.var(:ern), TF.var(:frn)
    τ = 0.0; M = intmass(b + 1)
    pr("   top-zone profile: zone  q_ext  r/Rsun  rho  T  P_gas  P_rad  kappa  tau  L/Lsun  grad  grad_ad  Fconv/F\n")
    for k in b:-1:b-nz+1
        τ += chif[k] * dn[k] * (rn[k+1] - rn[k])
        q = (M - intmass(k)) / M
        pr("   %4d q=%9.2e r=%8.5f rho=%9.3e T=%9.2f Pg=%9.3e Prad=%9.3e kappa=%9.3e tau=%9.3e Lrad=%9.3e u[km/s]=%9.3e\n", k - a + 1, q, rn[k] / Rsun, dn[k], tn[k], pgn[k], ern[k] / 3, chif[k], τ, 4π * rn[k+1]^2 * frn[k+1] / Lsun, TF.var(:un)[k+1] / 1e5)
    end
end
get(ENV, "PROFILE", "0") == "1" && profile_top(envi("PROFN", 30))
if get(ENV, "RELAX", "1") == "1" && !haskey(ENV, "RESTART")
    # relaxation of the initial model (as S5): no accretion (the surface boundary conditions stay: P_ram, incident radiation of the precursor), burning off, dt grows by 1.5 per step up to
    # 3 t_dyn (backward Euler damps the dynamical modes), until t > RELAXT t_dyn and max|u|/c_s < RELAXU in the star
    acc_on = false; phase[] = :relax; mat.nuc = false; TD.ALE_MDOT[] = 0.0
    pr("== relaxation (accretion and burning off), t_dyn = %.3e s\n", Tdyn())
    trel0 = TF.getscalar(:timen)
    for k in 1:envi("RELAXSTEPS", 600)
        xme_start = copy(TF.var(:xmen)[a:b+1]); Mz_start = 4π .* [dn[kk] * TF.var(:dvoln)[kk] for kk in a:b]
        update_bc!(); robust_step!(); nstep[] += 1
        sd = TX.step_diag(TF.getscalar(:timen) - TF.getscalar(:timeo)); LASTSD[] = sd
        shift_xmext!(xmext0)
        un = TF.var(:un); asn = TF.var(:asn)
        umax = maximum(abs(un[kk]) / asn[min(kk, b)] for kk in a:b-15)
        nstep[] % 20 == 0 && pr("relax %5d t=%.4e s (%.1f t_dyn) dt=%.3e s R=%.4f Rsun Tc=%.5e max|u|/cs(core)=%.2e it=%d\n", nstep[], TF.getscalar(:timen) - trel0, (TF.getscalar(:timen) - trel0) / Tdyn(), TF.getscalar(:timen) - TF.getscalar(:timeo), rn[b+1] / Rsun, tn[a], umax, Int(TF.getscalar(:iter)))
        if TF.getscalar(:timen) - trel0 > envf("RELAXT", 20.0) * Tdyn() && umax < envf("RELAXU", 2e-3)
            pr("relaxed at step %d: t = %.3e s = %.1f t_dyn, max|u|/c_s = %.2e\n", nstep[], TF.getscalar(:timen) - trel0, (TF.getscalar(:timen) - trel0) / Tdyn(), umax); break
        end
    end
    acc_on = acc_final; phase[] = :evol; mat.nuc = get(ENV, "NUC", "1") == "1"
    mat.acc.Einj .= 0.0; led.E_prev = energies().Etot; led.E0 = led.E_prev
    TF.setscalar!(:dtime, envf("DT1", 100.0)); update_bc!(); tstart0[] = TF.getscalar(:timen); Mstart0[] = intmass(b + 1)
end
maxsteps = envi("MAXSTEPS", 100); tend = envf("TEND", 1e12) * YR; mend = envf("MEND", 1e9) * Msun
try
    for k in 1:maxsteps
        xme_start = copy(TF.var(:xmen)[a:b+1]); Mz_start = 4π .* [dn[kk] * TF.var(:dvoln)[kk] for kk in a:b]
        update_bc!(); prep_xmext!(TF.getscalar(:dtime))
        xme_start = copy(TF.var(:xmen)[a:b+1])
        robust_step!()
        nstep[] += 1
        d = after_step(nstep[], xme_start, Mz_start)
        (nstep[] % envi("CKEVERY", 500) == 0) && haskey(ENV, "SAVE") && checkpoint(ENV["SAVE"] * "_auto" * string((nstep[] ÷ envi("CKEVERY", 500)) % 2))
        nstep[] % envi("CSVEVERY", 200) == 0 && write_csv()
        while !isempty(profm) && d.M / Msun >= profm[1]; dump_profile(d); popfirst!(profm); end
        if d.M >= mend; pr("MEND reached\n"); break; end
        if TF.getscalar(:timen) > tend; pr("TEND reached\n"); break; end
    end
catch e
    pr("STOPPED at step %d: %s\n", nstep[], sprint(showerror, e, catch_backtrace())[1:min(end, 1500)])
end
write_csv()
haskey(ENV, "SAVE") && checkpoint(ENV["SAVE"])
let d = diag(TX.step_diag(1.0)), t0 = TF.getscalar(:timen)
    dM = d.M - Mstart0[]
    pr("summary: steps %d (failed attempts %d), t = %.5e yr, M = %.6f Msun, R = %.3f Rsun\n", nstep[], FAILS[], d.t / YR, d.M / Msun, d.R / Rsun)
    pr("mass accretion check: (M - M0)/(Mdot t) = %.8f (expected 1)\n", dM / (Mdot * (d.t - tstart0[])))
    pr("composition advection: max |change of species mass - inflow|/mass = %.2e; accreted-gas abundances kept in the top zone: X_H(top) = %.5f (inflow %.5f)\n", MAXCONS[], SM.hmass(mat, b), SM.hmass(mat, b) * 0 + Y0[SM.IH] * NN.AA[SM.IH] + Y0[SM.ID] * NN.AA[SM.ID])
    pr("energy ledger: cumulative sum|err|/sum|terms| = %.3e (accretion off: ~1e-3)\n", led.errabs / max(led.Eabs, 1e-300))
end
pr("wall %.1f s\n", time() - t0wall)
