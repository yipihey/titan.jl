# S6a (production): Pop III protostar accreting at constant Mdot on a self-similar prescribed mass grid with SPS/OP03-type surface boundary conditions
# (see s6b_setup.jl). Phases: (1) relaxation of the isentropic initial star with accretion and burning off (dt grows 1.5x per step up to 3 t_dyn, until max|u|/c_s < 2e-3),
# (2) accretion with composition advection, burning (StarGas) and convective mixing; step control on max|d ln T, rho, r| of the zones with q > 1e-5. Output lib/popiii/s6b_<OUT>.log / .csv.
#   TITAN_NOLIB=1 julia -t 4 --project=env s6b_run.jl        options (environment, defaults):
#   star/BC:   MDOT [1.1e-3] Msun/yr, RHOC [5e-4] and TC [1e5] (central density/temperature of the embryo isentrope, ends where P = P_ram), NIN/NENV [70/50] zones, MFLOOR [3e-7], QSPLIT [0.05],
#              QC [1.12], ALPHA [1.5] MLT, OPAC [stellar|z0], KLSTAR [15] (L* from this many interfaces below the surface), ESW [0.5] / ESLIM [1.15] (damping of E_s), BCHOOK [1] (E_s(R), P_ram(R) in the
#              Jacobian), BCNR/BCSPAN, TFLOOROPAC [3000], ABSTRAD [1], WOLD [1], HSAD [1], NLIM [0.3], CONV [1e-7], XMEXT [1e-5]
#   steps:     TOL [0.03] (use 0.05-0.15), GROW [1.6], FMAX [3e-3] (cap of Mdot dt/M per step), NSKIP [20], MEND [Msun], TEND [yr], MAXSTEPS, RELAX [1], RELAXSTEPS, RELAXT [20], RELAXU [2e-3]
#   run:       EVERY, OUT, SAVE/RESTART/CKEVERY (checkpoints in /tmp/s6), PROFM [list of Msun: structure dumps], CSVEVERY, ACC [1], PSI [1] (hybrid mesh: Lagrangian interior, mass grid following the accreted mass from q = QHI [0.9] smoothly to q = QLO [1e-4]; PSI=0 = all interfaces self-similar), DTWHY, TRACE, SCAN, LEDGER, MASSDBG
#   S6a-fix:   GRAVWORK [auto|1|0] (energy-conserving gravity work on the mass grid, TitanResidFast.GRAVWORK; auto = on with accretion), CONVSTALL [1e-5] (accept Newton 2-cycles of a switch below
#              this, TitanDriver.CONV_STALL), ABSSMOOTH [0.05] (C¹ blend in ABSORB_TRAD), LSTARQ (L* for the BC at this exterior mass fraction instead of KLSTAR), LTAU [0 yr] (relaxation of that L*),
#              DTR (dt after RESTART), NITER [25], ITSLOW [20], GWDBG (gravity-work correction per interface), DECTRACE (with TRACE: branch decisions that flip between Newton iterations).
#              Production grid of S6a-fix: NIN=80 QC=1.083 QSPLIT=0.3 NENV=70 (150 zones, zone masses continuous at the split; see notes)
include(joinpath(@__DIR__, "s4_common.jl"))
include(joinpath(@__DIR__, "s6_star.jl"))
include(joinpath(@__DIR__, "s6_setup.jl"))
include(joinpath(@__DIR__, "s6b_setup.jl"))
include(joinpath(@__DIR__, "s6c_precursor.jl"))
include(joinpath(@__DIR__, "e2e_handoff.jl"))
envf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
envi(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
# E2E-A: performance recording (TitanDriver.Perf): always accumulated, written to lib/popiii/perf/<PERF>_* when PERF=<name> is set (PERFSTAGE labels the stage)
const Perf = TD.Perf; const PERFNAME = get(ENV, "PERF", "")
Perf.reset!(isempty(PERFNAME) ? "s6b" : PERFNAME; stage = get(ENV, "PERFSTAGE", "protostar"), driver = "s6b_run.jl", out = get(ENV, "OUT", "run"), mdotfile = get(ENV, "MDOTFILE", ""), mdot = get(ENV, "MDOT", ""),
            prec = get(ENV, "PRECTHICK", "0"), nin = get(ENV, "NIN", "70"), nenv = get(ENV, "NENV", "50"), tol = get(ENV, "TOL", "0.03"))
# S6b Part 3: MDOTFILE = CSV (header, then M [Msun], Mdot [Msun/yr]): time-dependent accretion rate Mdot(M*) (log-log interpolation, constant outside), e.g. from the collapse (popiii_m2c.jl SINK=1)
const MDOTTAB = Ref{Any}(nothing)
# S7b (opt-in, S7B_STAR=1; s7b_star.jl): hooks of the coupled run — S7BH[](d) after every accepted step, S7BDT[](dt) cap of the next time step, S7BMDOT[]() accretion rate (or nothing)
const S7BH = Ref{Any}(nothing); const S7BDT = Ref{Any}(nothing); const S7BMDOT = Ref{Any}(nothing)
if haskey(ENV, "MDOTFILE")
    rows = [parse.(Float64, split(l, ",")[1:2]) for l in readlines(ENV["MDOTFILE"])[2:end] if !isempty(strip(l))]
    MDOTTAB[] = (log.([r[1] for r in rows]), log.([r[2] for r in rows]))
end
function mdot_of(M)
    lx, ly = MDOTTAB[]; x = log(M / Msun)
    x <= lx[1] && return exp(ly[1]) * Msun / YR; x >= lx[end] && return exp(ly[end]) * Msun / YR
    i = searchsortedlast(lx, x); w = (x - lx[i]) / (lx[i+1] - lx[i]); exp((1 - w) * ly[i] + w * ly[i+1]) * Msun / YR
end
# MDOTRAMP [Msun/yr] (default off): the rate is capped at MDOTRAMP (M/M0GUESS)^MDOTRAMPP — as OHY10, who "start from a smaller rate and increase it with protostellar mass until the rate
# of the collapse is reached" because their protostar models do not converge at the collapse rate for the first 0.1 Msun
mdot_cap(M) = haskey(ENV, "MDOTRAMP") ? envf("MDOTRAMP", 4.41e-3) * Msun / YR * (M / (envf("M0GUESS", 0.22) * Msun))^envf("MDOTRAMPP", 2.5) : Inf
# E2E-A: MDOTSCALE [1] multiplies the table rate (a few-per-cent jitter steps the protostar over a precursor stall that depends chaotically on Mdot(M*), see notes E2E-A)
mdot_now(M) = envf("MDOTSCALE", 1.0) * min(mdot_of(M), mdot_cap(M))
Mdot = haskey(ENV, "MDOTFILE") ? mdot_now(envf("M0GUESS", 0.22) * Msun) : envf("MDOT", 1.1e-3) * Msun / YR      # (MDOTFILE: the rate at M0GUESS sets the ram pressure of the embryo)
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
# S6b Part 2: PRECTHICK=1 (default 0 = S6a): optically thick precursor with advection of the radiation and the enthalpy sink of the accreted gas (s6c_precursor.jl);
# PRECADV / PRECH [1] switch the two terms separately
const HOOKOFF = Ref(0)
const PRECTHICK = get(ENV, "PRECTHICK", "0") == "1"; const PRECADV = get(ENV, "PRECADV", "1") == "1"; const PRECH = get(ENV, "PRECH", "1") == "1"; const H0 = Ref{Any}(nothing)
const LBC = Ref(0.0)
"S6a-fix: interface of L* for the boundary conditions: KLSTAR [15] interfaces below the surface, or (LSTARQ > 0) the interface at exterior mass fraction LSTARQ."
const PROT = Ref(false)      # S7c-1 (PROTGAMMA > 0): protections of the surface condition while the star is near its Eddington limit (see update_bc!)
function kLstar()
    (envf("LSTARQ", 0.0) > 0 || PROT[]) || return b + 1 - envi("KLSTAR", 15)      # S7c-1: LSTARQ = 0 (default) as before; with PROTGAMMA the interface at q = 1e-3 while PROT (KLSTAR = 15 lies inside the convective/H-ionisation skin of a luminous, extended star)
    qL = envf("LSTARQ", 1e-3); M_ = intmass(b + 1)
    for k in b+1:-1:a+1; (M_ - intmass(k)) / M_ >= qL && return k; end
    a + 1
end
"Surface boundary conditions for the next step from the current state (lagged): P_ext, inflow, incident radiation; ALE rate."
const LASTSD = Ref{Any}(nothing); const ES_OLD = Ref(0.0)
function update_bc!()
    HOOKOFF[] -= 1
    R = rn[b+1]; M = intmass(b + 1); frn = TF.var(:frn)
    MDOTTAB[] === nothing || (global Mdot = mdot_now(M))
    S7BMDOT[] === nothing || (let m_ = S7BMDOT[](); m_ === nothing || (global Mdot = m_); end)
    if envf("PROTGAMMA", 0.0) > 0       # S7c-1: Gamma = (L* + L_acc)/L_Edd (kappa_es = 0.34) of the previous step; protections on above PROTGAMMA, off below 0.8 PROTGAMMA
        Lacc_ = G * M * Mdot / R; Γ_ = (LBC[] + Lacc_) / (4π * Z0Physics.c * G * M / 0.34)
        PROT[] = Γ_ > (PROT[] ? 0.8 : 1.0) * envf("PROTGAMMA", 0.0)
        TD.SAMEDT_RETRY[] = PROT[]
    end
    kL = kLstar()         # L* from a few interfaces below the surface (the surface flux itself is the stiff response to the boundary condition)
    Lstar = 4π * rn[kL]^2 * (frn[kL] + (LASTSD[] === nothing ? 0.0 : LASTSD[].fcv[kL]))
    # S6a-fix: L* of the precursor relaxed over LTAU [yr] (0 = off): the luminosity of the surface layers flickers by 10 when the precursor sits at the ionisation cliff
    if envf("LTAU", 0.0) > 0 && acc_on && LBC[] > 0
        w_ = min(1.0, (TF.getscalar(:timen) - TF.getscalar(:timeo)) / (envf("LTAU", 0.0) * YR)); Lstar = exp((1 - w_) * log(LBC[]) + w_ * log(max(Lstar, 1e-3 * LBC[])))
    end
    LBC[] = max(Lstar, 1e-6 * Lsun)
    vff = sqrt(2G * M / R); Pram = Mdot * vff / (4π * R^2)
    Lacc = G * M * Mdot / R
    pc = PRECTHICK ? precursor_thick(R, M, Mdot, max(Lstar, 0.0); X = SM.hmass(mat, b), advect = PRECADV, hterm = PRECH, h0 = H0[]) :
                     precursor(R, M, Mdot, max(Lstar, 0.0) + Lacc; X = SM.hmass(mat, b))
    PRECTHICK && (H0[] = pc.hR)
    TX.PEXT[] = Pram
    TF.setscalar!(:phir0, acc_on ? -Mdot / (4π * R^2) : 0.0)
    # S7c-1: the damping of E_s (weight ESW, limit ESLIM per step) is a rate per UNIT TIME: the precursor is a steady flow that adjusts on its flow time τ_E = ESTAU R/v_ff, so over a step dt < τ_E the surface
    # radiation moves only the fraction fdt = 1 − exp(−dt/τ_E) of the damped amount (fdt = 1 for dt ≫ τ_E: the old per-step rule). With a per-step rule the boundary kept moving by ESLIM per step while the
    # time-step control, which sees the response of the surface zones to that move, cut dt by 2–3 per step (dt → seconds in ~20 steps, "the S6b stall"): the move per step never became smaller.
    estau_ = PROT[] ? envf("ESTAU", 1.0) : envf("ESTAU", 0.0)
    fdt = estau_ > 0 ? -expm1(-TF.getscalar(:dtime) / (estau_ * R / vff)) : 1.0
    Es = ES_OLD[] > 0 ? exp((1 - envf("ESW", 0.5) * fdt) * log(ES_OLD[]) + envf("ESW", 0.5) * fdt * log(pc.Es)) : pc.Es
    ES_OLD[] > 0 && (Es = clamp(Es, ES_OLD[] / envf("ESLIM", 1.15)^fdt, ES_OLD[] * envf("ESLIM", 1.15)^fdt))       # the precursor flips to an optically thick, ionised state (T_s ×4): follow it in steps of at most 15 % in E_s
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
        if get(ENV, "BCHOOK", "1") == "1" && acc_on && HOOKOFF[] <= 0
            # E_s(R) and P_ram(R) as functions of the radius of the iterate: precursor integrals on a grid of radii around R, interpolated in ln R (linear, differentiable);
            # the iterate then sees the dependence on the radius (steep: opacity cliff) instead of the value of the previous step
            nR = envi("BCNR", 13); span = envf("BCSPAN", 0.06)
            lrs = collect(range(log(R) - span, log(R) + span; length = nR)); lEs = zeros(nR)
            for (i, lr) in enumerate(lrs)
                Ri = exp(lr)
                pci = PRECTHICK ? precursor_thick(Ri, M, Mdot, max(Lstar, 0.0); X = X_, n = envi("BCNPTS", 120), advect = PRECADV, hterm = PRECH, h0 = pc.hR) :
                                  precursor(Ri, M, Mdot, max(Lstar, 0.0) + G * M * Mdot / Ri; X = X_, n = envi("BCNPTS", 120))
                lEs[i] = log(pci.Es)
            end
            # S6b: BCSMOOTH = σ (in ln R; default 0 = off): Gaussian smoothing of ln E_s(R) over the nodes (the thin → thick transition of the precursor is a steep, kinked function of R
            # near r_ph = R; the Newton iteration of the top zones 2-cycles on it)
            if envf("BCSMOOTH", 0.0) > 0
                σ_ = envf("BCSMOOTH", 0.0); dl = lrs .- (log(R)); ws = [exp.(-0.5 .* ((lrs .- lrs[i]) ./ σ_) .^ 2) for i in 1:nR]
                lEs = [sum(ws[i] .* lEs) / sum(ws[i]) for i in 1:nR]
            end
            # S6b: BCSLOPE = s (default Inf = off): |d ln E_s / d ln R| ≤ s in the hook (increments limited from the centre node outwards): a steep E_s(R) (the knee where the precursor turns thick)
            # in the Jacobian of the top zones made Newton fail at every dt (the same step passes with BCHOOK = 0)
            if isfinite(envf("BCSLOPE", Inf))
                sl = envf("BCSLOPE", Inf); i0 = (nR + 1) ÷ 2
                for i in i0+1:nR; lEs[i] = clamp(lEs[i], lEs[i-1] - sl * (lrs[i] - lrs[i-1]), lEs[i-1] + sl * (lrs[i] - lrs[i-1])); end
                for i in i0-1:-1:1; lEs[i] = clamp(lEs[i], lEs[i+1] - sl * (lrs[i+1] - lrs[i]), lEs[i+1] + sl * (lrs[i+1] - lrs[i])); end
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
    kL = kLstar()
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
    (; t = TF.getscalar(:timen), M, R, Lstar, Lacc, Ltot, Teff, rph, Ts = pc.Ts, tau = pc.tau, Tc = tn[a], ρc = dn[a], Tmax, mTmax, Lpp, Ld, Lcno, Lhe, fD, conv = join(conv, ","), tKH, tacc, Sion,
     Linf = hasproperty(pc, :Linf) ? pc.Linf : Ltot, hRf = hasproperty(pc, :hR) ? pc.hR / (G * M / R) : 0.0, Mdot = Mdot * YR / Msun)
end

# ----------------------------------------------------------------------------------- energy ledger
"Total energy components of the current state [erg]: internal + radiation, kinetic, gravitational (as S5)."
function energies()
    egrn, un, dvn = TF.var(:egrn), TF.var(:un), TF.var(:dvoln)
    Eint = 0.0; Ekin = 0.0; Eg = 0.0; Egd = 0.0; mcum = 0.0
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
        Mk = k < b ? 0.5 * (dmz[j] + dmz[j+1]) : 0.5 * dmz[j]
        Egd -= G * mcum * Mk / rn[k+1]          # S6a-fix: W_d = −Σ G m_k m̄_k/r_k, the potential whose rate is the work of the momentum rows' gravity (TitanResidFast.GRAVWORK)
        fine || (Eg -= G * mcum * Mk / rn[k+1])
    end
    (; Eint, Ekin, Eg, Egd, Etot = Eint + Ekin + Eg, Etotd = Eint + Ekin + Egd)
end
mutable struct Ledger; Ed_prev::Float64; errd::Float64; errabsd::Float64; Eint_prev::Float64; Ekin_prev::Float64; Eg_prev::Float64; E_prev::Float64; E0::Float64; Ein::Float64; Eout::Float64; Ebind::Float64; Einj::Float64; Ebur::Float64; err::Float64; errabs::Float64; Eabs::Float64; end
led = Ledger(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
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
    S7BDT[] === nothing || (dtn = S7BDT[](dtn))
    haskey(ENV, "DTMAX") && (dtn = min(dtn, parse(Float64, ENV["DTMAX"])))        # accreted mass per step ≤ FMAX × M
    dtn
end
TC.DT_RULE[] = get(ENV, "NORULE", "0") == "1" ? nothing : dt_rule
TD.FAIL_TFAC[] = envf("FAILTFAC", 0.5)
TX.MLT_RCAP[] = true; TX.MLT_HSDISC[] = true; TX.MLT_WOLD[] = get(ENV, "WOLD", "1") == "1"; TX.MLT_HS[] = get(ENV, "MLTHS", "1") == "1"; TX.TAU_MLT[] = envf("TAUMLT", 1.0); TX.MLT_HSAD[] = get(ENV, "HSAD", "1") == "1"
TC.EXTRAP_FRAC[] = envf("EXTRAP", 0.3)
TD.SAMEDT_RETRY[] = get(ENV, "SAMEDTRETRY", "0") == "1"      # S7c-1: after the first Newton failure of a step retry at the same dt from the old state (see TitanDriver.SAMEDT_RETRY)
TC.PRED_LAGRHO[] = false; SM.TFLOOR[] = envf("TFLOOR", 0.0); TX.ABSORB_TRAD[] = get(ENV, "ABSTRAD", "1") == "1"; SM.TFLOOR_OPAC[] = envf("TFLOOROPAC", 3000.0); TC.NEWTON_LIM[] = envf("NLIM", 0.3)
SM.REGION_CAP[] = true; SM.install_stage_hook!(mat); SM.install_dt_hook!(mat); get(ENV, "PREDICT", "1") == "1" && SM.install_state_predictor!()
TF.setscalar!(:dtime, envf("DT0", 1.0))
get(ENV, "NOVIS", "0") == "1" && TF.setscalar!(:cqvis, 0.0)
get(ENV, "NOCONV", "0") == "1" && (TX.CONVECTION[] = false)
TD.STEPFN[] = TD.step!
# S7c-1 (opt-in): DAMPED = 1: Newton with a backtracking line search on the row-equilibrated residual (s7c1_newton.jl = DampedNewton + the star driver's hooks); DTRIG [0.05] engages it when the full correction exceeds this, DALPHA [0.5] max relative change per iterate
if get(ENV, "DAMPED", "0") == "1"
    include(joinpath(@__DIR__, "s7c1_newton.jl")); DampedNewtonS.TFAC[] = envf("FAILTFAC", 0.5); DampedNewtonS.DTRIG[] = envf("DTRIG", 0.05); DampedNewtonS.ALPHA[] = envf("DALPHA", 0.5); DampedNewtonS.LMIN[] = envf("DLMIN", 1 / 32); TD.STEPFN[] = DampedNewtonS.step!
end
TX.GHOST_RHO[] = get(ENV, "GHOSTRHO", "1") == "1"      # S6a-fix: ghost densities of the free surface explicit in the residual (quadratic Newton convergence at the surface)
TX.ABSORB_SMOOTH[] = envf("ABSSMOOTH", 0.05)      # S6a-fix: C¹ blend of max(T, T_rad) in ABSORB_TRAD (the hard switch 2-cycled Newton in optically thick zones)
TD.CONV_STALL[] = envf("CONVSTALL", 1e-5)       # S6a-fix: accept Newton 2-cycles of a switch at ≤ 1e-5 (see TitanDriver.CONV_STALL)
haskey(ENV, "GRAVWORK") && (TX.GRAVWORK[] = ENV["GRAVWORK"] == "1" ? true : ENV["GRAVWORK"] == "0" ? false : :auto)      # S6a-fix: energy-conserving gravity work on the mass grid (default :auto = on with ALE_MDOT)

const DECPREV = Dict{Tuple{Symbol,Int},Vector{Bool}}()
get(ENV, "DECTRACE", "0") == "1" && (TX.FREEZE[] = 1)
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
        if get(ENV, "DECTRACE", "0") == "1"      # branch decisions (limiter, upwind, viscosity) that changed since the last iteration
            TX.FREEZE[] = 1
            ch = String[]
            for (key, v) in TX.DEC
                old_ = get(DECPREV, key, nothing)
                old_ === nothing || for k in a-1:b+2; v[k] != old_[k] && push!(ch, "$(key[1])/$(key[2])@$(k-a+1)"); end
                DECPREV[key] = copy(v)
            end
            pr("      decisions changed: %s\n", join(ch, " "))
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
        # S6b: BCFALLBACK = n (default off): from attempt n on the surface-condition hook XIMR_HOOK (E_s(R) in the Jacobian) is dropped for that step: at the thin <-> thick transitions of
        # the precursor its Jacobian fails Newton at every dt, and the step converges at once without it; update_bc! installs the hook again for the next step
        attempt >= envi("BCFALLBACK", 1000) && attempt > 0 && (TX.XIMR_HOOK[] = nothing; HOOKOFF[] = envi("BCFALLSTEPS", 100))      # ... and for the next BCFALLSTEPS [100] steps
        res = nothing
        try
            res = Perf.@section(:solve_step, SM.solve_step!(mat, TD.STEPFN[]; tol = envf("SPLITTOL", 3e-3), maxit = envi("SPLITIT", 4), ztol = envf("ZSPLITTOL", 3e-3)))
        catch e
            e isa InterruptException && rethrow()
            FAILS[] += 1; Perf.count!(:step_fail); pr("   step failed (%s): retry with dt x 0.3\n", sprint(showerror, e)[1:min(end, 100)])
            continue
        end
        if res[1] > envf("SPLITREJ", 0.02) || res[3] > envf("ZSPLITREJ", 1e-2)
            Perf.count!(:step_reject); pr("   step rejected: source/network mismatch %.3f (zone %.1e), retry with dt x 0.3\n", res[1], res[3]); continue
        end
        # consistency of the zone masses: Δx_me (prescribed grid + mass row) against ρΔV, top zones
        xm_ = TF.var(:xmen); mmis = 0.0
        for kk in b-4:b; mmis = max(mmis, abs((xm_[kk] - xm_[kk+1]) / (dn[kk] * TF.var(:dvoln)[kk]) - 1)); end
        get(ENV, "MASSDBG", "0") == "1" && pr("      mass mismatch per zone (b-6..b): %s; xmen top: %s\n", join([@sprintf("%.2e", (xm_[kk] - xm_[kk+1]) / (dn[kk] * TF.var(:dvoln)[kk]) - 1) for kk in b-6:b], " "), join([@sprintf("%.10e", xm_[kk]) for kk in b-2:b+1], " "))
        if !(mmis < envf("MASSMIS", 5e-2)) || any(!isfinite, view(tn, a:b)) || any(<=(0), view(dn, a:b))
            Perf.count!(:step_reject); pr("   step rejected: zone-mass mismatch %.2e (or non-finite state) at dt = %.3e s, retry with dt x 0.3\n", mmis, TF.getscalar(:timen) - TF.getscalar(:timeo)); continue
        end
        dmax = 0.0
        for kk in a:b-envi("NSKIP", 20); dmax = max(dmax, abs(log(tn[kk] / to[kk])), abs(log(dn[kk] / do_[kk]))); end
        if dmax > envf("REJECT", 0.3) && attempt < 3
            Perf.count!(:step_reject); pr("   step rejected: max relative change %.3f at dt = %.3e s, retry with dt x 0.3\n", dmax, TF.getscalar(:timen) - TF.getscalar(:timeo)); continue
        end
        attempt > 0 && (COOL[] = 12)
        return res
    end
    error("robust_step!: 8 attempts failed")
end

# E2E-A: EDDCAP = Γ_max [0 = off]: soft cap of the opacity at κ_cap = Γ_max 4π c G M / L* (L* of the boundary condition, updated every step), as s5_run.jl: a star at L/L_Edd,es ≈ 0.9 has Γ > 1 in
# the He II opacity bump and its 1D envelope inflates (phase IIIb: R* 50 -> 1000 R☉, dt -> seconds) when mass loss is off
const EDDCAP = envf("EDDCAP", 0.0)
update_eddcap!() = EDDCAP > 0 && (mat.kcap = EDDCAP * 4π * Z0Physics.c * G * intmass(b + 1) / max(LBC[], 1e-6 * Lsun))
XME_START = Ref{Any}(nothing)
function after_step(k, xme_start, Mz_start)
    XME_START[] = xme_start
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    sd = TX.step_diag(dt); LASTSD[] = sd
    Eh = Perf.@section(:diag, energies()); acc = mat.acc
    Einj_h = sum(acc.Einj[a:b]); Fout_h = acc.Fout
    Ebind = acc_on ? G * (intmass(b + 1) - 0.5 * Mdot * dt) * Mdot * dt / (0.5 * (rn[b+1] + ro[b+1])) : 0.0       # gravitational energy released by the added mass (radiated at the shock): G M δm/R
    err = (Eh.Etot - led.E_prev) - (Einj_h - Fout_h - Ebind)
    errd = (Eh.Etotd - led.Ed_prev) - (Einj_h - Fout_h - Ebind)         # the same with W_d
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
    if get(ENV, "GWDBG", "0") == "1"       # S6a-fix: gravitational-work correction per interface (TitanResidFast.GRAVWORK) at the converged state
        c_ = TX.context(; dt); R_ = zeros(c_.neqn, c_.nb); ws_ = TX.getws(Float64, c_.mgr); TX.residual!(R_, ws_, c_, zeros(c_.neqn * c_.nb)); W_ = ws_.a
        tot = 0.0
        for j in a+1:b
            cj = TX.gravwork(j, c_, W_, Float64) * 4π; tot += cj
            (j > b - 25 || abs(cj) > 0.02 * abs(4π * rn[b+1]^2 * TF.var(:frn)[b+1])) && pr("   GW j=%3d q=%.2e r=%.3e Rsun C=%+.3e Lsun rho_b/rho_bar=%.3f rho=%.3e T=%.3e Mdot*e=%.3e Lsun\n", j - a, (intmass(b + 1) - intmass(j)) / intmass(b + 1), rn[j] / Rsun, cj / Lsun,
                W_.db[j] / (0.5 * (W_.d[j] + W_.d[j-1])), dn[j], tn[j], Mdot * TF.var(:egn)[j] / Lsun)
        end
        pr("   GW total %.4e Lsun (L* %.4e)\n", tot / Lsun, 4π * rn[b+1-envi("KLSTAR", 15)]^2 * TF.var(:frn)[b+1-envi("KLSTAR", 15)] / Lsun)
    end
    r = Perf.@section :burn SM.burn!(mat; sd)
    Perf.@section :mix SM.mix!(mat; sd, capmix = false, bridge = envi("BRIDGE", 4), tmax = envf("MIXTMAX", 0.02))
    led.Ein += 0.0; led.Eout += Fout_h; led.Ebind += Ebind; led.Einj += r.Einj; led.Ebur += r.Ebur; led.err += err; led.errabs += abs(err); led.errd += errd; led.errabsd += abs(errd); led.Eabs += abs(Einj_h) + abs(Fout_h) + abs(Ebind)
    adv = acc_on ? Perf.@section(:advect, SM.advect!(mat, xme_start, Mz_start; xin = Y0, order = envi("ADVORDER", 2), fac = Mdot * dt / ((TF.getscalar(:xmtot) - xme_start[end]) * 4π))) : (; cons = 0.0, cfl = 0.0, dTmax = 0.0, nsub = 0)
    MAXCONS[] = max(MAXCONS[], adv.cons)
    shift_xmext!(xmext0)
    let Eaft = energies(); led.Ed_prev = Eaft.Etotd; led.E_prev = Eaft.Etot; led.Eint_prev = Eaft.Eint; led.Ekin_prev = Eaft.Ekin; led.Eg_prev = Eaft.Eg; end
    d = Perf.@section :diag diag(sd)
    get(ENV, "DTWHY", "0") == "1" && pr("   DTWHY %d dt=%.3e s -> next %.3e s: max rel change %.4f (%s zone %d), iter %d, caps: dtcap %.3e predict %.3e; adv cfl %.2g; BC: Ts_used=%.1f Lstar=%.2f Pram=%.4e ein=%.3e rho_in=%.3e Ttop=%.1f rho_top=%.3e u_top=%.3e\n", k, dt, TF.getscalar(:dtime), DTW[][1], DTW[][3], DTW[][2], Int(TF.getscalar(:iter)), mat.dtcap, SM.predict_cap(mat), adv.cfl, BC[].Tsh, BC[].Lstar / Lsun, BC[].Pram, BC[].ein, BC[].ρin, tn[b], dn[b], TF.var(:un)[b+1])
    push!(series, (d.t / YR, dt / YR, d.M / Msun, d.R / Rsun, d.Lstar / Lsun, d.Lacc / Lsun, d.Teff, d.Tc, d.ρc, d.Tmax, d.mTmax, d.Ld / Lsun, d.Lpp / Lsun, d.Lcno / Lsun, d.fD, d.tKH / YR, d.tacc / YR, d.Sion, d.rph / Rsun, d.Ts, d.tau, d.Lhe / Lsun, d.Linf / Lsun, d.hRf, d.Mdot))
    if k % envi("EVERY", 10) == 0 || k <= 5
        pr("%5d t=%.4e yr dt=%.3e s M=%.4f R=%.2f Rsun L*=%.3f Lacc=%.2f Teff=%.0f Ts=%.0f tau=%.2f Tc=%.4e rho_c=%.3e Tmax=%.3e@%.3f LD=%.2e fD=%.4f it=%d adv(cons %.1e cfl %.1e dT %.1e) ledger(err/|terms| step %.1e cum %.1e; W_d: step %.1e cum %.1e signed %.1e) lim=%.3f %s@%d conv=%s\n", k, d.t / YR, dt, d.M / Msun, d.R / Rsun, d.Lstar / Lsun, d.Lacc / Lsun,
           d.Teff, d.Ts, d.tau, d.Tc, d.ρc, d.Tmax, d.mTmax, d.Ld / Lsun, d.fD, Int(TF.getscalar(:iter)), adv.cons, adv.cfl, adv.dTmax, err / max(abs(Einj_h) + abs(Fout_h) + abs(Ebind), 1e-300), led.errabs / max(led.Eabs, 1e-300), errd / max(abs(Einj_h) + abs(Fout_h) + abs(Ebind), 1e-300), led.errabsd / max(led.Eabs, 1e-300), led.errd / max(led.Eabs, 1e-300), DTW[][1], DTW[][3], DTW[][2], d.conv)
    end
    d
end

function write_csv()
    io = open(joinpath(@__DIR__, "..", "lib", "popiii", "s6b_$(tag).csv"), "w")
    println(io, "t_yr,dt_yr,M,R_Rsun,Lstar_Lsun,Lacc_Lsun,Teff,Tc,rho_c,Tmax,mTmax,LD,Lpp,Lcno,fD,tKH_yr,tacc_yr,Sion,rph_Rsun,Ts,tau_prec,Lhe,Linf_Lsun,hR_over_GMR,Mdot")
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
    serialize(joinpath(ckdir, "s6b_$(name).aux"), (; Y = mat.Y, Yprev = mat.Yprev, fscale = mat.fscale, fslope = mat.fslope, Dmix = mat.Dmix, ledger = mat.ledger, dtcap = mat.dtcap, capE = mat.capE, series, nstep = nstep[], pext = TX.PEXT[], xmext0, es_old = ES_OLD[], h0 = H0[], lbc = LBC[]))
end
function restart!(name)
    TF.load_state!(joinpath(ckdir, "s6b_$(name).tstate.gz"))
    aux = deserialize(joinpath(ckdir, "s6b_$(name).aux"))
    mat.Y .= aux.Y; mat.Yprev .= aux.Yprev; mat.fscale .= aux.fscale; mat.fslope .= aux.fslope; mat.Dmix = aux.Dmix; mat.ledger = aux.ledger; mat.dtcap = aux.dtcap; mat.capE = aux.capE
    append!(series, aux.series); nstep[] = aux.nstep; TX.PEXT[] = aux.pext
    hasproperty(aux, :es_old) && (ES_OLD[] = aux.es_old; H0[] = aux.h0; LBC[] = aux.lbc)      # S7c-1: boundary-condition memory (E_s damping, enthalpy guess, L*) restored with the state
    global dn, tn, rn, to, do_, ro = TF.var(:dn), TF.var(:tn), TF.var(:rn), TF.var(:to), TF.var(:do), TF.var(:ro)
    TF.setscalar!(:thet, 1.0); TF.setscalar!(:stol, 0.3); TF.setscalar!(:tmax, 1e30); TF.setscalar!(:niter, envi("NITER", 25))
    TX.REGULAR_CENTRE[] = true; TX.XMINT[] = Float64[]; TX.ALPHA_MLT[] = envf("ALPHA", 1.5); TC.FREE_OUTER[] = true
    TD.EOS[](); TD.OPAC[](); SM.refresh!()
    LASTSD[] = TX.step_diag(1.0)      # S7c-1: the convective flux of the restored state (without it L* of the boundary condition lacks the convective part for one step)
end
haskey(ENV, "RESTART") && (restart!(ENV["RESTART"]); pr("restarted from %s at step %d\n", ENV["RESTART"], nstep[]); haskey(ENV, "DTR") && TF.setscalar!(:dtime, envf("DTR", 1.0)))

# ----------------------------------------------------------------------------------- run
# hybrid mesh: interfaces deep in the star are Lagrangian (psi = 0), those in the envelope follow the accreted mass (psi = 1), smooth ramp (smoothstep in ln q) between QHI and QLO.
# The spurious source proportional to Mdot G M/R_core of the first version (both meshes) was the gravity work of the momentum rows on the mass that crosses the
# interfaces (upwinded density at steep gradients); removed by TitanResidFast.GRAVWORK (S6a-fix).
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
Perf.@section(:bc, update_bc!())
let E0_ = energies(); led.Ed_prev = E0_.Etotd; led.E_prev = E0_.Etot; led.Eint_prev = E0_.Eint; led.Ekin_prev = E0_.Ekin; led.Eg_prev = E0_.Eg; end; led.E0 = led.E_prev; const tstart0 = Ref(TF.getscalar(:timen)); const Mstart0 = Ref(intmass(b + 1))
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
        pr("   %4d q=%9.2e r=%8.5f rho=%9.3e T=%9.2f Pg=%9.3e Prad=%9.3e kappa=%9.3e tau=%9.3e Lrad=%9.3e Lconv=%9.3e u[km/s]=%9.3e\n", k - a + 1, q, rn[k] / Rsun, dn[k], tn[k], pgn[k], ern[k] / 3, chif[k], τ, 4π * rn[k+1]^2 * frn[k+1] / Lsun, LASTSD[] === nothing ? 0.0 : 4π * rn[k+1]^2 * LASTSD[].fcv[k+1] / Lsun, TF.var(:un)[k+1] / 1e5)
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
        Perf.@section(:bc, update_bc!()); robust_step!(); nstep[] += 1
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
    mat.acc.Einj .= 0.0; led.E_prev = energies().Etot; led.Ed_prev = energies().Etotd; led.E0 = led.E_prev
    TF.setscalar!(:dtime, envf("DT1", 100.0)); Perf.@section(:bc, update_bc!()); tstart0[] = TF.getscalar(:timen); Mstart0[] = intmass(b + 1)
end
if get(ENV, "S7B_STAR", "0") == "1"
    include(joinpath(@__DIR__, "s7b_exchange.jl")); include(joinpath(@__DIR__, "s7b_star.jl"))
end
maxsteps = envi("MAXSTEPS", 100); tend = envf("TEND", 1e12) * YR; mend = envf("MEND", 1e9) * Msun
try
    for k in 1:maxsteps
        xme_start = copy(TF.var(:xmen)[a:b+1]); Mz_start = 4π .* [dn[kk] * TF.var(:dvoln)[kk] for kk in a:b]
        Perf.@section(:bc, update_bc!()); prep_xmext!(TF.getscalar(:dtime))
        update_eddcap!()
        xme_start = copy(TF.var(:xmen)[a:b+1])
        Perf.@section :step begin
            Perf.@section :robust_step robust_step!()
            nstep[] += 1
            d = after_step(nstep[], xme_start, Mz_start)
            S7BH[] === nothing || S7BH[](d)
        end
        Perf.count!(:steps)
        haskey(ENV, "PROFSTEP") && nstep[] == envi("PROFSTEP", 0) && profile_top(envi("PROFN", 30))      # S7c-1 diagnostic: top-zone profile after this step
        nstep[] % 50 == 0 && Perf.sample!(step = nstep[], t_yr = d.t / YR, dt_s = TF.getscalar(:dtime), M_Msun = d.M / Msun, zones = b - a + 1, newton = Int(TF.getscalar(:iter)), wall_s = Perf.elapsed())
        Perf.@section :io begin
            (nstep[] % envi("CKEVERY", 500) == 0) && haskey(ENV, "SAVE") && checkpoint(ENV["SAVE"] * "_auto" * string((nstep[] ÷ envi("CKEVERY", 500)) % 2))
            nstep[] % envi("CSVEVERY", 200) == 0 && write_csv()
            while !isempty(profm) && d.M / Msun >= profm[1]; dump_profile(d); popfirst!(profm); end
        end
        if d.M >= mend; pr("MEND reached\n"); break; end
        if TF.getscalar(:timen) > tend; pr("TEND reached\n"); break; end
    end
catch e
    pr("STOPPED at step %d: %s\n", nstep[], sprint(showerror, e, catch_backtrace())[1:min(end, 1500)])
end
write_csv()
S7BH[] === nothing || Main.S7bStar.finish()
haskey(ENV, "SAVE") && checkpoint(ENV["SAVE"])
# E2E-A: HANDOFF_OUT=<file>: the star at the end of this run (grid, profiles, composition) for s5_run.jl (HANDOFF=<file>), see e2e_handoff.jl
haskey(ENV, "HANDOFF_OUT") && (export_star(ENV["HANDOFF_OUT"]; mat, extra = (; Mdot = Mdot * YR / Msun, tacc = (TF.getscalar(:timen) - tstart0[]) / YR, M0 = Mstart0[] / Msun)))
let d = diag(TX.step_diag(1.0)), t0 = TF.getscalar(:timen)
    dM = d.M - Mstart0[]
    pr("summary: steps %d (failed attempts %d), t = %.5e yr, M = %.6f Msun, R = %.3f Rsun\n", nstep[], FAILS[], d.t / YR, d.M / Msun, d.R / Rsun)
    pr("mass accretion check: (M - M0)/(Mdot t) = %.8f (expected 1)\n", dM / (Mdot * (d.t - tstart0[])))
    pr("composition advection: max |change of species mass - inflow|/mass = %.2e; accreted-gas abundances kept in the top zone: X_H(top) = %.5f (inflow %.5f)\n", MAXCONS[], SM.hmass(mat, b), SM.hmass(mat, b) * 0 + Y0[SM.IH] * NN.AA[SM.IH] + Y0[SM.ID] * NN.AA[SM.ID])
    pr("energy ledger: cumulative sum|err|/sum|terms| = %.3e (accretion off: ~1e-3); with W_d: sum|err|/sum|terms| = %.3e, sum err/sum|terms| = %.3e\n", led.errabs / max(led.Eabs, 1e-300), led.errabsd / max(led.Eabs, 1e-300), led.errd / max(led.Eabs, 1e-300))
end
pr("wall %.1f s\n", time() - t0wall)
if !isempty(PERFNAME)
    Perf.write(); Perf.summary()
end
