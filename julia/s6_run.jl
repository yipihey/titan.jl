# S6a, approach A: accreting Pop III protostar with the accretion shock in the domain (free-fall inflow at r_out, adaptive grid lgrid = 1),
# radiation hydrodynamics, StarGas burning + mixing, composition advection through the moving grid.
#   TITAN_NOLIB=1 julia -t 4 --project=env s6_run.jl        options (environment): M [0.1] Msun (the isentrope stops where P = P_ram), MDOT [1.1e-3] Msun/yr, TC [1e5] K,
#   NIN/NENV/NINF [60/40/70], ROUT [4] R*, DR0 [1e-3], TAU [1e7] grid time constant, GW [1,1,1] weights of ln m, ln ρ, ln T, TOL [0.03], GROW [1.6], MAXSTEPS, TEND [yr],
#   MEND [Msun] stop mass, OUT [name] log/csv tag, SAVE/RESTART checkpoints, EVERY [10]
include(joinpath(@__DIR__, "s4_common.jl"))
include(joinpath(@__DIR__, "s6_star.jl"))
include(joinpath(@__DIR__, "s6_setup.jl"))
envf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
envi(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
Mdot = envf("MDOT", 1.1e-3) * Msun / YR
tag = get(ENV, "OUT", "run"); logio = open(joinpath(@__DIR__, "..", "lib", "popiii", "s6_$(tag).log"), "w")
pr(f, args...) = (s = Printf.format(Printf.Format(f), args...); print(s); print(logio, s); flush(logio); flush(stdout))

# isentrope: central density from RHOC (default: pick the isentrope that reaches M at the contact)
setup = build_accreting(; rhoc = envf("RHOC", 5e-4), Mdot, Tc = envf("TC", 1e5), nin = envi("NIN", 60), nenv = envi("NENV", 40), ninf = envi("NINF", 70), rout = envf("ROUT", 4.0), dr0 = envf("DR0", 1e-3),
                        Lstar = envf("LSTAR", 3.0) * Lsun, gridtau = envf("TAU", 1e7), gridw = Tuple(parse.(Float64, split(get(ENV, "GW", "1,1,1"), ","))), nuc = get(ENV, "NUC", "1") == "1",
                        xmext_fac = envf("XMEXT", 100.0), efrac = envf("EFRAC", 0.5), dXmax = envf("DXMAX", 4e-3), alpha = envf("ALPHA", 1.5))
mat = setup.mat
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
TF.setscalar!(:thet, 1.0); TF.setscalar!(:stol, 0.3); TF.setscalar!(:tmax, 1e30); TF.setscalar!(:niter, envi("NITER", 25))
haskey(ENV, "CONV") && TF.setscalar!(:conv, parse(Float64, ENV["CONV"]))
dn, tn, rn, to, do_, ro = TF.var(:dn), TF.var(:tn), TF.var(:rn), TF.var(:to), TF.var(:do), TF.var(:ro)
tol = envf("TOL", 0.03)
const G = Z0Physics.G

# ----------------------------------------------------------------------------------- diagnostics
"Interior mass (g) of interface k: xmtot − xmen (per sr) × 4π."
intmass(k) = (TF.getscalar(:xmtot) - TF.var(:xmen)[k]) * 4π
function zone_masses()
    [4π * dn[k] * TF.var(:dvoln)[k] for k in a:b]
end
"Shock (stellar surface) radius: outermost interface where |u| < ½ v_ff(r) going inwards from the boundary; returns (R, k_interface)."
function shock_radius()
    un = TF.var(:un)
    for k in b+1:-1:a+1
        vff = sqrt(2G * intmass(k) / rn[k])
        if abs(un[k]) < 0.5 * vff && k < b + 1
            # linear interpolation between k and k+1
            f1 = abs(un[k]) / vff - 0.5; f2 = abs(un[k+1]) / sqrt(2G * intmass(k + 1) / rn[k+1]) - 0.5
            w = f1 / (f1 - f2)
            return rn[k] + w * (rn[k+1] - rn[k]), k
        end
    end
    rn[a+1], a + 1
end
"Photosphere of the infall (Rosseland optical depth 2/3 from the outer boundary): radius and the index."
function photosphere()
    chif = TF.var(:chifn); τ = 0.0
    for k in b:-1:a
        dτ = chif[k] * (rn[k+1] - rn[k])
        if τ + dτ >= 2 / 3
            w = (2 / 3 - τ) / dτ
            return rn[k+1] - w * (rn[k+1] - rn[k]), k
        end
        τ += dτ
    end
    rn[a], a
end
function diag(sd)
    frn, ern = TF.var(:frn), TF.var(:ern)
    M = intmass(b + 1)
    Rsh, ksh = shock_radius()
    # base of the settling layer: innermost interface (searching inwards from the shock) with |u| < 2 % v_ff... use the interface 1 zone inside the shock
    kb_ = max(ksh - 2, a + 2)
    Lint(k) = 4π * rn[k]^2 * (frn[k] + sd.fcv[k])
    Lstar = Lint(kb_)
    Mst = intmass(kb_)
    Lacc = G * Mst * Mdot / Rsh
    rph, kph = photosphere()
    Lph = 4π * rn[kph+1]^2 * frn[kph+1]
    Teff = (Lph / (4π * rph^2 * Z0Physics.σSB))^0.25
    Ltot_out = 4π * rn[b+1]^2 * frn[b+1]
    # T max and its mass coordinate (star only, below the shock)
    zm = zone_masses(); mcum = cumsum(zm)
    kT = a - 1 + argmax(tn[a:kb_]); Tmax = tn[kT]; mTmax = (mcum[kT-a+1] - 0.5zm[kT-a+1]) / Mst
    # nuclear luminosities by channel
    Lpp = 0.0; Ld = 0.0; Lcno = 0.0; Lhe = 0.0; YD = 0.0; YD0 = mat.Y[SM.ID, b]
    for k in a:kb_-1
        Yk = view(mat.Y, :, k)
        eg = NN.energy(Yk, NN.rates_k(tn[k], dn[k], NN.Plasma(dn[k], Yk), NN.ye_of(Yk)))[3]
        mk = zm[k-a+1]; Lpp += mk * eg.pp; Ld += mk * eg.d; Lcno += mk * eg.cno; Lhe += mk * eg.he
        YD += mk * Yk[SM.ID]
    end
    Mstar_zones = sum(zm[1:kb_-a])
    fD = YD / (Mstar_zones * setup.Y0[SM.ID])
    # convective zones (mass fractions, interfaces with F_conv > 0)
    conv = String[]; inreg = false; m0 = 0.0
    for k in a+1:kb_
        c = sd.fcv[k] > 0
        mk_ = (mcum[k-a] ) / Mst
        if c && !inreg; inreg = true; m0 = mk_; end
        if !c && inreg; inreg = false; push!(conv, @sprintf("%.3f-%.3f", m0, mk_)); end
    end
    inreg && push!(conv, @sprintf("%.3f-%.3f", m0, 1.0))
    tKH = G * Mst^2 / (Rsh * max(Lstar, 1e-30)); tacc = Mst / Mdot
    # ionising photon rate of a black body of T_eff and radius r_ph
    x0 = 13.6 * 1.602177e-12 / (Z0Physics.kB * Teff)
    I = sum(exp(-n * x0) * (x0^2 / n + 2x0 / n^2 + 2 / n^3) for n in 1:200)
    Sion = 4π * rph^2 * 2π * (Z0Physics.kB * Teff)^3 / (6.62607e-27^3 * Z0Physics.c^2) * I
    (; t = TF.getscalar(:timen), M = Mst, Mtot = M, R = Rsh, Lstar, Lacc, Ltot = Lstar + Lacc, Lph, Teff, rph, Tc = tn[a], ρc = dn[a], Tmax, mTmax, Lpp, Ld, Lcno, Lhe, fD, conv = join(conv, ","), tKH, tacc, Sion, Lout = Ltot_out, kb = kb_, ksh)
end

# ----------------------------------------------------------------------------------- step control
const COOL = Ref(0); const DTW = Ref((0.0, 0, ""))
function dt_rule(dt0)
    d = 0.0; dk = 0; dv = ""
    @inbounds for k in a:b
        x = abs(log(tn[k] / to[k])); x > d && (d = x; dk = k - a + 1; dv = "T")
        x = abs(log(dn[k] / do_[k])); x > d && (d = x; dk = k - a + 1; dv = "rho")
    end
    DTW[] = (d, dk, dv)
    gmax = COOL[] > 0 ? 1.15 : envf("GROW", 1.6); COOL[] = max(COOL[] - 1, 0)
    f = clamp(0.9 * tol / max(d, 1e-12), 0.5, gmax)
    TF.getscalar(:iter) > envi("ITSLOW", 20) && (f = min(f, 0.8))
    f * dt0
end
TC.DT_RULE[] = get(ENV, "NORULE", "0") == "1" ? nothing : dt_rule
TD.FAIL_TFAC[] = envf("FAILTFAC", 0.5)
TX.MLT_RCAP[] = true; TX.MLT_HSDISC[] = true
TC.EXTRAP_FRAC[] = envf("EXTRAP", 0.3)
SM.REGION_CAP[] = true; SM.install_stage_hook!(mat); SM.install_dt_hook!(mat); get(ENV, "PREDICT", "1") == "1" && SM.install_state_predictor!()
TF.setscalar!(:dtime, envf("DT0", 1.0))
TD.STEPFN[] = TD.step!

if get(ENV, "TRACE", "0") == "1"       # Newton trace: per iteration the largest correction of each variable (zone index counted from the centre)
    oldmd = TD.ROUTINES[:maxdel]
    ixn = Dict(Int(TF.getscalar(s)) => string(s) for s in (:ir, :im, :id, :iu, :it, :ie, :if))
    TD.ROUTINES[:maxdel] = () -> begin
        oldmd()
        kx = TF.var(:kx); dx = TF.var(:dx)
        pr("   it %2d dt=%.3e dmax=%.3e  %s\n", Int(TF.getscalar(:iter)), TF.getscalar(:dtime), TF.getscalar(:dmax),
           join([@sprintf("%s:%+.2e@%d", ixn[i], dx[i], kx[i]) for i in 1:Int(TF.getscalar(:neqn))], " "))
    end
end
nstep = Ref(0); t0wall = time(); series = Tuple[]
Y0 = setup.Y0
function robust_step!()
    snap = TD.snapshot(); dt0 = TF.getscalar(:dtime)
    for attempt in 0:7
        TD.restore!(snap); TF.setscalar!(:dtime, dt0 * 0.3^attempt)
        res = nothing
        try
            res = SM.solve_step!(mat, TD.STEPFN[]; tol = envf("SPLITTOL", 3e-3), maxit = envi("SPLITIT", 4), ztol = envf("ZSPLITTOL", 3e-3))
        catch e
            e isa InterruptException && rethrow()
            pr("   step failed (%s): retry with dt x 0.3\n", sprint(showerror, e)[1:min(end, 100)])
            continue
        end
        if res[1] > envf("SPLITREJ", 0.02) || res[3] > envf("ZSPLITREJ", 1e-2)
            pr("   step rejected: source/network mismatch %.3f (zone %.1e), retry with dt x 0.3\n", res[1], res[3]); continue
        end
        attempt > 0 && (COOL[] = 12)
        return res
    end
    error("robust_step!: 8 attempts failed")
end

function after_step(k, xme_start, Mz_start)
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    sd = TX.step_diag(dt)
    r = SM.burn!(mat; sd)
    SM.mix!(mat; sd, capmix = false, bridge = envi("BRIDGE", 4), tmax = envf("MIXTMAX", 0.02))
    adv = SM.advect!(mat, xme_start, Mz_start; xin = Y0, order = envi("ADVORDER", 2))
    d = diag(sd)
    get(ENV, "DTWHY", "0") == "1" && pr("   DTWHY %d dt=%.3e s -> next %.3e s: max rel change %.4f (%s zone %d), iter %d, caps: dtcap %.3e predict %.3e; adv cfl %.2f x%d\n", k, dt, TF.getscalar(:dtime), DTW[][1], DTW[][3], DTW[][2], Int(TF.getscalar(:iter)), mat.dtcap, SM.predict_cap(mat), adv.cfl, adv.nsub)
    push!(series, (d.t / YR, dt / YR, d.M / Msun, d.R / Rsun, d.Lstar / Lsun, d.Lacc / Lsun, d.Teff, d.Tc, d.ρc, d.Tmax, d.mTmax, d.Ld / Lsun, d.Lpp / Lsun, d.Lcno / Lsun, d.fD, d.tKH / YR, d.tacc / YR, d.Sion, d.rph / Rsun, d.Lph / Lsun))
    if k % envi("EVERY", 10) == 0 || k <= 5
        pr("%5d t=%.4e yr dt=%.3e s M=%.4f R=%.2f Rsun L*=%.3f Lacc=%.2f Lph=%.2f Teff=%.0f Tc=%.4e rho_c=%.3e Tmax=%.3e@%.3f LD=%.2e fD=%.4f it=%d adv(cons %.1e cfl %.2f x%d dT %.1e) lim=%.3f %s@%d ksh=%d conv=%s\n", k, d.t / YR, dt, d.M / Msun, d.R / Rsun, d.Lstar / Lsun, d.Lacc / Lsun,
           d.Lph / Lsun, d.Teff, d.Tc, d.ρc, d.Tmax, d.mTmax, d.Ld / Lsun, d.fD, Int(TF.getscalar(:iter)), adv.cons, adv.cfl, adv.nsub, adv.dTmax, DTW[][1], DTW[][3], DTW[][2], d.ksh - a + 1, d.conv)
    end
    d
end

# ----------------------------------------------------------------------------------- checkpoints
using Serialization
ckdir = get(ENV, "CKDIR", "/tmp/s6")
function checkpoint(name)
    mkpath(ckdir)
    TF.save_state(joinpath(ckdir, "s6_$(name).tstate.gz"))
    serialize(joinpath(ckdir, "s6_$(name).aux"), (; Y = mat.Y, Yprev = mat.Yprev, fscale = mat.fscale, fslope = mat.fslope, Dmix = mat.Dmix, ledger = mat.ledger, dtcap = mat.dtcap, capE = mat.capE, series, nstep = nstep[]))
end
function restart!(name)
    TF.load_state!(joinpath(ckdir, "s6_$(name).tstate.gz"))
    aux = deserialize(joinpath(ckdir, "s6_$(name).aux"))
    mat.Y .= aux.Y; mat.Yprev .= aux.Yprev; mat.fscale .= aux.fscale; mat.fslope .= aux.fslope; mat.Dmix = aux.Dmix; mat.ledger = aux.ledger; mat.dtcap = aux.dtcap; mat.capE = aux.capE
    append!(series, aux.series); nstep[] = aux.nstep
    global dn, tn, rn, to, do_, ro = TF.var(:dn), TF.var(:tn), TF.var(:rn), TF.var(:to), TF.var(:do), TF.var(:ro)
    TF.setscalar!(:thet, 1.0); TF.setscalar!(:stol, 0.3); TF.setscalar!(:tmax, 1e30); TF.setscalar!(:niter, envi("NITER", 25))
    TX.REGULAR_CENTRE[] = true; TX.ALPHA_MLT[] = envf("ALPHA", 1.5)
    TD.EOS[](); TD.OPAC[](); SM.refresh!()
end
# ----------------------------------------------------------------------------------- run
haskey(ENV, "RESTART") && (restart!(ENV["RESTART"]); pr("restarted from %s at step %d\n", ENV["RESTART"], nstep[]))
let d = diag(TX.step_diag(1.0))
    pr("S6 start: M* = %.4f Msun, R* = %.2f Rsun, Mdot = %.3e Msun/yr, n = %d zones (star %d), tau_grid = %.2e s\n", d.M / Msun, d.R / Rsun, Mdot * YR / Msun, setup.n, setup.nst, envf("TAU", 1e7))
end
maxsteps = envi("MAXSTEPS", 100); tend = envf("TEND", 1e12) * YR; mend = envf("MEND", 1e9) * Msun
try
    for k in 1:maxsteps
        xme_start = copy(TF.var(:xmen)[a:b+1]); Mz_start = 4π .* [dn[kk] * TF.var(:dvoln)[kk] for kk in a:b]
        robust_step!()
        nstep[] += 1
        d = after_step(nstep[], xme_start, Mz_start)
        (nstep[] % envi("CKEVERY", 500) == 0) && haskey(ENV, "SAVE") && checkpoint(ENV["SAVE"] * "_auto" * string((nstep[] ÷ envi("CKEVERY", 500)) % 2))
        if d.M >= mend; pr("MEND reached\n"); break; end
        if TF.getscalar(:timen) > tend; pr("TEND reached\n"); break; end
    end
catch e
    pr("STOPPED at step %d: %s\n", nstep[], sprint(showerror, e, catch_backtrace())[1:min(end, 1200)])
end
let io = open(joinpath(@__DIR__, "..", "lib", "popiii", "s6_$(tag).csv"), "w")
    println(io, "t_yr,dt_yr,M,R_Rsun,Lstar_Lsun,Lacc_Lsun,Teff,Tc,rho_c,Tmax,mTmax,LD,Lpp,Lcno,fD,tKH_yr,tacc_yr,Sion,rph_Rsun,Lph_Lsun")
    for row in series; println(io, join(row, ",")); end
    close(io)
end
haskey(ENV, "SAVE") && checkpoint(ENV["SAVE"])
pr("wall %.1f s\n", time() - t0wall)
