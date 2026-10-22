# S5: isolated Pop III star (Z = 0) from an n = 3 polytrope (Eddington standard model) through the pre-main sequence and the main sequence,
# with the energy ledger. Phases: (1) relaxation of the initial model with burning off (backward Euler, dt ≤ 3 t_dyn, until max|u|/c_s < 1e-3
# after ≥ 30 t_dyn), (2) evolution with burning from the first step (solve_step!: consistent source/burn coupling, robust_step!: step rejection,
# DT_RULE: step control on the relative change of T, ρ, r), stopping at X_c < XSTOP. Output: lib/popiii/s5_<OUT>.log and .csv.
#   TITAN_NOLIB=1 julia -t 4 --project=env s5_run.jl        configuration through environment variables (default):
#   star:       M [100] Msun, X [0.75], R0 [6.5] Rsun (initial polytrope radius; ≈ 1.3-1.5 × ZAMS), N [140] zones, NENV [60] log-spaced envelope zones,
#               QSPLIT [0.05] exterior mass fraction where they start, MFLOOR [1e-10] exterior mass fraction of the surface zone, ALPHA [1.5] MLT,
#               ATM [1] grey-atmosphere offset of the initial T (see s4_star.jl), KCAP [Inf] soft κ cap (not needed)
#   physics:    FCOUPLING [1] level factor f on ε in the residual (0 = S4 coupling), BURNT [path|end] T history of the network integration,
#               EFRAC [0.5], DXMAX [2e-3] caps of dt on the energy / hydrogen change per step (S5 runs: 5, 4e-3)
#   numerics:   INT [trbdf2|theta] integrator of the evolution phase, TOL [0.03] target of max |Δ ln T, ρ, r| per step, GROW [1.6] dt growth per step,
#               CONV [1e-5] Newton tolerance (1e-7 … 1e-8 for smooth near-adiabatic ∇), NITER [25], SPLITTOL [3e-3] / SPLITIT [4] solve_step!,
#               SPLITREJ [0.02] / REJECT [0.15] step rejection, FAILTFAC [0.5] dt factor after a Newton failure, EXTRAP [0.3] max grid
#               extrapolation in zone widths, PREDICT [1] T/E_r time extrapolation, LAGRHO [1], RCAP [1] / HSDISC [1] MLT options (TitanResidFast)
#   run:        MAXSTEPS, TEND [yr], XSTOP [1e-3], EVERY [20] log interval, PROFILE_X [list of X_c] profile dumps, RELAX [1], STOPZAMS [0],
#               SAVE [name] checkpoint in /tmp/s5 at the end (+ every CKEVERY [400] steps as <name>_auto), RESTART [name]
#   debugging:  DETAIL, DTWHY, NORULE, TRACE, ZONEDUMP [n], RESDUMP, JACTEST, ITERTEST, DUMPNOCORE, NANREPORT, VSPLIT
include(joinpath(@__DIR__, "s4_star.jl"))
include(joinpath(@__DIR__, "Trbdf2.jl"))
envf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
envi(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
# E2E-A: performance recording (TitanDriver.Perf): always accumulated, written to lib/popiii/perf/<PERF>_* when PERF=<name> is set (PERFSTAGE labels the stage)
const Perf = TD.Perf; const PERFNAME = get(ENV, "PERF", "")
Perf.reset!(isempty(PERFNAME) ? "s5" : PERFNAME; stage = get(ENV, "PERFSTAGE", "ms"), driver = "s5_run.jl", out = get(ENV, "OUT", "run"), M = get(ENV, "M", "100"), n = get(ENV, "N", "140"),
            integ = get(ENV, "INT", "trbdf2"), tol = get(ENV, "TOL", "0.03"), handoff = get(ENV, "HANDOFF", ""))
Mstar = envf("M", 100.0) * Msun; Rs = envf("R0", 6.5) * Rsun; X0 = envf("X", 0.75)
n = envi("N", 140); nenv = envi("NENV", 60); yr = 3.15576e7
integ = get(ENV, "INT", "trbdf2")
tol = envf("TOL", 0.03)
# E2E-A: HANDOFF=<file written by s6b_run.jl HANDOFF_OUT> starts from the accreted star (its grid, profiles and composition) instead of the polytrope
include(joinpath(@__DIR__, "e2e_handoff.jl"))
const HANDOFF = get(ENV, "HANDOFF", "")
st = !isempty(HANDOFF) ? star_from_handoff(HANDOFF; alpha = envf("ALPHA", 1.5), efrac = envf("EFRAC", 0.5), dXmax = envf("DXMAX", 2e-3), burnT = Symbol(get(ENV, "BURNT", "path")),
                                           kcap = envf("KCAP", Inf), fcoupling = get(ENV, "FCOUPLING", "1") == "1") :
     build_star(; M = Mstar, R = Rs, X = X0, n, nenv, nuc = true, alpha = envf("ALPHA", 1.5), efrac = envf("EFRAC", 0.5), dXmax = envf("DXMAX", 2e-3),
                mfloor = envf("MFLOOR", 1e-10), qsplit = envf("QSPLIT", 0.05), kcap = envf("KCAP", Inf), burnT = Symbol(get(ENV, "BURNT", "path")),
                fcoupling = get(ENV, "FCOUPLING", "1") == "1", atm = get(ENV, "ATM", "1") == "1")
mat = st.mat
mat.dXmax = envf("DXMAX", 2e-3)
# E2E-A: EDDCAP = Γ_max (default 0 = off): soft cap of the opacity at κ_cap = Γ_max 4π c G M / L_surface, updated every step (see notes E2E-A: a very massive star at L/L_Edd,es ≈ 0.8 has
# Γ > 1 in the He II opacity bump and its 1D hydrostatic envelope is unbound; mass loss is off, so the opacity is limited instead)
const EDDCAP = envf("EDDCAP", 0.0)
update_eddcap!(L) = EDDCAP > 0 && (mat.kcap = EDDCAP * 4π * Z0.c * Z0.G * st.M / L)
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
tag = get(ENV, "OUT", "run"); logio = open(joinpath(@__DIR__, "..", "lib", "popiii", "s5_$(tag).log"), "w")
pr(f, args...) = (s = Printf.format(Printf.Format(f), args...); print(s); print(logio, s); flush(logio))
TF.setscalar!(:thet, 1.0)
TD.STEPFN[] = () -> (phase[] === :evol && integ == "trbdf2" ? Trbdf2.step!() : TD.step!())    # relaxation always with backward Euler (damps the dynamical modes)
function apply_settings!()      # (a restart overwrites the scalars of the state file)
    TF.setscalar!(:stol, envf("STOL", 0.3)); TF.setscalar!(:tmax, 1e30); TF.setscalar!(:niter, envi("NITER", 25)); TF.setscalar!(:thet, 1.0)
    haskey(ENV, "CONV") && TF.setscalar!(:conv, parse(Float64, ENV["CONV"]))
end
apply_settings!()
dn, tn, rn, to, do_, ro = TF.var(:dn), TF.var(:tn), TF.var(:rn), TF.var(:to), TF.var(:do), TF.var(:ro)
pr("S5 star: M = %.1f Msun, X = %.3f, R0 = %.3f Rsun, n = %d (env %d, qsplit %.3g), integrator %s, theta %.2f, alpha_MLT %.2f, tol %.3g, dXmax %.3g, efrac %.3g, fcoupling %s, burnT %s, kcap %.3g\n",
   st.M / Msun, X0, Rs / Rsun, n, nenv, envf("QSPLIT", 0.05), integ, TF.getscalar(:thet), TX.ALPHA_MLT[], tol, mat.dXmax, mat.efrac, mat.fcoupling, mat.burnT, mat.kcap)
pr("polytrope: beta = %.4f, rho_c = %.3f, T_c = %.4e K\n", st.β, st.ρc, st.Tc)

phase = Ref(:relax)
update_eddcap!(star_diag(st).Ls)
# ----------------------------------------------------------------------------------- energy ledger
G = Z0.G
"Total energy components of the current state [erg]: internal+radiation, kinetic, gravitational."
function energies()
    egrn, un, dvn = TF.var(:egrn), TF.var(:un), TF.var(:dvoln)
    Eint = 0.0; Ekin = 0.0; Eg = 0.0
    mcum = 0.0
    dmz = [4π * dn[k] * dvn[k] for k in a:b]
    for k in a:b
        j = k - a + 1
        Eint += dmz[j] * egrn[k]; Ekin += dmz[j] * 0.25 * (un[k]^2 + un[k+1]^2)
        mcum += dmz[j]                                  # mass interior to interface k+1
        Mk = k < b ? 0.5 * (dmz[j] + dmz[j+1]) : 0.5 * dmz[j]
        Eg -= G * mcum * Mk / rn[k+1]
    end
    (; Eint, Ekin, Eg, Etot = Eint + Ekin + Eg)
end
mutable struct Ledger
    E_prev::Float64; Ebur::Float64; Fout::Float64; Einj::Float64; err_hyd::Float64; err_hyd_abs::Float64; E0::Float64; t0::Float64; Ebur_abs::Float64
    corr::Float64
end
led = Ledger(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
function conv_core()      # mass fraction of the convective core: outermost convective interface of the first convective region (gaps of ≤ 2 interfaces tolerated)
    D = mat.Dmix; j = a + 1; last = a; gap = 0
    while j <= b && gap <= 2
        D[j] > 0 ? (last = j; gap = 0) : (gap += 1); j += 1
    end
    sum((4π * dn[k] * TF.var(:dvoln)[k] for k in a:last-1); init = 0.0) / st.M
end
Xc() = mat.Y[SM.IH, a]
Tdyn() = sqrt(rn[b+1]^3 / (G * st.M))

# ----------------------------------------------------------------------------------- time-step control
nmv = ["r", "m", "rho", "u", "T", "E", "F"]
const DTWHAT = Ref{Any}((0.0, ("", 0))); const COOL = Ref(0)
function dt_rule(dt0)
    if phase[] === :relax
        return min(1.5 * dt0, 3 * Tdyn())
    end
    d = 0.0; dwhat = ("", 0)
    @inbounds for k in a:b
        x = abs(log(tn[k] / to[k])); x > d && (d = x; dwhat = ("T", k - a + 1))
        x = abs(log(dn[k] / do_[k])); x > d && (d = x; dwhat = ("rho", k - a + 1))
    end
    @inbounds for k in a+1:b+1; x = abs(log(rn[k] / ro[k])); x > d && (d = x; dwhat = ("r", k - a + 1)); end
    DTWHAT[] = (d, dwhat)
    gmax = COOL[] > 0 ? 1.15 : envf("GROW", 1.6); COOL[] = max(COOL[] - 1, 0)       # after a rejected / failed step grow cautiously (AIMD)
    f = clamp(0.9 * tol / max(d, 1e-12), 0.5, gmax)
    TF.getscalar(:iter) > envi("ITSLOW", 20) && (f = min(f, 0.8))
    f * dt0
end
TC.DT_RULE[] = get(ENV, "NORULE", "0") == "1" ? nothing : dt_rule
get(ENV, "NOCONV", "0") == "1" && (TX.CONVECTION[] = false)
TD.FAIL_TFAC[] = envf("FAILTFAC", 0.5); Trbdf2.TFAC[] = envf("FAILTFAC", 0.5)
TC.PRED_LAGRHO[] = get(ENV, "LAGRHO", "1") == "1"
TX.MLT_RCAP[] = get(ENV, "RCAP", "1") == "1"; TX.MLT_HSDISC[] = get(ENV, "HSDISC", "1") == "1"
SM.VERBOSE_SPLIT[] = get(ENV, "VSPLIT", "0") == "1"
TC.EXTRAP_FRAC[] = envf("EXTRAP", 0.3)
SM.REGION_CAP[] = get(ENV, "REGIONCAP", "1") == "1"; SM.install_stage_hook!(mat); SM.install_dt_hook!(mat); get(ENV, "PREDICT", "1") == "1" && SM.install_state_predictor!()
TF.setscalar!(:dtime, envf("DT0", 1e2))
mat.nuc = false

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
# ----------------------------------------------------------------------------------- per-step bookkeeping
nstep = Ref(0); t0wall = time(); tstart = Ref(0.0); lastsd = Ref{Any}(nothing)
series = Tuple[]; win = Tuple{Float64,Float64,Float64}[]; zams = Ref{Any}(nothing)
lastmism = (NaN, 0, NaN)
pxs = haskey(ENV, "PROFILE_X") ? sort(parse.(Float64, split(ENV["PROFILE_X"], ",")); rev = true) : Float64[]
function profile(io = nothing)
    rnv, un, frn, ern, pgn, xkn = (TF.var(s) for s in (:rn, :un, :frn, :ern, :pgn, :chifn))
    sd = lastsd[]
    pr("%4s %10s %10s %10s %10s %9s %8s %9s %10s %10s %8s %8s %8s %9s\n", "zone", "m_ext/M", "r/Rsun", "rho", "T", "kappa_R", "beta", "u [km/s]", "Lrad/Lsun", "Lconv/Lsun", "grad", "grad_ad", "Gamma", "X_H")
    msum = cumsum([4π * dn[k] * TF.var(:dvoln)[k] for k in a:b]); Mt = msum[end]
    for k in unique(sort([a:4:b; b-30:b]))
        k < a && continue
        j = k - a + 1; Pr = ern[k] / 3; L = 4π * rnv[k+1]^2 * frn[k+1]; mi = msum[j]
        Lc = 4π * rnv[k+1]^2 * sd.fcv[k+1]
        Γ = xkn[k] * (L + Lc) / (4π * Z0.c * Z0.G * mi)
        Pm = pgn[k] + Pr; ∇ = k < b ? log(tn[k+1] / tn[k]) / log((pgn[k+1] + ern[k+1] / 3) / Pm) : NaN
        TX.ZONE[] = k; ∇ad = TX.material_thermo(mat, dn[k], tn[k])[1]
        pr("%4d %10.3e %10.4f %10.3e %10.4e %9.4f %8.4f %9.3f %10.3e %10.3e %8.4f %8.4f %8.4f %9.5f\n", j, 1 - mi / Mt, rnv[k+1] / Rsun, dn[k], tn[k], xkn[k], pgn[k] / (pgn[k] + Pr), un[k+1] / 1e5, L / Lsun, Lc / Lsun, ∇, ∇ad, Γ, SM.hmass(mat, k))
    end
end

function after_step(k)
    nstep[] = k
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    acc = mat.acc
    if haskey(ENV, "ITERTEST") && k == envi("ITERSTEP", 1)      # extra Newton iterations on the converged state: does the step satisfy its equations?
        for it in 1:envi("ITERTEST", 4)
            TD.r(:matgen); TD.r(:matslv); TD.r(:maxdel)
            kx = TF.var(:kx); dx = TF.var(:dx)
            pr("extra Newton it %d: dmax = %.3e  %s\n", it, TF.getscalar(:dmax), join([@sprintf("%.2e@%d", dx[i], kx[i]) for i in 1:Int(TF.getscalar(:neqn))], " "))
            TD.r(:update)
        end
    end
    Eh = energies()                                   # after the hydro step, before burn!/mix!
    sd = TX.step_diag(dt); lastsd[] = sd
    Tc_h = tn[a]; fs_old = mat.fcur[a]; Ls_ = star_diag(st).Ls
    if get(ENV, "DETAIL", "0") == "1" && phase[] === :evol      # central-zone comparison of the three estimates of the energy release
        Y0c = copy(mat.Y[:, a]); mkc = 4π * dn[a] * TF.var(:dvoln)[a]
        e_res = SM.eps_nuc(dn[a], tn[a], Y0c) * dt * mkc               # residual-like: eps(T_n, Y_old) dt
        rb = NN.burn(Y0c, dn[a], tn[a], dt); e_bconst = rb.ε * dt * mkc   # network at constant T_n
        pr("   CENTRE dt=%.3e yr  T_o=%.6e T_n=%.6e  Einj0=%.5e (booked/f)  eps(T_n,Yold)dt m=%.5e  burn(T_n) dt m=%.5e  ratio burn/eps=%.4f  Einj0/eps=%.4f\n", dt / yr, to[a], tn[a], acc.Einj[a] / fs_old, e_res, e_bconst, e_bconst / e_res, acc.Einj[a] / fs_old / e_res)
    end
    if phase[] === :evol
        r = Perf.@section :burn SM.burn!(mat; sd)
        mixr = Perf.@section :mix SM.mix!(mat; sd, capmix = false, bridge = envi("BRIDGE", 4), tmax = envf("MIXTMAX", 0.02))
        Ea = energies()
        get(ENV, "NORULE", "0") == "1" && pr("   TITAN timstp step %d dt=%.3e yr -> next %.3e yr: smax %.3e (stol %.2f) limited by %s\n", k, dt / yr, TF.getscalar(:dtime) / yr, TF.getscalar(:smax), TF.getscalar(:stol),
            join(["$(nmv[i]):zone $(TF.var(:kx)[i]) dx=$(round(TF.var(:dx)[i]; sigdigits = 2))" for i in 1:Int(TF.getscalar(:neqn))], ", "))
        get(ENV, "DTWHY", "0") == "1" && pr("   DTWHY step %d dt=%.3e yr: max rel change %.3e at %s zone %d, newton iter %d; caps [yr]: dX %.3e  E(spec) %.3e  E(eps) %.3e  predict %.3e\n", k, dt / yr, DTWHAT[][1], DTWHAT[][2][1], DTWHAT[][2][2], TF.getscalar(:iter),
            dt * mat.dXmax / max(r.dXH, 1e-300) / yr, dt * mat.efrac / max(r.dEsp, 1e-300) / yr, dt * mat.efrac / max(r.dEfrac, 1e-300) / yr, SM.predict_cap(mat) / yr)
        get(ENV, "DETAIL", "0") == "1" && pr("   DETAIL step %d dt=%.3e yr kind=%s it=%d split(%d it, mism %.1e) Tc_hyd=%.6e Tc_end=%.6e f_c(used)=%.4f f_c(next)=%.4f Einj/dt/L=%.4f Ebur/dt/L=%.4f Fout/dt/L=%.4f Lnuc_end/L=%.4f\n", k, dt / yr, string(acc.last), TF.getscalar(:iter), lastmism[2], lastmism[1],
            Tc_h, tn[a], fs_old, mat.fslope[a] * dt, sum(acc.Einj[a:b]) / dt / Ls_, r.Ebur / dt / Ls_, acc.Fout / dt / Ls_, SM.nuclear_power(mat) / Ls_)
        # ledger of this step: hydro step vs booked injection and boundary loss; burn! adds Ebur − Einj
        Einj = sum(acc.Einj[a:b]); Fout = acc.Fout
        err = (Eh.Etot - led.E_prev) - (Einj - Fout)
        led.err_hyd += err; led.err_hyd_abs += abs(err)
        led.Ebur += r.Ebur; led.Fout += Fout; led.Einj += r.Einj; led.Ebur_abs += abs(r.Ebur); led.corr += abs(r.Ebur - r.Einj)
        led.E_prev = Ea.Etot
        step_info = (; Ebur = r.Ebur, Fout, Einj = r.Einj, err, dEtot = Ea.Etot - Eh.Etot)
    else
        Ea = Eh; led.E_prev = Eh.Etot; step_info = (; Ebur = 0.0, Fout = acc.Fout, Einj = 0.0, err = 0.0, dEtot = 0.0)
    end
    d = star_diag(st)
    update_eddcap!(d.Ls)
    ts = (d.t - tstart[]) / yr
    Yc = view(mat.Y, :, a)
    if get(ENV, "DUMPNOCORE", "0") == "1" && phase[] === :evol && sum(mat.Dmix[a:a+30] .> 0) < 20      # the convective core is not detected at this step
        pr("   NOCORE step %d (hydro state before burn!/mix!): j sup fcv vcv lcv T(hyd)\n", k)
        for j in 1:12
            kk = a + j; pr("      %2d  %11.4e  %10.3e  %10.3e  %10.3e  %.8e\n", j, sd.nab[kk] - sd.nad[kk], sd.fcv[kk], sd.vcv[kk], sd.lcv[kk], tn[kk-1])
        end
    end
    if phase[] === :evol
        # ledger closure over the evolution so far: Δ(E_int+E_kin+E_grav) = ΣE_bur − ∫Φ_out dt
        dEt = Ea.Etot - led.E0
        clos = (dEt - (led.Ebur - led.Fout)) / max(led.Ebur_abs, 1e-300)
        eg = NN.energy(Yc, NN.rates_k(tn[a], dn[a], NN.Plasma(dn[a], Yc), NN.ye_of(Yc)))[3]; etot = max(eg.pp + eg.d + eg.cno + eg.he, 1e-300)
        push!(series, (ts, dt / yr, d.R / Rsun, log10(d.Ls / Lsun), log10(d.Teff), d.Tc, d.ρc, Xc(), conv_core(), SM.xyz(Yc)[3], eg.cno / etot,
                       step_info.Ebur / max(dt, 1e-30) / Lsun, step_info.Fout / max(dt, 1e-30) / Lsun, clos, lastmism[1], TF.getscalar(:iter)))
        push!(win, (step_info.Ebur, step_info.Fout, dt)); length(win) > 10 && popfirst!(win)
        wr = sum(w[2] for w in win) / max(sum(w[1] for w in win), 1e-300)          # L_surface / <L_nuc> over the last ≤ 10 steps
        if zams[] === nothing && length(win) >= 10 && abs(wr - 1) < 0.02 && ts > 100 * 3.15e3 / yr
            zams[] = (; ts, R = d.R / Rsun, logL = log10(d.Ls / Lsun), logT = log10(d.Teff), Tc = d.Tc, ρc = d.ρc, Xc = Xc(), Mcc = conv_core(), step = k, XCNO = SM.xyz(Yc)[3])
            pr("-- ZAMS (L_surface/<L_nuc> = %.4f over the last 10 steps): t = %.4e yr, R = %.4f Rsun, log L = %.4f, log Teff = %.5f, T_c = %.4e K, rho_c = %.3f, X_c = %.5f, M_cc/M = %.4f, X_CNO,c = %.2e\n",
               wr, ts, zams[].R, zams[].logL, zams[].logT, zams[].Tc, zams[].ρc, zams[].Xc, zams[].Mcc, zams[].XCNO)
        end
    end
    if k % envi("EVERY", 20) == 0 || k <= 3
        eg = NN.energy(Yc, NN.rates_k(tn[a], dn[a], NN.Plasma(dn[a], Yc), NN.ye_of(Yc)))[3]; etot = max(eg.pp + eg.d + eg.cno + eg.he, 1e-300)
        if phase[] === :relax
            umax = maximum(abs.(TF.var(:un)[a:b+1]) ./ [sqrt(5 / 3 * TF.var(:pgn)[min(kk, b)] / dn[min(kk, b)]) for kk in a:b+1])
            pr("relax %5d t=%.4e s dt=%.3e s R=%.4f logL=%.4f logTeff=%.5f Tc=%.4e rho_c=%.3f max|u|/cs=%.2e iter=%d\n", k, d.t, dt, d.R / Rsun, log10(d.Ls / Lsun), log10(d.Teff), d.Tc, d.ρc, umax, TF.getscalar(:iter))
        else
            dEt = Ea.Etot - led.E0
            wr = sum(w[2] for w in win) / max(sum(w[1] for w in win), 1e-300)
            pr("%6d %10.4e %9.3e %8.4f %8.4f %9.5f %9.5f %9.4f %8.4f %9.5f %8.5f %9.2e %3d %d | %9.1e %9.2e %9.2e %7.4f %9.3e\n", k, ts, dt / yr, d.R / Rsun, log10(d.Ls / Lsun), log10(d.Teff), d.Tc / 1e8, d.ρc,
               wr, Xc(), conv_core(), SM.xyz(Yc)[3], TF.getscalar(:iter), lastmism[2], (dEt - (led.Ebur - led.Fout)) / max(led.Ebur_abs, 1e-300),
               led.corr / max(led.Ebur_abs, 1e-300), led.err_hyd_abs / max(led.Ebur_abs, 1e-300), eg.cno / etot, mat.dtcap / yr)
        end
    end
    if phase[] === :evol && !isempty(pxs) && Xc() <= pxs[1]
        pr("-- profile at X_c = %.4f (t = %.4e yr, step %d)\n", Xc(), ts, k); profile(); popfirst!(pxs)
    end
end

using Serialization
ckdir = get(ENV, "CKDIR", "/tmp/s5")
function checkpoint(name)
    mkpath(ckdir)
    TF.save_state(joinpath(ckdir, "s5_$(name).tstate.gz"))
    serialize(joinpath(ckdir, "s5_$(name).aux"), (; Y = mat.Y, Yprev = mat.Yprev, fscale = mat.fscale, fslope = mat.fslope, Dmix = mat.Dmix, ledger = mat.ledger, dtcap = mat.dtcap, capE = mat.capE,
              led = (led.E_prev, led.Ebur, led.Fout, led.Einj, led.err_hyd, led.err_hyd_abs, led.E0, led.t0, led.Ebur_abs, led.corr), tstart = tstart[], series, nstep = nstep[], nrelax = nrelax[], pxs, zams = zams[], win))
end
function restart!(name)
    TF.load_state!(joinpath(ckdir, "s5_$(name).tstate.gz"))
    aux = deserialize(joinpath(ckdir, "s5_$(name).aux"))
    mat.Y .= aux.Y; mat.Yprev .= aux.Yprev; mat.fscale .= aux.fscale; hasproperty(aux, :fslope) && (mat.fslope .= aux.fslope); mat.Dmix = aux.Dmix; mat.ledger = aux.ledger; mat.dtcap = aux.dtcap; mat.capE = aux.capE
    (led.E_prev, led.Ebur, led.Fout, led.Einj, led.err_hyd, led.err_hyd_abs, led.E0, led.t0, led.Ebur_abs, led.corr) = aux.led
    tstart[] = aux.tstart; append!(series, aux.series); nstep[] = aux.nstep; nrelax[] = aux.nrelax; empty!(pxs); append!(pxs, aux.pxs); hasproperty(aux, :zams) && (zams[] = aux.zams; append!(win, aux.win))
    global dn, tn, rn, to, do_, ro = TF.var(:dn), TF.var(:tn), TF.var(:rn), TF.var(:to), TF.var(:do), TF.var(:ro)      # load_state! replaces the arrays
    apply_settings!(); TD.EOS[](); TD.OPAC[](); SM.refresh!()
end
"""
One evolution step with rejection: the step is repeated with Δ × 0.3 from the saved state if Newton fails (TITAN's own retries cut Δ by 20 per
attempt and give up after 6) or if the converged state is not the neighbour of the old one (max |Δ ln T|, |Δ ln ρ| over the zones with
exterior mass fraction > 1e-8 above REJECT, default 0.15: a spurious Newton solution).
"""
function robust_step!()
    snap = TD.snapshot(); dt0 = TF.getscalar(:dtime); tol_ = envf("SPLITTOL", 3e-3); maxit_ = envi("SPLITIT", 4); rej = envf("REJECT", 0.15)
    dprev = Inf
    for attempt in 0:7
        TD.restore!(snap); TF.setscalar!(:dtime, dt0 * 0.3^attempt)
        res = nothing
        try
            res = Perf.@section(:solve_step, SM.solve_step!(mat, TD.STEPFN[]; tol = tol_, maxit = maxit_, ztol = envf("ZSPLITTOL", 3e-3)))
        catch e
            e isa InterruptException && rethrow()
            Perf.count!(:step_fail)
            pr("   step failed (%s): retry with dt x 0.3\n", sprint(showerror, e)[1:min(end, 80)])
            continue
        end
        if res[1] > envf("SPLITREJ", 0.02) || res[3] > envf("ZSPLITREJ", 1e-2)       # the heating the residual used does not match the network's release: do not hide it in the post-step correction
            Perf.count!(:step_reject)
            pr("   step rejected: source/network mismatch %.3f (zone-wise %.1e of the zone energy) after %d iterations at dt = %.3e yr, retry with dt x 0.3\n", res[1], res[3], res[2], (TF.getscalar(:timen) - TF.getscalar(:timeo)) / yr)
            continue
        end
        get(ENV, "DBGTO", "0") == "1" && pr("   DBGTO attempt %d: to[a]=%.5e tn[a]=%.5e to[a+1]=%.5e tn[a+1]=%.5e timeo=%.6e timen=%.6e\n", attempt, to[a], tn[a], to[a+1], tn[a+1], TF.getscalar(:timeo), TF.getscalar(:timen))
        d = 0.0; dk = 0; dv = ""
        for kk in a:b-12
            x = abs(log(tn[kk] / to[kk])); x > d && (d = x; dk = kk - a + 1; dv = "T")
            x = abs(log(dn[kk] / do_[kk])); x > d && (d = x; dk = kk - a + 1; dv = "rho")
        end
        if d > rej && !(attempt >= 2 && abs(d - dprev) < 0.05 * d)     # a change that does not shrink with dt is a property of the state, not of the step
            Perf.count!(:step_reject)
            pr("   step rejected: max relative change %.3f (%s, zone %d) > %.2f at dt = %.3e yr, retry with dt x 0.3\n", d, dv, dk, rej, (TF.getscalar(:timen) - TF.getscalar(:timeo)) / yr)
            if get(ENV, "DUMPREJ", "0") == "1"
                sdx = TX.step_diag(TF.getscalar(:timen) - TF.getscalar(:timeo))
                pr("      zone  T_old  T_new  rho_old  rho_new  sup_new  Fconv_new  X_H  (Dmix last)\n")
                for jj in max(dk - 5, 1):min(dk + 5, b - a)
                    kk = a + jj - 1
                    pr("      %4d %.5e %.5e %.4e %.4e %10.3e %10.3e %.5f %.3e\n", jj, to[kk], tn[kk], do_[kk], dn[kk], sdx.nab[kk+1] - sdx.nad[kk+1], sdx.fcv[kk+1], SM.hmass(mat, kk), mat.Dmix[kk+1])
                end
            end
            dprev = d
            continue
        end
        d > rej && pr("   step accepted although the relative change is %.3f (%s, zone %d): independent of dt (bifurcation / state-intrinsic)\n", d, dv, dk)
        attempt > 0 && (COOL[] = 12)
        return res
    end
    error("robust_step!: 8 attempts failed")
end
# ----------------------------------------------------------------------------------- relaxation (burning off), then evolution
nrelax = Ref(0)
restarted = haskey(ENV, "RESTART")
if restarted
    restart!(ENV["RESTART"]); phase[] = :evol; mat.nuc = true; pr("restarted from %s at step %d, t = %.5e yr\n", ENV["RESTART"], nstep[], (TF.getscalar(:timen) - tstart[]) / yr)
elseif get(ENV, "RELAX", "1") == "1"
    pr("== relaxation of the initial model (burning off, dt = 1.3x per step up to 3 t_dyn = %.2e s, t_dyn = %.2e s)\n", 3Tdyn(), Tdyn())
    TF.setscalar!(:dtime, 1.0)
    phase[] = :relax
    uprev = Inf
    try
        for k in 1:envi("RELAXSTEPS", 400)
            TD.run_julia!(maxsteps = 1, output = nothing); after_step(k); nrelax[] += 1
            umax = maximum(abs.(TF.var(:un)[a:b+1]) ./ [sqrt(5 / 3 * TF.var(:pgn)[min(kk, b)] / dn[min(kk, b)]) for kk in a:b+1])
            if TF.getscalar(:timen) > 30 * Tdyn() && umax < envf("UTOL", 1e-3)
                pr("relaxed at step %d: t = %.3e s = %.1f t_dyn, max|u|/c_s = %.2e\n", k, TF.getscalar(:timen), TF.getscalar(:timen) / Tdyn(), umax); break
            end
        end
    catch e
        pr("relaxation STOPPED at step %d: %s\n", nstep[], sprint(showerror, e)[1:min(end, 300)])
        rethrow()
    end
end
restarted || get(ENV, "PROFILE_RELAX", "0") != "1" || profile()
let d = star_diag(st)
    pr("start of evolution: R = %.4f Rsun, log L = %.4f, log Teff = %.5f, Tc = %.4e, rho_c = %.3f, L_nuc(X0)/L = %.4f\n", d.R / Rsun, log10(d.Ls / Lsun), log10(d.Teff), d.Tc, d.ρc, d.Lnuc / d.Ls)
end
if !restarted
    tstart[] = TF.getscalar(:timen); phase[] = :evol; mat.nuc = true
    let Ea = energies(); led.E0 = Ea.Etot; led.E_prev = Ea.Etot; end
    TF.setscalar!(:dtime, envf("DT1", 1e3 * 3.15e0))      # first evolution step: ~3e3 s, grows with the controller
end
pr("%6s %10s %9s %8s %8s %9s %9s %9s %8s %9s %8s %9s %3s %s | %9s %9s %9s %7s %9s\n", "step", "t [yr]", "dt [yr]", "R [Rsun]", "logL", "logTeff", "Tc [1e8]", "rho_c", "Ls/<Ln>", "X_c", "Mcc/M", "X_CNO,c", "it", "s", "closure", "corr", "err_hyd", "eCNO/e", "cap [yr]")
xstop = envf("XSTOP", 1e-3); tend = envf("TEND", 1e12) * yr; maxsteps = envi("MAXSTEPS", 3000)
try
    for k in 1:maxsteps
        Perf.@section :step begin
            global lastmism = Perf.@section(:robust_step, robust_step!()); Perf.@section(:after_step, after_step(nrelax[] + k))
        end
        Perf.count!(:steps)
        k % 50 == 0 && Perf.sample!(step = nrelax[] + k, t_yr = (TF.getscalar(:timen) - tstart[]) / yr, dt_yr = TF.getscalar(:dtime) / yr, zones = b - a + 1, newton = Int(TF.getscalar(:iter)), Xc = Xc(), wall_s = Perf.elapsed())
        k == envi("SAVE_AT_STEP", -1) && haskey(ENV, "SAVE") && checkpoint(ENV["SAVE"])
        (k % envi("CKEVERY", 400) == 0) && haskey(ENV, "SAVE") && checkpoint(ENV["SAVE"] * "_auto" * string((k ÷ envi("CKEVERY", 400)) % 2))
        if k == envi("DEBUGSTEP", -1)
            TF.setscalar!(:dtime, envf("DEBUGDT", 1.6e9))
            TD.r(:peq)
            for s in (:rn, :dn, :tn, :ern, :frn, :un, :egn, :xmen, :drdt, :dmdt, :pgn, :chifn, :fedd)
                v = TF.var(s); bad = findall(!isfinite, v[a-2:b+3]); isempty(bad) || pr("  after peq: non-finite %s at zones %s\n", s, string(bad[1:min(end, 8)] .+ (a - 3)))
            end
            pr("after peq: dtime %.4e, max|drdt| %.3e\n", TF.getscalar(:dtime), maximum(abs.(TF.var(:drdt)[a:b+1])))
            TD.r(:matgen)
            rhs = TF.var(:rhs); pr("rhs finite: %s\n", all(isfinite, rhs[:, a-2:b+3]))
            for s in (:fedd, :fcv)
                v = TF.var(s); bad = findall(!isfinite, v[a-2:b+3]); isempty(bad) || pr("  after matgen: non-finite %s at zones %s\n", s, string(bad[1:min(end, 8)] .+ (a - 3)))
            end
            break
        end
        if get(ENV, "STOPZAMS", "0") == "1" && zams[] !== nothing && nstep[] >= zams[].step + 10; pr("stopping 10 steps after the ZAMS\n"); break; end
        if Xc() < xstop; pr("central hydrogen exhausted (X_c = %.2e) at step %d\n", Xc(), k); break; end
        if TF.getscalar(:timen) - tstart[] > tend; pr("TEND reached\n"); break; end
    end
catch e
    pr("STOPPED at step %d: %s\n", nstep[], sprint(showerror, e)[1:min(end, 300)])
    if get(ENV, "NANREPORT", "1") == "1"
        for s in (:rn, :dn, :tn, :ern, :frn, :un, :egn, :ro, :do, :to, :ero, :uo, :fro)
            v = TF.var(s); bad = findall(!isfinite, v[a-2:b+3]); isempty(bad) || pr("  non-finite %s at zones %s\n", s, string(bad[1:min(end, 8)] .+ (a - 3)))
        end
        for k in a:b
            TX.ZONE[] = k; ρk, Tk = do_[k], to[k]
            for (nm, f) in (("source", T -> TX.material_source(mat, ρk, T)), ("gas", T -> TX.material_gas(mat, ρk, T)[1]), ("e", T -> TX.material_gas(mat, ρk, T)[2]),
                            ("kappa", T -> TX.material_opac(mat, ρk, T)[1]), ("kappaP", T -> TX.material_opac(mat, ρk, T)[2]), ("nad", T -> TX.material_thermo(mat, ρk, T)[1]),
                            ("cp", T -> TX.material_thermo(mat, ρk, T)[2]))
                v = try ForwardDiff.derivative(f, Tk) catch ex; NaN end
                isfinite(v) || pr("  zone %d: d(%s)/dT = %s at rho=%.4e T=%.5e fscale=%.4g\n", k - a + 1, nm, string(v), ρk, Tk, mat.fscale[k])
            end
        end
        bad = findall(!isfinite, mat.fscale); isempty(bad) || pr("  non-finite fscale at %s\n", string(bad[1:min(end, 8)]))
        bad = findall(!isfinite, mat.Y); isempty(bad) || pr("  non-finite Y at %s\n", string(bad[1:min(end, 8)]))
    end
end
wall = time() - t0wall
d = star_diag(st)
if !isempty(PERFNAME)
    Perf.write(); Perf.summary()
end
pr("final: t = %.4e yr, R = %.4f Rsun, log L = %.4f, log Teff = %.5f, T_c = %.5e K, rho_c = %.4f, L_nuc/L = %.4f, X_c = %.5f, M_cc/M = %.4f; %d steps (%d relax), %.1f s wall\n",
   (d.t - tstart[]) / yr, d.R / Rsun, log10(d.Ls / Lsun), log10(d.Teff), d.Tc, d.ρc, d.Lnuc / d.Ls, Xc(), conv_core(), nstep[], nrelax[], wall)
pr("ledger: sum E_bur = %.4e erg, sum Fout dt = %.4e, E_tot change = %.4e; closure (dE - (Ebur - Fout))/sum|Ebur| = %.3e; splitting correction sum|Ebur-Einj|/sum|Ebur| = %.3e; hydro-step energy error sum|err|/sum|Ebur| = %.3e\n",
   led.Ebur, led.Fout, energies().Etot - led.E0, (energies().Etot - led.E0 - (led.Ebur - led.Fout)) / max(led.Ebur_abs, 1e-300), led.corr / max(led.Ebur_abs, 1e-300), led.err_hyd_abs / max(led.Ebur_abs, 1e-300))
if zams[] !== nothing
    z = zams[]; th = (d.t - tstart[]) / yr - z.ts
    pr("ZAMS: t = %.4e yr, R = %.4f Rsun, log L = %.4f, log Teff = %.5f, T_c = %.4e K, rho_c = %.3f, X_c = %.5f, M_cc/M = %.4f, X_CNO,c = %.2e (step %d)\n", z.ts, z.R, z.logL, z.logT, z.Tc, z.ρc, z.Xc, z.Mcc, z.XCNO, z.step)
    pr("H-burning lifetime from the ZAMS to X_c = %.1e: tau_H = %.4e yr (%.4e yr from the start of the evolution)\n", Xc(), th, (d.t - tstart[]) / yr)
end
let io = open(joinpath(@__DIR__, "..", "lib", "popiii", "s5_$(tag).csv"), "w")
    println(io, "t_yr,dt_yr,R_Rsun,logL,logTeff,Tc,rho_c,Xc,Mcc_over_M,XCNO_c,eCNO_over_e,Lnuc_avg_Lsun,Lsurf_avg_Lsun,closure,split_mismatch,newton_iter")
    for row in series; println(io, join(row, ",")); end
    close(io)
end
let sd = lastsd[]
    Lin = 0.0; vals = String[]
    for j in 1:4
        k = a + j - 1
        Lin += 4π * dn[k] * TF.var(:dvoln)[k] * SM.eps_nuc(dn[k], tn[k], view(mat.Y, :, k))
        push!(vals, @sprintf("%.3e/%.3e", 4π * rn[k+1]^2 * sd.fcv[k+1] / Lsun, Lin / Lsun))
    end
    pr("centre MLT: L_conv/L_nuc(<r) [Lsun] at the first 4 interfaces: %s;  L_conv(1)/L_conv(2) = %.3f (r^3 law: %.3f)\n", join(vals, "  "), sd.fcv[a+1] * rn[a+1]^2 / (sd.fcv[a+2] * rn[a+2]^2 + 1e-300), (rn[a+1] / rn[a+2])^3)
end
pr("integrator stats: %s\n", integ == "trbdf2" ? string(Trbdf2.STATS) : "theta")
get(ENV, "PROFILE", "0") == "1" && profile()
haskey(ENV, "SAVE") && checkpoint(ENV["SAVE"])
if haskey(ENV, "ZONEDUMP")       # diagnostic: zone table of the centre (convective flux, gradients, nuclear luminosity inside)
    nz = envi("ZONEDUMP", 30); sd = lastsd[]
    frn_, ern_, pgn_ = TF.var(:frn), TF.var(:ern), TF.var(:pgn)
    pr("zone dump (centre outwards): j  r[Rsun]  rho  T  lnP  grad  grad_ad  sup  Fconv  Frad  Lconv  Lrad  eps  Lnuc_inside  fscale\n")
    global Lin = 0.0
    for k in a:a+nz-1
        Pk = pgn_[k] + ern_[k] / 3; ε = SM.eps_nuc(dn[k], tn[k], view(mat.Y, :, k)); global Lin += 4π * dn[k] * TF.var(:dvoln)[k] * ε
        ∇ = sd.nab[k+1]; ∇ad = sd.nad[k+1]
        pr("%3d %8.4f %9.4e %10.5e %9.5f %9.6f %9.6f %10.3e %10.3e %10.3e %10.3e %10.3e %9.3e %10.3e %6.3f\n", k - a + 1, rn[k+1] / Rsun, dn[k], tn[k], log(Pk), ∇, ∇ad, ∇ - ∇ad, sd.fcv[k+1], frn_[k+1],
           4π * rn[k+1]^2 * sd.fcv[k+1] / Lsun, 4π * rn[k+1]^2 * frn_[k+1] / Lsun, ε, Lin / Lsun, mat.fscale[k])
    end
end

if haskey(ENV, "RESDUMP")      # residual rows of the energy equation at the converged state (zone, residual, source, convective flux divergence)
    dt_ = TF.getscalar(:timen) - TF.getscalar(:timeo)
    c = TX.context(; dt = dt_); R = zeros(c.neqn, c.nb); ws = TX.getws(Float64, c.mgr)
    TX.residual!(R, ws, c, zeros(c.neqn * c.nb)); W = ws.a
    it = c.ix[5]; ie_ = c.ix[6]; sd = lastsd[]
    pr("residual rows at the converged state (dt = %.4e s): zone  R_energy  R_rad  src  rmu*fcv(k+1)  rmu*fcv(k)  d(e rho dV)/dt\n", dt_)
    for k in a:a+7
        kbk = k - c.k0 + 1
        pr("%3d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n", k - a + 1, R[it, kbk], R[ie_, kbk], W.srcw[k], W.rmu[k+1] * W.fcv[k+1], W.rmu[k] * W.fcv[k],
           (W.dn[k] * W.egn[k] * W.dvoln[k] - c.do_[k] * c.ego[k] * c.dvolo[k]) / dt_)
    end
end

if haskey(ENV, "JACTEST")      # AD Jacobian vs central finite differences of the same residual, at the converged state of the last step
    dt_ = TF.getscalar(:timen) - TF.getscalar(:timeo)
    TD.r(:update)   # (idempotent) make sure derived arrays are current
    c = TX.context(; dt = dt_); neqn = c.neqn; nb = c.nb; k0 = c.k0; N = neqn * nb
    TX.assemble!(c)
    Jad = zeros(N, N)
    em2, em1, e00, ep1, ep2 = (TF.var(s) for s in (:em2, :em1, :e00, :ep1, :ep2))
    for kb in 1:nb, i in 1:neqn, o in -2:2
        cb = kb + o; (1 <= cb <= nb) || continue
        M = o == -2 ? em2 : o == -1 ? em1 : o == 0 ? e00 : o == 1 ? ep1 : ep2
        for j in 1:neqn; Jad[(kb-1)*neqn+i, (cb-1)*neqn+j] = M[i, j, k0+kb-1]; end
    end
    ws = TX.getws(Float64, c.mgr)
    function Rvec(x)
        R = zeros(neqn, nb); TX.residual!(R, ws, c, x); vec(R)
    end
    x0 = zeros(N); h = envf("JACH", 1e-7)
    Jfd = zeros(N, N)
    for col in 1:N
        xp = copy(x0); xm = copy(x0); xp[col] += h; xm[col] -= h
        Jfd[:, col] = (Rvec(xp) - Rvec(xm)) / (2h)
    end
    for kz in (1, 2, 3, 50)
        kb = kz + a - k0
        iT = ix_t = Int(TF.getscalar(:it))
        row = (kb - 1) * neqn + iT
        pr("   zone %d energy row: dR/dlnT_{k-1,k,k+1}: AD %s   FD %s\n", kz, join([@sprintf("%.4e", Jad[row, (kb-1+o)*neqn+iT]) for o in -1:1 if 1 <= kb+o <= nb], " "), join([@sprintf("%.4e", Jfd[row, (kb-1+o)*neqn+iT]) for o in -1:1 if 1 <= kb+o <= nb], " "))
    end
    names = ["r", "m", "rho", "u", "T", "E", "F"]
    ix = [Int(TF.getscalar(s)) for s in (:ir, :im, :id, :iu, :it, :ie, :if)]
    nm(i) = names[findfirst(==(i), ix)]
    pr("Jacobian check at the converged state (dt = %.4e s): AD vs central FD (h = 1e-7); entries with |ΔJ| > 1e-3 max|J row|:\n", dt_)
    bad = Tuple[]
    for row in 1:N
        scale = maximum(abs.(Jfd[row, :])) + 1e-300
        for col in 1:N
            d = abs(Jad[row, col] - Jfd[row, col])
            d > 1e-3 * scale && push!(bad, (d / scale, row, col, Jad[row, col], Jfd[row, col]))
        end
    end
    sort!(bad; rev = true)
    pr("%d entries differ by more than 1e-3 of their row scale\n", length(bad))
    for (r_, row, col, ja, jf) in bad[1:min(end, 25)]
        kr = (row - 1) ÷ neqn + 1; ir_ = (row - 1) % neqn + 1; kc = (col - 1) ÷ neqn + 1; ic = (col - 1) % neqn + 1
        pr("   row eq %-3s zone %3d  col var %-3s zone %3d:  AD %+.4e  FD %+.4e  rel.diff of row scale %.2e\n", nm(ir_), kr + k0 - 1 - a + 1, nm(ic), kc + k0 - 1 - a + 1, ja, jf, r_)
    end
end
