# S4 test d: a 100 Msun Eddington-model star in TITAN with radiation, MLT, nuclear burning and mixing, relaxed in time.
#   TITAN_NOLIB=1 julia -t 4 --project=env s4_d.jl           env: M (Msun), RSUN, N, NENV, TEND (s), MAXSTEPS, NUC, MIX, OUT (log tag), DTMAX
include(joinpath(@__DIR__, "s4_star.jl"))
using Serialization
Mstar = parse(Float64, get(ENV, "M", "100")) * Msun; Rs = parse(Float64, get(ENV, "RSUN", "5")) * Rsun
n = parse(Int, get(ENV, "N", "100")); nenv = parse(Int, get(ENV, "NENV", "24"))
nuc = get(ENV, "NUC", "1") == "1"; mixon = get(ENV, "MIX", "1") == "1"
st = build_star(; M = Mstar, R = Rs, n, nenv, nuc, alpha = parse(Float64, get(ENV, "ALPHA", "1.5")),
                efrac = parse(Float64, get(ENV, "EFRAC", "0.05")), mfloor = parse(Float64, get(ENV, "MFLOOR", "1e-10")),
                kcap = haskey(ENV, "KCAP") ? parse(Float64, ENV["KCAP"]) : Inf, theta = haskey(ENV, "THETA") ? parse(Float64, ENV["THETA"]) : nothing)
mat = st.mat; SM.install_dt_hook!(mat)
a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
tag = get(ENV, "OUT", "d"); logio = open(joinpath(@__DIR__, "..", "lib", "popiii", "s4_$(tag).log"), "w")
pr(f, args...) = (s = Printf.format(Printf.Format(f), args...); print(s); print(logio, s); flush(logio))
pr("Eddington model: M = %.1f Msun, R = %.3f Rsun, mu = %.4f, beta = %.4f, rho_c = %.3f, T_c = %.4e K (polytrope), M_poly/M = %.6f\n",
   st.M / Msun, st.R / Rsun, st.μmol, st.β, st.ρc, st.Tc, st.Mpoly / st.M)
pr("zones %d (env %d), theta = %.2f, stol = %.3g, nuc %s, mixing %s, alpha_MLT %.2f, efrac %.3g\n", n, nenv, TF.getscalar(:thet), TF.getscalar(:stol), nuc, mixon, TX.ALPHA_MLT[], mat.efrac)
TF.setscalar!(:dtime, parse(Float64, get(ENV, "DT0", "1e3"))); TF.setscalar!(:tmax, parse(Float64, get(ENV, "TEND", "1e12")))
haskey(ENV, "STOL") && TF.setscalar!(:stol, parse(Float64, ENV["STOL"]))
TC.STATIC_FALLBACK[] = true
dtmax = parse(Float64, get(ENV, "DTMAX", "1e30")); grow = parse(Float64, get(ENV, "GROW", "1.3")); smaxf = parse(Float64, get(ENV, "SMAXF", "0.3")); itmax = parse(Int, get(ENV, "ITMAX", "4"))
dn, tn, rn = TF.var(:dn), TF.var(:tn), TF.var(:rn)
function conv_core()          # mass fraction of the central convective region (contiguous from the centre), from the last mix diagnostics
    D = mat.Dmix; j = a + 1
    while j <= b && D[j] > 0; j += 1; end
    sum(4π * dn[k] * TF.var(:dvoln)[k] for k in a:j-1) / st.M
end
Xc() = mat.Y[SM.IH, a]
pr("%7s %10s %9s %9s %8s %9s %9s %9s %8s %9s %9s %7s %9s %9s %10s %7s %7s\n", "step", "t [yr]", "dt [yr]", "R [Rsun]", "logL", "logTeff", "Tc [1e8]", "rho_c", "Lnuc/L", "X_c", "Mcc/M", "iter", "cap_post", "cap_pred", "Zc(CNO F)", "eCNO/e", "e3a/e")
nstep = Ref(0); t0 = time()
function profile(io = nothing)
    rn, un, frn, ern, pgn, xkn = (TF.var(s) for s in (:rn, :un, :frn, :ern, :pgn, :chifn))
    sd = lastsd[]
    pr("%4s %10s %10s %10s %10s %9s %8s %9s %10s %10s %8s %8s %8s %9s\n", "zone", "m_ext/M", "r/Rsun", "rho", "T", "kappa_R", "beta", "u [km/s]", "Lrad/Lsun", "Lconv/Lsun", "grad", "grad_ad", "Gamma", "X_H")
    msum = cumsum([4π * dn[k] * TF.var(:dvoln)[k] for k in a:b]); Mt = msum[end]
    for k in unique(sort([a:6:b; b-20:b]))
        k < a && continue
        j = k - a + 1; Pr = ern[k] / 3; L = 4π * rn[k+1]^2 * frn[k+1]; mi = msum[j]
        Γ = xkn[k] * (L + 4π * rn[k+1]^2 * sd.fcv[k+1]) / (4π * Z0.c * Z0.G * mi)
        Pm = pgn[k] + Pr; ∇ = k < b ? log(tn[k+1] / tn[k]) / log((pgn[k+1] + ern[k+1] / 3) / Pm) : NaN
        TX.ZONE[] = k; ∇ad = TX.material_thermo(mat, dn[k], tn[k])[1]
        pr("%4d %10.3e %10.4f %10.3e %10.4e %9.4f %8.4f %9.3f %10.3e %10.3e %8.4f %8.4f %8.4f %9.5f\n", j, 1 - mi / Mt, rn[k+1] / Rsun, dn[k], tn[k], xkn[k], pgn[k] / (pgn[k] + Pr), un[k+1] / 1e5, L / Lsun, 4π * rn[k+1]^2 * sd.fcv[k+1] / Lsun, ∇, ∇ad, Γ, SM.hmass(mat, k))
    end
end
const lastsd = Ref{Any}(nothing)
pat = haskey(ENV, "PROFILE_AT") ? parse.(Float64, split(ENV["PROFILE_AT"], ",")) : Float64[]
nucstep = parse(Int, get(ENV, "NUCSTEP", "150"))      # the relaxation of the initial model's dynamical imbalance runs without burning
mat.nuc = false
function every(k)
    k == nucstep && (mat.nuc = nuc; pr("-- nuclear burning on at step %d (t = %.3e yr)\n", k, TF.getscalar(:timen) / 3.15576e7))
    sd = TX.step_diag(TF.getscalar(:timen) - TF.getscalar(:timeo))      # MLT/source diagnostics at the converged state, before burn!/mix! touch it
    lastsd[] = sd
    SM.burn!(mat; sd)
    mixon && SM.mix!(mat; sd)
    nstep[] = k
    if !isempty(pat) && TF.getscalar(:timen) / 3.15576e7 >= pat[1]
        pr("-- profile at t = %.4e yr (step %d)\n", TF.getscalar(:timen) / 3.15576e7, k); profile(); popfirst!(pat)
    end
    if k % parse(Int, get(ENV, "EVERY", "20")) == 0 || k <= 5
        d = star_diag(st)
        eg = let Yc = view(mat.Y, :, a); NN.energy(Yc, NN.rates_k(tn[a], dn[a], NN.Plasma(dn[a], Yc), NN.ye_of(Yc)))[3]; end
        pr("%7d %10.4e %9.3e %9.4f %8.4f %9.5f %9.5f %9.4f %8.4f %9.5f %9.5f %7d %9.2e %9.2e %10.3e %7.3f %7.4f\n", k, d.t / 3.15576e7, TF.getscalar(:dtime) / 3.15576e7, d.R / Rsun,
           log10(d.Ls / Lsun), log10(d.Teff), d.Tc / 1e8, d.ρc, d.Lnuc / d.Ls, Xc(), conv_core(), TF.getscalar(:iter), mat.dtcap / 3.15576e7, SM.predict_cap(mat; verbose = true) / 3.15576e7, SM.xyz(view(mat.Y, :, a))[3], eg.cno / max(eg.pp + eg.d + eg.cno + eg.he, 1e-300), eg.he / max(eg.pp + eg.d + eg.cno + eg.he, 1e-300))
        get(ENV, "CAPDIAG", "0") == "1" && pr("      predict_cap limiter: %s\n", string(SM.WHICH[]))
    end
    # quasi-static relaxation: the star is far from thermal equilibrium on the dynamical scale only briefly; allow dt to grow
    # beyond TITAN's stol rule while the step changed the state by less than SMAXF (implicit θ = 1 damps the dynamical modes)
    dtl = TF.getscalar(:timen) - TF.getscalar(:timeo)
    if grow > 1 && TF.getscalar(:smax) < smaxf && TF.getscalar(:iter) <= itmax && mat.dtcap > grow * dtl
        TF.setscalar!(:dtime, max(TF.getscalar(:dtime), grow * dtl))
    end
    TF.setscalar!(:dtime, min(TF.getscalar(:dtime), dtmax))
end
try
    n_, why = TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "400")), output = every)
    pr("done: %d steps (%s), %.1f s wall\n", n_, why, time() - t0)
catch e
    pr("STOPPED at step %d: %s\n", nstep[], sprint(showerror, e)[1:min(end, 300)])
end
d = star_diag(st)
pr("final: t = %.4e yr, R = %.4f Rsun, log L = %.4f, log Teff = %.5f, T_c = %.5e K, rho_c = %.4f, L_nuc/L = %.4f, X_c = %.5f, M_cc/M = %.5f\n",
   d.t / 3.15576e7, d.R / Rsun, log10(d.Ls / Lsun), log10(d.Teff), d.Tc, d.ρc, d.Lnuc / d.Ls, Xc(), conv_core())
get(ENV, "PROFILE", "0") == "1" && profile()
if get(ENV, "DTDIAG", "0") == "1"
    kx = TF.var(:kx); dx = TF.var(:dx)
    pr("timstp limiter (var: zone, relative change): %s\n", join(["$(i): k=$(kx[i]) dx=$(round(dx[i]; sigdigits=3))" for i in 1:Int(TF.getscalar(:neqn))], "; "))
    pr("index map ir im id iu it ie iF: %s, ngrs=%d\n", string([Int(TF.getscalar(s)) for s in (:ir, :im, :id, :iu, :it, :ie, :if)]), Int(TF.getscalar(:ngrs)))
end

# ---- composition of the centre and zone-by-zone energy / luminosity summary
let
    Yc = mat.Y[:, a]; Xs = SM.xyz(Yc)
    pr("centre: X_H = %.5f, Y_He = %.5f, Z(C,N,O,F) = %.3e;  X(12C) = %.3e, X(14N) = %.3e, X(16O) = %.3e, X(3He) = %.3e\n", Xs..., 12Yc[NN.IDX[:c12]], 14Yc[NN.IDX[:n14]], 16Yc[NN.IDX[:o16]], 3Yc[NN.IDX[:he3]])
    d = star_diag(st)
    pr("L_surface = %.4e erg/s (%.4e Lsun), L_nuc = ∫ε dm = %.4e (ratio L_nuc/L = %.4f), tracked ledger: ΣEinj %.4e ΣEbur %.4e, splitting correction %.3e, T re-solve error %.1e\n",
       d.Ls, d.Ls / Lsun, d.Lnuc, d.Lnuc / d.Ls, mat.ledger.Einj, mat.ledger.Ebur, mat.ledger.dabs / max(mat.ledger.Eabs, 1e-300), mat.ledger.resid)
    pr("convective core: M_cc/M = %.4f (contiguous from the centre, D > 0), Dmax = %.3e cm^2/s\n", conv_core(), maximum(mat.Dmix))
end
function mixtest()
    sd = lastsd[]; dvn = TF.var(:dvoln)
    jc = a + 1; while jc <= b && sd.fcv[jc] > 0 && sd.vcv[jc] > 0; jc += 1; end      # interfaces a+1 .. jc-1 convective; zones a .. jc-1
    ncz = jc - a
    if ncz >= 6
        tcross = sum((rn[k+1] - rn[k]) / sd.vcv[k+1] for k in a+1:jc-2)         # core crossing time ∫ dr / v
        vmean = (rn[jc] - rn[a+1]) / tcross
        pr("c2: convective core zones %d..%d (M_cc/M %.4f, r_cc = %.3f Rsun), mean v = %.3e cm/s, crossing time ∫dr/v = %.3e s (%.3e yr), l/v at the middle = %.3e s\n",
           a, jc - 1, conv_core(), rn[jc] / Rsun, vmean, tcross, tcross / 3.15576e7, sd.lcv[a+ncz÷2] / sd.vcv[a+ncz÷2])
        half = a:a+ncz÷2-1; step = 0.1
        for k in half; mat.Y[SM.IH, k] -= step; mat.Y[NN.IDX[:he4], k] += step / 4; end
        sm0 = SM.species_mass(mat); mcore = [4π * dn[k] * dvn[k] for k in a:jc-1]
        jo = jc:min(b, jc + 40); outside0 = copy(mat.Y[:, jo]); e0 = copy(TF.var(:egn)[a:b])      # jo: the radiative layers above the core (the thin convective skin zones near the surface also mix, trivially)
        spread() = (x = [SM.hmass(mat, k) for k in a:jc-1]; (maximum(x) - minimum(x)) / step)
        xm = sum(SM.hmass(mat, k) * mcore[k-a+1] for k in a:jc-1) / sum(mcore)
        pr("c2: initial step in X_H of %.2f over the inner half of the core; mass-mean X_H %.6f; spread/step %.4f\n", step, xm, spread())
        tdone = 0.0; nsub = 40; cons = 0.0
        for (nt, ttarget) in enumerate((1.0, 3.0, 10.0, 30.0))
            while tdone < ttarget * tcross * (1 - 1e-12)
                dts = tcross / nsub * (nt == 1 ? 1 : 4)
                r = SM.mix!(mat; dt = dts, sd = sd); cons = max(cons, r.cons); tdone += dts
            end
            sm1 = SM.species_mass(mat)
            pr("c2: after %5.1f crossing times: spread/step = %.3e, mass-mean X_H %.8f, max conservation error of species masses (whole star) %.2e\n",
               tdone / tcross, spread(), sum(SM.hmass(mat, k) * mcore[k-a+1] for k in a:jc-1) / sum(mcore), maximum(abs.(sm1 .- sm0) ./ max.(sm0, 1e-300)))
        end
        pr("c2: the 40 radiative zones above the core: max |ΔY| = %.1e (exact zero expected); max |Δe|/e of any zone = %.2e (e held fixed)\n",
           maximum(abs.(mat.Y[:, jo] .- outside0)), maximum(abs.(TF.var(:egn)[a:b] ./ e0 .- 1)))
        pr("c2: %s\n", spread() < 1e-3 && maximum(abs.(mat.Y[:, jo] .- outside0)) == 0.0 ? "PASS" : "FAIL")
    else
        pr("c2: no convective core of ≥ 6 zones, mixing test skipped\n")
    end
end
get(ENV, "MIXTEST", "1") == "1" && mixon && mixtest()
