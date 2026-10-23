# S8b acceptance tests: the 21 cm hyperfine diagnostic (HI21.jl).
#   julia --project=env s8b_test.jl [a|b|c|d|e|f|all]
# a  κ10 tables against the published values (FOB06 Tables 3–4 text values, FF07b Table 1), interpolation, ordering of the three species
# b  limits of T_s: no coupling → T_γ (also with a radio excess); strong collisions → T_k; strong Lyα → T_c; explicit rates vs the linearised formula; detailed balance; relaxation
# c  Wouthuysen–Field: critical flux J_αC (Pritchard & Loeb 2012 eq. 13 ff.), P_α ↔ Hirata's x̃_α, S̃_α and T_c^eff fits, fixed point
# d  IGM sanity check: mean-density gas at z = 69 (T_s, δT_b ≈ −40 mK: FOB06 Fig. 6, PL12 Fig. 3), z = 20 without / with Lyα coupling against PL12 eq. 7
# e  ray tracing: uniform static sphere vs the analytic chord integral; disc-mean of the thin sphere; Sobolev limit of an expanding sphere; mirror symmetry
# f  snapshots of the z = 19 collapse (if present): sane T_s, τ, equilibration times
using Printf, Serialization
const J = @__DIR__
include(joinpath(J, "HI21.jl")); include(joinpath(J, "Perf.jl"))
const H = HI21
const RESULTS = Tuple{String,Bool}[]
function check(name, ok, detail)
    push!(RESULTS, (name, ok)); @printf("%s  %-28s %s\n", ok ? "PASS" : "FAIL", name, detail); flush(stdout)
end
which = isempty(ARGS) ? "all" : ARGS[1]
Perf.reset!("s8b_test"; stage = "test")
relerr(a, b) = abs(a / b - 1)

# ====================================================================== a. κ10 tables
function test_a()
    # values quoted in the text of FOB06 / FF07 (independent of the tables): κHH(10 K) = 2.88e-12, κHH(1000 K) = 3.31e-10, κeH(100 K) = 2.26e-9, κeH(1e4) = 8.37e-9, κpH(10 K) = 0.3699e-9,
    # the H–H fall-off at low T (1.38e-13 at 1 K)
    pub = [("HH 10 K", H.κ_HH(10.0), 2.88e-12), ("HH 100 K", H.κ_HH(100.0), 1.19e-10), ("HH 1000 K", H.κ_HH(1000.0), 3.31e-10), ("HH 1 K", H.κ_HH(1.0), 1.38e-13),
           ("eH 100 K", H.κ_eH(100.0), 2.26e-9), ("eH 1e4 K", H.κ_eH(1e4), 8.37e-9), ("eH 1 K", H.κ_eH(1.0), 2.39e-10),
           ("pH 10 K", H.κ_pH(10.0), 0.3699e-9), ("pH 1e4 K", H.κ_pH(1e4), 1.695e-9), ("HH 70 K", H.κ_HH(70.0), 9.25e-11)]
    worst = maximum(relerr(a, b) for (_, a, b) in pub)
    check("a1 published values", worst < 0.02, @sprintf("worst relative deviation %.1e over %d values (κ_HH(88 K) = %.3e, used in FOB06 eq. 67)", worst, length(pub), H.κ_HH(88.0)))
    # table integrity: all rows positive, the two electron/proton rows share their temperatures; monotonic H–H and H–e; interpolation inside the bracketing rows
    mono = all(diff(H.K_HH) .> 0) && all(diff(H.K_EH[1:15]) .> 0)
    ins = all(H.K_HH[j] < H.κ_HH(sqrt(H.T_HH[j] * H.T_HH[j+1])) < H.K_HH[j+1] for j in 1:length(H.T_HH)-1)
    check("a2 table interpolation", mono && ins && H.κ_HH(3e4) == H.K_HH[end] && H.κ_eH(0.1) == H.K_EH[1], "monotonic; log-log interpolation inside the rows; clamped outside (HH above 1e4 K, eH/pH above 2e4 K, all below 1 K)")
    rat = [H.κ_eH(T) / H.κ_HH(T) for T in (10.0, 100.0, 1000.0, 1e4)]; pr = H.κ_pH(1e4) / H.κ_HH(1e4)
    # FF07a: κ_eH/κ_HH ≫ 1 at low T and ~ sqrt(m_p/m_e)/few at 1e4 K; FF07b: p–H marginally above H–H (factor ~2) at high T, vastly above it at low T
    check("a3 species ordering", rat[1] > 200 && rat[4] > 5 && 1.5 < pr < 3.5 && H.κ_pH(10.0) / H.κ_HH(10.0) > 100, @sprintf("κeH/κHH = %.0f, %.0f, %.0f, %.0f at 10, 100, 1e3, 1e4 K; κpH/κHH = %.2f at 1e4 K, %.0f at 10 K", rat..., pr, H.κ_pH(10.0) / H.κ_HH(10.0)))
end

# ====================================================================== b. limits
function test_b()
    Tg = 54.5; Tk = 300.0
    Ts0, = H.spin_temperature(Tg, Tk, 0.0)
    Tgr = H.Tgamma(2.7255 * 20, 19.0; Ar = 0.1)
    Tsr, = H.spin_temperature(Tgr, Tk, 0.0)
    check("b1 no coupling", relerr(Ts0, Tg) < 1e-14 && relerr(Tsr, Tgr) < 1e-14 && Tgr > Tg, @sprintf("T_s = T_γ = %.4f K; with a radio excess Ar = 0.1 K at 78 MHz: T_γ = %.3f K and T_s = %.3f K", Ts0, Tgr, Tsr))
    xc, = H.xcoll(Tk, Tg; nH = 1e8)
    Tsc, = H.spin_temperature(Tg, Tk, xc)
    check("b2 collisions dominate", relerr(Tsc, Tk) < 1e-3 && xc > 1e3, @sprintf("n_H = 1e8 cm⁻³: x_c = %.2e, T_s = %.5f K (T_k = %.1f)", xc, Tsc, Tk))
    Tc = 20.0; P = 1e-4
    Tsa, xa, Tcv, S = H.spin_temperature(Tg, Tk, 0.0; P = P, Tc = Tc, kind = :actual)
    check("b3 Lyα dominates", relerr(Tsa, Tc) < 1e-3 && S == 1.0, @sprintf("P_α = 1e-4 s⁻¹ (actual line intensity), T_c = 20 K, T_k = 300 K: x_α = %.2e, T_s = %.5f K", xa, Tsa))
    Tsq, xq, Tcq, Sq = H.spin_temperature(Tg, Tk, 0.0; P = P, kind = :continuum, τGP = 1e6)
    check("b4 Lyα continuum", relerr(Tsq, Tk) < 1e-2 && abs(Tcq - Tk) < 0.01 * Tk && Sq < 1, @sprintf("continuum Lyα: T_s = %.3f K (T_k), T_c^eff = %.2f K, S̃_α = %.4f", Tsq, Tcq, Sq))
    # linearised vs explicit rates (the linearisation is O(T*/T_γ) ~ 1e-3); detailed balance exactly
    d = 0.0
    for (tg, tk, xc, xa) in ((54.5, 7.0, 1e-2, 0.0), (54.5, 7.0, 0.5, 2.0), (190.0, 88.0, 0.9, 0.0), (54.5, 3000.0, 10.0, 1e-2))
        a = H.ts_lin(tg, tk, xc, xa, tk); b = H.ts_rates(tg, tk, xc, xa, tk); d = max(d, relerr(a, b))
    end
    e1 = relerr(H.ts_rates(54.5, 54.5, 0.0, 0.0, 54.5), 54.5); e2 = relerr(H.ts_rates(54.5, 123.0, 1e14, 0.0, 123.0), 123.0)
    check("b5 explicit rates", d < 3e-3 && e1 < 1e-9 && e2 < 1e-6, @sprintf("linearised vs explicit: worst %.1e (T*/T_γ = %.1e); detailed balance: CMB only %.1e, collisions only %.1e", d, H.TSTAR / 54.5, e1, e2))
    # relaxation: t → ∞ equilibrium, inverse maps, CMB-only t_eq = 1/(A10 (1 + 4 n̄))
    up, dn = H.rates(54.5, 7.0, 0.5, 0.0, 7.0)
    y, Tn, teq = H.relax(H.y_from_ts(54.5), 1e30, 54.5, 7.0, 0.5, 0.0, 7.0)
    y0, = H.relax(0.123, 0.0, 54.5, 7.0, 0.5, 0.0, 7.0); yh, Th, _ = H.relax(H.y_from_ts(54.5), 1 / (up + dn), 54.5, 7.0, 0.5, 0.0, 7.0)
    nbar = 1 / expm1(H.TSTAR / 54.5); _, _, teq0 = H.relax(0.1, 1.0, 54.5, 54.5, 0.0, 0.0, 54.5)
    ye = up / (up + dn)
    check("b6 relaxation", relerr(Tn, H.ts_lin(54.5, 7.0, 0.5, 0.0, 7.0)) < 5e-3 && y0 == 0.123 && abs((yh - ye) / (H.y_from_ts(54.5) - ye) - exp(-1)) < 1e-12 && relerr(H.ts_from_y(H.y_from_ts(33.0)), 33.0) < 1e-12 && relerr(teq0, 1 / (H.A10 * (1 + 4nbar))) < 1e-12,
          @sprintf("T_s(∞) = %.4f K (linear %.4f); one e-fold after t_eq = %.3e s; CMB-only t_eq = %.3e s = %.0f yr (z = 19)", Tn, H.ts_lin(54.5, 7.0, 0.5, 0.0, 7.0), teq, teq0, teq0 / 3.15576e7))
end

# ====================================================================== c. Wouthuysen–Field
function test_c()
    Tg = 2.7255 * 20; JC = 27 * H.A10 * Tg / (16π * H.SIGINT_A * H.TSTAR)          # x_α = 1 for S = 1; photons cm⁻² s⁻¹ Hz⁻¹ sr⁻¹
    # PL12: J_αC = 1.165e-10 [(1+z)/20] (their eq. 13 with rounded constants; direct evaluation of eq. 13 with e, m_e, c, f_α, T* = 0.0682, A10 = 2.85e-15 gives 1.107e-10)
    check("c1 critical flux", relerr(JC, 1.165e-10) < 0.07, @sprintf("J_αC(z = 19) = %.4e (PL12: 1.165e-10; ∫σ dν = %.4e cm² Hz = π e² f/m_e c with f = 0.4162)", JC, H.SIGINT_A))
    # x̃_α = 8π λ² γ T*/(9 A10 T_γ) S̃ J_α (Hirata eq. 38), J_α in photons: same number
    Jn = 3e-11; P = H.P_from_J(Jn * H.hP * H.NU_A)
    xh = 8π * H.LAM_A^2 * H.GAM_A * H.TSTAR / (9 * H.A10 * Tg) * Jn
    check("c2 P_α ↔ Hirata x̃_α", relerr(H.xalpha(P, Tg), xh) < 2e-3, @sprintf("x_α = %.5e from P_α = %.3e s⁻¹; Hirata eq. 38: %.5e (λ, γ = A_α/4π = %.3e Hz)", H.xalpha(P, Tg), P, xh, H.GAM_A))
    # S̃_α against the older fit exp(-0.803 T^-2/3 (1e-6 τ)^1/3) (Chen & Miralda-Escudé 2004; PL12) and Hirata's statement that the effects are largest at low T_k
    s100 = H.S_tilde(100.0, 100.0, 2e6); old100 = exp(-0.803 * 100^(-2 / 3) * 2^(1 / 3)); s5 = H.S_tilde(5.0, 5.0, 2e6); old5 = exp(-0.803 * 5^(-2 / 3) * 2^(1 / 3))
    check("c3 S̃_α fit", abs(s100 / old100 - 1) < 0.06 && s5 < s100 < 1 && H.S_tilde(1e4, 1e4, 2e6) > 0.99, @sprintf("S̃_α(T_k = T_s = 100 K, τ_GP = 2e6) = %.4f (CM04 fit %.4f); at 5 K %.4f (%.4f); at 1e4 K %.4f", s100, old100, s5, old5, H.S_tilde(1e4, 1e4, 2e6)))
    # T_c^eff between T_s and T_k, equal to T_k when T_s = T_k; implicit solution is a fixed point
    tc = H.Tc_eff(100.0, 20.0); Ts, xa, Tcv, S = H.spin_temperature(54.5, 20.0, 0.1; P = 2e-12, kind = :continuum, τGP = 2e6)
    res = abs(Ts - H.ts_lin(54.5, 20.0, 0.1, H.xalpha(2e-12, 54.5, H.S_tilde(20.0, Ts, 2e6)), H.Tc_eff(20.0, Ts)))
    check("c4 T_c, fixed point", H.Tc_eff(100.0, 100.0) == 100.0 && 20.0 < tc < 100.0 && res < 1e-8 * Ts, @sprintf("T_c^eff(T_k = 100, T_s = 20) = %.2f K; implicit T_s = %.5f K (x_α = %.3e, S̃ = %.4f), residual %.1e", tc, Ts, xa, S, res))
end

# ====================================================================== d. IGM
function test_d()
    Ωm = 0.31; h = 0.68; Ωbh2 = 0.0223
    nH0 = (1 - 0.245) * Ωbh2 * 1.8788e-29 / H.mH                      # n_H today [cm⁻³] (Y_p = 0.245)
    igm(z; xe = 2e-4, Tk = 2.7255 * (1 + z)^2 / 151, Lyα = false) = begin
        Tg = 2.7255 * (1 + z); nH = nH0 * (1 + z)^3
        xc, = H.xcoll(Tk, Tg; nH = nH * (1 - xe), ne = nH * xe, np = nH * xe)
        Ts = Lyα ? Tk : H.ts_lin(Tg, Tk, xc, 0.0, Tk)
        τ = H.tau_sobolev(nH * (1 - xe), Ts, H.Hz(z; Ωm = Ωm, h = h))
        (; z, Tg, Tk, nH, xc, Ts, τ, dTb = (Ts - Tg) * (-expm1(-τ)) / (1 + z) * 1e3)
    end
    a = igm(69.0)
    check("d1 dark-age absorption", -45 < a.dTb < -30 && 0.5 < a.xc < 1.5, @sprintf("z = 69 (1+z = 70), T_k = %.1f K (FOB06: 88 K), n_H = %.3e: x_c = %.3f, T_s = %.1f K, T_γ = %.1f K, τ = %.4f, δT_b = %.1f mK (published peak ≈ −40 mK, FOB06 Fig. 6b / PL12 Fig. 3)", a.Tk, a.nH, a.xc, a.Ts, a.Tg, a.τ, a.dTb))
    b = igm(19.0); c = igm(19.0; Lyα = true)
    # PL12 eq. 7 (neglecting the velocity-gradient term): δT_b = 27 x_HI (1+δ) (Ωbh²/0.023) sqrt(0.15 (1+z)/(Ωm h² 10)) (T_s − T_R)/T_s mK
    dlin = (c.Ts - c.Tg) / (1 + 19.0) * c.τ * 1e3                          # optically thin form (PL12 eq. 6), the form of eq. 7
    pl = 27 * (1 - 2e-4) * (Ωbh2 / 0.023) * sqrt(0.15 / (Ωm * h^2) * (1 + 19.0) / 10) * (c.Ts - c.Tg) / c.Ts
    check("d2 z = 19 IGM", abs(b.dTb) < 1 && relerr(dlin, pl) < 0.03, @sprintf("z = 19, T_k = %.2f K, x_c = %.2e: without Lyα T_s = %.2f K, δT_b = %.2f mK (collisions decoupled); fully Lyα-coupled T_s = T_k: δT_b = %.1f mK (τ = %.3f; thin form %.1f mK) against PL12 eq. 7 (27 mK coefficient): %.1f mK", b.Tk, b.xc, b.Ts, b.dTb, c.dTb, c.τ, dlin, pl))
end

# ====================================================================== e. rays
function uniform_sphere(n; R = 1e17, T = 100.0, nHI = 1e-2, Ts = 20.0, u = 0.0)
    r = collect(range(0, R, length = n + 1)); (; r, u = u .* r ./ R, Tk = fill(T, n), Ts = fill(Ts, n), nHI = fill(nHI, n), R)
end
function test_e()
    sph = uniform_sphere(100; nHI = 100.0); Tg = 54.5
    b = H.bth(100.0); vg = collect(range(-6b, 6b, length = 241))
    k0 = H.κ0(sph.nHI[1], sph.Ts[1])
    worst = 0.0; tmax = 0.0
    for p in (0.0, 0.3sph.R, 0.6sph.R, 0.9sph.R)
        Tb = H.spectrum(sph.r, sph.u, sph.Tk, sph.Ts, sph.nHI, p, vg; Tbg = Tg)
        τ = k0 .* exp.(-(vg ./ b) .^ 2) ./ (sqrt(π) * b) .* 2 * sqrt(sph.R^2 - p^2)
        ana = Tg .* exp.(-τ) .+ sph.Ts[1] .* (1 .- exp.(-τ))
        dd = abs.(ana .- Tg); m = dd .> 1e-4 * maximum(dd); worst = max(worst, maximum(abs.(Tb[m] .- ana[m]) ./ dd[m])); tmax = max(tmax, maximum(τ))
    end
    check("e1 static sphere", worst < 1e-3, @sprintf("T_b(v) of the rays p/R = 0, 0.3, 0.6, 0.9 vs T_γ e^-τ + T_s (1 - e^-τ), τ = κ0 φ_v 2√(R²-p²) (max τ %.2f): worst relative deviation of δT %.1e", tmax, worst))
    # thin sphere: disc-mean of δT = (T_s − T_γ) <τ>, <τ> = κ0 φ_v (4R/3) (volume / disc area)
    th = uniform_sphere(300; nHI = 1e-8)
    o = (; Tk = th.Tk, Ts_use = th.Ts, nHI = th.nHI, Tg = fill(Tg, 300))
    s = (; r = th.r, u = th.u)
    sp = H.cloud_spectra(s, o; z = 19.0)
    k0t = H.κ0(1e-8, 20.0); ana = (20.0 - Tg) .* k0t .* exp.(-(sp.v ./ b) .^ 2) ./ (sqrt(π) * b) .* 4th.R / 3
    i0 = argmin(abs.(sp.v)); e = relerr(sp.mean[i0], ana[i0]); ew = relerr(H.eqwidth(sp.v, sp.mean), H.eqwidth(sp.v, ana))
    check("e2 thin sphere, disc mean", e < 1e-2 && ew < 5e-3, @sprintf("area-weighted δT(v=0) = %.4e K vs analytic %.4e (rel. %.1e); equivalent width rel. %.1e (%d rays, %d velocity bins)", sp.mean[i0], ana[i0], e, ew, length(sp.p), length(sp.v)))
    # Sobolev: uniform expanding sphere u = H r, R ≫ b/H; central ray at v = 0 has τ = κ0/H; check optically thin and thick
    Hh = 1e-14; R = 3e21; b100 = H.bth(100.0)
    res = String[]; ok = true
    for nHI in (0.3, 12.0)
        ex = uniform_sphere(2000; R = R, nHI = nHI, u = Hh * R); vg2 = collect(range(-4b100, 4b100, length = 161))
        Tb = H.spectrum(ex.r, ex.u, ex.Tk, ex.Ts, ex.nHI, 0.0, vg2; Tbg = Tg)
        τS = H.tau_sobolev(nHI, 20.0, Hh); ana = Tg * exp(-τS) + 20.0 * (1 - exp(-τS)); v0 = Tb[81]
        ok &= abs((v0 - Tg) / (ana - Tg) - 1) < 0.03; push!(res, @sprintf("τ_S = %.3f: T_b(v=0) = %.4f vs %.4f", τS, v0, ana))
    end
    check("e3 Sobolev limit", ok, join(res, "; ") * @sprintf(" (b/H = %.1e cm = %.1e R)", b100 / Hh, b100 / Hh / R))
    # mirror symmetry: infall ↔ outflow with v → −v
    inf = uniform_sphere(60; nHI = 1.0, u = -3e5); out = uniform_sphere(60; nHI = 1.0, u = 3e5); vg3 = collect(range(-8e5, 8e5, length = 321))
    d = maximum(abs.(H.spectrum(inf.r, inf.u, inf.Tk, inf.Ts, inf.nHI, 0.4e17, vg3; Tbg = Tg) .- reverse(H.spectrum(out.r, out.u, out.Tk, out.Ts, out.nHI, 0.4e17, vg3; Tbg = Tg))))
    check("e4 mirror symmetry", d < 1e-9, @sprintf("infall vs outflow (3 km/s at the edge) with v → −v: max |ΔT_b| = %.1e K", d))
end

# ====================================================================== f. snapshot
function test_f()
    using_file = joinpath(J, "..", "lib", "popiii", "snaps_abn02z19_lines_n200_final2.jls")
    alt = "/Users/tabel/Projects/titan.jl/lib/popiii/snaps_abn02z19_lines_n200_final2.jls"
    f = isfile(using_file) ? using_file : alt
    if !isfile(f)
        check("f1 snapshot (skipped)", true, "no snapshot file in lib/popiii (not committed): skipped"); return
    end
    snaps = Serialization.deserialize(f); s = snaps[end]
    o = H.zones(s; Tcmb = 54.5, z = 19.0)
    lo = min.(o.Tk, o.Tg); hi = max.(o.Tk, o.Tg)
    ok = all(isfinite, o.Ts) && all(o.Ts .>= lo .* (1 - 1e-9)) && all(o.Ts .<= hi .* (1 + 1e-9)) && all(o.tau .>= 0) && all(isfinite, o.teq)
    j = argmin(abs.(o.nHI .- 1e6)); k = argmax(o.nHI)
    check("f1 snapshot z = 19", ok, @sprintf("final2 snapshot: T_s between min/max(T_k, T_γ) in all %d zones; n_HI max %.2e → T_s = %.1f K (T_k = %.1f), outer zone T_s = %.2f K (T_k = %.1f, T_γ = 54.5), max τ21 = %.2e, t_eq/t_ff: %.1e (outer) … %.1e (centre)",
          length(o.Ts), o.nHI[k], o.Ts[k], o.Tk[k], o.Ts[end], o.Tk[end], maximum(o.tau), o.teq[end] / o.tff[end], o.teq[1] / o.tff[1]))
end

which in ("a", "all") && Perf.@section :a test_a()
which in ("b", "all") && Perf.@section :b test_b()
which in ("c", "all") && Perf.@section :c test_c()
which in ("d", "all") && Perf.@section :d test_d()
which in ("e", "all") && Perf.@section :e test_e()
which in ("f", "all") && Perf.@section :f test_f()
Perf.write()
println("S8B: ", all(last, RESULTS) ? "ALL PASS" : "$(count(!last, RESULTS)) FAIL")
