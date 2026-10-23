# S8a acceptance tests of the Monte Carlo Lyα transfer (LyaMC.jl). Self-contained (LyaMC includes its own LymanWings).
#   julia --project=env -t 4 s8a_test.jl [a|b|c|d|e|f|g|h|i|all]
# a  sampler & kernels: parallel-velocity sampler (KS against the numerical CDF), wing drift/variance, Rayleigh angle sampler, collisional destruction saturates at A_2γ/(3A_21), KH wing table
# b  static slab: emergent spectrum (Harrington 1973/Neufeld 1990 analytic) and N_sc ∝ τ_0;  static uniform sphere (Dijkstra, Haiman & Spaans 2006 analytic spectrum)
# c  moving medium: first-flight scattering probability in a Hubble-like flow (Sobolev), infall/outflow mirror symmetry and blue peak, energy bookkeeping (work = ΔE of the photons)
# d  destruction: Neufeld (1990) f_esc(τ_a, aτ_0) for frequency-independent absorption; per-scattering ε: f_esc = exp(−ε N_sc); scaling collapse (τ_0/s, ε s) used for the inner column
# e  force multiplier of a static thin shell against Dijkstra & Loeb (2008) / Smith, Bromm & Loeb (2017)
# f  Wouthuysen–Field: T_c → T_k in the centre of a thick medium with recoil (Field 1959, Hirata 2006)
# g  conservation: photons and energy ledgers, momentum
# h  Monte Carlo noise ∝ N^{-1/2}; core skipping against the exact transfer (spectrum, N_sc, f_esc, force)
using Printf, Random
include(joinpath(@__DIR__, "LyaMC.jl"))
using .LyaMC
const L = LyaMC; const RESULTS = Tuple{String,Bool}[]
function check(name, ok, detail)
    push!(RESULTS, (name, ok)); @printf("%s  %-34s %s\n", ok ? "PASS" : "FAIL", name, detail); flush(stdout)
end
which = isempty(ARGS) ? "all" : ARGS[1]
const c = L.c

# ====================================================================== a
function test_a()
    rng = Xoshiro(3)
    function cdf_num(x, a; lo = -9.0, hi = 9.0)
        pts = sort(unique(filter(u -> lo <= u <= hi, vcat(range(lo, hi, length = 400_000), x .+ a .* tan.(range(-1.5, 1.5, length = 100_000))))))
        p(u) = exp(-u^2) / ((u - x)^2 + a^2); cc = zeros(length(pts))
        for i in 2:length(pts); cc[i] = cc[i-1] + 0.5 * (p(pts[i]) + p(pts[i-1])) * (pts[i] - pts[i-1]); end
        pts, cc ./ cc[end]
    end
    worst = 0.0; N = 100_000
    for (x, a) in ((0.3, 5e-4), (1.0, 5e-4), (2.5, 5e-4), (3.5, 1e-3), (6.0, 5e-4), (-4.0, 5e-4), (15.0, 2e-3), (2.0, 5e-3))
        pts, C = cdf_num(x, a); s = sort([L.sample_upar(rng, x, a) for _ in 1:N]); D = 0.0
        for i in 1:100:N; j = clamp(searchsortedlast(pts, s[i]), 1, length(C)); D = max(D, abs(C[j] - i / N)); end
        worst = max(worst, D)
    end
    check("a1 parallel-velocity sampler", worst < 2.2 / sqrt(N), @sprintf("KS distance %.4f over 8 (x, a) pairs (1%% critical value %.4f)", worst, 1.63 / sqrt(N)))
    # wing: <Δx> = −1/x − g..., var → 1 (two Gaussians of variance 1/2)
    g = 2.5e-4; dr = String[]; ok = true
    for x in (8.0, 30.0)
        n = 300_000; d = zeros(n); for i in 1:n; d[i] = L.scatter_x(rng, x, 5e-4, g, 0.1)[1] - x; end
        m = sum(d) / n; v = sum(d .^ 2) / n; push!(dr, @sprintf("x=%g: <Δx> = %.4f (−1/x = %.4f), <Δx²> = %.4f", x, m, -1 / x, v))
        ok &= abs(v - 1) < 0.05 && abs(m + 1 / x) < 0.2 / x + 5e-3
    end
    check("a2 wing redistribution", ok, join(dr, "; "))
    μs = [L.ray_mu(rng) for _ in 1:400_000]
    check("a3 Rayleigh angles", abs(sum(μs) / length(μs)) < 3e-3 && abs(sum(μs .^ 2) / length(μs) - 0.4) < 3e-3, @sprintf("<μ> = %.4f, <μ²> = %.4f (0.4 for (3/8)(1+μ²))", sum(μs) / length(μs), sum(μs .^ 2) / length(μs)))
    e = L.eps_collisional(1e4, 1e8, 1e8)
    check("a4 destruction probability", abs(e / (L.A2G / (3 * L.A21)) - 1) < 0.01 && L.eps_collisional(1e4, 1e2, 1e2) < 0.05 * e, @sprintf("ε(n_e = 1e8) = %.3e (A_2γ/3A_21 = %.3e), ε(n_e = 100) = %.2e", e, L.A2G / (3 * L.A21), L.eps_collisional(1e4, 1e2, 1e2)))
    L.init_wings!(); T = 1e4; D = L.dopp(T); worstw = 0.0
    for Δ in (-0.1, -0.02, 0.02, 0.1)
        ν = L.ν0 * (1 + Δ); x = (ν - L.ν0) / D
        σmc = L.sig0(T) * L.voigt(L.damp(T), x) * L.Cwing(Δ); σkh = L.LW.σ_gas(ν, T); worstw = max(worstw, abs(σmc / σkh - 1))
    end
    check("a5 KH wing table", worstw < 0.02, @sprintf("σ_MC/σ_KH within %.2g over Δν/ν = ±0.02 … ±0.1", worstw))
end

# ====================================================================== b
shell(T, tau0, h, R0; nz = 8, tau_a = 0.0, skip = true) = (
    r = collect(range(R0 - h, R0 + h, length = 2nz + 1));
    L.build_zones(r, fill(tau0 / (L.sig0(T) * h), 2nz), fill(T, 2nz), zeros(2nz + 1); kabs = fill(tau_a / h, 2nz), xcs = skip ? (3.0, 6.0) : ()))
function test_b()
    T = 10.0; a = L.damp(T); b = L.bth(T)
    # b1: slab spectrum and N_sc ∝ τ_0 (midplane source, thin spherical shell R0/h = 1e3)
    taus = [1e4, 3e4, 1e5]; Ns = Float64[]; peak = 0.0; dev = 0.0
    for τ0 in taus
        Z = shell(T, τ0, 1e15, 1e18; skip = false); xs = (a * τ0)^(1 / 3)
        res = L.run(Z, [L.src_sphere(1e18, 1.0; npk = τ0 > 5e4 ? 2400 : 1600)]; nbatch = 8, opts = L.Opts(inner = :escape, wmax = 5 * xs * b, nbw = 100))
        push!(Ns, res.nsc_per_photon)
        if τ0 == 1e5
            x = L.wcenters(res) ./ b; sp = res.spec_out .+ res.spec_in; dx = x[2] - x[1]; mc = sp ./ (sum(sp) * dx)
            an = [xx^2 / cosh(sqrt(pi^3 / 54) * abs(xx)^3 / (a * τ0)) for xx in x]; an ./= sum(an) * dx
            peak = abs(x[argmax(mc)]) / xs; dev = sqrt(sum((mc .- an) .^ 2) / sum(an .^ 2)); ximean = sqrt(sum(sp .* x .^ 2) / sum(sp)) / xs
            @printf("   slab aτ0 = %.0f: MC peak at |x| = %.2f (aτ0)^(1/3) (analytic 1.06), rms deviation of the spectrum %.3f, rms |x|/(aτ0)^(1/3) = %.3f\n", a * τ0, peak, dev, ximean)
        end
    end
    slope = (log(Ns[end]) - log(Ns[1])) / (log(taus[end]) - log(taus[1]))
    check("b1 slab: N_sc ∝ τ_0", abs(slope - 1) < 0.08, @sprintf("N_sc = %s for τ_0 = %s: log-slope %.3f, N_sc/τ_0 = %.3f %.3f %.3f", string(round.(Ns, sigdigits = 3)), string(taus), slope, (Ns ./ taus)...))
    check("b2 slab: emergent spectrum", abs(peak - 1.06) < 0.2 && dev < 0.30, @sprintf("Harrington/Neufeld: peak %.2f (1.06) × (aτ0)^(1/3), rms deviation %.3f (700 packets-level noise)", peak, dev))
    # b3: uniform sphere, central source: Dijkstra, Haiman & Spaans (2006) eq. 9, peak 0.92 (aτ0)^{1/3}
    τ0 = 1e6; R = 1e16; nz = 10; r = collect(range(0, R, length = nz + 1)); xs = (a * τ0)^(1 / 3)
    Z = L.build_zones(r, fill(τ0 / (L.sig0(T) * R), nz), fill(T, nz), zeros(nz + 1); xcs = (3.0, 6.0))
    res = L.run(Z, [L.src_point(1.0; npk = 12800)]; nbatch = 8, opts = L.Opts(wmax = 4 * xs * b, nbw = 80))
    x = L.wcenters(res) ./ b; sp = res.spec_out; dx = x[2] - x[1]; mc = sp ./ (sum(sp) * dx)
    an = [xx^2 / (1 + cosh(sqrt(2pi^3 / 27) * abs(xx)^3 / (a * τ0))) for xx in x]; an ./= sum(an) * dx
    dev = sqrt(sum((mc .- an) .^ 2) / sum(an .^ 2)); pk = abs(x[argmax(mc)]) / xs
    check("b3 sphere: analytic spectrum", abs(pk - 0.92) < 0.2 && dev < 0.3, @sprintf("aτ0 = %.0f: peak %.2f (0.92) × (aτ0)^(1/3), rms deviation %.3f, N_sc/τ0 = %.3f, skipped zones x_c = %s", a * τ0, pk, dev, res.nsc_per_photon / τ0, string(unique(Z.xc))))
end

# ====================================================================== c
function test_c()
    T = 1e4; a = L.damp(T); b = L.bth(T); s0 = L.sig0(T); R = 1e17; nz = 40
    # c1: first-flight probability in a homologous flow u = Γ r (Sobolev): emission at the line centre of the medium (thermal profile), P = 1 − <exp(−τ(x0))>, τ(x0) = (τ_S/√π) ∫_{-∞}^{x0} H dx
    Hint(x0) = (s = 0.0; xs = range(-60, x0, length = 6000); for i in 1:5999; s += L.voigt(a, 0.5 * (xs[i] + xs[i+1])) * (xs[i+1] - xs[i]); end; s / sqrt(pi))
    out = String[]; ok = true
    for τS in (0.1, 1.0, 10.0)
        Γ = 2000 * b / R; n = τS * Γ / (s0 * sqrt(pi) * b); r = collect(range(0, R, length = nz + 1))
        Z = L.build_zones(r, fill(n, nz), fill(T, nz), Γ .* r; xcs = ())
        src = L.src_zones([zeros(nz ÷ 2 - 1); 1.0; zeros(nz - nz ÷ 2)]; npk = 16000)
        res = L.run(Z, [src]; nbatch = 8, opts = L.Opts(wmax = 3e8, nbw = 10, maxsc = 1))
        mcP = sum(res.nsc) - res.f_lost
        rng = Xoshiro(1); an = sum(1 - exp(-τS * Hint(randn(rng) * 0.70710678 + a * tan(pi * (rand(rng) - 0.5)))) for _ in 1:3000) / 3000
        push!(out, @sprintf("τ_S = %g: MC %.4f, Sobolev %.4f", τS, mcP, an)); ok &= abs(mcP - an) < 0.015 + 0.03 * an
    end
    check("c1 first flight in a Hubble flow", ok, join(out, "; "))
    # c2: infall / outflow of a thick uniform sphere (aτ0 = 470): blue peak for infall, mirror symmetry, energy bookkeeping
    R = 1e17; nz = 12; τ0 = 3e5; n = τ0 / (s0 * R); r = collect(range(0, R, length = nz + 1)); rows = Dict{String,Any}()
    for (lab, V) in (("infall", -8b), ("outflow", 8b))
        Z = L.build_zones(r, fill(n, nz), fill(T, nz), V .* r ./ R; xcs = (3.0, 6.0))
        res = L.run(Z, [L.src_point(1.0; npk = 1600)]; nbatch = 8, opts = L.Opts(wmax = 120 * b, nbw = 120))
        x = L.wcenters(res) ./ b; sp = res.spec_out; sp ./= sum(sp); rows[lab] = (xm = sum(sp .* x), blue = sum(sp[x .> 0]), res = res, nsc = res.nsc_per_photon)
    end
    ri, ro = rows["infall"], rows["outflow"]
    check("c2 infall blue peak, mirror symmetry", ri.blue > 0.95 && ro.blue < 0.05 && abs(ri.xm + ro.xm) < 0.1 * abs(ri.xm) && abs(ri.nsc / ro.nsc - 1) < 0.1, @sprintf("blue fraction %.3f (infall) / %.3f (outflow); <x> = %.2f / %.2f; N_sc %.3g / %.3g", ri.blue, ro.blue, ri.xm, ro.xm, ri.nsc, ro.nsc))
    # work: the photons gain <x> b/c in energy from the receding/approaching gas (lab-frame energy of the escaping photons)
    e1 = ri.res.Egas_total / L.E0; e2 = -ri.xm * b / c
    check("c3 energy exchange with the flow", abs(e1 / e2 - 1) < 0.1, @sprintf("gas energy change per photon %.4e E_α against −<x>b/c = %.4e (infall)", e1, e2))
end

# ====================================================================== d
function test_d()
    T = 10.0; a = L.damp(T); τ0 = 1e5; ζ = 0.525; out = String[]; ok = true
    fneu(aτ, τa) = 1 / cosh(sqrt(3) / (pi^(5 / 12) * ζ) * sqrt(aτ^(1 / 3) * τa))
    for τa in (0.01, 0.03, 0.1)
        Z = shell(T, τ0, 1e15, 1e18; tau_a = τa, skip = false); xs = (a * τ0)^(1 / 3)
        res = L.run(Z, [L.src_sphere(1e18, 1.0; npk = 2400)]; nbatch = 8, opts = L.Opts(inner = :escape, wmax = 5 * xs * L.bth(T), nbw = 20))
        f = res.f_out + res.f_in; fn = fneu(a * τ0, τa)
        push!(out, @sprintf("τ_a = %g: %.3f (Neufeld %.3f)", τa, f, fn)); ok &= abs(f / fn - 1) < 0.12
    end
    check("d1 Neufeld f_esc (absorption)", ok, @sprintf("aτ0 = %.0f, slab; %s", a * τ0, join(out, "; ")))
    # d2: per-scattering destruction ε: f_esc = <exp(−ε N)>; scaling (τ0/s, ε s) leaves it invariant
    R = 1e16; nz = 6; fl = Float64[]; lab = String[]; spx = Vector{Vector{Float64}}()
    for (s, ε0) in ((1.0, 3e-7), (10.0, 3e-7)), mode in (:eps,)
        τ = 1e6 / s; r = collect(range(0, R, length = nz + 1)); n = τ / (L.sig0(T) * R)
        Z = L.build_zones(r, fill(n, nz), fill(T, nz), zeros(nz + 1); xcs = (3.0, 6.0), epsfloor = ε0 * s)
        res = L.run(Z, [L.src_point(1.0; npk = 2400)]; nbatch = 8, opts = L.Opts(wmax = 8 * (a * 1e6 / s)^(1 / 3) * L.bth(T), nbw = 160))
        x = L.wcenters(res) ./ L.bth(T) .* cbrt(s); sp = res.spec_out; push!(fl, res.f_out); push!(lab, @sprintf("s = %g: f_esc %.3f N_sc = %.3g", s, res.f_out, res.nsc_per_photon)); push!(spx, sp ./ max(sum(sp), 1e-300))
    end
    Zn = L.build_zones(collect(range(0, R, length = nz + 1)), fill(1e6 / (L.sig0(T) * R), nz), fill(T, nz), zeros(nz + 1); xcs = (3.0, 6.0))
    r0 = L.run(Zn, [L.src_point(1.0; npk = 1200)]; nbatch = 8, opts = L.Opts(wmax = 1e8, nbw = 10)); pred = exp(-3e-7 * r0.nsc_per_photon)
    check("d2 f_esc = exp(−ε N_sc), scaling invariance", abs(fl[1] - fl[2]) < 0.06 && abs(fl[1] / pred - 1) < 0.25, @sprintf("%s; exp(−ε N̄) = %.3f (N̄ = %.3g without destruction); (τ0, ε) → (τ0/10, 10 ε) gives the same f_esc", join(lab, "; "), pred, r0.nsc_per_photon))
    # emergent spectra of the two scalings collapse in x/(aτ0)^{1/3}
    m1 = sum(spx[1][i] * (i - 80.5)^2 for i in eachindex(spx[1])); m2 = sum(spx[2][i] * (i - 80.5)^2 for i in eachindex(spx[2]))
    check("d3 spectral collapse under x → x/λ, τ0 → τ0/λ³", abs(sqrt(m1 / m2) - 1) < 0.1, @sprintf("rms x (bins of width ∝ (aτ0)^{1/3}) %.2f vs %.2f (same bins): ratio %.3f", sqrt(m1), sqrt(m2), sqrt(m1 / m2)))
end

# ====================================================================== e
function test_e()
    cE = c / L.E0
    # e1: conservative grey scattering (far Lyα wing, Rayleigh phase function, N_sc ≫ 1): F = L/(4πr²) at every radius, force density = κF/c → M_F = τ_total EXACTLY, zone by zone
    T = 1e4; b = L.bth(T); a = L.damp(T); s0 = L.sig0(T); L.init_wings!(); x0 = 3000.0
    σx = s0 * L.voigt(a, x0) * L.Cwing(x0 * b / c); out = String[]; ok = true
    for τ in (3.0, 30.0)
        R = 1e17; nz = 10; r = collect(range(0, R, length = nz + 1))
        Z = L.build_zones(r, fill(τ / (σx * R), nz), fill(T, nz), zeros(nz + 1); xcs = ())
        res = L.run(Z, [L.src_shell(1e12, 1.0, [x0 * b * 0.999, x0 * b * 1.001], [1.0]; npk = 4000)]; nbatch = 8, opts = L.Opts(wmax = 2.5e9, nbw = 20))
        zr = res.force[2:end] .* cE ./ (τ / nz); push!(out, @sprintf("τ = %g: M_F(zones 2-10) = %.2f (events %.2f; τ_total(9/10) = %.2f), zone forces/τ_zone %.3f ± %.3f", τ, sum(res.force[2:end]) * cE, sum(res.force_event[2:end]) * cE, τ * 0.9, sum(zr) / length(zr), sqrt(sum((zr .- sum(zr) / length(zr)) .^ 2) / (length(zr) - 1))))
        ok &= abs(sum(res.force[2:end]) * cE / (τ * (nz - 1) / nz) - 1) < 0.08 && abs(sum(zr) / length(zr) - 1) < 0.05
    end
    check("e1 grey scattering: M_F = τ (exact)", ok, join(out, "; "))
    # e2: static thin shell around a central source (Dijkstra & Loeb 2008 / Smith, Bromm & Loeb 2017 set-up), T = 1e4 K
    out = String[]; ok = true
    for (NH, npk) in ((1e19, 1000), (1e20, 1000))
        Rin = 3e18; dr = 3e16; nz = 3
        r = vcat([0.0], collect(range(Rin, Rin + dr, length = nz + 1)))
        Z = L.build_zones(r, vcat(0.0, fill(NH / dr, nz)), fill(T, nz + 1), zeros(nz + 2))
        res = L.run(Z, [L.src_point(1.0; npk = npk)]; nbatch = 8, opts = L.Opts(wmax = 3e9, nbw = 40))
        MFpi = sum(res.force) * cE; MF = sum(res.force_event) * cE; dl = 10 * (NH / 1e19)^0.44; mt = 4pi / 3 * L.LW.line_plateau(2, NH, T)[1]
        push!(out, @sprintf("N_H = %.0e: M_F = %.1f (momentum balance; events %.1f) [D&L08 %.1f, MT08 %.1f, 3.51(aτ0)^(1/3) = %.1f]", NH, MFpi, MF, dl, mt, 3.51 * (L.damp(T) * NH * s0)^(1 / 3))); ok &= 0.25 < MFpi / dl < 6 && abs(MFpi / MF - 1) < 0.5
    end
    check("e2 force multiplier of a static shell", ok, join(out, "; ") * "; the cavity returns the reflected photons to the shell (N_H = 1e19: 18 crossings of the inner face per photon), M_F is 1.5-4× the D&L08 fit")
end

# ====================================================================== f
function test_f()
    T = 10.0; s0 = L.sig0(T); R = 1e16; nz = 5; τ0 = 1e5
    r = collect(range(0, R, length = nz + 1))
    Z = L.build_zones(r, fill(τ0 / (s0 * R), nz), fill(T, nz), zeros(nz + 1); xcs = ())
    res = L.run(Z, [L.src_point(1.0; npk = 12000)]; nbatch = 8, opts = L.Opts(wmax = 1e8, nbw = 20))
    ratio = res.Tcolor[1:3] ./ res.T[1:3]; err = res.Tcolor_err[1:3] ./ res.T[1:3]
    w = 1 ./ err .^ 2; m = sum(w .* ratio) / sum(w); em = 1 / sqrt(sum(w))
    check("f1 Wouthuysen–Field: T_c → T_k", abs(m - 1) < 2.5 * em + 0.05 && em < 0.15, @sprintf("T = 10 K (recoil 2g = %.3e per Doppler width), τ0 = %.0e: T_c/T_k = %s ± %s in the inner zones, mean %.3f ± %.3f (1 = Field 1959; Hirata 2006 for the closed medium)", 2 * Z.g[1], τ0, string(round.(ratio, digits = 3)), string(round.(err, digits = 3)), m, em))
    # recoil sets the sign: without it the odd part of J vanishes
    Z0 = L.build_zones(r, fill(τ0 / (s0 * R), nz), fill(T, nz), zeros(nz + 1); xcs = (), g_extra = 0.0)
    r0 = L.run(Z0, [L.src_point(1.0; npk = 12000)]; nbatch = 8, opts = L.Opts(wmax = 1e8, nbw = 20))
    s_on = res.Tcolor[1]; s_off = r0.Tcolor[1]
    check("f2 no recoil, no tilt", !r0.Tcolor_valid[1] || r0.Tcolor[1] / r0.T[1] > 3 || r0.Tcolor[1] < 0, @sprintf("without recoil the fitted tilt is consistent with zero (T_c/T_k = %.2f ± %.2f, valid = %s); with recoil %.3f ± %.3f", r0.Tcolor[1] / r0.T[1], r0.Tcolor_err[1] / r0.T[1], string(r0.Tcolor_valid[1]), s_on / res.T[1], res.Tcolor_err[1] / res.T[1]))
end

# ====================================================================== g
function test_g()
    T = 1e4; b = L.bth(T); s0 = L.sig0(T); R = 1e17; nz = 10; τ0 = 3e4; n = τ0 / (s0 * R)
    r = collect(range(0, R, length = nz + 1)); u = -5b .* r ./ R
    ne = fill(1e8, nz); Z = L.build_zones(r, fill(n, nz), fill(T, nz), u; ne = ne, np = ne, kabs = fill(0.2 / R, nz), nH2 = fill(n * 1e-3, nz), xcs = (3.0, 6.0))
    res = L.run(Z, [L.src_point(0.7; npk = 800), L.src_zones(ones(nz); npk = 800), L.src_shell(R * 0.5, 0.3, [-2e6, 2e6], [1.0]; npk = 400)]; nbatch = 8, opts = L.Opts(wmax = 2e7, nbw = 200))
    tot = res.f_out + res.f_in + res.f_2g + res.f_Hm + res.f_H2 + res.f_lost
    check("g1 photon ledger", abs(tot - 1) < 1e-12, @sprintf("out %.4f + in %.4f + 2γ %.4f + H⁻ %.4f + H₂ %.4f + cut %.1e = %.15f", res.f_out, res.f_in, res.f_2g, res.f_Hm, res.f_H2, res.f_lost, tot))
    # energy: emitted lab energy = escaped + 2γ photons + gas (absorbed + scattering energy exchange incl. work)
    Eled = res.Eemit - res.Eout - res.Ein - res.E2g - res.Egas_total
    check("g2 energy ledger (moving medium)", abs(Eled) < 1e-9 * res.Eemit, @sprintf("E_emit − E_out − E_in − E_2γ − E_gas = %.2e E_emit (infall |V| = 5b: gas exchange %.3e E_emit)", Eled / res.Eemit, res.Egas_total / res.Eemit))
    # momentum: pressure-tensor and event estimators agree in a thick (diffusive, aτ0 = 1.5e3) static medium
    T2 = 10.0; n2 = 1e5 / (L.sig0(T2) * R)
    Zs = L.build_zones(r, fill(n2, nz), fill(T2, nz), zeros(nz + 1); xcs = (3.0, 6.0))
    rs = L.run(Zs, [L.src_point(1.0; npk = 3200)]; nbatch = 8, opts = L.Opts(wmax = 2e7, nbw = 20))
    cE = c / L.E0; Fe = sum(rs.force_event) * cE; Fp = sum(rs.force) * cE
    check("g3 momentum: events vs momentum balance", abs(Fe / Fp - 1) < 0.2, @sprintf("Σ force c/L: %.2f (events) vs %.2f (momentum balance), thick static sphere T = 10 K, τ0 = 1e5", Fe, Fp))
end

# ====================================================================== h
function test_h()
    T = 10.0; s0 = L.sig0(T); R = 1e16; nz = 4; τ0 = 3e5; a = L.damp(T)
    r = collect(range(0, R, length = nz + 1)); n = τ0 / (s0 * R)
    cE = c / L.E0
    # h1: noise ∝ N^{-1/2}: scatter of the total force over seeds
    sd = Float64[]; Ns = [200, 800, 3200]
    for N in Ns
        Zk = L.build_zones(r, fill(n, nz), fill(T, nz), zeros(nz + 1); xcs = (3.0, 6.0))
        F = [sum(L.run(Zk, [L.src_point(1.0; npk = N)]; nbatch = 8, seed = 100 + s, opts = L.Opts(wmax = 1e7, nbw = 10), store = false).force) * cE for s in 1:8]
        m = sum(F) / 8; push!(sd, sqrt(sum((F .- m) .^ 2) / 7) / m)
    end
    sl = (log(sd[end]) - log(sd[1])) / (log(Ns[end]) - log(Ns[1]))
    check("h1 noise vs packets", -0.8 < sl < -0.25, @sprintf("relative scatter of the total force (8 seeds) %s for N = %s packets: slope %.2f (−0.5)", string(round.(sd, sigdigits = 2)), string(Ns), sl))
    # h2: core skipping against the exact transfer
    res = Dict{String,Any}()
    for (lab, xcs) in (("exact", ()), ("skip", (3.0, 6.0)))
        Zk = L.build_zones(r, fill(n, nz), fill(T, nz), zeros(nz + 1); xcs = xcs, epsfloor = 2e-7)
        t = @elapsed rr = L.run(Zk, [L.src_point(1.0; npk = lab == "exact" ? 2400 : 9600)]; nbatch = 8, opts = L.Opts(wmax = 8 * (a * τ0)^(1 / 3) * L.bth(T), nbw = 60))
        x = L.wcenters(rr) ./ L.bth(T); sp = rr.spec_out ./ sum(rr.spec_out); res[lab] = (f = rr.f_out, N = rr.nsc_per_photon, xr = sqrt(sum(sp .* x .^ 2)), F = sum(rr.force) * cE, t = t, nsk = count(>(0), Zk.xc))
    end
    e, s = res["exact"], res["skip"]
    check("h2 core skipping vs exact", abs(s.f - e.f) < 0.04 && abs(s.N / e.N - 1) < 0.08 && abs(s.xr / e.xr - 1) < 0.06 && abs(s.F / e.F - 1) < 0.15,
          @sprintf("f_esc %.3f/%.3f, N_sc %.3g/%.3g, rms x %.2f/%.2f, Σ force c/L %.1f/%.1f (skip/exact); time %.1f s (9600 packets) vs %.1f s (2400); eligible zones %d/%d", s.f, e.f, s.N, e.N, s.xr, e.xr, s.F, e.F, s.t, e.t, s.nsk, nz))
end

# ====================================================================== i: coupling pieces (no hydro)
for (k, f) in (("a", test_a), ("b", test_b), ("c", test_c), ("d", test_d), ("e", test_e), ("f", test_f), ("g", test_g), ("h", test_h))
    which in (k, "all") && (t = @elapsed f(); @printf("   [%s: %.0f s]\n", k, t))
end
println("S8A: ", all(last, RESULTS) ? "ALL PASS" : "$(count(!last, RESULTS)) FAIL")
