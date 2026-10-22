# S7c-2 acceptance tests: Lyman-series scattering of atomic H (LymanWings.jl in RadTransfer) and the trapped-Lyα radiation force (RadLya.jl).
#   TITAN_NOLIB=1 julia --project=env -t 4 s7c2_test.jl [a|b|c|d|e|all]
# a  cross section: sum rules, Rayleigh limit against Lee & Kim (2005) / Gavrila (1967), Lyα peak and line integral, damping wings and the region between Lyα and Lyβ against the Lee (2003) expansion, Voigt cores
# b  force of a blackbody star on a neutral shell (TITAN grid, beam) against the direct numerical integral over frequency of the exact cross section
# c  trapped-Lyα force: force multiplier against Dijkstra & Loeb (2008) / Smith, Bromm & Loeb (2017), the flux-conserving closure in the coherent regime, thin limit, bookkeeping
# d  regression: default options (lyman = lya = false) reproduce the S7b/S7a numbers to 12 digits
# e  LW band: overlap of the H I Lyman wings and the H₂ self-shielding fit (no double counting)
using Printf
const J = @__DIR__
include(joinpath(J, "PrimordialChem.jl")); include(joinpath(J, "ContinuumOpacity.jl")); include(joinpath(J, "ChemNet.jl")); include(joinpath(J, "ChemRates.jl"))
include(joinpath(J, "AtomicHHe.jl")); include(joinpath(J, "RadTransfer.jl")); include(joinpath(J, "Perf.jl")); include(joinpath(J, "s7b_exchange.jl")); include(joinpath(J, "s7b_feedback.jl"))
const CN = ChemNet; const RT = RadTransfer; const LW = RT.LymanWings; const S = S7b
const RESULTS = Tuple{String,Bool}[]
function check(name, ok, detail)
    push!(RESULTS, (name, ok)); @printf("%s  %-30s %s\n", ok ? "PASS" : "FAIL", name, detail); flush(stdout)
end
const AU = 1.495978707e13; const Msun = 1.98847e33; const Lsun = 3.828e33; const c = CN.c; const eV = CN.eV; const hP = CN.hP
which = isempty(ARGS) ? "all" : ARGS[1]
Perf.reset!("s7c2_test"; stage = "test")
net, photos = AtomicHHe.atomic()
GR = RT.default_grid(dIR = 0.3, dopt = 0.25, dLW = 0.02, dUV = 0.05, lyman = true)         # the production grid with the lines Ly-2 … Ly-12 resolved

# ====================================================================== a. cross section
function test_a()
    @printf("a: ν_Lyα = %.5e Hz (%.3f Å, %.4f eV), Lyman limit %.5e Hz; f_2 = %.5f, f_3 = %.5f, A_21 = %.4e s^-1; %d bins in the production grid (%d without lines)\n",
            LW.ν_line(2), c / LW.ν_line(2) * 1e8, hP * LW.ν_line(2) / eV, LW.νRy, LW.FN[2], LW.FN[3], LW.A_np1s(2), RT.nfreq(GR), RT.nfreq(RT.default_grid(dIR = 0.3, dopt = 0.25, dLW = 0.02, dUV = 0.05, lyman = false)))
    # a1: sum rule and static polarizability
    cont = 0.0
    for i in eachindex(LW.GLX); y = 0.5 * (LW.GLX[i] + 1); cont += 0.5 * LW.GLW[i] * LW.dfde(LW.E0 / y) * LW.E0 / y^2; end
    tot = sum(LW.FN[2:end]) + LW.dfde(LW.E0) / (2 * (LW.NLINES + 0.5)^2) + cont
    S0 = real(LW.S_total(0.0))
    check("a1 sum rule, polarizability", abs(tot - 1) < 2e-3 && abs(S0 - 4.5) < 5e-3, @sprintf("Σ f_bound = %.5f, continuum %.5f, total %.5f (1); S(0) = %.5f (4.5 a.u.)", sum(LW.FN[2:end]), cont, tot, S0))
    # a2: Rayleigh limit against Lee & Kim (2005, MNRAS 358, 1472) low-energy expansion (their Table 1, c_p), which they show to coincide with Gavrila (1967)
    cp = [1.26537, 3.73766, 8.8127, 19.1515, 39.919, 81.1018, 161.896, 319.001, 622.229, 1203.82]
    worst = 0.0; line = ""
    for rr in (0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
        s = LW.σ_coh(rr * LW.νRy); lee = LW.σT * rr^4 * sum(cp[p+1] * rr^(2p) for p in 0:9)
        worst = max(worst, abs(s / lee - 1)); line *= @sprintf(" %.2f: %.4e/%.4e", rr, s, lee)
    end
    s6 = LW.σ_coh(0.6 * LW.νRy); l6 = LW.σT * 0.6^4 * sum(cp[p+1] * 0.6^(2p) for p in 0:9)
    @printf("   σ(ω/ω_l)  [this/Lee-Kim]:%s\n   ω/ω_l = 0.6: %.4e vs series (truncated at p = 9, 'accurate to 4 %%') %.4e: ratio %.4f\n", line, s6, l6, s6 / l6)
    check("a2 Rayleigh vs Lee & Kim", worst < 2e-3 && abs(s6 / l6 - 1) < 0.06, @sprintf("max |ratio − 1| = %.1e for ω/ω_l ≤ 0.5; leading coefficient σ/σ_T (ω_l/ω)⁴ at 0.02 = %.4f (c_0 = 1.26537)", worst, LW.σ_coh(0.02 * LW.νRy) / (LW.σT * 0.02^4)))
    # a3: Ferland (2001)'s polynomial fit to Gavrila (1967) (Lee & Kim eq. 16) and Dalgarno's polynomial used before (ContinuumOpacity.σ_ray_H)
    wf = 0.0; wd = 0.0; ld = ""
    for rr in (0.2, 0.3, 0.4, 0.5, 0.6)
        s = LW.σ_coh(rr * LW.νRy)
        fer = 8.41e-25 * rr^4 + 3.37e-24 * rr^6 + 4.71e-22 * rr^14
        dal = ContinuumOpacity.σ_ray_H(c / (rr * LW.νRy) * 1e8)
        rr <= 0.6 && (wf = max(wf, abs(s / fer - 1))); wd = max(wd, abs(s / dal - 1)); ld *= @sprintf(" %.1f: %.2f", rr, dal / s)
    end
    check("a3 Ferland/Gavrila fit", wf < 0.07, @sprintf("max |σ/σ_Ferland − 1| = %.3f (ω/ω_l = 0.2-0.6); Dalgarno's polynomial / this (ω/ω_l:ratio):%s", wf, ld))
    # a4: resonance: peak 3λ²/2π, line integral π e² f/(m c)
    ν0 = LW.ν_line(2); pk = LW.σ_coh(ν0); an = 3 * (c / ν0)^2 / (2π)
    D = LW.A_tot(2) / (4π)
    Iw = LW.geo_integral(ν -> LW.σ_res(ν, 2), ν0 * (1 - 1e-3), ν0 * (1 + 1e-3), ν0, 1e3); tail = 2 * LW.πe2mc * LW.FN[2] / π * D / (1e-3 * ν0)
    check("a4 Lyα resonance", abs(pk / an - 1) < 1e-4 && abs((Iw + tail) / (LW.πe2mc * LW.FN[2]) - 1) < 5e-3,
          @sprintf("peak σ = %.4e (3λ²/2π = %.4e); ∫σ dν over ±1e-3 ν₀ + Lorentz tail = %.4e (π e² f/mc = %.4e)", pk, an, Iw + tail, LW.πe2mc * LW.FN[2]))
    # a5: wings: Lee (2003) fifth-order expansion about Lyα (Lee & Kim eqs. 17-19), both sides of the line (red wing and the region between Lyα and Lyβ)
    A = [1.0, -0.8961, -12.22, -52.52, -243.8, -1210.0]
    # eq. (19) as printed: polynomial of the square truncated at the fifth order
    lee19(ωt) = 0.0433056 / ωt^2 * (1 - 1.792ωt - 23.637ωt^2 - 83.1393ωt^3 - 244.1453ωt^4 - 699.473ωt^5)
    ws = 0.0; ln = ""
    for ωt in (-0.1, -0.05, -0.02, -0.01, 0.01, 0.02, 0.05)
        s = LW.σ_coh(ν0 * (1 + ωt)) / LW.σT; l = lee19(ωt)
        a_ = abs(ωt) <= 0.05 && (ws = max(ws, abs(s / l - 1))); ln *= @sprintf(" %+.2f: %.3e/%.3e", ωt, s, l)
    end
    @printf("   σ/σ_T near Lyα [this/Lee eq. 19]:%s\n", ln)
    # between Lyα and Lyβ: the sum over lines, Lorentz wings of Lyα and Lyβ and the continuum: spot values
    nu_mid = 0.5 * (LW.ν_line(2) + LW.ν_line(3))
    wing_a = LW.πe2mc * LW.FN[2] * (LW.A_tot(2) / (4π^2)) / (nu_mid - LW.ν_line(2))^2; wing_b = LW.πe2mc * LW.FN[3] * (LW.A_tot(3) / (4π^2)) / (LW.ν_line(3) - nu_mid)^2
    @printf("   midway Lyα–Lyβ (%.4f eV): σ = %.4e cm² (Lorentz wings of Lyα %.3e + Lyβ %.3e; the dispersion terms carry the rest)\n", hP * nu_mid / eV, LW.σ_coh(nu_mid), wing_a, wing_b)
    check("a5 wings vs Lee expansion", ws < 0.03, @sprintf("max |σ/σ_Lee − 1| = %.3f for |Δω/ω₂₁| ≤ 0.05 (fifth order; the series itself is good to a few %% there)", ws))
    # a6: Voigt function and thermal cores
    function voigt_ref(a, x)
        # H(a, x) = (1/π) ∫_{−π/2}^{π/2} exp(−(x + a tan φ)²) dφ  (y = x + a tan φ removes the Lorentzian peak), composite Gauss–Legendre
        s = 0.0; nsub = 400; h = π / nsub
        for i in 1:nsub, k in 1:5
            φ = -π / 2 + (i - 1) * h + 0.5 * h * (LW.G5X[k] + 1); s += 0.5 * h * LW.G5W[k] * exp(-(x + a * tan(φ))^2)
        end
        s / π
    end
    wv = 0.0
    for a_ in (1e-3, 0.05, 0.5, 2.0), x in (0.0, 0.7, 2.0)
        wv = max(wv, abs(LW.voigt_H(a_, x) - voigt_ref(a_, x)) / voigt_ref(a_, x))
    end
    for a_ in (0.5, 2.0); wv = max(wv, abs(LW.voigt_H(a_, 5.0) - voigt_ref(a_, 5.0)) / voigt_ref(a_, 5.0)); end
    for a_ in (1e-3, 0.05, 0.5, 2.0), x in (30.0, 300.0)          # the damping wing: H → a/(√π x²) (1 + (3/2 − a²)/x² + …)
        w_ = a_ / (sqrt(π) * x^2) * (1 + (1.5 - a_^2) / x^2); wv = max(wv, abs(LW.voigt_H(a_, x) - w_) / w_)
    end
    # area conservation of the thermal convolution over a window about Lyα (the same window for both), and the core peak
    T = 1e4; w = 3e-3
    I1 = LW.geo_integral(ν -> LW.σ_gas(ν, T), ν0 * (1 - w), ν0 * (1 + w), ν0, 0.25 * LW.Δν_D(2, T))
    I0 = LW.geo_integral(ν -> LW.σ_coh(ν), ν0 * (1 - w), ν0 * (1 + w), ν0, 0.5 * LW.A_tot(2) / (4π))      # (the atom's Lorentz core is 5e7 Hz wide)
    pkT = LW.σ_gas(ν0, T); pkan = LW.πe2mc * LW.FN[2] / (sqrt(π) * LW.Δν_D(2, T)) * LW.voigt_H(LW.A_tot(2) / (4π * LW.Δν_D(2, T)), 0.0)
    check("a6 Voigt, thermal cores", wv < 2e-3 && abs(I1 / I0 - 1) < 2e-3 && abs(pkT / pkan - 1) < 2e-3,
          @sprintf("Voigt H(a,x) vs quadrature: max rel. diff %.1e; ∫σ(T = 1e4 K) / ∫σ(atom) over ±0.3 %% = %.5f; core peak %.4e cm² (σ₀ H(a,0) = %.4e)", wv, I1 / I0, pkT, pkan))
    # a7: bin averages: table against brute-force sampling, temperature interpolation
    tab = RT.lytable(GR); worst = 0.0; nb = 0; worstbin = 0; sabs = 0.0; stot = 0.0
    for i in 1:RT.nfreq(GR)
        a_, b_ = GR.edges[i], GR.edges[i+1]
        (a_ > 0.7 * ν0 && b_ < LW.ν_line(12) * 1.001) || continue         # (above Ly-12 the bins hold the unresolved lines n > 12)
        Tt = 1500.0
        s_tab = LW.sigma_bin(tab, i, Tt)
        pts = sort(unique(vcat(range(a_, b_, length = 600), [LW.ν_line(n) + d * LW.ν_line(n) * q for n in 2:12 for d in (-3e-7, -1e-6, -1e-5, 0.0, 1e-6, 1e-5, 3e-7) for q in (1.0,) if a_ < LW.ν_line(n) * (1 + d) < b_])))
        s_ref = sum(LW.σ_gas(0.5 * (pts[q] + pts[q+1]), Tt) * (pts[q+1] - pts[q]) for q in 1:length(pts)-1) / (b_ - a_)
        rel = abs(s_tab / s_ref - 1); nb += 1; sabs += abs(s_tab - s_ref) * (b_ - a_); stot += s_ref * (b_ - a_)
        if rel > worst; worst = rel; worstbin = i; end
    end
    @printf("   table vs brute-force bin average at T = 1500 K (between nodes): %d bins, worst %.3f in bin %d (%.4f eV, width %.1e)\n", nb, worst, worstbin, GR.E[worstbin], (GR.edges[worstbin+1] - GR.edges[worstbin]) / GR.ν[worstbin])
    check("a7 tabulated bin averages", worst < 0.25 && sabs / stot < 0.02, @sprintf("max rel. difference %.3f (bins in the Gaussian-to-Lorentz transition of a core, 4-7 Doppler widths, where ln σ is not linear in ln T between the nodes 10^(1/8) apart); width-weighted mean %.1e", worst, sabs / stot))
end

# ====================================================================== b. force on a neutral shell
function shell_force(g, N, T; Teff = 5e4, L = 1e6 * Lsun, r0 = 1e16, o = RT.Opts(diffuse = false, lyman = true, lya = false))
    ab = RT.absorbers(photos, g)
    nH = 1e5; dr = N / nH; r = [r0, r0 + dr]
    X = Dict{Symbol,Vector{Float64}}(:H => [1.0], :H2 => [0.0], :He => [0.0], :e => [1e-8], :Hp => [1e-8])
    med = RT.Medium(r, [T], [nH], X)
    src = RT.blackbody_source(g, Teff, L)
    f = RT.transfer(g, ab, src, med, o; rho = [nH * CN.μH()])
    sum((f.Ldir[1, i] - f.Ldir[2, i]) for i in 1:RT.nfreq(g) if g.edges[i+1] <= 13.59 * eV / hP) / L       # momentum (units L/c) removed from the beam below the Lyman limit
end
function direct_force(N, T; Teff = 5e4, L = 1e6 * Lsun, nν = 20000, Emin = 0.1, Emax = 13.59)
    νs = exp.(range(log(Emin * eV / hP), log(Emax * eV / hP), length = nν))
    νs = sort(unique(vcat(νs, filter(x -> x > νs[1] && x < νs[end], LW.line_edges(2, 12; dmin = 3e-7, per_decade = 60, dmax = 0.1)))))
    Lν(ν) = L * π * RT.Bν(ν, Teff) / (RT.σSB * Teff^4)
    s = 0.0
    for i in 1:length(νs)-1
        ν = 0.5 * (νs[i] + νs[i+1]); s += Lν(ν) * (-expm1(-N * LW.σ_gas(ν, T))) * (νs[i+1] - νs[i])
    end
    s / L
end
function test_b()
    worst = 0.0; lines = String[]
    for N in (1e18, 1e20, 1e22, 1e24)
        fg = shell_force(GR, N, 1000.0); fd = direct_force(N, 1000.0)
        push!(lines, @sprintf("   N_H = %.0e cm^-2: grid %.5e, direct integral %.5e (fraction of L/c removed below 13.59 eV), ratio %.4f", N, fg, fd, fg / fd))
        N >= 1e20 && (worst = max(worst, abs(fg / fd - 1)))
    end
    foreach(println, lines)
    f0 = shell_force(GR, 1e22, 1000.0; o = RT.Opts(diffuse = false, lyman = false, lya = false))
    check("b1 shell force vs integral", worst < 0.05, @sprintf("5e4 K blackbody, neutral shell at T = 1000 K: max |grid/direct − 1| = %.3f for N_H ≥ 1e20 (Voigt cores and Lorentz wings; the grid resolves Ly-2 … Ly-12); with Dalgarno's polynomial instead the force at 1e22 is %.4e (%.1f times smaller)", worst, f0, shell_force(GR, 1e22, 1000.0) / f0))
end

# ====================================================================== c. trapped Lyα force
function test_c()
    # c1: force multiplier of a thin shell, M_F = (4π/3) J/F = (4π/3) 0.411 (aτ̄)^{1/3} (sat, no destruction), T = 1e4 K (Δv_D = 12.9 km/s), against the fit of Dijkstra & Loeb (2008) Fig. 6 quoted by Smith, Bromm & Loeb (2017)
    mf(N) = 4π / 3 * LW.line_plateau(2, N, 1e4)[1]
    fit(N) = 10 * (N / 1e19)^0.44
    ratios = [mf(N) / fit(N) for N in (1e19, 1e20, 1e21)]
    slope = log(mf(1e21) / mf(1e19)) / log(100)
    τ0(N) = LW.πe2mc * LW.FN[2] / (sqrt(π) * LW.Δν_D(2, 1e4)) * N
    slab(N) = 15 * (τ0(N) / 10^5.5)^(1 / 3)
    @printf("   M_F(N_H = 1e19, 1e20, 1e21, 1e22) = %.1f %.1f %.1f %.1f; fit 10 N19^0.44 = %.1f %.1f %.1f %.1f; slab estimate 15 (τ_0/10^5.5)^{1/3} = %.1f %.1f %.1f\n", mf(1e19), mf(1e20), mf(1e21), mf(1e22), fit(1e19), fit(1e20), fit(1e21), fit(1e22), slab(1e19), slab(1e20), slab(1e21))
    check("c1 force multiplier", all(0.7 .< ratios .< 1.5) && 0.30 < slope < 0.46, @sprintf("M_F / (10 N19^0.44) = %.2f %.2f %.2f at N_H = 1e19, 1e20, 1e21; log-slope %.3f (analytic 1/3; their fit 0.44)", ratios..., slope))
    # c2: the flux-conserving closure in the coherent regime: monochromatic source at 8 eV, ρ ∝ r^-1.5 envelope, force per mass σ n_H F/(ρ c) with F = L/(4π r²), no saturation of the beam
    E = 8.0; ν = E * eV / hP
    g = RT.edge_grid(range(7.5, 8.5, length = 6)); tabg = RT.lytable(g)
    n = 60; r0 = 1e16; r = r0 .* exp.(range(0, log(100.0), length = n + 1)); rc = sqrt.(r[1:end-1] .* r[2:end])
    σ8 = LW.σ_gas(ν, 300.0)
    τtot = 40.0; n0 = τtot / σ8 / (2 * r0)                       # N_out(r0) = 2 n0 r0 (1 − (r0/r1)^½)
    nH = n0 .* (rc ./ r0) .^ -1.5
    X = Dict{Symbol,Vector{Float64}}(:H => ones(n), :e => fill(1e-9, n), :Hp => fill(1e-9, n))
    med = RT.Medium(r, fill(300.0, n), nH, X); rho = nH .* CN.μH()
    ab = RT.absorbers(photos, g)
    L0 = 1e38; src = RT.mono_source(g, E, L0 / (hP * ν))
    fon = RT.transfer(g, ab, src, med, RT.Opts(diffuse = false, lyman = true, lya = true); rho = rho)
    foff = RT.transfer(g, ab, src, med, RT.Opts(diffuse = false, lyman = true, lya = false); rho = rho)
    σb = LW.sigma_bin(tabg, RT.ibin(g, E), 300.0)
    # the flux that is conserved is the one attenuated by true absorption only (here the weak H₂⁺/H⁻ opacities of the continuum module: τ_abs ≈ 0.25 over the envelope)
    op = RT.new_opacity(n, RT.nfreq(g), length(ab)); RT.fill_opacity!(op, g, ab, med, RT.Opts(diffuse = false, lyman = true, lya = true))
    fb = RT.ibin(g, E); Lj = zeros(n + 1); Lj[1] = L0; for k in 1:n; Lj[k+1] = Lj[k] * exp(-op.kabs[k, fb] * (r[k+1] - r[k])); end
    expect = [σb * nH[k] * sqrt(Lj[k] * Lj[k+1]) / (4π * rc[k]^2) / (rho[k] * c) for k in 1:n]
    ks = 5:n-5
    ratio = fon.arad[ks] ./ expect[ks]
    sat = foff.arad[ks] ./ expect[ks]
    @printf("   coherent regime (τ_total(r_0) = %.0f at 8 eV): a_lya / (σ n F/ρc) at zones 5, 20, 40, 55 = %.3f %.3f %.3f %.3f; the beam alone (saturating): %.3g %.3g %.3g %.3g\n", τtot,
            fon.arad[5] / expect[5], fon.arad[20] / expect[20], fon.arad[40] / expect[40], fon.arad[55] / expect[55], foff.arad[5] / expect[5], foff.arad[20] / expect[20], foff.arad[40] / expect[40], foff.arad[55] / expect[55])
    check("c2 flux-conserving closure", all(0.95 .< ratio .< 1.05) && sat[1] < 0.2, @sprintf("force per mass / (σ n_H F/ρc): mean %.3f, range %.3f-%.3f over zones 5-55 (conserved flux, N_eff from the profile); the beam (scattered photons removed) gives %.3f at zone 5", sum(ratio) / length(ratio), minimum(ratio), maximum(ratio), sat[1]))
    # c3: thin limit: the component vanishes when nothing is trapped
    medt = RT.Medium(r, fill(300.0, n), nH .* 1e-9, X); rhot = nH .* 1e-9 .* CN.μH()
    ft_on = RT.transfer(g, ab, src, medt, RT.Opts(diffuse = false, lyman = true, lya = true); rho = rhot)
    ft_off = RT.transfer(g, ab, src, medt, RT.Opts(diffuse = false, lyman = true, lya = false); rho = rhot)
    dth = maximum(abs.(ft_on.arad .- ft_off.arad) ./ ft_off.arad)
    check("c3 thin limit", dth < 0.15, @sprintf("τ_total = %.1e: max |a(lya) − a(beam)|/a(beam) = %.3f (bookkeeping of N_eff only)", τtot * 1e-9, dth))
    # c4: bookkeeping of the pressure of an injected line: P(1) = (4π/3c) (J/F) F, F = L/(4π r²), N_eff = N_out (k_ρ − 1)/(k_ρ + 1)
    Linj = 1e39; o = RT.Opts(diffuse = false, lyman = true, lya = true)
    Tz = fill(1000.0, n); nHl = 1e8 .* (rc ./ r0) .^ -1.5
    X2 = Dict{Symbol,Vector{Float64}}(:H => ones(n), :H2 => zeros(n), :He => zeros(n), :e => fill(1e-8, n), :Hp => fill(1e-8, n))
    medl = RT.Medium(r, Tz, nHl, X2)
    La = zeros(n + 1, RT.nfreq(GR))
    P, parts = RT.lya_pressures(GR, RT.lytable(GR), La, medl, o; Linject = Linj)
    Nout1 = sum(nHl[k] * (r[k+1] - r[k]) for k in 1:n); Neff = parts.Neff[1]
    new_ = sum(nHl[k] * (r[k+1] - r[k]) * 1e-8 * nHl[k] for k in 1:n) / Nout1          # column-weighted n_e of the neutral gas
    Jl = LW.line_plateau(2, Neff, 1000.0; eps = 4.4e-9 / (1 + 1.55e4 / new_))
    pan = 4π / (3c) * Jl[1] * Linj / (4π * r[1]^2)
    @printf("   injected Lyα (L = 1e39 erg/s) at r = %.0f AU in an n_H = 1e8 (r/r0)^-1.5 envelope: N_out = %.3e, N_eff = %.3e (k_ρ = %.2f), J/F = %.1f, M_F = %.1f, P_iso = %.4e dyn/cm² (formula %.4e)\n", r[1] / AU, Nout1, Neff, parts.kρ[1], Jl[1], 4π / 3 * Jl[1], P[1], pan)
    check("c4 injected line", abs(P[1] / pan - 1) < 1e-9 && P[n+1] == 0 && all(diff(P[1:n+1]) .<= 1e-30 * P[1]), @sprintf("P_iso(r_1) / formula = %.12f; P = 0 at the outer edge, decreasing outwards", P[1] / pan))
end

# ====================================================================== d. regression
function test_d()
    M, R, Ls, Linf, Tph = 100.0, 41.0, 2.5e6, 3.4e6, 8200.0
    g = RT.default_grid(dIR = 0.3, dopt = 0.25, dLW = 0.02, dUV = 0.05, lyman = false)
    ia = RT.ion_atoms(photos, g; dr = AtomicHHe.k_di_HeII)
    rec = Dict("Lstar" => Ls * S.Lsun, "R" => R * S.Rsun, "Linf" => Linf * S.Lsun, "Tph" => Tph)
    src0 = S.star_source(g, rec)
    inn = S.inner_column(ia, src0, M * S.Msun, R * S.Rsun, 500 * S.AU, 4.7e-16, 1.9e6; reemit = (g, Tph))
    src = RT.Source(inn.Lout, 0.0, sum(inn.Lout))
    ab = RT.absorbers(photos, g)
    n = 80; r = 500 * S.AU .* exp.(range(0, log(400.0), length = n + 1)); rc = sqrt.(r[1:end-1] .* r[2:end])
    nH = 1e7 .* (rc ./ (3000 * S.AU)) .^ -1.5
    X = Dict{Symbol,Vector{Float64}}(:H => ones(n), :H2 => fill(1e-3, n), :He => fill(CN.yHe, n), :e => fill(1e-5, n), :Hp => fill(1e-5, n))
    med = RT.Medium(r, fill(1000.0, n), nH, X); rho = nH .* CN.μH()
    # reference values of the S7b code (RadTransfer.jl at commit 1da6d40), 12 digits
    ref = Dict(false => (1.727178251030e-04, 8.376734081802e-06, 3.002954401616e-08, 5.942762103282e-11, 5.456841205430e-04, 2.715295247994e-05),
               true => (-2.415254042752e-01, 1.801703955696e-01, 3.203935376045e-07, 4.321786963215e-10, 6.684575261536e-01, 1.858380851446e-04))
    worst = 0.0
    for diff in (false, true)
        f = RT.transfer(g, ab, src, med, RT.Opts(diffuse = diff); rho = rho)
        got = (f.arad[1], f.arad[10], f.arad[40], f.arad[80], sum(f.arad), f.kLW[1])
        worst = max(worst, maximum(abs.(got .- ref[diff]) ./ abs.(ref[diff])))
    end
    check("d1 defaults = S7b", worst < 1e-10 && !RT.LYMAN[] && !RT.LYA[] && all(RT.ARAD_LYA[] .== 0), @sprintf("a_rad (zones 1, 10, 40, 80, Σ) and k_LW of the 100 Msun star in a 1e7 cm^-3 (r/3000 AU)^-1.5 envelope, direct beam and with the diffuse field: max relative difference to the S7b numbers %.1e; ARAD_LYA = 0", worst))
    # with lyman = lya = false but the grid refined by the lines: the same physics (Dalgarno), only more bins
    o_on = RT.Opts(diffuse = false, lyman = true, lya = false)
    g1 = RT.default_grid(dIR = 0.3, dopt = 0.25, dLW = 0.02, dUV = 0.05, lyman = true); ia1 = RT.ion_atoms(photos, g1; dr = AtomicHHe.k_di_HeII)
    src1 = S.star_source(g1, rec); inn1 = S.inner_column(ia1, src1, M * S.Msun, R * S.Rsun, 500 * S.AU, 4.7e-16, 1.9e6; reemit = (g1, Tph))
    ab1 = RT.absorbers(photos, g1)
    f0 = RT.transfer(g1, ab1, RT.Source(inn1.Lout, 0.0, sum(inn1.Lout)), med, RT.Opts(diffuse = false, lyman = false, lya = false); rho = rho)
    f1 = RT.transfer(g1, ab1, RT.Source(inn1.Lout, 0.0, sum(inn1.Lout)), med, o_on; rho = rho)
    f2 = RT.transfer(g1, ab1, RT.Source(inn1.Lout, 0.0, sum(inn1.Lout)), med, RT.Opts(diffuse = false, lyman = true, lya = true); rho = rho)
    f0g = RT.transfer(g, ab, src, med, RT.Opts(diffuse = false); rho = rho)
    gl = S.G * M * S.Msun ./ rc .^ 2
    @printf("   a_rad/g at zones 1, 10, 20, 40, 80: S7b (coarse grid) %s | refined grid, Dalgarno %s | + Lyman wings %s | + trapped Lyα %s\n",
            join([@sprintf("%.3g", f0g.arad[k] / gl[k]) for k in (1, 10, 20, 40, 80)], " "), join([@sprintf("%.3g", f0.arad[k] / gl[k]) for k in (1, 10, 20, 40, 80)], " "),
            join([@sprintf("%.3g", f1.arad[k] / gl[k]) for k in (1, 10, 20, 40, 80)], " "), join([@sprintf("%.3g", f2.arad[k] / gl[k]) for k in (1, 10, 20, 40, 80)], " "))
    check("d2 grid refinement alone", maximum(abs.(f0.arad .- f0g.arad) ./ f0g.arad) < 0.08, @sprintf("lines resolved, Dalgarno opacity: max |a − a(S7b grid)|/a = %.1e (bin splitting only: the coarse 25 %% bins of S7b sample the spectrum and the opacities at their centres)", maximum(abs.(f0.arad .- f0g.arad) ./ f0g.arad)))
end

# ====================================================================== e. LW band
function test_e()
    # H I wings and cores sit on top of the Lyman-Werner band (11.2-13.6 eV). The H₂ lines of the band are NOT in the opacity (kLW, self-shielding fit f_sh(N_H₂, T)); the H I scattering is. The two are
    # independent attenuations of the same photons: J_LW(zone) = J(beam attenuated by H I wings, Rayleigh, …) and k_LW = 1.38e9 J_LW f_sh(N_H₂): no photon is counted twice.
    g = GR; ab = RT.absorbers(photos, g)
    out = String[]; ratios = Float64[]
    for N in (1e20, 1e22, 1e23)
        nH = 1e5; r = [1e16, 1e16 + N / nH]
        X = Dict{Symbol,Vector{Float64}}(:H => [1.0], :H2 => [1e-4], :He => [0.0], :e => [1e-8], :Hp => [1e-8])
        med = RT.Medium(r, [1000.0], [nH], X)
        src = RT.blackbody_source(g, 3e4, 1e6 * Lsun)
        k0 = RT.transfer(g, ab, src, med, RT.Opts(diffuse = false, lyman = false); rho = [nH * CN.μH()]).kLW[1]
        k1 = RT.transfer(g, ab, src, med, RT.Opts(diffuse = false, lyman = true); rho = [nH * CN.μH()]).kLW[1]
        push!(ratios, k1 / k0); push!(out, @sprintf("N_H = %.0e: k_LW %.3e (Dalgarno) -> %.3e (Lyman wings), ratio %.3f", N, k0, k1, k1 / k0))
    end
    foreach(l -> println("   ", l), out)
    check("e1 LW overlap", all(0.5 .< ratios .< 2.0), "k_LW changes by the difference of the two H I scattering cross sections in 11.2-13.6 eV (Dalgarno's polynomial is an extrapolation there, the KH sum is not): the H₂ Lyman/Werner lines are not in the opacity (kLW = 1.38e9 J_LW f_sh), the H I wings attenuate the same photons independently; no photon is counted twice")
end

which in ("a", "all") && Perf.@section :a test_a()
which in ("b", "all") && Perf.@section :b test_b()
which in ("c", "all") && Perf.@section :c test_c()
which in ("d", "all") && Perf.@section :d test_d()
which in ("e", "all") && Perf.@section :e test_e()
Perf.write()
println("S7C2: ", all(last, RESULTS) ? "ALL PASS" : "$(count(!last, RESULTS)) FAIL")
