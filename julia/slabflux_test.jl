# SLABFLUX unit test (TitanResidFast.sf_mean / sf_avchi!, docs/POPIII_SLABFLUX.md):
#   1. limits and ordering of the interface opacity χ_n (equal χ, symmetry, H ≤ χ_∞ ≤ G = χ_2 ≤ χ_1 = logarithmic mean, monotone in n,
#      bounded by max χ±) and smoothness of the series branch of ℰ;
#   2. (a) against a fine-grid numerical solution of the two-zone slab: a grey LTE layer at constant flux with χ ∝ T^q across a
#      transparency front (χ falls by up to 1e5 between the centres), shooting on 2e5 graded RK4 cells; (b) a sharp opacity step at
#      a position φ between the centres, to show which front position each closure implies;
#   3. AD vs central finite differences of avchi on a 4-zone strip (elasticities in χ), including the series branch.
#   TITAN_NOLIB=1 julia --project=<env> slabflux_test.jl
using Printf, ForwardDiff
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl"))
const TX = TitanResidFast
nbad = Ref(0)
check(ok, msg) = (ok || (nbad[] += 1); println(ok ? "  ok   " : "  FAIL ", msg))
H(a, b) = 2 / (1 / a + 1 / b); G(a, b) = sqrt(a * b); Lm(a, b) = a == b ? a : (a - b) / log(a / b)
χn(a, b, n) = TX.sf_mean(a, b, n)

println("1. limits and ordering")
for χ in (1e-7, 3.1e-5, 1.0, 21.7, 2e4)
    check(all(abs(χn(χ, χ, n) / χ - 1) < 4e-16 for n in (1.0, 2.0, 3.25, 10.0, Inf)), @sprintf("equal chi %.2e: chi_n = chi for n = 1, 2, 3.25, 10, Inf", χ))
end
let worst = 0.0, ordok = true, symok = true, mono = true, bound = true
    for a in 10.0 .^ (-6:0.5:4), b in 10.0 .^ (-6:0.5:4)
        c1 = χn(a, b, 1.0); c2 = χn(a, b, 2.0); ci = χn(a, b, Inf)
        worst = max(worst, abs(c1 / Lm(a, b) - 1), abs(c2 / G(a, b) - 1))
        ordok &= H(a, b) * (1 - 1e-13) <= ci <= c2 * (1 + 1e-13) <= c1 * (1 + 2e-13)
        symok &= all(abs(χn(a, b, n) / χn(b, a, n) - 1) < 1e-13 for n in (1.0, 1.7, 3.25, Inf))
        vs = [χn(a, b, n) for n in (1.0, 1.5, 2.0, 3.0, 5.0, 20.0, 1e3, Inf)]
        mono &= all(vs[i+1] <= vs[i] * (1 + 1e-12) for i in 1:length(vs)-1)
        bound &= all(min(a, b) * (1 - 1e-13) <= v <= max(a, b) * (1 + 1e-13) for v in vs)
    end
    check(worst < 1e-12, @sprintf("n = 1 is the logarithmic mean (exponential in column), n = 2 the geometric mean: worst rel. difference %.1e", worst))
    check(ordok, "H <= chi_Inf <= G = chi_2 <= chi_1 = L on a 21x21 grid of chi- , chi+ in 1e-6..1e4")
    check(symok, "chi_n(chi-, chi+) = chi_n(chi+, chi-)")
    check(mono, "chi_n decreasing in n (n = 1 ... Inf)")
    check(bound, "min chi <= chi_n <= max chi")
end
# the series branch of ℰ (|y| < 1e-2): value and first two derivatives continuous across the switch
let e0 = 0.0, e1 = 0.0, e2 = 0.0
    for y0 in (1e-2, -1e-2)
        yi = y0 * (1 - 1e-12); yo = y0 * (1 + 1e-12)
        d1(y) = ForwardDiff.derivative(TX.sf_E, y); d2(y) = ForwardDiff.derivative(d1, y)
        ex(y) = expm1(y) / y; ex1(y) = (exp(y) * y - expm1(y)) / y^2
        Δ = yo - yi       # jumps net of the slope across the 2e-12 relative gap
        e0 = max(e0, abs(TX.sf_E(yo) - TX.sf_E(yi) - d1(yi) * Δ), abs(TX.sf_E(yi) - ex(yi)))
        e1 = max(e1, abs(d1(yo) - d1(yi) - d2(yi) * Δ), abs(d1(yi) - ex1(yi)))
        e2 = max(e2, abs(d2(yi) - d2(yo)))
    end
    check(e0 < 1e-15 && e1 < 1e-12 && e2 < 1e-8, @sprintf("E(y) series branch at |y| = 1e-2: jumps of value %.1e, 1st derivative %.1e, 2nd %.1e", e0, e1, e2))
end

println("2a. fine-grid grey LTE layer at constant flux, chi = k0 T^q (front: chi falls by R between the centres)")
# column ΔM between the centres, Eddington f = 1/3 (thick layer): dE/dm = −3 χ(E) F/c, χ(E) = χ₊ (E/E₊)^(q/4).
# Given E₋, E₊ the exact F follows from shooting; each closure predicts F = c (E₋ − E₊)/(3 ΔM χ_eff) from the centre values alone.
const cc = 2.998e10
function fine_flux(E₋, E₊, χ₊, p, ΔM; N = 200_000, k = 6)
    χ(E) = χ₊ * (E / E₊)^p
    rhs(E, F) = -3χ(max(E, 1e-300)) * F / cc
    # RK4 on a mesh graded towards the opaque centre (m_i = ΔM (i/N)^k): the opaque part of the layer is a tiny column
    shoot(F) = (E = E₋; for i in 1:N
        h = ΔM * ((i / N)^k - ((i - 1) / N)^k)
        k1 = rhs(E, F); k2 = rhs(E + h / 2 * k1, F); k3 = rhs(E + h / 2 * k2, F); k4 = rhs(E + h * k3, F)
        E = max(E + h / 6 * (k1 + 2k2 + 2k3 + k4), 1e-300); end; E)
    lo, hi = 0.0, cc * (E₋ - E₊) / (3ΔM * χ₊)        # F with the transparent χ everywhere is an upper bound
    for _ in 1:80; mid = (lo + hi) / 2; shoot(mid) > E₊ ? (lo = mid) : (hi = mid); end
    (lo + hi) / 2
end
for (q, R, b) in ((13.0, 1e3, 0.3), (13.0, 1e5, 0.3), (10.0, 1e5, 1.0), (8.0, 1e4, 0.05))
    p = q / 4; χ₊ = 3.1e-5; χ₋ = R * χ₊; E₊ = 1.0; E₋ = E₊ * R^(1 / p); ΔM = 2b / χ₊   # b = χ₊ ΔM/2: the thin half-zone's depth
    Ff = fine_flux(E₋, E₊, χ₊, p, ΔM)
    Fk = cc * (E₋^(1 - p) - E₊^(1 - p)) / (χ₊ * E₊^(-p) * (1 - p)) / (3ΔM)     # Kirchhoff (closed form of the same problem)
    pred(χe) = cc * (E₋ - E₊) / (3ΔM * χe)
    Fs = pred(χn(χ₋, χ₊, p)); Fi = pred(χn(χ₋, χ₊, Inf)); Fh = pred(H(χ₋, χ₊)); Fg = pred(G(χ₋, χ₊)); Fl = pred(Lm(χ₋, χ₊))
    check(abs(Ff / Fk - 1) < 1e-6, @sprintf("q %4.1f R %.0e b %.2f: fine grid vs Kirchhoff closed form %.1e", q, R, b, Ff / Fk - 1))
    check(abs(Fs / Ff - 1) < 1e-6, @sprintf("   SLABFLUX n = q/4 = %.2f exact: F_SF/F_fine - 1 = %.1e", p, Fs / Ff - 1))
    @printf("         F/F_fine: harmonic %.3g  SLABFLUX n=Inf %.3g  geometric %.3g  log mean %.3g   (dtau_fine %.3g, harmonic %.3g, n=Inf %.3g)\n",
            Fh / Ff, Fi / Ff, Fg / Ff, Fl / Ff, cc * (E₋ - E₊) / (3Ff), ΔM * H(χ₋, χ₊), ΔM * χn(χ₋, χ₊, Inf))
end

println("2b. sharp opacity step (zdiagB 67|68 values) at fraction phi of the slab from the opaque centre: dtau = dM (phi chi- + (1 - phi) chi+)")
let χ₋ = 21.7, χ₊ = 3.1e-5, ΔM = 2 * 0.42 / 3.1e-5
    for φ in (0.0, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5)
        @printf("  phi %.0e: dtau %.4g\n", φ, ΔM * (φ * χ₋ + (1 - φ) * χ₊))
    end
    φof(χe) = (χe - χ₊) / (χ₋ - χ₊)
    @printf("  closures: harmonic dtau %.4g (phi %.1e)  SLABFLUX n=Inf %.4g (phi %.1e)  n=3 %.4g (phi %.1e)  geometric %.4g (phi %.1e)\n",
            ΔM * H(χ₋, χ₊), φof(H(χ₋, χ₊)), ΔM * χn(χ₋, χ₊, Inf), φof(χn(χ₋, χ₊, Inf)), ΔM * χn(χ₋, χ₊, 3.0), φof(χn(χ₋, χ₊, 3.0)), ΔM * G(χ₋, χ₊), φof(G(χ₋, χ₊)))
end

println("3. AD vs FD of sf_avchi! on a 4-zone strip (elasticities d ln avchi_k / d ln chi_j)")
TX.SLABFLUX[] = true
function strip(χ, n)
    T = eltype(χ); TX.SFN[] = n
    avf = vcat(zero(T), T[2 / (1 / χ[k-1] + 1 / χ[k]) for k in 2:4], zero(T))
    TX.sf_avchi!(avf, χ, 1, 4); avf[2:4]
end
cases = (("front (zdiagB 67|68: hot thick, cold thin)", [543.0, 21.7, 3.1e-5, 2.2e-5]),
         ("handover (cold zone at tau ~ 1)",            [543.0, 21.7, 3.0e-4, 2.2e-5]),
         ("thick interior",                             [2e4, 1.4e4, 1.0e4, 7.0e3]),
         ("series branch (|ln chi ratio| < 1e-2)",      [1.0, 1.004, 1.0, 0.997]),
         ("thin envelope",                              [2.6e-5, 2.1e-5, 1.9e-5, 1.8e-5]))
for n in (Inf, 3.0), (name, p) in cases
    f(x) = log.(strip(exp.(x), n)); x0 = log.(p)
    J = ForwardDiff.jacobian(f, x0); worst = 0.0; h = 1e-6
    for j in eachindex(x0)
        xp = copy(x0); xm = copy(x0); xp[j] += h; xm[j] -= h
        worst = max(worst, maximum(abs.((f(xp) .- f(xm)) ./ 2h .- J[:, j])))
    end
    check(worst < 1e-7, @sprintf("n %-4s %-44s worst |AD - FD| %.1e (max |elasticity| %.3f)", string(n), name, worst, maximum(abs, J)))
end
# sensitivity of the front interface to the transparent zone's opacity (TITAN's harmonic mean: 1.0)
for χ in (1e-5, 3.1e-5, 1e-4, 1e-3)
    g(n) = ForwardDiff.derivative(x -> log(χn(21.7, exp(x), n)), log(χ))
    @printf("  front 21.7 | %.1e: chi_Inf/H %.3g  dln chi_Inf/dln chi+ %.3f   chi_3/H %.3g  dln chi_3/dln chi+ %.3f\n", χ, χn(21.7, χ, Inf) / H(21.7, χ), g(Inf), χn(21.7, χ, 3.0) / H(21.7, χ), g(3.0))
end
TX.SLABFLUX[] = false; TX.SFN[] = Inf

println("4. SLABTIE: depth of the donor centre below the interface and the donor factor h")
for (χm, Mm, Mp) in ((21.7, 1e4, 1e4), (2e4, 3e3, 7e3), (3e-5, 2e4, 1e4))
    τi = TX.st_tauint(χm, χm, Mm, Mp)
    check(abs(τi / (χm * Mm) - 1) < 1e-12, @sprintf("equal chi %.2e: tau_i = chi M- (%.1e)", χm, τi / (χm * Mm) - 1))
    τn = TX.st_tauint(χm, χm * (1 + 1e-8), Mm, Mp)
    check(abs(τn / τi - 1) < 1e-7, @sprintf("   series branch |ln R| < 1e-6 continuous: %.1e", τn / τi - 1))
end
# dark-state front of t_sf8 step 5805 (zones 22|23: chi 15.7 | 7.27e-6, cold half-zone depth 0.0042, columns ~ equal)
let χm = 15.67, χp = 7.268e-6, Mp = 0.0042 / 7.268e-6, Mm = Mp
    τi = TX.st_tauint(χm, χp, Mm, Mp); hfs = 1 / (3 * χm * Mm)
    @printf("  dark front: tau_i %.3g (piecewise-constant donor half-depth %.3g), h_FS ~ %.2e -> h %.3f\n", τi, χm * Mm, hfs, TX.st_hdonor(hfs, 1 / 3, χm, χp, Mm, Mp))
    h(x) = TX.st_hdonor(hfs, 1 / 3, exp(x[1]), exp(x[2]), exp(x[3]), exp(x[4]))
    x0 = log.([χm, χp, Mm, Mp]); g = ForwardDiff.gradient(h, x0); worst = 0.0
    for j in 1:4; xp = copy(x0); xm = copy(x0); xp[j] += 1e-6; xm[j] -= 1e-6; worst = max(worst, abs((h(xp) - h(xm)) / 2e-6 - g[j])); end
    check(worst < 1e-8, @sprintf("   AD vs FD of h in ln(chi-, chi+, M-, M+): worst %.1e", worst))
end
let hfs = 0.7, f = 0.8   # no front: the formal solution's h is kept exactly; thin envelope gradient: nearly
    check(TX.st_hdonor(hfs, f, 1e-6, 1e-6, 100.0, 100.0) == hfs, "no front (chi- = chi+): h = h_FS exactly")
    check(abs(TX.st_hdonor(hfs, f, 1.3e-6, 1e-6, 100.0, 100.0) / hfs - 1) < 1e-4, @sprintf("thin envelope (chi ratio 1.3, tau ~ 1e-4): h/h_FS - 1 = %.1e", TX.st_hdonor(hfs, f, 1.3e-6, 1e-6, 100.0, 100.0) / hfs - 1))
end
println(nbad[] == 0 ? "SLABFLUX unit test: all passed" : "SLABFLUX unit test: $(nbad[]) FAILED")
