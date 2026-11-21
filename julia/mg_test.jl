# MGROUPS unit test (TitanResidFast MGRad.jl, G2Opacity.group_fast; docs/POPIII_MULTIGROUP.md):
#   1. Planck group fractions b_g(T): sum to 1, against direct quadrature of B_ν over each group, and AD vs central differences in T;
#   2. group means of the G2 continuum on a set of states (H⁻-dominated photosphere gas, cool molecular gas, warm ionised gas): the
#      identities Σ_g b_g(T) κ_P,g = κ_P, Σ_g b_g(T_w) κ_E,g = κ_E, Σ_g b'_g(T)/χ_g = 1/χ (grey means of the same call), and the
#      ionising group's H(1s) part.
#   TITAN_NOLIB=1 julia --project=<env> mg_test.jl
using Printf, ForwardDiff
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "Z0Physics.jl"))
include(joinpath(@__DIR__, "ContinuumOpacity.jl")); include(joinpath(@__DIR__, "G2Opacity.jl"))
const TX = TitanResidFast; const G2O = G2Opacity; const COp = ContinuumOpacity; const Z0 = Z0Physics
nbad = Ref(0)
check(ok, msg) = (ok || (nbad[] += 1); println(ok ? "  ok   " : "  FAIL ", msg))
const EV = 1.602176634e-12; const KB = 1.380649e-16; const H = 6.62607015e-27; const C = 2.99792458e10
eHm = H * C * (1e4 / 1.6419) / EV; eH = Z0.χH / EV
for edges in ([eHm, 11.2, eH], [eHm, eH / 9, eH / 4, 11.2, eH, Z0.χHe1 / EV, Z0.χHe2 / EV])
    TX.mg_setup!(edges, EV / KB, 54.5, 9, 10)
    G = length(edges) + 1
    println("1. Planck fractions, $G groups, edges ", join((@sprintf("%.4f", e) for e in edges), ", "), " eV")
    Bν(ν, T) = (x = H * ν / (KB * T); x > 700 ? 0.0 : 2H * ν^3 / C^2 / expm1(x))
    worst = 0.0; wsum = 0.0; wad = 0.0; wloc = (0.0, 0)
    for T in (30.0, 54.5, 300.0, 1e3, 3e3, 5300.0, 1e4, 3e4, 1e5, 3e5)
        b = [TX.mg_bfrac(g, T) for g in 1:G]
        bs = sum(b[g] for g in 1:G if b[g] > TX.MGBFLOOR); wsum = max(wsum, abs(bs - 1))
        lo = [0.0; edges .* (EV / H)]; hi = [edges .* (EV / H); Inf]
        for g in 1:G                                   # direct quadrature in ln ν (5000 points per group, Simpson)
            a = log(max(lo[g], 1e-6 * KB * T / H)); z = log(min(hi[g], 200 * KB * T / H)); z > a || continue
            n = 5000; h = (z - a) / n
            s = sum((i == 0 || i == n ? 1 : iseven(i) ? 2 : 4) * (ν = exp(a + i * h); Bν(ν, T) * ν) for i in 0:n) * h / 3
            q = s / (5.670374419e-5 * T^4 / π)
            q > 1e-150 && (worst = max(worst, abs(b[g] / q - 1)))
        end
        for g in 1:G                                   # AD vs central difference in T
            # (4th-order stencil, h = 1e-3 T / max(1, x) with x the group's largest finite edge in kT: d ln b/d ln T ~ x, so a fixed relative step
            #  is truncation-limited at large x; a derivative below 1e-12 b/T (b ≈ 1 of the IR group in cold gas) is round-off in the difference)
            xm = maximum(e * (EV / KB) / T for e in edges[max(g-1, 1):min(g, G-1)])
            d = ForwardDiff.derivative(t -> TX.mg_bfrac(g, t), T); hh = 1e-3T / max(1.0, xm); f(t) = TX.mg_bfrac(g, t)
            fd = (8 * (f(T + hh) - f(T - hh)) - (f(T + 2hh) - f(T - 2hh))) / 12hh
            abs(d) * T < 1e-12 * b[g] && continue
            if b[g] > 0.5               # (the dominant group: its difference quotient is round-off of 1 − b; Σ_g b_g = 1 instead)
                fd = -sum(ForwardDiff.derivative(t -> TX.mg_bfrac(g2, t), T) for g2 in 1:G if g2 != g)
            end
            b[g] > 1e-100 && abs(fd) > 1e-300 && abs(d / fd - 1) > wad && (wad = abs(d / fd - 1); wloc = (T, g))
        end
    end
    check(wsum < 1e-12, @sprintf("Σ_g b_g = 1 (30 K … 3e5 K): worst %.1e", wsum))
    check(worst < 1e-7, @sprintf("b_g against direct quadrature (b_g > 1e-150): worst relative difference %.1e", worst))
    check(wad < 1e-6, @sprintf("d b_g/dT: AD vs central difference worst %.1e (T %g, group %d)", wad, wloc...))
end

println("2. group means of the G2 continuum (4 groups) and the grey identities")
edges = [eHm, 11.2, eH]; G = 4; TX.mg_setup!(edges, EV / KB, 54.5, 9, 10)
gidx = G2O.group_index(edges .* (EV / H))
dbfrac(g, T) = ForwardDiff.derivative(t -> TX.mg_bfrac(g, t) * t^4, T) / (4T^3)        # fraction of ∂B/∂T in group g
states = [  # (label, ρ, T, Tw, nH fractions: x_H, x_e(=x_p), x_H2, x_H-)
    ("photosphere gas 5000 K", 1e-9, 5000.0, 5300.0, 0.999, 1e-4, 1e-4, 1e-8),
    ("molecular 1500 K", 1e-11, 1500.0, 5300.0, 0.2, 1e-7, 0.4, 1e-12),
    ("warm ionised 1.2e4 K", 1e-10, 1.2e4, 1.2e4, 0.05, 0.95, 0.0, 1e-12),
    ("neutral 8000 K", 1e-8, 8000.0, 8000.0, 0.99, 0.01, 0.0, 1e-8)]
for (lab, ρ, T, Tw, xH, xe, xH2, xHm) in states
    nH = ρ / (1.4 * 1.6735575e-24); n1 = xH * nH; ne = xe * nH; np = ne; nH2 = xH2 * nH; nHm = xHm * nH; nHe = 0.0789 * nH
    n2 = G2O.n2_lte(T, n1, np, n1 + nHe + nH2)
    κE, κP, χ, k1, gE, gP, gχ = G2O.group_fast(ρ, T, gidx, G; n1, n2, nH2, nHe, nHp = np, ne, nHm, hbfP = false, Tw, Tmin = 1.0)
    sP = sum(TX.mg_bfrac(g, T) * κP[g] for g in 1:G) / gP; sE = sum(TX.mg_bfrac(g, Tw) * κE[g] for g in 1:G) / gE
    sR = sum(dbfrac(g, T) / χ[g] for g in 1:G if χ[g] > 0) * gχ
    check(max(abs(sP - 1), abs(sE - 1), abs(sR - 1)) < 2e-2, @sprintf("%-24s Σ b κP/κP %.4f, Σ b(Tw) κE/κE %.4f, Σ b'/χ·χ %.4f | rP %s | rR %s | κ1s ion %.3e (κE grey %.3e)",
          lab, sP, sE, sR, join((@sprintf("%.3g", κP[g] / gP) for g in 1:G), " "), join((@sprintf("%.3g", χ[g] / gχ) for g in 1:G), " "), k1[4], gE))
    # the fast grey means of the same physics
    m = G2O.means_fast(ρ, T; n1, n2, nH2, nHe, nHp = np, ne, nHm, hbfP = false, Tw)
    check(maximum(abs.([gE, gP, gχ] ./ collect(m) .- 1)) < 1e-12, @sprintf("%-24s grey means of group_fast = means_fast (%.1e)", lab, maximum(abs.([gE, gP, gχ] ./ collect(m) .- 1))))
end
println(nbad[] == 0 ? "all passed" : "$(nbad[]) FAILED")
