# Tests of H2Lines (M3): thin limit, homogeneous sphere, convergence.  julia -t 8 --project=env h2lines_test.jl
using Printf
include(joinpath(@__DIR__, "PrimordialChem.jl")); include(joinpath(@__DIR__, "H2Lines.jl"))
const PC = PrimordialChem; const HL = H2Lines
L = HL.load()
@printf("RHFC02 subset: %d lines\n", length(L.A))

# 1. thin LTE emission per molecule vs Hollenbach & McKee's LTE fit (RA05 eq. 1.23; HM79 used all lines, ours v ≤ 2, J ≤ 20)
println("\n1. LTE emission per H2 molecule [erg/s]: ExoMol v<=2, J<=20 vs HM79 fit")
for T in (300.0, 1000.0, 2000.0, 3000.0, 5000.0)
    a = HL.thin_per_molecule(L, T); b = PC.Λ_H2_LTE_HM(T)
    @printf("  T = %6.0f K: lines %.4e  HM79 %.4e  ratio %.3f\n", T, a, b, a / b)
end

# 2. thin sphere: transfer result vs n_H2 × thin emission (low n_H2: τ ≪ 1)
function sphere(n; R = 1e16, T0 = 1000.0, nH2 = 1e-3)
    r = collect(range(0.0, R, length = n + 1)); u = zeros(n + 1)
    (r, u, fill(1.0, n), fill(T0, n), fill(nH2, n))
end
println("\n2. thin uniform sphere (R = 1e16 cm, T = 1000 K, n_H2 = 1e-3): L_transfer / (n_H2 V Σ f_u A hν)")
for (n, nb) in ((20, 30), (40, 30), (40, 60))
    r, u, ρ, T, nH2 = sphere(n)
    Λ, Lr, _ = HL.cooling(L, r, u, ρ, T, nH2; nbin = nb)
    Lthin = nH2[1] * 4π / 3 * r[end]^3 * HL.thin_per_molecule(L, T[1])
    @printf("  n = %3d, nbin = %3d: ratio %.6f   (Λ zone 1 / thin %.6f)\n", n, nb, Lr[end] / Lthin, Λ[1] / (nH2[1] * HL.thin_per_molecule(L, T[1])))
end

# 3. thick homogeneous sphere, one frequency bin at a time: emergent flux F = πS [1 − (1 − (1 + 2τ) e^{−2τ})/(2τ²)],
#    τ = radial optical depth. Compare line by line, frequency-integrated with the same bins.
println("\n3. homogeneous sphere, line luminosities vs the analytic emergent flux (same frequency bins)")
for (nH2, n, nb) in ((1e8, 40, 30), (1e8, 80, 30), (1e12, 80, 30), (1e14, 80, 30), (1e14, 160, 30), (1e16, 160, 30), (1e16, 160, 60))
    R = 1e16; T0 = 1000.0
    r, u, ρ, T, nn = sphere(n; R, T0, nH2)
    _, Lr, Ll = HL.cooling(L, r, u, ρ, T, nn; nbin = nb)
    f = HL.populations(L, T0); vth = sqrt(2HL.kB * T0 / HL.mH2); vspan = 3vth
    xs = range(-vspan, vspan, length = nb); dv = step(xs)
    worst = 0.0; τmax = 0.0; Lan_tot = 0.0
    for i in eachindex(L.A)
        ν0 = L.ν[i]; κ0 = HL.h * ν0 / (4π) * HL.c^2 / (2HL.h * ν0^3) * L.A[i] * nH2 * (L.g[L.iu[i]] / L.g[L.il[i]] * f[L.il[i]] - f[L.iu[i]])
        S = HL.Bν(ν0, T0); Lan = 0.0
        for vb in xs
            φ = exp(-(vb / vth)^2) / (sqrt(π) * vth / HL.c * ν0); τ = κ0 * φ * R; τmax = max(τmax, τ)
            Fb = τ < 1e-3 ? π * S * (4τ / 3) : π * S * (1 - (1 - (1 + 2τ) * exp(-2τ)) / (2τ^2))
            Lan += 4π * R^2 * Fb * (dv / HL.c * ν0)
        end
        Lan_tot += Lan
        Lan > 1e-6 * maximum(Ll[:, end]) && (worst = max(worst, abs(Ll[i, end] / Lan - 1)))
    end
    @printf("  n_H2 = %.0e, n = %3d, nbin = %3d: max τ_line = %.2e, total L ratio %.5f, worst line |ΔL/L| %.3e\n",
            nH2, n, nb, τmax, Lr[end] / Lan_tot, worst)
end

# 4. diffusion limit, grey: χ uniform, S = B(T(r)) at one frequency, T = T0 (1 − 0.5 (r/R)²): F → −(4π/3χ) dS/dr when
#    χR ≫ 1 (formal_flux, used for the continuum in popiii_regimes and the same scheme as the lines)
println("\n4. grey diffusion limit: F/F_D at r = 0.3, 0.5, 0.7 R (χR = optical radius)")
for (χR, n) in ((1e1, 80), (1e2, 80), (1e3, 80), (1e3, 160), (1e5, 80), (1e5, 160))
    R = 1e16; χ = χR / R; ν = 1e14; T0 = 2000.0
    r = collect(range(0.0, R, length = n + 1)); rc = 0.5 .* (r[1:end-1] .+ r[2:end])
    S = reshape([HL.Bν(ν, T0 * (1 - 0.5 * (x / R)^2)) for x in rc], n, 1)
    F = HL.formal_flux(r, fill(χ, n, 1), S)
    out = String[]
    for xr in (0.3, 0.5, 0.7)
        k = round(Int, xr * n) + 1; dSdr = (HL.Bν(ν, T0 * (1 - 0.5 * ((r[k] + 1e-6R) / R)^2)) - HL.Bν(ν, T0 * (1 - 0.5 * ((r[k] - 1e-6R) / R)^2))) / 2e-6R
        push!(out, @sprintf("%.4f", F[k] / (-(4π / (3χ)) * dSdr)))
    end
    @printf("  chi R = %.0e (Δτ per zone %.1e), n = %3d: F/F_D = %s\n", χR, χR / n, n, join(out, ", "))
end

# 5. per-line frequency grids (G6a): hot, mostly dissociated core (2e4 K, 1 % of the H₂ density) inside a cold molecular
#    envelope (200–2000 K) falling in at up to 20 km/s. Reference: all zones in every line's span, up to 4000 bins.
println("\n5. per-line frequency grids vs a fully resolved global grid (hot core + cold infalling envelope)")
let n = 80, R = 1e15
    r = collect(range(0.0, R, length = n + 1)); rc = 0.5 .* (r[1:end-1] .+ r[2:end])
    T = [x < 0.1R ? 2e4 : 2000 * (0.1R / x) for x in rc]
    nH2 = [x < 0.1R ? 1e10 : 1e12 * (0.1R / x)^2 for x in rc]
    u = [-2e6 * min(1.0, x / (0.2R)) for x in r]
    t0 = time(); Λa, La, _ = HL.cooling(L, r, u, ones(n), T, nH2); ta = time() - t0; nba = HL.LAST[].nbins
    t0 = time(); Λb, Lb, _ = HL.cooling(L, r, u, ones(n), T, nH2; cover = 1.0, nbinmax = 4000, nbin = 4000); tb = time() - t0
    @printf("  L(surface): per-line %.5e (≤ %d bins, %.1f s), reference %.5e (4000 bins, %.1f s): ratio %.5f; max |ΔΛ/Λ| over zones with ≥ 1e-3 of the peak cooling: %.3e\n",
            La[end], nba, ta, Lb[end], tb, La[end] / Lb[end],
            maximum(abs(Λa[k] / Λb[k] - 1) for k in 1:n if abs(Λb[k]) > 1e-3 * maximum(abs, Λb)))
end
