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
