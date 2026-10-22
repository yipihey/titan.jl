# continuum_ff_test.jl — checks of ContinuumFF.jl and its wiring into ContinuumOpacity / Z0Physics:
# table nodes reproduced, Gaunt factors in range, every new absorption ≥ 0 and continuous at the table
# edges, and the effect on the Planck/Rosseland means of the dense primordial gas (EXT_FF on/off).
#   julia --project=env continuum_ff_test.jl
using Printf
include(joinpath(@__DIR__, "ContinuumOpacity.jl")); const CO = ContinuumOpacity; const FF = CO.FF
include(joinpath(@__DIR__, "Z0Physics.jl")); const Z0 = Z0Physics
pass = Ref(0); fail = Ref(0)
chk(name, ok) = (ok ? (pass[] += 1) : (fail[] += 1); @printf("%-78s %s\n", name, ok ? "PASS" : "FAIL"))
const c = 2.99792458e10

# ---- table nodes
chk("H₂⁻ ff node (1.0003 μm, 5000 K) = 0.380e-26", isapprox(FF.k_H2m_ff(10002.74, 5000.0), 0.380e-26; rtol = 1e-12))
chk("H₂⁻ ff node (2.5007 μm, 1000 K) = 3.92e-26", isapprox(FF.k_H2m_ff(25006.82, 1000.0), 3.92e-26; rtol = 1e-12))
chk("He⁻ ff node (λ = 0.7594 μm, θ = 1) = 0.103e-26", isapprox(FF.k_Hem_ff(7594.0, 5040.0), 0.103e-26; rtol = 1e-12))
chk("He⁻ ff node (λ = 1.51878 μm, θ = 0.5) = 27.979e-26", isapprox(FF.k_Hem_ff(151878.0, 10080.0), 27.979e-26; rtol = 1e-12))
chk("H₂⁺ K(8400 K) = 2.6532e21", isapprox(FF.K_H2p(8400.0), 2.6532e21; rtol = 1e-12))
bf, ff = FF.α_H2p(c / 1e-4, 12600.0, 1.0, 1.0)   # λ = 1 μm (a node of both λ grids), T node 12600 K (column 6)
stim = -expm1(-6.62607015e-27 * c / 1e-4 / (1.380649e-16 * 12600.0))
ibf = findfirst(==(10000.0), FF.H2P_LAM_BF); iff = findfirst(==(10000.0), FF.H2P_LAM_FF)
chk("H₂⁺ bf node (1 μm, 12600 K) reproduces σ_bf(table) × n_H n_H⁺/K × stim", isapprox(bf, FF.H2P_BF[ibf, 6] / FF.K_H2p(12600.0) * stim; rtol = 1e-12))
chk("H₂⁺ ff node (1 μm, 12600 K) reproduces the table entry × stim", isapprox(ff, FF.H2P_FF[iff, 6] * stim; rtol = 1e-12))
chk("H₂⁺ bf/ff magnitudes at that node (Stancil: σ_bf ~1e-18 cm², ff ~1e-39 cm⁵)", 1e-19 < FF.H2P_BF[ibf, 6] < 1e-17 && 1e-40 < FF.H2P_FF[iff, 6] < 1e-38)
chk("K_H2p Arrhenius below the table: K(2000 K) < K(3150 K)/10", FF.K_H2p(2000.0) < FF.K_H2p(3150.0) / 10)
# ---- Gaunt factors: within the known range, 1 at the classical point, Z-scaling exact (γ² depends on Z²/T only)
gs = [FF.gaunt_ff(ν, T, 1) for ν in (1e12, 1e13, 1e14, 1e15, 3e15), T in (3e3, 1e4, 3e4, 1e5, 1e6)]
chk("ḡ_ff(Z = 1) within [1.0, 6] over 1e12–3e15 Hz, 3e3–1e6 K (van Hoof table range)", all(1.0 .<= gs .<= 6.0))
chk("ḡ_ff(ν, T, Z = 2) = ḡ_ff(ν, T/4, Z = 1) with u scaled: table coordinates consistent",
    isapprox(FF.gaunt_ff(1e15, 1e5, 2), FF.gaunt_ff(1e15 / 4, 1e5 / 4, 1); rtol = 1e-12))
chk("ḡ_ff(1e4 K, 1 μm, Z = 1) ≈ 1.16 (Karzas–Latter regime 1.1–1.3)", 1.05 < FF.gaunt_ff(c / 1e-4, 1e4, 1) < 1.35)
# ---- relativistic merged tables, integrated Gaunt factor
chk("merged (2015) and non-relativistic (2014) Gaunt factors agree < 0.5 % at ≤ 1e7 K",
    all(isapprox(FF.gaunt_ff(ν, T, Z), FF.gaunt_ff(ν, T, Z; relativistic = false); rtol = 5e-3) for ν in (1e13, 1e15, 1e17), T in (1e4, 1e6, 1e7), Z in (1, 2)))
chk("merged table differs from non-relativistic at 1e9 K (relativistic regime)",
    !isapprox(FF.gaunt_ff(1e18, 1e9, 1), FF.gaunt_ff(1e18, 1e9, 1; relativistic = false); rtol = 1e-2))
cen(T) = 1.1 + 0.34 * exp(-(5.5 - log10(T))^2 / 3)
chk("⟨g_ff⟩ (van Hoof) vs Cen 1992 integrated fit ≤ 2 % for Z = 1, 2 over 3e3–3e8 K",
    all(isapprox(FF.gaunt_ff_total(T, Z), cen(T / Z^2); rtol = 2e-2) for T in 10 .^ (3.5:0.25:8.5), Z in (1, 2)))
# ---- continuity at the clamps / extrapolations
chk("He⁻ ff continuous at 10 080 K and positive to 3e4 K", isapprox(FF.k_Hem_ff(1e4, 10080.0), FF.k_Hem_ff(1e4, 10081.0); rtol = 1e-3) &&
    all(FF.k_Hem_ff(λ, T) > 0 for λ in (6e3, 1e4, 3e4), T in (12000.0, 20000.0, 30000.0)))
chk("He⁻ ff continuous across the red table edge", isapprox(FF.k_Hem_ff(151878.0 * (1 + 1e-9), 5040.0), FF.k_Hem_ff(151878.0, 5040.0); rtol = 1e-6))
chk("H₂⁻ ff continuous across the T clamp (15 000 K)", isapprox(FF.k_H2m_ff(1e4, 15000.0), FF.k_H2m_ff(1e4, 15001.0); rtol = 1e-6))
chk("H₂⁺ continuous across 2 μm", all(isapprox.(FF.α_H2p(c / 2.0e-4, 8400.0, 1.0, 1.0), FF.α_H2p(c / 2.0000001e-4, 8400.0, 1.0, 1.0); rtol = 1e-5)))
# ---- alpha: all parts ≥ 0 at representative dense states; means with/without the new processes
function state(nH, T, fH2, xHp, xe)
    nH2 = fH2 * nH / 2; nHI = nH * (1 - fH2 - xHp); nHe = nH / 12; ne = xe * nH
    # H⁻ in Saha with the local electrons: n_H⁻ = n_H n_e e^{0.754 eV/kT}/(4 λ_e⁻³)
    kB = 1.380649e-16; λe3 = (6.62607015e-27^2 / (2π * 9.1093837e-28 * kB * T))^1.5
    nHm = nHI * ne * exp(0.754 * 1.602176634e-12 / (kB * T)) * λe3 / 4
    (nH = nHI, nH2 = nH2, nHe = nHe, nHp = xHp * nH, ne = ne, nHm = nHm)
end
println("\nPlanck / Rosseland means [cm² g⁻¹], EXT_FF off → on (n_H, T, f_H2, x_H⁺, x_e; He neutral unless given):")
cases = ((1e14, 2000.0, 0.9, 1e-12, 4e-10), (1e16, 3000.0, 0.5, 1e-10, 4e-10), (1e18, 3500.0, 0.3, 1e-8, 5e-10),
         (1e20, 8000.0, 0.05, 2e-4, 2e-4), (1e21, 2.0e4, 0.0, 0.3, 0.3))
okpos = true
for (nH, T, fH2, xHp, xe) in cases
    global okpos
    s = state(nH, T, fH2, xHp, xe); ρ = nH * 1.6735575e-24 * (1 + 4 / 12)
    for ν in (3e13, 1e14, 3e14, 1e15)
        a, sc, p = CO.alpha(ν, T; s...)
        okpos &= all(v >= 0 for v in values(p)) && a >= 0
    end
    CO.EXT_FF[] = false; m0 = CO.means(ρ, T; s...)
    CO.EXT_FF[] = true;  m1 = CO.means(ρ, T; s...)
    @printf("  %7.0e %7.0f %4.2f %7.0e %7.0e : κ_P %.3e → %.3e (×%.2f)   κ_R %.3e → %.3e (×%.2f)\n", nH, T, fH2, xHp, xe,
            m0.κP, m1.κP, m1.κP / m0.κP, m0.κR, m1.κR, m1.κR / m0.κR)
end
chk("all absorption parts ≥ 0 at the dense states", okpos)
# hot ionised gas with He ions: He⁺/He²⁺ ff present and the Z-scaling visible
s = state(1e20, 3.0e4, 0.0, 0.9, 0.9)
a0, _, p0 = CO.alpha(1e14, 3.0e4; s...); a1, _, p1 = CO.alpha(1e14, 3.0e4; s..., nHep = 0.5e20 / 12, nHepp = 0.3e20 / 12)
chk("He⁺/He²⁺ free–free adds absorption in ionised gas", p1.he_ff > 0 && a1 > a0)
# ---- Z0Physics Gaunt switch: ratio of the ff term in the expected range, default off unchanged
ρz, Tz = 1e-8, 3.0e4
sz = Z0.saha(ρz, Tz); pH, pHe = Z0.level_pops(Tz, sz)
Z0.GAUNT_FF[] = false; k0 = Z0.kappa_nu(1e14, Tz, sz, pH, pHe)[1]
Z0.GAUNT_FF[] = true;  k1 = Z0.kappa_nu(1e14, Tz, sz, pH, pHe)[1]; Z0.GAUNT_FF[] = false
@printf("Z0 kappa_nu(1e14 Hz, 3e4 K, 1e-8 g/cc): g = 1 %.4e, van Hoof %.4e (×%.3f)\n", k0, k1, k1 / k0)
chk("Z0 Gaunt switch raises the ff-dominated κ_ν by 1.0–2.5×", 1.0 < k1 / k0 < 2.5)
@printf("\n%d passed, %d failed\n", pass[], fail[])
