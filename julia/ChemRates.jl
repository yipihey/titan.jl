"""
    ChemRates

Rate sets for ChemNet (G3/G4). `rhfc02()` reproduces PrimordialChem.rates_rhfc02 (RHFC02 Table A1, no reverses beyond
those listed). `modern()`: AAZN97/GA08 reactions with detailed-balance reverses (`ChemNet.with_db`), H⁻ photodetachment
from John (1988) σ with the Milne inverse, H ground-state photoionisation (hydrogenic σ) with the Milne inverse plus
case-B recombination to excited levels; Li when `li = true` (rates: see rates_li_extraction.md, G4). Sources of every
rate: docs/popiii_rates_extraction.md.
"""
module ChemRates

import ..ChemNet, ..PrimordialChem, ..ContinuumOpacity
const CN = ChemNet; const PC = PrimordialChem; const CO = ContinuumOpacity
R(name, Rs, Ps, k) = CN.Reaction(name, Rs, Ps, k)

"RHFC02 Table A1 exactly as PrimordialChem.rates_rhfc02 (PSS83 k5 or the detailed-balance k5)."
function rhfc02(; k5 = PC.k5_pss)
    rx = [
        R("H+ + e → H (AAZN97 k2, case A)", [:Hp, :e], [:H], (T, env) -> PC.k2_aazn(T)),
        R("H + e → H- (AAZN97 k7)", [:H, :e], [:Hm], (T, env) -> PC.k7_aazn(T)),
        R("H + H- → H2 + e (AAZN97 k8)", [:H, :Hm], [:H2, :e], (T, env) -> PC.k8_aazn(T)),
        R("H + e → H+ + 2e (AAZN97 k1)", [:H, :e], [:Hp, :e, :e], (T, env) -> PC.k1_aazn(T)),
        R("3H → H2 + H (PSS83 k4)", [:H, :H, :H], [:H2, :H], (T, env) -> PC.k4_pss(T)),
        R("H2 + H → 3H (PSS83 k5)", [:H2, :H], [:H, :H, :H], (T, env) -> k5(T)),
        R("2H + H2 → 2H2 (PSS83 k6 = k4/8)", [:H, :H, :H2], [:H2, :H2], (T, env) -> PC.k4_pss(T) / 8),
        R("2H2 → 2H + H2 (PSS83 k7 = k5/8)", [:H2, :H2], [:H, :H, :H2], (T, env) -> k5(T) / 8),
        R("H + H → H+ + e + H (PSS83 k9)", [:H, :H], [:Hp, :e, :H], (T, env) -> PC.k9_pss(T)),
    ]
    CN.Network([:H, :Hp, :Hm, :H2, :H2p, :e, :He], [:Hp, :Hm, :H2, :H2p], rx)
end

"""
Photo-ionisation of atomic H from its levels n ≤ ChemNet.NH_LEVELS (Boltzmann populations, ChemNet.Photo): 1s from
Verner et al. 1996, n ≥ 2 Kramers σ_n = 7.91e-18 n (ν_n/ν)³ cm² (Gaunt factor 1) for ν ≥ ν_n = ν₁/n². Its Milne inverse
is the total radiative recombination (≈ case A); with J = B the pair gives Saha exactly (the partition function of H
includes the same levels). Order-of-magnitude cross sections for n ≥ 2; to refine with exact hydrogenic σ.
"""
levels_H() = [(CN.χH * (1 - 1 / n^2), 2.0n^2, CN.χH / n^2,
               n == 1 ? (ν -> σ_verner(ν, 13.598, VERNER_HI)) : (ν -> 7.91e-18 * n * (CN.χH / n^2 / (6.62607015e-27 * ν))^3))
              for n in 1:CN.NH_LEVELS]
"""
Verner et al. (1996, ApJ 465, 487) outer-shell fit, as coded in D. A. Verner's phfit2.f (version 2, 1996; copy in FAC
`ionis/phfit2.f`): σ = σ₀ [(x−1)² + y_w²] z^{0.5P−5.5} (1 + √(z/y_a))^{−P} Mb, x = E/E₀ − y₀, z = √(x² + y₁²), for E ≥ E_th.
Parameters (E₀ [eV], σ₀ [Mb], y_a, P, y_w, y₀, y₁).
"""
function σ_verner(ν, Eth, p)
    E = 6.62607015e-27 * ν / 1.602176634e-12; E < Eth && return 0.0
    E0, σ0, ya, P, yw, y0, y1 = p
    x = E / E0 - y0; z = sqrt(x^2 + y1^2)
    σ0 * ((x - 1)^2 + yw^2) * z^(0.5P - 5.5) * (1 + sqrt(z / ya))^(-P) * 1e-18
end
const VERNER_HI = (4.298e-1, 5.475e4, 3.288e1, 2.963, 0.0, 0.0, 0.0)       # PH2(1,1): H I 1s
const VERNER_LiI = (3.107, 6.245e1, 1.501e1, 4.895, 0.0, 0.0, 0.0)        # PH2(3,3): Li I 2s, threshold 5.392 eV
"""
Li photo-ionisation from its levels (ChemNet.LI_LEVELS, Boltzmann populations): ground 2s from Verner et al. 1996; excited
levels hydrogenic (Kramers, σ = 7.91e-18 n* (ν_th/ν)³ cm², n* the effective quantum number) — approximate. The Milne
inverse is then the total radiative recombination (≈ case A), checked against GP98 Li1 (Verner & Ferland 1996).
"""
levels_Li() = [(Ex, g, CN.χLi - Ex,
                i == 1 ? (ν -> σ_verner(ν, 5.392, VERNER_LiI)) : (ν -> 7.91e-18 * ns * ((CN.χLi - Ex) / (6.62607015e-27 * ν))^3))
               for (i, (Ex, g, ns)) in enumerate(CN.LI_LEVELS)]
"""
Voronov (1997, ADNDT 65, 1) collisional ionisation by electrons: k = A (1 + P√U) U^K e^{−U}/(X + U) cm³ s⁻¹, U = ΔE/kT
(coefficients as in the AMUSE/SPHRAY and FAC copies of Voronov's table). Li I: ΔE 5.4 eV, P 0, A 1.39e-7, X 0.438, K 0.41.
"""
voronov(T, dE, P, A, X, K) = (U = dE * 11604.518 / T; A * (1 + P * sqrt(U)) * U^K * exp(-U) / (X + U))
"H⁻ photodetachment (John 1988), ν in Hz."
σ_Hm(ν) = CO.σ_Hm_bf(2.99792458e14 / ν)

"""
Modern set: detailed balance for every reversible channel. Returns (net, photos).
"""
function modern(; li = false)
    k(f) = (T, env) -> f(T)
    lst = [
        # collisional ionisation of H by e (Janev, AAZN97 k1) ⇌ three-body recombination (detailed balance)
        (R("H + e → H+ + 2e (AAZN97 k1)", [:H, :e], [:Hp, :e, :e], k(PC.k1_aazn)), true),
        # H⁻ channels
        (R("H + H- → H2 + e (AAZN97 k8)", [:H, :Hm], [:H2, :e], k(PC.k8_aazn)), true),
        (R("H- + H+ → 2H (GA08 k5, Croft et al. 1999)", [:Hm, :Hp], [:H, :H], k(T -> 2.4e-6 / sqrt(T) * (1 + T / 20000))), true),
        (R("H- + e → H + 2e (AAZN97 k14)", [:Hm, :e], [:H, :e, :e], k(T -> exp(sum(c * log(T / 11605)^(i - 1) for (i, c) in enumerate((-18.01849334, 2.3608522, -0.28274430, 1.62331664e-2, -3.36501203e-2, 1.17832978e-2, -1.65619470e-3, 1.06827520e-4, -2.63128581e-6)))))), true),
        # H₂⁺ channels (AAZN97 numbering)
        # radiative association: its inverse is H₂⁺ photodissociation (needs σ(ν) — not yet available): one-way for now
        (R("H + H+ → H2+ + γ (GA08 k3)", [:H, :Hp], [:H2p], k(T -> (l = log10(clamp(T, 1.0, 32000.0)); 10^(-19.38 - 1.523l + 1.118l^2 - 0.1269l^3)))), false),
        (R("H2+ + H → H2 + H+ (AAZN97 k10)", [:H2p, :H], [:H2, :Hp], k(T -> 6.4e-10 + zero(T))), true),
        (R("H2+ + e → 2H (AAZN97 k18)", [:H2p, :e], [:H, :H], k(T -> T < 617 ? 1e-8 + zero(T) : 1.32e-6 * T^-0.76)), true),
        (R("H- + H+ → H2+ + e (GA08 k16)", [:Hm, :Hp], [:H2p, :e], k(T -> T <= 8000 ? 6.9e-9 * T^-0.35 : 9.6e-7 * T^-0.90)), true),
        # three-body H₂ formation (PSS83 k4; H₂ as third body k4/8) ⇌ collisional dissociation (detailed balance)
        (R("3H → H2 + H (PSS83 k4)", [:H, :H, :H], [:H2, :H], k(PC.k4_pss)), true),
        (R("2H + H2 → 2H2 (k4/8)", [:H, :H, :H2], [:H2, :H2], k(T -> PC.k4_pss(T) / 8)), true),
        (R("2H + He → H2 + He (GA08 k32)", [:H, :H, :He], [:H2, :He], k(T -> 6.9e-32 * T^-0.4)), true),
    ]
    rx = CN.with_db(lst)
    photos = [CN.Photo("H- + γ → H + e (John 1988)", :Hm, [:H, :e], σ_Hm, 0.754 * 1.602176634e-12, true),
              CN.Photo("H(n ≤ $(CN.NH_LEVELS)) + γ → H+ + e (Verner 1s + Kramers, Boltzmann levels)", :H, [:Hp, :e], levels_H(), CN.ZH_el, true)]
    species = [:H, :Hp, :Hm, :H2, :H2p, :e, :He]; free = [:Hp, :Hm, :H2, :H2p]
    if li
        # G4: lithium (Li, Li⁺; Li⁻ and LiH not yet). Sources: rates_li_extraction.md (GP98 Table 4), Verner 1996,
        # Voronov 1997. GP98's Li7 (charge exchange Li + H⁺ ⇌ Li⁺ + H, Kimura et al. 1994) is left out: its direction /
        # validity as printed in GP98 is doubtful (flagged in the extraction).
        push!(species, :Li, :Lip); push!(free, :Lip)
        lst_li = [
            (R("Li + e → Li+ + 2e (Voronov 1997)", [:Li, :e], [:Lip, :e, :e], k(T -> voronov(T, 5.4, 0.0, 1.39e-7, 0.438, 0.41))), true),
            (R("Li+ + H- → Li + H (GP98 Li3, Peart & Hayton 1994)", [:Lip, :Hm], [:Li, :H], k(T -> max(6.3e-6 / sqrt(T) - 7.6e-9 + 2.6e-10 * sqrt(T) + 2.7e-14 * T, 0.0))), true),
            (R("Li + H+ → Li+ + H + γ (GP98 Li8, Stancil & Zygelman 1996)", [:Li, :Hp], [:Lip, :H], k(T -> 1.7e-13 * T^-0.051 * exp(-T / 282000))), false),
        ]
        append!(rx, CN.with_db(lst_li))
        push!(photos, CN.Photo("Li(levels) + γ → Li+ + e (Verner 1996 + hydrogenic)", :Li, [:Lip, :e], levels_Li(), CN.ZLi_el, true))
    end
    CN.Network(species, free, rx), photos
end

end # module
