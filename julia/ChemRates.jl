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
Photo-ionisation cross section of atomic H with its levels n ≤ ChemNet.NH_LEVELS populated in Boltzmann ratio to the
ground state at T: σ(ν, T) = Σ_n p_n(T) σ_n(ν), Kramers σ_n = 7.91e-18 n (ν_n/ν)³ cm² (g = 1) for ν ≥ ν_n = ν₁/n². Its Milne
inverse is the total radiative recombination (≈ case A); with J = B the pair gives Saha exactly (the partition function
of H includes the same levels). Order-of-magnitude cross sections (Gaunt factors 1); to refine with exact hydrogenic σ.
"""
function σ_H(ν, T)
    ν1 = 13.598 * 1.602176634e-12 / 6.62607015e-27; Z = CN.ZH_el(T); s = 0.0
    for n in 1:CN.NH_LEVELS
        νn = ν1 / n^2; ν < νn && continue
        pn = 2n^2 * exp(-13.598 * 1.602176634e-12 * (1 - 1 / n^2) / (1.380649e-16 * T)) / Z
        s += pn * 7.91e-18 * n * (νn / ν)^3
    end
    s
end
"H⁻ photodetachment (John 1988), ν in Hz."
σ_Hm(ν, T) = CO.σ_Hm_bf(2.99792458e14 / ν)

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
              CN.Photo("H(n ≤ 8) + γ → H+ + e (hydrogenic, Boltzmann levels)", :H, [:Hp, :e], σ_H, 13.598 * 1.602176634e-12 / CN.NH_LEVELS^2, true)]
    species = [:H, :Hp, :Hm, :H2, :H2p, :e, :He]; free = [:Hp, :Hm, :H2, :H2p]
    if li
        push!(species, :Li, :Lip); push!(free, :Lip)
        # (G4: Li reactions from rates_li_extraction.md are added here)
    end
    CN.Network(species, free, rx), photos
end

end # module
