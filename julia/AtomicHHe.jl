"""
    AtomicHHe

S6c: atomic H and He chemistry and cooling for ChemNet (opt-in, `AtomicHHe.atomic()`). Species He, He⁺, He²⁺ (term-resolved He I and hydrogenic
He II partition functions in ChemNet: Keq gives Saha exactly); collisional ionisation by electrons (Voronov 1997 / AAZN97-Janev) with detailed-balance
three-body recombination; radiative recombination of H⁺, He⁺, He²⁺ through Milne photo pairs (level lists: Verner et al. 1996 ground states,
Johnson 1972 hydrogenic / Kramers excited levels); dielectronic recombination of He⁺ (Aldrovandi & Péquignot 1973, one-way); charge exchange
He⁺ + H ⇌ He + H⁺. Cooling (`cooling`): net radiative loss of every photo pair from the same integrals (recombination photons minus absorbed
J_ν), collisional excitation of H, He I, He II (multi-level), dielectronic recombination, free–free, Compton; escape probabilities of the lines
and on-the-spot absorption of the Lyman-type recombination continua. Sources: docs/POPIII_NOTES.md S6c.
"""
module AtomicHHe

import ..ChemNet, ..ChemRates, ..PrimordialChem
using ForwardDiff
const CN = ChemNet; const CR = ChemRates; const PC = PrimordialChem
const kB = CN.kB; const hP = CN.hP; const c = CN.c; const eV = CN.eV; const me = CN.me
const NL = CN.NH_LEVELS

# ------------------------------------------------------------------ photo-ionisation cross sections and level lists
"Johnson (1972) bound–free Gaunt coefficients (g0, g1, g2) for level n ≥ 2 (as NLTEH)."
function gcoef(n)
    n == 2 && return (1.0785, -0.2319, 0.02947)
    (0.9935 + 0.2328 / n - 0.1296 / n^2, -(0.6282 - 0.5598 / n + 0.5299 / n^2) / n, (0.3887 - 1.181 / n + 1.470 / n^2) / n^2)
end
"Hydrogenic photo-ionisation cross section of level n of a nucleus of charge Z at frequency ν [cm²] (Johnson 1972 × Z⁻²σ_H(ν/Z²)); 0 below threshold."
function σ_hyd(n, Z, ν)
    ν0 = Z^2 * CN.χH / n^2 / hP; ν < ν0 && return 0.0
    x = ν / ν0; g0, g1, g2 = gcoef(n); νH = ν / Z^2
    2.815e29 / (n^5 * νH^3) * (g0 + g1 / x + g2 / x^2) / Z^2
end
const VERNER_HeI = (1.361e1, 9.492e2, 1.469, 3.188, 2.039, 4.434e-1, 2.136)     # PH2(2,2) of phfit2.f: He I 1s², threshold 24.59 eV
const VERNER_HeII = (1.720, 1.369e4, 3.288e1, 2.963, 0.0, 0.0, 0.0)             # PH2(2,1): He II 1s, threshold 54.42 eV
"Scale of the Kramers cross sections of the excited He I levels relative to 7.91e-18 n*(ν_th/ν)³ (calibrated to Verner & Ferland 1996 case A, see the test)."
const HE_EXC_SCALE = Ref(1.15)

"H levels as ChemRates.levels_H but with the Johnson (1972) Gaunt factors for n ≥ 2 (the 1s from Verner et al. 1996)."
levels_H_johnson() = [(CN.χH * (1 - 1 / n^2), 2.0n^2, CN.χH / n^2,
                       n == 1 ? (ν -> CR.σ_verner(ν, 13.598, CR.VERNER_HI)) : (ν -> σ_hyd(n, 1, ν))) for n in 1:NL]
"He II levels: hydrogenic Z = 2 (ground: Verner et al. 1996, n ≥ 2: Johnson 1972)."
levels_HeII() = [(CN.χHe2 * (1 - 1 / n^2), 2.0n^2, CN.χHe2 / n^2,
                  n == 1 ? (ν -> CR.σ_verner(ν, 54.418, VERNER_HeII)) : (ν -> σ_hyd(n, 2, ν))) for n in 1:NL]
"He I levels (CN.HE_LEVELS): ground 1 ¹S Verner et al. 1996; excited Kramers σ = scale 7.91e-18 n*(ν_th/ν)³, n* = (Ry/E_th)^½."
levels_HeI() = [(Ex, g, CN.χHe1 - Ex,
                 i == 1 ? (ν -> CR.σ_verner(ν, 24.587, VERNER_HeI)) :
                 (Eth = CN.χHe1 - Ex; ns = sqrt(CN.Ry_eV * eV / Eth); ν -> ν < Eth / hP ? 0.0 : HE_EXC_SCALE[] * 7.91e-18 * ns * (Eth / (hP * ν))^3))
                for (i, (Ex, g, n)) in enumerate(CN.HE_LEVELS)]

# ------------------------------------------------------------------ reaction rates
const KEV = 11604.518
"Janev et al. (1987) fits as AAZN97 Table 1 (3), (5) / GA08 k17, k18, for comparison with Voronov; T in K."
function janev_He0(T)
    l = log(T / KEV); co = (-44.09864886, 23.91596563, -10.7532302, 3.05803875, -0.56851189, 6.79539123e-2, -5.00905610e-3, 2.06723616e-4, -3.64916141e-6)
    T < 4000 ? zero(T) : exp(sum(a * l^(i - 1) for (i, a) in enumerate(co)))
end
function janev_Hep(T)
    l = log(T / KEV); co = (-68.71040990, 43.93347633, -18.4806699, 4.70162649, -0.76924663, 8.113042e-2, -5.32402063e-3, 1.97570531e-4, -3.16558106e-6)
    exp(sum(a * l^(i - 1) for (i, a) in enumerate(co)))
end
"Voronov (1997) He I (24.6 eV; A 1.75e-8, X 0.180, K 0.35) and He II (54.4 eV; P 1, A 2.05e-9, X 0.265, K 0.25)."
voronov_He0(T) = CR.voronov(T, 24.6, 0.0, 1.75e-8, 0.180, 0.35)
voronov_Hep(T) = CR.voronov(T, 54.4, 1.0, 2.05e-9, 0.265, 0.25)
"Voronov (1997) H I (13.6 eV; A 2.91e-8, X 0.232, K 0.39), for comparison with the Janev fit used by ChemRates."
voronov_H(T) = CR.voronov(T, 13.6, 0.0, 2.91e-8, 0.232, 0.39)
"Dielectronic recombination He⁺ + e → He, Aldrovandi & Péquignot (1973) = Black (1981) = GA08 k19,di = AAZN97 k4d [cm³ s⁻¹]."
k_di_HeII(T) = 1.9e-3 * T^-1.5 * exp(-473421 / T) * (1 + 0.3 * exp(-94684 / T))
"Charge exchange He⁺ + H → He + H⁺ (GA08 k26; Kingdon & Ferland 1996 form), and GA08's He + H⁺ → He⁺ + H (k27)."
k26_CX(T) = 1.2e-15 * (T / 300)^0.25
k27_CX_GA08(T) = T <= 1e4 ? 1.26e-9 * T^-0.75 * exp(-127500 / T) : 4.0e-37 * T^4.74
"Radiative recombination fits for comparison: Verner & Ferland (1996) total (case A) a/[√(T/T0)(1+√(T/T0))^{1−b}(1+√(T/T1))^{1+b}]."
vf96(a, b, T0, T1, T) = a / (sqrt(T / T0) * (1 + sqrt(T / T0))^(1 - b) * (1 + sqrt(T / T1))^(1 + b))
vf96_H(T) = vf96(7.982e-11, 0.7480, 3.148, 7.036e5, T)
vf96_HeI(T) = vf96(3.294e-11, 0.6910, 1.554e1, 3.676e7, T)
vf96_HeII(T) = vf96(1.891e-10, 0.7524, 9.370, 2.774e6, T)

"""
Forward fit continued below its validity limit Tx as k(Tx) (T/Tx)^½ exp[−(ΔE/k)(1/T − 1/Tx)] (Arrhenius with the threshold energy ΔE/k [K]): physically
irrelevant there (≲ e^{−13} of the value at Tx) but it keeps the detailed-balance reverse k_f/K exact, which `ChemNet.reverse_k` would
otherwise extrapolate (error ≈ 1e-5 of Saha at n_H = 1e20, T = 3e3 K).
"""
arrh_ext(kf, Tx, dEK) = T -> T >= Tx ? kf(T) : kf(Tx) * sqrt(T / Tx) * exp(-dEK * (1 / T - 1 / Tx))

# ------------------------------------------------------------------ on-the-spot absorption (escape of the recombination continua) and line escape
"""
Zone state for the escape probabilities (set by the driver per zone before the chemistry step, `set_zone!`): radius r [cm] (0 = optically thin, the default;
r is the length of the continuum path), Sobolev length L_v = 1/max(|du/dr|, |u|/r) [s] (velocity-gradient scale; Inf = static), H-nuclei density, composition at the
start of the step (neutral fractions of H, He, He⁺).
Line escape length of an absorber with Doppler parameter b: ℓ = min(√π b L_v, r), i.e. the Sobolev optical depth τ = (πe²/mc) f λ n_l L_v capped by the static column σ₀ n r.
"""
mutable struct ZoneState
    r::Float64; Lv::Float64; nH::Float64; xH::Float64; xHe::Float64; xHep::Float64
end
const ZS = ZoneState(0.0, Inf, 1.0, 0.0, 0.0, 0.0)
const ONSPOT = Ref(true)
function set_zone!(r, Lv, nH, d)
    ZS.r = r; ZS.Lv = Lv; ZS.nH = nH
    ZS.xH = get(d, :H, 0.0); ZS.xHe = get(d, :He, 0.0); ZS.xHep = get(d, :Hep, 0.0)
    nothing
end
βesc(τ) = τ < 1e-8 ? one(τ) - τ / 2 : -expm1(-τ) / τ
σ_abs(E_eV) = begin                    # photo-absorption cross sections [cm²] at E: (H I, He I, He II)
    ν = E_eV * eV / hP
    (CR.σ_verner(ν, 13.598, CR.VERNER_HI), CR.σ_verner(ν, 24.587, VERNER_HeI), CR.σ_verner(ν, 54.418, VERNER_HeII))
end
"Continuum optical depth over r at photon energy E [eV] from the zone state."
function τ_cont(E_eV)
    sH, sHe, sHep = σ_abs(E_eV)
    ZS.r * ZS.nH * (sH * ZS.xH + sHe * ZS.xHe + sHep * ZS.xHep)
end
"Escape fraction of the recombination continuum to the ground level of species R (H, He, Hep): 1 when thin, (1 − e^{−τ})/τ otherwise (on-the-spot: case B for τ ≫ 1)."
function esc_ground(R::Symbol)
    Eth = R === :H ? 13.598 : R === :He ? 24.587 : 54.418
    (T, env) -> (ONSPOT[] && ZS.r > 0) ? βesc(τ_cont(Eth)) : 1.0
end

"Level lists split into the ground level (continuum subject to on-the-spot absorption) and the excited levels."
function photo_pair(name, R, ion, levels, Zf)
    g = CN.Photo(name * " ground", R, [ion, :e], levels[1:1], Zf, true, esc_ground(R))
    e = CN.Photo(name * " excited", R, [ion, :e], levels[2:end], Zf, true, (T, env) -> 1.0)
    [g, e]
end

# ------------------------------------------------------------------ the network
"""
    atomic(; li = false, cx = :db, ionfit = :voronov) -> (net, photos)

ChemRates.modern plus He, He⁺, He²⁺: Voronov (or Janev, `ionfit = :janev`) collisional ionisation of He and He⁺ with detailed-balance three-body
recombination, charge exchange He⁺ + H ⇌ He + H⁺ (`cx = :db`: GA08 k26 with the reverse by detailed balance; `:ga08`: k26 and the independent GA08 k27;
`:none`), dielectronic recombination of He⁺ (one-way: the detailed-balance reverse k/K would be a unimolecular thermal ionisation of He of 2e-9 s⁻¹ at 1.5e4 K — it only exists in a thermal photon bath — so it is added only with `jb_inverse` or `dr_db`), and the photo pairs of H, He I, He II (levels above, ground levels subject to the on-the-spot
escape fraction). `jb_inverse` adds the inverses of the one-way radiative processes (H₂⁺ association, Li + H⁺ → Li⁺ + H + γ, dielectronic recombination, mutual neutralisation) by detailed balance — only meaningful for J = B (test a). H⁻, H₂, H₂⁺ and (li) Li as in `modern`. HeH⁺ is not included (see the notes).
"""
function atomic(; li = false, cx = :db, ionfit = :voronov, dr_db = false, jb_inverse = false)
    net0, ph0 = CR.modern(; li)
    k(f) = (T, env) -> f(T)
    kHe0 = ionfit === :janev ? janev_He0 : voronov_He0; kHep = ionfit === :janev ? janev_Hep : voronov_Hep
    # H ionisation: the Janev fit of modern() continued below 1e4 K (Arrhenius), replacing its reaction and reverse
    filter!(r -> !startswith(r.name, "H + e → H+ + 2e") && !startswith(r.name, "H+ + 2e → H + e"), net0.reactions)
    kH = arrh_ext(PC.k1_aazn, 1e4, 157807.0)
    lst = [
        (CR.R("H + e → H+ + 2e (AAZN97 k1 = Janev, Arrhenius below 1e4 K)", [:H, :e], [:Hp, :e, :e], k(kH)), true),
        (CR.R("He + e → He+ + 2e (Voronov 1997)", [:He, :e], [:Hep, :e, :e], k(arrh_ext(kHe0, 5e3, 285335.0))), true),
        (CR.R("He+ + e → He2+ + 2e (Voronov 1997)", [:Hep, :e], [:Hepp, :e, :e], k(arrh_ext(kHep, 1e4, 631515.0))), true),
        (CR.R("He+ + e → He + γ (dielectronic, Aldrovandi & Péquignot 1973)", [:Hep, :e], [:He], k(k_di_HeII)), dr_db || jb_inverse),
    ]
    cx === :none || push!(lst, (CR.R("He+ + H → He + H+ (GA08 k26)", [:Hep, :H], [:He, :Hp], k(k26_CX)), cx === :db))
    cx === :ga08 && push!(lst, (CR.R("He + H+ → He+ + H (GA08 k27)", [:He, :Hp], [:Hep, :H], k(k27_CX_GA08)), false))
    rx = vcat(net0.reactions, CN.with_db(lst))
    # Mutual neutralisation H⁻ + H⁺ → 2H (GA08 k5): its products are mostly excited H (n = 2, 3), so the detailed-balance reverse to ground-state H + H
    # (collisional ion-pair formation, 4.6e-15 cm³ s⁻¹ at 1e4 K by k_f/K) is far too large (130× the electron-impact ionisation at x_e = 0.05, spurious
    # ionisation of warm atomic gas): one-way here, likewise H₂⁺ + e → 2H (associative ionisation 2H → H₂⁺ + e, 11 eV endothermic, same argument), its reverse only for the J = B test (`jb_inverse`).
    if !jb_inverse
        filter!(r -> !((startswith(r.name, "H- + H+ → 2H") || startswith(r.name, "H2+ + e → 2H")) && endswith(r.name, "(reverse, detailed balance)")), rx)
    end
    # H⁻ + e → H + 2e (AAZN97 k14) with the Arrhenius continuation below its 3e3 K limit, for the exact reverse
    filter!(r -> !startswith(r.name, "H- + e → H + 2e"), rx)
    k14(T) = exp(sum(c * log(T / 11605)^(i - 1) for (i, c) in enumerate((-18.01849334, 2.3608522, -0.28274430, 1.62331664e-2, -3.36501203e-2, 1.17832978e-2, -1.65619470e-3, 1.06827520e-4, -2.63128581e-6))))
    append!(rx, CN.with_db([(CR.R("H- + e → H + 2e (AAZN97 k14, Arrhenius below 3e3 K)", [:Hm, :e], [:H, :e, :e], k(arrh_ext(k14, 3e3, 8750.0))), true)]))
    if li       # Li ionisation continued below 1e4 K like H (exact detailed balance of the reverse)
        filter!(r -> !startswith(r.name, "Li + e → Li+ + 2e"), rx)
        append!(rx, CN.with_db([(CR.R("Li + e → Li+ + 2e (Voronov 1997, Arrhenius below 1e4 K)", [:Li, :e], [:Lip, :e, :e], k(arrh_ext(T -> CR.voronov(T, 5.4, 0.0, 1.39e-7, 0.438, 0.41), 1e4, 62600.0))), true)]))
    end
    if jb_inverse      # inverses of the one-way radiative processes (H₂⁺ association, Li + H⁺ → Li⁺ + H + γ) at J = B: k/K; valid only for a thermal radiation field (test a)
        for r in copy(rx)
            (startswith(r.name, "H + H+ → H2+ + γ") || startswith(r.name, "Li + H+ → Li+ + H + γ")) &&
                push!(rx, CN.Reaction(r.name * " (inverse at J = B)", r.P, r.R, CN.reverse_k(r)))
        end
    end
    photos = CN.Photo[p for p in ph0 if p.R !== :H]
    append!(photos, photo_pair("H(n ≤ $NL) + γ → H+ + e (Verner 1s + Kramers)", :H, :Hp, CR.levels_H(), CN.ZH_el))
    append!(photos, photo_pair("He I + γ → He+ + e (Verner 1996 + Kramers)", :He, :Hep, levels_HeI(), CN.ZHe_el))
    append!(photos, photo_pair("He II + γ → He2+ + e (Verner 1996 + Johnson 1972)", :Hep, :Hepp, levels_HeII(), CN.ZHeII_el))
    net = CN.Network(vcat(net0.species, [:Hep, :Hepp]), vcat(net0.free, [:Hep, :Hepp]), rx)
    net, photos
end

# ------------------------------------------------------------------ cooling
const πe2mc = 0.026540                                   # πe²/(m_e c) [cm² s⁻¹]
const EXCFRAC = Ref(1.0)                                 # sensitivity switch: scale of the collisional excitation cooling
"e^y E₁(y) (Abramowitz & Stegun 5.1.53 for y ≤ 1, 5.1.56 above)."
function expE1(y)
    if y <= 1
        E1 = -log(y) - 0.57721566 + 0.99999193y - 0.24991055y^2 + 0.05519968y^3 - 0.00976004y^4 + 0.00107857y^5
        return exp(y) * E1
    end
    (y^2 + 2.334733y + 0.250621) / (y^2 + 3.330657y + 1.681534) / y
end
"Absorption oscillator strength 1 → n of a hydrogenic ion (Kramers × Johnson 1972 Gaunt factor); Lyα 0.4162."
function f1n(n)
    g0, g1, g2 = (1.1330, -0.4059, 0.07014); x = 1 - 1 / n^2
    32 / (3π * sqrt(3)) / n^3 * (1 - 1 / n^2)^(-3) * (g0 + g1 / x + g2 / x^2)
end
struct Line
    sp::Symbol          # lower-level species (H, He, Hep)
    E::Float64          # [erg]
    f::Float64
    m::Float64          # mass of the absorber [g]
    q::Function         # T -> excitation rate coefficient by electrons [cm³ s⁻¹]
end
const NEXC = 8
"""
Collisional-excitation lines (electrons; ground → n):
H I and He II: n = 2 from Black (1981) / Cen (1992) [q₂ = 4.59e-8 e^{−118348/T}/(1 + √(T/1e5)), He II: 8.47e-7 T^{−0.397} e^{−473638/T}/(1 + √T₅)], n = 3…8 from the
Bethe/Coulomb scaling of the collision strength, Υ₁ₙ = Υ₁₂ (f₁ₙ/f₁₂)(E₁₂/E₁ₙ) with the same high-T correction;
He I: singlet 1¹S → n¹P, n = 2–4 (f = 0.2762, 0.0734, 0.0299), van Regemorter (1962) with ḡ = max(0.276 e^y E₁(y), 0.05) (floor from the H Lyα comparison with
Cen); the spin-forbidden and metastable channels (2³S, 2¹S, 2³P) are omitted.
"""
function exc_lines()
    L = Line[]
    mH = CN.mH; mHe = 4.0026 * 1.66053907e-24
    f2 = f1n(2)
    for n in 2:NEXC
        En = CN.χH * (1 - 1 / n^2); sc = f1n(n) / f2 * (1 - 1 / 4) / (1 - 1 / n^2)
        qH = T -> 4.59e-8 / (1 + sqrt(T / 1e5)) * exp(-118348.0 / T) * (n == 2 ? 1.0 : sc * exp(-(En - 0.75 * CN.χH) / (kB * T)))
        push!(L, Line(:H, En, f1n(n), mH, qH))
        E2 = 4En
        qHe2 = T -> 8.47e-7 * T^-0.397 / (1 + sqrt(T / 1e5)) * exp(-473638.0 / T) * (n == 2 ? 1.0 : sc * exp(-(4En - 3CN.χH) / (kB * T)))
        push!(L, Line(:Hep, E2, f1n(n), mHe, qHe2))
    end
    for (f, EeV) in ((0.2762, 21.218), (0.0734, 23.087), (0.0299, 23.742))
        E = EeV * eV
        q = T -> (y = E / (kB * T); 1.2521e-4 / sqrt(T) * (CN.χH / E) * f * max(0.276 * expE1(y), 0.05) * exp(-y))
        push!(L, Line(:He, E, f, mHe, q))
    end
    L
end
const LINES = exc_lines()
"Line escape probability (1 − e^{−τ})/τ for a transition of oscillator strength f and photon energy E from a level of density nl, Doppler width from T; 1 for r = 0 (thin)."
function β_line(f, E, nl, T, m)
    ZS.r > 0 || return one(T)
    ν0 = E / hP; b = sqrt(2kB * T / m)
    ℓ = min(sqrt(π) * b * ZS.Lv, ZS.r)
    βesc(πe2mc * f * nl * ℓ / (sqrt(π) * ν0 * b / c))
end
"Net factor 1 − J(ν₀)/B_ν₀(T) (stimulated absorption of the radiation field in the line; 1 for J = 0)."
function fJ(E, T, env)
    ν = E / hP; J = env.J(ν); J == 0 && return one(T)
    1 - J / (2hP * ν^3 / c^2 / expm1(E / (kB * T)))
end
"Resonance-line energy and f of the cascade end for each species (H Lyα, He I 584 Å, He II 304 Å)."
const CASC = Dict(:H => (10.2 * eV, 0.4162, CN.mH), :He => (21.218 * eV, 0.2762, 4.0026 * 1.66053907e-24), :Hep => (40.81 * eV, 0.4162, 4.0026 * 1.66053907e-24))
const E_DR = 40.81 * eV + CN.χHe1                         # photon energy radiated per dielectronic recombination (see the notes)
"Gaunt factor of the free–free emission (Cen 1992 / Katz et al. 1996 fit)."
gff(T) = 1.1 + 0.34 * exp(-(5.5 - log10(T))^2 / 3)

"""
    parts(ρ, T, d; photos, env, eqcomp = nothing) -> NamedTuple (erg cm⁻³ s⁻¹)

Terms of the atomic cooling at density ρ, temperature T, composition d (Dict of number fractions per H nucleus):
`rad` net radiative loss of the photo pairs (recombination photons incl. the cascade lines × β, minus absorbed J_ν; ground-level continua times their
on-the-spot escape fraction), `exc` collisional excitation (H, He I, He II) × escape probability × (1 − J/B), `dr` dielectronic recombination, `ff`, `comp` (Compton).
Collisional ionisation and three-body recombination carry no term here: their energy is in E0 of the species (ChemNet.eos).
"""
function parts(ρ, T, d; photos, env)
    nH = ρ / CN.μH(); g(s) = get(d, s, zero(T)); ne = nH * g(:e); Tv = ForwardDiff.value(T)
    rad = zero(T)
    for p in photos
        RE, AE, AX = CN.photoenergy_t(p, T, env); es = p.esc(Tv, env)
        xion = g(p.P[1]); xX = g(p.R)
        βc = one(T); fj = one(T)
        if haskey(CASC, p.R)
            Ec, fc, mc = CASC[p.R]; βc = β_line(fc, Ec, nH * xX, T, mc); fj = fJ(Ec, T, env)
        end
        rad += nH^2 * xion * g(:e) * es * (AE + βc * fj * AX) - nH * xX * RE
    end
    exc = zero(T)
    for l in LINES
        nl = nH * g(l.sp); nl == 0 && continue
        exc += ne * nl * l.q(T) * l.E * β_line(l.f, l.E, nl, T, l.m) * fJ(l.E, T, env)
    end
    exc *= EXCFRAC[]
    dr = nH^2 * g(:e) * g(:Hep) * k_di_HeII(T) * E_DR
    ff = 1.42e-27 * gff(T) * sqrt(T) * ne * nH * (g(:Hp) + g(:Hep) + 4g(:Hepp))
    comp = env.Trad > 0 ? 1.0178e-37 * env.Trad^4 * (T - env.Trad) * ne : zero(T)
    (rad = rad, exc = exc, dr = dr, ff = ff, comp = comp)
end
"""
    cooling(photos, env; eqcomp = nothing) -> (ρ, T, d) -> Λ

Total atomic cooling [erg cm⁻³ s⁻¹] for ChemNet.step (differentiable in T and the composition). `d === nothing` (the equilibrium-table regime of
PopIIIMaterial.NetGas) uses `eqcomp(ρ, T)` → Dict.
"""
function cooling(photos, env; eqcomp = nothing)
    (ρ, T, d) -> begin
        d === nothing && (d = eqcomp(ρ, T))
        p = parts(ρ, T, d; photos, env)
        p.rad + p.exc + p.dr + p.ff + p.comp
    end
end

end # module
