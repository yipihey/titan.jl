"""
    ChemNet

Population III project, G3: a general reaction-network engine for metal-free gas (replaces the hand-coded network of
PrimordialChem for new work; PrimordialChem's RHFC02 network stays as the reference it must reproduce).

- **Species** carry mass, charge, element content, an energy zero E0 (relative to free nuclei and electrons at rest —
  only differences within an element matter) and an internal partition function Z_int(T) (electronic × nuclear-spin ×
  rotational × vibrational; nuclear spins kept consistently so that equilibrium constants are exact).
- **Equilibrium constants** for any reaction from the species: K = Π_P (Z_tr Z_int e^{−E0/kT}) / Π_R (…), Z_tr =
  (2π m k T/h²)^{3/2}.
- **Reactions** as data: reactants, products, rate function k(T, env) (env: n_H, J_ν access, …). `db = true` adds the
  reverse reaction with k_rev = k_fwd / K (detailed balance), so the network's own steady state is the Saha/K equilibrium.
- **Photo-processes** from cross sections σ(ν) and the mean intensity J_ν supplied per zone (`env.J(ν)`; default: a
  blackbody at `env.Trad`): photo-rate 4π∫σJ/(hν)dν; the inverse (radiative attachment/recombination) from the Milne
  relation with stimulated emission, α = (n_A/(n_B n_e))_eq 4π∫σ/(hν)(2hν³/c² + J_ν)e^{−hν/kT}dν — so J = B gives the
  equilibrium exactly.
- **Unknowns**: number fractions x_i = n_i/n_H of all species except the eliminated ones — one per element (H, He, Li:
  from nuclei conservation) and e⁻ (charge). Conservation holds to round-off by construction.
- **Energy**: u = Σ x_i (c_v,i T + E0_i)/μ_H per gram (c_v: 3/2 k for atoms/ions/e⁻, RHFC02 A9 for H₂, 5/2 k for H₂⁺) —
  the generalisation of PrimordialChem.eos_rhfc02.
- **Step**: backward Euler in x with the rates at the new T plus u(ρ, T, x) = u − dt Λ/ρ (chemistry and T together),
  Newton with the exact Jacobian (ForwardDiff), damping that keeps every abundance (incl. the eliminated ones) ≥ 0,
  recursive halving.
"""
module ChemNet

using ForwardDiff, Printf
import LinearAlgebra

const kB = 1.380649e-16; const hP = 6.62607015e-27; const c = 2.99792458e10; const eV = 1.602176634e-12
const mH = 1.6735575e-24; const me = 9.1093837e-28; const mp = 1.67262192e-24
const yHe = 1 / 12; const yLi = 4.15e-10

# ------------------------------------------------------------------ species
struct Species
    name::Symbol
    m::Float64            # mass [g]
    q::Int                # charge
    el::NTuple{3,Int}     # (H, He, Li) nuclei
    E0::Float64           # energy zero [erg]
    Zint::Function        # T -> internal partition function (incl. nuclear spin)
    cv::Function          # T -> heat capacity per particle in units of k (thermal energy = ∫ c_v dT ≈ c_v T, see u_th)
end
const θrot_H2 = 85.3; const θvib_H2 = 6332.0
function Zrot_H2(T)
    Zr = zero(T); for J in 0:60; Zr += (iseven(J) ? 1 : 3) * (2J + 1) * exp(-J * (J + 1) * θrot_H2 / T); end; Zr
end
Zvib(T, θ) = 1 / (1 - exp(-θ / T))
"RHFC02 A9 (ON98): 1/(γ_H₂ − 1) = 5/2 + x² eˣ/(eˣ − 1)², x = 6100 K/T — energy per H₂ (thermal) = this × kT."
u_H2(T) = (x = 6100 / T; x < 300 ? 5 / 2 + x^2 * exp(x) / expm1(x)^2 : 5 / 2 + x^2 * exp(-x) / (1 - exp(-x))^2)     # (S6c: the second form, identical, below ≈ 20 K where exp(x)² overflows — Newton overshoots there in the stiff atomic-cooling regime)
const χH = 13.598eV; const EA_H = 0.754eV; const D0_H2 = 4.478eV; const D0_H2p = 2.65eV
const χLi = 5.3917eV
const χHe1 = 24.5874eV; const χHe2 = 54.4178eV; const Ry_eV = 13.605693
"Levels of atomic H kept in its partition function and photo-ionisation (n ≤ NH_LEVELS; the sum diverges without a
cut-off — at the densities where excited levels matter they are pressure-dissolved and SCVH95 takes over)."
const NH_LEVELS = 15
"""
Li I levels for the partition function and photo-ionisation: ground 2s (g = 2) plus quantum-defect hydrogenic levels
n ≤ NH_LEVELS, l ≤ 3, E_nl = χ_Li − Ry/(n − δ_l)², δ_s = 0.40, δ_p = 0.05, δ_d = δ_f = 0, g = 2(2l + 1)
[approximate: 2p at 1.85 eV, 3s 3.37 eV — the measured 1.848 and 3.373 eV]. Returns (E_exc [erg], g, n*) per level.
"""
const LI_LEVELS = let lv = [(0.0, 2.0, 2 - 0.40)]
    δ = (0.40, 0.05, 0.0, 0.0); Ry = 13.6057 * 1.602176634e-12
    for n in 2:NH_LEVELS, l in 0:min(n - 1, 3)
        (n == 2 && l == 0) && continue
        ns = n - δ[l+1]; push!(lv, (χLi - Ry / ns^2, 2.0 * (2l + 1), ns))
    end
    lv
end
ZLi_el(T) = (z = zero(T); for (E, g, _) in LI_LEVELS; z += g * exp(-E / (kB * T)); end; z)
"Electronic partition function of H, levels n ≤ NH_LEVELS (g = 2n², E_n = χ(1 − 1/n²))."
ZH_el(T) = (z = zero(T); for n in 1:NH_LEVELS; z += 2n^2 * exp(-χH * (1 - 1 / n^2) / (kB * T)); end; z)
"H₂⁺ rotation (B = 29.9 cm⁻¹ → θ = 43 K) and vibration (2297 cm⁻¹ → 3305 K), ortho/para weights 3/1, electron spin 2 [approximate]."
Zint_H2p(T) = (Zr = zero(T); for J in 0:60; Zr += (iseven(J) ? 1 : 3) * (2J + 1) * exp(-J * (J + 1) * 43.0 / T); end; 2Zr * Zvib(T, 3305.0))

"""
He I levels (S6c) for the partition function and the photo pairs: (E_exc [erg], g) — n = 1 ground (1 ¹S), n = 2 and 3 term-resolved
(NIST: 2³S 19.820, 2¹S 20.616, 2³P 20.964, 2¹P 21.218; 3³S 22.718, 3¹S 22.920, 3³P 23.007, 3³D 23.074, 3¹D 23.074, 3¹P 23.087 eV),
n = 4…NH_LEVELS lumped (g = 4n², E = χ_He − Ry/n²).
"""
const HE_LEVELS = let lv = Tuple{Float64,Float64,Int}[(0.0, 1.0, 1)]      # (E, g, n)
    for (E, g) in ((19.8196, 3.0), (20.6158, 1.0), (20.9642, 9.0), (21.2180, 3.0)); push!(lv, (E * eV, g, 2)); end
    for (E, g) in ((22.7185, 3.0), (22.9203, 1.0), (23.0070, 9.0), (23.0739, 15.0), (23.0743, 5.0), (23.0870, 3.0)); push!(lv, (E * eV, g, 3)); end
    for n in 4:NH_LEVELS; push!(lv, (χHe1 - Ry_eV * eV / n^2, 4.0n^2, n)); end
    lv
end
ZHe_el(T) = (z = zero(T); for (E, g, _) in HE_LEVELS; z += g * exp(-E / (kB * T)); end; z)
"Electronic partition function of He II (hydrogenic Z = 2, g = 2n², E_n = χ_He2 (1 − 1/n²), n ≤ NH_LEVELS)."
ZHeII_el(T) = (z = zero(T); for n in 1:NH_LEVELS; z += 2n^2 * exp(-χHe2 * (1 - 1 / n^2) / (kB * T)); end; z)

const SPECIES = Dict{Symbol,Species}(
    # nuclear spin: H nucleus 2 states; Li-7 nucleus 4; He-4 1 (only ratios within an element matter)
    :H   => Species(:H, mH, 0, (1, 0, 0), -χH, T -> 2 * ZH_el(T), T -> 1.5),     # × nuclear spin 2
    :Hp  => Species(:Hp, mp, 1, (1, 0, 0), 0.0, T -> 2 + zero(T), T -> 1.5),
    :Hm  => Species(:Hm, mH + me, -1, (1, 0, 0), -(χH + EA_H), T -> 2 + zero(T), T -> 1.5),
    :H2  => Species(:H2, 2mH, 0, (2, 0, 0), -(2χH + D0_H2), T -> Zrot_H2(T) * Zvib(T, θvib_H2), T -> NaN),   # c_v via u_H2
    :H2p => Species(:H2p, 2mH - me, 1, (2, 0, 0), -(χH + D0_H2p), Zint_H2p, T -> 2.5),
    :e   => Species(:e, me, -1, (0, 0, 0), 0.0, T -> 2 + zero(T), T -> 1.5),
    # S6c: He with its term-resolved electronic partition function (cancels wherever He is a spectator on both sides); E0 relative to neutral He
    :He  => Species(:He, 4.0026 * 1.66053907e-24, 0, (0, 1, 0), 0.0, ZHe_el, T -> 1.5),
    :Hep => Species(:Hep, 4.0026 * 1.66053907e-24 - me, 1, (0, 1, 0), χHe1, ZHeII_el, T -> 1.5),
    :Hepp => Species(:Hepp, 4.0026 * 1.66053907e-24 - 2me, 2, (0, 1, 0), χHe1 + χHe2, T -> 1 + zero(T), T -> 1.5),
    :Li  => Species(:Li, 7.016 * 1.66053907e-24, 0, (0, 0, 1), -χLi, T -> 4 * ZLi_el(T), T -> 1.5),     # × nuclear spin 4
    :Lip => Species(:Lip, 7.016 * 1.66053907e-24 - me, 1, (0, 0, 1), 0.0, T -> 4 + zero(T), T -> 1.5),
)
Ztr(s::Species, T) = (2π * s.m * kB * T / hP^2)^1.5
"Equilibrium constant Π_P n / Π_R n [cm^{-3 Δν}] of R → P (vectors of species names, repeats allowed)."
Keq(R, P, T) = exp(lnKeq(R, P, T))
"ln Keq — use this where K itself under- or overflows (ionisation equilibria below a few hundred K)."
function lnKeq(R, P, T)
    lnK = zero(T)
    for p in P; s = SPECIES[p]; lnK += log(Ztr(s, T) * s.Zint(T)) - s.E0 / (kB * T); end
    for r in R; s = SPECIES[r]; lnK -= log(Ztr(s, T) * s.Zint(T)) - s.E0 / (kB * T); end
    lnK
end

# ------------------------------------------------------------------ reactions
struct Reaction
    name::String
    R::Vector{Symbol}; P::Vector{Symbol}
    k::Function           # (T, env) -> rate coefficient [cm^{3(n_R−1)} s⁻¹]
end
"""
Photo-process X + γ → products from the levels of X: `levels` = [(E_x, g, E_th, σ(ν)), …] (excitation energy above X's
ground state, statistical weight, threshold, cross section), populated in Boltzmann ratio g e^{−E_x/kT}/Z(T) with Z the
electronic partition function of X (the same as in X's Zint, so that the Milne pair gives Saha exactly).
"""
struct Photo
    name::String
    R::Symbol; P::Vector{Symbol}
    levels::Vector{Tuple{Float64,Float64,Float64,Function}}
    Z::Function
    milne::Bool           # add the inverse (radiative attachment/recombination) from the Milne relation
    esc::Function         # S6c: (T, env) -> escape fraction of the recombination photons (on-the-spot: 0 = absorbed locally, case B); 1 = thin
end
Photo(name, R, P, levels, Z, milne::Bool) = Photo(name, R, P, levels, Z, milne, (T, env) -> 1.0)
"Single-level process with cross section σ(ν) above Eth."
Photo(name, R, P, σ::Function, Eth::Real, milne::Bool) = Photo(name, R, P, [(0.0, 1.0, Float64(Eth), σ)], T -> 1.0, milne, (T, env) -> 1.0)
struct Network
    species::Vector{Symbol}          # all species in the network (incl. eliminated; He as inert spectator)
    free::Vector{Symbol}             # the unknowns
    reactions::Vector{Reaction}
end

"""
Reaction list with detailed-balance reverses added for the entries flagged `db`: entries (r, db) or (r, db, Tmin), Tmin
the lower validity limit of r's fit (below it the reverse is extrapolated, see `reverse_k`).
"""
function with_db(list)
    out = Reaction[]
    for e in list
        r, db = e[1], e[2]; Tmin = length(e) > 2 ? e[3] : 0.0
        push!(out, r)
        db && push!(out, Reaction(r.name * " (reverse, detailed balance)", r.P, r.R, reverse_k(r; Tmin)))
    end
    out
end
"""
k_r = k_f / K, formed in logarithms. Below T_u = max(Tmin, T where k_f underflows or is cut off) k_r — the exothermic
direction, close to a power law in T — is extrapolated from [T_u, 2T_u] as a power law. Tmin matters: k_f/K amplifies any
error of an endothermic fit by e^{ΔE/kT} (AAZN97 k1 below 1e4 K gave three-body recombination of H⁺ in 1e-9 s at 200 K).
"""
function reverse_k(r::Reaction; Tmin = 0.0)
    Tu = Ref(NaN)
    lnkr(T, env) = log(r.k(T, env)) - lnKeq(r.R, r.P, T)
    function k(T, env)
        kf = r.k(T, env)
        T >= Tmin && kf > 1e-290 && return exp(log(kf) - lnKeq(r.R, r.P, T))
        iszero(kf) || kf > 0 || return zero(T)                 # a fit that is ≤ 0 (clipped) stays 0
        if isnan(Tu[])
            t = max(Float64(Tmin), ForwardDiff.value(T)); while r.k(t, env) <= 1e-290 && t < 1e6; t *= 1.05; end
            Tu[] = t
        end
        Tu[] >= 1e6 && return zero(T)
        a = lnkr(Tu[], env); s = (lnkr(2Tu[], env) - a) / log(2.0)
        exp(a + s * log(T / Tu[]))
    end
end

"Photo-rate [s⁻¹] 4π∫σJ/(hν)dν and the Milne inverse coefficient α for A + e → X + γ (stimulated emission included)."
function photorates(p::Photo, T, env; n = 200)
    # α = Σ_levels (n_l/Π n_P)_eq × 4π∫σ_l/(hν)(2hν³/c² + J) e^{−hν/kT} dν, (n_l/Π n_P)_eq = g_l e^{−E_x/kT}/(Z K(X → P));
    # the Boltzmann factors and 1/K are combined in one exponent (each alone under/overflows below a few hundred K)
    lnK = lnKeq([p.R], p.P, T); Z = p.Z(T); R = 0.0; A = 0.0
    for (Ex, g, Eth, σf) in p.levels
        ν0 = Eth / hP; νs = [ν0 * exp(x) for x in range(0.0, log(1 + 60kB * max(T, env.Trad) / Eth + 30.0), length = n + 1)]
        w = g / Z; Rl = 0.0; Al = 0.0
        for i in 1:n
            ν = 0.5 * (νs[i] + νs[i+1]); dν = νs[i+1] - νs[i]; σ = σf(ν); J = env.J(ν)
            Rl += 4π * σ * J / (hP * ν) * dν
            Al += 4π * σ / (hP * ν) * (2hP * ν^3 / c^2 + J) * exp(-(Ex + hP * ν) / (kB * T) - lnK) * dν
        end
        R += w * exp(-Ex / (kB * T)) * Rl; A += w * Al
    end
    R, A
end

"""
Energy moments of the photo pair (S6c): (RE, AE, AX) with RE = Σ_l w_l e^{-Ex/kT} 4π∫σ_l J dν [erg s⁻¹ per atom, absorbed from the radiation field]
and AE = Σ_l w_l 4π∫σ_l/(hν)(2hν³/c² + J) e^{−(Ex+hν)/kT − lnK} hν dν [erg cm³ s⁻¹ per ion·electron, energy of the emitted recombination photons
incl. stimulated emission], AX = Σ_l Ex_l α_l [erg cm³ s⁻¹: the excitation energy of the recombined level, radiated by its cascade]. Net radiative loss per volume of the process: n_ion n_e AE − n_X RE (= 0 for J = B: the Milne pair).
"""
function photoenergy(p::Photo, T, env; n = 200)
    lnK = lnKeq([p.R], p.P, T); Z = p.Z(T); R = 0.0; A = 0.0; AX = 0.0
    for (Ex, g, Eth, σf) in p.levels
        ν0 = Eth / hP; νs = [ν0 * exp(x) for x in range(0.0, log(1 + 60kB * max(T, env.Trad) / Eth + 30.0), length = n + 1)]
        w = g / Z; Rl = 0.0; Al = 0.0; Ar = 0.0
        for i in 1:n
            ν = 0.5 * (νs[i] + νs[i+1]); dν = νs[i+1] - νs[i]; σ = σf(ν); J = env.J(ν)
            Rl += 4π * σ * J * dν
            a = 4π * σ / (hP * ν) * (2hP * ν^3 / c^2 + J) * exp(-(Ex + hP * ν) / (kB * T) - lnK) * dν
            Ar += a; Al += a * hP * ν
        end
        R += w * exp(-Ex / (kB * T)) * Rl; A += w * Al; AX += w * Ex * Ar
    end
    R, A, AX
end
const PETAB = IdDict{Any,Any}()
"photoenergy through a table in log T (static radiation field only), Catmull–Rom in ln, differentiable in T."
function photoenergy_t(p::Photo, T, env)
    hasproperty(env, :energy_fn) && return env.energy_fn(p, ForwardDiff.value(T))      # S7a: per-zone tabulated Milne sums (RadTransfer.ZoneEnv)
    (hasproperty(env, :static) && env.static) || return photoenergy(p, ForwardDiff.value(T), env)
    tab = get!(PETAB, (p, env)) do
        lT = collect(0.0:0.01:6.0); v = [photoenergy(p, 10^x, env) for x in lT]
        (lT, [log(max(a[1], 1e-300)) for a in v], [log(max(a[2], 1e-300)) for a in v], [log(max(a[3], 1e-300)) for a in v])
    end
    lT, lR, lA, lX = tab; n = length(lT); x = clamp(log10(T), 0.0, 6.0); i = clamp(floor(Int, ForwardDiff.value(x) / 0.01) + 1, 2, n - 2); w = (x - lT[i]) / 0.01
    cr(f) = (let p0 = f[i-1], p1 = f[i], p2 = f[i+1], p3 = f[i+2]; p1 + 0.5w * (p2 - p0 + w * (2p0 - 5p1 + 4p2 - p3 + w * (3(p1 - p2) + p3 - p0))) end)
    R = exp(cr(lR)); A = exp(cr(lA)); AX = exp(cr(lX))
    (R < 1e-290 ? zero(R) : R), A, (AX < 1e-290 ? zero(AX) : AX)
end

"""
photorates through a table in log T (1–1e6 K, 0.005 dex, Catmull–Rom in ln R and ln α) when the radiation field is the same in
every zone (`env.static == true`, e.g. the CMB or none); otherwise direct. Table vs direct ≤ 3e-3 at 20 K–9e5 K (H, H⁻, Li; J = B(54 K)); the larger differences, H⁻ α above 7e4 K, are the direct quadrature's own (n = 200 vs 800: 3 %).
"""
const PHTAB = IdDict{Any,Any}()
function photorates_t(p::Photo, T, env)
    hasproperty(env, :rates_fn) && return env.rates_fn(p, T)                          # S7a: per-zone tabulated rates (RadTransfer.ZoneEnv)
    (hasproperty(env, :static) && env.static) || return photorates(p, T, env)
    tab = get!(PHTAB, (p, env)) do
        lT = collect(0.0:0.005:6.0)
        v = [photorates(p, 10^x, env) for x in lT]
        (lT, [log(max(a[1], 1e-300)) for a in v], [log(max(a[2], 1e-300)) for a in v])
    end
    lT, lR, lA = tab; n = length(lT); x = clamp(log10(T), 0.0, 6.0); i = clamp(Int(floor(x / 0.005)) + 1, 2, n - 2); w = (x - lT[i]) / 0.005
    function cr(f)                                                 # Catmull–Rom
        p0, p1, p2, p3 = f[i-1], f[i], f[i+1], f[i+2]
        p1 + 0.5w * (p2 - p0 + w * (2p0 - 5p1 + 4p2 - p3 + w * (3(p1 - p2) + p3 - p0)))
    end
    R = exp(cr(lR)); A = exp(cr(lA))
    (R < 1e-290 ? 0.0 : R), A
end

"(p, R, α·esc) for the dxdt photo list: the recombination (inverse) coefficient scaled by the on-the-spot escape fraction of p."
function photoset(p::Photo, T, env, f)
    R, A = f(p, T, env); (p, R, A * p.esc(T, env))
end

# ------------------------------------------------------------------ composition
"Eliminated species: H (H nuclei), Li (Li nuclei), e (charge); He fixed at y_He."
function full(x, net)
    d = Dict{Symbol,eltype(x)}()
    for (i, s) in enumerate(net.free); d[s] = x[i]; end
    # S6c: with He ionisation states among the free species He itself is eliminated (He nuclei conservation)
    heE = :Hep in net.species || :Hepp in net.species
    heE ? (d[:He] = yHe - sum((SPECIES[s].el[2] * d[s] for s in net.free if SPECIES[s].el[2] > 0); init = zero(eltype(x)))) : get!(d, :He, yHe + zero(eltype(x)))
    hasLi = :Li in net.species
    H = 1 - sum((SPECIES[s].el[1] * d[s] for s in net.free if s != :H); init = zero(eltype(x)))
    d[:H] = H
    hasLi && (d[:Li] = yLi - sum((SPECIES[s].el[3] * d[s] for s in net.free if SPECIES[s].el[3] > 0); init = zero(eltype(x))))
    d[:e] = sum((SPECIES[s].q * d[s] for s in keys(d) if s != :e); init = zero(eltype(x)))
    d
end
μH() = mH * (1 + 4yHe)    # mass per H nucleus (Li negligible)

"""
Index form of a network (built once per Network, cached): species vector = net.species ∪ {He}; free-species positions;
reactant/product index lists; element and charge coefficients for the eliminated H, Li, e.
"""
struct Compiled
    sp::Vector{Symbol}; idx::Dict{Symbol,Int}; free::Vector{Int}
    R::Vector{Vector{Int}}; P::Vector{Vector{Int}}
    iH::Int; iLi::Int; ie::Int; iHe::Int
    elH::Vector{Float64}; elLi::Vector{Float64}; q::Vector{Float64}       # per free species
    elHe::Vector{Float64}; heE::Bool
    cv::Vector{Function}; E0::Vector{Float64}; isH2::Vector{Bool}
end
const COMPILED = IdDict{Any,Compiled}()
compiled(net) = get!(COMPILED, net) do
    sp = unique(vcat(net.species, [:He])); idx = Dict(s => i for (i, s) in enumerate(sp))
    fr = [idx[s] for s in net.free]
    Compiled(sp, idx, fr, [[idx[s] for s in r.R] for r in net.reactions], [[idx[s] for s in r.P] for r in net.reactions],
             get(idx, :H, 0), get(idx, :Li, 0), get(idx, :e, 0), idx[:He],
             [s == :H ? 0.0 : Float64(SPECIES[s].el[1]) for s in net.free], [Float64(SPECIES[s].el[3]) for s in net.free],
             [Float64(SPECIES[s].q) for s in net.free],
             [Float64(SPECIES[s].el[2]) for s in net.free], (:Hep in net.species || :Hepp in net.species),
             [SPECIES[s].cv for s in sp], [SPECIES[s].E0 for s in sp], [s == :H2 for s in sp])
end
"Full composition as a vector over compiled(net).sp (same elimination as `full`)."
function fullvec(x, c::Compiled)
    v = zeros(eltype(x), length(c.sp))
    for (i, j) in enumerate(c.free); v[j] = x[i]; end
    v[c.iHe] = c.heE ? yHe - sum(c.elHe[i] * x[i] for i in eachindex(x)) : yHe
    c.iH > 0 && (v[c.iH] = 1 - sum(c.elH[i] * x[i] for i in eachindex(x)))
    c.iLi > 0 && (v[c.iLi] = yLi - sum(c.elLi[i] * x[i] for i in eachindex(x)))
    c.ie > 0 && (v[c.ie] = sum(c.q[i] * x[i] for i in eachindex(x)))
    v
end

"Internal energy per gram and pressure."
function eos(ρ, T, x, net)
    c = compiled(net); v = fullvec(x, c); nH = ρ / μH()
    u = zero(T); N = zero(T)
    for j in eachindex(v)
        u += v[j] * ((c.isH2[j] ? u_H2(T) : c.cv[j](T)) * kB * T + c.E0[j])
        N += v[j]
    end
    (u / μH(), N * nH * kB * T)
end

"dx/dt for the free species at H-nucleus density nH."
function dxdt(x, nH, T, net, env, photo)
    c = compiled(net); v = fullvec(x, c); dd = zeros(promote_type(eltype(x), typeof(T)), length(v))
    for (ir, r) in enumerate(net.reactions)
        Ri = c.R[ir]; rate = r.k(T, env) * nH^(length(Ri) - 1)
        for j in Ri; rate *= v[j]; end
        for j in Ri; dd[j] -= rate; end
        for j in c.P[ir]; dd[j] += rate; end
    end
    for (p, R, α) in photo                          # precomputed per zone: photo-rate R, inverse coefficient α
        jr = c.idx[p.R]; jp = [c.idx[s] for s in p.P]
        rate = R * v[jr]; for j in jp; dd[j] += rate; end; dd[jr] -= rate
        inv_ = α * nH^(length(jp) - 1); for j in jp; inv_ *= v[j]; end    # per H nucleus: α n_H^(k−1) Π x
        for j in jp; dd[j] -= inv_; end; dd[jr] += inv_
    end
    dd[c.free]
end

# S6b: floor of the abundance in the relative convergence test of `step` (default 1e-20 = as before). Species at 1e-18 (H₂⁺ in cold dense gas) have round-off of the
# Newton solve far above 1e-10 of their value: the step then "failed to converge" at every substep (the sink runs of popiii_m2c.jl use 1e-12).
const XFLOOR = Ref(1e-20)
const NESTED = Ref(false)              # S6c (opt-in): when the coupled Newton of `step` fails, solve T by bisection with the chemistry solved at each trial T (robust on the cliff of a steep cooling function) before halving dt
const NANHOOK = Ref{Any}(nothing)      # debugging: called with (z, R, J) when a Newton attempt meets a non-finite residual or Jacobian
"""
    step(x0, ρ, u, T0, dt, net; env, photos = Photo[], cooling = nothing) -> (x, T)

Chemistry and temperature together over dt at fixed ρ and u (u − dt Λ/ρ with `cooling(ρ, T, d)`, d the full composition).
"""
function step(x0, ρ, u, T0, dt, net; env, photos = Photo[], cooling = nothing, tol = 1e-10, maxit = 60)
    nH = ρ / μH(); nf = length(net.free)
    phT = Ref(NaN); phc = Ref{Any}(nothing)
    function resid(z, x0, h)
        x = z[1:nf]; T = exp(z[nf+1])
        Tv = ForwardDiff.value(T)                    # photo-rates lagged in T (Float), cached per T value
        Tv == phT[] || (phc[] = [photoset(p, Tv, env, photorates_t) for p in photos]; phT[] = Tv)
        ph = phc[]
        f = dxdt(x, nH, T, net, env, ph)
        r = x .- x0 .- h .* f
        ue = eos(ρ, T, x, net)[1] - u + (cooling === nothing ? zero(T) : h * cooling(ρ, T, full(x, net)) / ρ)
        vcat(r, ue / (kB * T / μH()))
    end
    function attempt(x0, T0, h)
        z = vcat(collect(Float64, x0), log(T0))
        for it in 1:maxit
            R = resid(z, x0, h)
            Jm = ForwardDiff.jacobian(zz -> resid(zz, x0, h), z)
            if !(all(isfinite, R) && all(isfinite, Jm))      # S6c: a non-finite residual/Jacobian (state outside a fit's range) fails the attempt, the step then halves
                NANHOOK[] === nothing || NANHOOK[](z, R, Jm)
                return (z, false)
            end
            δ = try -(Jm \ R) catch e; e isa LinearAlgebra.SingularException ? nothing : rethrow() end      # E2E-A: a singular Jacobian fails the attempt (the step halves / bisects in T) instead of aborting the run
            (δ === nothing || !all(isfinite, δ)) && return (z, false)
            λ = 1.0
            # keep every abundance ≥ 0, incl. the eliminated H, Li, e (linear in x)
            for i in 1:nf; (XFLOOR[] <= 1e-20 || z[i] > XFLOOR[]) && z[i] + λ * δ[i] < 0 && (λ = min(λ, 0.9 * z[i] / -δ[i])); end     # (XFLOOR > 1e-20: species below it are noise, not limited, clamped at 0 below)
            d0 = full(z[1:nf], net); d1 = full(z[1:nf] .+ δ[1:nf], net)
            for s in (:H, :Li, :e, :He)
                haskey(d0, s) || continue
                dv = d1[s] - d0[s]; d0[s] + λ * dv < 0 && (λ = min(λ, 0.9 * d0[s] / -dv))
            end
            abs(δ[nf+1]) * λ > 0.5 && (λ = 0.5 / abs(δ[nf+1]))
            z .+= λ .* δ
            XFLOOR[] > 1e-20 && for i in 1:nf; z[i] < 0 && (z[i] = 0.0); end
            λ == 1.0 && maximum(abs.(δ[1:nf]) ./ max.(abs.(z[1:nf]), XFLOOR[])) < tol && abs(δ[nf+1]) < tol && return (z, true)   # XFLOOR: abundance floor of the convergence test
        end
        (z, false)
    end
    function chem_fixedT(x0, T, h)          # backward-Euler chemistry at fixed T (Newton on x), nothing if it fails
        ph = [photoset(p, T, env, photorates_t) for p in photos]
        x = collect(Float64, x0)
        for it in 1:60
            F(xx) = xx .- x0 .- h .* dxdt(xx, nH, T, net, env, ph)
            Fv = F(x); Jm = ForwardDiff.jacobian(F, x); δ = try -(Jm \ Fv) catch e; e isa LinearAlgebra.SingularException ? nothing : rethrow() end; (δ === nothing || !all(isfinite, δ)) && return nothing
            λ = 1.0
            for i in 1:nf; (XFLOOR[] <= 1e-20 || x[i] > XFLOOR[]) && x[i] + λ * δ[i] < 0 && (λ = min(λ, 0.9 * x[i] / -δ[i])); end
            d0 = full(x, net); d1 = full(x .+ δ, net)
            for s in (:H, :Li, :e, :He); haskey(d0, s) || continue; dv = d1[s] - d0[s]; d0[s] + λ * dv < 0 && (λ = min(λ, 0.9 * d0[s] / -dv)); end
            x .+= λ .* δ; XFLOOR[] > 1e-20 && for i in 1:nf; x[i] < 0 && (x[i] = 0.0); end
            λ == 1.0 && maximum(abs.(δ) ./ max.(abs.(x), XFLOOR[])) < tol && return x
        end
        nothing
    end
    function nested(x0, T0, h)
        function Fe(lT)
            T = exp(lT); x = chem_fixedT(x0, T, h); x === nothing && return (NaN, nothing)
            ((eos(ρ, T, x, net)[1] - u + (cooling === nothing ? 0.0 : h * cooling(ρ, T, full(x, net)) / ρ)) / (kB * T / μH()), x)
        end
        lo = log(T0); f0, _ = Fe(lo); isfinite(f0) || return nothing
        hi = lo; dl = f0 > 0 ? -0.35 : 0.35
        for k in 1:40                             # bracket: walk in ln T until the sign changes
            hi = lo + dl; fh, _ = Fe(hi)
            if isfinite(fh) && sign(fh) != sign(f0); break; end
            isfinite(fh) && (lo = hi)
            k == 40 && return nothing
        end
        a, b = lo, hi
        for _ in 1:80
            m = 0.5 * (a + b); fm, xm = Fe(m)
            isfinite(fm) || return nothing
            sign(fm) == sign(f0) ? (a = m) : (b = m)
            abs(b - a) < 1e-13 && break
        end
        m = 0.5 * (a + b); fm, xm = Fe(m); xm === nothing && return nothing
        (xm, exp(m))
    end
    function adv(x0, T0, h, depth)
        z, ok = attempt(x0, T0, h)
        ok && return (z[1:nf], exp(z[nf+1]))
        if NESTED[]
            r = try nested(x0, T0, h) catch; nothing end
            r === nothing || return r
        end
        depth >= 12 && error("ChemNet.step: Newton failed (nH = $nH, T0 = $T0, dt = $dt)")
        x1, T1 = adv(x0, T0, h / 2, depth + 1)
        adv(x1, T1, h / 2, depth + 1)
    end
    adv(x0, T0, dt, 0)
end

"""
Steady state of the network at fixed (ρ, T) (for tests): pseudo-transient continuation — backward-Euler steps at fixed
T; a step is accepted only when its Newton converges (then dt × 3), otherwise it is retried with dt / 10; ends when an
accepted step with dt ≥ 1e30 s leaves the state unchanged.
"""
function steady(x0, ρ, T, net; env, photos = Photo[], maxsteps = 2000, xfloor = 1e-30, tol = 1e-13, xclamp = 0.0)
    nH = ρ / μH(); nf = length(net.free)
    ph = [photoset(p, T, env, photorates) for p in photos]
    x = collect(Float64, x0); dt = 1e-8
    function newton(xo, dt)
        x = copy(xo)
        F(xx) = xx .- xo .- dt .* dxdt(xx, nH, T, net, env, ph)
        for it in 1:60
            R = F(x); Jm = ForwardDiff.jacobian(F, x)
            δ = -(Jm \ R); all(isfinite, δ) || return (x, false)
            λ = 1.0
            for i in 1:nf; (xclamp == 0.0 || x[i] > xclamp) && x[i] + λ * δ[i] < 0 && (λ = min(λ, 0.9 * x[i] / -δ[i])); end    # xclamp > 0: species below it are noise, not limited, clamped at 0 below
            d0 = full(x, net); d1 = full(x .+ δ, net)
            for s in (:H, :Li, :e, :He); haskey(d0, s) || continue; dv = d1[s] - d0[s]; d0[s] + λ * dv < 0 && (λ = min(λ, 0.9 * d0[s] / -dv)); end
            x .+= λ .* δ
            xclamp > 0 && for i in 1:nf; x[i] < 0 && (x[i] = 0.0); end
            λ == 1.0 && maximum(abs.(δ) ./ max.(abs.(x), xfloor)) < tol && return (x, true)
        end
        (x, false)
    end
    for k in 1:maxsteps
        xn, ok = newton(x, dt)
        if ok
            done = dt >= 1e30 && maximum(abs.(xn .- x) ./ max.(abs.(xn), xfloor)) < 1e-10
            x = xn; dt *= 3
            done && return x
        else
            dt /= 10; dt < 1e-30 && error("ChemNet.steady: no convergence (nH = $nH, T = $T)")
        end
    end
    x
end

end # module
