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

const kB = 1.380649e-16; const hP = 6.62607015e-27; const c = 2.99792458e10; const eV = 1.602176634e-12
const mH = 1.6735575e-24; const me = 9.1093837e-28; const mp = 1.67262192e-24
const yHe = parse(Float64, get(ENV, "CHEM_YHE", string(1 / 12))); const yLi = 4.15e-10   # (CHEM_YHE=0: the pure-H benchmarks)

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
"""
RHFC02 A9 (ON98): 1/(γ_H₂ − 1) = 5/2 + x² eˣ/(eˣ − 1)², x = 6100 K/T — energy per H₂ (thermal) = this × kT.
For x ≥ 100 (T ≤ 61 K) the vibrational term is written x² e⁻ˣ/(1 − e⁻ˣ)², the same function. The first form overflows there:
its automatic derivative is NaN below 17.3 K (eˣ·(eˣ − 1)·∂x exceeds 1.8e308; a NaN Jacobian stopped the 146-zone run at step
415553 on 11 Nov 2026, when a zone at 2.6 pc reached 17.35 K) and its value is NaN below 8.6 K. For x < 100 nothing changes,
bit for bit; for x ≥ 100 the term is below 1e-39 of the 5/2.
"""
function u_H2(T)
    x = 6100 / T
    ForwardDiff.value(x) < 100 && return 5 / 2 + x^2 * exp(x) / expm1(x)^2
    em = exp(-x)
    5 / 2 + x^2 * em / (1 - em)^2
end
const χH = 13.598eV; const EA_H = 0.754eV; const D0_H2 = 4.478eV; const D0_H2p = 2.65eV
const χLi = 5.3917eV
"Levels of atomic H kept in its partition function and photo-ionisation (n ≤ NH_LEVELS; the sum diverges without a
cut-off — at the densities where excited levels matter they are pressure-dissolved and SCVH95 takes over)."
const NH_LEVELS = 15
const χHe1 = 24.587387eV; const χHe2 = 54.417763eV     # NIST ASD ionisation energies of He I and He II
"""
`HECHEM=1` (default off): helium is a chemical element of the network — He, He⁺, He²⁺ with He eliminated by nuclei
conservation (`full`, `fullvec`) and the electrons following from the charge sum. With it off, He stays the inert neutral
spectator of before (Z_int = 1, x_He = y_He), bit for bit. Energy zero: neutral He at 0, He⁺ at +χHe1, He²⁺ at +χHe1 + χHe2
(only differences within an element matter; this keeps u of every existing state unchanged).
"""
const HECHEM = get(ENV, "HECHEM", "0") == "1"
"""
He I levels for the partition function and photo-ionisation: 1¹S (g 1), the four n = 2 terms at their NIST energies
(2³S 19.8196, 2¹S 20.6158, 2³P 20.9641, 2¹P 21.2180 eV), and n = 3 … NH_LEVELS hydrogenic for the outer electron,
E_n = χHe1 − Ry/n², g = 4n² (singlets and triplets, all l) [approximate for n ≥ 3: 3³S at 22.72 eV against 23.08 hydrogenic].
Returns (E_exc [erg], g, n) per level.
"""
const HEI_LEVELS = let Ry = 13.605693 * 1.602176634e-12
    lv = [(0.0, 1.0, 1.0), (19.8196eV, 3.0, 2.0), (20.6158eV, 1.0, 2.0), (20.9641eV, 9.0, 2.0), (21.2180eV, 3.0, 2.0)]
    for n in 3:NH_LEVELS; push!(lv, (χHe1 - Ry / n^2, 4.0 * n^2, Float64(n))); end
    lv
end
ZHe_el(T) = (z = zero(T); for (E, g, _) in HEI_LEVELS; z += g * exp(-E / (kB * T)); end; z)
"He II (hydrogenic, Z = 2): g = 2n², E_n = χHe2 (1 − 1/n²), n ≤ NH_LEVELS."
ZHep_el(T) = (z = zero(T); for n in 1:NH_LEVELS; z += 2n^2 * exp(-χHe2 * (1 - 1 / n^2) / (kB * T)); end; z)
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

const SPECIES = Dict{Symbol,Species}(
    # nuclear spin: H nucleus 2 states; Li-7 nucleus 4; He-4 1 (only ratios within an element matter)
    :H   => Species(:H, mH, 0, (1, 0, 0), -χH, T -> 2 * ZH_el(T), T -> 1.5),     # × nuclear spin 2
    :Hp  => Species(:Hp, mp, 1, (1, 0, 0), 0.0, T -> 2 + zero(T), T -> 1.5),
    :Hm  => Species(:Hm, mH + me, -1, (1, 0, 0), -(χH + EA_H), T -> 2 + zero(T), T -> 1.5),
    :H2  => Species(:H2, 2mH, 0, (2, 0, 0), -(2χH + D0_H2), T -> Zrot_H2(T) * Zvib(T, θvib_H2), T -> NaN),   # c_v via u_H2
    :H2p => Species(:H2p, 2mH - me, 1, (2, 0, 0), -(χH + D0_H2p), Zint_H2p, T -> 2.5),
    :e   => Species(:e, me, -1, (0, 0, 0), 0.0, T -> 2 + zero(T), T -> 1.5),
    :He  => Species(:He, 4.0026 * 1.66053907e-24, 0, (0, 1, 0), 0.0, HECHEM ? ZHe_el : (T -> 1 + zero(T)), T -> 1.5),
    :Hep  => Species(:Hep, 4.0026 * 1.66053907e-24 - me, 1, (0, 1, 0), χHe1, ZHep_el, T -> 1.5),
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
end
"Single-level process with cross section σ(ν) above Eth."
Photo(name, R, P, σ::Function, Eth::Real, milne::Bool) = Photo(name, R, P, [(0.0, 1.0, Float64(Eth), σ)], T -> 1.0, milne)
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
    Tu = Ref(NaN); Ta = Ref(NaN); Ts = Ref(NaN)
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
        # the power law fixed at T_u: computed once (it was recomputed on every call; same values, bit-identical). The forward
        # fits of all detailed-balance pairs ignore env (the three-level atom's env-dependent factor wraps a forward reaction only).
        if isnan(Ta[])
            a = lnkr(Tu[], env); Ts[] = (lnkr(2Tu[], env) - a) / log(2.0); Ta[] = a
        end
        exp(Ta[] + Ts[] * log(T / Tu[]))
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
        # env.caseB (on the spot): recombinations directly to a level with Ex = 0 of a multi-level atom make a photon that is reabsorbed nearby, so they do not count
        caseB = hasproperty(env, :caseB) && env.caseB && Ex == 0 && length(p.levels) > 1
        R += w * exp(-Ex / (kB * T)) * Rl; caseB || (A += w * Al)
    end
    R, A
end

"""
Photo-rate R of `photorates` alone: the same arithmetic in the same order (R is bit-identical), without the integrand of the
Milne coefficient. For callers that use only R (the star-light rates of every zone and step: the unused α was half of their cost).
"""
function photorate(p::Photo, T, env; n = 200)
    Z = p.Z(T); R = 0.0
    for (Ex, g, Eth, σf) in p.levels
        Rl = _photoR(σf, env.J, Eth / hP, range(0.0, log(1 + 60kB * max(T, env.Trad) / Eth + 30.0), length = n + 1), n)
        w = g / Z
        R += w * exp(-Ex / (kB * T)) * Rl
    end
    R
end
# function barrier: compiled for the concrete cross-section and J closures (σf is an untyped field of Photo)
function _photoR(σf::F, Jf::G, ν0, xs, n) where {F,G}
    Rl = 0.0; νa = ν0 * exp(xs[1])
    for i in 1:n
        νb = ν0 * exp(xs[i+1]); ν = 0.5 * (νa + νb); dν = νb - νa; σ = σf(ν); J = Jf(ν)
        Rl += 4π * σ * J / (hP * ν) * dν; νa = νb
    end
    Rl
end

"""
photorates through a table in log T (1–1e6 K, 0.005 dex, Catmull–Rom in ln R and ln α) when the radiation field is the same in
every zone (`env.static == true`, e.g. the CMB or none); otherwise direct. Table vs direct ≤ 3e-3 at 20 K–9e5 K (H, H⁻, Li; J = B(54 K)); the larger differences, H⁻ α above 7e4 K, are the direct quadrature's own (n = 200 vs 800: 3 %).
"""
const PHTAB = IdDict{Any,Any}(); const PHLOCK = ReentrantLock()
function photorates_t(p::Photo, T, env)
    (hasproperty(env, :static) && env.static) || return photorates(p, T, env)
    # env.tabkey: a per-zone env (extra fields for reactions) that shares the static J of `tabkey`, and so its table
    key = (p, hasproperty(env, :tabkey) ? env.tabkey : env)
    tab = get(PHTAB, key, nothing)
    if tab === nothing
        # built once, under a lock: the chemistry threads all arrive here in the first step (they used to build the table 16
        # times over, 5 s each, and to insert into the IdDict concurrently)
        tab = lock(PHLOCK) do
            get!(PHTAB, key) do
                lT = collect(0.0:0.005:6.0)
                v = [photorates(p, 10^x, env) for x in lT]
                (lT, [log(max(a[1], 1e-300)) for a in v], [log(max(a[2], 1e-300)) for a in v])
            end
        end
    end
    lT, lR, lA = tab; n = length(lT); x = clamp(log10(T), 0.0, 6.0); i = clamp(Int(floor(x / 0.005)) + 1, 2, n - 2); w = (x - lT[i]) / 0.005
    function cr(f)                                                 # Catmull–Rom
        p0, p1, p2, p3 = f[i-1], f[i], f[i+1], f[i+2]
        p1 + 0.5w * (p2 - p0 + w * (2p0 - 5p1 + 4p2 - p3 + w * (3(p1 - p2) + p3 - p0)))
    end
    R = exp(cr(lR)); A = exp(cr(lA))
    (R < 1e-290 ? 0.0 : R), A
end

# ------------------------------------------------------------------ composition
"Eliminated species: H (H nuclei), Li (Li nuclei), e (charge); He fixed at y_He, or from He nuclei when He⁺/He²⁺ are free (HECHEM)."
function full(x, net)
    d = Dict{Symbol,eltype(x)}()
    for (i, s) in enumerate(net.free); d[s] = x[i]; end
    if any(s -> SPECIES[s].el[2] > 0, net.free)
        d[:He] = yHe - sum((SPECIES[s].el[2] * d[s] for s in net.free if SPECIES[s].el[2] > 0); init = zero(eltype(x)))
    else
        get!(d, :He, yHe + zero(eltype(x)))
    end
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
    elHe::Vector{Float64}; hasHe::Bool                                     # He nuclei per free species; He eliminated (HECHEM)
    cv::Vector{Function}; E0::Vector{Float64}; isH2::Vector{Bool}
end
const COMPILED = IdDict{Any,Compiled}()
compiled(net) = get!(COMPILED, net) do
    sp = unique(vcat(net.species, [:He])); idx = Dict(s => i for (i, s) in enumerate(sp))
    fr = [idx[s] for s in net.free]
    Compiled(sp, idx, fr, [[idx[s] for s in r.R] for r in net.reactions], [[idx[s] for s in r.P] for r in net.reactions],
             get(idx, :H, 0), get(idx, :Li, 0), get(idx, :e, 0), idx[:He],
             [s == :H ? 0.0 : Float64(SPECIES[s].el[1]) for s in net.free], [Float64(SPECIES[s].el[3]) for s in net.free],
             [Float64(SPECIES[s].q) for s in net.free], [Float64(SPECIES[s].el[2]) for s in net.free], any(s -> SPECIES[s].el[2] > 0, net.free),
             [SPECIES[s].cv for s in sp], [SPECIES[s].E0 for s in sp], [s == :H2 for s in sp])
end
"Full composition as a vector over compiled(net).sp (same elimination as `full`); entries of x behind the free species are ignored."
function fullvec(x, c::Compiled)
    v = zeros(eltype(x), length(c.sp))
    for (i, j) in enumerate(c.free); v[j] = x[i]; end
    v[c.iHe] = c.hasHe ? yHe - sum(c.elHe[i] * x[i] for i in eachindex(c.free)) : yHe
    c.iH > 0 && (v[c.iH] = 1 - sum(c.elH[i] * x[i] for i in eachindex(c.free)))
    c.iLi > 0 && (v[c.iLi] = yLi - sum(c.elLi[i] * x[i] for i in eachindex(c.free)))
    c.ie > 0 && (v[c.ie] = sum(c.q[i] * x[i] for i in eachindex(c.free)))
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

"""
From this recursion depth, a failed Newton is retried with each free abundance clipped on its own (x → 0.1 x) instead of scaling the
whole correction. A trace species whose linearised update is strongly negative (H₂⁺ at n ≈ 1e16, T ≈ 3300 K, with the H₂ row far from
balance) otherwise sets λ ≈ 0.9 x/|δ| for all unknowns: x falls by 10 % per iteration, λ → 1e-60, and Newton freezes for every h.
Runs that never reach this depth are unchanged bit for bit.
"""
const CLIPDEPTH = 4
"""Non-finite Jacobian = failed attempt (halve h), and the clipped retry clips trace Li on its own species (ENV CHEMROBUST=1, default off;
off keeps the old behaviour: `\\` throws, and Li → 0 sets λ for every unknown)."""
const ROBUST = Ref(get(ENV, "CHEMROBUST", "0") != "0")
"""
CHEMTRCLIP=f (default 0, off): positivity of trace species per species instead of through the common step length. A free species
below f per H nucleus whose Newton update would make it negative has its own update cut to −0.9 of its value, and an eliminated Li
nucleus going negative scales only the Li species' updates; the common λ is then set by the major species (H, e, He) and T alone.
Without it a trace species (He⁺ ~1e-10 with HECHEM=1, neutral Li) whose linear update overshoots by ×1e4 sets λ ~ 1e-6, every
other unknown stays put, λ falls ×10 per iteration and the attempt fails (13 Nov 2026, the pre-shock zones of t_b224c8: 42
failed attempts in one step at h ≥ dt/16).
"""
const TRCLIP = Ref(parse(Float64, get(ENV, "CHEMTRCLIP", "0")))
"""
CHEMLOGX=1 (default off): Newton of `step` in w = ln x for the free species (and ln T), so that their positivity holds by
construction and the exponential dependence of trace species in fast equilibrium (H⁻, H₂⁺, He⁺, Li⁺ …) on T is linear in the
unknowns. The eliminated H, e, He keep the common step length (halved until non-negative); the eliminated Li scales only the Li
species' updates. Species below 1e-25 per H nucleus do not enter the convergence test (floor 1e-60). On failure the halving and
the CLIPDEPTH retries of the linear form follow as before.
"""
const LOGX = Ref(get(ENV, "CHEMLOGX", "0") == "1")
"""
CHEMSUBCYC=ε (default 0, off): `step` without `extra` and `Tfix` integrates by sub-cycling instead of one backward-Euler Newton
solve (the approach of ChemistryKernels.jl's sub-cycler; Anninos et al. 1997, NewA 2, 209). Per sub-step h, every free species
from its formation C and destruction rate D at the start of the sub-step, x' = (x + C h)/(1 + D h) (positive for any h; species in
fast equilibrium, D h ≫ 1, take x = C/D: H⁻, H₂⁺, He⁺, Li⁺ …), then T from u(ρ, T, x') = u − h Λ(ρ, T, x')/ρ (one variable,
cooling implicit in T). h is limited to ε of x/|dx/dt| for x_e and every free species above 1e-8 per H that is not in fast
equilibrium, |Δ ln T| ≤ max(ε, 0.02), growth ×2 per sub-step; a sub-step that leaves an eliminated abundance negative or a T
solve that fails is retried with h/4.
Motivation (13 Nov 2026, the pre-shock zones of t_b224c8 at n ~ 3e12, T 5500-7000 K, dt ~ 1.3e6 s): the single backward-Euler
step needs 160-500 Newton iterations in 13-31 halvings per call; from the start state Newton is outside the basin of the root
(positivity damping by trace species, or a spurious x_e → 0 attractor in ln x), and the CHEMSUB halving tree costs 0.5-1 s per zone.
"""
const SUBCYC = Ref(parse(Float64, get(ENV, "CHEMSUBCYC", "0"))); const NSUBCYC = Threads.Atomic{Int}(0)
"Formation C and destruction L per H nucleus and second for every species of compiled(net).sp (dx/dt = C − L), with the full vector v."
function cdrates(x, nH, T, net, env, ph)
    c = compiled(net); v = fullvec(x, c); C = zeros(length(v)); L = zeros(length(v))
    for (ir, r) in enumerate(net.reactions)
        Ri = c.R[ir]; rate = r.k(T, env) * nH^(length(Ri) - 1)
        for j in Ri; rate *= v[j]; end
        for j in Ri; L[j] += rate; end
        for j in c.P[ir]; C[j] += rate; end
    end
    for (p, R, α) in ph
        jr = c.idx[p.R]; jp = [c.idx[s] for s in p.P]
        rate = R * v[jr]; for j in jp; C[j] += rate; end; L[jr] += rate
        inv_ = α * nH^(length(jp) - 1); for j in jp; inv_ *= v[j]; end
        for j in jp; L[j] += inv_; end; C[jr] += inv_
    end
    v, C, L
end
"The sub-cycled integration of CHEMSUBCYC (see there); returns (x, T) as `step`."
function step_sub(x0, ρ, u, T0, dt, net; env, photos = Photo[], cooling = nothing, ε = SUBCYC[], maxsub = 200000)
    nH = ρ / μH(); c = compiled(net); nf = length(net.free)
    x = collect(Float64, x0); T = Float64(T0); ucur = Float64(u)
    t = 0.0; h = dt; nsub = 0; εT = max(ε, 0.02)
    while t < dt * (1 - 1e-12)
        nsub += 1; nsub > maxsub && error("ChemNet.step_sub: more than $maxsub sub-steps (nH = $nH, T0 = $T0, dt = $dt, u = $u)")
        CHEMSTATS && Threads.atomic_add!(NADV, 1)
        ph = [(p, photorates_t(p, T, env)...) for p in photos]
        v, C, L = cdrates(x, nH, T, net, env, ph)
        hh = min(h, dt - t)
        for (i, j) in enumerate(c.free)
            D = L[j] / max(v[j], 1e-300)
            (x[i] > 1e-8 && D * hh < 1) || continue
            r = abs(C[j] - L[j]) / x[i]; r > 0 && (hh = min(hh, ε / r))
        end
        if c.ie > 0 && v[c.ie] > 1e-8
            de = sum(c.q[i] * (C[c.free[i]] - L[c.free[i]]) for i in 1:nf); de != 0 && (hh = min(hh, ε * v[c.ie] / abs(de)))
        end
        hh = max(hh, 1e-12 * dt)
        xn = similar(x)
        for (i, j) in enumerate(c.free)
            xn[i] = (x[i] + C[j] * hh) / (1 + L[j] / max(v[j], 1e-300) * hh)
        end
        vn = fullvec(xn, c)
        if any(j -> j > 0 && vn[j] < 0, (c.iH, c.ie, c.iHe, c.iLi))
            hh <= 1e-12 * dt && error("ChemNet.step_sub: negative eliminated abundance at the smallest sub-step (nH = $nH, T = $T)")
            h = hh / 4; continue
        end
        dn = full(xn, net)
        g(lT) = (Tq = exp(lT); (eos(ρ, Tq, xn, net)[1] - ucur + (cooling === nothing ? zero(Tq) : hh * cooling(ρ, Tq, dn) / ρ)) / (kB * Tq / μH()))
        lT = log(T); ok = false
        for it in 1:50
            gv = g(lT); dg = ForwardDiff.derivative(g, lT)
            (isfinite(gv) && isfinite(dg) && dg != 0) || break
            δ = clamp(-gv / dg, -0.2, 0.2); lT += δ
            abs(δ) < 1e-11 && (ok = true; break)
        end
        if !ok || abs(lT - log(T)) > εT
            hh <= 1e-12 * dt && error("ChemNet.step_sub: no temperature at the smallest sub-step (nH = $nH, T = $T)")
            h = hh / (ok ? 2 : 4); continue
        end
        x = xn; T = exp(lT); ucur = Float64(eos(ρ, T, x, net)[1]); t += hh; h = 2hh
    end
    CHEMSTATS && Threads.atomic_add!(NSUBCYC, nsub)
    x, T
end
"Print the first N non-finite Jacobians of `step` (ENV CHEMNANDBG; diagnostics)."
const NANDBG = Ref(parse(Int, get(ENV, "CHEMNANDBG", "0")))
"""
`CHEMSTATS=1` (diagnostic, default off): count calls of `step`, sub-step attempts (`adv` calls: halvings included) and Jacobians
(= Newton iterations); totals printed at exit.
"""
const CHEMSTATS = get(ENV, "CHEMSTATS", "0") == "1"
"""
`CHEMSTAG=n` (default 0 = off): accept a Newton attempt that has stagnated at round-off (full step λ = 1, every abundance changed by
< 1e-6 relative and ln T by < 1e-8, the acceptance test applied after `maxit` iterations) already from iteration n on. Measured
12 Nov: in b224's network zones the step took 24 Newton iterations on average against 3.2 in c128, many running to maxit = 60 with
a trace species flipping in its last digits. Results differ from the default at the round-off level of those species.
"""
const CHEMSTAG = parse(Int, get(ENV, "CHEMSTAG", "0"))
"""
`CHEMGIVEUP=n` (default 0 = off): an attempt that at iteration n is still far from converged (damped step, or a relative change
above 1e-3) is abandoned as failed, so the step is halved at once instead of after `maxit` = 60 iterations. Measured 12 Nov in b224:
~17 000 of ~56 000 attempts failed, each after 60 iterations, ≈ all of the ~0.9 million Jacobians.
"""
const CHEMGIVEUP = parse(Int, get(ENV, "CHEMGIVEUP", "0"))
const NSTEP = Threads.Atomic{Int}(0); const NADV = Threads.Atomic{Int}(0); const NJAC = Threads.Atomic{Int}(0)
CHEMSTATS && atexit(() -> println("CHEMSTATS: ChemNet.step calls ", NSTEP[], ", sub-step attempts ", NADV[], ", Jacobians (Newton iterations) ", NJAC[]))
"Print the Newton iterations of the first N failed steps (ENV CHEMTRACE; diagnostics)."
const TRACEDBG = Ref(parse(Int, get(ENV, "CHEMTRACE", "0")))
"Diagnostics (default false): trace every attempt of `step` (set from code, e.g. PopIIIMaterial.chemreplay)."
const TRACEALL = Ref(false)

"""
    step(x0, ρ, u, T0, dt, net; env, photos = Photo[], cooling = nothing, extra = nothing, Tfix = nothing) -> (x, T)

Chemistry and temperature together over dt at fixed ρ and u (u − dt Λ/ρ with `cooling(ρ, T, d)`, d the full composition).
`extra(x)`: added to dx/dt, a function of the new abundances (a photon-conserving ray, IonRay). `Tfix`: hold T at this value instead of solving for it (tests).
`accepted(x, h)`: called after each accepted sub-step, so that a ray can count the photons used through the halvings.
"""
function step(x0, ρ, u, T0, dt, net; env, photos = Photo[], cooling = nothing, extra = nothing, Tfix = nothing, accepted = nothing, tol = 1e-10, maxit = 60)
    nH = ρ / μH(); nf = length(net.free)
    CHEMSTATS && Threads.atomic_add!(NSTEP, 1)
    SUBCYC[] > 0 && extra === nothing && Tfix === nothing && return step_sub(x0, ρ, u, T0, dt, net; env, photos, cooling)
    phT = Ref(NaN); phc = Ref{Any}(nothing)
    ucur = Ref(Float64(u))                           # u at the start of the current sub-step: what cooling took in the earlier halves stays gone
    function resid(z, x0, h)
        x = z[1:nf]; T = exp(z[nf+1])
        Tv = ForwardDiff.value(T)                    # photo-rates lagged in T (Float), cached per T value
        Tv == phT[] || (phc[] = [(p, photorates_t(p, Tv, env)...) for p in photos]; phT[] = Tv)
        ph = phc[]
        f = dxdt(x, nH, T, net, env, ph)
        extra === nothing || (f = f .+ extra(x))
        r = x .- x0 .- h .* f
        Tfix === nothing || return vcat(r, z[nf+1] - log(Tfix))
        ue = eos(ρ, T, x, net)[1] - ucur[] + (cooling === nothing ? zero(T) : h * cooling(ρ, T, full(x, net)) / ρ)
        vcat(r, ue / (kB * T / μH()))
    end
    function attempt(x0, T0, h; clip = false, trace = false)
        z = vcat(collect(Float64, x0), log(T0)); rel = Inf; dT = Inf; λ = 0.0
        trace = trace || TRACEALL[]; lim = :none
        TRACEALL[] && println("  attempt h ", h, " T0 ", T0, " clip ", clip)
        for it in 1:maxit
            R = resid(z, x0, h)
            trace && !TRACEALL[] && println("  trace it ", it, " T ", exp(z[nf+1]), " x ", z[1:nf], " R ", R)
            CHEMSTATS && Threads.atomic_add!(NJAC, 1)
            Jm = ForwardDiff.jacobian(zz -> resid(zz, x0, h), z)
            # a non-finite Jacobian is a failed attempt like a non-finite δ (halve h): `\` would throw an ArgumentError past the retry
            if ROBUST[] && !(all(isfinite, Jm) && all(isfinite, R))
                NANDBG[] > 0 && (NANDBG[] -= 1; println("ChemNet nonfinite: h ", h, " T ", exp(z[nf+1]), " x ", z[1:nf], " R ", R, " bad J cols ", [j for j in axes(Jm, 2) if !all(isfinite, Jm[:, j])], " bad J rows ", [i for i in axes(Jm, 1) if !all(isfinite, Jm[i, :])]))
                return (z, false)
            end
            δ = -(Jm \ R); all(isfinite, δ) || return (z, false)
            λ = 1.0; lim = :none
            clip && for i in 1:nf; z[i] + δ[i] < 0 && (δ[i] = -0.9 * z[i]); end
            TRCLIP[] > 0 && for i in 1:nf; z[i] < TRCLIP[] && z[i] + δ[i] < 0 && (δ[i] = -0.9 * z[i]); end
            # keep every abundance ≥ 0, incl. the eliminated H, Li, e (linear in x)
            for i in 1:nf; z[i] + λ * δ[i] < 0 && 0.9 * z[i] / -δ[i] < λ && (λ = 0.9 * z[i] / -δ[i]; lim = net.free[i]); end
            d0 = full(z[1:nf], net); d1 = full(z[1:nf] .+ δ[1:nf], net)
            for s in (:H, :Li, :e, :He)
                haskey(d0, s) || continue
                dv = d1[s] - d0[s]
                if ((ROBUST[] && clip) || TRCLIP[] > 0) && s === :Li && d0[s] + dv < 0
                    # Li is a trace element: clip only its own species (Li fully ionized by star light, Li → 0 froze every unknown)
                    cc = compiled(net); f = 0.9 * d0[s] / -dv
                    for i in 1:nf; cc.elLi[i] != 0 && (δ[i] *= f); end
                    continue
                end
                d0[s] + λ * dv < 0 && 0.9 * d0[s] / -dv < λ && (λ = 0.9 * d0[s] / -dv; lim = s)
            end
            abs(δ[nf+1]) * λ > 0.5 && (λ = 0.5 / abs(δ[nf+1]); lim = :T)
            z .+= λ .* δ
            rel = maximum(abs.(δ[1:nf]) ./ max.(abs.(z[1:nf]), 1e-20)); dT = abs(δ[nf+1])                                # 1e-20: abundance floor
            TRACEALL[] && (ir = argmax(abs.(δ[1:nf]) ./ max.(abs.(z[1:nf]), 1e-20));
                           println("    it ", it, " λ ", round(λ, sigdigits = 3), " lim ", lim, " δlnT ", round(δ[nf+1], sigdigits = 3), " T ", round(exp(z[nf+1]), digits = 1),
                                   " rel ", round(rel, sigdigits = 3), " (", net.free[ir], ") |R| ", round(maximum(abs, R), sigdigits = 3), " x_e ", round(full(z[1:nf], net)[:e], sigdigits = 4)))
            λ == 1.0 && rel < tol && dT < tol && return (z, true)
            CHEMSTAG > 0 && it >= CHEMSTAG && λ == 1.0 && rel < 1e-6 && dT < 1e-8 && return (z, true)
            CHEMGIVEUP > 0 && it == CHEMGIVEUP && !(λ == 1.0 && rel < 1e-3) && return (z, false)
        end
        # Stagnated at round-off: a trace species in fast equilibrium (H⁺ ≈ 4e-12 at n = 6.5e16, T = 2850 K) flips in its 10th digit
        # from one iteration to the next while T and all other abundances are converged to the last digit; the test above (1e-10 per
        # species) is then never met for any h, and the halving ended in "Newton failed" at depth 12. Accepted only after maxit
        # full Newton steps, so steps that converged before are unchanged bit for bit.
        λ == 1.0 && rel < 1e-6 && dT < 1e-8 && return (z, true)
        (z, false)
    end
    function attempt_log(x0, T0, h)
        w = vcat(log.(max.(collect(Float64, x0), 1e-60)), log(T0)); rel = Inf; dT = Inf; λ = 0.0
        rw(ww) = resid(vcat(exp.(ww[1:nf]), ww[nf+1:nf+1]), x0, h)
        cc = compiled(net)
        TRACEALL[] && println("  log attempt h ", h, " T0 ", T0)
        xof(ww) = exp.(ww[1:nf])
        for it in 1:maxit
            R = rw(w)
            CHEMSTATS && Threads.atomic_add!(NJAC, 1)
            Jm = ForwardDiff.jacobian(rw, w)
            (all(isfinite, Jm) && all(isfinite, R)) || return (vcat(xof(w), w[nf+1]), false)
            δ = -(Jm \ R); all(isfinite, δ) || return (vcat(xof(w), w[nf+1]), false)
            λ = 1.0; lim = :none
            d0 = full(xof(w), net)
            if haskey(d0, :Li)
                for _ in 1:60
                    full(exp.(w[1:nf] .+ δ[1:nf]), net)[:Li] >= 0.1 * d0[:Li] && break
                    for i in 1:nf; cc.elLi[i] != 0 && (δ[i] *= 0.5); end
                end
            end
            for _ in 1:60
                d1 = full(exp.(w[1:nf] .+ λ .* δ[1:nf]), net)
                bad = false
                for s in (:H, :e, :He)
                    haskey(d1, s) && d1[s] < 0.1 * d0[s] && d0[s] > 0 && (bad = true; lim = s)
                end
                bad || break
                λ *= 0.5
            end
            abs(δ[nf+1]) * λ > 0.5 && (λ = 0.5 / abs(δ[nf+1]); lim = :T)
            w .+= λ .* δ
            for i in 1:nf; w[i] = max(w[i], log(1e-60)); end
            act = [i for i in 1:nf if w[i] > log(1e-25)]
            rel = isempty(act) ? 0.0 : maximum(abs(δ[i]) for i in act); dT = abs(δ[nf+1])
            TRACEALL[] && println("    lit ", it, " λ ", round(λ, sigdigits = 3), " lim ", lim, " δlnT ", round(δ[nf+1], sigdigits = 3), " T ", round(exp(w[nf+1]), digits = 1),
                                  " max|δlnx| ", round(rel, sigdigits = 3), " |R| ", round(maximum(abs, R), sigdigits = 3), " x_e ", round(full(xof(w), net)[:e], sigdigits = 4))
            λ == 1.0 && rel < tol && dT < tol && return (vcat(xof(w), w[nf+1]), true)
            CHEMSTAG > 0 && it >= CHEMSTAG && λ == 1.0 && rel < 1e-6 && dT < 1e-8 && return (vcat(xof(w), w[nf+1]), true)
            CHEMGIVEUP > 0 && it == CHEMGIVEUP && !(λ == 1.0 && rel < 1e-3) && return (vcat(xof(w), w[nf+1]), false)
        end
        λ == 1.0 && rel < 1e-6 && dT < 1e-8 && return (vcat(xof(w), w[nf+1]), true)
        (vcat(xof(w), w[nf+1]), false)
    end
    function adv(x0, T0, h, depth)
        CHEMSTATS && Threads.atomic_add!(NADV, 1)
        z, ok = LOGX[] ? attempt_log(x0, T0, h) : attempt(x0, T0, h)
        ok && (accepted === nothing || accepted(z[1:nf], h); return (z[1:nf], exp(z[nf+1])))
        if depth >= CLIPDEPTH                        # (see the note at CLIPDEPTH)
            z, ok = attempt(x0, T0, h; clip = true)
            ok && (accepted === nothing || accepted(z[1:nf], h); return (z[1:nf], exp(z[nf+1])))
        end
        depth >= 12 && TRACEDBG[] > 0 && (TRACEDBG[] -= 1; println("ChemNet trace: h ", h, " T0 ", T0, " x0 ", x0, " ρ ", ρ, " u ", ucur[], " u(T0,x0) ", eos(ρ, T0, x0, net)[1]); println("  free ", net.free); attempt(x0, T0, h; trace = true); println("  clip:"); attempt(x0, T0, h; clip = true, trace = true))
        depth >= 12 && error("ChemNet.step: Newton failed (nH = $nH, T0 = $T0, dt = $dt, u = $u)")
        x1, T1 = adv(x0, T0, h / 2, depth + 1)
        cooling === nothing || Tfix !== nothing || (ucur[] = Float64(eos(ρ, T1, x1, net)[1]))
        adv(x1, T1, h / 2, depth + 1)
    end
    adv(x0, T0, dt, 0)
end

"""
Steady state of the network at fixed (ρ, T) (for tests): pseudo-transient continuation — backward-Euler steps at fixed
T; a step is accepted only when its Newton converges (then dt × 3), otherwise it is retried with dt / 10; ends when an
accepted step with dt ≥ 1e30 s leaves the state unchanged.
"""
function steady(x0, ρ, T, net; env, photos = Photo[], maxsteps = 2000)
    nH = ρ / μH(); nf = length(net.free)
    ph = [(p, photorates(p, T, env)...) for p in photos]
    x = collect(Float64, x0); dt = 1e-8
    function newton(xo, dt)
        x = copy(xo)
        F(xx) = xx .- xo .- dt .* dxdt(xx, nH, T, net, env, ph)
        for it in 1:60
            R = F(x); Jm = ForwardDiff.jacobian(F, x)
            δ = -(Jm \ R); all(isfinite, δ) || return (x, false)
            λ = 1.0
            for i in 1:nf; x[i] + λ * δ[i] < 0 && (λ = min(λ, 0.9 * x[i] / -δ[i])); end
            d0 = full(x, net); d1 = full(x .+ δ, net)
            for s in (:H, :Li, :e); haskey(d0, s) || continue; dv = d1[s] - d0[s]; d0[s] + λ * dv < 0 && (λ = min(λ, 0.9 * d0[s] / -dv)); end
            x .+= λ .* δ
            λ == 1.0 && maximum(abs.(δ) ./ max.(abs.(x), 1e-30)) < 1e-13 && return (x, true)
        end
        (x, false)
    end
    for k in 1:maxsteps
        xn, ok = newton(x, dt)
        if ok
            done = dt >= 1e30 && maximum(abs.(xn .- x) ./ max.(abs.(xn), 1e-30)) < 1e-10
            x = xn; dt *= 3
            done && return x
        else
            dt /= 10; dt < 1e-30 && error("ChemNet.steady: no convergence (nH = $nH, T = $T)")
        end
    end
    x
end

end # module
