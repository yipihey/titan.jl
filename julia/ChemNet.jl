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
    H2PF && return u_H2_g2(T)
    u_H2_rhfc02(T)
end
function u_H2_rhfc02(T)
    x = 6100 / T
    ForwardDiff.value(x) < 100 && return 5 / 2 + x^2 * exp(x) / expm1(x)^2
    em = exp(-x)
    5 / 2 + x^2 * em / (1 - em)^2
end
const χH = 13.598eV; const EA_H = 0.754eV; const D0_H2 = 4.478eV; const D0_H2p = 2.65eV
"""
`H2PF=1` (default off; off = bit-identical): the internal partition function of H₂ as the sum over all 301 bound rovibrational levels of
the ground electronic state, ExoMol RACPPK (Roueff et al. 2019, A&A 630, A58; `POPIII_LIT/H2lines/1H2__RACPPK.states`: energies and
g_tot with nuclear-spin weights 1 (para) / 3 (ortho), the convention of Zrot_H2), instead of the rigid rotor × harmonic oscillator
Zrot_H2 Zvib; and the H₂ thermal energy consistently from it, u_H2 = 3/2 + ⟨E_int⟩/kT (instead of RHFC02 A9). Enters every equilibrium
constant with H₂ (dissociation) and the EOS. G2 audit (docs/POPIII_G2_MATERIAL.md §7): the rigid-rotor/harmonic form keeps H₂ bound
too long at 3000–5000 K (SCVH95, whose H₂ partition function is a level sum, dissociates earlier; network vs table 0.6–0.8 kT per particle
in u, 4 % in p at log ρ −5.5, log T 3.58). Applied above 1500 K only, blended in from 800 K (`h2pf_w`): below, the cold gas keeps RHFC02.
"""
const H2PF = get(ENV, "H2PF", "0") == "1"
const H2LEV = H2PF ? let f = joinpath(get(ENV, "POPIII_LIT", joinpath(homedir(), "Projects", "popiii_lit")), "H2lines", "1H2__RACPPK.states")
    rows = [split(l) for l in readlines(f) if !isempty(strip(l))]
    (E = [parse(Float64, r[2]) * hP * c for r in rows], g = [parse(Float64, r[3]) for r in rows])
end : (E = Float64[], g = Float64[])
"(H2PF) weight of the level sum: 0 below 800 K (cold gas unchanged: the ortho/para ratio there is not in equilibrium, and the run's
convention 5/2 kT stays), 1 above 1500 K, C¹ smoothstep in ln T between."
h2pf_w(T) = (x = clamp((log(T) - log(800.0)) / (log(1500.0) - log(800.0)), zero(T), one(T)); x * x * (3 - 2x))
"(H2PF) Z_int of H₂: rigid rotor × harmonic oscillator below 800 K, the RACPPK level sum above 1500 K, geometric blend between."
Zint_H2_g2(T) = (w = h2pf_w(T); w == 0 ? Zrot_H2(T) * Zvib(T, θvib_H2) : exp((1 - w) * log(Zrot_H2(T) * Zvib(T, θvib_H2)) + w * log(Zint_H2_pf(T))))
"(H2PF) thermal energy of H₂ per molecule / kT: RHFC02 A9 below 800 K, the level sum above 1500 K, blended between."
u_H2_g2(T) = (w = h2pf_w(T); w == 0 ? u_H2_rhfc02(T) : (1 - w) * u_H2_rhfc02(T) + w * u_H2_pf(T))
"H₂ internal partition function from the RACPPK levels (H2PF)."
Zint_H2_pf(T) = (z = zero(T); for i in eachindex(H2LEV.E); z += H2LEV.g[i] * exp(-H2LEV.E[i] / (kB * T)); end; z)
"H₂ thermal energy per molecule in units of kT from the RACPPK levels: 3/2 + ⟨E_int⟩/kT (H2PF)."
function u_H2_pf(T)
    z = zero(T); e = zero(T)
    for i in eachindex(H2LEV.E); t = H2LEV.g[i] * exp(-H2LEV.E[i] / (kB * T)); z += t; e += t * H2LEV.E[i]; end
    3 / 2 + e / (z * kB * T)
end
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
    :H2  => Species(:H2, 2mH, 0, (2, 0, 0), -(2χH + D0_H2), T -> H2PF ? Zint_H2_g2(T) : Zrot_H2(T) * Zvib(T, θvib_H2), T -> NaN),   # c_v via u_H2
    :H2p => Species(:H2p, 2mH - me, 1, (2, 0, 0), -(χH + D0_H2p), Zint_H2p, T -> 2.5),
    :e   => Species(:e, me, -1, (0, 0, 0), 0.0, T -> 2 + zero(T), T -> 1.5),
    :He  => Species(:He, 4.0026 * 1.66053907e-24, 0, (0, 1, 0), 0.0, HECHEM ? ZHe_el : (T -> 1 + zero(T)), T -> 1.5),
    :Hep  => Species(:Hep, 4.0026 * 1.66053907e-24 - me, 1, (0, 1, 0), χHe1, ZHep_el, T -> 1.5),
    :Hepp => Species(:Hepp, 4.0026 * 1.66053907e-24 - 2me, 2, (0, 1, 0), χHe1 + χHe2, T -> 1 + zero(T), T -> 1.5),
    :Li  => Species(:Li, 7.016 * 1.66053907e-24, 0, (0, 0, 1), -χLi, T -> 4 * ZLi_el(T), T -> 1.5),     # × nuclear spin 4
    :Lip => Species(:Lip, 7.016 * 1.66053907e-24 - me, 1, (0, 0, 1), 0.0, T -> 4 + zero(T), T -> 1.5),
)

"""
(FIX-RACE) Lazily filled caches read from threaded zone loops (MATTHREADS sweeps, chemistry threads): copy-on-write. Readers look up
the current dictionary, which is never mutated after it is published; a miss builds the value under the lock, inserts it into a copy and
publishes the copy. (A plain IdDict read while another task inserts — rehashes — can see undefined slots: UndefRefError, z100s50b.)
The values are the ones the plain dictionaries held: results unchanged.
"""
mutable struct CowCache{V}
    @atomic d::IdDict{Any,V}
    lk::ReentrantLock
end
CowCache{V}() where {V} = CowCache{V}(IdDict{Any,V}(), ReentrantLock())
@inline cowget(c::CowCache, key, default) = get((@atomic :acquire c.d), key, default)
function cowget!(f, c::CowCache{V}, key) where {V}
    v = get((@atomic :acquire c.d), key, nothing); v === nothing || return v
    lock(c.lk) do
        d0 = @atomic :acquire c.d
        v2 = get(d0, key, nothing); v2 === nothing || return v2
        val = f(); d1 = copy(d0); d1[key] = val
        @atomic :release c.d = d1
        val
    end
end
Ztr(s::Species, T) = (2π * s.m * kB * T / hP^2)^1.5
@inline function _zint(s::Species, T)
    ET = typeof(one(T)); n = s.name
    n === :H ? 2 * ZH_el(T) :
    (n === :Hp || n === :Hm || n === :e) ? 2 + zero(T) :
    n === :H2 ? (H2PF ? Zint_H2_g2(T) : Zrot_H2(T) * Zvib(T, θvib_H2)) :
    n === :H2p ? Zint_H2p(T) :
    n === :He ? (HECHEM ? ZHe_el(T) : 1 + zero(T)) :
    n === :Hep ? ZHep_el(T) :
    n === :Hepp ? 1 + zero(T) :
    n === :Li ? 4 * ZLi_el(T) :
    n === :Lip ? 4 + zero(T) : (ET(s.Zint(T)))::ET
end
"Equilibrium constant Π_P n / Π_R n [cm^{-3 Δν}] of R → P (vectors of species names, repeats allowed)."
Keq(R, P, T) = exp(lnKeq(R, P, T))
"ln Keq — use this where K itself under- or overflows (ionisation equilibria below a few hundred K)."
function lnKeq(R, P, T)
    lnK = zero(T)
    for p in P; s = SPECIES[p]; lnK += log(Ztr(s, T) * _zint(s, T)) - s.E0 / (kB * T); end
    for r in R; s = SPECIES[r]; lnK -= log(Ztr(s, T) * _zint(s, T)) - s.E0 / (kB * T); end
    lnK
end
"""
PERF512: Float64 T: the species terms ln(Z_tr Z_int) − E0/kT (partition functions summed over up to 61 levels each) are kept per thread for the
last T, so that the reverse rates of all reactions of a zone evaluated at one T compute each species' term once. The same term, added in
the same order: bit-identical to the generic method. (No yield point lies between reading the thread's memo and using it, so a task cannot
change threads in between.)
"""
mutable struct LnKMemo
    T::Float64
    v::Dict{Symbol,Float64}
end
const LNKMEMOON = get(ENV, "LNKMEMO", "1") != "0"     # (FIX-RACE) LNKMEMO=0: the direct sum (bit-identical)
const LNKMEMO = [LnKMemo(NaN, Dict{Symbol,Float64}()) for _ in 1:Threads.maxthreadid()]
__init__() = (resize!(LNKMEMO, Threads.maxthreadid()); for i in eachindex(LNKMEMO); isassigned(LNKMEMO, i) || (LNKMEMO[i] = LnKMemo(NaN, Dict{Symbol,Float64}())); end)   # (again when a precompiled TitanStack image is loaded: the build process has one thread)
@inline _spterm(s::Species, T) = log(Ztr(s, T) * _zint(s, T)) - s.E0 / (kB * T)
function _spterm_memo(m::LnKMemo, p::Symbol, T::Float64)
    t = get(m.v, p, NaN)
    isnan(t) || return t
    t = _spterm(SPECIES[p], T)::Float64; m.v[p] = t
    t
end
function lnKeq(R, P, T::Float64)
    tid = Threads.threadid()
    (LNKMEMOON && tid <= length(LNKMEMO)) || return _lnKeq_direct(R, P, T)
    m = LNKMEMO[tid]
    m.T === T || (m.T = T; empty!(m.v))
    lnK = 0.0
    for p in P; lnK += _spterm_memo(m, p, T); end
    for r in R; lnK -= _spterm_memo(m, r, T); end
    lnK
end
function _lnKeq_direct(R, P, T::Float64)
    lnK = zero(T)
    for p in P; s = SPECIES[p]; lnK += log(Ztr(s, T) * _zint(s, T)) - s.E0 / (kB * T); end
    for r in R; s = SPECIES[r]; lnK -= log(Ztr(s, T) * _zint(s, T)) - s.E0 / (kB * T); end
    lnK
end

# ------------------------------------------------------------------ reactions
struct Reaction
    name::String
    R::Vector{Symbol}; P::Vector{Symbol}
    k::Function           # (T, env) -> rate coefficient [cm^{3(n_R−1)} s⁻¹]
end
function hm_cs_hook end
struct VernerCS <: Function
    is_hm::Bool
    eth::Float64
    p::NTuple{7,Float64}
end
VernerCS(eth::Real, p::NTuple{7,Float64}) = VernerCS(false, Float64(eth), p)
VernerCS() = VernerCS(true, 0.0, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
@inline function (v::VernerCS)(ν)
    v.is_hm && return hm_cs_hook(ν)
    E = 6.62607015e-27 * ν / 1.602176634e-12; E < v.eth && return 0.0
    E0, σ0, ya, P, yw, y0, y1 = v.p
    x = E / E0 - y0; z = sqrt(x^2 + y1^2)
    σ0 * ((x - 1)^2 + yw^2) * z^(0.5P - 5.5) * (1 + sqrt(z / ya))^(-P) * 1e-18
end
struct KramersCS <: Function
    n::Float64
    div4::Bool
    eth::Float64
end
@inline (k::KramersCS)(ν) = (k.div4 ? 7.91e-18 * k.n / 4 : 7.91e-18 * k.n) * (k.eth / (6.62607015e-27 * ν))^3

struct PhotoFast{F1,F2}
    Ex::Vector{Float64}; g::Vector{Float64}; Eth::Vector{Float64}
    s1::F1; srest::Vector{F2}
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
    Ex::Vector{Float64}; g::Vector{Float64}; Eth::Vector{Float64}
    pf::Any
    function Photo(name, R, P, levels, Z, milne)
        lvs = Vector{Tuple{Float64,Float64,Float64,Function}}(levels)
        Ex = [lv[1] for lv in lvs]; g = [lv[2] for lv in lvs]; Eth = [lv[3] for lv in lvs]
        s1 = lvs[1][4]
        srest = if length(lvs) == 1 && s1 isa VernerCS
            KramersCS[]
        elseif length(lvs) > 1 && all(i -> lvs[i][4] isa typeof(lvs[2][4]), 2:length(lvs))
            [lvs[i][4]::typeof(lvs[2][4]) for i in 2:length(lvs)]
        else
            Function[lvs[i][4] for i in 2:length(lvs)]
        end
        new(String(name), Symbol(R), Vector{Symbol}(P), lvs, Z, Bool(milne), Ex, g, Eth, PhotoFast(Ex, g, Eth, s1, srest))
    end
end
"Single-level process with cross section σ(ν) above Eth."
Photo(name, R, P, σ::Function, Eth::Real, milne::Bool) = Photo(name, R, P, [(0.0, 1.0, Float64(Eth), σ)], T -> 1.0, milne)
@inline function _photoZ(p::Photo, T)
    ET = typeof(one(T)); f = p.Z
    f === ZH_el ? ZH_el(T) : f === ZHe_el ? ZHe_el(T) : f === ZHep_el ? ZHep_el(T) : f === ZLi_el ? ZLi_el(T) :
    (ET(f(T)))::ET
end
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
    rk = r.k; rR = r.R; rP = r.P
    lnkr(T, env) = log(rk(T, env)) - lnKeq(rR, rP, T)
    function k(T, env)
        kf = rk(T, env)
        dmu_shift = (FREE_EOS_ON[] && hasproperty(env, :dmu) && env.dmu !== nothing) ? _dmu_rxn(rR, rP, env.dmu) : zero(T)
        T >= Tmin && kf > 1e-290 && return exp(log(kf) - lnKeq(rR, rP, T) + dmu_shift)
        iszero(kf) || kf > 0 || return zero(T)                 # a fit that is ≤ 0 (clipped) stays 0
        if isnan(Tu[])
            t = max(Float64(Tmin), ForwardDiff.value(T)); while rk(t, env) <= 1e-290 && t < 1e6; t *= 1.05; end
            Tu[] = t
        end
        Tu[] >= 1e6 && return zero(T)
        # the power law fixed at T_u: computed once (it was recomputed on every call; same values, bit-identical). The forward
        # fits of all detailed-balance pairs ignore env (the three-level atom's env-dependent factor wraps a forward reaction only).
        if isnan(Ta[])
            a = lnkr(Tu[], env); Ts[] = (lnkr(2Tu[], env) - a) / log(2.0); Ta[] = a
        end
        exp(Ta[] + Ts[] * log(T / Tu[]) + dmu_shift)
    end
end

"The n + 1 nodes 0 … xe of the log-frequency quadrature of the photo-rates (a Dual end point, T_c a Dual: the nodes carry its derivative)."
@inline _xgrid(xe, n) = range(0.0, xe, length = n + 1)
@inline _xgrid(xe::ForwardDiff.Dual, n) = [xe * (i / n) for i in 0:n]
"Photo-rate [s⁻¹] 4π∫σJ/(hν)dν and the Milne inverse coefficient α for A + e → X + γ (stimulated emission included)."
function photorates(p::Photo, T, env; n = 200)
    pf = p.pf
    pf isa PhotoFast{VernerCS,KramersCS} ? _photorates_fast(p, pf, T, env, n) :
    pf isa PhotoFast{KramersCS,KramersCS} ? _photorates_fast(p, pf, T, env, n) :
    _photorates_fast(p, pf, T, env, n)
end
function _photorates_fast(p::Photo, pf::PhotoFast{F1,F2}, T, env::E, n::Int) where {F1,F2,E}
    # α = Σ_levels (n_l/Π n_P)_eq × 4π∫σ_l/(hν)(2hν³/c² + J) e^{−hν/kT} dν, (n_l/Π n_P)_eq = g_l e^{−E_x/kT}/(Z K(X → P));
    # the Boltzmann factors and 1/K are combined in one exponent (each alone under/overflows below a few hundred K)
    ET = typeof(one(T))
    lnK = lnKeq((p.R,), p.P, T)::ET; Z = _photoZ(p, T); R = zero(ET); A = zero(ET)
    caseB0 = hasproperty(env, :caseB) && env.caseB && length(pf.Ex) > 1
    @inbounds for l in eachindex(pf.Ex)
        Ex = pf.Ex[l]; g = pf.g[l]; Eth = pf.Eth[l]
        ν0 = Eth / hP; xs = _xgrid(log(1 + 60kB * max(T, env.Trad) / Eth + 30.0), n)
        w = g / Z
        Rl, Al = l == 1 ? _photoRA_fast(pf.s1, env.J, ν0, xs, n, Ex, T, lnK) : _photoRA_fast(pf.srest[l-1], env.J, ν0, xs, n, Ex, T, lnK)
        # env.caseB (on the spot): recombinations directly to a level with Ex = 0 of a multi-level atom make a photon that is reabsorbed nearby, so they do not count
        caseB = caseB0 && Ex == 0
        R += w * exp(-Ex / (kB * T)) * Rl; caseB || (A += w * Al)
    end
    R, A
end

# (PERF-RED) the quadrature of `photorates` behind a function barrier: σf (an untyped field of Photo) and env.J are concrete here, so the
# loop no longer boxes every intermediate (it allocated ~150 GB per step in the LOCALJ run, whose rates fall back to this direct integral
# where T_c is outside photorates_lj's tables). The same arithmetic in the same order: bit-identical.
function _photoRA(σf::F, Jf::G, νs, n, Ex, T, lnK) where {F,G}
    Rl = 0.0; Al = 0.0
    for i in 1:n
        ν = 0.5 * (νs[i] + νs[i+1]); dν = νs[i+1] - νs[i]; σ = σf(ν); J = Jf(ν)
        Rl += 4π * σ * J / (hP * ν) * dν
        Al += 4π * σ / (hP * ν) * (2hP * ν^3 / c^2 + J) * exp(-(Ex + hP * ν) / (kB * T) - lnK) * dν
    end
    Rl, Al
end
function _photoRA_fast(σf::F, Jf::G, ν0::Float64, xs, n::Int, Ex::Float64, T, lnK) where {F,G}
    νa = ν0 * exp(xs[1])
    ET = typeof(one(νa) * one(T) * one(lnK))
    Rl = zero(ET); Al = zero(ET)
    for i in 1:n
        νb = ν0 * exp(xs[i+1]); ν = 0.5 * (νa + νb); dν = νb - νa; σ = σf(ν); J = Jf(ν)
        Rl += 4π * σ * J / (hP * ν) * dν
        Al += 4π * σ / (hP * ν) * (2hP * ν^3 / c^2 + J) * exp(-(Ex + hP * ν) / (kB * T) - lnK) * dν
        νa = νb
    end
    Rl, Al
end
"""
Photo-rate R of `photorates` alone: the same arithmetic in the same order (R is bit-identical), without the integrand of the
Milne coefficient. For callers that use only R (the star-light rates of every zone and step: the unused α was half of their cost).
"""
function photorate(p::Photo, T, env; n = 200)
    pf = p.pf
    pf isa PhotoFast{VernerCS,KramersCS} ? _photorate_fast(p, pf, T, env, n) :
    pf isa PhotoFast{KramersCS,KramersCS} ? _photorate_fast(p, pf, T, env, n) :
    _photorate_fast(p, pf, T, env, n)
end
function _photorate_fast(p::Photo, pf::PhotoFast{F1,F2}, T, env::E, n::Int) where {F1,F2,E}
    ET = typeof(one(T))
    Z = _photoZ(p, T); R = zero(ET)
    @inbounds for l in eachindex(pf.Ex)
        Ex = pf.Ex[l]; g = pf.g[l]; Eth = pf.Eth[l]
        xs = _xgrid(log(1 + 60kB * max(T, env.Trad) / Eth + 30.0), n)
        Rl = l == 1 ? _photoR(pf.s1, env.J, Eth / hP, xs, n) : _photoR(pf.srest[l-1], env.J, Eth / hP, xs, n)
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
mutable struct PhTab
    const p::Photo
    const tabenv::Any
    const lT::Vector{Float64}
    const lR::Vector{Float64}
    const lA::Vector{Float64}
end
const PHTAB = CowCache{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}(); const PHLOCK = PHTAB.lk
const PHTAB_LIST = Ref(PhTab[])
@inline function _get_phtab(p::Photo, tabenv, env)
    list = PHTAB_LIST[]
    @inbounds for i in eachindex(list)
        e = list[i]
        (e.p === p && e.tabenv === tabenv) && return e
    end
    lock(PHLOCK) do
        for e in PHTAB_LIST[]
            (e.p === p && e.tabenv === tabenv) && return e
        end
        lT = collect(0.0:0.005:6.0)
        v = [photorates(p, 10^x, env) for x in lT]
        tb = PhTab(p, tabenv, lT, [log(max(a[1], 1e-300)) for a in v], [log(max(a[2], 1e-300)) for a in v])
        PHTAB_LIST[] = vcat(PHTAB_LIST[], [tb])
        tb
    end
end
@inline function _ph_et(T, env)
    ET = typeof(one(T))
    (HM88 && !FREE_EOS_ON[] && hasproperty(env, :hm)) && (ET = promote_type(ET, typeof(env.hm[1])))
    (FREE_EOS_ON[] && hasproperty(env, :dmu) && env.dmu !== nothing) && (ET = promote_type(ET, typeof(env.dmu[1])))
    ET
end
"(HM88 / FREE_EOS) the Milne inverse of photoionization × non-ideal activity factor; otherwise `_photorates_t`."
function photorates_t(p::Photo, T, env)
    ET = _ph_et(T, env); r = _photorates_t(p, T, env)
    a = if FREE_EOS_ON[] && hasproperty(env, :dmu) && env.dmu !== nothing
        r[2] * exp(_dmu_rxn((p.R,), p.P, env.dmu))
    elseif HM88 && !FREE_EOS_ON[] && p.R === :H && hasproperty(env, :hm)
        r[2] * hm_phi(T, env)
    else
        r[2]
    end
    (ET(r[1])::ET, ET(a)::ET)
end
function _photorates_t(p::Photo, T, env)
    hasproperty(env, :mgq) && return photorates_mgq(p, T, env)      # MGROUPS: the static field plus the group fields (photorates_mgq)
    ET = typeof(one(T))
    if hasproperty(env, :ljW)                          # J = J_static + W B(T_c): the table of the diluted blackbody plus the static field's rates (the spontaneous α once)
        if hasproperty(env, :ljs)
            r1, r2, a0 = photorates_lj_a0(p, T, env)
            s_ = _photorates_t(p, T, env.ljs)
            return (ET(r1 + s_[1])::ET, ET(r2 + s_[2] - a0)::ET)
        end
        r = photorates_lj(p, T, env)
        return (ET(r[1])::ET, ET(r[2])::ET)
    end
    if !(hasproperty(env, :static) && env.static)
        r0 = photorates(p, T, env)
        return (ET(r0[1])::ET, ET(r0[2])::ET)
    end
    # env.tabkey: a per-zone env (extra fields for reactions) that shares the static J of `tabkey`, and so its table
    tabenv = hasproperty(env, :tabkey) ? env.tabkey : env
    tab = _get_phtab(p, tabenv, env)
    lT = tab.lT; lR = tab.lR; lA = tab.lA; n = length(lT); x = clamp(log10(T), 0.0, 6.0); i = clamp(Int(floor((_ljv(x)) / 0.005)) + 1, 2, n - 2); w = (x - lT[i]) / 0.005
    R = exp(_cr4(lR[i-1], lR[i], lR[i+1], lR[i+2], w)); A = exp(_cr4(lA[i-1], lA[i], lA[i+1], lA[i+2], w))
    (ET(R < 1e-290 ? zero(R) : R)::ET, ET(A)::ET)
end

"""
Photo-rates in a diluted blackbody J = W B(T_c) from tables (PopIIIMaterial LOCALJ with LOCALJTAB=1: the env carries `ljW` = W and
`Trad` = T_c). Exact decompositions of `photorates`, tabulated once per process:
- R = W Σ_l (g_l/Z(T)) e^{−E_x,l/kT} R_l(T_c), R_l = 4π∫σ_l B(T_c)/(hν) dν per level: 1-D tables of ln R_l in log T_c (1.5–6.5, 0.005 dex);
- α = α₀(T) + W α_st(T, T_c) (α is linear in J): ln α₀ in log T (0.5–6.5, 0.005 dex, J = 0) and ln(α_st/α₀) on a (log T, log T_c)
  grid (log T 1.0–6.5, log T_c 1.5–6.5, 0.05 dex); the ratio matters only where it is not small, and there it is smooth. The α parts
  are integrated on the gas-T grid with n = 1000, R_l on `photorates`' T_c grid with n = 1000 (both more accurate than the direct path).
Catmull–Rom in each direction (as the static tables). Outside the ranges, or with `caseB`, the direct integrals.
"""
mutable struct LJTab
    const p::Photo
    const lR::Matrix{Float64}
    const lA0::Vector{Float64}
    const lr::Matrix{Float64}
end
const LJTAB = CowCache{LJTab}(); const LJLOCK = LJTAB.lk
const LJTAB_LIST = Ref(LJTab[])
# quadrature points of the table integrals (LJTABN, default 1000). Against n = 4000 integrals at n = 1000: Li ≤ 4e-4, H⁻ ≤ 5e-3 (flow 1e-5),
# He I ≤ 2e-3 (thick) / 1e-2 (flow α), He II ≤ 1e-2 / 5e-2; the direct path (n = 200) 5–30× larger. Build time ∝ n.
const LJTABN = parse(Int, get(ENV, "LJTABN", "1000"))
const LJ_T = (0.5, 0.005, 1201); const LJ_C = (1.5, 0.005, 1001); const LJ_2T = (1.0, 0.05, 111); const LJ_2C = (1.5, 0.05, 101)
_ljv(x::ForwardDiff.Dual) = _ljv(ForwardDiff.value(x)); _ljv(x::Real) = Float64(x)
_bnu_lj(ν, T) = (x = hP * ν / (kB * T); x > 700 ? 0.0 : 2hP * ν^3 / c^2 / expm1(x))
@inline function lj_tables(p::Photo)
    list = LJTAB_LIST[]
    @inbounds for i in eachindex(list)
        e = list[i]
        e.p === p && return e
    end
    lock(LJLOCK) do
        for e in LJTAB_LIST[]
            e.p === p && return e
        end
        z0 = (Trad = 0.0, J = ν -> 0.0)
        lc = [LJ_C[1] + LJ_C[2] * (i - 1) for i in 1:LJ_C[3]]
        # R_l with n = 1000 on the T_c grid (n = 200, as `photorates`, is under-resolved by up to 7 % (He I) and 30 % (He II levels))
        lvs = p.levels; lR = zeros(length(lvs), length(lc))
        Threads.@threads for ij in 1:length(lR)
            l, j = Tuple(CartesianIndices(lR)[ij]); (Ex, g, Eth, σf) = lvs[l]; Tc = 10^lc[j]; Jb = ν -> _bnu_lj(ν, Tc); nR = LJTABN
            lR[l, j] = log(max(_photoR(σf, Jb, Eth / hP, range(0.0, log(1 + 60kB * Tc / Eth + 30.0), length = nR + 1), nR), 1e-300))
        end
        # the Milne coefficients on a frequency grid set by the gas T (Trad = 0: the integrand lives within kT of each threshold; with
        # the grid of max(T, T_c) and n = 200, as `photorates` with this J, it is under-resolved by up to 30–60 % where T ≪ T_c), n = 1000
        lt = [LJ_T[1] + LJ_T[2] * (i - 1) for i in 1:LJ_T[3]]; lA0 = zeros(length(lt))
        Threads.@threads for i in eachindex(lt); lA0[i] = log(max(photorates(p, 10^lt[i], z0; n = LJTABN)[2], 1e-300)); end
        l2t = [LJ_2T[1] + LJ_2T[2] * (i - 1) for i in 1:LJ_2T[3]]; l2c = [LJ_2C[1] + LJ_2C[2] * (j - 1) for j in 1:LJ_2C[3]]
        lr = zeros(length(l2t), length(l2c)); a0s = [photorates(p, 10^x, z0; n = LJTABN)[2] for x in l2t]
        Threads.@threads for ij in 1:length(lr)
            i, j = Tuple(CartesianIndices(lr)[ij]); T = 10^l2t[i]; Tc = 10^l2c[j]
            a1 = photorates(p, T, (Trad = 0.0, J = ν -> _bnu_lj(ν, Tc)); n = LJTABN)[2]
            lr[i, j] = log(max(a1 / max(a0s[i], 1e-300) - 1, 1e-300))
        end
        tb = LJTab(p, lR, lA0, lr)
        LJTAB_LIST[] = vcat(LJTAB_LIST[], [tb])
        tb
    end
end
@inline function _cr4(p0, p1, p2, p3, w)
    y = p1 + 0.5w * (p2 - p0 + w * (2p0 - 5p1 + 4p2 - p3 + w * (3(p1 - p2) + p3 - p0)))
    min(p0, p1, p2, p3) <= -30.0 ? clamp(y, min(p1, p2), max(p1, p2)) : y
end
@inline _ljidx(x, g) = (i = clamp(Int(floor((_ljv(x) - g[1]) / g[2])) + 1, 2, g[3] - 2); (i, (x - (g[1] + g[2] * (i - 1))) / g[2]))
function photorates_lj_a0(p::Photo, T, env)
    Tc = env.Trad; W = Float64(env.ljW); lT = log10(T); lc = log10(Tc); lTv = _ljv(lT)
    okT = !(hasproperty(env, :caseB) && env.caseB) && LJ_T[1] + 2LJ_T[2] <= lTv <= LJ_T[1] + LJ_T[2] * (LJ_T[3] - 3)
    if okT && W == 0.0
        tb = lj_tables(p)
        it, wt = _ljidx(lT, LJ_T); a0 = exp(_cr4(tb.lA0[it-1], tb.lA0[it], tb.lA0[it+1], tb.lA0[it+2], wt))
        return zero(T), a0, a0
    end
    ok = okT && LJ_C[1] + 2LJ_C[2] <= lc <= LJ_C[1] + LJ_C[2] * (LJ_C[3] - 3) &&
         LJ_2C[1] + 2LJ_2C[2] <= lc <= LJ_2C[1] + LJ_2C[2] * (LJ_2C[3] - 3) &&
         LJ_2T[1] + 2LJ_2T[2] <= lTv <= LJ_2T[1] + LJ_2T[2] * (LJ_2T[3] - 3)
    if !ok
        r = photorates(p, T, env)
        a0 = if okT
            tb = lj_tables(p)
            it, wt = _ljidx(lT, LJ_T); exp(_cr4(tb.lA0[it-1], tb.lA0[it], tb.lA0[it+1], tb.lA0[it+2], wt))
        else
            photorates(p, T, merge(env, (; ljW = 0.0)))[2]
        end
        return r[1], r[2], a0
    end
    tb = lj_tables(p)
    ic, wc = _ljidx(lc, LJ_C); Z = _photoZ(p, T); R = zero(T)
    Exs = p.Ex; gs = p.g
    @inbounds for l in eachindex(Exs)
        Ex = Exs[l]; g = gs[l]
        Rl = exp(_cr4(tb.lR[l, ic-1], tb.lR[l, ic], tb.lR[l, ic+1], tb.lR[l, ic+2], wc))
        R += g / Z * exp(-Ex / (kB * T)) * Rl
    end
    it, wt = _ljidx(lT, LJ_T); a0 = exp(_cr4(tb.lA0[it-1], tb.lA0[it], tb.lA0[it+1], tb.lA0[it+2], wt))
    jt, vt = _ljidx(lT, LJ_2T); jc, vc = _ljidx(lc, LJ_2C)
    row(j) = _cr4(tb.lr[jt-1, j], tb.lr[jt, j], tb.lr[jt+1, j], tb.lr[jt+2, j], vt)
    r = exp(_cr4(row(jc - 1), row(jc), row(jc + 1), row(jc + 2), vc))
    R = W * R
    (R < 1e-290 ? zero(R) : R), a0 * (1 + W * r), a0
end
function photorates_lj(p::Photo, T, env)
    r1, r2, _ = photorates_lj_a0(p, T, env)
    r1, r2
end

# ------------------------------------------------------------------ group-resolved photo-rates (MGROUPS stage 2c)
"""
Photo-rates from multigroup radiation fields. In group g (ν_g … ν_g+1) the field is c u_g = 4π∫_g J_ν dν [erg cm⁻² s⁻¹] with the in-group
shape of B_ν(T_c); the rate of level l is then c u_g Q_{l,g}(T_c),

  Q_{l,g}(T_c) = ∫_g σ_l B_ν(T_c)/(hν) dν / ∫_g B_ν(T_c) dν   [cm² erg⁻¹],

tabulated here as ln Q on log₁₀T_c = 1 … 6.5 (0.01 dex; linear interpolation, clamped). Both integrals carry the factor e^{hν_g/kT_c}
(the Planck function relative to the group's lower edge), so the ratio is finite where the group lies far in the Wien tail; log-spaced
midpoint rule (LJTABN points) from max(ν_g, ν_threshold) to min(ν_g+1, ν_g + 60 kT_c/h). A group entirely below a level's threshold has
no entry (−Inf). Identity (checked in mb/p2c/gtab.jl): Σ_g b_g(T_c) 4σT_c⁴ Q_{l,g} = 4π∫σ_l B_ν/(hν) dν, the rate of `lj_tables` at W = 1.
"""
const LJG_C = (1.0, 0.01, 551)
const LJGTAB = CowCache{Any}()
function _gquad(wt::F, f::H, a, b, n) where {F,H}
    s = 0.0; la = log(a); h = (log(b) - la) / n; νa = a
    for i in 1:n
        νb = exp(la + i * h); ν = 0.5 * (νa + νb); s += wt(ν) * f(ν) * (νb - νa); νa = νb
    end
    s
end
function lj_gtables(p::Photo, ehz::Vector{Float64})
    key = (p, ehz)
    t = cowget(LJGTAB, key, nothing); t === nothing || return t
    cowget!(LJGTAB, key) do
        G = length(ehz) + 1; lvs = p.levels; nc = LJG_C[3]; n = LJTABN
        lQ = fill(-Inf, G, length(lvs), nc)
        Threads.@threads for j in 1:nc
            kt = kB * 10^(LJG_C[1] + LJG_C[2] * (j - 1)) / hP             # kT_c/h [Hz]
            for g in 1:G
                ν1 = g == 1 ? 0.0 : ehz[g-1]; ν2 = g == G ? Inf : ehz[g]
                wt = ν -> 2hP * ν^3 / c^2 * exp(-(ν - ν1) / kt) / (-expm1(-ν / kt))
                lo = max(ν1, 1e-6 * kt); hi = min(ν2, ν1 + 60 * kt)
                den = _gquad(wt, ν -> 1.0, lo, hi, n)
                for (l, lv) in enumerate(lvs)
                    νt = lv[3] / hP; σf = lv[4]
                    νt < ν2 || continue
                    νa = max(lo, νt)
                    num = νa < hi ? _gquad(wt, ν -> σf(ν) / (hP * ν), νa, hi, n) : 0.0
                    lQ[g, l, j] = log(max(num / den, 1e-300))
                end
            end
        end
        lQ
    end
end
"Σ_levels (g_l/Z) e^{−E_x/kT} Σ_g cu_g Q_{l,g}(T_c) [s⁻¹]: the photo-rate of process p in the group fields cu (see `lj_gtables`); T, T_c, cu may be Duals."
function photorates_g(p::Photo, T, Tc, cu, ehz::Vector{Float64})
    tb = lj_gtables(p, ehz); G = length(ehz) + 1
    x = clamp((log10(Tc) - LJG_C[1]) / LJG_C[2], 0.0, LJG_C[3] - 1.0); j = min(floor(Int, _ljv(x)) + 1, LJG_C[3] - 1); w = x - (j - 1)
    Z = p.Z(T); R = zero(T) * zero(cu[1]) * zero(Tc)
    for (l, lv) in enumerate(p.levels)
        s = zero(R)
        for g in 1:G
            a = tb[g, l, j]; b = tb[g, l, j+1]; isfinite(a) || continue
            q = a + w * (b - a); q > -600 && (s += cu[g] * exp(q))
        end
        R += lv[2] / Z * exp(-lv[1] / (kB * T)) * s
    end
    R
end
"Stimulated part of the Milne coefficient in J = W B_ν(T_c) per unit W: α = α₀(T)(1 + W r(T, T_c)) of `photorates_lj`; 0 outside its tables."
function lj_stim(p::Photo, T, Tc)
    lT = log10(T); lc = log10(Tc); lTv = _ljv(lT)
    ok = LJ_T[1] + 2LJ_T[2] <= lTv <= LJ_T[1] + LJ_T[2] * (LJ_T[3] - 3) && LJ_2C[1] + 2LJ_2C[2] <= lc <= LJ_2C[1] + LJ_2C[2] * (LJ_2C[3] - 3) &&
         LJ_2T[1] + 2LJ_2T[2] <= lTv <= LJ_2T[1] + LJ_2T[2] * (LJ_2T[3] - 3)
    ok || return zero(T)
    tb = lj_tables(p)
    it, wt = _ljidx(lT, LJ_T); a0 = exp(_cr4(tb.lA0[it-1], tb.lA0[it], tb.lA0[it+1], tb.lA0[it+2], wt))
    jt, vt = _ljidx(lT, LJ_2T); jc, vc = _ljidx(lc, LJ_2C)
    row(j) = _cr4(tb.lr[jt-1, j], tb.lr[jt, j], tb.lr[jt+1, j], tb.lr[jt+2, j], vt)
    a0 * exp(_cr4(row(jc - 1), row(jc), row(jc + 1), row(jc + 2), vc))
end
"""
(R, α) of process p in a zone whose radiation field is the static field of `env.tabkey` (the CMB, tabulated) plus the multigroup star light
`env.mgq = (; Tc, cu, W, ehz)`: R = R_static + Σ_g c u_g Q_g(T_c) (`photorates_g`), α = α_static + W_g α₀ r(T, T_c) with the dilution W_g of the
group that holds the ground level's threshold (`lj_stim`).
"""
function photorates_mgq(p::Photo, T, env)
    q = env.mgq; R0, A0 = _photorates_t(p, T, env.tabkey)
    g = searchsortedlast(q.ehz, p.levels[1][3] / hP * 1.01) + 1      # (the group just above the threshold; 1 %: the H⁻ edge 0.7551 eV against the pair's 0.754 eV)
    (R0 + photorates_g(p, T, q.Tc, q.cu, q.ehz), p.milne ? A0 + q.W[g] * lj_stim(p, T, q.Tc) : A0 + zero(q.W[g]))
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

struct CKWS{T}
    x::Vector{T}; xt::Vector{T}; v::Vector{T}; veos::Vector{T}; ks::Vector{T}; dd::Vector{T}; fsl::Vector{T}
    ud1::Vector{T}; ud2::Vector{T}; w::Vector{T}
end
function CKWS{T}(nf::Int, nsp::Int, nrx::Int, nsl::Int) where {T}
    CKWS{T}(zeros(T, nf), zeros(T, nf), zeros(T, nsp), zeros(T, nsp), zeros(T, nrx), zeros(T, nsp), zeros(T, nsl),
            zeros(T, 4), zeros(T, 4), zeros(T, 5))
end

"""
Stack-allocated immutable composition dictionary (`full_fast(x, net)`), bit-identical to `full(x, net)` on every species.
"""
struct FullComp{T} <: AbstractDict{Symbol,T}
    H::T; Hp::T; Hm::T; H2::T; H2p::T; e::T; He::T; Hep::T; Hepp::T; Li::T; Lip::T
    hasHe::Bool; hasLi::Bool
end
@inline function Base.getindex(d::FullComp{T}, s::Symbol) where {T}
    s === :H ? d.H : s === :Hp ? d.Hp : s === :Hm ? d.Hm : s === :H2 ? d.H2 : s === :H2p ? d.H2p :
    s === :e ? d.e : s === :He ? d.He : s === :Hep ? d.Hep : s === :Hepp ? d.Hepp :
    s === :Li ? d.Li : s === :Lip ? d.Lip : throw(KeyError(s))
end
@inline function Base.get(d::FullComp{T}, s::Symbol, def) where {T}
    s === :H ? d.H : s === :Hp ? d.Hp : s === :Hm ? d.Hm : s === :H2 ? d.H2 : s === :H2p ? d.H2p :
    s === :e ? d.e : s === :He ? d.He : s === :Hep ? (d.hasHe ? d.Hep : def) : s === :Hepp ? (d.hasHe ? d.Hepp : def) :
    s === :Li ? (d.hasLi ? d.Li : def) : s === :Lip ? (d.hasLi ? d.Lip : def) : def
end
@inline Base.haskey(d::FullComp, s::Symbol) = s in (:H, :Hp, :Hm, :H2, :H2p, :e, :He) || (d.hasHe && s in (:Hep, :Hepp)) || (d.hasLi && s in (:Li, :Lip))

"""
Index form of a network (built once per Network, cached): species vector = net.species ∪ {He}; free-species positions;
reactant/product index lists; element and charge coefficients for the eliminated H, Li, e.
"""
mutable struct PIdx
    const p::Photo
    const jr::Int
    const jp::Vector{Int}
end
struct Compiled
    sp::Vector{Symbol}; idx::Dict{Symbol,Int}; free::Vector{Int}
    R::Vector{Vector{Int}}; P::Vector{Vector{Int}}
    iH::Int; iLi::Int; ie::Int; iHe::Int
    elH::Vector{Float64}; elLi::Vector{Float64}; q::Vector{Float64}       # per free species
    elHe::Vector{Float64}; hasHe::Bool                                     # He nuclei per free species; He eliminated (HECHEM)
    cv::Vector{Function}; E0::Vector{Float64}; isH2::Vector{Bool}
    cv0::Vector{Float64}; ktup::Any; nR1::Vector{Int}; nslow::Int
    pidx::Base.RefValue{Vector{PIdx}}; pidx_lk::ReentrantLock
    ifree::NTuple{7,Int}; hasLi::Bool; ekeys::Vector{Symbol}
    sp_kind::Vector{UInt8}; sp_q::Vector{Float64}                          # per compiled species: 0=e, 1=H, 2=He, 3=Li; charge
    fe_idx::NTuple{11,Int}                                                 # indices of (H2, H, Hp, He, Hep, Hepp, Hm, H2p, Li, Lip, e)
    ws_f64::Vector{CKWS{Float64}}; ws_pool::Vector{Dict{DataType,Any}}
end
const COMPILED = CowCache{Compiled}()
compiled(net) = cowget!(COMPILED, net) do
    sp = unique(vcat(net.species, [:He])); idx = Dict(s => i for (i, s) in enumerate(sp))
    fr = [idx[s] for s in net.free]
    fp(s) = (i = findfirst(==(s), net.free); i === nothing ? 0 : i)
    nslow = count(i -> !(net.free[i] in (:Hm, :H2p, :Hep, :Hepp, :Lip)), eachindex(net.free))
    hasHe = any(s -> SPECIES[s].el[2] > 0, net.free); hasLi = :Li in net.species
    dtmp = Dict{Symbol,Float64}()
    for s in net.free; dtmp[s] = 0.0; end
    dtmp[:He] = 0.0; dtmp[:H] = 0.0; hasLi && (dtmp[:Li] = 0.0)
    ekeys = [s for s in keys(dtmp) if s != :e]
    sp_kind = UInt8[s === :e ? 0x00 : (SPECIES[s].el[2] > 0 ? 0x02 : (SPECIES[s].el[3] > 0 ? 0x03 : 0x01)) for s in sp]
    sp_q = Float64[Float64(SPECIES[s].q) for s in sp]
    fe_idx = (get(idx, :H2, 0), get(idx, :H, 0), get(idx, :Hp, 0), get(idx, :He, 0),
              get(idx, :Hep, 0), get(idx, :Hepp, 0), get(idx, :Hm, 0), get(idx, :H2p, 0),
              get(idx, :Li, 0), get(idx, :Lip, 0), get(idx, :e, 0))
    nt = Threads.maxthreadid()
    Compiled(sp, idx, fr, [[idx[s] for s in r.R] for r in net.reactions], [[idx[s] for s in r.P] for r in net.reactions],
             get(idx, :H, 0), get(idx, :Li, 0), get(idx, :e, 0), idx[:He],
             [s == :H ? 0.0 : Float64(SPECIES[s].el[1]) for s in net.free], [Float64(SPECIES[s].el[3]) for s in net.free],
             [Float64(SPECIES[s].q) for s in net.free], [Float64(SPECIES[s].el[2]) for s in net.free], hasHe,
             [SPECIES[s].cv for s in sp], [SPECIES[s].E0 for s in sp], [s == :H2 for s in sp],
             [s == :H2 ? NaN : Float64(SPECIES[s].cv(300.0)) for s in sp], Tuple(r.k for r in net.reactions),
             [length(r.R) - 1 for r in net.reactions], nslow, Ref(PIdx[]), ReentrantLock(),
             (fp(:Hp), fp(:Hm), fp(:H2), fp(:H2p), fp(:Lip), fp(:Hep), fp(:Hepp)), hasLi, ekeys,
             sp_kind, sp_q, fe_idx,
             [CKWS{Float64}(length(fr), length(sp), length(net.reactions), nslow) for _ in 1:nt],
             [Dict{DataType,Any}() for _ in 1:nt])
end
@inline function _pidx(c::Compiled, p::Photo)
    list = c.pidx[]
    @inbounds for i in eachindex(list)
        e = list[i]
        e.p === p && return (e.jr, e.jp)
    end
    lock(c.pidx_lk) do
        for e in c.pidx[]
            e.p === p && return (e.jr, e.jp)
        end
        e = PIdx(p, c.idx[p.R], [c.idx[s] for s in p.P])
        c.pidx[] = vcat(c.pidx[], [e])
        (e.jr, e.jp)
    end
end
@inline get_ckws(c::Compiled, ::Type{Float64}) = @inbounds c.ws_f64[Threads.threadid()]
@inline function get_ckws(c::Compiled, ::Type{T}) where {T}
    d = @inbounds c.ws_pool[Threads.threadid()]
    w = get(d, T, nothing)
    w === nothing ? (w2 = CKWS{T}(length(c.free), length(c.sp), length(c.R), c.nslow); d[T] = w2; w2) : (w::CKWS{T})
end
function full_fast(x, net; xn = nothing)
    c = compiled(net); T = xn === nothing ? eltype(x) : promote_type(eltype(x), typeof(xn)); z = zero(T)
    iHp, iHm, iH2, iH2p, iLip, iHep, iHepp = c.ifree
    @inbounds begin
        Hp = iHp > 0 ? x[iHp] : z; Hm = iHm > 0 ? x[iHm] : z
        H2 = iH2 > 0 ? x[iH2] : z; H2p = iH2p > 0 ? x[iH2p] : z
        Lip = iLip > 0 ? x[iLip] : z; Hep = iHep > 0 ? x[iHep] : z; Hepp = iHepp > 0 ? x[iHepp] : z
    end
    He = c.hasHe ? yHe - (z + Hep + Hepp) : yHe + z
    sh = z
    if xn === nothing
        @inbounds for i in eachindex(c.free)
            el1 = Int(c.elH[i]); el1 > 0 && (sh += el1 * x[i])
        end
        H = 1 - sh
    else                                   # xn: the neutral hydrogen of the slow species, 1 − x_H⁺ − 2x_H₂ given directly (no cancellation); the fast ones are subtracted
        @inbounds for i in eachindex(c.free)
            (i == iHp || i == iH2) && continue
            el1 = Int(c.elH[i]); el1 > 0 && (sh += el1 * x[i])
        end
        H = xn - sh
    end
    Li = c.hasLi ? yLi - (z + Lip) : z
    d0 = FullComp{T}(H, Hp, Hm, H2, H2p, z, He, Hep, Hepp, Li, Lip, c.hasHe, c.hasLi)
    e = z
    @inbounds for s in c.ekeys
        e += SPECIES[s].q * d0[s]
    end
    FullComp{T}(H, Hp, Hm, H2, H2p, e, He, Hep, Hepp, Li, Lip, c.hasHe, c.hasLi)
end
"Full composition as a vector over compiled(net).sp (same elimination as `full`); entries of x behind the free species are ignored."
function fullvec!(v::Vector{T}, x, c::Compiled; xn = nothing) where {T}
    fill!(v, zero(T))
    @inbounds for (i, j) in enumerate(c.free); v[j] = x[i]; end
    v[c.iHe] = c.hasHe ? yHe - sum(c.elHe[i] * x[i] for i in eachindex(c.free)) : yHe
    if c.iH > 0
        if xn === nothing
            v[c.iH] = 1 - sum(c.elH[i] * x[i] for i in eachindex(c.free))
        else                               # xn given: H = xn − (the fast species' H nuclei), the slow ones (H⁺, H₂) are inside xn
            iHp, iH2 = c.ifree[1], c.ifree[3]
            v[c.iH] = xn - sum((c.elH[i] * x[i] for i in eachindex(c.free) if i != iHp && i != iH2); init = zero(T))
        end
    end
    c.iLi > 0 && (v[c.iLi] = yLi - sum(c.elLi[i] * x[i] for i in eachindex(c.free)))
    c.ie > 0 && (v[c.ie] = sum(c.q[i] * x[i] for i in eachindex(c.free)))
    v
end
fullvec(x, c::Compiled) = fullvec!(zeros(eltype(x), length(c.sp)), x, c)

"Chemical energy per gram Σ x_j E0_j/μ_H of composition x (the T-independent part of `eos`'s u; 0 for fully ionized gas)."
function echem(x, net)
    c = compiled(net); v = fullvec!(get_ckws(c, eltype(x)).veos, x, c)
    sum(v[j] * c.E0[j] for j in eachindex(v)) / μH()
end
"""
`UNIFIED_EOS=1` (default off; off = bit-identical): unified smooth C² microphysics EOS removing the fixed-primordial,
ideal-electron, no-pair, no-Coulomb assumptions from ChemNet. Uses per-zone elemental composition `(X, Y, Z)` (from
`UE_COMP[]` / NucNet, defaulting to primordial `X = 0.75, Y = 0.25, Z = 0` with real atomic masses matching SCVH95 and
StellarEOS), exact Fermi–Dirac + relativistic + e± pair electrons (`UE_ELECTRON[]` = `StellarEOS.electron_pe`), and
YS89 Coulomb corrections (`UE_COULOMB[]` = `StellarEOS.coulomb_pe`). Energy zero matches ChemNet & PopIIIEqEOS
(`E0(H⁺) = 0, E0(H) = -χH, E0(He) = 0, E0(He²⁺) = χHe1 + χHe2`, bare metals at `+χ_m`), eliminating any Δ0 offset.
"""
const UNIFIED_EOS = get(ENV, "UNIFIED_EOS", "0") == "1"
const UNIFIED_EOS_ON = Ref(UNIFIED_EOS)
const FREE_EOS = get(ENV, "FREE_EOS", "0") == "1"
const FREE_EOS_ON = Ref(FREE_EOS)
const UE_ELECTRON = Ref{Any}(nothing)      # (ne, T; pairs) -> (Pe, Ee)
const UE_COULOMB = Ref{Any}(nothing)       # (ni, zbar, T) -> (Pc, Ec)
const UE_COMP = Ref{Any}(nothing)          # () -> comp of current zone
const UE_COMPK = Ref{Any}(nothing)         # (k) -> comp of zone k
const UE_PAIRS = Ref(true)
const UE_COULOMB_ON = Ref(true)
const FE_EOS_NONLTE = Ref{Any}(nothing)    # FreeEnergyEOS.eos_nonlte_dens
const FE_EOS_FULLION = Ref{Any}(nothing)   # FreeEnergyEOS.eos_full_ion
const FE_DELTA_MU = Ref{Any}(nothing)      # FreeEnergyEOS.delta_mu_ex_7
const FE_DMU_NH_MIN = Ref(5e11)            # n_H threshold (~1e-12 g cm⁻³) above which Δμ_i^ex is evaluated in ChemNet
@inline function _dmu_sp(s::Symbol, dmu)
    s === :H2 ? dmu[1] : s === :H ? dmu[2] : s === :Hp ? dmu[3] :
    s === :He ? dmu[4] : s === :Hep ? dmu[5] : s === :Hepp ? dmu[6] :
    s === :e ? dmu[7] : zero(dmu[1])
end
@inline function _dmu_rxn(R, P, dmu)
    s = zero(dmu[1])
    for p in P; s += _dmu_sp(p, dmu); end
    for r in R; s -= _dmu_sp(r, dmu); end
    s
end
@inline function _fe_dmu_from_v(nH, T, v, c::Compiled)
    (FREE_EOS_ON[] && FE_DELTA_MU[] !== nothing && ForwardDiff.value(nH) >= FE_DMU_NH_MIN[]) || return nothing
    iH2, iH, iHp, iHe, iHep, iHepp, _, _, _, _, ie = c.fe_idx
    z = zero(nH * eltype(v)(0))
    pos(idx) = idx > 0 ? max(nH * v[idx], z) : z
    # Scale He by real mass ratio so nHe matches FreeEnergyEOS primordial (X=0.75, Y=0.25) exactly
    he_fac = (0.25 / _MHE) / (0.75 / mH) / yHe
    nH2 = pos(iH2); nH1 = pos(iH); nHp = pos(iHp)
    nHe = he_fac * pos(iHe); nHep = he_fac * pos(iHep); nHepp = he_fac * pos(iHepp)
    ne = ie > 0 ? max(nH * v[ie] + (he_fac - 1.0) * (pos(iHep) + 2pos(iHepp)), z) : (nHp + nHep + 2nHepp)
    FE_DELTA_MU[](T, nH2, nH1, nHp, nHe, nHep, nHepp, ne; pairs = UE_PAIRS[], coulomb = UE_COULOMB_ON[])
end
const _AMU = 1.66053906660e-24
const _MHE = 4.002602 * _AMU
const PRIM_COMP = (X = 0.75, Y = 0.25, Z = 0.0,
                   x = (0.75, 0.25, 0.0, 0.0, 0.0),
                   m = (mH, _MHE, 12.011 * _AMU, 14.007 * _AMU, 15.999 * _AMU),
                   zc = (1.0, 2.0, 6.0, 7.0, 8.0),
                   χ = (χH, χHe1 + χHe2, 489.99eV, 667.05eV, 871.41eV),
                   νion = 0.75 / mH + 0.25 / _MHE,
                   νe = 0.75 / mH + 2 * 0.25 / _MHE)
@inline zone_comp() = UE_COMP[] === nothing ? PRIM_COMP : UE_COMP[]()
@inline zone_comp_k(k) = UE_COMPK[] === nothing ? PRIM_COMP : UE_COMPK[](k)
@inline _ue_smoother5(x) = (x = clamp(x, zero(x), one(x)); x^3 * (10 - 15x + 6x^2))
"C² metal ionisation weight: 0 below 10⁵ K (neutral CNO), 1 above 10^{5.6} K (bare CNO nuclei)."
@inline _ue_metal_w(T) = _ue_smoother5((log10(T) - 5.0) / 0.6)

"Unified microphysics EOS (before optional table blend) for species vector `v` and elemental composition `comp`."
function eos_unified_raw(ρ, T, v, c::Compiled, comp; pairs::Bool = UE_PAIRS[], coulomb::Bool = UE_COULOMB_ON[])
    if FREE_EOS_ON[] && FE_EOS_NONLTE[] !== nothing
        nH_tot = ρ * comp.x[1] / comp.m[1]
        nHe_scale = (ρ * comp.x[2] / comp.m[2]) / yHe
        iH2, iH, iHp, iHe, iHep, iHepp, iHm, iH2p, iLi, iLip, _ = c.fe_idx
        z = zero(ρ * T * eltype(v)(0))
        gH(idx)  = idx > 0 ? max(nH_tot * v[idx], z) : z
        gHe(idx) = idx > 0 ? max(nHe_scale * v[idx], z) : z
        nZ = (comp.x[3] > 0 ? ρ * comp.x[3] / comp.m[3] : zero(ρ),
              comp.x[4] > 0 ? ρ * comp.x[4] / comp.m[4] : zero(ρ),
              comp.x[5] > 0 ? ρ * comp.x[5] / comp.m[5] : zero(ρ))
        return FE_EOS_NONLTE[](ρ, T, gH(iH2), gH(iH), gH(iHp), gHe(iHe), gHe(iHep), gHe(iHepp);
                               nHm = gH(iHm), nH2p = gH(iH2p), nLi = gH(iLi), nLip = gH(iLip),
                               nZ, pairs, coulomb, hs = true)
    end
    νH = comp.x[1] / comp.m[1]
    νHe_tot = comp.x[2] / comp.m[2]
    he_scale = νHe_tot / yHe
    ET = promote_type(typeof(ρ), typeof(T), eltype(v), Float64)
    u_heavy = zero(ET); ν_heavy = zero(ET); ν_e = zero(ET); ν_ion = zero(ET); ν_qpos = zero(ET)
    @inbounds for j in eachindex(v)
        knd = c.sp_kind[j]
        knd == 0x00 && continue
        νj = (knd == 0x02 ? he_scale : νH) * v[j]
        u_heavy += νj * ((c.isH2[j] ? u_H2(T) : c.cv0[j]) * kB * T + c.E0[j])
        ν_heavy += νj
        q = c.sp_q[j]
        if q != 0.0
            ν_e += q * νj
            if q > 0.0
                ν_ion += νj
                ν_qpos += q * νj
            end
        end
    end
    if comp.x[3] > 0 || comp.x[4] > 0 || comp.x[5] > 0
        wZ = _ue_metal_w(T)
        @inbounds for i in 3:5
            comp.x[i] > 0 || continue
            νm = comp.x[i] / comp.m[i]
            ν_heavy += νm
            u_heavy += νm * (1.5 * kB * T + wZ * comp.χ[i])
            ν_e += wZ * comp.zc[i] * νm
            ν_ion += wZ * νm
            ν_qpos += wZ * comp.zc[i] * νm
        end
    end
    if HM88 && c.iH > 0 && νH > 0 && ForwardDiff.value(v[c.iH]) > 0
        Nneu = ρ * (νH * v[c.iH] + he_scale * v[c.iHe])
        @inbounds for j in eachindex(v); c.isH2[j] && (Nneu += ρ * νH * v[j]); end
        Nion = ρ * (ForwardDiff.value(ν_ion) > 0 ? ν_ion : zero(ν_ion))
        eex, pex = hm_exc(T, Nneu, Nion)
        u_heavy += νH * v[c.iH] * eex
        ν_heavy += νH * v[c.iH] * pex
    end
    p = ρ * ν_heavy * kB * T
    u = u_heavy
    ne = ForwardDiff.value(ν_e) > 0 ? ρ * ν_e : zero(ρ * ν_e)
    Pe, Ee = UE_ELECTRON[](ne, T; pairs)
    p += Pe
    u += Ee / ρ
    if coulomb
        ni = ρ * ν_ion
        if ForwardDiff.value(ni) > 1e-30
            zbar = ν_qpos / ν_ion
            Pc, Ec = UE_COULOMB[](ni, zbar, T)
            p += Pc
            u += Ec / ρ
        end
    end
    return (u, p)
end

"Exact fully-ionised limit of `eos_unified_raw` (H⁺, He²⁺, bare CNO, Fermi–Dirac e±, Coulomb; ChemNet energy zero)."
function eos_full_ion(ρ, T, comp; pairs::Bool = UE_PAIRS[], coulomb::Bool = UE_COULOMB_ON[])
    if FREE_EOS_ON[] && FE_EOS_FULLION[] !== nothing
        return FE_EOS_FULLION[](ρ, T, comp; pairs, coulomb)
    end
    νH = comp.x[1] / comp.m[1]
    νHe = comp.x[2] / comp.m[2]
    ET = promote_type(typeof(ρ), typeof(T), Float64)
    ν_heavy = ET(νH + νHe)
    ν_ion = ET(νH + νHe)
    ν_e = ET(νH + 2νHe)
    u_heavy = (1.5 * (νH + νHe)) * kB * T + νHe * (χHe1 + χHe2)
    if comp.x[3] > 0 || comp.x[4] > 0 || comp.x[5] > 0
        wZ = _ue_metal_w(T)
        @inbounds for i in 3:5
            comp.x[i] > 0 || continue
            νm = comp.x[i] / comp.m[i]
            ν_heavy += νm
            u_heavy += νm * (1.5 * kB * T + wZ * comp.χ[i])
            ν_e += wZ * comp.zc[i] * νm
            ν_ion += wZ * νm
        end
    end
    p = ρ * ν_heavy * kB * T
    u = u_heavy
    Pe, Ee = UE_ELECTRON[](ρ * ν_e, T; pairs)
    p += Pe
    u += Ee / ρ
    if coulomb
        ni = ρ * ν_ion
        if ForwardDiff.value(ni) > 1e-30
            zbar = ν_e / ν_ion
            Pc, Ec = UE_COULOMB[](ni, zbar, T)
            p += Pc
            u += Ec / ρ
        end
    end
    return (u, p)
end

"Internal energy per gram and pressure."
function eos(ρ, T, x, net; xn = nothing, comp = nothing)
    c = compiled(net); v = fullvec!(get_ckws(c, xn === nothing ? eltype(x) : promote_type(eltype(x), typeof(xn))).veos, x, c; xn)
    if UNIFIED_EOS_ON[]
        cp = comp === nothing ? zone_comp() : comp
        u, p = eos_unified_raw(ρ, T, v, c, cp)
        EOSMIX[] === nothing && return (u, p)
        return EOSMIX[](ρ, T, u, p)
    end
    nH = ρ / μH()
    u = zero(T); N = zero(T)
    @inbounds for j in eachindex(v)
        u += v[j] * ((c.isH2[j] ? u_H2(T) : c.cv0[j]) * kB * T + c.E0[j])
        N += v[j]
    end
    if HM88 && c.iH > 0                    # HM88: excitation energy and occupation-probability pressure of H I (see HM88)
        Nneu, Nion = hm_perturbers(v, c, nH)
        eex, pex = hm_exc(T, Nneu, Nion)
        EOSMIX[] === nothing && return ((u + v[c.iH] * eex) / μH(), (N + v[c.iH] * pex) * nH * kB * T)
        return EOSMIX[](ρ, T, (u + v[c.iH] * eex) / μH(), (N + v[c.iH] * pex) * nH * kB * T)
    end
    EOSMIX[] === nothing && return (u / μH(), N * nH * kB * T)
    EOSMIX[](ρ, T, u / μH(), N * nH * kB * T)
end
"""
EOSMIX[] = f (default nothing = off): `eos` returns f(ρ, T, u, p) of its own (u, p) — set by the material (PopIIIMaterial G2BLEND: the smooth
network/table blend in the non-ideal regime), so that the chemistry's T(u) and the hydro's EOS are the same function.
"""
const EOSMIX = Ref{Any}(nothing)

# ------------------------------------------------------------------ HM88: occupation probabilities (G2, docs/POPIII_G2_MATERIAL.md)
"""
`HM88=1` (default off; off = bit-identical): atomic H with the occupation-probability formalism of Hummer & Mihalas (1988, ApJ 331, 794;
w_n of their eq. 4.71 with Z_p = 1, as Z0Physics.occupation): partition function Z_w = Σ_{n ≤ $(80)} 2n² (w_n/w_1) e^{−E_n/kT} instead of
ZH_el's hard cut at n = 15 without w (which overpopulates neutral H by Z_H/2 ≈ 10–500 above 3e4 K and has no excitation energy).
- EOS (`eos`, exact in ρ, T, x): u += x_H Σ_n P_n E_n (excitation energy, P_n = level fractions of Z_w), p += n_H x_H kT Σ_n P_n (−ln(w_n/w_1))
  (−∂F_int/∂V at fixed particle numbers, F_int = −kT N_H ln Z_w; the cross terms through the perturbers' chemical potentials are omitted).
- Chemistry: the equilibrium constant of H ⇌ H⁺ + e becomes K Z_el/Z_w; every H recombination channel (reactions H⁺ + e (+ …) → H (+ …),
  incl. the detailed-balance reverses and the three-level C α_B, and the Milne inverse of H photoionization) is multiplied by
  φ = Z_w/Z_el ≤ 1, so detailed balance and J = B give the HM88 Saha equilibrium; ionization rates are unchanged (ionization from dissolved
  levels is lumped into the recombination deficit). φ takes the perturber densities from `env.hm = (N_neu, N_ion)` [cm⁻³], set per zone and
  step by the material from the composition at the start of the step (lagged within the step; exact in T); without `env.hm`, φ = 1.
Networks: wrap with `hm88_network(net)`.
"""
const HM88 = get(ENV, "HM88", "0") == "1"
const NHM = 80
const A0CM = 5.29177211e-9; const ECGS = 4.80320471e-10; const RYDINF = 13.605693eV
"HM88 occupation probability of level n of a hydrogenic ion of charge Z in neutral/ion perturbers (cm⁻³), eq. 4.71 (Z_p = 1)."
@inline function w_hm(n, Z, Nneu, Nion)
    exp(_lw_hm(n, Z, Nneu, Nion))
end
"ln w_n of `w_hm`, its exponent (the same expression: `w_hm` is bit-identical to before)."
@inline function _lw_hm(n, Z, Nneu, Nion)
    χn = Z^2 * RYDINF / n^2; K = n <= 3 ? 1.0 : (16 / 3) * n / (n + 1)^2
    q = sqrt(Float64(Z)) * ECGS^2 / (sqrt(K) * χn); rn = n^2 * A0CM / Z
    -(4π / 3) * (16 * Nion * q^3 + Nneu * (rn + A0CM)^3)
end
"Z_w/2 of H I (ground weight 1, w_n/w_1)."
function ZH_hm(T, Nneu, Nion)
    w1 = w_hm(1, 1, Nneu, Nion); z = zero(T) + zero(Nneu) + zero(Nion)
    for n in 1:NHM; z += n^2 * (w_hm(n, 1, Nneu, Nion) / w1) * exp(-χH * (1 - 1 / n^2) / (kB * T)); end
    z
end
"(⟨E_exc⟩, Σ P_n (−ln(w_n/w_1))) per H atom [erg, dimensionless] of the HM88 level distribution."
function hm_exc(T, Nneu, Nion)
    w1 = w_hm(1, 1, Nneu, Nion); lw1 = _lw_hm(1, 1, Nneu, Nion); z = zero(T) + zero(Nneu) + zero(Nion); e = z; q = z
    for n in 1:NHM
        wr = w_hm(n, 1, Nneu, Nion) / w1; En = χH * (1 - 1 / n^2); t = n^2 * wr * exp(-En / (kB * T))
        # ln(w_n/w_1): log(wr) as ever where w_n/w_1 ≥ 1e-150; below, where w underflows towards 0 (the weight t ≤ 1e-143 there: far below one rounding of q) the exponents' difference, which
        # stays finite to every derivative order (the second derivative of log(w) for w ≲ 1e-154 is Inf/NaN: w'² underflows). M6: exact nested ForwardDiff through the EOS.
        lwr = wr >= 1e-150 ? log(wr) : _lw_hm(n, 1, Nneu, Nion) - lw1
        z += t; e += t * En; q -= t * lwr
    end
    if ForwardDiff.value(z) == 0.0
        return zero(e), zero(q)
    end
    e / z, q / z
end
"Perturber densities of composition v (fullvec) at n_H: neutrals H, He, H₂ (and Li); ions as the electrons."
function hm_perturbers(v, c, nH)
    Nneu = v[c.iH] + v[c.iHe]
    for j in eachindex(v); c.isH2[j] && (Nneu += v[j]); end
    Nion = c.ie > 0 ? max(v[c.ie], zero(v[c.ie])) : zero(Nneu)
    nH * Nneu, nH * Nion
end
"φ = Z_w/Z_el of the zone's env (1 without env.hm, with HM88 off, or when FREE_EOS supplies full Potekhin unbounding)."
hm_phi(T, env) = (HM88 && !FREE_EOS_ON[] && hasproperty(env, :hm)) ? 2ZH_hm(T, env.hm[1], env.hm[2]) / ZH_el(T) : one(T)
"The network with every H recombination channel × φ (HM88)."
function hm88_network(net)
    isrec(r) = (:Hp in r.R) && (:e in r.R) && (:H in r.P) && !(:H in r.R)
    rx = [isrec(r) ? Reaction(r.name * " [HM88 φ]", r.R, r.P, let k = r.k; (T, env) -> k(T, env) * hm_phi(T, env); end) : r for r in net.reactions]
    Network(net.species, net.free, rx)
end

"dx/dt for the free species at H-nucleus density nH."
function dxdt(x, nH, T, net, env, photo)
    c = compiled(net); v = fullvec(x, c)
    dmu = _fe_dmu_from_v(nH, T, v, c)
    ev = dmu === nothing ? env : merge(env, (; dmu))
    ET = dmu === nothing ? promote_type(eltype(x), typeof(T)) : promote_type(eltype(x), typeof(T), typeof(dmu[1]))
    dd = zeros(ET, length(v))
    for (ir, r) in enumerate(net.reactions)
        Ri = c.R[ir]; rate = r.k(T, ev) * nH^(length(Ri) - 1)
        for j in Ri; rate *= v[j]; end
        for j in Ri; dd[j] -= rate; end
        for j in c.P[ir]; dd[j] += rate; end
    end
    for (p, R, α) in photo                          # precomputed per zone: photo-rate R, inverse coefficient α
        jr = c.idx[p.R]; jp = [c.idx[s] for s in p.P]
        rate = R * v[jr]; for j in jp; dd[j] += rate; end; dd[jr] -= rate
        α_eff = dmu === nothing ? α : α * exp(_dmu_rxn((p.R,), p.P, dmu))
        inv_ = α_eff * nH^(length(jp) - 1); for j in jp; inv_ *= v[j]; end    # per H nucleus: α n_H^(k−1) Π x
        for j in jp; dd[j] -= inv_; end; dd[jr] += inv_
    end
    dd[c.free]
end

"dxdt with the rate coefficients precomputed: ks[ir] = k_ir(T, env) n_H^(n_R − 1) (as in `ck_expand`); same result as `dxdt`."
function dxdt_k!(ws::CKWS{ET}, x, nH, net, c::Compiled, ph, ks; xn = nothing) where {ET}
    v = fullvec!(ws.v, x, c; xn); dd = fill!(ws.dd, zero(ET))
    @inbounds for ir in eachindex(net.reactions)
        Ri = c.R[ir]; rate = ks[ir]
        for j in Ri; rate *= v[j]; end
        for j in Ri; dd[j] -= rate; end
        for j in c.P[ir]; dd[j] += rate; end
    end
    @inbounds for (p, R, α) in ph
        jr, jp = _pidx(c, p)
        rate = R * v[jr]; for j in jp; dd[j] += rate; end; dd[jr] -= rate
        inv_ = α * nH^(length(jp) - 1); for j in jp; inv_ *= v[j]; end
        for j in jp; dd[j] -= inv_; end; dd[jr] += inv_
    end
    dd
end
function dxdt_k(x, nH, net, ph, ks)
    c = compiled(net); ET = promote_type(eltype(x), eltype(ks)); ws = get_ckws(c, ET)
    dd = dxdt_k!(ws, x, nH, net, c, ph, ks)
    dd[c.free]
end

"""
From this recursion depth, a failed Newton is retried with each free abundance clipped on its own (x → 0.1 x) instead of scaling the
whole correction. A trace species whose linearised update is strongly negative (H₂⁺ at n ≈ 1e16, T ≈ 3300 K, with the H₂ row far from
balance) otherwise sets λ ≈ 0.9 x/|δ| for all unknowns: x falls by 10 % per iteration, λ → 1e-60, and Newton freezes for every h.
Runs that never reach this depth are unchanged bit for bit.
"""
const CLIPDEPTH = 4
"""
`CHEMSTATS=1` (diagnostic, default off): count calls of `step`, sub-step attempts (`adv` calls: halvings included) and Jacobians
(= Newton iterations); totals printed at exit.
"""
const CHEMSTATS = get(ENV, "CHEMSTATS", "0") == "1"
const NSTEP = Threads.Atomic{Int}(0); const NADV = Threads.Atomic{Int}(0); const NJAC = Threads.Atomic{Int}(0)
CHEMSTATS && atexit(() -> println("CHEMSTATS: ChemNet.step calls ", NSTEP[], ", sub-step attempts ", NADV[], ", Jacobians (Newton iterations) ", NJAC[]))
"Print the Newton iterations of the first N failed steps (ENV CHEMTRACE; diagnostics)."
const TRACEDBG = Ref(parse(Int, get(ENV, "CHEMTRACE", "0")))


# ------------------------------------------------------------------ CHEMCK: reduced network (ChemistryKernels approach)
"""
`CHEMCK=1` (default off): the approach of ChemistryKernels.jl. Only H₂, H⁺ (the slow species) and T are Newton unknowns of the
backward-Euler step; the fast species are algebraic functions of them: H⁻ and H₂⁺ in chemical steady state (creation = destruction,
`CHEMCKSWEEPS` Gauss–Seidel sweeps, default 3, as CK evaluates them twice), He, He⁺, He²⁺ and Li, Li⁺ in ionisation equilibrium
(every process of the network between adjacent stages balanced: collisional ionisation ⇌ three-body recombination, photo-ionisation ⇌
radiative recombination, charge transfer), with n_e from charge conservation. The fast species cannot go negative and no longer limit
the Newton step (the stall of 13 Nov: λ falling ×10 per iteration on trace He⁺/He²⁺ and the eliminated Li). Assumes the fast species
relax faster than the step. `extra` (star light: H⁻ photodetachment, LW, H photoionisation; a ray) enters the slow rows directly and
the H⁻/H₂⁺ balance as a source and a loss linear in the species' own abundance (extra with the fast singles at 0 and at 1, evaluated
once per call at the slow state); not the He/Li stages.
"""
const CK = Ref(get(ENV, "CHEMCK", "0") == "1"); const CKSWEEPS = Ref(parse(Int, get(ENV, "CHEMCKSWEEPS", "3")))
"""
`CHEMCKHSS=f` (with CHEMCK=1; default 0 = off): H⁺ in steady state where it relaxes fast. In a three-level zone (env carries `tl3`) n_e in
the Peebles C is taken from the current iterate (x_H⁺) instead of the start-of-step value, and where |∂f_H⁺/∂x_H⁺| dt > f at the start of
the step the H⁺ row is the steady state f_H⁺ = 0 instead of backward Euler. The pre-shock zones (n 3–9e12, T 5300–6700 K) ionize by an
electron avalanche e-folding in seconds within steps of 1e5–1e6 s, which the material's CHEMSUB tree (lagged n_e) bisected to ~10 s
(48 % of the chemistry time in 4 zones, t_cka, 13 Nov). The material skips CHEMSUB and takes the same n_e in tl3_heat when this is on.
"""
const CKHSS = Ref(parse(Float64, get(ENV, "CHEMCKHSS", "0")))
const CKFAST = (:Hm, :H2p, :Hep, :Hepp, :Lip)
struct CKInfo
    fast::Vector{Int}; slow::Vector{Int}                 # positions in net.free
    single::Vector{Tuple{Int,Int}}                       # H⁻, H₂⁺: (position in net.free, index in compiled sp)
    chains::Vector{Vector{Tuple{Int,Int}}}               # element ionisation stages from the neutral up: (free position or 0, sp index)
    kHp::Int                                             # index of :Hp in slow (0 if not present)
end
const CKINFO = CowCache{CKInfo}()
ckinfo(net) = cowget!(CKINFO, net) do
    c = compiled(net); fp(s) = (i = findfirst(==(s), net.free); i === nothing ? 0 : i)
    fast = [fp(s) for s in CKFAST if fp(s) > 0]; slow = [i for i in eachindex(net.free) if !(i in fast)]
    single = [(fp(s), c.idx[s]) for s in (:Hm, :H2p) if fp(s) > 0]
    chains = Vector{Tuple{Int,Int}}[]
    for ch in ((:He, :Hep, :Hepp), (:Li, :Lip))
        all(s -> haskey(c.idx, s), ch) && any(s -> fp(s) > 0, ch[2:end]) || continue
        push!(chains, [(fp(s), c.idx[s]) for s in ch])
    end
    kHp = findfirst(==(:Hp), net.free[slow])
    CKInfo(fast, slow, single, chains, kHp === nothing ? 0 : kHp)
end
@generated function _fill_ks!(ks::Vector{ET}, ktup::KT, nR1::Vector{Int}, T, ev, nH) where {ET,KT<:Tuple}
    M = fieldcount(KT)
    ex = Expr(:block)
    for i in 1:M
        push!(ex.args, :(@inbounds ks[$i] = ktup[$i](T, ev) * nH^nR1[$i]))
    end
    push!(ex.args, :ks)
    ex
end
"k · Π v[Ri] with one occurrence of index j left out."
@inline function _rwo(r0, v, Ri, j)
    r = r0; sk = true
    for i in Ri; if sk && i == j; sk = false; else r *= v[i]; end; end
    r
end
"Creation rate of sp index j and its loss rate per unit x_j (per H nucleus), ks = k n_H^(n−1) per reaction, ph per-zone photo-rates."
function _ck_cl_p(v, j, nH, net, c::Compiled, ph, ks)
    C = zero(eltype(v)); L = zero(eltype(v))
    @inbounds for ir in eachindex(net.reactions)
        Ri = c.R[ir]; Pi = c.P[ir]; nr = count(==(j), Ri); np = count(==(j), Pi)
        nr == 0 && np == 0 && continue
        nr > 0 && (L += nr * _rwo(ks[ir], v, Ri, j))
        np > 0 && (C += np * _rwo(ks[ir], v, Ri, 0))
    end
    @inbounds for (p, R, α) in ph
        jr, jp = _pidx(c, p)
        inv0 = α * nH^(length(jp) - 1)
        if jr == j; L += R; C += _rwo(inv0, v, jp, 0); end
        if j in jp; C += R * v[jr]; L += _rwo(inv0, v, jp, j); end
    end
    C, L
end
"Rates per unit x_a of a → b and per unit x_b of b → a (a, b sp indices of adjacent ionisation stages)."
function _ck_ud_p(v, a, b, nH, net, c::Compiled, ph, ks)
    up = zero(eltype(v)); dn = zero(eltype(v))
    @inbounds for ir in eachindex(net.reactions)
        Ri = c.R[ir]; Pi = c.P[ir]
        a in Ri && b in Pi && (up += _rwo(ks[ir], v, Ri, a))
        b in Ri && a in Pi && (dn += _rwo(ks[ir], v, Ri, b))
    end
    @inbounds for (p, R, α) in ph
        jr, jp = _pidx(c, p)
        jr == a && b in jp && (up += R; dn += _rwo(α * nH^(length(jp) - 1), v, jp, b))
    end
    up, dn
end
function ck_expand!(ws::CKWS{ET}, xs, T, nH, net, c::Compiled, ci::CKInfo, ph, ks::Vector{ET}, extra = nothing; xn = nothing) where {ET}
    x = fill!(ws.x, zero(ET))
    @inbounds for (k, i) in enumerate(ci.slow); x[i] = xs[k]; end
    # extra (star light, ray) in the H⁻/H₂⁺ rows, taken as a + b x_j: two calls, all fast singles at 0 and at 1
    # (xn: the neutral hydrogen of the slow species, given directly; an `extra` that takes it is called extra(x, xn))
    xex(f_, x_) = xn === nothing ? f_(x_) : f_(x_, xn)
    ea = eb = nothing
    if extra !== nothing && !isempty(ci.single)
        xt = copyto!(ws.xt, x); ea = xex(extra, xt); for (f, _) in ci.single; xt[f] = one(ET); end; eb = xex(extra, xt) .- ea
    end
    for _ in 1:CKSWEEPS[]
        for (f, j) in ci.single
            v = fullvec!(ws.v, x, c; xn); C, L = _ck_cl_p(v, j, nH, net, c, ph, ks)
            ea === nothing || (ea[f] > 0 && (C += ea[f]); eb[f] < 0 && (L -= eb[f]))
            x[f] = L > 0 ? C / L : zero(C)
        end
        for ch in ci.chains
            v = fullvec!(ws.v, x, c; xn); m = length(ch)
            ud1 = ws.ud1; ud2 = ws.ud2; w = ws.w
            @inbounds for k in 1:m-1
                u_, d_ = _ck_ud_p(v, ch[k][2], ch[k+1][2], nH, net, c, ph, ks)
                ud1[k] = u_; ud2[k] = d_
            end
            @inbounds for k in 0:m-1
                p1 = one(ET); for i in 1:k; p1 *= ud1[i]; end
                p2 = one(ET); for i in k+1:m-1; p2 *= ud2[i]; end
                w[k+1] = p1 * p2
            end
            sw = zero(ET); @inbounds for k in 1:m; sw += w[k]; end
            tot = ch[1][2] == c.iHe ? yHe : yLi
            @inbounds for k in 1:m
                ch[k][1] > 0 || continue
                x[ch[k][1]] = sw > 0 ? tot * w[k] / sw : zero(sw)
            end
        end
    end
    x
end
"Free abundances with the slow ones xs (over ckinfo(net).slow) and the fast ones from equilibrium at T."
function ck_expand(xs, T, nH, net, ph, ks, extra = nothing; xn = nothing)
    c = compiled(net); ci = ckinfo(net)
    ET = promote_type(eltype(xs), typeof(T), eltype(ks))
    ws = get_ckws(c, ET)
    copy(ck_expand!(ws, xs, T, nH, net, c, ci, ph, ET === eltype(ks) ? ks : copyto!(ws.ks, ks), extra; xn))
end
"""
    ck_eval(xs, ρ, T, net; env, phc, extra) -> (x, f)

COUPLED (TitanResidFast): the full free-species vector x of the reduced network (slow abundances `xs` over `ckinfo(net).slow`, fast ones
from `ck_expand`) and dx/dt of the slow species, at (ρ, T), generic in the number types (the hydro Newton differentiates through it). In a
three-level zone (env carries `tl3`) n_e of the Peebles C is n_H x_H⁺ of `xs` (self-consistent, as CHEMCKHSS). `phc(Tv)` returns the
per-zone photo-rate tuples at the Float temperature Tv (lagged in T as in `step`).
"""
function ck_eval(xs, ρ, T, net; env, phc, extra, xn = nothing)
    ph = phc(T)                                   # (a Dual T: the photo-rates carry ∂/∂T, see PopIIIMaterial.photo_cache)
    _ck_eval_barrier(xs, ρ, T, net, env, ph, extra, xn)
end
function _ck_eval_barrier(xs, ρ, T, net, env, ph::Vector{Tuple{Photo,TP,TP}}, extra, xn) where {TP}
    nH = ρ / μH(); c = compiled(net); ci = ckinfo(net); sl = ci.slow
    kHp = ci.kHp
    tl3 = kHp > 0 && hasproperty(env, :tl3) && env.tl3 !== nothing
    ev = if tl3
        t3 = merge(env.tl3, (; ne = nH * xs[kHp]))
        (HM88 && !FREE_EOS_ON[] && hasproperty(env, :hm)) ? (; tl3 = t3, hm = env.hm) : (; tl3 = t3)
    else
        env
    end
    ET = promote_type(eltype(xs), typeof(ρ), typeof(T), TP)
    ws = get_ckws(c, ET)
    ks = _fill_ks!(ws.ks, c.ktup, c.nR1, T, ev, nH)
    x = ck_expand!(ws, xs, T, nH, net, c, ci, ph, ks, extra; xn)
    ph_use = ph
    if FREE_EOS_ON[] && FE_DELTA_MU[] !== nothing && ForwardDiff.value(nH) >= FE_DMU_NH_MIN[]
        v0 = fullvec!(ws.v, x, c; xn)
        dmu = _fe_dmu_from_v(nH, T, v0, c)
        if dmu !== nothing
            ev_fe = merge(ev, (; dmu))
            ks = _fill_ks!(ws.ks, c.ktup, c.nR1, T, ev_fe, nH)
            ph_use = Tuple{Photo,ET,ET}[(p, ET(R), ET(α * exp(_dmu_rxn((p.R,), p.P, dmu)))) for (p, R, α) in ph]
            x = ck_expand!(ws, xs, T, nH, net, c, ci, ph_use, ks, extra; xn)
        end
    end
    dd = dxdt_k!(ws, x, nH, net, c, ph_use, ks; xn)
    fsl = ws.fsl
    if extra === nothing
        @inbounds for (k, i) in enumerate(sl); fsl[k] = dd[c.free[i]]; end
    else
        ex = xn === nothing ? extra(x) : extra(x, xn)
        @inbounds for (k, i) in enumerate(sl); fsl[k] = dd[c.free[i]] + ex[i]; end
    end
    x, fsl
end
function _ck_eval_barrier(xs, ρ, T, net, env, ph, extra, xn)
    nH = ρ / μH(); c = compiled(net); ci = ckinfo(net); sl = ci.slow
    kHp = ci.kHp
    tl3 = kHp > 0 && hasproperty(env, :tl3) && env.tl3 !== nothing
    ev = if tl3
        t3 = merge(env.tl3, (; ne = nH * xs[kHp]))
        (HM88 && !FREE_EOS_ON[] && hasproperty(env, :hm)) ? (; tl3 = t3, hm = env.hm) : (; tl3 = t3)
    else
        env
    end
    ET = promote_type(eltype(xs), typeof(ρ), typeof(T), eltype(ph) <: Tuple ? fieldtype(eltype(ph), 2) : typeof(T))
    ws = get_ckws(c, ET)
    ks = _fill_ks!(ws.ks, c.ktup, c.nR1, T, ev, nH)
    x = ck_expand!(ws, xs, T, nH, net, c, ci, ph, ks, extra; xn)
    ph_use = ph
    if FREE_EOS_ON[] && FE_DELTA_MU[] !== nothing && ForwardDiff.value(nH) >= FE_DMU_NH_MIN[]
        v0 = fullvec!(ws.v, x, c; xn)
        dmu = _fe_dmu_from_v(nH, T, v0, c)
        if dmu !== nothing
            ev_fe = merge(ev, (; dmu))
            ks = _fill_ks!(ws.ks, c.ktup, c.nR1, T, ev_fe, nH)
            ph_use = Tuple{Photo,ET,ET}[(p, ET(R), ET(α * exp(_dmu_rxn((p.R,), p.P, dmu)))) for (p, R, α) in ph]
            x = ck_expand!(ws, xs, T, nH, net, c, ci, ph_use, ks, extra; xn)
        end
    end
    dd = dxdt_k!(ws, x, nH, net, c, ph_use, ks; xn)
    fsl = ws.fsl
    if extra === nothing
        @inbounds for (k, i) in enumerate(sl); fsl[k] = dd[c.free[i]]; end
    else
        ex = xn === nothing ? extra(x) : extra(x, xn)
        @inbounds for (k, i) in enumerate(sl); fsl[k] = dd[c.free[i]] + ex[i]; end
    end
    x, fsl
end
function step_ck(x0, ρ, u, T0, dt, net; env, photos, cooling, extra, accepted, tol, maxit)
    nH = ρ / μH(); ci = ckinfo(net); sl = ci.slow; ns = length(sl); cc = compiled(net)
    phT = Ref(NaN); phc = Ref{Any}(nothing); ucur = Ref(Float64(u))
    kHp = findfirst(==(:Hp), net.free[sl]); tl3sc = CKHSS[] > 0 && kHp !== nothing && hasproperty(env, :tl3) && env.tl3 !== nothing
    envz(z) = tl3sc ? merge(env, (; tl3 = merge(env.tl3, (; ne = nH * z[kHp])))) : env
    hss = Ref(false)
    function xT(z)
        T = exp(z[ns+1]); Tv = ForwardDiff.value(T); ev = envz(z)
        Tv == phT[] || (phc[] = [(p, photorates_t(p, Tv, env)...) for p in photos]; phT[] = Tv)
        ph = phc[]; ks = [r.k(T, ev) * nH^(length(r.R) - 1) for r in net.reactions]
        x = ck_expand(z[1:ns], T, nH, net, ph, ks, extra)
        (x, T, ph, ev)
    end
    function fslow(z)
        x, T, ph, ev = xT(z)
        f = dxdt(x, nH, T, net, ev, ph); extra === nothing || (f = f .+ extra(x))
        x, T, f
    end
    function resid(z, xs0, h)
        x, T, f = fslow(z)
        r = x[sl] .- xs0 .- h .* f[sl]
        hss[] && (r = [i == kHp ? -h * f[sl[i]] : r[i] for i in eachindex(r)])
        ue = eos(ρ, T, x, net)[1] - ucur[] + (cooling === nothing ? zero(T) : h * cooling(ρ, T, full(x, net)) / ρ)
        vcat(r, ue / (kB * T / μH()))
    end
    function attempt(xs0, T0, h)
        z = vcat(collect(Float64, xs0), log(T0)); rel = Inf; dT = Inf; λ = 0.0
        for it in 1:maxit
            R = resid(z, xs0, h)
            CHEMSTATS && Threads.atomic_add!(NJAC, 1)
            Jm = ForwardDiff.jacobian(zz -> resid(zz, xs0, h), z)
            (all(isfinite, Jm) && all(isfinite, R)) || return (z, false)
            δ = -(Jm \ R); all(isfinite, δ) || return (z, false)
            λ = 1.0
            # a trace slow species (H₂ in ionized gas, H⁺ in molecular gas) is clipped on its own; only a major one limits the whole step
            for i in 1:ns; z[i] < 1e-6 && z[i] + δ[i] < 0 && (δ[i] = -0.9 * z[i]); end
            for i in 1:ns; z[i] + λ * δ[i] < 0 && 0.9 * z[i] / -δ[i] < λ && (λ = 0.9 * z[i] / -δ[i]); end
            # the eliminated H, linear in the slow species (the fast ones are trace)
            H0 = 1 - sum(cc.elH[i] * xi for (i, xi) in enumerate(xT(z)[1])); dH = -sum(cc.elH[sl[k]] * δ[k] for k in 1:ns)
            H0 + λ * dH < 0 && 0.9 * H0 / -dH < λ && (λ = 0.9 * H0 / -dH)
            abs(δ[ns+1]) * λ > 0.5 && (λ = 0.5 / abs(δ[ns+1]))
            z .+= λ .* δ
            rel = maximum(abs.(δ[1:ns]) ./ max.(abs.(z[1:ns]), 1e-20)); dT = abs(δ[ns+1])
            λ == 1.0 && rel < tol && dT < tol && return (z, true)
        end
        λ == 1.0 && rel < 1e-6 && dT < 1e-8 && return (z, true)
        (z, false)
    end
    function adv(xs0, T0, h, depth)
        CHEMSTATS && Threads.atomic_add!(NADV, 1)
        z, ok = attempt(xs0, T0, h)
        if ok
            accepted === nothing || accepted(xT(z)[1], h)
            return (z[1:ns], exp(z[ns+1]))
        end
        depth >= 12 && error("ChemNet.step_ck: Newton failed (nH = $nH, T0 = $T0, dt = $dt, u = $u)")
        xs1, T1 = adv(xs0, T0, h / 2, depth + 1)
        cooling === nothing || (ucur[] = Float64(eos(ρ, T1, xT(vcat(xs1, log(T1)))[1], net)[1]))
        adv(xs1, T1, h / 2, depth + 1)
    end
    if tl3sc
        z0 = vcat(collect(Float64, x0[sl]), log(T0))
        dfdx = ForwardDiff.derivative(y -> fslow([i == kHp ? y + zero(y) * z0[i] : z0[i] + zero(y) for i in eachindex(z0)])[3][sl[kHp]], z0[kHp])
        hss[] = abs(dfdx) * dt > CKHSS[]
    end
    zs, T = adv(collect(Float64, x0[sl]), T0, dt, 0)
    (xT(vcat(zs, log(T)))[1], T)
end

"""
    step(x0, ρ, u, T0, dt, net; env, photos = Photo[], cooling = nothing, extra = nothing, Tfix = nothing) -> (x, T)

Chemistry and temperature together over dt at fixed ρ and u (u − dt Λ/ρ with `cooling(ρ, T, d)`, d the full composition).
`extra(x)`: added to dx/dt, a function of the new abundances (a photon-conserving ray, IonRay). `Tfix`: hold T at this value instead of solving for it (tests).
`accepted(x, h)`: called after each accepted sub-step, so that a ray can count the photons used through the halvings.
"""
function step(x0, ρ, u, T0, dt, net; env, photos = Photo[], cooling = nothing, extra = nothing, Tfix = nothing, accepted = nothing, tol = 1e-10, maxit = 60)
    nH = ρ / μH(); nf = length(net.free)
    CHEMSTATS && Threads.atomic_add!(NSTEP, 1)
    CK[] && Tfix === nothing && !isempty(ckinfo(net).fast) && return step_ck(x0, ρ, u, T0, dt, net; env, photos, cooling, extra, accepted, tol, maxit)
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
        lim = :none
        for it in 1:maxit
            R = resid(z, x0, h)
            trace && println("  trace it ", it, " T ", exp(z[nf+1]), " x ", z[1:nf], " R ", R)
            CHEMSTATS && Threads.atomic_add!(NJAC, 1)
            Jm = ForwardDiff.jacobian(zz -> resid(zz, x0, h), z)
            δ = -(Jm \ R); all(isfinite, δ) || return (z, false)
            λ = 1.0; lim = :none
            clip && for i in 1:nf; z[i] + δ[i] < 0 && (δ[i] = -0.9 * z[i]); end
            # keep every abundance ≥ 0, incl. the eliminated H, Li, e (linear in x)
            for i in 1:nf; z[i] + λ * δ[i] < 0 && 0.9 * z[i] / -δ[i] < λ && (λ = 0.9 * z[i] / -δ[i]; lim = net.free[i]); end
            d0 = full(z[1:nf], net); d1 = full(z[1:nf] .+ δ[1:nf], net)
            for s in (:H, :Li, :e, :He)
                haskey(d0, s) || continue
                dv = d1[s] - d0[s]
                d0[s] + λ * dv < 0 && 0.9 * d0[s] / -dv < λ && (λ = 0.9 * d0[s] / -dv; lim = s)
            end
            abs(δ[nf+1]) * λ > 0.5 && (λ = 0.5 / abs(δ[nf+1]); lim = :T)
            z .+= λ .* δ
            rel = maximum(abs.(δ[1:nf]) ./ max.(abs.(z[1:nf]), 1e-20)); dT = abs(δ[nf+1])                                # 1e-20: abundance floor
            λ == 1.0 && rel < tol && dT < tol && return (z, true)
        end
        # Stagnated at round-off: a trace species in fast equilibrium (H⁺ ≈ 4e-12 at n = 6.5e16, T = 2850 K) flips in its 10th digit
        # from one iteration to the next while T and all other abundances are converged to the last digit; the test above (1e-10 per
        # species) is then never met for any h, and the halving ended in "Newton failed" at depth 12. Accepted only after maxit
        # full Newton steps, so steps that converged before are unchanged bit for bit.
        λ == 1.0 && rel < 1e-6 && dT < 1e-8 && return (z, true)
        (z, false)
    end
    function adv(x0, T0, h, depth)
        CHEMSTATS && Threads.atomic_add!(NADV, 1)
        z, ok = attempt(x0, T0, h)
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
