"""
    PopIIIMaterial

Population III project, M2c: metal-free gas with non-equilibrium chemistry as a TITAN material. Composition per zone
(`y[k]` = (x_H⁺, x_H⁻, x_H₂, x_H₂⁺), PrimordialChem), RHFC02 EOS (`eos_rhfc02`: e_g includes the chemical binding energy).
The hydro step (TITAN's Newton) sees the composition frozen; after each converged step `chem!` advances chemistry,
cooling and T together in every zone at fixed (ρ, e_g) over the step just taken (`PrimordialChem.step_coupled`,
backward Euler) — operator splitting, option (A) of M0.2. Needs `TitanResidFast.ZONE` (zone index of a material call).
Hydrodynamics only for now (no opacities: neqn = 5).
"""
module PopIIIMaterial

import ..TitanDriver, ..TitanResidFast, ..TitanMaterials, ..TitanCore, ..PrimordialChem, ..PopIIIEqEOS
using ForwardDiff, Printf
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const TM = TitanMaterials; const TC = TitanCore
const PC = PrimordialChem; const EQ = PopIIIEqEOS

mutable struct PrimordialGas <: TM.Material
    y::Vector{NTuple{4,Float64}}          # per zone index k (length mgr)
    rates::Function                       # T -> named tuple of rate coefficients
    cooling::Union{Nothing,Function}      # (ρ, T, y) -> Λ [erg cm⁻³ s⁻¹]
    ρeq::Float64                          # equilibrium chemistry above this density (RHFC02: 5e-11 g cm⁻³); Inf = never
    eqtable::Bool                         # true: equilibrium EOS from PopIIIEqEOS (Saha + SCVH95, RHFC02); false: Saha only
    Λext::Vector{Float64}                 # per zone k: extra cooling [erg cm⁻³ s⁻¹] (line transfer), evaluated at T = Tref
    Tref::Vector{Float64}                 # G7: over the next step Λext × (T/Tref)^slope, implicit in T (see ext_cooling)
    slope::Vector{Float64}                #     slope = d ln Λ/d ln T (0: Λext held fixed, explicit)
end
PrimordialGas(y, rates, cooling) = PrimordialGas(y, rates, cooling, 5e-11, false)
PrimordialGas(y, rates, cooling, ρeq) = PrimordialGas(y, rates, cooling, ρeq, false)
PrimordialGas(y, rates, cooling, ρeq, eqtable) = PrimordialGas(y, rates, cooling, ρeq, eqtable, zeros(length(y)), ones(length(y)), zeros(length(y)))

"""
G7: the transfer cooling Λ_k of zone k as a function of its own temperature over the step, Λ_k (T/T_ref)^s with
s = d ln Λ_k/d ln T from the transfer (all zones perturbed together: the local response including the coupling to the
neighbours). Inside the implicit chemistry–temperature step this makes the line cooling implicit (backward Euler in T,
linearised in ln T) instead of explicit with Δt ≤ 0.1 uρ/Λ.
"""
ext_cooling(m, k) = (Λx = m.Λext[k]; Λx == 0 ? nothing :
                     m.slope[k] == 0 ? ((T) -> Λx + zero(T)) : (Tr = m.Tref[k]; s = m.slope[k]; (T) -> Λx * (T / Tr)^s))
"Cooling function for zone k: the material's own plus the transfer term."
function zone_cooling(m, k)
    f = ext_cooling(m, k)
    c = f === nothing ? m.cooling : m.cooling === nothing ? ((ρ, T, y) -> f(T)) : ((ρ, T, y) -> m.cooling(ρ, T, y) + f(T))
    (HDCHEM[] && k <= length(HDX[]) && HDX[][k] > 0) || return c
    xhd = HDX[][k]                                   # HDCHEM=1: HD cooling with the x_HD of the start of the step
    c === nothing ? ((ρ, T, y) -> hd_cool(ρ, T, y, xhd)) : ((ρ, T, y) -> c(ρ, T, y) + hd_cool(ρ, T, y, xhd))
end

@inline function TX.material_gas(m::PrimordialGas, ρ, T)
    if m.eqtable && ForwardDiff.value(ρ) >= m.ρeq
        u, p = EQ.eos(ρ, T)                 # composition in equilibrium at (ρ, T): differentiated through by the Newton
        return p, u, sqrt(5.0 / 3.0 * p / ρ)
    end
    y = m.y[TX.ZONE[]]
    u, p = PC.eos_rhfc02(ρ, T, y)
    # a_s enters only the artificial viscosity: γ ≈ 5/3 estimate
    return p, u, sqrt(5.0 / 3.0 * p / ρ)
end

function TM.fill_eos!(m::PrimordialGas)
    dn, tn, pgn, egn, asn = TM.v(:dn), TM.v(:tn), TM.v(:pgn), TM.v(:egn), TM.v(:asn)
    @inbounds for k in TM.krange()
        (dn[k] > 0 && tn[k] > 0) || continue
        TX.ZONE[] = k
        pg, eg, as_ = TX.material_gas(m, dn[k], tn[k])
        pgn[k] = pg; egn[k] = eg; asn[k] = as_
    end
    nothing
end
TM.fill_opac!(::PrimordialGas) = nothing

"Copy the composition of the first/last real zone into the ghost zones."
function ghosts!(m::PrimordialGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    for k in a-2:a-1; k >= 1 && (m.y[k] = m.y[a]); end
    for k in b+1:b+2; m.y[k] = m.y[b]; end
end

"""
After a converged step: chemistry + cooling + T in every zone over the step just taken, at fixed ρ and e_g; then the
EOS fill (p, e consistent with the new T and y). Returns the largest relative change of T (for diagnostics).
"""
function chem!(m::PrimordialGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    dn, tn, egn = TM.v(:dn), TM.v(:tn), TM.v(:egn)
    dTmax = 0.0
    for k in a:b
        cool = zone_cooling(m, k)
        y1, T1 = dn[k] >= m.ρeq ? (m.eqtable ? EQ.step(dn[k], egn[k], tn[k], dt; cooling = cool) :
                                             PC.step_equilibrium(dn[k], egn[k], tn[k], dt; cooling = cool)) :
                                  PC.step_coupled(m.y[k], dn[k], egn[k], tn[k], dt, m.rates; cooling = cool)
        dTmax = max(dTmax, abs(T1 / tn[k] - 1))
        m.y[k] = y1; tn[k] = T1
    end
    ghosts!(m)
    # ghost temperatures as TITAN's boundary rules (reflecting centre, copy outside), then the EOS fill
    tn[a-1] = tn[a]; tn[b+1] = tn[b]; tn[b+2] = tn[b]
    TM.fill_eos!(m)
    dTmax
end


# ====================================================================== G5: network gas (ChemNet) with an EOS-validity switch
"""
    NetGas(x, net, photos, cooling; env, pswitch = 1e4)

Metal-free gas with the general network (ChemNet/ChemRates; composition `x[k]` = free-species abundances per zone). The
switch to the equilibrium EOS table (PopIIIEqEOS: Saha + SCVH95) is made on **EOS validity, not density**: where the
ideal-gas pressure of the zone is ≥ `pswitch` (RHFC02: non-ideal effects from 1e4 dyn cm⁻², SCVH95 from 1e6 with a
blend between) — the network (detailed balance, implicit, L-stable) already relaxes to equilibrium by itself wherever
the chemistry is fast, so no ρ_eq is needed. `env` supplies J_ν for photo-processes (`ChemNet.photorates`); default
the CMB/zero field given by the caller.
"""
mutable struct NetGas <: TM.Material
    x::Vector{Vector{Float64}}           # per zone index k
    net::Any; photos::Any
    cooling::Union{Nothing,Function}     # (ρ, T, d::Dict) -> Λ
    env::Any                             # (Trad, J)
    pswitch::Float64
    eqzone::Vector{Bool}                 # per zone: in the equilibrium (table) regime after the last chem! step
    Λext::Vector{Float64}; Tref::Vector{Float64}; slope::Vector{Float64}     # as PrimordialGas (G7)
end
NetGas(x, net, photos, cooling; env, pswitch = 1e4) =
    NetGas(x, net, photos, cooling, env, pswitch, falses(length(x)), zeros(length(x)), ones(length(x)), zeros(length(x)))

import ..ChemNet
const CNn = ChemNet

"""
S6-C1: hot branch of `NetGas` for the protostar (popiii_m2c HOT=1). `HOT[] = (; T1, T2, Δ0, gas, thermo)` with
`gas(ρ, T, k) -> (p, e)` the stellar EOS of zone k's composition, `thermo(ρ, T, k) -> (∇_ad, c_P, δ)`. Between T1 and T2 p and e are
blended with a C² weight in ln T, e_gas = e_table + w (e_star − Δ0 − e_table), Δ0 the zero-point difference of the two EOS
(neutral atoms vs PopIIIEqEOS's zero), measured in the fully ionised dilute limit. The blend is not thermodynamically consistent
inside the band, but e(ρ, T) and p(ρ, T) are continuous and differentiable, so a zone crossing it conserves energy.
`nothing` (default): no change, the runs below T1 are identical.
"""
const HOT = Ref{Any}(nothing)
@inline function hotw(T, T1, T2)
    x = log(T / T1) / log(T2 / T1)
    x * x * x * (x * (6x - 15) + 10)
end
@inline function TX.material_gas(m::NetGas, ρ, T)
    k = TX.ZONE[]
    h = HOT[]
    if h !== nothing && ForwardDiff.value(T) > h.T1
        ps, es = h.gas(ρ, T, k); es -= h.Δ0
        if ForwardDiff.value(T) >= h.T2
            return ps, es, sqrt(5.0 / 3.0 * ps / ρ)
        end
        u, p = EQ.eos(ρ, T); w = hotw(T, h.T1, h.T2)
        pb = p + w * (ps - p)
        return pb, u + w * (es - u), sqrt(5.0 / 3.0 * pb / ρ)
    end
    if m.eqzone[k]
        u, p = EQ.eos(ρ, T)
    else
        u, p = CNn.eos(ρ, T, m.x[k], m.net)
    end
    return p, u, sqrt(5.0 / 3.0 * p / ρ)
end
"""
S6: nuclear source of the hot branch in the implicit energy row (S5's coupling: `h.eps(ρ, T, k)` = level factor × ε(ρ, T; Y lagged), generic in the
number type, so ∂ε/∂T is in the Jacobian; 0 below T1).
"""
TX.has_source(::NetGas) = HOT[] !== nothing
TX.material_source(::NetGas, ρ, T) = (h = HOT[]; ForwardDiff.value(T) > h.T1 ? h.eps(ρ, T, TX.ZONE[]) : zero(T))
"""
S6-C1d (experiment): MLT also on the table EOS, in equilibrium zones with T > `COLDMLT[]` (0 = off, the default). `cold_thermo` is ∇_ad, c_P, δ of
PopIIIEqEOS + equilibrium radiation, from derivatives of the table, generic in the number type.
"""
const COLDMLT = Ref(0.0)
function cold_thermo(ρ, T)
    arad = 7.565723e-15
    pf(r, t) = EQ.eos(r, t)[2] + arad * t^4 / 3
    ef(r, t) = EQ.eos(r, t)[1] + arad * t^4 / r
    p = pf(ρ, T)
    pT = ForwardDiff.derivative(t -> pf(ρ, t), T); pρ = ForwardDiff.derivative(r -> pf(r, T), ρ)
    eT = ForwardDiff.derivative(t -> ef(ρ, t), T); eρ = ForwardDiff.derivative(r -> ef(r, T), ρ)
    cP = eT - (eρ - p / ρ^2) * pT / pρ
    δ = T / ρ * pT / pρ
    return p * δ / (T * ρ * cP), cP, δ
end
"""MLT thermodynamics: only the hot branch convects (and, with `COLDMLT[]` > 0, the hot equilibrium zones of the table EOS). Below T1 the huge ∇_ad switches the flux off in the chemistry gas: it must be
out of reach (it was 10, and at a T jump inside the protostar ∇ = Δln T/Δln P crossed it: an MLT flux with this placeholder c_P
of up to 1e40 erg/s between two zones, and Newton diverged)."""
function TX.material_thermo(m::NetGas, ρ, T)
    h = HOT[]
    (h !== nothing && T > h.T1) && return h.thermo(ρ, T, TX.ZONE[])
    (COLDMLT[] > 0 && T > COLDMLT[] && m.eqzone[TX.ZONE[]]) && return cold_thermo(ρ, T)
    return (1e30, 1e9, 1.0)
end
function TM.fill_eos!(m::NetGas)
    dn, tn, pgn, egn, asn = TM.v(:dn), TM.v(:tn), TM.v(:pgn), TM.v(:egn), TM.v(:asn)
    @inbounds for k in TM.krange()
        (dn[k] > 0 && tn[k] > 0) || continue
        TX.ZONE[] = k
        pg, eg, as_ = TX.material_gas(m, dn[k], tn[k])
        pgn[k] = pg; egn[k] = eg; asn[k] = as_
    end
    nothing
end
"""
S6: continuum opacity of `NetGas` when TITAN's radiation is on (popiii_m2c RAD=1): `OPAC[] = (ρ, T) -> (κ_R, κ_P)`
[cm² g⁻¹], e.g. StellarOpacity (own CIA/H⁻ means in chemical equilibrium below 10^3.8 K, OPAL above). `nothing` (default):
hydrodynamics only, no opacity calls. The H₂ lines stay outside the continuum (H2Lines, as Λext).
"""
const OPAC = Ref{Any}(nothing)
"""
Non-LTE opacity of the network zones of a three-level run (SL3=1, `NLOPAC[] = true`). `OPAC[]` is a table in chemical *and* excitation
equilibrium: its electrons are Saha's at the gas T and its H levels Boltzmann's. Outside the photosphere neither holds. At 36 M☉ just outside
τ = 1 (T ≈ 5700 K, ρ ≈ 2e-12) 83 % of κ_P is H bound–free (n = 1 and n ≥ 2 from Boltzmann populations ∝ e^{−10.2 eV/kT}), so κ_P ∝ T²²
while the network zone's heat capacity at frozen composition is 9× below the table's: Newton solved an emission ∝ T²⁶ there and diverged
every few steps (Δt 1e3–5e4 s, `_sl3b`). In the three-level atom n = 2 is set by recombination and the radiation field (4 n_e q21 ≪ β):
it hardly depends on the gas T at fixed composition, and its bound–free energy exchange is already in the chemistry (tl3_heat, GASQ/RAYQ).

So in every network zone with three-level parameters, from the end of `chem!` (T0, final composition), for the next hydro step
(TITAN's exchange term is c κ_E E − 4π κ_P B):
- κ_E (absorption): continuum Planck mean of the zone's own composition: the network's H⁻ (photodetachment, associative detachment),
  H⁻ ff, H ff, CIA, no H bound–free (case B: Lyman continuum in detailed balance; n ≥ 2 in the three-level atom);
- κ_P (emission): the same with H⁻ bound–free from H⁻ in Saha with the zone's own H and electrons (Kirchhoff–Milne for radiative
  attachment with Maxwellian electrons: emission ∝ n_H n_e, whatever n_H⁻ is). With the network's E0 the gas keeps −kT per attachment;
  the table, with equilibrium (molecular, electron-poor) gas, had no such emission in the cold infall while the network released 0.754 eV
  into heat per attachment;
- κ_R: Rosseland mean of the extinction with the network's H⁻ and H bound–free from the non-LTE n = 2 of the three-level atom (n ≥ 3
  with the same departure coefficient as n = 2: an approximation, these levels matter only through the Paschen continuum), Thomson on the
  zone's electrons;
- each as a ratio to the table at (ρ0, T0) below 6000 K (where the table is the same continuum code in equilibrium composition, agreement
  ≤ 3 %), above to the continuum code in equilibrium composition (so the table's calibration, OPAL and lines above 6300 K, is kept); then
  κ = κ_table(ρ, T0) × ratio × (T/T0)^q with q the frozen-composition slope of the continuum code, so Newton sees frozen composition in
  the EOS and in the opacity alike. No threshold in T: a T switch of κ made zones bistable (sawtooth 2900/4000 K at 3000–3500 K).
Deep in (τ ≳ 1, W = 1, T_c = T) the three-level n = 2 is Boltzmann's (b2 ≈ 0.8–1.1 at 36 M☉), the electrons are Saha's, and only the
dropped H bound–free of κ_P differs from the table zones: at zone 129 (τ 0.2) the three-level bound–free heating, 2.1e-5 erg cm⁻³ s⁻¹,
and the LTE κ_P,Hbf (B − J) it replaces, 5.2e-5, agree in sign and order; further out the diluted 6000 K starlight heats by n = 2
photoionization where the gray LTE term (T > T_r) cooled.
"""
const NLOPAC = Ref(false)
const NLOPC = Ref{Vector{Any}}(Any[])          # per zone: (; T0, RE, qE, RP, qP, RR, qR) or nothing
"""
NLODAMP=w (0 < w < 1; default 1 = off): under-relaxation in time of the lagged non-LTE opacity factors of the network zones (NLOPAC),
ln R ← (1 − w) ln R_prev + w ln R_new for RE, RP, RR, where the zone had factors in the previous step too. The factors are frozen inside
Newton and set by the chemistry after it; at 6.5 M_sun the pre-shock gas (T 4000–6000 K, x_e ≈ Saha, d ln x_e/d ln T ≈ 16) flips between
a hot opaque and a cool transparent state every 20–40 steps through that one-step lag. Test switch.
"""
const NLODAMP = parse(Float64, get(ENV, "NLODAMP", "1"))
const NLOPREV = Ref{Vector{Any}}(Any[])
import ..ContinuumOpacity
const COp = ContinuumOpacity
const NLGRID = Ref{Vector{Float64}}(Float64[])
"""
Frequency-only quantities of the nlte_means grid (midpoint ν, dν, σ = ν/c, λ in μm, ν³, σ_H⁻bf(λ), Rayleigh σ of H₂ and H), computed once:
the same functions of the same arguments as in the loop, so the means are bit-identical; they were recomputed 2–3 × 1000 times per
network zone per step.
"""
const NLFREQ = Ref{Any}(nothing)
function nlfreq()
    NLFREQ[] === nothing || return NLFREQ[]
    isempty(NLGRID[]) && (NLGRID[] = COp.grid(; n = 1000))
    νs = NLGRID[]; n = length(νs) - 1
    ν = [0.5 * (νs[i] + νs[i+1]) for i in 1:n]; dν = [νs[i+1] - νs[i] for i in 1:n]; σ = [x / COp.c for x in ν]; λμ = [1e4 / x for x in σ]
    NLFREQ[] = (; ν, dν, σ, λμ, ν3 = [x^3 for x in ν], σhm = [COp.σ_Hm_bf(x) for x in λμ], ray2 = [COp.σ_ray_H2(x * 1e3) for x in λμ],
                ray1 = [COp.σ_ray_H(x * 1e4) for x in λμ])
end
"H bound–free absorption [cm⁻¹] (Kramers, as ContinuumOpacity.α_H_kramers) from n1 and the population of n = 2; n ≥ 3 Boltzmann relative to n = 2."
function hbf_nlte(ν, T, n1, n2)
    χ = 13.598 * EVERG; σ0 = 7.906e-18; hν = COp.h * ν; α = 0.0
    hν >= χ && (α += n1 * σ0 * (χ / hν)^3)
    for n in 2:30
        hν >= χ / n^2 || continue
        nn = n2 * (n^2 / 4) * exp(-χ * (0.25 - 1 / n^2) / (CNn.kB * T))
        α += nn * σ0 * n * (χ / n^2 / hν)^3
    end
    α
end
"The same with the n = 2…30 populations precomputed (`nn[n-1]`, independent of ν): bit-identical to `hbf_nlte`, without 29 exp per frequency."
function hbf_nlte_pre(ν, n1, nn)
    χ = 13.598 * EVERG; σ0 = 7.906e-18; hν = COp.h * ν; α = 0.0
    hν >= χ && (α += n1 * σ0 * (χ / hν)^3)
    for n in 2:30
        hν >= χ / n^2 || continue
        α += nn[n-1] * σ0 * n * (χ / n^2 / hν)^3
    end
    α
end
hbf_levels(T, n2) = ntuple(i -> (n = i + 1; n2 * (n^2 / 4) * exp(-13.598 * EVERG * (0.25 - 1 / n^2) / (CNn.kB * T))), 29)
"""Planck means of absorption (κ_E, H⁻ bf from nHm) and emission (κ_P, H⁻ bf from nHmS), H bound–free in both or neither (`hbfP`), and
Rosseland mean of the extinction (H⁻ from nHm, H bound–free always) [cm² g⁻¹]; number densities in cm⁻³."""
function nlte_means(ρ, T; n1, n2, nH2, nHe, nHp, ne, nHm, nHmS, hbfP, Tw = T)
    F = nlfreq(); t = COp.tables(); sE = 0.0; sP = 0.0; sB = 0.0; sR = 0.0; sD = 0.0; sW = 0.0; wJ = Tw != T
    cia = nH2 > 1e-3 * n1; acia = 0.0; nnl = hbf_levels(T, n2)
    @inbounds for i in eachindex(F.ν)
        ν = F.ν[i]; dν = F.dν[i]; σ = F.σ[i]; λμ = F.λμ[i]
        stim = -expm1(-COp.h * ν / (COp.kB * T))
        # CIA bands are broad (≳ 100 cm⁻¹; the grid step is 1 %): evaluated on every 8th point, held in between (cost: 30 ms per call otherwise)
        cia && i % 8 == 1 && (acia = COp.k_h2h2(t, σ, T) * nH2^2 + COp.kcia(t.H2He, σ, T) * nH2 * nHe + COp.kcia(t.H2H, σ, T) * nH2 * n1)
        a0 = acia + COp.k_Hm_ff(λμ, T) * n1 * ne * COp.kB * T + 3.692e8 * ne * nHp / (sqrt(T) * F.ν3[i]) * stim
        σhm = F.σhm[i] * stim
        ab = hbf_nlte_pre(ν, n1, nnl) * stim
        sca = F.ray2[i] * nH2 + F.ray1[i] * n1 + COp.σT * ne
        B = COp.Bν(ν, T); D = COp.dBdT(ν, T); h = hbfP ? ab : 0.0
        BJ = wJ ? COp.Bν(ν, Tw) : B
        sE += (a0 + σhm * nHm + h) * BJ * dν; sW += BJ * dν; sP += (a0 + σhm * nHmS + h) * B * dν; sB += B * dν
        sR += D / max(a0 + σhm * nHm + ab + sca, 1e-300) * dν; sD += D * dν
    end
    (sE / sW / ρ, sP / sB / ρ, sD / sR / ρ)
end
"n = 2 population [cm⁻³] of the three-level atom (2s + 2p, 4 n_2s) in steady state for the zone's n1, n_e, n_p and lagged radiation p3."
function n2_tl3(T, n1, ne, np, p)
    D = p.RL + LAM2S + 4ne * q21_H(T)
    4 * (n1 * (ne * q12_H(T) + p.P12) + alphaB_recfast(T) * ne * np) / (D + p.β)
end
"(κ_E, κ_P, κ_R) of the continuum code for zone composition d (full species fractions) or, with `d = nothing`, equilibrium composition."
function nlte_kappas(ρ, T, d, p)
    nH = ρ / CNn.μH(); nHe = PC.yHe * nH
    if d === nothing
        y = EQ.composition(ρ, T); xp = y[1]; xh2 = y[3]; n1 = nH * max(1 - xp - 2xh2, 0.0); ne = nH * xp
        if CNn.HECHEM                                 # helium electrons; neutral He for the CIA and Rayleigh terms
            z1, z2 = EQ.composition_he(ρ, T); ne = nH * (xp + z1 + 2z2); nHe = nH * max(PC.yHe - z1 - z2, 0.0)
        end
        n2 = 4n1 * exp(-E21 / (CNn.kB * T)); nHm = n1 * ne / CNn.Keq([:Hm], [:H, :e], T)
        return nlte_means(ρ, T; n1, n2, nH2 = nH * xh2, nHe, nHp = ne, ne, nHm, nHmS = nHm, hbfP = true)
    end
    n1 = nH * max(d[:H], 0.0); ne = nH * max(d[:e], 0.0); np = nH * max(d[:Hp], 0.0)
    CNn.HECHEM && (nHe = nH * max(d[:He], 0.0))
    nHmS = n1 * ne / CNn.Keq([:Hm], [:H, :e], T)
    nHm = haskey(d, :Hm) ? nH * max(d[:Hm], 0.0) : nHmS
    Tw = NLOPACJ[] && hasproperty(p, :Tc) ? p.Tc : T
    nlte_means(ρ, T; n1, n2 = n2_tl3(T, n1, ne, np, p), nH2 = nH * max(d[:H2], 0.0), nHe, nHp = np, ne, nHm, nHmS, hbfP = false, Tw)
end
"Lagged non-LTE opacity factors of zone k at T0 (see NLOPAC)."
function nlopc(ρ, T0, d, p)
    OPAC[] === nothing && return nothing     # no opacity table yet (the cold collapse, before RAD's opacity is set): LTE tables later
    f = 1.02
    E0, P0, R0 = nlte_kappas(ρ, T0, d, p); E1, P1, R1 = nlte_kappas(ρ, T0 * f, d, p)
    if T0 < 6000
        Rt, Pt = OPAC[](ρ, T0); Pe, Re = Pt, Rt
    else
        _, Pe, Re = nlte_kappas(ρ, T0, nothing, p)
    end
    lf = log(f)
    (; T0, RE = E0 / Pe, qE = log(E1 / E0) / lf, RP = P0 / Pe, qP = log(P1 / P0) / lf, RR = R0 / Re, qR = log(R1 / R0) / lf)
end
"""
Diagnostic switch, off by default: `CHICAP=κ[,Tmax]` caps the flux-mean opacity χ at κ [cm²/g] where T < Tmax (default 5e4 K). The energy
means κ_E, κ_P are not changed. With κ = electron scattering it removes the H/H⁻ opacity peak from the radiation force and the diffusion,
to test what drives the envelope cycle at L ≈ L_Edd. Not physics.
"""
const CHICAP = let s = split(get(ENV, "CHICAP", "0"), ","); (parse(Float64, s[1]), length(s) > 1 ? parse(Float64, s[2]) : 5e4); end
@inline chicap(κ, T) = (CHICAP[1] > 0 && T < CHICAP[2]) ? min(κ, CHICAP[1]) : κ
@inline function TX.material_opac(m::NetGas, ρ, T)
    if NLOPAC[]
        k = TX.ZONE[]; c = NLOPC[]
        if k <= length(c) && !m.eqzone[k] && (o = c[k]) !== nothing
            κR0, κP0 = OPAC[](ρ, o.T0); s = T / o.T0
            return chicap(κR0 * o.RR * s^o.qR, T), κP0 * o.RE * s^o.qE, κP0 * o.RP * s^o.qP
        end
    end
    κR, κP = OPAC[](ρ, T)
    return chicap(κR, T), κP, κP
end
function TM.fill_opac!(m::NetGas)
    OPAC[] === nothing && return nothing
    dn, tn = TM.v(:dn), TM.v(:tn)
    chifn, xken, xkpn, xnen = TM.v(:chifn), TM.v(:xken), TM.v(:xkpn), TM.v(:xnen)
    @inbounds for k in TM.krange()
        (dn[k] > 0 && tn[k] > 0) || continue
        TX.ZONE[] = k
        χ, κE, κP = TX.material_opac(m, dn[k], tn[k])
        chifn[k] = χ; xken[k] = κE; xkpn[k] = κP
        # electrons for the scattering term of the formal solution (the template state's xnen is not ours)
        xnen[k] = (CNn.HECHEM ? (m.eqzone[k] ? (Tq = min(tn[k], 1e6); z = EQ.composition_he(dn[k], Tq); EQ.composition(dn[k], Tq)[1] + z[1] + 2z[2]) :
                                 CNn.full(m.x[k], m.net)[:e]) :
                   (m.eqzone[k] ? EQ.composition(dn[k], min(tn[k], 1e6))[1] : CNn.full(m.x[k], m.net)[:Hp])) * dn[k] / CNn.μH()
    end
    nothing
end

"Free-species vector from the equilibrium table's composition (H₂, H⁺; H⁻ and Li⁺ in Saha with the table's electrons)."
function x_from_table(m::NetGas, ρ, T)
    y = EQ.composition(ρ, T); nH = ρ / CNn.μH()
    xH2 = y[3]; xHp = y[1]; xH = 1 - xHp - 2xH2
    ne = xHp * nH
    z1, z2 = CNn.HECHEM ? EQ.composition_he(ρ, T) : (0.0, 0.0)
    CNn.HECHEM && (ne = (xHp + z1 + 2z2) * nH)
    xHm = xH * ne / CNn.Keq([:Hm], [:H, :e], T)
    out = zeros(length(m.net.free))
    for (i, s) in enumerate(m.net.free)
        out[i] = s === :Hp ? xHp : s === :H2 ? xH2 : s === :Hm ? xHm : s === :Hep ? z1 : s === :Hepp ? z2 :
                 s === :Lip ? (K = CNn.Keq([:Li], [:Lip, :e], T); CNn.yLi * K / (K + max(ne, 1e-300))) : 0.0
    end
    out
end

"""
S6: hysteresis of the EOS switch: a zone in the table regime stays there until p < PSWHYS × pswitch (1 = the sharp switch of G5). Warm gas
of the accretion layer sitting on pswitch otherwise flips every step between table and network, which differ by ~10 % in T at fixed e.
"""
const PSWHYS = Ref(1.0)
const PSWTMIN = Ref(0.0)
"""
Table to network, falling pressure: PSWTEX[] > 0 keeps a table zone on the table until T < PSWTEX[], whatever its pressure. The network has no photoionization by the local radiation field, so
at 8500 K it holds x(H+) well below Saha at fixed e, and the recombination energy went into heat: T 8509 -> 13938 K in one step at the surface of a 28 M_sun protostar, then Newton failed at every dt.
E_ion/E_th at rho = 1e-11..1e-10 is <= 4 % at 5000 K, 15-48 % at 6000 K. Zones the ray reaches go to the network before this test.
"""
const PSWTEX = Ref(0.0)
"""
PSWSTAY[] = false drops the PSWTEX "stay" rule above, so a table zone below the pressure switch stays on the table only by the τ rule (PSWTAU).
With the three-level network (SL3) the network carries the star-light ionization the rule was written to stand in for; the rule otherwise holds table
zones out to τ ≈ 0.4 while T ≥ PSWTEX, so the boundary becomes a temperature switch in the zones where table and network disagree most.
"""
const PSWSTAY = Ref(true)
"""
PSWTAU[] > 0: a zone with T >= PSWTEX[] and Rosseland depth from the surface >= PSWTAU[] (half that to stay) uses the table whatever its pressure. Where radiation is thermalised, J = B holds H at Saha.
The network has no radiation field: at the outer edge of an opaque layer (p ~ 6, T ~ 1e4 K, 39 M_sun) it sits at coronal x(H+) ~ 0.05 against 0.99, Newton heats the zone by a factor of 2 at frozen
composition and chem! brings it back each step, until Newton fails. On the table the ionization heat capacity is inside Newton. The switch at fixed e lowers T once; radiation reheats the zone.
"""
const PSWTAU = Ref(0.0)
"""
PSWAGREE[] = ε > 0 (three-level runs): below the pressure switch the table is used only where it agrees with the network. The τ rule (PSWTAU)
flipped bands of zones in the extended photosphere of the non-LTE opacity (τ 1 → 0.4 over nine zones): a zone left the table at τ < 0.5, the
three-level atom ionized it in one step (x_e 3e-3 → 4.5e-2: the diluted 6000 K star light is hotter than the gas), its κ took τ back over 1,
the table reset x_e to Saha, T jumped, and Δt fell to ~1e5 s. Now, at T ≥ PSWTEX[] and with three-level parameters for the zone:
- a network zone goes to the table only if the table's T at its u is within ε of the network T, *and* relaxing the table composition to the
  three-level steady state (x_ss, H only, `xss_tl3`) at fixed u would move T by less than ε (it would not leave again at once);
- a table zone stays while x_ss moves T by less than 2ε (hysteresis); leaving is seamless (the network starts from the table composition);
- at τ ≥ PSWTAUDEEP[] the table always (thermalised radiation, see there).
The switch then happens only where both models give the same state to ε in T. The pressure rule (EOS validity) is unchanged.
"""
const PSWAGREE = Ref(0.0)
# At τ ≥ PSWTAUDEEP[] (3, where SLBLEND reaches T_c = T_r, W = 1) the table always: there a network zone, whose n = 2 bound–free coupling to
# the radiation is in the explicit three-level terms (κ_P without H bound–free), fell 1800 K below T_r at τ ≈ 20, so x_ss ≠ Saha(T) kept it off the
# table and the disagreement ran inward (`_nlo4`, first try).
const PSWTAUDEEP = Ref(3.0)
"SLBLEND[]: star-light colour temperature and dilution blended smoothly in τ (see starlight_setup); false: switch at the τ = 2/3 zone."
const SLBLEND = Ref(false)
"""
SLEB[] (with SLBLEND): the star-light colour temperature of zone k is local, T_r at τ_k + 2/3 (Eddington–Barbier: the outward light at depth τ_k
comes from about Δτ = 2/3 further in), interpolated in log τ, in place of the photosphere's T_eff. In the extended precursor τ is nearly flat
(0.8 → 0.4 over 400–1400 R_sun at 36 M_sun), the τ = 2/3 zone hopped between 122 and 130 from step to step, and T_eff = (L/4πr²σ)^¼ there
swung 5580 ↔ 7100 K at constant L: every lit zone's spectrum jumped by 1500 K, faster than its ionization relaxes (C α_B n_e ~ 1e-5 s⁻¹).
"""
const SLEB = Ref(false)
"""
NLOPACJ[] (with NLOPAC): κ_E of a network zone, the absorption mean of the gray E_r, weighted with the star light's spectrum B_ν(T_c) (T_c of
the three-level parameters) in place of B_ν(T_gas). In the precursor T_gas ≈ 4800 K under T_c ≈ 6000 K light; H⁻ bound–free stops at 1.64 μm,
where B(4800 K) has more of its energy than B(6000 K).
"""
const NLOPACJ = Ref(false)
const NSWIN = Ref(0); const NSWOUT = Ref(0)
"""
PSWRAD=ε[,phard] (default off): a network zone ENTERS the table only if its gas and radiation temperatures agree, |T_rad/T − 1| ≤ ε with
T_rad = (E_r/a)^¼, unless its pressure is at least phard × pswitch (default 1: only the entries below the pressure switch, i.e. by the τ rule
or PSWAGREE, are held to the test). Staying on the table and leaving it are unchanged.
Why: the table is LTE, composition and κ_P at Saha of T_gas, which is the state of the gas only where the radiation is thermalised at T_gas.
The τ rule assumes that from τ ≥ PSWTAU alone. In a dense shell piled up at the photosphere a zone is optically thick by itself while shock
or compression heat still holds T_gas above T_rad: production zone 144 at 101.3 M☉ entered the table at T 7600 K with T_rad 5416 K; the LTE
κ_P of 9 cm²/g then equilibrates gas and radiation in 0.5 s, the step had no solution near the old state for any Δt above that, and the run
stopped on the retry limit (docs/POPIII_RUNS.md). Held on the network the zone couples through the opacity of its actual composition
(NLOPAC), n_e and κ_P rising together, and joins the table once the two temperatures have met, where the three-level steady state (W = 1,
T_c = T_r at τ ≥ 3) is Saha.
"""
const PSWRAD = Ref{Union{Nothing,NTuple{2,Float64}}}(let s = get(ENV, "PSWRAD", "0"); s == "0" ? nothing : (v = parse.(Float64, split(s, ",")); (v[1], length(v) > 1 ? v[2] : 1.0)) end)
const NSWRAD = Ref(0)
"H⁺ fraction of the three-level atom in steady state at (ρ, T) for the lagged radiation p (H only, n_e = n_H⁺)."
function xss_tl3(ρ, T, p)
    nH = ρ / CNn.μH(); q = q12_H(T); a = alphaB_recfast(T)
    f(y) = (ne = nH * y; C = tl3_C(T, (; RL = p.RL, ne, β = p.β)); (1 - C) * (1 - y) * (ne * q + p.P12) - C * a * ne * y)
    lo, hi = 1e-14, 1 - 1e-12
    for _ in 1:200; mid = sqrt(lo * hi); f(mid) > 0 ? (lo = mid) : (hi = mid); hi / lo - 1 < 1e-6 && break; end
    sqrt(lo * hi)
end
"""T at fixed (ρ, u) of the three-level atom in steady state: u(ρ, T, x(T)) = u with x(T) the table composition and H⁺ = x_ss(T)
(the three-level model's own equilibrium, ionization heat capacity included; bisection in ln T), NaN without a bracket."""
function T_ss_at_u(m::NetGas, ρ, u, p, T0, iHp)
    xT(T) = (x = x_from_table(m, ρ, T); x[iHp] = xss_tl3(ρ, T, p); x)
    T_at_u(ρ, u, T -> xT(T), m.net, T0)
end
"T at fixed (ρ, u) for network composition x (a vector, or a function of T), bisection in ln T, NaN without a bracket."
function T_at_u(ρ, u, x, net, T0)
    F(T) = CNn.eos(ρ, T, x isa Function ? x(T) : x, net)[1] - u
    lo, hi = T0 / 1.5, T0 * 1.5
    while F(lo) > 0 && lo > 20; lo /= 1.5; end
    while F(hi) < 0 && hi < 1e7; hi *= 1.5; end
    (F(lo) > 0 || F(hi) < 0) && return NaN
    for _ in 1:100; mid = sqrt(lo * hi); F(mid) > 0 ? (hi = mid) : (lo = mid); hi / lo - 1 < 1e-9 && break; end
    sqrt(lo * hi)
end
"After a converged hydro step: network step per zone in the ideal regime, equilibrium-table step where p ≥ pswitch."
const CHEMFALLBACK = Ref(false); const NCHEMFB = Ref(0)
"""
CHEMFBDENSE=1 (default): a failed network step of a zone whose pressure is above the table switch (p ≥ pswitch; only its T < PSWTMIN keeps
it on the network) takes the table-equilibrium fallback too, in every phase, counted in NCHEMFB. Such gas is in chemical equilibrium to a
good approximation (n ≳ 1e15 cm⁻³), and it is where the backward-Euler Newton has failed in the collapse runs: the H₂ dissociation layer
just below PSWTMIN at Δt of hundreds of seconds (full3, core growth stage, step 11935: n_H 4.5e18, T 4346 K, Δt 515 s, twelve halvings
without convergence), which without the fallback ends the run. 0: the fallback only where CHEMFALLBACK is set (the accretion phase).
"""
const CHEMFBDENSE = Ref(get(ENV, "CHEMFBDENSE", "1") == "1")
"""
M3: the star's ionizing light (IonRay), `nothing` = off. `RAY[] = (; start, source, env, Ec, Q, xmin)`: `start(a, b)` -> (photosphere zone, photons per group) or `nothing`; `source` as IonRay.source;
`env`: ChemNet environment for lit zones (case B); `Q`: energy taken per zone, shared with TitanResidFast.RAYQ. Zones are done outward from the photosphere, each with the photons the last one left. A zone
where the photons of one step could ionize more than `xmin` of its hydrogen stays on the network whatever its pressure, because the table assumes equilibrium.
"""
const RAY = Ref{Any}(nothing)
"""
M3a: the star's light below the Lyman edge in the network zones (`STARLIGHT[] = (; pHm, pH, shield)`, `nothing` = off). Each zone sees
J_ν = W B_ν(T_col): outside the photosphere (τ_R = 2/3, from the outside) T_col = T_eff of the star and W = c(E_r − aT_CMB⁴)/(4σT_eff⁴), so the
dilution and attenuation are those of TITAN's own transfer; inside it T_col = (E_r/a)^¼ and W = 1. Rates, lagged over the step (start-of-step T):
- H⁻ + γ → H + e (`pHm`, John 1988), and optionally (`pH`, `nothing` = off, the default) H(n ≥ 2) + γ → H⁺ + e (the network's H levels without
  n = 1: Lyman-continuum photons are IonRay's, photon-conserving; a dilute Wien tail applied locally would ionize every zone). No extra heat for
  these: in RAD runs the network has no radiative cooling either, all radiative energy exchange is the gray κ_P term (LTE opacities), and ChemNet
  pays the binding energy out of u. That is harmless for H⁻ (x ~ 1e-8, 0.75 eV) but not for H(n ≥ 2): near the photosphere k_H ~ 1 s⁻¹ ionizes
  in ≪ dt what u cannot pay for (13.6 eV per H against ~1 eV of thermal energy), T collapses and ChemNet's Newton fails (zone at 7700 K, 2026-11-05).
  Doing it right needs the H bound-free heating and recombination losses in the network and the H bound-free part removed from κ_P there.
- Lyman-Werner H₂ dissociation (two-step Solomon), k = 1.1e8 F_ν(12.87 eV) s⁻¹ (AAZN97 k27, F_ν = 4πJ_ν), times the H₂ self-shielding
  f_sh(N_H2, T) of Wolcott-Green, Haiman & Bryan (2011, eq. 10, α = 1.1) and the shielding by H I Lyman lines, their eq. 13 (β = 1.6, γ = 0.15,
  N_HI/2.85e23 cm⁻²); fits valid to a factor 2 for 1e13 < N_H2 < 1e21, 1e22 < N_HI < 1e24, 500 < T < 5000 K, static gas. Columns from the
  photosphere out to the zone centre (start-of-step composition). `shield = :static` (that), `:sobolev` (the H₂ column of one Sobolev length
  n v_th/|du/dr| when smaller: the infall's Doppler shifts move the H₂ lines out of the shielding; N_HI stays static, since the H I Lyman lines
  shield through damping wings thousands of km/s wide at N_HI ~ 1e24), `:none` (no shielding), `:off` (no LW). Heat
  4.48 + 0.4 eV per dissociation (binding + fragment kinetic energy; UV-pumping heat in dense gas is left out), moved from E_r to the gas in the
  hydro step (TitanResidFast.GASQ), not in the chemistry: there, without the gas's radiative cooling, the H⁻ → H₂ → LW cycle heated a zone by
  ×2.5 in one step.
`SLDIAG[]`: per-step diagnostics (photosphere zone, T_eff, and per zone k_H⁻, k_H, k_LW, f_sh).
"""
const STARLIGHT = Ref{Any}(nothing); const SLDIAG = Ref{Any}(nothing)
const EVERG = 1.602176634e-12
_bnu(ν, T) = (x = CNn.hP * ν / (CNn.kB * T); x > 700 ? 0.0 : 2CNn.hP * ν^3 / CNn.c^2 / expm1(x))
"H₂ self-shielding, Wolcott-Green, Haiman & Bryan (2011) eq. (10) with α = 1.1; b = (2kT/m_H₂)^½."
fsh_h2(N, T) = (x = N / 5e14; b5 = sqrt(CNn.kB * T / 1.6735575e-24) / 1e5; 0.965 / (1 + x / b5)^1.1 + 0.035 / sqrt(1 + x) * exp(-8.5e-4 * sqrt(1 + x)))
"Shielding of H₂ by H I Lyman-series lines, Wolcott-Green & Haiman (2011) eq. (13)."
fsh_hi(N) = (x = N / 2.85e23; (1 + x)^-1.6 * exp(-0.15x))
"""
Three-level hydrogen (1s, n = 2, continuum) for the network zones under STARLIGHT (`STARLIGHT[].tl3`): the Peebles (1968) C factor as coded in
ChemistryKernels.jl (recombination.jl, `_peebles_k2_from_beta`, RECFAST constants), with the Hubble rate replaced by the Sobolev velocity
gradient of the flow and electron collisions added. Per 2s atom (2s and 2p in statistical equilibrium, n_2p = 3 n_2s):
  C = (R_Lyα + Λ_2s + 4 n_e q21) / (R_Lyα + Λ_2s + 4 n_e q21 + β),
R_Lyα = 3A_21 (1 − e^{−τ_S})/τ_S, τ_S = 3A_21 λ_α³ n_1s / (8π (dv/dr)) (→ 8π (dv/dr)/(λ_α³ n_1s), CK's 1/(K n_1s) with H → dv/dr), dv/dr the
angle mean (|∂u/∂r| + 2|u|/r)/3; Λ_2s = 8.2245809 s⁻¹; β = W α_B(T_c) (C_R T_c)^{3/2} e^{−χ_2/kT_c} (RECFAST's Milne β of a blackbody, diluted by W);
α_B the RECFAST fit (CK `recfast_alpha`); q12 = Λ_ceHI/E_21 from Black (1981, via Cen 1992, as EmissionKernels ceHI: 7.5e-19 e^{−118348/T}
/(1 + (T/1e5)^½) erg cm³ s⁻¹), q21 = q12 e^{E_21/kT}/4. Net rate per volume
  dn_p/dt = (1 − C) n_1s (n_e q12 + P12) − C α_B n_e n_p,   P12 = (R_Lyα + Λ_2s) n̄_α,  n̄_α = W/(e^{E_21/kT_c} − 1)
(pumping by the local field at Lyα). These replace the Boltzmann (LTE) excited levels of the H photo pair and, on the spot, the 1s recombination
(Lyman-continuum photons: IonRay). With T_c = T, W = 1 the pair gives Saha exactly (ground g = 2, the network's Z_H = 2 to 1e-6 below 1e4 K);
with collisions dominant (n_e q21 ≫ R_Lyα, Λ, β) n_2 is Boltzmann at T and the old LTE treatment is recovered.
Energy: ChemNet pays χ_H out of u per net ionization (E0). The chemistry also gets, per ionization through n = 2, the photon energy ⟨hν⟩ (Balmer
continuum mean, σ ∝ ν⁻³ in W B(T_c)) + E_21 P12/(n_e q12 + P12) (the pumped share of the excitation), and per net recombination loses
χ_H + 0.75 kT to radiation; that energy leaves (or enters) the gray field through RAYQ. Inside the chemistry so that the collisional-excitation
cost (10.2 eV per ionization, ≫ kT) limits the rate implicitly in T. The recombine/re-ionize cycle's heat (1 − C) α_B n_e n_p (⟨hν⟩ − χ_2 −
0.75 kT) goes to the gas in the hydro step (GASQ), where it meets the gas's radiative cooling. Left out: H–H collisional excitation (Drawin),
n ≥ 3 except through α_B, and the Lyα cooling of excitations that do not ionize (in the gray κ_P's lines).
"""
const TL3DBG = Ref(parse(Int, get(ENV, "TL3DBG", "0")))
"Print the chemistry/table step of this zone every step (ENV ZONEDBG; diagnostics)."
const ZONEDBG = Ref(parse(Int, get(ENV, "ZONEDBG", "0")))
const LYA_LAM = 1.215668e-5; const LAM2S = 8.2245809; const A21LYA = 6.2648e8
const E21 = 10.19881 * EVERG; const CHI2 = CNn.χH - E21; const CR_REC = 1.799920e10      # 2π m_e k/h² [cm⁻² K⁻¹]
alphaB_recfast(T) = (t = T / 1e4; 4.309e-13 * t^-0.6166 / (1 + 0.6703 * t^0.53))          # cm³ s⁻¹ (= CK recfast_alpha × 1e6)
beta2s_bb(Tc) = alphaB_recfast(Tc) * (CR_REC * Tc)^1.5 * exp(-CHI2 / (CNn.kB * Tc))         # s⁻¹ per 2s atom, J = B(T_c)
q12_H(T) = 7.5e-19 * exp(-118348.0 / T) / (1 + sqrt(T / 1e5)) / E21
q21_H(T) = 7.5e-19 / (1 + sqrt(T / 1e5)) / E21 / 4
tl3_C(T, p) = (D = p.RL + LAM2S + 4p.ne * q21_H(T); D / (D + p.β))
"Mean photon energy of n = 2 photoionization (σ ∝ ν⁻³) in J ∝ B(T_c)."
function balmer_mean(Tc)
    ν0 = CHI2 / CNn.hP; s0 = 0.0; s1 = 0.0
    νs = [ν0 * exp(x) for x in range(0.0, log(1 + 40CNn.kB * Tc / CHI2), length = 101)]
    for i in 1:100
        ν = 0.5(νs[i] + νs[i+1]); w = _bnu(ν, Tc) / ν^4 * (νs[i+1] - νs[i]); s0 += w; s1 += w * CNn.hP * ν
    end
    s0 > 0 ? s1 / s0 : CHI2
end
"The three reactions of the three-level H atom (zero unless the zone's env carries `tl3`)."
function tl3_reactions()
    on(env) = hasproperty(env, :tl3) && env.tl3 !== nothing
    [CNn.Reaction("H + e → H(n=2) + e, then n = 2 photoionization (three-level, (1−C) q12)", [:H, :e], [:Hp, :e, :e],
                  (T, env) -> on(env) ? (1 - tl3_C(T, env.tl3)) * q12_H(T) : zero(T)),
     CNn.Reaction("H + γ(Lyα) → H(n=2), then n = 2 photoionization (three-level, (1−C) P12)", [:H], [:Hp, :e],
                  (T, env) -> on(env) ? (1 - tl3_C(T, env.tl3)) * env.tl3.P12 + zero(T) : zero(T)),
     CNn.Reaction("H+ + e → H(n≥2) → H(1s) (three-level, C α_B)", [:Hp, :e], [:H],
                  (T, env) -> on(env) ? tl3_C(T, env.tl3) * alphaB_recfast(T) : zero(T))]
end
"""
The network with the three-level H atom: `tl3_reactions` added, and mutual neutralization H⁻ + H⁺ → 2H (GA08 k5) split by the same C: it
leaves one H in n ≥ 2 (mostly n = 3), which the star light ionizes again with probability 1 − C, so in lit zones only C of it recombines
(the rest is H⁻ + H⁺ → H + H⁺ + e). Without the split, at x_e ~ 0.4 it recombined ~1 % of the H per step and put 12.9 eV per event into
the gas (the network has no radiative losses in RAD runs). Its detailed-balance reverse is unchanged.
"""
function tl3_network(net)
    on(env) = hasproperty(env, :tl3) && env.tl3 !== nothing
    rx = CNn.Reaction[]
    for r in net.reactions
        if startswith(r.name, "H- + H+ → 2H") && !occursin("reverse", r.name)
            k0 = r.k
            push!(rx, CNn.Reaction(r.name * " × C (three-level)", r.R, r.P, (T, env) -> on(env) ? tl3_C(T, env.tl3) * k0(T, env) : k0(T, env)))
            push!(rx, CNn.Reaction("H- + H+ → H(n≥2) + H+, re-ionized: H + H+ + e ((1−C) × GA08 k5)", [:Hm, :Hp], [:H, :Hp, :e],
                                   (T, env) -> on(env) ? (1 - tl3_C(T, env.tl3)) * k0(T, env) : zero(T)))
        else
            push!(rx, r)
        end
    end
    CNn.Network(net.species, net.free, vcat(rx, tl3_reactions()))
end
"Energy terms of the three-level atom per volume: (chemistry heat, cycle heat for GASQ) at (ρ, T, full composition d)."
function tl3_heat(ρ, T, d, p)
    nH = ρ / CNn.μH(); C = tl3_C(T, p); xH = max(d[:H], zero(T)); xe = max(d[:e], zero(T)); xp = max(d[:Hp], zero(T))
    Q12 = nH * xe * q12_H(T); I = (1 - C) * xH * (Q12 + p.P12); rec = alphaB_recfast(T) * nH * xe * xp
    Ein = p.hν + E21 * p.P12 / (Q12 + p.P12 + 1e-300)
    chemq = nH * (I * Ein - C * rec * (CNn.χH + 0.75CNn.kB * T))
    cycq = nH * (1 - C) * rec * (p.hν - CHI2 - 0.75CNn.kB * T)
    LYABB[] && (chemq += lya_bb(ρ, T, d, p))
    chemq, cycq
end
"Bound–bound n = 1 ↔ 2 heat of the three-level atom per volume [erg cm⁻³ s⁻¹] (see LYABB)."
function lya_bb(ρ, T, d, p)
    nH = ρ / CNn.μH(); xH = max(d[:H], zero(T)); xe = max(d[:e], zero(T)); xp = max(d[:Hp], zero(T))
    # branching of n = 2 as in the reactions and in chemq: C and D with the lagged p.ne, so the de-excited share 4 p.ne q21/(D + β) ≤ 1 of
    # every n = 2 atom is the one the network's C assumes (with the current n_e in the collisional term and p.ne in D, a zone ionizing
    # within the step de-excited 41× its n = 2 atoms: energy from nothing, 3654 K → 5e5 K, l2 step 99005); collisional excitation with
    # the current n_e (as the network's H + e → H(n=2) rate), its radiative share (R + Λ)/(D + β) ≤ 1
    rec = alphaB_recfast(T) * nH * xe * xp; D = p.RL + LAM2S + 4p.ne * q21_H(T); C = D / (D + p.β)
    nH * E21 * C * ((xH * p.P12 + rec) * 4p.ne * q21_H(T) - xH * nH * xe * q12_H(T) * (p.RL + LAM2S)) / D
end
"""
LYABB=1: the bound–bound energy exchange of the three-level atom with the gas, missing from tl3_heat's ionization bookkeeping. The atom's
exact radiative balance with the gas is E21 [n1 P12 − n2s (R_Lyα + Λ_2s)] + β n2s ⟨hν⟩ − α_B n_e n_p (χ_2 + 0.75 kT) (n2s in steady state);
tl3_heat's chemq + cycq differ from it by
  E21 C [(n1 P12 + α_B n_e n_p) 4 n_e q21 − n1 n_e q12 (R_Lyα + Λ_2s)] / D,   D = R_Lyα + Λ_2s + 4 n_e q21:
n = 2 atoms made by Lyα/two-photon pumping (excitation temperature of the diluted starlight T_c) or by recombination that are collisionally
de-excited heat the gas; collisional excitations that decay radiatively (Lyα escape, two-photon continuum) cool it. With W = 1, T_c = T the
sum vanishes (detailed balance; without the term, LTE gas lost E21 β n2s 4 n_e q21/(D+β)). The energy goes to/from the gray field via RAYQ
with the rest of chemq, implicit in T in the chemistry. Left out: H–H (Drawin) collisions.
"""
const LYABB = Ref(get(ENV, "LYABB", "0") == "1")
"Per-step star-light setup: photosphere, T_eff, and the LW columns; returns a function k -> (kHm, kH, kLW, fsh) for a network zone k."
function starlight_setup(m::NetGas, a, b, tauz_)
    sl = STARLIGHT[]; dn, tn, rn = TM.v(:dn), TM.v(:tn), TM.v(:rn); er, fr, un = TM.v(:ern), TM.v(:frn), TM.v(:un)
    kph = a; for k in b:-1:a; tauz_[k] >= 2 / 3 && (kph = k; break); end
    L = 4π * rn[kph+1]^2 * fr[kph+1]; Te = L > 0 ? (L / (4π * rn[kph+1]^2 * 5.670374419e-5))^0.25 : 0.0
    arad = 7.565723e-15; Ecmb = arad * sl.tcmb^4
    nz = length(dn); NH2 = zeros(nz); NHI = zeros(nz)
    s2 = 0.0; s1 = 0.0
    for k in kph+1:b
        d = CNn.full(m.x[k], m.net); nH = dn[k] / CNn.μH(); dr = rn[k+1] - rn[k]
        nh2 = nH * max(d[:H2], 0.0); nhi = nH * max(d[:H], 0.0)
        NH2[k] = s2 + 0.5nh2 * dr; NHI[k] = s1 + 0.5nhi * dr
        if sl.shield === :sobolev
            dudr = abs(un[k+1] - un[k]) / dr
            if dudr > 0
                vth = sqrt(2CNn.kB * tn[k] / 1.6735575e-24)
                NH2[k] = min(NH2[k], nh2 * vth / sqrt(2) / dudr)     # (N_HI stays static: its Lyman lines shield through Lorentz wings ≫ the infall's Δv)
            end
        end
        s2 += nh2 * dr; s1 += nhi * dr
    end
    νLW = 12.87 * EVERG / CNn.hP
    diag = (; kph, Te, L, kHm = zeros(nz), kH = zeros(nz), kLW = zeros(nz), fsh = ones(nz), tl3 = Vector{Any}(nothing, nz))
    tl3on = hasproperty(sl, :tl3) && sl.tl3
    SLDIAG[] = diag
    # T_r at optical depth τ (log–log interpolation in the zones a..b; τ beyond the innermost zone: its T_r), for SLEB
    function tr_at_tau(τ)
        τ >= tauz_[a] && return (max(er[a], Ecmb) / arad)^0.25
        for j in a:b-1
            if tauz_[j] >= τ > tauz_[j+1]
                t1, t2 = log(max(tauz_[j], 1e-30)), log(max(tauz_[j+1], 1e-30)); w = t1 > t2 ? (log(τ) - t2) / (t1 - t2) : 0.0
                l1, l2 = log(max(er[j], Ecmb, 1e-300)), log(max(er[j+1], Ecmb, 1e-300))
                return (exp(l2 + w * (l1 - l2)) / arad)^0.25
            end
        end
        (max(er[b], Ecmb) / arad)^0.25
    end
    function rates(k)
        if SLBLEND[] && Te > 0
            # smooth in τ (no jump at the moving τ = 2/3 zone): colour temperature from T_eff (τ ≤ 0.3) to the local T_r (τ ≥ 3), log blend,
            # dilution W from the local energy density (W = 1 at T_c = T_r). The sharp switch changed x_ss of the zone just outside τ = 2/3
            # by ×2 at equal energy density (diluted 6500 K vs a 4750 K blackbody).
            Tr = (max(er[k], Ecmb) / arad)^0.25
            Tsrc = SLEB[] ? tr_at_tau(tauz_[k] + 2 / 3) : Te
            wq = clamp((log10(max(tauz_[k], 1e-30)) + log10(1 / 0.3)) / 1.0, 0.0, 1.0); wq = wq * wq * (3 - 2wq)
            Tc = Tsrc^(1 - wq) * Tr^wq
            W = wq >= 1 ? 1.0 : min(max(CNn.c * (er[k] - Ecmb), 0.0) / (4 * 5.670374419e-5 * Tc^4), 1.0)
        elseif k > kph && Te > 0
            Tc = Te; W = max(CNn.c * (er[k] - Ecmb), 0.0) / (4 * 5.670374419e-5 * Te^4)
        else
            Tc = (max(er[k], Ecmb) / arad)^0.25; W = 1.0
        end
        (W > 0 && Tc > 0) || return (0.0, 0.0, 0.0, 1.0)
        env = (Trad = Tc, J = ν -> W * _bnu(ν, Tc))
        kHm = sl.pHm === nothing ? 0.0 : CNn.photorate(sl.pHm, tn[k], env); kH = sl.pH === nothing ? 0.0 : CNn.photorate(sl.pH, tn[k], env)
        f = 1.0; kLW = 0.0
        if k > kph && sl.shield !== :off
            f = sl.shield === :none ? 1.0 : fsh_h2(NH2[k], tn[k]) * fsh_hi(NHI[k])
            kLW = 1.1e8 * 4π * W * _bnu(νLW, Tc) * f
        end
        diag.kHm[k] = kHm; diag.kH[k] = kH; diag.kLW[k] = kLW; diag.fsh[k] = f
        if tl3on                                        # three-level H parameters (start-of-step state), see tl3_reactions
            d = CNn.full(m.x[k], m.net); nH = dn[k] / CNn.μH(); n1 = nH * max(d[:H], 1e-30)
            rc = 0.5(rn[k] + rn[k+1]); ub = 0.5(un[k] + un[k+1])
            g = (abs(un[k+1] - un[k]) / (rn[k+1] - rn[k]) + 2abs(ub) / rc) / 3
            τS = 3A21LYA * LYA_LAM^3 * n1 / (8π * max(g, 1e-30))
            RL = 3A21LYA * (τS < 1e-6 ? 1.0 : -expm1(-τS) / τS)
            nb = W / expm1(E21 / (CNn.kB * Tc))
            diag.tl3[k] = (; β = W * beta2s_bb(Tc), RL, ne = nH * max(d[:e], 0.0), P12 = (RL + LAM2S) * nb, hν = balmer_mean(Tc), Tc)
            ov_ = LYAOVR[]; (k <= length(ov_) && ov_[k] !== nothing) && (diag.tl3[k] = merge(diag.tl3[k], ov_[k]))   # LYART: transferred Lyα
        end
        (kHm, kH, kLW, f)
    end
    tl3on && LYART[] && lya_transfer!(m, a, b, tauz_, kph, rates, diag)
    rates
end
include(joinpath(@__DIR__, "LyaTransfer.jl"))
"dx/dt of the star-light reactions (rates lagged) for ChemNet.step's `extra`. No heat here: the LW energy goes to TitanResidFast.GASQ (see `lwheat`)."
function starlight_terms(m::NetGas, ρ, kHm, kH, kLW)
    fr_ = m.net.free; iHp = findfirst(==(:Hp), fr_); iHm = findfirst(==(:Hm), fr_); iH2 = findfirst(==(:H2), fr_)
    nH = ρ / CNn.μH()
    extra = x -> begin
        f = zeros(eltype(x), length(x)); d = CNn.full(x, m.net)
        iHm === nothing || (f[iHm] -= kHm * x[iHm])
        f[iHp] += kH * max(d[:H], zero(d[:H]))
        iH2 === nothing || (f[iH2] -= kLW * x[iH2])
        f
    end
    extra
end
"CHEMTHREADS=1: the zones of `chem!(::NetGas)` in parallel (Threads.@threads) when no ray is set (the ray's photons pass the zones in order). Default off."
const CHEMTHREADS = Ref(get(ENV, "CHEMTHREADS", "0") == "1")
const CHEMLOCK = ReentrantLock()
"""
CHEMSCHED=zone (default): the chemistry tasks take zones one at a time from a shared counter. `chunk`: `Threads.@threads`,
which gives every task one contiguous block of zones; the expensive zones (sub-cycled, fallback) are neighbours, so a few
tasks then carry most of the work. The zones are independent: the results are the same.
"""
const CHEMSCHED = Ref(Symbol(get(ENV, "CHEMSCHED", "zone")))
"""
CHEMSUB=fT,fx (default off): a three-level (tl3) network zone whose single backward-Euler chemistry step changes T by more than the factor
fT, or n_e by more than the factor fx (with x_e > 1e-4 at either end), is integrated again in halved sub-steps, each held to the same test
(to depth CHEMSUBMAX), with n_e in the Peebles C refreshed from the current composition between sub-steps.
Why: C (and with it the share 1 − C of collisionally excited n = 2 atoms that the star light ionizes) is frozen at the start-of-step n_e.
A dense zone just out of the table (n 9e13, T 7000 K, x_e 3e-5, β 7e4 /s: C 2e-3) then has an electron avalanche with an e-folding time
of 6 s inside a step of 1e4 s; backward Euler on a growing mode with h λ ≫ 1 has spurious roots, and Newton ended on one at 2.6e6 K, fully
ionized, the recombine/re-ionize loop feeding ⟨hν⟩ per cycle: 5.5e14 erg/g from nothing (full3 step 34704, zone 157, the photosphere zone
in the accretion shock; the radiation spike that followed stopped the run). With the sub-steps the same zones
(n 5e13–2.4e14, T 5300–5500 K at the start of the step) end at x_e 2e-4–1.6e-3 with T 4–5 % lower, and no step heats (test t_cs1, steps
34700–34900); where a single step still converges on the right branch it overshoots x_e by 2–5×.
"""
const CHEMSUB = Ref{Union{Nothing,NTuple{2,Float64}}}(let s = get(ENV, "CHEMSUB", "0"); s == "0" ? nothing : (v = parse.(Float64, split(s, ",")); (v[1], length(v) > 1 ? v[2] : 3.0)) end)
const CHEMSUBMAX = Ref(parse(Int, get(ENV, "CHEMSUBMAX", "16"))); const NCHEMSUB = Ref(0)
"""
CHEMREPLAY=n (diagnostic, default 0): the first n zone steps whose CHEMSUB sub-stepping reaches CHEMSUBMAX are solved again in place
under solver variants (reference rule; T criterion only; depth caps; n_e iterated to self-consistency in the three-level C and P),
printing per variant the ChemNet.step calls, the wall time and the end state (T, x_e). The zone's result is the production one.
"""
const CHEMREPLAY = Ref(parse(Int, get(ENV, "CHEMREPLAY", "0"))); const NREPLAY = Ref(0)
function chemreplay(k, m, ρ, u0, T0, x0, dt, phk, exs, c0, p3; xref, Tref)
    nHk = ρ / CNn.μH(); xel(x) = max(CNn.full(x, m.net)[:e], 0.0)
    ncall = Ref(0)
    function st(xa, Ta, ua, hh, pa; npic = 0)
        pc = pa; xb = xa; Tb = Ta
        for it in 0:npic
            ca = (ρ_, T, d) -> (c0 === nothing ? zero(T) : c0(ρ_, T, d)) - tl3_heat(ρ_, T, d, pc)[1]
            ncall[] += 1
            xb, Tb = CNn.step(xa, ρ, ua, Ta, hh, m.net; env = merge(m.env, (; tabkey = m.env, tl3 = pc)), photos = phk, cooling = ca, extra = exs)
            pn = merge(pc, (; ne = nHk * xel(xb)))
            it < npic && abs(pn.ne / pc.ne - 1) < 1e-3 && (pc = pn; break)
            pc = pn
        end
        xb, Tb, pc
    end
    function run(fT, fx, dmaxc; npic = 0)
        bad(xa, Ta, xb, Tb) = (ea = xel(xa); eb = xel(xb); !(isfinite(Tb) && all(isfinite, xb)) || abs(log(Tb / Ta)) > log(fT) ||
                               (max(ea, eb) > 1e-4 && max(ea, eb) > fx * min(ea, eb)))
        function sub(xa, Ta, ua, hh, depth, pa)
            xb, Tb, pb = st(xa, Ta, ua, hh, pa; npic)
            (depth >= dmaxc || !bad(xa, Ta, xb, Tb)) && return (xb, Tb, pb)
            xc, Tc_, pc = sub(xa, Ta, ua, hh / 2, depth + 1, pa)
            sub(xc, Tc_, Float64(CNn.eos(ρ, Tc_, xc, m.net)[1]), hh / 2, depth + 1, pc)
        end
        ncall[] = 0; t0 = time(); j0 = CNn.NJAC[]; a0 = CNn.NADV[]
        r = try
            xb, Tb, pb = st(x0, T0, u0, dt, p3; npic)
            if dmaxc > 0 && bad(x0, T0, xb, Tb)
                xm, Tm, pm = sub(x0, T0, u0, dt / 2, 1, p3)
                xb, Tb, pb = sub(xm, Tm, Float64(CNn.eos(ρ, Tm, xm, m.net)[1]), dt / 2, 1, pm)
            end
            (Tb, xel(xb))
        catch e
            (NaN, NaN)
        end
        (ncall[], time() - t0, r..., CNn.NJAC[] - j0, CNn.NADV[] - a0)
    end
    @printf("chemreplay k %d step %d n %.2e dt %.2e T0 %.1f xe0 %.2e | production: T %.1f xe %.3e\n", k, TF.getscalar(:jstep), nHk, dt, T0, xel(x0), Tref, xel(xref))
    for (lab, fT, fx, dm, np) in (("ref 1.3,3 d16", 1.3, 3.0, 16, 0), ("T-only 1.3 d16", 1.3, Inf, 16, 0), ("1.3,3 d8", 1.3, 3.0, 8, 0),
                                  ("1.3,3 d4", 1.3, 3.0, 4, 0), ("single step", 1.3, 3.0, 0, 0), ("single + ne-picard", 1.3, 3.0, 0, 8),
                                  ("ne-picard T-only d16", 1.3, Inf, 16, 4), ("ne-picard 1.3,3 d16", 1.3, 3.0, 16, 4), ("1.3,10 d16", 1.3, 10.0, 16, 0))
        nc, w, T, xe, nj, na = run(fT, fx, dm; npic = np)
        @printf("    %-22s calls %6d  wall %8.3f s  T %8.1f  xe %.3e  Newton its %7d  attempts %6d  (%.2f ms/its)\n", lab, nc, w, T, xe, nj, na, 1e3w / max(nj, 1))
    end
    # cost of the residual's pieces at the start state (Float64 and with 9-wide dual numbers, as in the ForwardDiff Jacobian)
    d0 = CNn.full(x0, m.net); nr = 200
    ca = (ρ_, T, d) -> (c0 === nothing ? zero(T) : c0(ρ_, T, d))
    Td = ForwardDiff.Dual(T0, ntuple(_ -> 1.0, 9)...)
    tim(f) = (f(); t = time(); for _ in 1:nr; f(); end; 1e6 * (time() - t) / nr)
    @printf("    pieces [us/eval, Float64 / dual]: cooling %.1f / %.1f  tl3_heat %.1f / %.1f  eos %.1f / %.1f\n",
            tim(() -> ca(ρ, T0, d0)), tim(() -> ca(ρ, Td, d0)), tim(() -> tl3_heat(ρ, T0, d0, p3)), tim(() -> tl3_heat(ρ, Td, d0, p3)),
            tim(() -> CNn.eos(ρ, T0, x0, m.net)), tim(() -> CNn.eos(ρ, Td, x0, m.net)))
end
"""
`CHEMTMIN=T` (K, default 0 = off): a chemistry step (single or sub-stepped) that returns T below this floor is treated as failed
and the zone takes the table fallback for that step (CHEMFALLBACK). Added 12 Nov after f128 stopped at 46.6 M☉: zone 69 (795 R☉,
n 8.8e11) went 1516 K -> 0.3 K in the sub-stepped chemistry (single step 57 K) and the next hydro solve was NaN. Below the CMB
(55 K at z = 19) such a T is a spurious root of the implicit step, not physics.
"""
const CHEMTMIN = Ref(parse(Float64, get(ENV, "CHEMTMIN", "0"))); const NCHEMTMIN = Ref(0)
"""
CHEMZPROF=1 (diagnostic, default off): wall time of the chemistry per zone; every 500 steps the twelve zones that took most of it.
"""
const CHEMZPROF = get(ENV, "CHEMZPROF", "0") == "1"
const ZPT = Ref(Float64[]); const ZPN = Ref(0)
"The chemistry–temperature step of the H network (`chem!` adds the deuterium trace network with HDCHEM=1)."
function chem_core!(m::NetGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    dn, tn, egn, pgn = TM.v(:dn), TM.v(:tn), TM.v(:egn), TM.v(:pgn)
    dTmax = 0.0
    h = HOT[]
    ray = RAY[]; kph = typemax(Int); Nin = nothing
    if ray !== nothing
        fill!(ray.Q, 0.0)
        st = ray.start(a, b); st === nothing || ((kph, Nin) = st)
    end
    rn = TM.v(:rn)
    tauz = zeros(length(dn))
    if PSWTAU[] > 0 || STARLIGHT[] !== nothing
        chi_ = TM.v(:chifn); tt = 0.0
        for k in b:-1:a; tt += chi_[k] * dn[k] * (rn[k+1] - rn[k]); tauz[k] = tt; end
    end
    slrates = STARLIGHT[] === nothing ? nothing : starlight_setup(m, a, b, tauz)
    # the star-light terms of zone k for ChemNet.step: (extra, cooling) composed with the ray's / the zone's own
    function slterms(k, extra0, cool0)
        slrates === nothing && return extra0, cool0
        kHm, kH, kLW, _ = slrates(k)
        all(isfinite, (kHm, kH, kLW)) || error("starlight: non-finite rates in zone $k: kHm $kHm kH $kH kLW $kLW, T $(tn[k]), ρ $(dn[k]), E_r $(TM.v(:ern)[k])")
        (kHm == 0 && kH == 0 && kLW == 0) && return extra0, cool0
        ex = starlight_terms(m, dn[k], kHm, kH, kLW)
        ex2 = extra0 === nothing ? ex : (x -> extra0(x) .+ ex(x))
        ex2, cool0
    end
    # LW energy absorbed in zone k over the step, from the end-of-step H₂ (steady under LW): 4.48 eV binding + 0.4 eV fragment kinetic energy
    # per dissociation (UV pumping heat left out); the chemistry is energy-neutral (it takes the binding energy from u and gives it back on re-formation)
    lwheat!(k, x1) = (slrates === nothing || isempty(TX.GASQ[])) ? nothing :
        (TX.GASQ[][k] = SLDIAG[].kLW[k] * max(CNn.full(x1, m.net)[:H2], 0.0) * dn[k] / CNn.μH() * (4.48 + 0.4) * EVERG; nothing)
    slrates === nothing || isempty(TX.GASQ[]) || fill!(TX.GASQ[], 0.0)
    ray === nothing && !isempty(TX.RAYQ[]) && fill!(TX.RAYQ[], 0.0)            # (RAYQ also carries the three-level H's photon energy)
    isempty(TX.CHEMB[]) || fill!(TX.CHEMBON[], false)
    ph3 = filter(p -> p.R !== :H, m.photos)                                 # three-level zones: no Boltzmann-level H photo pair
    tl3sw = slrates !== nothing && hasproperty(STARLIGHT[], :tl3) && STARLIGHT[].tl3; iHpf = findfirst(==(:Hp), m.net.free)
    if NLOPAC[]
        length(NLOPC[]) == length(dn) || (NLOPC[] = Any[nothing for _ in dn])
        NLODAMP < 1 && (NLOPREV[] = copy(NLOPC[]))
        fill!(NLOPC[], nothing)
    end
    dTz = zeros(length(tn))
    function zone!(k)
        if h !== nothing && tn[k] > h.T1             # hot branch: no chemistry, no re-solve of T with the table (driver: burn and mix)
            m.x[k] = x_from_table(m, dn[k], min(tn[k], 1e6)); m.eqzone[k] = true
            return nothing
        end
        cool = zone_cooling(m, k)
        if Nin !== nothing && k > kph
            V = 4π / 3 * (rn[k+1]^3 - rn[k]^3)
            if sum(Nin) * dt / (dn[k] / CNn.μH() * V) > ray.xmin
                src = ray.source(dn[k], rn[k], rn[k+1], Nin)
                cool2 = (ρ, T, d) -> (cool === nothing ? zero(T) : cool(ρ, T, d)) - src.heat(d)
                m.eqzone[k] = false
                ex3, cool3 = slterms(k, src.extra, cool2)
                x1, T1 = try
                    CNn.step(m.x[k], dn[k], egn[k], tn[k], dt, m.net; env = ray.env, photos = m.photos, cooling = cool3, extra = ex3, accepted = src.accepted)
                catch e
                    println("ray zone failed: k ", k, " x ", m.x[k], " free ", m.net.free, " ρ ", dn[k], " u ", egn[k], " T ", tn[k], " dt ", dt, " Nin ", sum(Nin), " r ", rn[k], " ", rn[k+1])
                    rethrow()
                end
                lwheat!(k, x1)
                out = src.out(x1, dt)
                ray.Q[k] = sum((Nin .- out) .* ray.Ec) / V
                m.x[k] = x1; Nin = out
                dTz[k] = abs(T1 / tn[k] - 1); tn[k] = T1
                return nothing
            end
            Nin = nothing                                # the rest of the light cannot change any zone further out
        end
        # PSWTMIN[] > 0: also require T above it (5 % hysteresis). Below ~5000 K ionization holds negligible energy, so the network and the table agree there
        # to ~1 % in T at fixed e; at 7000-9500 K the network's H ionization lags the table's Saha by ×5, and switching there drops T by up to 16 %.
        ptab = (pgn[k] >= m.pswitch && tn[k] >= PSWTMIN[]) || (m.eqzone[k] && pgn[k] >= PSWHYS[] * m.pswitch && tn[k] >= 0.95 * PSWTMIN[])
        ag = nothing                                  # PSWAGREE: table only where it agrees with the network (see PSWAGREE)
        if !ptab && PSWAGREE[] > 0 && tl3sw && tn[k] >= PSWTEX[]
            slrates(k); p = SLDIAG[].tl3[k]
            if tauz[k] >= PSWTAUDEEP[]
                ag = true                     # thermalised (J = B): the table's implicit LTE coupling; the three-level n = 2 exchange is explicit
            elseif p !== nothing
                # the two equilibrium models at the zone's u: table (Saha) and three-level steady state, both with their ionization heat
                # capacity; u from the network EOS at the table composition (offset to the table EOS ≤ 0.3 % in T)
                ε = PSWAGREE[]; _, Ttab = EQ.step(dn[k], egn[k], tn[k], dt)
                u0 = CNn.eos(dn[k], Ttab, x_from_table(m, dn[k], Ttab), m.net)[1]
                Tq = T_ss_at_u(m, dn[k], u0, p, Ttab, iHpf)
                if m.eqzone[k]
                    ag = isfinite(Tq) && abs(Tq / Ttab - 1) < 2ε
                    ag || lock(() -> (NSWOUT[] += 1), CHEMLOCK)
                else
                    ag = abs(Ttab / tn[k] - 1) < ε && isfinite(Tq) && abs(Tq / Ttab - 1) < ε
                    ag && lock(() -> (NSWIN[] += 1), CHEMLOCK)
                end
            end
        end
        tab = ptab || (ag === nothing ? ((PSWSTAY[] && m.eqzone[k] && PSWTEX[] > 0 && tn[k] >= PSWTEX[]) || (PSWTAU[] > 0 && tn[k] >= PSWTEX[] && tauz[k] >= (m.eqzone[k] ? 0.5 : 1.0) * PSWTAU[])) : ag)
        if tab && !m.eqzone[k] && PSWRAD[] !== nothing && pgn[k] < PSWRAD[][2] * m.pswitch       # entry only with T_gas = T_rad (see PSWRAD)
            Tr = (max(TM.v(:ern)[k], 0.0) / 7.565723e-15)^0.25
            if abs(Tr / tn[k] - 1) > PSWRAD[][1]
                tab = false; nsw = lock(() -> (NSWRAD[] += 1), CHEMLOCK)
                (nsw <= 40 || nsw % 500 == 0 || TL3DBG[] > 0) && @printf("pswrad k %d step %d (%d): stays on the network, T %.1f T_rad %.1f tau %.2f p/pswitch %.3g\n", k, TF.getscalar(:jstep), nsw, tn[k], Tr, tauz[k], pgn[k] / m.pswitch)
            end
        end
        if tab
            coolt = cool === nothing ? nothing : ((ρ, T, y) -> cool(ρ, T, nothing))
            _, T1 = EQ.step(dn[k], egn[k], tn[k], dt; cooling = coolt)
            m.x[k] = x_from_table(m, dn[k], T1); m.eqzone[k] = true
        else
            eqprev = m.eqzone[k]
            # (back to the ideal regime: the network continues from the table composition). With the three-level atom, from the table composition
            # at the zone's current T, not m.x: m.x of a table zone is still advected, and took in the ionized network gas upstream (x_e 0.3–0.4)
            # while the table's e had radiated that energy away; the network then started 3e12 erg/g below u(10 K, x) (zones at 1800–1900 R_sun).
            eqprev && STARLIGHT[] !== nothing && hasproperty(STARLIGHT[], :tl3) && STARLIGHT[].tl3 && (m.x[k] = x_from_table(m, dn[k], tn[k]))
            m.eqzone[k] && (m.eqzone[k] = false)
            exs, cools = slterms(k, nothing, cool)
            p3 = (slrates === nothing || m.eqzone[k]) ? nothing : SLDIAG[].tl3[k]
            envk, phk = m.env, m.photos; fb = false
            p3e = p3                                 # the three-level parameters at the end of the chemistry (n_e refreshed if CHEMSUB sub-stepped)
            if p3 !== nothing
                envk = merge(m.env, (; tabkey = m.env, tl3 = p3)); phk = ph3
                c0 = cools
                cools = (ρ, T, d) -> (c0 === nothing ? zero(T) : c0(ρ, T, d)) - tl3_heat(ρ, T, d, p3)[1]
            end
            if !isempty(TX.CHEMB[])                  # BALSPLIT: the hydro step already gave the zone c̄ (TitanResidFast.CHEMB)
                cb = eqprev ? 0.0 : TX.CHEMB[][k]; TX.CHEMBON[][k] = !eqprev
                if cb != 0
                    c1 = cools; cools = (ρ, T, d) -> (c1 === nothing ? zero(T) : c1(ρ, T, d)) + ρ * cb
                end
            end
            x1, T1 = try
                # u below the zone's own composition at 10 K (x and e advected by different schemes; with x_e ~ 0.1-0.5 from the three-level
                # atom the mismatch reaches ~1 eV per H): no solution for this x, take the fallback below
                p3 !== nothing && egn[k] < CNn.eos(dn[k], 10.0, m.x[k], m.net)[1] && error("tl3: u below u(10 K, x) in zone $k")
                xs1, Ts1 = CNn.step(m.x[k], dn[k], egn[k], tn[k], dt, m.net; env = envk, photos = phk, cooling = cools, extra = exs)
                if p3 !== nothing && CHEMSUB[] !== nothing && isempty(TX.CHEMB[])
                    fT, fx = CHEMSUB[]; nHk = dn[k] / CNn.μH(); xel(x) = max(CNn.full(x, m.net)[:e], 0.0)
                    bad(xa, Ta, xb, Tb) = (ea = xel(xa); eb = xel(xb); !(isfinite(Tb) && all(isfinite, xb)) || abs(log(Tb / Ta)) > log(fT) ||
                                           (max(ea, eb) > 1e-4 && max(ea, eb) > fx * min(ea, eb)))
                    if bad(m.x[k], tn[k], xs1, Ts1)
                        dmax = Ref(0)
                        function sub(xa, Ta, ua, hh, depth, pa)
                            dmax[] = max(dmax[], depth)
                            ca = (ρ, T, d) -> (c0 === nothing ? zero(T) : c0(ρ, T, d)) - tl3_heat(ρ, T, d, pa)[1]
                            xb, Tb = CNn.step(xa, dn[k], ua, Ta, hh, m.net; env = merge(m.env, (; tabkey = m.env, tl3 = pa)), photos = phk, cooling = ca, extra = exs)
                            (depth >= CHEMSUBMAX[] || !bad(xa, Ta, xb, Tb)) && return (xb, Tb, merge(pa, (; ne = nHk * xel(xb))))
                            xc, Tc_, pc = sub(xa, Ta, ua, hh / 2, depth + 1, pa)
                            sub(xc, Tc_, Float64(CNn.eos(dn[k], Tc_, xc, m.net)[1]), hh / 2, depth + 1, pc)
                        end
                        xs0, Ts0 = xs1, Ts1
                        xm, Tm, pm = sub(m.x[k], tn[k], Float64(egn[k]), dt / 2, 1, p3)
                        xs1, Ts1, p3e = sub(xm, Tm, Float64(CNn.eos(dn[k], Tm, xm, m.net)[1]), dt / 2, 1, pm)
                        lock(() -> (NCHEMSUB[] += 1), CHEMLOCK)
                        dmax[] >= CHEMSUBMAX[] && CHEMREPLAY[] > 0 && lock(() -> (NREPLAY[] += 1) <= CHEMREPLAY[], CHEMLOCK) &&
                            lock(() -> chemreplay(k, m, dn[k], Float64(egn[k]), tn[k], m.x[k], dt, phk, exs, c0, p3; xref = xs1, Tref = Ts1), CHEMLOCK)
                        (dmax[] >= 4 || TL3DBG[] > 0) && @printf("chemsub k %d step %d r %.3e R_sun n %.2e dt %.2e depth %d: T %.1f -> %.1f (one step: %.1f), x_e %.2e -> %.2e (one step: %.2e)\n", k, TF.getscalar(:jstep), rn[k] / 6.957e10, nHk, dt, dmax[], tn[k], Ts1, Ts0, xel(m.x[k]), xel(xs1), xel(xs0))
                    end
                end
                CHEMTMIN[] > 0 && Ts1 < CHEMTMIN[] && (lock(() -> (NCHEMTMIN[] += 1), CHEMLOCK); error("chem: T $(Ts1) below CHEMTMIN in zone $k (from $(tn[k]))"))
                (STARLIGHT[] !== nothing && !(all(isfinite, xs1) && isfinite(Ts1))) && error("starlight: non-finite chemistry in zone $k: x $xs1 T $Ts1 from x $(m.x[k]) T $(tn[k]) ρ $(dn[k]) rates $(slrates(k))")
                (xs1, Ts1)
            catch e
                # Fixed mesh: e and x are advected by different schemes, and e may then have no solution for the zone's x. Hot dense gas: table equilibrium, same e.
                STARLIGHT[] !== nothing && !(e isa ErrorException) && println("starlight zone failed: k ", k, " x ", m.x[k], " ρ ", dn[k], " u ", egn[k], " T ", tn[k], " dt ", dt, " rates ", slrates(k), " eq ", m.eqzone[k], " E_r ", TM.v(:ern)[k])
                ((CHEMFALLBACK[] || (CHEMFBDENSE[] && pgn[k] >= m.pswitch)) && (e isa ErrorException || (p3 !== nothing && e isa ArgumentError))) || rethrow()
                lock(() -> (NCHEMFB[] += 1), CHEMLOCK); fb = true
                TL3DBG[] > 0 && @printf("tl3fb k %d r %.3e T %.0f tauz %.2e eqprev %s xe %.2e u %.4e u(10K,x) %.4e u(T,x) %.4e: %s\n", k, rn[k] / 6.957e10, tn[k], tauz[k], eqprev, CNn.full(m.x[k], m.net)[:e], egn[k], CNn.eos(dn[k], 10.0, m.x[k], m.net)[1], CNn.eos(dn[k], tn[k], m.x[k], m.net)[1], first(sprint(showerror, e), 60))
                coolt = cool === nothing ? nothing : ((ρ, T, y) -> cool(ρ, T, nothing))
                _, T2 = EQ.step(dn[k], egn[k], tn[k], dt; cooling = coolt)
                (x_from_table(m, dn[k], T2), T2)
            end
            fb && !isempty(TX.CHEMB[]) && (TX.CHEMBON[][k] = false)
            m.eqzone[k] || lwheat!(k, x1)
            if p3 !== nothing && !m.eqzone[k] && !fb
                cq, yq = tl3_heat(dn[k], T1, CNn.full(x1, m.net), p3e)
                isempty(TX.RAYQ[]) || (TX.RAYQ[][k] += cq)       # the photon energy the chemistry gave the gas leaves the gray field
                isempty(TX.GASQ[]) || (TX.GASQ[][k] += yq)
                SLDIAG[].tl3[k] = merge(p3e, (; chemq = cq, cycq = yq, qbb = lya_bb(dn[k], T1, CNn.full(x1, m.net), p3e)))
                NLOPAC[] && (NLOPC[][k] = nlopc(dn[k], T1, CNn.full(x1, m.net), p3e))
                if NLODAMP < 1 && NLOPAC[] && k <= length(NLOPREV[]) && (o0 = NLOPREV[][k]) !== nothing && (o1 = NLOPC[][k]) !== nothing
                    w_ = NLODAMP; bl(a0, a1) = exp((1 - w_) * log(a0) + w_ * log(a1))
                    NLOPC[][k] = merge(o1, (; RE = bl(o0.RE, o1.RE), RP = bl(o0.RP, o1.RP), RR = bl(o0.RR, o1.RR)))
                end
                if TL3DBG[] > 0 && TF.getscalar(:jstep) % TL3DBG[] == 0 && k % 6 == 0
                    d0 = CNn.full(m.x[k], m.net); d1 = CNn.full(x1, m.net)
                    @printf("tl3dbg k %d r %.3e tauz %.2e T %.0f -> %.0f  xe %.2e -> %.2e  C %.2e β %.2e RL %.2e  chemq %.2e cycq %.2e  u %.3e  u(T0,x0) %.3e  eqprev %s\n", k, rn[k] / 6.957e10, tauz[k], tn[k], T1, d0[:e], d1[:e], tl3_C(T1, p3), p3.β, p3.RL, cq, yq, egn[k], CNn.eos(dn[k], tn[k], m.x[k], m.net)[1], eqprev)
                end
            end
            m.x[k] = x1
        end
        k == ZONEDBG[] && @printf("zonedbg step %d k %d p %.4e T %.1f -> %.1f eq %s tauz %.3e u %.5e dt %.3e\n", TF.getscalar(:jstep), k, pgn[k], tn[k], T1, m.eqzone[k], tauz[k], egn[k], dt)
        dTz[k] = abs(T1 / tn[k] - 1); tn[k] = T1
        nothing
    end
    CHEMZPROF && length(ZPT[]) != length(tn) && (ZPT[] = zeros(length(tn)))
    zrun!(k) = CHEMZPROF ? (t0_ = time_ns(); zone!(k); ZPT[][k] += (time_ns() - t0_) * 1e-9; nothing) : zone!(k)   # CHEMZPROF: time per zone (one task per zone, no race)
    if CHEMTHREADS[] && ray === nothing && Threads.nthreads() > 1    # zones are independent (start-of-step state, per-zone outputs)
        if CHEMSCHED[] === :zone
            next = Threads.Atomic{Int}(a)
            @sync for _ in 1:Threads.nthreads()
                Threads.@spawn while true
                    kz = Threads.atomic_add!(next, 1)
                    kz > b && break
                    zrun!(kz)
                end
            end
        else
            Threads.@threads :dynamic for k in a:b; zrun!(k); end
        end
    else
        for k in a:b; zrun!(k); end
    end
    dTmax = maximum(view(dTz, a:b))
    if CHEMZPROF && (ZPN[] += 1) % 500 == 0           # the zones that took most of the network time over the last 500 steps
        zt = ZPT[]; tot = sum(view(zt, a:b)); js = sort(collect(a:b); by = k -> -zt[k])[1:min(12, b - a + 1)]
        @printf("CHEMZPROF step %d: %.1f s in zones over 500 steps (%d table, %d network now); top zones:\n", TF.getscalar(:jstep), tot, count(k -> m.eqzone[k], a:b), count(k -> !m.eqzone[k], a:b))
        for k in js
            @printf("    zone %3d  %5.1f %%  r %.3e R_sun  n %.3e  T %8.1f  p/pswitch %.3g  %s\n", k - a + 1, 100 * zt[k] / tot, TM.v(:rn)[k] / 6.957e10, dn[k] / CNn.μH(), tn[k], pgn[k] / m.pswitch, m.eqzone[k] ? "table" : "network")
        end
        fill!(zt, 0.0)
    end
    for k in a-2:a-1; k >= 1 && (m.x[k] = m.x[a]; m.eqzone[k] = m.eqzone[a]); end
    for k in b+1:b+2; m.x[k] = m.x[b]; m.eqzone[k] = m.eqzone[b]; end
    tn[a-1] = tn[a]; tn[b+1] = tn[b]; tn[b+2] = tn[b]
    TM.fill_eos!(m)
    dTmax
end

# ------------------------------------------------------------------------------------------ deuterium and HD (HDCHEM=1)
"""
HDCHEM=1 (default off): D, D⁺ and HD as a trace network carried with the H network, and optically thin HD cooling.

State per zone: (x_D⁺, x_HD) per H nucleus, appended to the zone's `m.x[k]` behind the free species of the network (the advection
and every remap treat `m.x[k]` as a vector; `ChemNet.full` and `fullvec` read only the free species). x_D = y_D − x_D⁺ − x_HD with
y_D = D/H by number (ENV `DH`, the variable the nuclear network uses; 2.5e-5).
After the chemistry–temperature step of a zone the three-state linear system takes one backward-Euler step at the zone's new T and
composition: rates k50–k56 of ChemistryKernels.jl (Savin 2002; Galli & Palla 2002; Clark et al. 2011; Kreckel et al. 2010 form for
D + H⁻), D ionisation and recombination with the H rates, HD photodissociation with the zone's Lyman–Werner rate. The step is written
as the stationary state of the chain with a restart at rate 1/Δt to the old state: sums of positive products only, exact for any
rate × Δt, and y_D is conserved. Table and hot zones: the stationary state.
Cooling: n_HD [W(T) − W(T_CMB)], W = Λ_LTE/(1 + Λ_LTE/(n_HI Λ_low)) (Coppola et al. 2011; Wrathmall, Gusdorf & Flower 2007: HD–H
collisions only), inside the implicit chemistry–temperature step, with the x_HD of the start of the step.
Zones without the two entries (restart from a run without HD, a rebuilt mesh): table zones get the stationary state; network zones
the stationary state with `HDINIT=eq` (default), or with `HDINIT=<x_HD>` that HD abundance and D⁺/D = H⁺/H.
Not included: HD in the EOS (binding energy, particle number: 2.5e-5 of H), HD–H₂ and HD–He collisions, HD line optical depths,
pumping of HD by the star, the nuclear network's deuterium (the star's zones are fully ionised, so their D is D⁺ here either way).
"""
const HDCHEM = Ref(get(ENV, "HDCHEM", "0") == "1")
const HDYD = Ref(parse(Float64, get(ENV, "DH", "2.5e-5")))
const HDINIT = Ref(get(ENV, "HDINIT", "eq"))
const HDTCMB = Ref(parse(Float64, get(ENV, "TCMB", "0")))
const HDX = Ref(Float64[])                 # per zone: x_HD at the start of the step (the cooling term of zone_cooling)
const HDLAM = Ref(Float64[])               # per zone: HD cooling at the end of the step [erg cm⁻³ s⁻¹] (diagnostic)
const HDNPAD = Ref(0)
hd_k50(T) = T <= 2e5 ? max(0.0, 2.0e-10 * T^0.402 * exp(-37.1 / T) - 3.31e-17 * T^1.48) : 2.5e-8 * (T / 2e5)^0.402      # H⁺ + D → H + D⁺
hd_k51(T) = 2.06e-10 * T^0.396 * exp(-33.0 / T) + 2.03e-9 * T^-0.332                                                    # H + D⁺ → H⁺ + D
hd_k52(T) = T <= 1e4 ? (l = log10(T); max(0.0, 1.0e-9 * (0.417 + 0.846l - 0.137l^2))) : 1.609e-9                        # H₂ + D⁺ → HD + H⁺
hd_k53(T) = 1.1e-9 * exp(-488.0 / T)                                                                                    # HD + H⁺ → H₂ + D⁺
hd_k54(T) = T <= 2e3 ? (l = log10(T); exp10(-56.4737 + l * (5.88886 + l * (7.19692 + l * (2.25069 + l * (-2.16903 + l * 0.317887)))))) :
                       3.17e-10 * exp(-5207.0 / T)                                                                      # H₂ + D → HD + H
# (the published fit turns upward below 78.5 K: below 200 K an Arrhenius tail matched in value and slope, as ChemistryKernels.jl)
hd_k55(T) = T <= 200 ? exp(log(5.25e-11) - 4430.0 / 200 + 173900.0 / 200^2 - (4430.0 - 2 * 173900.0 / 200) * (1 / T - 1 / 200)) :
                       5.25e-11 * exp(-4430.0 / T + 173900.0 / T^2)                                                     # HD + H → H₂ + D
hd_k56(T) = 1.35e-9 * (T^0.098493 + 0.32852 * T^0.5561 + 2.771e-7 * T^2.1826) /
            (1 + 6.191e-3 * T^1.0461 + 8.9712e-11 * T^3.0424 + 3.2576e-14 * T^3.7741)                                   # D + H⁻ → HD + e
"HD cooling per molecule in LTE [erg/s] (Coppola et al. 2011) and the low-density coefficient for HD–H [erg cm³/s] (Wrathmall et al. 2007)."
hd_lte(T) = (t = min(max(T, 10.0), 3.0e4); l = log10(t); exp10(min(-55.5725 + l * (56.649 + l * (-37.9102 + l * (12.698 + l * (-2.02424 + l * 0.122393)))), 0.0)))
hd_low(T) = (t = min(max(T, 10.0), 6.0e3); l = log10(t / 1e3); exp10(-23.175780 + l * (1.5035261 + l * (0.40871403 + l * (0.17849311 + l * (-0.077291388 + l * 0.10031326))))))
hd_w(T, nHI) = (l = hd_lte(T); lo = nHI * hd_low(T); (l > 0 && lo > 0) ? l / (1 + l / lo) : zero(l * lo))
"HD cooling [erg cm⁻³ s⁻¹] of gas with composition d (ChemNet.full; nothing: table zone, no HD) and x_HD; < 0 below T_CMB."
function hd_cool(ρ, T, d, xHD)
    d === nothing && return zero(T)
    nH = ρ / CNn.μH(); nHI = nH * max(d[:H], 0.0); Tc = HDTCMB[]
    nH * xHD * (hd_w(T, nHI) - (Tc > 0 ? hd_w(Tc, nHI) : 0.0))
end
"Transition rates [1/s] between the states (D, D⁺, HD): (D→D⁺, D⁺→D, D⁺→HD, HD→D⁺, D→HD, HD→D)."
function hd_rates(T, nH, d, kLW)
    t = min(max(Float64(T), 10.0), 1e6)
    g(s) = max(Float64(get(d, s, 0.0)), 0.0)
    xH, xp, xe, xH2, xm = g(:H), g(:Hp), g(:e), g(:H2), g(:Hm)
    (nH * (PC.k1_aazn(t) * xe + hd_k50(t) * xp), nH * (PC.k2_aazn(t) * xe + hd_k51(t) * xH), nH * hd_k52(t) * xH2,
     nH * hd_k53(t) * xp, nH * (hd_k54(t) * xH2 + hd_k56(t) * xm), nH * hd_k55(t) * xH + kLW)
end
"""
One backward-Euler step of the chain from p0 = fractions (D, D⁺, HD) of y_D, idt = 1/Δt (0: the stationary state). p(I − Δt Q) = p0 is
the stationary state of the chain with every state restarting at rate 1/Δt into p0; for three states that is a sum over spanning trees.
"""
function hd_step(q, p0, idt)
    q12, q21, q23, q32, q13, q31 = q
    if idt > 0
        isfinite(idt) || return p0
        q21 += idt * p0[1]; q31 += idt * p0[1]; q12 += idt * p0[2]; q32 += idt * p0[2]; q13 += idt * p0[3]; q23 += idt * p0[3]
    end
    w1 = q21 * q31 + q21 * q32 + q23 * q31
    w2 = q12 * q32 + q12 * q31 + q32 * q13
    w3 = q13 * q23 + q13 * q21 + q23 * q12
    s = w1 + w2 + w3
    (isfinite(s) && s > 0) || return p0
    (w1 / s, w2 / s, w3 / s)
end
hd_fracs(xDp, xHD) = (yD = HDYD[]; p2 = clamp(xDp / yD, 0.0, 1.0); p3 = clamp(xHD / yD, 0.0, 1.0 - p2); (1.0 - p2 - p3, p2, p3))
function hd_default(m, k, x)
    dn, tn = TM.v(:dn), TM.v(:tn); yD = HDYD[]
    d = CNn.full(x, m.net); nH = dn[k] / CNn.μH()
    if m.eqzone[k] || HDINIT[] == "eq"
        p = hd_step(hd_rates(tn[k], nH, d, 0.0), (1.0, 0.0, 0.0), 0.0)
        return (yD * p[2], yD * p[3])
    end
    xhd = min(parse(Float64, HDINIT[]), yD); xH = max(Float64(d[:H]), 0.0); xp = max(Float64(d[:Hp]), 0.0)
    ((yD - xhd) * (xH + xp > 0 ? xp / (xH + xp) : 0.0), xhd)
end
"Give every zone its (x_D⁺, x_HD) entries (no-op where they are present)."
function hd_pad!(m::NetGas)
    nf = length(m.net.free); n = 0
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    for k in a:b
        length(m.x[k]) == nf || continue
        p = hd_default(m, k, m.x[k]); m.x[k] = vcat(m.x[k], p[1], p[2]); n += 1
    end
    for k in eachindex(m.x)                                        # ghosts and the unused tail of the arrays: as the nearest zone
        (a <= k <= b || length(m.x[k]) != nf) && continue
        j = k < a ? a : b; m.x[k] = vcat(m.x[k], m.x[j][nf+1], m.x[j][nf+2])
    end
    n > 0 && (HDNPAD[] += n; @printf("HDCHEM: (x_D+, x_HD) set in %d zones (HDINIT=%s, D/H = %.3g)\n", n, HDINIT[], HDYD[]))
    m
end
"""
Restart from a checkpoint written without `HECHEM`: insert the He⁺ and He²⁺ fractions (Saha at the zone's ρ, T from the
equilibrium table) in front of the HD tail of every zone's vector (no-op where they are present). Ghost entries use the nearest zone.
"""
function he_pad!(m::NetGas)
    CNn.HECHEM || return m
    i1 = findfirst(==(:Hep), m.net.free); i1 === nothing && return m
    nf = length(m.net.free); want = nf + (HDCHEM[] ? 2 : 0); a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dn, tn = TM.v(:dn), TM.v(:tn); n = 0
    for k in eachindex(m.x)
        length(m.x[k]) == want - 2 || continue
        j = clamp(k, a, b); z = EQ.composition_he(dn[j], min(tn[j], 1e6))
        m.x[k] = vcat(m.x[k][1:i1-1], z[1], z[2], m.x[k][i1:end]); a <= k <= b && (n += 1)
    end
    n > 0 && @printf("HECHEM: (x_He+, x_He++) set from Saha in %d zones\n", n)
    m
end
hd_strip!(m::NetGas) = (nf = length(m.net.free); for k in eachindex(m.x); length(m.x[k]) > nf && (m.x[k] = m.x[k][1:nf]); end; m)
"Per-zone HD record for the snapshots (zones a:b): x_D⁺, x_HD and the HD cooling."
function hd_snap(m::NetGas, a, b)
    hd_pad!(m); nf = length(m.net.free); lam = HDLAM[]
    (; Dp = [m.x[k][nf+1] for k in a:b], HD = [m.x[k][nf+2] for k in a:b], Λ = length(lam) >= b ? lam[a:b] : zeros(b - a + 1))
end
function chem!(m::NetGas)
    nf = length(m.net.free)
    if !HDCHEM[]
        any(x -> length(x) != nf, m.x) && hd_strip!(m)             # (restart from a run with HDCHEM=1)
        return chem_core!(m)
    end
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    dn, tn = TM.v(:dn), TM.v(:tn)
    hd_pad!(m)
    pas = [(x[nf+1], x[nf+2]) for x in m.x]
    hd_strip!(m)                                                   # the network step sees its own free species only
    HDX[] = [max(p[2], 0.0) for p in pas]
    dT = chem_core!(m)
    yD = HDYD[]; lam = zeros(length(m.x)); sd = STARLIGHT[] === nothing ? nothing : SLDIAG[]
    for k in a:b
        d = CNn.full(m.x[k], m.net); nH = dn[k] / CNn.μH()
        kLW = (sd === nothing || k > length(sd.kLW)) ? 0.0 : sd.kLW[k]
        p0 = hd_fracs(pas[k]...)
        p = hd_step(hd_rates(tn[k], nH, d, kLW), p0, m.eqzone[k] ? 0.0 : (dt > 0 ? 1 / dt : Inf))
        m.x[k] = vcat(m.x[k], yD * p[2], yD * p[3])
        lam[k] = Float64(hd_cool(dn[k], tn[k], d, yD * p[3]))
    end
    for k in 1:a-1; m.x[k] = copy(m.x[a]); end
    for k in b+1:length(m.x); m.x[k] = copy(m.x[b]); end
    HDLAM[] = lam
    dT
end

end # module
