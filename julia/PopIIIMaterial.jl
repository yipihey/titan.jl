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
    f === nothing && return m.cooling
    m.cooling === nothing ? ((ρ, T, y) -> f(T)) : ((ρ, T, y) -> m.cooling(ρ, T, y) + f(T))
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
@inline function TX.material_opac(m::NetGas, ρ, T)
    κR, κP = OPAC[](ρ, T)
    return κR, κP, κP
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
        xnen[k] = (m.eqzone[k] ? EQ.composition(dn[k], min(tn[k], 1e6))[1] : CNn.full(m.x[k], m.net)[:Hp]) * dn[k] / CNn.μH()
    end
    nothing
end

"Free-species vector from the equilibrium table's composition (H₂, H⁺; H⁻ and Li⁺ in Saha with the table's electrons)."
function x_from_table(m::NetGas, ρ, T)
    y = EQ.composition(ρ, T); nH = ρ / CNn.μH()
    xH2 = y[3]; xHp = y[1]; xH = 1 - xHp - 2xH2
    ne = xHp * nH
    xHm = xH * ne / CNn.Keq([:Hm], [:H, :e], T)
    out = zeros(length(m.net.free))
    for (i, s) in enumerate(m.net.free)
        out[i] = s === :Hp ? xHp : s === :H2 ? xH2 : s === :Hm ? xHm :
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
PSWTAU[] > 0: a zone with T >= PSWTEX[] and Rosseland depth from the surface >= PSWTAU[] (half that to stay) uses the table whatever its pressure. Where radiation is thermalised, J = B holds H at Saha.
The network has no radiation field: at the outer edge of an opaque layer (p ~ 6, T ~ 1e4 K, 39 M_sun) it sits at coronal x(H+) ~ 0.05 against 0.99, Newton heats the zone by a factor of 2 at frozen
composition and chem! brings it back each step, until Newton fails. On the table the ionization heat capacity is inside Newton. The switch at fixed e lowers T once; radiation reheats the zone.
"""
const PSWTAU = Ref(0.0)
"After a converged hydro step: network step per zone in the ideal regime, equilibrium-table step where p ≥ pswitch."
const CHEMFALLBACK = Ref(false); const NCHEMFB = Ref(0)
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
    chemq, cycq
end
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
    function rates(k)
        if k > kph && Te > 0
            Tc = Te; W = max(CNn.c * (er[k] - Ecmb), 0.0) / (4 * 5.670374419e-5 * Te^4)
        else
            Tc = (max(er[k], Ecmb) / arad)^0.25; W = 1.0
        end
        (W > 0 && Tc > 0) || return (0.0, 0.0, 0.0, 1.0)
        env = (Trad = Tc, J = ν -> W * _bnu(ν, Tc))
        kHm = sl.pHm === nothing ? 0.0 : CNn.photorates(sl.pHm, tn[k], env)[1]; kH = sl.pH === nothing ? 0.0 : CNn.photorates(sl.pH, tn[k], env)[1]
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
            diag.tl3[k] = (; β = W * beta2s_bb(Tc), RL, ne = nH * max(d[:e], 0.0), P12 = (RL + LAM2S) * nb, hν = balmer_mean(Tc))
        end
        (kHm, kH, kLW, f)
    end
    rates
end
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
function chem!(m::NetGas)
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
    ph3 = filter(p -> p.R !== :H, m.photos)                                 # three-level zones: no Boltzmann-level H photo pair
    for k in a:b
        if h !== nothing && tn[k] > h.T1             # hot branch: no chemistry, no re-solve of T with the table (driver: burn and mix)
            m.x[k] = x_from_table(m, dn[k], min(tn[k], 1e6)); m.eqzone[k] = true
            continue
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
                dTmax = max(dTmax, abs(T1 / tn[k] - 1)); tn[k] = T1
                continue
            end
            Nin = nothing                                # the rest of the light cannot change any zone further out
        end
        # PSWTMIN[] > 0: also require T above it (5 % hysteresis). Below ~5000 K ionization holds negligible energy, so the network and the table agree there
        # to ~1 % in T at fixed e; at 7000-9500 K the network's H ionization lags the table's Saha by ×5, and switching there drops T by up to 16 %.
        if (pgn[k] >= m.pswitch && tn[k] >= PSWTMIN[]) || (m.eqzone[k] && pgn[k] >= PSWHYS[] * m.pswitch && tn[k] >= 0.95 * PSWTMIN[]) || (m.eqzone[k] && PSWTEX[] > 0 && tn[k] >= PSWTEX[]) || (PSWTAU[] > 0 && tn[k] >= PSWTEX[] && tauz[k] >= (m.eqzone[k] ? 0.5 : 1.0) * PSWTAU[])
            coolt = cool === nothing ? nothing : ((ρ, T, y) -> cool(ρ, T, nothing))
            _, T1 = EQ.step(dn[k], egn[k], tn[k], dt; cooling = coolt)
            m.x[k] = x_from_table(m, dn[k], T1); m.eqzone[k] = true
        else
            eqprev = m.eqzone[k]
            m.eqzone[k] && (m.eqzone[k] = false)     # (back to the ideal regime: the network continues from the table composition)
            exs, cools = slterms(k, nothing, cool)
            p3 = (slrates === nothing || m.eqzone[k]) ? nothing : SLDIAG[].tl3[k]
            envk, phk = m.env, m.photos; fb = false
            if p3 !== nothing
                envk = merge(m.env, (; tabkey = m.env, tl3 = p3)); phk = ph3
                c0 = cools
                cools = (ρ, T, d) -> (c0 === nothing ? zero(T) : c0(ρ, T, d)) - tl3_heat(ρ, T, d, p3)[1]
            end
            x1, T1 = try
                # u below the zone's own composition at 10 K (x and e advected by different schemes; with x_e ~ 0.1-0.5 from the three-level
                # atom the mismatch reaches ~1 eV per H): no solution for this x, take the fallback below
                p3 !== nothing && egn[k] < CNn.eos(dn[k], 10.0, m.x[k], m.net)[1] && error("tl3: u below u(10 K, x) in zone $k")
                xs1, Ts1 = CNn.step(m.x[k], dn[k], egn[k], tn[k], dt, m.net; env = envk, photos = phk, cooling = cools, extra = exs)
                (STARLIGHT[] !== nothing && !(all(isfinite, xs1) && isfinite(Ts1))) && error("starlight: non-finite chemistry in zone $k: x $xs1 T $Ts1 from x $(m.x[k]) T $(tn[k]) ρ $(dn[k]) rates $(slrates(k))")
                (xs1, Ts1)
            catch e
                # Fixed mesh: e and x are advected by different schemes, and e may then have no solution for the zone's x. Hot dense gas: table equilibrium, same e.
                STARLIGHT[] !== nothing && !(e isa ErrorException) && println("starlight zone failed: k ", k, " x ", m.x[k], " ρ ", dn[k], " u ", egn[k], " T ", tn[k], " dt ", dt, " rates ", slrates(k), " eq ", m.eqzone[k], " E_r ", TM.v(:ern)[k])
                (CHEMFALLBACK[] && (e isa ErrorException || (p3 !== nothing && e isa ArgumentError))) || rethrow()
                NCHEMFB[] += 1; fb = true
                coolt = cool === nothing ? nothing : ((ρ, T, y) -> cool(ρ, T, nothing))
                _, T2 = EQ.step(dn[k], egn[k], tn[k], dt; cooling = coolt)
                (x_from_table(m, dn[k], T2), T2)
            end
            m.eqzone[k] || lwheat!(k, x1)
            if p3 !== nothing && !m.eqzone[k] && !fb
                cq, yq = tl3_heat(dn[k], T1, CNn.full(x1, m.net), p3)
                isempty(TX.RAYQ[]) || (TX.RAYQ[][k] += cq)       # the photon energy the chemistry gave the gas leaves the gray field
                isempty(TX.GASQ[]) || (TX.GASQ[][k] += yq)
                SLDIAG[].tl3[k] = merge(p3, (; chemq = cq, cycq = yq))
                if TL3DBG[] > 0 && TF.getscalar(:jstep) % TL3DBG[] == 0 && k % 6 == 0
                    d0 = CNn.full(m.x[k], m.net); d1 = CNn.full(x1, m.net)
                    @printf("tl3dbg k %d r %.3e tauz %.2e T %.0f -> %.0f  xe %.2e -> %.2e  C %.2e β %.2e RL %.2e  chemq %.2e cycq %.2e  u %.3e  u(T0,x0) %.3e  eqprev %s\n", k, rn[k] / 6.957e10, tauz[k], tn[k], T1, d0[:e], d1[:e], tl3_C(T1, p3), p3.β, p3.RL, cq, yq, egn[k], CNn.eos(dn[k], tn[k], m.x[k], m.net)[1], eqprev)
                end
            end
            m.x[k] = x1
        end
        dTmax = max(dTmax, abs(T1 / tn[k] - 1)); tn[k] = T1
    end
    for k in a-2:a-1; k >= 1 && (m.x[k] = m.x[a]; m.eqzone[k] = m.eqzone[a]); end
    for k in b+1:b+2; m.x[k] = m.x[b]; m.eqzone[k] = m.eqzone[b]; end
    tn[a-1] = tn[a]; tn[b+1] = tn[b]; tn[b+2] = tn[b]
    TM.fill_eos!(m)
    dTmax
end

end # module
