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
using ForwardDiff
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
"After a converged hydro step: network step per zone in the ideal regime, equilibrium-table step where p ≥ pswitch."
function chem!(m::NetGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    dn, tn, egn, pgn = TM.v(:dn), TM.v(:tn), TM.v(:egn), TM.v(:pgn)
    dTmax = 0.0
    h = HOT[]
    for k in a:b
        if h !== nothing && tn[k] > h.T1             # hot branch: no chemistry, no re-solve of T with the table (driver: burn and mix)
            m.x[k] = x_from_table(m, dn[k], min(tn[k], 1e6)); m.eqzone[k] = true
            continue
        end
        cool = zone_cooling(m, k)
        if pgn[k] >= m.pswitch || (m.eqzone[k] && pgn[k] >= PSWHYS[] * m.pswitch)
            coolt = cool === nothing ? nothing : ((ρ, T, y) -> cool(ρ, T, nothing))
            _, T1 = EQ.step(dn[k], egn[k], tn[k], dt; cooling = coolt)
            m.x[k] = x_from_table(m, dn[k], T1); m.eqzone[k] = true
        else
            m.eqzone[k] && (m.eqzone[k] = false)     # (back to the ideal regime: the network continues from the table composition)
            x1, T1 = CNn.step(m.x[k], dn[k], egn[k], tn[k], dt, m.net; env = m.env, photos = m.photos, cooling = cool)
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
