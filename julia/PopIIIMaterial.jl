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

@inline function TX.material_gas(m::NetGas, ρ, T)
    k = TX.ZONE[]
    if m.eqzone[k]
        u, p = EQ.eos(ρ, T)
    else
        u, p = CNn.eos(ρ, T, m.x[k], m.net)
    end
    return p, u, sqrt(5.0 / 3.0 * p / ρ)
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
TM.fill_opac!(::NetGas) = nothing

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

"After a converged hydro step: network step per zone in the ideal regime, equilibrium-table step where p ≥ pswitch."
function chem!(m::NetGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
    dn, tn, egn, pgn = TM.v(:dn), TM.v(:tn), TM.v(:egn), TM.v(:pgn)
    dTmax = 0.0
    for k in a:b
        cool = zone_cooling(m, k)
        if pgn[k] >= m.pswitch
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
