"""
    StellarMaterial

S2: `StellarGas`, the TITAN material built on `StellarEOS` (degenerate/relativistic electrons, Saha join, composition per
zone). Per-zone mass fractions X[k], Y[k], Z[k] (k = zone index, as PopIIIMaterial's `y`); the residual reads the zone from
`TitanResidFast.ZONE[]` (set by `_material!` before every material call; also set in the MLT coefficient loops).
The gas is gas-only (radiation is TITAN's E_r). Opacity: PLACEHOLDER = the Z = 0 tables of Z0Physics (composition-independent)
until S3 supplies stellar opacities; `material_thermo` returns ∇_ad, c_P, δ of gas + equilibrium radiation.
Use: `m = StellarGas(mgr; X, Y, Z)`; `TitanMaterials.use!(m)`; `setcomp!(m, k, X, Y, Z)` for a zone;
`ghosts!(m)` copies the first/last real zone into the ghosts.
"""
module StellarMaterial

import ..TitanDriver, ..TitanResidFast, ..TitanMaterials, ..Z0Physics, ..StellarEOS
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const TM = TitanMaterials
const SE = StellarEOS; const Z0 = Z0Physics

mutable struct StellarGas{T} <: TM.Material
    X::Vector{Float64}; Y::Vector{Float64}; Z::Vector{Float64}
    zmix::NamedTuple{(:C, :N, :O),NTuple{3,Float64}}
    tab::T                  # opacity table (placeholder: Z0Physics)
    grey::Bool
    coulomb::Bool; pairs::Bool
end
function StellarGas(n::Integer; X = 0.76, Y = 0.24, Z = 0.0, zmix = SE.ZMIX, tab = Z0.table(), grey = false,
                    coulomb = true, pairs = true)
    StellarGas(fill(Float64(X), n), fill(Float64(Y), n), fill(Float64(Z), n), zmix, tab, grey, coulomb, pairs)
end
setcomp!(m::StellarGas, k, X, Y, Z) = (m.X[k] = X; m.Y[k] = Y; m.Z[k] = Z; m)
function ghosts!(m::StellarGas)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    for k in max(a - 2, 1):a-1; setcomp!(m, k, m.X[a], m.Y[a], m.Z[a]); end
    for k in b+1:min(b + 2, length(m.X)); setcomp!(m, k, m.X[b], m.Y[b], m.Z[b]); end
end
@inline zonecomp(m::StellarGas) = (k = TX.ZONE[]; SE.composition(m.X[k], m.Y[k], m.Z[k], m.zmix))

@inline function TX.material_gas(m::StellarGas, ρ, T)
    g = SE.gas_state(ρ, T, zonecomp(m); pairs = m.pairs, coulomb = m.coulomb)
    # a_s enters only the artificial viscosity: γ = 5/3 estimate (as the other materials)
    return g.p, g.e, sqrt(5.0 / 3.0 * g.p / ρ)
end
@inline function TX.material_gasne(m::StellarGas, ρ, T)
    g = SE.gas_state(ρ, T, zonecomp(m); pairs = m.pairs, coulomb = m.coulomb)
    return g.p, g.e, sqrt(5.0 / 3.0 * g.p / ρ), g.ne
end
@inline function TX.material_opac(m::StellarGas, ρ, T)
    κR = Z0.kap(m.tab, ρ, T)
    m.grey && return κR, κR, κR
    κP = Z0.kapP(m.tab, ρ, T)
    return κR, κP, κP
end
"∇_ad, c_P, δ of gas + equilibrium radiation (Float64) for the MLT flux; composition of zone `ZONE[]`."
function TX.material_thermo(m::StellarGas, ρ, T)
    st = SE.eos(ρ, T; comp = zonecomp(m), pairs = m.pairs, coulomb = m.coulomb, radiation = true)
    return st.∇ad, st.cP, st.δ
end

function TM.fill_eos!(m::StellarGas)
    dn, tn = TM.v(:dn), TM.v(:tn)
    pgn, egn, asn, xnen = TM.v(:pgn), TM.v(:egn), TM.v(:asn), TM.v(:xnen)
    @inbounds for k in TM.krange()
        ρ, T = dn[k], tn[k]
        (ρ > 0 && T > 0) || continue
        TX.ZONE[] = k
        pg, eg, as_, ne = TX.material_gasne(m, ρ, T)
        pgn[k] = pg; egn[k] = eg; asn[k] = as_; xnen[k] = ne
    end
    nothing
end
function TM.fill_opac!(m::StellarGas)
    dn, tn = TM.v(:dn), TM.v(:tn)
    chifn, xken, xkpn, xnen = TM.v(:chifn), TM.v(:xken), TM.v(:xkpn), TM.v(:xnen)
    @inbounds for k in TM.krange()
        (dn[k] > 0 && tn[k] > 0) || continue
        TX.ZONE[] = k
        χ, κE, κP = TX.material_opac(m, dn[k], tn[k])
        chifn[k] = χ; xken[k] = κE; xkpn[k] = κP
        xnen[k] = TX.material_gasne(m, dn[k], tn[k])[4]
    end
    nothing
end

end # module
