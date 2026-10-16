"""
    TitanMaterials

Pluggable gas + opacity physics for the Julia TITAN, used in two places with one definition:
* the residual (`TitanResidFast.material_gas/material_opac`, called with ForwardDiff Duals), and
* the n-level EOS/opacity fills used by peq!/update! (`TitanDriver.EOS[]/OPAC[]`, Float64 on the state arrays).

Materials:
* `Builtin`  — TITAN's own ideal gas (leos = 2) with lopac 3/4 (the benchmarks); `use!(Builtin)` keeps the
  original inline code paths, so benchmark results are unchanged.
* `Z0Gas`    — zero-metallicity LTE Saha gas (H, He, e⁻) with tabulated Rosseland (flux mean, incl. scattering)
  and Planck (absorption) means from Z0Physics. Radiation is NOT in the gas EOS (TITAN evolves E_r separately).

`use!(m)` installs a material everywhere.
"""
module TitanMaterials

import ..TitanDriver, ..TitanResidFast, ..JuliaPhysics, ..Z0Physics
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const JP = JuliaPhysics; const Z0 = Z0Physics

abstract type Material end

struct Builtin <: Material end

"Zero-metallicity Saha gas with tabulated mean opacities."
struct Z0Gas{T} <: Material
    tab::T
end
Z0Gas() = Z0Gas(Z0.table())

# ---------------- residual interface (generic in the number type)
@inline function TX.material_gas(m::Z0Gas, ρ, T)
    s = Z0.saha(ρ, T)
    pg = (s.nH + s.nHe + s.ne) * Z0.kB * T
    eg = Z0._eint(s, ρ, T) - Z0.arad * T^4 / ρ          # gas only: translational + ionization
    # a_s enters only the artificial viscosity: monatomic-gas estimate (see notes)
    return pg, eg, sqrt(5.0 / 3.0 * pg / ρ)
end
@inline function TX.material_opac(m::Z0Gas, ρ, T)
    κR = Z0.kap(m.tab, ρ, T); κP = Z0.kapP(m.tab, ρ, T)
    return κR, κP, κP
end

# ---------------- state-array fills (peq!/update!)
const V = Dict{Symbol,Vector{Float64}}()
v(s::Symbol) = get!(() -> TF.var(s), V, s)::Vector{Float64}
push!(TF.ON_SWITCH, () -> empty!(V))
krange() = (Int(TD.g(:ngrs))::Int - 1):(Int(TD.g(:ngre))::Int + 2)

fill_eos!(::Builtin) = JP.eos!()
fill_opac!(::Builtin) = JP.opac!()

function fill_eos!(m::Z0Gas)
    dn, tn = v(:dn), v(:tn)
    pgn, egn, asn, xnen = v(:pgn), v(:egn), v(:asn), v(:xnen)
    @inbounds for k in krange()
        ρ, T = dn[k], tn[k]
        (ρ > 0 && T > 0) || continue
        pg, eg, as_ = TX.material_gas(m, ρ, T)
        pgn[k] = pg; egn[k] = eg; asn[k] = as_
        xnen[k] = Z0.saha(ρ, T).ne
    end
    nothing
end
function fill_opac!(m::Z0Gas)
    dn, tn = v(:dn), v(:tn)
    chifn, xken, xkpn, xnen = v(:chifn), v(:xken), v(:xkpn), v(:xnen)
    @inbounds for k in krange()
        (dn[k] > 0 && tn[k] > 0) || continue
        χ, κE, κP = TX.material_opac(m, dn[k], tn[k])
        chifn[k] = χ; xken[k] = κE; xkpn[k] = κP
        xnen[k] = Z0.saha(dn[k], tn[k]).ne               # electron density for the scattering source term
    end
    nothing
end

"Install a material in the residual and in the EOS/opacity fills."
function use!(m::Material)
    TX.MATERIAL[] = m isa Builtin ? nothing : m
    TD.EOS[] = () -> fill_eos!(m)
    TD.OPAC[] = () -> fill_opac!(m)
    m
end

end # module
