module QuasiStarInfall
using Serialization
include("constants.jl"); include("eos.jl"); include("opacity.jl"); include("envelope.jl"); include("boundary.jl")
const TABFILE = joinpath(@__DIR__, "..", "output", "kaptable_Z0_v2.jls")
const TABFILE_V1 = joinpath(@__DIR__, "..", "output", "kaptable_Z0.jls")
load_table() = isfile(TABFILE) ? KapTable(deserialize(TABFILE)...) : KapTable(deserialize(TABFILE_V1)...)
export QSParams, load_table, solve_R, outer_bc, integrate_inward, shoot_sign, LEdd_es, Msun, Lsun, Rsun, yr, G, σSB,
       kap, kapP, kap_bilinear, thermo, saha, eos_rhoT, eint, entropy, rho_PT, Ptot, kappa_means, in_table, vff, rho_infall
end
