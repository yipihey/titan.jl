# PROTO for grid.jl: the network's means with the G2 additions (G2Opacity.means, the in-tree code), LTE populations.
include(joinpath(@__DIR__, "load.jl")); include(joinpath(@__DIR__, "mats.jl"))
include(joinpath(JD, "G2Opacity.jl")); const G2O = G2Opacity
const PROTO = (ρ, T, s) -> begin
    n2 = G2O.n2_lte(T, s.nH, s.ne, s.nH + s.nHe + s.nH2)
    E, P, R = G2O.means(ρ, T; n1 = s.nH, n2, nH2 = s.nH2, nHe = s.nHe, nHep = s.nHep, nHepp = s.nHepp, nHp = s.nHp, ne = s.ne, nHm = s.nHm, nHmS = s.nHm, hbfP = true)
    (; κR = R, κP = P)
end
G2O.freq()
include(joinpath(@__DIR__, "grid.jl"))
