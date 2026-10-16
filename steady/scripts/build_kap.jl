push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; const Q = QuasiStarInfall
using Serialization
@time tab = Q.build_table(lT = 3.3:0.005:8.0, lR = -12.0:0.025:1.0)
serialize(Q.TABFILE, (tab.lT, tab.lR, tab.lk, tab.lkP))
println("table ", size(tab.lk), " log κR ", extrema(tab.lk), " log κP ", extrema(tab.lkP))
