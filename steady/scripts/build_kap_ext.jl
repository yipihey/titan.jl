# Extended table (v3): log T 3.0–8.0, log R −20–1 (same spacing as v2, so v2's nodes are a subset).
# Needed by the TITAN infall runs: cold (< 2000 K) and rarefied (log R < −12) inflow gas.
# No molecules (H2 CIA) are included: below ~2000 K this table is a lower bound.
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; const Q = QuasiStarInfall
using Serialization
@time tab = Q.build_table(lT = 3.0:0.005:8.0, lR = -20.0:0.025:1.0)
bad = count(!isfinite, tab.lk) + count(!isfinite, tab.lkP)
serialize(joinpath(@__DIR__, "..", "output", "kaptable_Z0_v3ext.jls"), (tab.lT, tab.lR, tab.lk, tab.lkP))
println("table ", size(tab.lk), " log κR ", extrema(filter(isfinite, tab.lk)), " log κP ", extrema(filter(isfinite, tab.lkP)), " nonfinite ", bad)
