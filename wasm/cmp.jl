include(joinpath(@__DIR__, "qsw.jl")); include(joinpath(@__DIR__, "probe.jl"))
using Serialization
const Q = QSW
lT, lR, lk, lkP = deserialize("/home/claude/titan.jl/steady/output/kaptable_Z0_v2.jls")
Q.table_init(Int32(941), Int32(521), 3.3, 0.0050000000000003375, -12.0, 0.025000000000000355)
for i in 1:941, j in 1:521; Q.table_set(Int32(i), Int32(j), lk[i, j], lkP[i, j]); end
p = Q.params(0.1, 0.8, 1e-3); Md = 5.0 * Q.Msun / Q.yr
for R in (8.26785524971e15, 8.26785524972e15, 8.26785524973e15, 8.26785524974e15, 8e15, 3e15)
    bc = Q.outer_bc(R, p; Mdot = Md)
    st, m = Probe3.integ_dbg(R, bc.Pb, bc.Tb, p)
    println(repr(R), " ", repr(bc.Pb), " ", repr(bc.Tb), " ", st, " ", repr(m))
end
