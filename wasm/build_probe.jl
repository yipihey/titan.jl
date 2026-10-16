using WasmTarget, ForwardDiff, Serialization
include(joinpath(@__DIR__, "qsw.jl")); include(joinpath(@__DIR__, "probe.jl")); const Q = QSW; const P = Probe; const P2 = Probe2; const P3 = Probe3; const R2 = Roots2
lT, lR, lk, lkP = deserialize("/home/claude/titan.jl/steady/output/kaptable_Z0_v2.jls")
mkpath("out"); open(io -> write(io, vec(permutedims(lk)), vec(permutedims(lkP))), "out/kap64.bin", "w")   # Float64 table for Node, [i][j]
Q.table_init(Int32(941), Int32(521), 3.3, 0.0050000000000003375, -12.0, 0.025000000000000355)
for i in 1:941, j in 1:521; Q.table_set(Int32(i), Int32(j), lk[i, j], lkP[i, j]); end
println("native thermo_nad=", repr(P.thermo_nad(10.0, 7000.0)), " mlt=", repr(P.mlt_nabla(1e4 + 0.0, 7.5e4, 7.5e15, 1.9e39)),
        " Pb=", repr(P.bc_Pb(8e15, 0.0)), " tau=", repr(P.bc_tau(8e15, 1.0)), " integ=", repr(P.integ(8e15, P.bc_Pb(8e15, 0.0), 7300.0)))
println("native rhs_dr=", repr(P2.rhs_dr(log(8.0), 8e15, 1.98e39, log(7300.0))), " rhs_rho=", repr(P2.rhs_rho(log(8.0), 8e15, 1.98e39, log(7300.0))), " mlt_full=", repr(P2.mlt_full(1e7, 2.6e5, 3.4e15, 1.6e39)))
println("native integd=", repr(P3.integd(8e15, P.bc_Pb(8e15, 0.0), 7300.0)))
for (f, md) in ((0.1, 0.0), (0.1, 5.0), (0.01, 10.0)); t = @elapsed n = R2.all_roots2(f, md, 0.8, 1e-3); t = @elapsed n = R2.all_roots2(f, md, 0.8, 1e-3); println("native all_roots2($f,$md) ", [repr(R2.root2(Int32(k))) for k in 1:n], " ", round(t, digits = 2), " s"); end
fill!(R2.ROOTS, 0.0); R2.NROOT[] = 0; fill!(P3.MARK, 0)
Q.TAB[] = Q.KapTable(zeros(1), zeros(1), zeros(1, 1), zeros(1, 1))
fs = [(Q.kap_R, (Float64, Float64)), (R2.all_roots2, (Float64, Float64, Float64, Float64)), (R2.root2, (Int32,)), (P3.integd, (Float64, Float64, Float64)), (P3.mark, (Int32,)), (P2.rhs_dr, (Float64, Float64, Float64, Float64)), (P2.rhs_rho, (Float64, Float64, Float64, Float64)), (P2.mlt_full, (Float64, Float64, Float64, Float64)), (P.thermo_nad, (Float64, Float64)), (P.mlt_nabla, (Float64, Float64, Float64, Float64)),
      (P.bc_Pb, (Float64, Float64)), (P.bc_tau, (Float64, Float64)), (P.integ, (Float64, Float64, Float64)),
      (Q.table_init, (Int32, Int32, Float64, Float64, Float64, Float64)), (Q.table_set, (Int32, Int32, Float64, Float64))]
b = compile_multi(fs); write("out/probe.wasm", b); println("OK ", length(b))
