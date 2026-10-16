# Item 4: BlockSolver vs Fortran matslv and JuliaSolver (gbtrf) on mid-Newton systems.
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaSolver.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
using .TitanDriver; const TD = TitanDriver; const TF = TD.TF; const JS = JuliaSolver; const BS = BlockSolver
name, nn, nstep = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3])
TD.bitwise!(true); TF.init(nn; rundir = joinpath(TF.ROOT, "lib", "p4_runs", "sc_$name"))
for i in 1:nstep-1; TD.titan_step!(); end
TD.r(:peq); TD.r(:matgen); TD.r(:matslv); TD.r(:maxdel); TD.r(:update); TD.r(:matgen)
S = TD.snapshot()
sol(f) = (TD.restore!(S); f(); copy(TF.var(:rhs)))
xf = sol(() -> TD.fcall(:matslv_)); xj = sol(JS.matslv!); xb = sol(BS.matslv!); xp = sol(BS.Penta.matslv!)
a = Int(TF.getscalar(:ngrs)) - 2; b = Int(TF.getscalar(:ngre)) + 2
d(x, y) = maximum(abs.(x[:, a:b] .- y[:, a:b])) / maximum(abs.(y[:, a:b]))
@printf("%s step %d: max|Δx|/max|x|  block vs Fortran %.2e   block vs gbtrf %.2e   gbtrf vs Fortran %.2e   (min scaled pivot %.2e)\n",
        name, nstep, d(xb, xf), d(xb, xj), d(xj, xf), BS.MINPIV[])
@printf("   penta vs Fortran %.2e   penta vs gbtrf %.2e\n", d(xp, xf), d(xp, xj))
for (lab, f) in (("Fortran matslv", () -> TD.fcall(:matslv_)), ("JuliaSolver gbtrf", JS.matslv!), ("Block-tri (pairs)", BS.matslv!), ("Block-penta static", BS.Penta.matslv!))
    TD.restore!(S); f(); TD.restore!(S); GC.gc(); al = @allocated f()
    t = 0.0; for i in 1:200; TD.restore!(S); t += @elapsed f(); end
    @printf("  %-18s %7.3f ms/call  %8d bytes/call\n", lab, 1e3t / 200, al)
end
TF.finalize()
