# Item 4: TitanResidFast vs Fortran matgen and vs TitanResid (reference) on mid-Newton states;
# Jacobian, rhs and side outputs; timing and allocations per call.
# usage: julia --project=env fast_check.jl NAME NN STEP
using Printf, ForwardDiff
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "TitanResid.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl"))
using .TitanDriver; const TD = TitanDriver; const TF = TD.TF; const TR = TitanResid; const TX = TitanResidFast
name, nn, nstep = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3])
TD.bitwise!(true)
TF.init(nn; rundir = joinpath(TF.ROOT, "lib", "p4_runs", "fc_$name"))
for i in 1:nstep-1; TD.titan_step!(); end
TD.r(:peq); TD.r(:matgen); TD.r(:matslv); TD.r(:maxdel); TD.r(:update)
S = TD.snapshot()
SIDE = (:dsn, :usn, :egrsn, :ersn, :frsn, :egrb, :erb, :frb, :egsn, :egb, :qf, :dudr, :db, :ub, :xnu)
grab() = ([copy(TF.var(s)) for s in (:em2, :em1, :e00, :ep1, :ep2)], copy(TF.var(:rhs)),
          Dict(s => copy(TF.var(s)) for s in SIDE))
TD.fcall(:matgen_); F = grab(); TD.restore!(S)
TR.matgen!(); A = grab(); TD.restore!(S)
TX.matgen!(ncol = 7); B7 = grab(); TD.restore!(S)
TX.matgen!(ncol = 6); B5 = grab(); TD.restore!(S)
ngrs = Int(TF.getscalar(:ngrs)); ngre = Int(TF.getscalar(:ngre)); neqn = Int(TF.getscalar(:neqn)); K = ngrs-2:ngre+2
function cmp(X, Y)
    jr = 0.0; rr = 0.0
    for i in 1:neqn, k in K
        rowmax = maximum(abs(Y[1][o+3][i, j, k]) for o in -2:2, j in 1:neqn)
        rr = max(rr, abs(X[2][i, k] - Y[2][i, k]) / max(abs(Y[2][i, k]), rowmax, floatmin()))
        for o in -2:2, j in 1:neqn
            jr = max(jr, abs(X[1][o+3][i, j, k] - Y[1][o+3][i, j, k]) / max(rowmax, floatmin()))
        end
    end
    sd = maximum(s -> (a = X[3][s][K]; b = Y[3][s][K]; m = maximum(abs, b); m == 0 ? maximum(abs, a) : maximum(abs.(a .- b)) / m), SIDE)
    (jr, rr, sd)
end
@printf("%s step %d mid-Newton:\n", name, nstep)
for (lab, X, Y) in (("Fast7 vs Ref", B7, A), ("Fast6 vs Ref", B5, A), ("Fast7 vs Fortran", B7, F), ("Fast6 vs Fortran", B5, F), ("Ref vs Fortran", A, F))
    jr, rr, sd = cmp(X, Y)
    @printf("  %-17s Jacobian max|Δ|/rowmax %.2e  rhs %.2e  side outputs (max rel) %.2e\n", lab, jr, rr, sd)
end
function tm(f, n)
    f(); GC.gc(); a = @allocated f()
    t = @elapsed for i in 1:n; f(); end
    (1e3t / n, a)
end
for (lab, f) in (("Fortran matgen", () -> TD.fcall(:matgen_)), ("TitanResid (ref)", () -> TR.matgen!()),
                 ("Fast ncol=7", () -> TX.matgen!(ncol = 7)), ("Fast ncol=6", () -> TX.matgen!(ncol = 6)))
    ms, a = tm(f, 50); TD.restore!(S)
    @printf("  %-18s %8.3f ms/call  %10d bytes allocated/call\n", lab, ms, a)
end

for s in SIDE
    a = A[3][s][K]; b = F[3][s][K]; c7 = B7[3][s][K]; m = max(maximum(abs, b), floatmin())
    @printf("  side %-6s Ref-Fortran %.2e   Fast-Fortran %.2e\n", s, maximum(abs.(a .- b)) / m, maximum(abs.(c7 .- b)) / m)
end
TF.finalize()
