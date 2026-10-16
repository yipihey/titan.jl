# Tier (a) for the full Julia residual (TitanResid) vs Fortran matgen on a mid-Newton
# benchmark state: residual rows and Jacobian blocks, per equation, incl. phantom rows.
# usage: julia --project=env resid_check.jl NAME NN STEP [dense]
using Printf, ForwardDiff
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "TitanResid.jl"))
using .TitanDriver, .TitanResid; const TD = TitanDriver; const TF = TD.TF; const TR = TitanResid
name, nn, nstep = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]); dense = length(ARGS) >= 4
TD.bitwise!(true)
TF.init(nn; rundir = joinpath(TF.ROOT, "lib", "p4_runs", "rc_$name"))
for i in 1:nstep-1; TD.titan_step!(); end
TD.r(:peq); TD.r(:matgen); TD.r(:matslv); TD.r(:maxdel); TD.r(:update)   # mid-Newton state
S = TD.snapshot()
TD.fcall(:matgen_)
Ef = [copy(TF.var(s)) for s in (:em2, :em1, :e00, :ep1, :ep2)]; Rf = copy(TF.var(:rhs))
TD.restore!(S)
c = TR.context()
TR.assemble!(c)
Ej = [copy(TF.var(s)) for s in (:em2, :em1, :e00, :ep1, :ep2)]; Rj = copy(TF.var(:rhs))
ngrs, ngre, neqn = c.ngrs, c.ngre, c.neqn; K = ngrs-2:ngre+2
names = Dict(c.ix[1] => "grid", c.ix[2] => "mass", c.ix[3] => "contin", c.ix[4] => "gasmom", c.ix[5] => "gasnrg")
c.lrad && (names[c.ix[6]] = "radnrg"; names[c.ix[7]] = "radmom")
@printf("%s step %d mid-Newton (t=%.5e, dt=%.4e): residual and Jacobian, Julia (AD) vs Fortran\n", name, nstep,
        TF.getscalar(:timen), TF.getscalar(:dtime))
worst = Dict{String,Any}()
for i in 1:neqn
    rr = 0.0; jr = 0.0; rk = 0; jk = (0, 0, 0)
    for k in K
        # residual: relative to the row's Jacobian scale × (largest variable change of order 1)
        rowmax = maximum(abs(Ef[o+3][i, j, k]) for o in -2:2, j in 1:neqn)
        sc = max(abs(Rf[i, k]), rowmax, floatmin())
        d = abs(Rj[i, k] - Rf[i, k]) / sc
        d > rr && (rr = d; rk = k)
        for o in -2:2, j in 1:neqn
            dj = abs(Ej[o+3][i, j, k] - Ef[o+3][i, j, k]) / max(rowmax, floatmin())
            dj > jr && (jr = dj; jk = (k, o, j))
        end
    end
    @printf("  %-7s residual max |Δ|/max(|R|,rowmax) = %.2e (k=%d)   Jacobian max |ΔJ|/rowmax = %.2e (k=%d, off=%+d, var=%s)\n",
            names[i], rr, rk, jr, jk[1], jk[2], jk[3] == 0 ? "-" : String(c.var[jk[3]]))
end
if dense let
    N = neqn * length(K)
    Jd = ForwardDiff.jacobian(x -> (R = zeros(eltype(x), neqn, length(K)); TR.residual!(R, c, x); vec(R)), zeros(N))
    # compare colored blocks with dense AD and count dependencies outside the ±2-block band
    out = 0.0; mism = 0.0
    for kb in 1:length(K), i in 1:neqn, cb in 1:length(K), j in 1:neqn
        v = Jd[(kb-1)*neqn+i, (cb-1)*neqn+j]
        if abs(cb - kb) > 2
            rowmax = maximum(abs(Ef[o+3][i, jj, K[kb]]) for o in -2:2, jj in 1:neqn)
            if abs(v) / max(rowmax, floatmin()) > out
                out = abs(v) / max(rowmax, floatmin())
                global OUTWHERE = (names[i], K[kb], cb - kb, String(c.var[j]))
            end
        else
            mism = max(mism, abs(v - Ej[cb-kb+3][i, j, K[kb]]))
        end
    end
    @printf("  dense AD check: max |entry|/rowmax outside ±2-block band = %.2e at %s ; max |colored − dense| (abs) = %.2e\n", out, string(OUTWHERE), mism)
end end
# entry-level detail for the largest Jacobian differences
let det = []
    for i in 1:neqn, k in K
        rowmax = maximum(abs(Ef[o+3][i, j, k]) for o in -2:2, j in 1:neqn)
        for o in -2:2, j in 1:neqn
            a = Ej[o+3][i, j, k]; f = Ef[o+3][i, j, k]
            dj = abs(a - f) / max(rowmax, floatmin())
            dj > 1e-9 && push!(det, (names[i], k, o, String(c.var[j]), a, f, dj))
        end
    end
    sort!(det, by = e -> -e[end])
    println("  Jacobian entries differing by > 1e-9 rowmax: ", length(det))
    for e in det[1:min(end, 12)]
        @printf("    %-7s k=%3d off=%+d d/d%-5s AD %+.6e  Fortran %+.6e  |Δ|/rowmax %.2e\n", e...)
    end
end
TF.finalize()
