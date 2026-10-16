# usage: julia --project=env phase3_run.jl NAME NN MODE [tierb]
#   MODE = fortran : Julia driver, all-Fortran routines (timing baseline)
#          jsolve  : matslv -> Julia pivoted banded LU
#          jall    : jsolve + continuity rows assembled in Julia (ForwardDiff)
#                    + Julia eos/opac (hooks)
#   tierb : additionally do tier-(b) one-step comparisons MODE vs fortran at checkpoints
# Writes snaps to lib/p2_runs/NAME_MODE/snaps.jls (for tierc.jl) and a log + timing table.
using Serialization, Printf
include(joinpath(@__DIR__, "TitanDriver.jl"))
include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "JuliaSolver.jl"))
include(joinpath(@__DIR__, "TitanResid.jl"))
include(joinpath(@__DIR__, "Metrics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); const TC = TitanCore
const TX = TitanResidFast; const BS = BlockSolver
using .TitanDriver, .JuliaPhysics, .JuliaSolver, .Metrics, .TitanResid
const TR = TitanResid
const TD = TitanDriver; const TF = TD.TF; const JP = JuliaPhysics; const JS = JuliaSolver

name, nn, mode = ARGS[1], parse(Int, ARGS[2]), ARGS[3]
dotierb = length(ARGS) >= 4 && ARGS[4] == "tierb"
dir = joinpath(TF.ROOT, "lib", "p2_runs", "$(name)_$mode" * (dotierb ? "_tb" : ""))
logio = open(joinpath(TF.ROOT, "lib", "p3_runs", "run_$(name)_$mode" * (dotierb ? "_tb" : "") * ".log"), "w")
say(x...) = (println(x...); println(logio, x...); flush(logio))

# ---------------------------------------------------------- routine sets
const FORT = Dict(s => TD.ROUTINES[s] for s in keys(TD.ROUTINES))
const TIMES = Dict{Symbol,Float64}(); const CALLS = Dict{Symbol,Int}()
const FIRST = Dict{Symbol,Float64}()     # first call of each routine (JIT compile), excluded
timed(s, f) = function ()
    t = time_ns(); f(); dtt = (time_ns() - t) * 1e-9
    if !haskey(FIRST, s); FIRST[s] = dtt
    else TIMES[s] = get(TIMES, s, 0.0) + dtt; CALLS[s] = get(CALLS, s, 0) + 1; end
    nothing
end
function use!(m)
    for s in keys(FORT); TD.ROUTINES[s] = timed(s, FORT[s]); end
    TD.install_hooks!()
    if m == "jfull"                          # ALL residuals + Jacobian (AD) + solver in Julia
        TD.ROUTINES[:matgen] = timed(:matgen, () -> TR.matgen!())
        TD.ROUTINES[:matslv] = timed(:matslv, JS.matslv!)
        TD.install_hooks!(eos = JP.eos!, opac = JP.opac!)
    end
    if m == "jstep"                          # Phase 5b: the whole time step in Julia (only I/O and init remain Fortran)
        TD.ROUTINES[:matgen] = timed(:matgen, () -> TX.matgen!())
        TD.ROUTINES[:matslv] = timed(:matslv, BS.Penta.matslv!)
        TD.install_hooks!(eos = JP.eos!, opac = JP.opac!)
        for (s, f) in ((:peq, TC.peq!), (:reset, TC.reset!), (:update, TC.update!), (:maxdel, TC.maxdel!),
                       (:timstp, TC.timstp!), (:totnrg, TC.totnrg!))
            TD.ROUTINES[s] = timed(s, f)
        end
        TC.EDDFAC[] = TC.eddfac!
    end
    if m == "jcore"                          # Phase 5: jfast + Julia reset/update/maxdel/timstp/totnrg (only peq, I/O in Fortran)
        TD.ROUTINES[:matgen] = timed(:matgen, () -> TX.matgen!())
        TD.ROUTINES[:matslv] = timed(:matslv, BS.Penta.matslv!)
        TD.install_hooks!(eos = JP.eos!, opac = JP.opac!)
        for (s, f) in ((:reset, TC.reset!), (:update, TC.update!), (:maxdel, TC.maxdel!), (:timstp, TC.timstp!), (:totnrg, TC.totnrg!))
            TD.ROUTINES[s] = timed(s, f)
        end
    end
    if m == "jfast"                          # Phase 4: allocation-free assembly + static block-penta solver
        TD.ROUTINES[:matgen] = timed(:matgen, () -> TX.matgen!())
        TD.ROUTINES[:matslv] = timed(:matslv, BS.Penta.matslv!)
        TD.install_hooks!(eos = JP.eos!, opac = JP.opac!)
    end
    if m in ("jsolve", "jall", "jcontin", "jeos")
        TD.ROUTINES[:matslv] = timed(:matslv, JS.matslv!)
    end
    if m in ("jall", "jcontin")
        TD.ROUTINES[:matgen] = timed(:matgen, () -> (FORT[:matgen](); JS.contin_rows!(); m == "jall" && JS.mass_rows!()))
    end
    m in ("jall", "jeos") && TD.install_hooks!(eos = JP.eos!, opac = JP.opac!)
end

use!(mode)                                # (hooks active during init for jall)
# bookkeeping fixes: ON by default; TITAN_BITWISE=1 for bitwise TITAN reproduction
TF.init(nn; rundir = dir)
# optional COMMON overrides after init, e.g. TITAN_SET="jstepe=2000"
for kv in split(get(ENV, "TITAN_SET", ""), ","; keepempty = false)
    k, val = split(kv, "="); T = eltype(TF.var(Symbol(k)))
    TF.setscalar!(Symbol(k), T <: Integer ? T(parse(Int, val)) : parse(Float64, val))
end
empty!(TIMES); empty!(CALLS)

# ---------------------------------------------------------- tier (b)
relmax(a, b) = maximum(abs(x - y) / max(abs(x), abs(y), floatmin()) for (x, y) in zip(a, b); init = 0.0)
function stepstate()
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    d = Dict{Symbol,Any}(n => TF.var(n)[a:b+1] for n in (:rn, :xmen, :un, :frn))
    for n in (:dn, :tn, :ern); d[n] = TF.var(n)[a:b]; end
    d[:iter] = Int(TF.getscalar(:iter)); d[:dtime] = TF.getscalar(:dtime); d
end
res_b = []
function tier_b(n)
    s0 = TD.snapshot()
    use!("fortran"); TD.step!(); A = stepstate()
    TD.restore!(s0)
    use!(mode); TD.step!(); B = stepstate()
    TD.restore!(s0)
    lg = maximum(relmax(A[k], B[k]) for k in (:rn, :xmen, :dn, :tn, :ern))
    nu = maximum(abs.(A[:un] .- B[:un])) / max(maximum(abs.(A[:un])), floatmin())
    nf = maximum(abs.(A[:frn] .- B[:frn])) / max(maximum(abs.(A[:frn])), floatmin())
    conv = TF.getscalar(:conv); ok = lg <= conv && nu <= conv && nf <= conv
    say(@sprintf("  tier(b) step %4d: iter F/J = %d/%d  max rel (r,m,ρ,T,E) %.2e  |Δu| %.2e  |ΔF| %.2e  Δdt/dt %.2e  %s",
                 n, A[:iter], B[:iter], lg, nu, nf, abs(B[:dtime] / A[:dtime] - 1), ok ? "PASS" : "FAIL"))
    push!(res_b, ok)
end

# ---------------------------------------------------------- run
snaps = Metrics.Snap[]; checkpoints = Set([1, 5, 25, 100, 200, 300, 400, 500, 600, 700])
t0 = time(); n = 0
while true
    global n += 1
    dotierb && n in checkpoints && tier_b(n)
    n == 3 && (global t3 = time(); empty!(TIMES); empty!(CALLS))   # steady state: steps 3..n
    tm_ = @timed TD.titan_step!(); st = tm_.value
    n <= 2 && say(@sprintf("  step %d wall %.2f s (compile %.2f s, gc %.2f s)", n, tm_.time, tm_.compile_time, tm_.gctime))
    push!(snaps, Metrics.snap(TF))
    st != 0 && break
end
wall = time() - t0
s = snaps[end]
say("$name/$mode: steps = $n, final jstep = $(s.jstep), t = $(s.time), wall = $(round(wall, digits = 2)) s, ",
    @sprintf("%.2f ms/step", 1e3wall / n))
dotierb && say("tier(b): $(count(res_b))/$(length(res_b)) pass")
say("  time by routine, steps 3..n (s, calls, ms/call):")
for (k, t) in sort(collect(TIMES), by = x -> -x[2])
    say(@sprintf("    %-8s %8.3f  %6d  %8.3f", k, t, CALLS[k], 1e3t / CALLS[k]))
end
wall3 = time() - t3
say(@sprintf("    %-8s %8.3f   (archive/writout/writsty/driver)", "other", wall3 - sum(values(TIMES))))
say(@sprintf("  steady state (steps 3..%d, excl. JIT): %.3f s, %.3f ms/step; first-call (compile) times: %s", n, wall3,
             1e3wall3 / (n - 2), join(["$k=$(round(v, digits = 2))s" for (k, v) in FIRST if v > 0.05], ", ")))
mkpath(dir); serialize(joinpath(dir, "snaps.jls"), snaps)
TF.finalize(); close(logio)
