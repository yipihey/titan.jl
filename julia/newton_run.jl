# Phase 4 item 5: TITAN's dt/20 retry vs damped/line-search Newton (DampedNewton), Julia stack
# (TitanResidFast + Penta solver + Julia EOS/opacity), with optional one-time 1e-14 relative
# T perturbations (seeds) to measure the step-count spread.
# usage: julia --project=env newton_run.jl NAME NN MODES SEEDS     e.g. SOD 3 titan,damped 0,1,2
#   seed 0 = unperturbed; its snaps go to lib/p2_runs/NAME_{jfast|damped}/ for tierc.jl
using Serialization, Printf, Random
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "DampedNewton.jl")); include(joinpath(@__DIR__, "Metrics.jl"))
using .TitanDriver, .Metrics
const TD = TitanDriver; const TF = TD.TF; const JP = JuliaPhysics; const TX = TitanResidFast
const BS = BlockSolver; const DN = DampedNewton
name, nn = ARGS[1], parse(Int, ARGS[2]); modes = split(ARGS[3], ","); seeds = parse.(Int, split(ARGS[4], ","))
for kv in split(get(ENV, "DN_SET", ""), ","; keepempty = false)      # e.g. DN_SET="ALPHA=0.3,TFAC=0.5"
    k, v = split(kv, "="); getfield(DN, Symbol(k))[] = parse(Float64, v)
end
const NMG = Ref(0)
TD.ROUTINES[:matgen] = () -> (NMG[] += 1; TX.matgen!())
TD.ROUTINES[:matslv] = BS.Penta.matslv!
TD.install_hooks!(eos = JP.eos!, opac = JP.opac!)
out = open(joinpath(TF.ROOT, "lib", "p4_runs", "newton_$(name).log"), "a")
say(x...) = (println(x...); println(out, x...); flush(out))
say("# $(name): ALPHA=$(DN.ALPHA[]) DTRIG=$(DN.DTRIG[]) LMIN=$(DN.LMIN[]) TFAC=$(DN.TFAC[])  TITAN_BITWISE=$(get(ENV, "TITAN_BITWISE", "0"))")
for mode in modes, seed in seeds
    TD.STEPFN[] = mode == "damped" ? DN.step! : TD.step!
    DN.resetstats!(); NMG[] = 0
    tag = mode == "damped" ? "damped" : "jfast"
    dir = joinpath(TF.ROOT, "lib", "p2_runs", seed == 0 ? "$(name)_$(tag)_nr" : "$(name)_$(tag)_s$seed"); mkpath(dir)
    TF.init(nn; rundir = dir)
    for kv in split(get(ENV, "TITAN_SET", ""), ","; keepempty = false)
        k, val = split(kv, "="); T = eltype(TF.var(Symbol(k)))
        TF.setscalar!(Symbol(k), T <: Integer ? T(parse(Int, val)) : parse(Float64, val))
    end
    if seed > 0
        rng = MersenneTwister(seed); tn = TF.var(:tn)
        a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
        for k in a:b; tn[k] *= 1 + 1e-14 * randn(rng); end
    end
    snaps = Metrics.Snap[]; n = 0; cuts = 0; t0 = time()
    while true
        n += 1
        st = TD.titan_step!()
        mode == "titan" && (cuts += Int(TF.getscalar(:jtry)) - 1)
        seed == 0 && push!(snaps, Metrics.snap(TF))
        st != 0 && break
    end
    wall = time() - t0; tend = TF.getscalar(:timen)
    TF.finalize()
    seed == 0 && serialize(joinpath(dir, "snaps.jls"), snaps)
    mode == "damped" && (cuts = DN.STATS[:dtcuts])
    say(@sprintf("%-6s %-6s seed %d: steps %4d  t_end %.6f  dt-cut retries %3d  Newton iterations %5d  damped iterations %4d  backtracks %4d  wall %.1f s",
                 name, mode, seed, n, tend, cuts, NMG[], DN.STATS[:damped], DN.STATS[:backtracks], wall))
end
close(out)
