# Phase 6: run TITAN problems with no Fortran at run time.
#   julia --project=env phase6_run.jl capture NAME NN    Fortran set-up (start*.f) once → reference_runs/init/NAME.jls
#   TITAN_NOLIB=1 julia --project=env phase6_run.jl pure NAME    load that state, run every routine in Julia
# The pure run never opens libtitan.so (TITAN_NOLIB=1): state lives in Julia arrays (TitanFortran :julia backend).
using Serialization, Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Metrics.jl"))
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
using .TitanDriver, .Metrics
const TD = TitanDriver; const TF = TD.TF; const JP = JuliaPhysics; const TX = TitanResidFast
const BS = BlockSolver; const TC = TitanCore; const TM = TitanMaterials

initfile(name) = joinpath(TF.ROOT, "reference_runs", "init", "$name.tstate.gz")
mode, name = ARGS[1], ARGS[2]

if mode == "capture"
    nn = parse(Int, ARGS[3])
    mkpath(dirname(initfile(name)))
    TF.init(nn; rundir = joinpath(TF.ROOT, "lib", "p6_runs", "$(name)_capture"))
    TF.save_state(initfile(name))
    TF.finalize()
    println("captured initial state of $name (deck titan.in", lpad(nn, 2, '0'), ") → ", initfile(name))
elseif mode == "pure"
    TF.load_state!(initfile(name))
    for kv in split(get(ENV, "TITAN_SET", ""), ","; keepempty = false)
        k, val = split(kv, "="); T = eltype(TF.var(Symbol(k)))
        TF.setscalar!(Symbol(k), T <: Integer ? T(parse(Int, val)) : parse(Float64, val))
    end
    # every routine of the time step in Julia
    TD.ROUTINES[:matgen] = () -> TX.matgen!()
    TD.ROUTINES[:matslv] = BS.Penta.matslv!
    for (s, f) in ((:peq, TC.peq!), (:reset, TC.reset!), (:update, TC.update!), (:maxdel, TC.maxdel!),
                   (:timstp, TC.timstp!), (:totnrg, TC.totnrg!))
        TD.ROUTINES[s] = f
    end
    TM.use!(TM.Builtin()); TC.EDDFAC[] = TC.eddfac!       # benchmarks: TITAN's own ideal gas + opacities
    TD.MSGIO[] = devnull
    snaps = Metrics.Snap[]; tstep = Float64[]
    t0 = time(); tl = Ref(time())
    n, why = TD.run_julia!(output = function (k)
        push!(snaps, Metrics.snap(TF)); push!(tstep, time() - tl[]); tl[] = time()
    end)
    wall = time() - t0
    dir = joinpath(TF.ROOT, "lib", "p2_runs", "$(name)_pure"); mkpath(dir)
    serialize(joinpath(dir, "snaps.jls"), snaps)
    steady = sum(tstep[3:end]) / max(1, length(tstep) - 2)
    @printf("%s/pure: steps = %d (%s), t = %.10g, wall = %.2f s, steady state %.3f ms/step, libtitan loaded: %s\n",
            name, n, why, TF.getscalar(:timen), wall, 1e3 * steady, TF.HANDLE[] != C_NULL)
else
    error("mode must be capture or pure")
end
