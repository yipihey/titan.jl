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
    # ECOMPAT=1: the compatible internal-energy rows (TX.ECOMPAT: first law, discrete adjoint of the momentum row, docs/POPIII_ENERGY.md §23); unset: TITAN's own rows
    get(ENV, "ECOMPAT", "0") == "1" && (TX.ECOMPAT[] = true; TX.GASROW[] = true)
    snaps = Metrics.Snap[]; tstep = Float64[]
    # ETOTMON=1 (diagnostic, default off; no influence on the solution): TX.etot_monitor after every converged implicit solve, summary below
    mon = Any[]
    get(ENV, "ETOTMON", "0") != "0" && (TD.STAGE_HOOK[] = (kind, dt) -> (push!(mon, merge((; dt, kind), TX.etot_monitor(dt))); nothing))
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
    if !isempty(mon) && haskey(ENV, "ETOTMON_OUT")      # one row per solve: the monitor's numbers
        open(ENV["ETOTMON_OUT"], "w") do io
            println(io, "solve,kind,dt,imb,thr,G,gnet,zmax,vmax,zshare,vshare,dE,fin,fout")
            for (i, m) in enumerate(mon); @printf(io, "%d,%s,%.6e,%.6e,%.6e,%.6e,%.6e,%d,%.6e,%d,%.6e,%.6e,%.6e,%.6e\n", i, m.kind, m.dt, m.imb, m.thr, m.G, m.gnet, m.zmax, m.vmax, m.zshare, m.vshare, m.dE, m.fin, m.fout); end
        end
    end
    if !isempty(mon)
        Gs = abs.([m.G for m in mon]); vm = [m.vmax for m in mon]; im = argmax(vm)
        @printf("ETOTMON %s: %d solves; G = (total-energy residual)/(energy throughput): median |G| %.3e, mean |G| %.3e, max |G| %.3e; sum imb / sum thr %.3e; |gnet| max %.3e; worst zone/own throughput %.3e (zone %d, solve %d)\n",
                name, length(mon), sort(Gs)[(length(Gs)+1)÷2], sum(Gs) / length(Gs), maximum(Gs), sum(m.imb for m in mon) / sum(m.thr for m in mon), maximum(abs(m.gnet) for m in mon), vm[im], mon[im].zmax, im)
    end
else
    error("mode must be capture or pure")
end
