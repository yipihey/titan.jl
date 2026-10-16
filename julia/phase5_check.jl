# Phase 5: in-flight tier-(a) check of the TitanCore translations on a full benchmark run.
# usage: julia --project=env phase5_check.jl NAME NN [routines]      e.g. SUPC 9 reset,update,maxdel,timstp,totnrg
# The run itself uses the full Julia physics stack (jfast: TitanResidFast + Penta solver + Julia eos/opac).
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "PortHarness.jl"))
using .TitanDriver
const TD = TitanDriver; const TF = TD.TF; const JP = JuliaPhysics; const TX = TitanResidFast
const BS = BlockSolver; const TC = TitanCore; const PH = PortHarness

name, nn = ARGS[1], parse(Int, ARGS[2])
which = length(ARGS) >= 3 ? Symbol.(split(ARGS[3], ",")) : [:reset, :update, :maxdel, :timstp, :totnrg]
TD.ROUTINES[:matgen] = () -> TX.matgen!()
TD.ROUTINES[:matslv] = BS.Penta.matslv!
TD.install_hooks!(eos = JP.eos!, opac = JP.opac!)
jl = Dict(:reset => TC.reset!, :update => TC.update!, :maxdel => TC.maxdel!, :timstp => TC.timstp!, :totnrg => TC.totnrg!)
for s in which; PH.wrap!(s, jl[s]); end
get(ENV, "FOLLOW_JULIA", "0") == "1" && (PH.FOLLOW_JULIA[] = true)
get(ENV, "PARITY", "0") == "1" && (TC.PARITY[] = true)

dir = joinpath(TF.ROOT, "lib", "p5_runs", "$(name)_check"); mkpath(dir)
TF.init(nn; rundir = dir)
for kv in split(get(ENV, "TITAN_SET", ""), ","; keepempty = false)
    k, val = split(kv, "="); T = eltype(TF.var(Symbol(k)))
    TF.setscalar!(Symbol(k), T <: Integer ? T(parse(Int, val)) : parse(Float64, val))
end
t0 = time()
n = TD.run!()
@printf("%s: %d steps, t = %.6e, wall %.1f s, FOLLOW_JULIA = %s\n", name, n, TF.getscalar(:timen), time() - t0, PH.FOLLOW_JULIA[])
PH.report()
TF.finalize()
