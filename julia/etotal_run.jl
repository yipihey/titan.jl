# Phase 4 item 2: run a benchmark with the Julia stack (TitanResidFast + Penta solver + Julia
# EOS/opacity) using either TITAN's internal-energy row (std) or the conservative total-energy
# row (etot[ηsw]); records per-step discrete total-energy imbalance, Newton iterations, retries,
# and Metrics snaps for tier (c) (lib/p2_runs/NAME_TAG/snaps.jls, readable by tierc.jl TAG).
# usage: julia --project=env etotal_run.jl NAME NN std|etot [ηsw]
using Serialization, Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "Metrics.jl"))
using .TitanDriver, .Metrics
const TD = TitanDriver; const TF = TD.TF; const JP = JuliaPhysics; const TX = TitanResidFast; const BS = BlockSolver
name, nn, form = ARGS[1], parse(Int, ARGS[2]), ARGS[3]
ηsw = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 0.0
tag = (form == "std" ? "std" : (ηsw > 0 ? "etot$(ARGS[4])" : "etot")) * get(ENV, "RUNTAG", "")
TX.ETOTAL[] = form == "etot"; TX.ETASW[] = ηsw
const STEPDT = Ref(0.0); const NIT = Ref(0)
TD.ROUTINES[:matgen] = () -> (NIT[] += 1; TX.matgen!())
TD.ROUTINES[:matslv] = BS.Penta.matslv!
TD.ROUTINES[:timstp] = () -> (STEPDT[] = TF.getscalar(:dtime); TD.fcall(:timstp_))
TD.install_hooks!(eos = JP.eos!, opac = JP.opac!)
dir = joinpath(TF.ROOT, "lib", "p2_runs", "$(name)_$tag"); mkpath(dir)
TF.init(nn; rundir = dir)
for kv in split(get(ENV, "TITAN_SET", ""), ","; keepempty = false)
    k, val = split(kv, "="); T = eltype(TF.var(Symbol(k)))
    TF.setscalar!(Symbol(k), T <: Integer ? T(parse(Int, val)) : parse(Float64, val))
end
snaps = Metrics.Snap[]; L = NamedTuple[]
n = 0; t0 = time()
while true
    global n += 1; NIT[] = 0
    st = TD.titan_step!()
    e = TX.energy_imbalance(STEPDT[])
    nsw = TX.ETASW[] > 0 ? count(TX.dualswitch(Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre)),
                                                Int(TF.getscalar(:lrad)) == 1, TX.ETASW[])) : 0
    push!(L, (; e..., t = TF.getscalar(:timen), dt = STEPDT[], nit = NIT[], nsw, cuts = Int(TF.getscalar(:jtry)) - 1))
    get(ENV, "ETOT_VERBOSE", "0") == "1" && (n <= 10 || n % 25 == 0 || L[end].cuts > 0) && @printf("  step %3d t %.4e dt %.3e nit %d cuts %d imb %.2e\n", n, L[end].t, L[end].dt, L[end].nit, L[end].cuts, L[end].imb)
    push!(snaps, Metrics.snap(TF))
    st != 0 && break
end
wall = time() - t0
TF.finalize()
serialize(joinpath(dir, "snaps.jls"), snaps); serialize(joinpath(dir, "ledger.jls"), L)
cum = cumsum([l.imb for l in L]); W = cumsum([l.fout - l.fin for l in L])
E0 = L[1].Enew - L[1].dE
scale = max(abs(E0), maximum(abs, W))
@printf("%s/%s: steps %d, t_end %.6e, Newton iterations %d (%.2f/step), wall %.1f s\n", name, tag, n, L[end].t,
        sum(l.nit for l in L), sum(l.nit for l in L) / n, wall)
@printf("  discrete total energy: E0 %.6e, scale %.6e; cumulative imbalance/scale %+.3e; max |per-step imbalance|/scale %.3e; max zones on internal-energy row %d\n",
        E0, scale, cum[end] / scale, maximum(abs(l.imb) for l in L) / scale, maximum(l.nsw for l in L))
