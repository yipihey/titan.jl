# Step a saved quasi-star state and report, per step, what limits it: Newton iterations, the variable and zone with
# the largest relative change (timstp), and the local profile there.
#   TITAN_NOLIB=1 julia --project=env stepdiag.jl STATE.tstate.gz NSTEPS [Z0Gas grey=0]
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl")); include(joinpath(@__DIR__, "Trbdf2.jl")); include(joinpath(@__DIR__, "DampedNewton.jl"))
using Printf
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const TM = TitanMaterials; const PS = ProblemSetup
TF.load_state!(ARGS[1]); PS.install_julia!(BlockSolver)
if haskey(ENV, "QS_FREEZE")     # DampedNewton with branch freezing from iteration QS_FREEZE on
    DampedNewton.KFREEZE[] = parse(Int, ENV["QS_FREEZE"]); DampedNewton.DTRIG[] = parse(Float64, get(ENV, "QS_DTRIG", "1e30"))
    TD.STEPFN[] = DampedNewton.step!
end
get(ENV, "QS_TRBDF2", "0") == "1" && (TD.STEPFN[] = Trbdf2.step!)
haskey(ENV, "QS_MLT_OLD") && (TX.MLT_AD[] = false; TX.MLT_HS[] = false)
TX.CONVECTION[] = get(ENV, "QS_CONV", "1") == "1"; TX.GASROW[] = get(ENV, "QS_GASROW", "0") == "1"; TX.COMPTON[] = get(ENV, "QS_COMPTON", "0") == "1"; TM.use!(TM.Z0Gas()); TD.MSGIO[] = stdout
for kv in split(get(ENV, "QS_SET", ""), ","; keepempty = false)
    k, val = split(kv, "="); T = eltype(TF.var(Symbol(k))); TF.setscalar!(Symbol(k), T <: Integer ? T(parse(Int, val)) : parse(Float64, val))
end
names = ("r", "m", "ρ", "u", "T", "E", "F")
a = Int(TF.getscalar(:ngrs)); R = 8.4352e15
v(s) = TF.var(s)::Vector{Float64}
nits = Ref(0); oldupd = TD.ROUTINES[:update]; TD.ROUTINES[:update] = () -> (nits[] += 1; oldupd())
for n in 1:parse(Int, ARGS[2])
    nits[] = 0
    TD.STEPFN[]()
    kx = TF.var(:kx); dx = TF.var(:dx)
    i = argmax(abs.(dx[[1, 2, 3, 5, 6]])); i = (1, 2, 3, 5, 6)[i]; k = a + kx[i] - 1
    @printf("step %d t=%.6e dt=%.3e updates=%d  limiter %s at k=%d (r/R=%.6f) Δ=%.3g | ", n, TF.getscalar(:timen), TF.getscalar(:dtime), nits[], names[i], k, v(:rn)[k] / R, dx[i])
    @printf("ρ %.3e %.3e %.3e  T %.4g %.4g %.4g  u %.3e %.3e %.3e  E %.4e %.4e %.4e\n", v(:dn)[k-1], v(:dn)[k], v(:dn)[k+1], v(:tn)[k-1], v(:tn)[k], v(:tn)[k+1], v(:un)[k-1], v(:un)[k], v(:un)[k+1], v(:ern)[k-1], v(:ern)[k], v(:ern)[k+1])
end
haskey(ENV, "QS_SAVE") && (TF.save_state(ENV["QS_SAVE"]); println("saved ", ENV["QS_SAVE"]))
