# Check of TitanResidFast.GASROW (gas-only energy row) on TITAN's radiating-shock benchmarks:
# the same run with TITAN's gas+radiation row and with the gas-only row must agree to truncation error.
#   TITAN_NOLIB=1 julia --project=env gasrow_check.jl SUBC|SUPC [maxsteps]
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TM = TitanMaterials
name = ARGS[1]; maxsteps = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 100000
function run(gasrow)
    TF.load_state!(joinpath(TF.ROOT, "reference_runs", "init", "$name.tstate.gz"))
    PS.install_julia!(BlockSolver); TM.use!(TM.Builtin()); TX.GASROW[] = gasrow; TD.MSGIO[] = devnull
    TF.setscalar!(:jstep, Int32(1))
    t0 = time(); n, why = TD.run_julia!(maxsteps = maxsteps)
    a, b = Int(TF.getscalar(:ngrs)), Int(TF.getscalar(:ngre))
    s = (n = n, why = why, t = TF.getscalar(:timen), sec = time() - t0, r = copy(TF.var(:rn)[a:b+1]),
         T = copy(TF.var(:tn)[a:b]), ρ = copy(TF.var(:dn)[a:b]), E = copy(TF.var(:ern)[a:b]), u = copy(TF.var(:un)[a:b+1]),
         etot = (TF.getscalar(:tes), TF.getscalar(:tew), TF.getscalar(:tel), TF.getscalar(:teq)))
    @printf("%s gasrow=%d: %d steps (%s), t = %.5e in %.1f s\n", name, gasrow, n, why, s.t, s.sec)
    s
end
A = run(false); B = run(true)
rel(x, y) = maximum(abs.(x ./ y .- 1))
if length(A.T) == length(B.T)
    @printf("max rel diff: T %.3e  ρ %.3e  E %.3e  r %.3e | max|Δu|/max|u| %.3e\n", rel(B.T, A.T), rel(B.ρ, A.ρ), rel(B.E, A.E),
            rel(B.r, A.r), maximum(abs.(B.u .- A.u)) / maximum(abs.(A.u)))
    # profile differences excluding the shock zone(s): median of the zone-wise relative T difference
    d = sort(abs.(B.T ./ A.T .- 1)); @printf("T rel diff: median %.3e, 90%% %.3e\n", d[end÷2], d[9end÷10])
end
println("energy ledger A ", A.etot, "\n              B ", B.etot)
