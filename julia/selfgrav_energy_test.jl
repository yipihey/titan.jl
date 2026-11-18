# Self-gravitating energy test on the fixed Eulerian mesh: a cold uniform sphere (M = 100 Msun, rho0 = 1e-19) falls in; ETOTMON's total-energy residual G
# (TitanResidFast.etot_monitor, which books gravity as the mass flux times the potential drop) against the rows that were solved.
#   TITAN_NOLIB=1 [ETOTAL=0] julia --project=env selfgrav_energy_test.jl [n] [fraction of the free-fall time]
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const Z0 = Z0Physics; const TC = TitanCore
n = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 100; fr = length(ARGS) > 1 ? parse(Float64, ARGS[2]) : 0.5
M = 100 * 1.98847e33; ρ0 = 1e-19; T0 = 300.0
R = cbrt(3M / (4π * ρ0)); rI = collect(range(0.02R, R, length = n + 1))
PS.build!(; rI, rho = fill(ρ0, n), T = fill(T0, n), u = zeros(n + 1), geom = 2, lgrid = 2, inner = :wall, outer = :wall, selfgrav = true)
PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
etot = get(ENV, "ETOTAL", "")
etot == "0" && (TX.ECOMPAT[] = true; TX.GASROW[] = true)
tff = sqrt(3π / (32Z0.G * ρ0)); TF.setscalar!(:dtime, 1e-4tff); TF.setscalar!(:tmax, fr * tff)
TC.DT_DYN[] = parse(Float64, get(ENV, "CDYN", "0.02"))      # dynamical-time step cap (as popiii_m1)
mon = Any[]; TD.STAGE_HOOK[] = (kind, dt) -> (push!(mon, TX.etot_monitor(dt)); nothing)
nst, why = TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "20000")))
Gs = sort(abs.([m.G for m in mon]))
@printf("selfgrav n=%d ETOTAL=%s: %d steps (%s), t/t_ff = %.4f; ETOTMON %d solves: median |G| %.3e, max |G| %.3e, sum imb/sum thr %.3e\n", n, etot == "" ? "unset" : etot, nst, why,
        TF.getscalar(:timen) / tff, length(mon), Gs[(length(Gs)+1)÷2], Gs[end], sum(m.imb for m in mon) / sum(m.thr for m in mon))
