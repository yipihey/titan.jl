# Check of TitanResidFast.COMPTON: hot, fully ionized gas in a strong, nearly uniform radiation field, optically
# thick to electron scattering (radiation leaks out over ~1e11 s), so gas and radiation exchange energy only by
# Compton scattering and (2–6 % of it) free-free emission. Without free-free and with E fixed:
#   T − T_r = (T0 − T_r) exp(−t/t_C),  t_C = 1.5 n m_e c² / (4 σ_T n_e c E).
# Needs Z0_TABLE=ext (and the log R extrapolation of Z0Physics.EXTRAP_LOWR: this gas lies below the table).
#   TITAN_NOLIB=1 julia --project=env compton_box.jl
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TM = TitanMaterials; const Z0 = Z0Physics
n = 10; L = 3e19; ρ0 = 1e-17; T0 = 1e8; E0 = 10.0
rI = collect(range(L, 2L, length = n + 1))
PS.build!(; rI, rho = fill(ρ0, n), T = fill(T0, n), u = zeros(n + 1), geom = 0, rad = true, Er = fill(E0, n),
          material = TM.Z0Gas(), lgrid = 2, xscale = L)
PS.install_julia!(BlockSolver); TX.GASROW[] = true; TX.COMPTON[] = get(ENV, "NOCOMPTON", "0") != "1"; TD.MSGIO[] = devnull
s = Z0.saha(ρ0, T0); ntot = s.nH + s.nHe + s.ne
tC = 1.5 * ntot * Z0.me * Z0.c^2 / (4 * 6.6524616e-25 * s.ne * Z0.c * E0)
Tr = (E0 / Z0.arad)^0.25
TF.setscalar!(:dtime, 1e-3tC); TF.setscalar!(:tmax, 5tC)
haskey(ENV, "STOL") && TF.setscalar!(:stol, parse(Float64, ENV["STOL"]))   # time-step control (default 0.1)
a = Int(TF.getscalar(:ngrs)); k = a + n ÷ 2
etot() = sum((TF.var(:dn)[j] * Z0.gas_e(TF.var(:dn)[j], TF.var(:tn)[j]) + TF.var(:ern)[j]) for j in a:a+n-1)
E_0 = etot()
@printf("t_C = %.4e s, T_r = %.1f K\n", tC, Tr)
worst = Ref(0.0); Tref = Ref(T0); told = Ref(0.0)
# reference: ρ c_v dT/dt = −Λ_C − κ_P ρ (4σT⁴ − cE), integrated with RK4 sub-steps, using the zone's own E(t)
# (E grows by ~2 % as the gas cools); c_v = 1.5 n k/ρ (fully ionized). Free-free adds 2–6 % to the Compton rate.
tab = Z0.table(); σT = 6.6524616e-25
rate(T, E) = -(4σT * s.ne * Z0.c * E * Z0.kB * (T - (E / Z0.arad)^0.25) / (Z0.me * Z0.c^2) +
               Z0.kapP(tab, ρ0, T) * ρ0 * (4Z0.σSB * T^4 - Z0.c * E)) / (1.5 * ntot * Z0.kB)
TD.run_julia!(output = function (st)
    t = TF.getscalar(:timen); T = TF.var(:tn)[k]; E = TF.var(:ern)[k]
    nsub = 50; h = (t - told[]) / nsub; y = Tref[]
    for _ in 1:nsub
        k1 = rate(y, E); k2 = rate(y + h / 2 * k1, E); k3 = rate(y + h / 2 * k2, E); k4 = rate(y + h * k3, E)
        y += h / 6 * (k1 + 2k2 + 2k3 + k4)
    end
    Tref[] = y; told[] = t; Tex = y
    worst[] = max(worst[], abs(T / Tex - 1))
    st % 20 == 0 && @printf("t/t_C = %.3f  T = %.5e  exact %.5e  rel err %.2e  E = %.6f\n", t / tC, T, Tex, T / Tex - 1, TF.var(:ern)[k])
end)
@printf("max rel error of T(t) in the middle zone: %.3e; total (gas + radiation) energy change %.3e (relative)\n",
        worst[], etot() / E_0 - 1)
