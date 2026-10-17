# Temporal order of the time integrators with fixed steps h = T_end/N, on a fixed grid: TITAN's θ-scheme with θ = 1
# (backward Euler) and θ = ½, and Trbdf2. Error at T_end against a Trbdf2 run with N = 4096 (same grid).
#   Z0_TABLE=ext TITAN_NOLIB=1 julia --project=env order_test.jl [compton|wave|radwave]
# compton: the box of compton_box.jl (gas at 1e8 K relaxing in a radiation field; gas-only row, Compton, free-free);
#          error in T of the middle zone.
# wave:    standing sound wave between two walls (ideal gas, 64 zones, amplitude 1e-3 c_s), half a period;
#          error in u (max norm over interfaces, relative to the amplitude).
# radwave: the same with radiation (TITAN grey opacity, optical depth 1 per zone; the outer boundary
#          is irradiated with I = B(T), so the uniform state does not leak, E = aT⁴, p_rad ≈ 1.8 p_gas).
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl")); include(joinpath(@__DIR__, "Trbdf2.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TM = TitanMaterials; const Z0 = Z0Physics
n = 10; L = 3e19; ρ0 = 1e-17; T0 = 1e8; E0 = 10.0
s = Z0.saha(ρ0, T0); tC = 1.5 * (s.nH + s.nHe + s.ne) * Z0.me * Z0.c^2 / (4 * 6.6524616e-25 * s.ne * Z0.c * E0)
Tend = 2tC
test = length(ARGS) > 0 ? Symbol(ARGS[1]) : :compton
nw = 64; Lw = 1e10; ρw = 1e-8; Tw = 1e5; γw = 5 / 3; gmu = 0.6
csw = sqrt(γw * 1.380658e-16 / (1.6605402e-24 * gmu) * Tw); Aw = parse(Float64, get(ENV, "AMP", "1e-3")) * csw
test == :compton || (Tend = Lw / csw)          # half a period of the fundamental standing wave
function setup!()
    if test == :compton
        PS.build!(; rI = collect(range(L, 2L, length = n + 1)), rho = fill(ρ0, n), T = fill(T0, n), u = zeros(n + 1), geom = 0,
                  rad = true, Er = fill(E0, n), material = TM.Z0Gas(), lgrid = 2, xscale = L)
        TX.GASROW[] = true; TX.COMPTON[] = true
    else
        rI = collect(range(Lw, 2Lw, length = nw + 1)); u = [Aw * sin(π * (x - Lw) / Lw) for x in rI]
        rad = test == :radwave
        Ew = nothing                       # radiation: E = aT⁴ (p_rad ≈ 1.8 p_gas at 1e5 K, 1e-8 g/cm³)
        PS.build!(; rI, rho = fill(ρw, nw), T = fill(Tw, nw), u, geom = 0, rad, Er = Ew, gam = γw, gmu, lgrid = 2, xscale = Lw,
                  opac = rad ? (lopac = 3, chif0 = nw / Lw) : nothing, av = (cq1 = 0.0, cq2 = 0.0, cadv = parse(Float64, get(ENV, "CADV", "1.0"))))   # no artificial viscosity (smooth flow); CADV=0: no van Leer slopes
        rad && TF.setscalar!(:ximr, TF.getscalar(:csigr) * Tw^4 / TF.getscalar(:cpi))   # incident I = B(T): no net leak
        TX.GASROW[] = false; TX.COMPTON[] = false
    end
end
err(q, ref) = test == :compton ? abs(q.T / ref.T - 1) : maximum(abs.(q.u .- ref.u)) / Aw
function run(scheme, N)
    setup!()
    PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
    TF.setscalar!(:conv, parse(Float64, get(ENV, "CONV", "1e-5")))     # Newton tolerance (TITAN default 1e-5)
    TD.STEPFN[] = scheme in (:be, :cn) ? TD.step! : Trbdf2.step!
    TF.setscalar!(:thet, scheme == :cn ? 0.5 : 1.0)
    h = Tend / N
    for i in 1:N
        TF.setscalar!(:dtime, h)
        TD.STEPFN[]()
    end
    a = Int(TF.getscalar(:ngrs)); k = a + n ÷ 2
    (t = TF.getscalar(:timen), T = TF.var(:tn)[k], E = TF.var(:ern)[k], u = copy(TF.var(:un)[a:a+nw]))
end
NREF = parse(Int, get(ENV, "NREF", "4096"))
ref = run(:trbdf2, NREF)
@printf("%s reference (Trbdf2, N = %d): t = %.6e (T_end %.6e)\n", test, NREF, ref.t, Tend)
for sch in Symbol.(split(get(ENV, "SCHEMES", "cn,be,trbdf2"), ","))
    prev = NaN
    for N in parse.(Int, split(get(ENV, "NLIST", "8,16,32,64,128"), ","))
        q = run(sch, N); e = err(q, ref)
        @printf("%-7s N = %4d  t = %.6e  error = %.3e  order %.2f\n", sch, N, q.t, e, isnan(prev) ? NaN : log2(prev / e))
        prev = e
    end
end
