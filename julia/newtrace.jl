# Trace the Newton iterations of a quasi-star step: per iteration dmax, the variable and zone of the largest
# correction, and the local state. Steps NSKIP steps quietly first.
#   TITAN_NOLIB=1 Z0_TABLE=ext QS_SET=... julia --project=env newtrace.jl STATE.tstate.gz NSKIP NTRACE
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl")); include(joinpath(@__DIR__, "Trbdf2.jl"))
using Printf
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const TM = TitanMaterials; const PS = ProblemSetup
TF.load_state!(ARGS[1]); PS.install_julia!(BlockSolver)
get(ENV, "QS_TRBDF2", "0") == "1" && (TD.STEPFN[] = Trbdf2.step!)
TX.CONVECTION[] = get(ENV, "QS_CONV", "1") == "1"; TX.GASROW[] = get(ENV, "QS_GASROW", "0") == "1"; TX.COMPTON[] = get(ENV, "QS_COMPTON", "0") == "1"; TM.use!(TM.Z0Gas()); TD.MSGIO[] = stdout
for kv in split(get(ENV, "QS_SET", ""), ","; keepempty = false)
    k, val = split(kv, "="); T = eltype(TF.var(Symbol(k))); TF.setscalar!(Symbol(k), T <: Integer ? T(parse(Int, val)) : parse(Float64, val))
end
names = ("r", "m", "ρ", "u", "T", "E", "F")
a = Int(TF.getscalar(:ngrs)); R = 8.4352e15
v(s) = TF.var(s)::Vector{Float64}
TRACE = Ref(false)
oldmd = TD.ROUTINES[:maxdel]
TD.ROUTINES[:maxdel] = () -> begin
    oldmd()
    if TRACE[]
        kx = TF.var(:kx); dx = TF.var(:dx)
        s = join([@sprintf("%s:%+.2e@%d", names[i], dx[i], a + kx[i] - 1) for i in 1:7], " ")
        @printf("   it %2d dt=%.3e dmax=%.3e  %s\n", Int(TF.getscalar(:iter)), TF.getscalar(:dtime), TF.getscalar(:dmax), s)
        if get(ENV, "SCAN", "0") == "1" && Int(TF.getscalar(:iter)) in (4, 5)
            c = TX.context(); neqn = c.neqn; nb = c.nb; k0 = c.ngrs - 2
            rhs = TF.var(:rhs); δ = zeros(neqn * nb)
            for kb in 1:nb, i in 1:neqn; δ[(kb-1)*neqn+i] = rhs[i, k0+kb-1]; end
            R = zeros(neqn, nb); ws = TX.getws(Float64, c.mgr)
            it = Int(TF.getscalar(:it)); kz = 153 - k0 + 1
            for λ in -1.0:0.25:2.0
                TX.residual!(R, ws, c, λ .* δ)
                @printf("      λ=%5.2f  R_T[152..155] = %s  R_E[153]=%.4e\n", λ, join([@sprintf("%+.6e", R[it, kz+j]) for j in -1:2], " "), R[Int(TF.getscalar(:ie)), kz])
            end
        end
    end
end
nskip = parse(Int, ARGS[2]); ntr = parse(Int, ARGS[3])
for n in 1:nskip+ntr
    TRACE[] = n > nskip
    TD.STEPFN[]()
    @printf("step %d t=%.6e dt=%.3e\n", n, TF.getscalar(:timen), TF.getscalar(:dtime))
end
