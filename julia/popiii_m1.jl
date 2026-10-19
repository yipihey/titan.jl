# Population III project, milestone M1: regular centre, self-gravity, Lagrangian mass grid, free outer boundary.
#   TITAN_NOLIB=1 julia --project=env popiii_m1.jl static      uniform sphere at rest, no gravity, p_ext = p: stays static
#   TITAN_NOLIB=1 julia --project=env popiii_m1.jl dust [n]    pressure-free collapse of a uniform sphere vs the parametric
#                                                              solution r = r0 cos²β, β + sin β cos β = t (8πGρ0/3)^½
# Hydrodynamics only (neqn = 5), TITAN's built-in ideal gas. Logs: lib/popiii/m1_*.log
using Printf
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl"))
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TC = TitanCore
const Z0 = Z0Physics
const Msun = 1.98847e33

"""
Interface masses of n shells, geometric towards the centre: Δm_{j+1}/Δm_j = q (j = 1 innermost), total M.
Returns the interior masses m_j (n+1 values, m_1 = 0).
"""
function shellmasses(M, n, q)
    w = [q^(j - 1) for j in 1:n]; dm = M .* w ./ sum(w)
    [0.0; cumsum(dm)]
end

"Uniform sphere (ρ0, T0) on the mass grid m; regular centre, Lagrangian, free outer surface at pressure pext."
function uniform!(; M, ρ0, T0, n, q, gam = 5 / 3, gmu = 1.0, selfgrav = true, pext = 0.0)
    m = shellmasses(M, n, q)
    rI = @. cbrt(3m / (4π * ρ0)); rI[1] = 0.0
    PS.build!(; rI, rho = fill(ρ0, n), T = fill(T0, n), u = zeros(n + 1), geom = 2, lgrid = 3,
              inner = :wall, outer = :free, pext, selfgrav, gam, gmu)
    TX.REGULAR_CENTRE[] = true
    TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})          # per steradian, fixed (Lagrangian)
    PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
    m, rI
end

state() = (a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre));
           (; t = TF.getscalar(:timen), r = TF.var(:rn)[a:b+1], u = TF.var(:un)[a:b+1], d = TF.var(:dn)[a:b],
              T = TF.var(:tn)[a:b]))

which = isempty(ARGS) ? "static" : ARGS[1]
logio = open(joinpath(TF.ROOT, "lib", "popiii", "m1_$(which)$(get(ENV, "TAG", "")).log"), "w")
out(f, args...) = (s = Printf.format(Printf.Format(f), args...); print(s); print(logio, s); flush(logio))

if which == "static"
    # no gravity, uniform p: every residual vanishes at rest; anything that moves is a centre/boundary defect
    n = 40; M = 100Msun; ρ0 = 1e-19; T0 = 300.0
    uniform!(; M, ρ0, T0, n, q = 1.3, selfgrav = false)
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    TX.PEXT[] = TF.var(:pgn)[b]                                # outside pressure = gas pressure
    cs = sqrt(5 / 3 * TF.var(:pgn)[b] / ρ0); R = TF.var(:rn)[b+1]
    TF.setscalar!(:dtime, 1e-3R / cs); TF.setscalar!(:tmax, 10R / cs)
    s0 = state(); worst = Ref((0.0, 0.0, 0.0)); nan = Ref(false)
    n_, why = TD.run_julia!(output = function (k)
        s = state()
        any(isnan, s.r) || any(isnan, s.u) || any(isnan, s.d) ? (nan[] = true) : nothing
        worst[] = max.(worst[], (maximum(abs.(s.r[2:end] ./ s0.r[2:end] .- 1)), maximum(abs.(s.u)) / cs,
                                 maximum(abs.(s.d ./ s0.d .- 1))))
    end)
    s = state()
    out("static uniform sphere, n = %d, %d steps (%s) to t = %.3g R/c_s; NaN: %s\n", n, n_, why, s.t * cs / R, nan[])
    out("max |Δr/r| = %.2e, max |u|/c_s = %.2e, max |Δρ/ρ| = %.2e, r(centre) = %g\n", worst[]..., s.r[1])
    ok = !nan[] && maximum(worst[]) < 1e-10 && s.r[1] == 0.0
    out("%s\n", ok ? "PASS" : "FAIL")
elseif which == "dust"
    n = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 40
    M = 100Msun; ρ0 = 1e-19; T0 = parse(Float64, get(ENV, "T0", "1e-3"))
    uniform!(; M, ρ0, T0, n, q = 1.2, selfgrav = true, pext = 0.0)
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    tff = sqrt(3π / (32Z0.G * ρ0)); ω = sqrt(8π * Z0.G * ρ0 / 3)
    TF.setscalar!(:dtime, 1e-4tff); TF.setscalar!(:tmax, tff)
    haskey(ENV, "STOL") && TF.setscalar!(:stol, parse(Float64, ENV["STOL"]))
    haskey(ENV, "THETA") && TF.setscalar!(:thet, parse(Float64, ENV["THETA"]))
    haskey(ENV, "CONV") && TF.setscalar!(:conv, parse(Float64, ENV["CONV"]))
    out("theta = %.3f, stol = %.3g, conv = %.1e\n", TF.getscalar(:thet), TF.getscalar(:stol), TF.getscalar(:conv))
    # β(t) from β + sin β cos β = ωt (Newton); r/r0 = cos²β for every shell
    function βof(τ)
        β = min(π / 2 - 1e-12, cbrt(1.5τ))
        for _ in 1:60; f = β + sin(β) * cos(β) - τ; β -= f / (2cos(β)^2 + 1e-300); β = clamp(β, 0.0, π / 2); end
        β
    end
    s0 = state(); rlim = parse(Float64, get(ENV, "RLIM", "1e-2"))
    worst = Ref(0.0); worstu = Ref(0.0); slast = Ref(s0); nsteps = Ref(0)
    out("pressure-free collapse, M = 100 Msun, rho0 = 1e-19, n = %d, q = 1.2, t_ff = %.6e s\n", n, tff)
    out("%6s %12s %12s %12s %12s %12s\n", "step", "t/t_ff", "R/R0", "exact", "max|Δr|/r", "max|Δu|/u")
    try
        TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "20000")), output = function (k)
            s = state(); nsteps[] = k; slast[] = s
            β = βof(ω * s.t); x = cos(β)^2
            # exact: r = r0 x, u = −(2Gm/r0)^½ (1/x − 1)^½ (energy conservation of the shell)
            rex = s0.r[2:end] .* x
            uex = @. -sqrt(2Z0.G * 4π * TX.XMINT[][a+1:b+1] / s0.r[2:end]) * sqrt(1 / x - 1)
            er = maximum(abs.(s.r[2:end] ./ rex .- 1)); eu = maximum(abs.(s.u[2:end] ./ uex .- 1))
            worst[] = max(worst[], er); worstu[] = max(worstu[], eu)
            (k % 20 == 0 || x < rlim) && out("%6d %12.6f %12.4e %12.4e %12.2e %12.2e\n", k, s.t / tff, s.r[end] / s0.r[end], x, er, eu)
            x < rlim && TF.setscalar!(:tmax, s.t)            # stop when the sphere has shrunk by 1/rlim
        end)
    catch e
        out("stopped: %s\n", sprint(showerror, e)[1:min(end, 200)])
    end
    s = slast[]
    out("steps %d, t/t_ff = %.6f, R/R0 = %.3e; max |r/r_exact − 1| = %.3e, max |u/u_exact − 1| = %.3e\n",
        nsteps[], s.t / tff, s.r[end] / s0.r[end], worst[], worstu[])
end
close(logio)
