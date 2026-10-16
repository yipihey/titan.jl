# Test ladder for the accretion problem: strong inflow onto a wall, from pure hydrodynamics up to the quasi-star.
#   TITAN_NOLIB=1 julia --project=env infall_tests.jl NAME [key=value …]
# Tests
#   noh1d      planar Noh problem: cold gas at −v0 onto a reflecting wall (ideal gas, γ = 5/3), exact solution known
#   noh_sph    spherical Noh (inner wall at small radius), exact solution known
# Options (key=value): n (zones), lgrid (2 Eulerian, 1 adaptive), gridtau, mach, ql (AV length in zones), cq1, cq2,
#   tfrac (end time in units of the time the shock needs to cross half the domain), every (snapshot cadence)
using Printf, Serialization
include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
include(joinpath(@__DIR__, "TitanCore.jl"))
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
include(joinpath(@__DIR__, "ProblemSetup.jl"))
using .TitanDriver
const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const BS = BlockSolver
const TC = TitanCore; const TM = TitanMaterials; const PS = ProblemSetup; const Z0 = Z0Physics

name = ARGS[1]
opt = Dict(Symbol(split(a, "=")[1]) => split(a, "=")[2] for a in ARGS[2:end])
o(k, d) = haskey(opt, k) ? (d isa Integer ? parse(Int, opt[k]) : d isa AbstractString ? opt[k] : parse(Float64, opt[k])) : d
outdir = joinpath(@__DIR__, "..", "lib", "tests"); mkpath(outdir)

"Exact Noh solution (γ, geometry index ν = 0 planar, 2 spherical) at time t for inflow speed v0, density ρ0."
function noh_exact(x, t, v0, ρ0, γ, ν)
    us = (γ - 1) * v0 / 2                       # shock speed
    xs = us * t
    if x < xs
        ρ = ρ0 * ((γ + 1) / (γ - 1))^(ν + 1); u = 0.0; p = (γ - 1) / 2 * ρ * v0^2 * 0 + ρ * v0^2 * (γ - 1) / 2   # p2 = ρ2 (γ−1) v0²/2
        return (ρ = ρ, u = 0.0, p = 0.5 * (γ - 1) * ρ * v0^2, xs = xs)
    else
        return (ρ = ρ0 * (1 + v0 * t / max(x, 1e-300))^ν, u = -v0, p = 0.0, xs = xs)
    end
end

function run_noh(ν)
    n = o(:n, 100); lgrid = o(:lgrid, 2); mach = o(:mach, 200.0)
    L = 1e12; v0 = 1e7; ρ0 = 1e-10; γ = 5 / 3; gmu = 0.6
    Rg = 1.380658e-16 / (1.6605402e-24 * gmu)
    T0 = (v0 / mach)^2 / (γ * Rg)
    rmin = ν == 0 ? L : 0.05L          # TITAN's unknowns are ln r: the wall cannot sit at r = 0 (planar: shift by L)
    rI = collect(range(rmin, rmin + L, length = n + 1))
    dx = L / n
    ql = o(:ql, 2.0)                         # AV length in zones
    av = ν == 0 ? (ql0 = ql * dx, ql1 = 0.0, cq1 = o(:cq1, 0.0), cq2 = o(:cq2, 1.0)) :
                  (ql0 = ql * dx, ql1 = 0.0, cq1 = o(:cq1, 0.0), cq2 = o(:cq2, 1.0))
    u = [k == 1 ? 0.0 : -v0 for k in 1:n+1]
    ρ = [ν == 0 ? ρ0 : ρ0 for _ in 1:n]
    PS.build!(; rI, rho = ρ, T = fill(T0, n), u, geom = ν, rad = false, inner = :wall, outer = :inflow,
              inflow = (ρ0 * (ν == 0 ? 1.0 : 1.0), T0, -v0), gam = γ, gmu, lgrid, gridtau = o(:gridtau, 1e-20),
              lady = (0, 2, 0, 0, 2, 2), xscale = L, av)
    # spherical: the inflow density at r_out grows with time as (1 + v0 t/r)²; for this test keep the
    # boundary state fixed (valid until the rarefaction from the boundary reaches the shock); compare inside
    PS.install_julia!(BS)
    TF.setscalar!(:dtime, 1e-3 * dx / v0); TF.setscalar!(:thet, o(:thet, 1.0))
    tend = o(:tfrac, 1.0) * (0.5L) / ((γ - 1) * v0 / 2)
    TF.setscalar!(:tmax, tend)
    TD.MSGIO[] = get(ENV, "QS_MSG", "0") == "1" ? stdout : devnull
    snaps = Any[PS.snap()]; every = o(:every, 10)
    t0 = time()
    nst, why = try TD.run_julia!(output = k -> (k % every == 0 && push!(snaps, PS.snap())))
    catch e
        println("STOPPED: ", first(sprint(showerror, e), 200)); (0, :failed)
    end
    push!(snaps, PS.snap())
    s = snaps[end]
    # compare with the exact solution at zone centres
    xc = 0.5 .* (s.r[1:end-1] .+ s.r[2:end])
    ex = [noh_exact(x - (ν == 0 ? rmin : 0.0), s.t, v0, ρ0, γ, ν) for x in xc]
    xs = ex[1].xs
    # error excluding 3 zones around the shock and the wall-heating zones (first 3)
    mask = [abs(x - (ν == 0 ? rmin : 0.0) - xs) > 3dx && x > rmin + 3dx for x in xc]
    eρ = maximum(abs.(s.rho[mask] ./ [e.ρ for e in ex][mask] .- 1))
    # numerical shock position: half-jump crossing of ρ
    ρ2 = ex[1].ρ; ρhalf = 0.5 * (ρ2 + ρ0 * (1 + v0 * s.t / xs)^ν)
    ish = findlast(i -> s.rho[i] > ρhalf, eachindex(s.rho))
    xs_num = ish === nothing ? NaN : xc[ish] - (ν == 0 ? rmin : 0.0)
    wall = s.rho[1] / ρ2
    ov = maximum(s.rho) / ρ2
    umax = maximum(s.u[2:end]) / v0     # positive velocities = backflow/oscillation
    @printf("%s: %s after %d steps, t = %.4g (target %.4g) in %.1f s | shock %.4g vs %.4g (Δ = %.2f zones) | max|ρ/ρ_exact−1| away from shock & wall = %.3g | wall ρ/ρ2 = %.3f | overshoot ρmax/ρ2 = %.3f | max u/v0 = %.3g\n",
            name, why, nst, s.t, tend, time() - t0, xs_num, xs, (xs_num - xs) / dx, eρ, wall, ov, umax)
    serialize(joinpath(outdir, "$(name)_$(join(["$k=$v" for (k, v) in opt], "_")).jls"), (snaps = snaps, v0, ρ0, γ, ν, L))
end

name == "noh1d" && run_noh(0)
name == "noh_sph" && run_noh(2)
