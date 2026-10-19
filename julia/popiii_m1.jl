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

"""
Energy ledger (per steradian × 4π = total) on the Lagrangian grid: K = Σ m̄_k u_k²/2 over interfaces (m̄ = half the masses
of the two neighbouring zones, as in the momentum row), U = Σ Δm e_g, W = −Σ G m_k m̄_k / r_k (the potential whose time
derivative is the work of the momentum row's gravity term), boundary term p_ext V (constant p_ext). Returns
(K, U, W, pV, E = K + U + W + pV).
"""
function energy()
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    r, u, d, e, dv = TF.var(:rn), TF.var(:un), TF.var(:dn), TF.var(:egn), TF.var(:dvoln); xm = TX.XMINT[]
    dm(k) = (a <= k <= b) ? d[k] * dv[k] : 0.0
    K = 0.0; U = 0.0; W = 0.0
    for k in a:b+1
        mb = 0.5 * (dm(k - 1) + dm(k)); K += 0.5mb * u[k]^2
        r[k] > 0 && (W -= TF.getscalar(:cgrav) * 4π * xm[k] * mb / r[k])
    end
    for k in a:b; U += dm(k) * e[k]; end
    pV = isnan(TX.PEXT[]) ? 0.0 : TX.PEXT[] * r[b+1]^3 / 3
    4π .* (K, U, W, pV, K + U + W + pV)
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
    TC.DT_DYN[] = parse(Float64, get(ENV, "CDYN", "0"))
    out("C_dyn = %g\n", TC.DT_DYN[])
    out("theta = %.3f, stol = %.3g, conv = %.1e\n", TF.getscalar(:thet), TF.getscalar(:stol), TF.getscalar(:conv))
    # β(t) from β + sin β cos β = ωt (Newton); r/r0 = cos²β for every shell
    function βof(τ)
        β = min(π / 2 - 1e-12, cbrt(1.5τ))
        for _ in 1:60; f = β + sin(β) * cos(β) - τ; β -= f / (2cos(β)^2 + 1e-300); β = clamp(β, 0.0, π / 2); end
        β
    end
    s0 = state(); rlim = parse(Float64, get(ENV, "RLIM", "1e-2")); E0 = energy()
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
    E1 = energy()
    out("energy: K %.4e U %.4e W %.4e pV %.4e; E %.6e -> %.6e, |dE|/|W| = %.3e\n", E1[1], E1[2], E1[3], E1[4], E0[5], E1[5],
        abs(E1[5] - E0[5]) / abs(E1[3]))
    s = slast[]
    out("steps %d, t/t_ff = %.6f, R/R0 = %.3e; max |r/r_exact − 1| = %.3e, max |u/u_exact − 1| = %.3e\n",
        nsteps[], s.t / tff, s.r[end] / s0.r[end], worst[], worstu[])
elseif which == "be"
    # Bonnor–Ebert sphere, ξ_max = 3 (stable; critical 6.451), nearly isothermal gas (γ = 1.001), free surface at the
    # pressure that holds the outer half shell. Variant "discrete" (default): radii marched outward from the centre
    # so that TITAN's momentum rows vanish at rest (exact discrete equilibrium; tests stability and round-off).
    # Variant "le": zone densities sampled from the Lane–Emden solution (tests the discretisation of hydrostatics).
    n = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 60; variant = length(ARGS) > 2 ? ARGS[3] : "discrete"
    ξmax = parse(Float64, get(ENV, "XIMAX", "3.0")); gam = 1.001; gmu = 1.0; T0 = 300.0; ρc = 1e-19
    TF.load_state!(joinpath(PS.INIT, "SOD.tstate.gz")); Rg = TF.getscalar(:ck) / (TF.getscalar(:cm0) * gmu)
    a_ = sqrt(Rg * T0 / (4π * Z0.G * ρc))
    # Lane–Emden ψ'' + 2ψ'/ξ = e^−ψ, RK4 from the series ψ = ξ²/6 at small ξ; μ(ξ) = ρc a³ ξ² ψ' per steradian
    function lane_emden(ξs)
        f(ξ, y) = (y[2], exp(-y[1]) - 2y[2] / ξ)
        ξ = 1e-4; y = (ξ^2 / 6, ξ / 3); out_ = Tuple{Float64,Float64}[]
        for ξt in ξs
            while ξ < ξt
                h = min(1e-3, ξt - ξ)
                k1 = f(ξ, y); k2 = f(ξ + h / 2, y .+ h / 2 .* k1); k3 = f(ξ + h / 2, y .+ h / 2 .* k2); k4 = f(ξ + h, y .+ h .* k3)
                y = y .+ h / 6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4); ξ += h
            end
            push!(out_, y)
        end
        out_
    end
    ξI = [ξmax * (j / n)^1.5 for j in 0:n]
    le = lane_emden(ξI[2:end]); μI = [0.0; [ρc * a_^3 * ξI[j+1]^2 * le[j][2] for j in 1:n]]; dμ = diff(μI)
    if variant == "le"
        rI = a_ .* ξI; rho = [dμ[j] / ((rI[j+1]^3 - rI[j]^3) / 3) for j in 1:n]
    else
        # march: zone 1 from the central density; interface k: r_k²(p_k − p_{k−1}) + 4πG μ_k/r_k² · (dμ_{k−1}+dμ_k)/2 = 0
        rI = zeros(n + 1); rho = zeros(n)
        rI[2] = a_ * ξI[2]; rho[1] = dμ[1] / (rI[2]^3 / 3)
        for k in 2:n
            pk = Rg * T0 * rho[k-1] - 4π * Z0.G * μI[k] / rI[k]^4 * 0.5 * (dμ[k-1] + dμ[k])
            rho[k] = pk / (Rg * T0); rI[k+1] = cbrt(rI[k]^3 + 3dμ[k] / rho[k])
        end
    end
    PS.build!(; rI, rho, T = fill(T0, n), u = zeros(n + 1), geom = 2, lgrid = 3, inner = :wall, outer = :free,
              pext = 0.0, selfgrav = true, gam, gmu)
    TX.REGULAR_CENTRE[] = true; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); xm = TX.XMINT[]
    TX.PEXT[] = TF.var(:pgn)[b] - 4π * Z0.G * xm[b+1] / rI[end]^4 * 0.5 * TF.var(:dn)[b] * TF.var(:dvoln)[b]
    cs = sqrt(Rg * T0); tsc = rI[end] / cs
    haskey(ENV, "THETA") && TF.setscalar!(:thet, parse(Float64, ENV["THETA"]))
    TC.DT_DYN[] = parse(Float64, get(ENV, "CDYN", "0.01"))
    TF.setscalar!(:dtime, 1e-3tsc); TF.setscalar!(:tmax, 10tsc)
    s0 = state(); wd = Ref(0.0); wu = Ref(0.0); nan = Ref(false)
    out("Bonnor-Ebert xi_max = %.3f, n = %d, variant %s, M = %.4g Msun, R = %.4e cm, rho_c/rho_edge = %.3f, p_ext = %.4e, theta = %.2f, C_dyn = %g\n",
        ξmax, n, variant, 4π * xm[b+1] / Msun, rI[end], rho[1] / rho[end], TX.PEXT[], TF.getscalar(:thet), TC.DT_DYN[])
    n_, why = TD.run_julia!(maxsteps = 100000, output = function (k)
        s = state(); (any(isnan, s.r) || any(isnan, s.d)) && (nan[] = true)
        wd[] = max(wd[], maximum(abs.(s.d ./ s0.d .- 1))); wu[] = max(wu[], maximum(abs.(s.u)) / cs)
        k % 200 == 0 && out("step %6d t/t_sc = %.3f  max|Δρ/ρ| = %.3e  max|u|/c_s = %.3e\n", k, s.t / tsc, wd[], wu[])
    end)
    s = state()
    out("%d steps (%s) to t = %.2f t_sc; max |Δρ/ρ| = %.3e, max |u|/c_s = %.3e, NaN %s\n", n_, why, s.t / tsc, wd[], wu[], nan[])
    out("%s (criterion: < 1e-3 both)\n", !nan[] && wd[] < 1e-3 && wu[] < 1e-3 ? "PASS" : "FAIL")
elseif which == "lp"
    # Larson–Penston: discrete-equilibrium Bonnor–Ebert sphere at ξ_max = 6.5 (just above critical) on a geometric
    # mass grid (n shells, Δm ratio q outward), T lowered by 10 % at t = 0. Nearly isothermal (γ = 1.001). Follow the
    # centre over DECADES decades of ρ_c; print the inner profiles of ρ r² 4πG/c_s² and u/c_s against r.
    n = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 200; q = parse(Float64, get(ENV, "Q", "1.1"))
    decades = parse(Float64, get(ENV, "DECADES", "8")); ξmax = 6.5; gam = 1.001; gmu = 1.0; T0 = 300.0; ρc = 1e-19
    TF.load_state!(joinpath(PS.INIT, "SOD.tstate.gz")); Rg = TF.getscalar(:ck) / (TF.getscalar(:cm0) * gmu)
    a_ = sqrt(Rg * T0 / (4π * Z0.G * ρc))
    # total mass of the BE sphere: μ = ρc a³ ξ²ψ'(ξ_max), ψ' from RK4 as in "be"
    ξ = 1e-4; y = (ξ^2 / 6, ξ / 3); fle(ξ, y) = (y[2], exp(-y[1]) - 2y[2] / ξ)
    while ξ < ξmax
        h = min(1e-3, ξmax - ξ); k1 = fle(ξ, y); k2 = fle(ξ + h / 2, y .+ h / 2 .* k1); k3 = fle(ξ + h / 2, y .+ h / 2 .* k2)
        k4 = fle(ξ + h, y .+ h .* k3); global y = y .+ h / 6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4); global ξ += h
    end
    μtot = ρc * a_^3 * ξmax^2 * y[2]
    μI = shellmasses(μtot, n, q); dμ = diff(μI)
    rI = zeros(n + 1); rho = zeros(n); rI[2] = cbrt(3dμ[1] / ρc); rho[1] = ρc
    for k in 2:n
        pk = Rg * T0 * rho[k-1] - 4π * Z0.G * μI[k] / rI[k]^4 * 0.5 * (dμ[k-1] + dμ[k])
        pk > 0 || error("march: p ≤ 0 at k = $k")
        rho[k] = pk / (Rg * T0); rI[k+1] = cbrt(rI[k]^3 + 3dμ[k] / rho[k])
    end
    PS.build!(; rI, rho, T = fill(T0, n), u = zeros(n + 1), geom = 2, lgrid = 3, inner = :wall, outer = :free,
              pext = 0.0, selfgrav = true, gam, gmu)
    TX.REGULAR_CENTRE[] = true; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre)); xm = TX.XMINT[]
    TX.PEXT[] = TF.var(:pgn)[b] - 4π * Z0.G * xm[b+1] / rI[end]^4 * 0.5 * TF.var(:dn)[b] * TF.var(:dvoln)[b]
    # start the collapse: T × 0.9 everywhere, external pressure × 0.9
    fT = parse(Float64, get(ENV, "FT", "0.9")); tn = TF.var(:tn)
    for k in a-1:b+2; tn[k] *= fT; end
    TX.PEXT[] *= fT          # the surface pressure drops with T: no compression wave driven from outside
    for (o, nn) in TC.PEQ_OLD1; copyto!(view(TF.var(o), a-1:b+2), view(TF.var(nn), a-1:b+2)); end
    TD.EOS[](); for (o, nn) in TC.PEQ_OLD2_H; copyto!(view(TF.var(o), a-1:b+1), view(TF.var(nn), a-1:b+1)); end
    cs = sqrt(Rg * T0 * fT); tff0 = sqrt(3π / (32Z0.G * ρc))
    haskey(ENV, "THETA") && TF.setscalar!(:thet, parse(Float64, ENV["THETA"]))
    TC.DT_DYN[] = parse(Float64, get(ENV, "CDYN", "0.01"))
    TF.setscalar!(:dtime, 1e-4tff0); TF.setscalar!(:tmax, 100tff0)
    out("Larson-Penston: BE xi_max = %.2f, M = %.4g Msun, n = %d, q = %.3f, dm_1/M = %.2e, T x %.2f, theta = %.2f, C_dyn = %g\n",
        ξmax, 4π * μtot / Msun, n, q, dμ[1] / μtot, fT, TF.getscalar(:thet), TC.DT_DYN[])
    function profile(s, tag)
        out("--- %s: t/t_ff0 = %.6f, rho_c/rho_c0 = %.4e\n", tag, s.t / tff0, s.d[1] / ρc)
        out("%12s %12s %14s %10s\n", "r [cm]", "m/M", "4piG rho r^2/cs^2", "u/c_s")
        for k in 1:4:n
            rz = 0.5 * (s.r[k] + s.r[k+1])
            out("%12.4e %12.4e %14.4f %10.4f\n", rz, μI[k+1] / μtot, 4π * Z0.G * s.d[k] * rz^2 / cs^2, s.u[k+1] / cs)
        end
    end
    next_dec = Ref(1.0); s0 = state(); slast = Ref(s0); E0 = energy()
    try
        TD.run_julia!(maxsteps = parse(Int, get(ENV, "MAXSTEPS", "200000")), output = function (k)
            s = state(); slast[] = s; dc = log10(s.d[1] / ρc)
            dbg = parse(Int, get(ENV, "DEBUG", "0"))
            if haskey(ENV, "DUMPRHO") && k in parse.(Int, split(ENV["DUMPRHO"], ","))
                kmx = argmax(s.d[1:n-20])
                out("  rho k=%d t/t_ff0 = %.5f: max rho/rho_1 = %.6f at k = %d; rho_k/rho_1 for k = 1:20:%d:", k, s.t / tff0, s.d[kmx] / s.d[1], kmx, n)
                for j in 1:20:n; out(" %.5f", s.d[j] / s.d[1]); end; out("\n")
            end
            if dbg > 0 && k % dbg == 0
                ku = argmin(s.u); du = diff(s.u); kj = argmax(abs.(diff(du))[1:n-20])
                out("  dbg %6d t/t_ff0 = %.6f dt = %.3e  min u/cs = %.4f at k = %d (m/M = %.2e)  max |Δ²u|/cs = %.3e at k = %d\n",
                    k, s.t / tff0, TF.getscalar(:dtime), s.u[ku] / cs, ku, μI[ku] / μtot, abs(diff(du)[kj]) / cs, kj + 1)
                if abs(diff(du)[kj]) / cs > parse(Float64, get(ENV, "DUMPAT", "Inf"))
                    q_ = TF.var(:qf); a0 = Int(TF.getscalar(:ngrs))
                    for j in max(1, kj - 8):min(n, kj + 8)
                        out("    k %4d r %.6e u/cs %9.5f rho %.6e T %.4f qf %.3e dm/dm- %.4f\n", j, s.r[j], s.u[j] / cs, s.d[j], s.T[j],
                            q_[a0+j-1], j > 1 ? dμ[j] / dμ[j-1] : 0.0)
                    end
                    TF.setscalar!(:tmax, s.t)
                end
            end
            if dc >= next_dec[]
                out("step %7d t/t_ff0 = %.8f  log rho_c/rho_c0 = %.2f  dt = %.3e s  T_c = %.2f\n", k, s.t / tff0, dc, TF.getscalar(:dtime), s.T[1])
                Ek = energy(); out("        energy K %.4e U %.4e W %.4e pV %.4e  (E - E0)/|W| = %.3e\n", Ek[1], Ek[2], Ek[3], Ek[4], (Ek[5] - E0[5]) / abs(Ek[3]))
                (round(Int, next_dec[]) % 2 == 0) && profile(s, @sprintf("%.0f decades", next_dec[]))
                next_dec[] += 1
            end
            dc >= decades && TF.setscalar!(:tmax, s.t)
        end)
    catch e
        out("stopped: %s\n", sprint(showerror, e)[1:min(end, 300)])
    end
    profile(slast[], "final")
end
close(logio)
