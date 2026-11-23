# Stand-alone harness for the TDC rows of julia/TDC.jl (M3a): implicit backward-Euler stepping of the w rows alone,
# prescribed background, Newton with a ForwardDiff (tridiagonal, 3-colour) Jacobian, plane (mu=0) or spherical (mu=2).
using ForwardDiff, LinearAlgebra, Printf
include(joinpath(@__DIR__, "..", "TDC.jl"))
using .TDC

const LOGDIR = normpath(joinpath(@__DIR__, "..", "..", "docs", "logs", "subconv", "m3a"))
const LOGIO = Ref{Union{IO,Nothing}}(nothing)
open_log(name) = (mkpath(LOGDIR); LOGIO[] = open(joinpath(LOGDIR, name * ".log"), "w"); nothing)
close_log() = (LOGIO[] === nothing || close(LOGIO[]); LOGIO[] = nothing)
function lp(args...)
    println(args...)
    LOGIO[] === nothing || (println(LOGIO[], args...); flush(LOGIO[]))
end
const NFAIL = Ref(0)
function check(name, ok::Bool, detail = "")
    lp(ok ? "PASS " : "FAIL ", name, "  ", detail)
    ok || (NFAIL[] += 1)
    return ok
end
write_csv(name, header, rows) = open(joinpath(LOGDIR, name), "w") do io
    println(io, header); for r in rows; println(io, join(r, ",")); end
end

"""
    make_grid(r; μ=0, rho=1.0, mdot=0.0, dVdt=0.0) -> TDCGeom

Interfaces r (nz+1 increasing). μ = 0: plane (area 1, dvol = Δr); μ = 2: spherical (area 4π r_c², dvol = 4π/3 Δ(r³)).
rho: zone density (scalar or vector); mdot: mass flux at the interfaces (scalar, vector, or function of r), TITAN sign (positive inward);
masses new = old (static grid).
"""
function make_grid(r; μ = 0, rho = 1.0, mdot = 0.0, dVdt = 0.0)
    nz = length(r) - 1
    rc = [0.5 * (r[j] + r[j+1]) for j in 1:nz]
    dvol = μ == 0 ? [r[j+1] - r[j] for j in 1:nz] : [4π / 3 * (r[j+1]^3 - r[j]^3) for j in 1:nz]
    Az = μ == 0 ? ones(nz) : [4π * rc[j]^2 for j in 1:nz]
    ρ = rho isa Number ? fill(float(rho), nz) : collect(float.(rho))
    md = mdot isa Number ? fill(float(mdot), nz + 1) : mdot isa Function ? [float(mdot(x)) for x in r] : collect(float.(mdot))
    dv = dVdt isa Number ? fill(float(dVdt), nz) : collect(float.(dVdt))
    return TDCGeom(nz, collect(float.(r)), Az, ρ .* dvol, ρ .* dvol, md, dv)
end
zone_rho(G::TDCGeom, rho) = rho isa Number ? fill(float(rho), G.nz) : collect(float.(rho))
_prof(x, r) = x isa Number ? fill(float(x), length(r)) : x isa Function ? [float(x(ri)) for ri in r] : collect(float.(x))
"background with a, U, sup, ell as scalars, vectors (nz+1) or functions of r; rho, ellz per zone (scalars/vectors)"
function make_back(G::TDCGeom; a = 1.0, U = 0.0, sup = 1.0, ell = 1.0, rho = 1.0, ellz = ell isa Number ? ell : 1.0)
    rc = [0.5 * (G.r[j] + G.r[j+1]) for j in 1:G.nz]
    return TDCBack(_prof(a, G.r), _prof(U, G.r), _prof(sup, G.r), _prof(ell, G.r), zone_rho(G, rho),
                   ellz isa Function ? [float(ellz(x)) for x in rc] : zone_rho(G, ellz))
end

"tridiagonal Jacobian of f!(φ, w) by ForwardDiff with 3 colours (row i depends on w_{i-1}, w_i, w_{i+1} only)"
function jac_tridiag(f!, w)
    n = length(w); D = ForwardDiff.Dual{Nothing,Float64,1}
    dl = zeros(n - 1); d = zeros(n); du = zeros(n - 1)
    for c in 0:2
        wd = [D(w[i], ForwardDiff.Partials((i % 3 == c ? 1.0 : 0.0,))) for i in 1:n]
        φ = zeros(D, n); f!(φ, wd)
        for i in 1:n, j in max(1, i - 1):min(n, i + 1)
            j % 3 == c || continue
            v = ForwardDiff.partials(φ[i], 1)
            j == i ? (d[i] = v) : j == i - 1 ? (dl[i-1] = v) : (du[i] = v)
        end
    end
    return Tridiagonal(dl, d, du)
end

"""
    (w, info) = newton_w(wo, G, B, dt, par; w0, cen, tol, maxit)

One implicit step of the TDC rows from wo (semismooth Newton on the Fischer–Burmeister rows, backtracking on ‖φ‖₂).
Convergence: ‖φ‖∞ ≤ tol·max(‖w‖∞, ‖wo‖∞) or the step below 1e-15 ‖w‖∞.
"""
function newton_w(wo, G, B, dt, par; w0 = copy(wo), cen = 0.0, tol = 1e-12, maxit = 80)
    f!(φ, w) = tdc_rows!(φ, w, wo, G, B, dt, par; cen = cen)
    n = length(wo); w = copy(w0); φ = zeros(n); f!(φ, w); nls = 0
    for it in 1:maxit
        sc = max(norm(w, Inf), norm(wo, Inf), 1e-300)
        if norm(φ, Inf) <= tol * sc
            return w, (iters = it - 1, ok = true, resid = norm(φ, Inf) / sc, ls = nls)
        end
        J = jac_tridiag(f!, w)
        dw = -(J \ φ)
        λ = 1.0; φt = similar(φ); wt = similar(w); nφ = norm(φ)
        while true
            wt .= w .+ λ .* dw; f!(φt, wt)
            (norm(φt) < (1 - 1e-4λ) * nφ || λ < 1e-10) && break
            λ /= 2; nls += 1
        end
        w .= wt; φ .= φt
        if λ * norm(dw, Inf) <= 1e-13 * max(norm(w, Inf), 1e-300) && norm(φ, Inf) <= 1e-8 * sc
            return w, (iters = it, ok = true, resid = norm(φ, Inf) / sc, ls = nls)
        end
    end
    sc = max(norm(w, Inf), norm(wo, Inf), 1e-300)
    return w, (iters = maxit, ok = false, resid = norm(φ, Inf) / sc, ls = nls)
end

"""
    evolve(w, G, B, dt, nsteps, par; cen, stop, cb) -> (w, hist)

Backtracking-free time loop; stops early when max|Δw| ≤ stop·max|w| (steady). cb(step, t, w, wo, info) is called after each step.
"""
function evolve(w, G, B, dt, nsteps, par; cen = 0.0, stop = 0.0, cb = nothing, tol = 1e-12)
    w = copy(w); t = 0.0; nbad = 0; maxit = 0
    for s in 1:nsteps
        wo = copy(w)
        w, info = newton_w(wo, G, B, dt, par; cen = cen, tol = tol)
        info.ok || (nbad += 1); maxit = max(maxit, info.iters)
        t += dt
        cb === nothing || cb(s, t, w, wo, info)
        if stop > 0 && norm(w - wo, Inf) <= stop * max(norm(w, Inf), 1e-300)
            return w, (steps = s, t = t, nbad = nbad, maxit = maxit, steady = true)
        end
    end
    return w, (steps = nsteps, t = t, nbad = nbad, maxit = maxit, steady = false)
end

"scalar implicit step of the local ODE dw/dt = R_loc(w): fb(w, w − wo − dt R) = 0 (Newton with ForwardDiff derivative)"
function local_step(wo, a, U, sup, ℓ, dt, par; tol = 1e-14, maxit = 100)
    f(w) = fb(w, (w - wo) - dt * tdc_local(w, a, U, sup, ℓ, par))
    w = wo
    for it in 1:maxit
        v = f(w); d = ForwardDiff.derivative(f, w)
        abs(v) <= tol * max(abs(w), abs(wo), 1e-300) && return w
        dw = -v / d; λ = 1.0
        while abs(f(w + λ * dw)) > (1 - 1e-4λ) * abs(v) && λ > 1e-10; λ /= 2; end
        w += λ * dw
        abs(λ * dw) <= 1e-13 * abs(w) && return w          # stagnation: the F noise floor ε dt|R| exceeds tol|w| for dt >> relaxation time
    end
    return w
end
