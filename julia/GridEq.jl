"""
GridEq: one global, scale-free equidistribution rule for the radial mesh (`GRIDEQ`, popiii_m2c.jl).

Every zone covers the same ∫ M d ln r with

    M = sqrt(1 + G²),  G² = Σ_q w_q (∂ ln q/∂ ln r)²,  q ∈ {ρ, P, T, E_r, m (interior mass), τ (clamped to [τlo, τhi])},

G² capped smoothly at gmax² (G² gmax²/(G² + gmax²): no zone finer than √(1 + gmax²) times the baseline, so a jump is spread over a
bounded number of zones), and neighbouring zones within the factor (1 + α)/α (the bound of Dorfi–Drury's spatial smoothing). All
terms are logarithmic derivatives: the rule has no scale and no threshold. A power law ρ ∝ r^−n gives a constant M (a logarithmic
mesh that does not move). The implicit form of the same rule is TITAN's grid equation (TitanResidFast.GRIDEQ); this module gives the
explicit target (the starting mesh and the prescribed mode). Plain vectors, no dependence on the rest of the code.
"""
module GridEq

"Cumulative zone count on the fine cells (q = max(λ f0, ·) with the growth limit e^lg per zone), λ such that Σ q = K."
function equi(f0::AbstractVector, K; lg = log(2.5 / 1.5))
    nf = length(f0); q = similar(f0, Float64)
    function fill!(lam)
        @inbounds for i in 1:nf; q[i] = lam * f0[i]; end
        for _ in 1:2
            @inbounds for i in 2:nf; q[i] = max(q[i], q[i-1] * exp(-lg * q[i-1])); end
            @inbounds for i in nf-1:-1:1; q[i] = max(q[i], q[i+1] * exp(-lg * q[i+1])); end
        end
        sum(q)
    end
    lo, hi = 1e-8, 1e8
    for _ in 1:70; lam = sqrt(lo * hi); fill!(lam) > K ? (hi = lam) : (lo = lam); end
    tot = fill!(lo); c = Vector{Float64}(undef, nf + 1); c[1] = 0.0
    @inbounds for i in 1:nf; c[i+1] = c[i] + q[i] * K / tot; end
    c[end] = K; c
end
function positions(xe, c, K)
    p = Vector{Float64}(undef, K + 1); p[1] = xe[1]; p[end] = xe[end]; i = 1
    @inbounds for z in 1:K-1
        while c[i+1] < z && i < length(c) - 1; i += 1; end
        p[z+1] = xe[i] + (xe[i+1] - xe[i]) * (z - c[i]) / max(c[i+1] - c[i], 1e-300)
    end
    p
end
"Linear interpolation in x (flat beyond the ends)."
function lin(xs, ys, x)
    x <= xs[1] && return ys[1]; x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x); ys[i] + (ys[i+1] - ys[i]) * (x - xs[i]) / (xs[i+1] - xs[i])
end
"""
    density(r, m, ρ, P, T, E, τ; w, gmax, τlo, τhi) -> (x, M)

M at the zone centres x = ln r_c (zones 2…n: the central sphere is left out) of the state: interfaces r (r[1] = 0), m (interior
mass), τ (optical depth from the outer edge); zones ρ, P, T, E. Derivatives centred in ln r between zone centres. X (optional, species ×
zones 1…n): further logarithmic fields, e.g. ln(x + x_floor) of the composition, with the common weight wx.
"""
function density(r, m, ρ, P, T, E, τ; w = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0), gmax = 100.0, τlo = 0.1, τhi = 1e4, X = nothing, wx = 0.0)
    n = length(ρ); x = [log(sqrt(r[j] * r[j+1])) for j in 2:n]
    lmc = [log(0.5 * (m[j] + m[j+1])) for j in 2:n]; ltc = [log(clamp(sqrt(max(τ[j], 1e-300) * max(τ[j+1], 1e-300)), τlo, τhi)) for j in 2:n]
    fs = (log.(ρ[2:n]), log.(P[2:n]), log.(T[2:n]), log.(max.(E[2:n], 1e-300)), lmc, ltc)
    M = similar(x)
    for i in eachindex(x)
        il, ih = max(i - 1, 1), min(i + 1, length(x))
        dl = i > il ? (x[i] - x[il]) : 1.0; dh = ih > i ? (x[ih] - x[i]) : 1.0
        w2 = (i > il && ih > i) ? 0.5 : 1.0
        G2 = 0.0
        for (wq, f) in zip(w, fs); G2 += wq * w2 * (((f[i] - f[il]) / dl)^2 + ((f[ih] - f[i]) / dh)^2); end
        (X === nothing || wx == 0) || (G2 += wx * w2 * sum(((X[q, i+1] - X[q, il+1]) / dl)^2 + ((X[q, ih+1] - X[q, i+1]) / dh)^2 for q in axes(X, 1)))   # (X: zones 1…n, x: zones 2…n)
        G2 = G2 * gmax^2 / (G2 + gmax^2); M[i] = sqrt(1 + G2)
    end
    (; x, M)
end
"""
    target(r, m, ρ, P, T, E, τ, N; w, gmax, alph, τlo, τhi, nf) -> N + 1 radii from r[2] (the first interface) to r[end]

The mesh of N zones between the fixed first interface and the outer edge that equidistributes ∫ M d ln r, with neighbour ratio
≤ (1 + α)/α.
"""
function target(r, m, ρ, P, T, E, τ, N; w = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0), gmax = 100.0, alph = 1.5, τlo = 0.1, τhi = 1e4, nf = 6000, X = nothing, wx = 0.0)
    D = density(r, m, ρ, P, T, E, τ; w, gmax, τlo, τhi, X, wx)
    xa, xb = log(r[2]), log(r[end]); dx = (xb - xa) / nf; xe = [xa + (i - 1) * dx for i in 1:nf+1]
    f0 = [lin(D.x, D.M, 0.5 * (xe[i] + xe[i+1])) * dx for i in 1:nf]
    c = equi(f0, N; lg = log((1 + alph) / alph))
    exp.(positions(xe, c, N))
end

end # module
