"""
GridMap: where the interfaces of a radial mesh with a fixed number of zones should be, given the enclosed mass m(r) and the mass of
the star M★. Every zone covers the same arc length in three coordinates,

    X  = ln(r + r_c) / hr          radius; r_c is the radius that encloses fc M★ (the uniform centre is not resolved in ln r),
    Q1 = ln(m + fc M★) / h1        enclosed mass: the interior of the star, dense shells, the envelope once it outweighs the star,
    Q2 = asinh((m − M★)/(ε M★)) / h2   mass counted from the surface, in both directions: the outer layers of the star and the infall,

and, optionally (`tau` given, GMTAU), the grey optical depth from the outer edge,

    Q3 = ln(clamp(τ, τlo, τhi)) / h3   the photosphere and the layers just below and above it (τ from τhi = 10 to τlo = 0.1),

so radius decides the spacing where the enclosed mass does not change and mass decides it where it does, and the place where one
hands over to the other follows M★(t). Two limits: no zone narrower than `dmin` in ln r, and neighbouring zones differ by at most
the factor `g`. The three widths hr, h1, h2 are scaled by one common factor λ so that the zones add up to N.

Used by the driver (`GRIDMAP`, popiii_m2c.jl) for the starting mesh at the hand-over and, every step, for the target the mesh
relaxes to. No dependence on the rest of the code; units cancel (only ratios of r and of m enter).
"""
module GridMap

"""
    anchor(m, u; frac = 0.2) -> M★

Mass inside the surface of the star: the fastest infall is just outside the accretion shock; going inwards from there, the surface is
where the infall has dropped to `frac` of that speed (linear in u between two interfaces, so M★ changes continuously in time). Motions
inside the star slower than `frac` × the infall do not move it. `m`, `u`: enclosed mass and velocity at the interfaces.
"""
function anchor(m, u; frac = 0.2)
    ku = argmin(u); u[ku] < 0 || return m[end]
    thr = frac * u[ku]; j = ku
    while j > 1 && u[j-1] < thr; j -= 1; end
    j == 1 && return m[1]
    m[j-1] + (thr - u[j-1]) / (u[j] - u[j-1]) * (m[j] - m[j-1])
end

"""
    target(r, m, Ms; N, hr, h1, h2, eps, fc, dmin, g, cr, rlo, nf, lam0) -> (; r, lam, nfloor)

Interfaces (N + 1 radii from `rlo` to r[end]) of the mesh described in the module text for the state `r`, `m` (interfaces with r > 0,
m > 0, both increasing) and star mass `Ms`. ln m is taken linear in ln r between the given interfaces and as r³ inside the first.
`cr = 0` switches the centre term off (X = ln r). `lam0`: the λ of the previous call (narrows the search). `nfloor`: zones at the
floor `dmin`.
"""
function target(r, m, Ms; N::Int, hr = 0.3, h1 = 0.3, h2 = 0.4, eps = 1e-5, fc = 1e-3, dmin = 0.02, g = 1.2, cr = 1.0, rlo = r[1],
                nf::Int = 12000, lam0 = NaN, tau = nothing, h3 = 0.4, tlo = 0.1, thi = 10.0, hw = 0.0, wlo = 0.0, whi = 0.0)
    x0 = log.(r); lm0 = log.(m); ns = length(x0)
    xa = log(rlo); xb = x0[end]; dx = (xb - xa) / (nf - 1)
    lnm = Vector{Float64}(undef, nf); j = 1
    @inbounds for i in 1:nf
        x = xa + (i - 1) * dx
        if x <= x0[1]
            lnm[i] = lm0[1] + 3 * (x - x0[1])
        else
            while j < ns - 1 && x > x0[j+1]; j += 1; end
            lnm[i] = lm0[j] + (lm0[j+1] - lm0[j]) * (x - x0[j]) / (x0[j+1] - x0[j])
        end
    end
    # radius that encloses fc Ms
    lrc = -Inf
    if cr > 0
        t = log(fc * Ms); i = searchsortedfirst(lnm, t)
        lrc = i <= 1 ? xa : i > nf ? xb : xa + (i - 2) * dx + dx * (t - lnm[i-1]) / (lnm[i] - lnm[i-1])
        lrc += log(cr)
    end
    rc = exp(lrc)
    # optical depth on the fine grid (`tau` at the interfaces r, decreasing outwards): ln τ linear in r inside a zone (an atmosphere of
    # constant scale height), clamped to [tlo, thi]; outside that range Q3 is constant and adds nothing
    lt = nothing
    if tau !== nothing
        lt = Vector{Float64}(undef, nf); j = 1; ltl, lth = log(tlo), log(thi)
        @inbounds for i in 1:nf
            x = xa + (i - 1) * dx; rr = exp(x)
            if x <= x0[1]
                t = tau[1]
            else
                while j < ns - 1 && x > x0[j+1]; j += 1; end
                f = (rr - r[j]) / (r[j+1] - r[j])
                t = tau[j] > 0 && tau[j+1] > 0 ? tau[j] * (tau[j+1] / tau[j])^f : tau[j] + (tau[j+1] - tau[j]) * f
            end
            lt[i] = clamp(log(max(t, 1e-300)), ltl, lth)
        end
    end
    Q3(i) = lt === nothing ? 0.0 : lt[i] / h3
    X(i) = (x = xa + (i - 1) * dx; cr > 0 ? log(exp(x) + rc) : x)
    Q1(i) = log(exp(lnm[i]) + fc * Ms); Q2(i) = asinh((exp(lnm[i]) - Ms) / (eps * Ms))
    # arc length per cell of the fine grid from the differences of the three coordinates: a jump of Q2 narrower than a cell is counted in full
    n0 = Vector{Float64}(undef, nf - 1); Xp = X(1); Q1p = Q1(1); Q2p = Q2(1); Q3p = Q3(1)
    @inbounds for i in 1:nf-1
        Xn = X(i + 1); Q1n = Q1(i + 1); Q2n = Q2(i + 1); Q3n = Q3(i + 1)
        n0[i] = lt === nothing ? sqrt(((Xn - Xp) / hr)^2 + ((Q1n - Q1p) / h1)^2 + ((Q2n - Q2p) / h2)^2) / dx :
                                 sqrt(((Xn - Xp) / hr)^2 + ((Q1n - Q1p) / h1)^2 + ((Q2n - Q2p) / h2)^2 + (Q3n - Q3p)^2) / dx
        if hw > 0                            # GMWIN: extra radius weight ln r / hw for wlo < r < whi (the infall above the star)
            xm = xa + (i - 0.5) * dx
            wlo < exp(xm) < whi && (n0[i] = sqrt(n0[i]^2 + (1 / hw)^2))
        end
        Xp = Xn; Q1p = Q1n; Q2p = Q2n; Q3p = Q3n
    end
    nmax = 1 / dmin; lg = log(g); q = similar(n0)
    function fill!(lam)                      # zones per unit ln r: floor, then the growth limit |d ln(zone width)/d(zone)| <= ln g
        @inbounds for i in 1:nf-1; q[i] = min(nmax, lam * n0[i]); end
        for _ in 1:2
            @inbounds for i in 2:nf-1; q[i] = max(q[i], q[i-1] * exp(-lg * q[i-1] * dx)); end
            @inbounds for i in nf-2:-1:1; q[i] = max(q[i], q[i+1] * exp(-lg * q[i+1] * dx)); end
        end
        s = 0.0
        @inbounds for i in 1:nf-1; s += q[i]; end
        s * dx
    end
    lo, hi = 1e-4, 1e4
    if isfinite(lam0) && lam0 > 0
        a_, b_ = lam0 / 1.1, lam0 * 1.1
        fill!(a_) < N <= fill!(b_) && ((lo, hi) = (a_, b_))
    end
    for _ in 1:(hi / lo < 2 ? 18 : 48)
        lam = sqrt(lo * hi); fill!(lam) > N ? (hi = lam) : (lo = lam)
    end
    lam = sqrt(lo * hi); tot = fill!(lam)
    xi = Vector{Float64}(undef, N + 1); xi[1] = xa; xi[N+1] = xb; i = 1; c0 = 0.0; c1 = q[1] * dx
    @inbounds for z in 1:N-1
        t = z * tot / N
        while c1 < t && i < nf - 1; i += 1; c0 = c1; c1 += q[i] * dx; end
        xi[z+1] = xa + (i - 1) * dx + dx * (t - c0) / (c1 - c0)
    end
    nfloor = count(z -> xi[z+1] - xi[z] < 1.05 * dmin, 1:N)
    (; r = exp.(xi), lam, nfloor)
end

end # module
