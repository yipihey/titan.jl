# Larson–Penston self-similar solution of isothermal collapse (t < 0, before the central singularity), computed here
# as the reference for popiii_m1.jl lp (M1).
#   r = c_s (−t) ξ... written in Shu's (1977) variables with x = r/(c_s t) < 0 for t < 0:
#   ρ = α(x)/(4πG t²),  u = c_s v(x),  m = c_s³ t m(x)/G (m < 0 for t < 0, so −m is the physical mass)
#   [(x−v)² − 1] v' = [α(x−v) − 2/x](x−v),   [(x−v)² − 1] α' = [α − 2(x−v)/x](x−v) α
# Regular centre: v = 2x/3, α = α0 (+ O(x²)). The Larson–Penston solution is the one α0 for which the integral curve
# passes regularly through the sonic line (x−v)² = 1 (numerator and denominator vanish together). Bisection on α0:
# for α0 too small/large the v-numerator vanishes after/before the sonic line is reached.
#   julia --project=env lp_selfsimilar.jl          → α0, the sonic point, and the asymptotes α x² and v at x → −∞
using Printf

rhs(x, v, α) = (D = (x - v)^2 - 1; ((α * (x - v) - 2 / x) * (x - v) / D, (α - 2(x - v) / x) * (x - v) * α / D))

"Integrate from the centre towards negative x until the sonic line or a v-numerator zero; return which came first."
function shoot(α0; h = -1e-5, xend = -10.0)
    x = -1e-4; v = 2x / 3; α = α0
    while x > xend
        D = (x - v)^2 - 1; Nv = (α * (x - v) - 2 / x)
        D > -1e-6 && return (:sonic, x, v, α)
        Nv < 0 && return (:numer, x, v, α)             # Nv changes sign from + (centre: 2/x → −∞ gives Nv > 0)
        k1 = rhs(x, v, α); k2 = rhs(x + h / 2, v + h / 2 * k1[1], α + h / 2 * k1[2])
        k3 = rhs(x + h / 2, v + h / 2 * k2[1], α + h / 2 * k2[2]); k4 = rhs(x + h, v + h * k3[1], α + h * k3[2])
        v += h / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]); α += h / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2]); x += h
    end
    (:none, x, v, α)
end

# bracket: classify α0 on a coarse scan, then bisect
cls(a) = shoot(a)[1]
lo, hi = 1.0, 3.0
println("scan: ", join((@sprintf("%.2f:%s", a, cls(a)) for a in 1.0:0.25:3.0), "  "))
cl0 = cls(lo); cls(hi) == cl0 && error("no bracket")
for _ in 1:50
    global lo, hi
    mid = (lo + hi) / 2; cls(mid) == cl0 ? (lo = mid) : (hi = mid)
end
α0 = (lo + hi) / 2
_, xs, vs, αs = shoot(α0)
@printf("alpha0 = %.6f (bracket width %.1e); sonic point x* = %.5f, v* = %.5f, alpha* = %.5f\n", α0, hi - lo, xs, vs, αs)

# continue across the sonic point: step over it on the straight line through the last two points of the two shots
# that bracket α0 (they leave on either side), then integrate outward from just beyond the sonic line to large |x|
function outer(α0; dx = 0.02, xend = -1e4)
    _, x1, v1, a1 = shoot(α0)
    # local slope from the regular branch: use finite difference of the curve just before the sonic point
    _, x0, v0, a0 = shoot(α0; h = -1e-5, xend = x1 + 0.01)
    sv = (v1 - v0) / (x1 - x0); sa = (a1 - a0) / (x1 - x0)
    x = x1 - dx; v = v1 - dx * sv; α = a1 - dx * sa
    h = -1e-4
    while x > xend
        k1 = rhs(x, v, α); k2 = rhs(x + h / 2, v + h / 2 * k1[1], α + h / 2 * k1[2])
        k3 = rhs(x + h / 2, v + h / 2 * k2[1], α + h / 2 * k2[2]); k4 = rhs(x + h, v + h * k3[1], α + h * k3[2])
        v += h / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]); α += h / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2]); x += h
        h = -1e-4 * max(1.0, abs(x))
        abs(x) in (10.0, 100.0) && nothing
    end
    (x, v, α)
end
for xe in (-10.0, -100.0, -1000.0, -10000.0)
    x, v, α = outer(α0; xend = xe)
    @printf("x = %9.1f: alpha x^2 = %.5f, v = %.5f\n", x, α * x^2, v)
end

# ---- check by the critical-point method: on the sonic line x − v = −1 the regular crossing needs α* = −2/x*
# (both numerators vanish: α(x−v) = 2/x and α = 2(x−v)/x); the slopes (p = v', s = α') there follow from l'Hôpital:
#   p (D_x + D_v p) = G_x + G_v p + G_α s,   s (D_x + D_v p) = H_x + H_v p + H_α s
# with D = (x−v)² − 1, G = [α(x−v) − 2/x](x−v), H = [α − 2(x−v)/x](x−v)α. Eliminate s, solve for p (two roots), step
# off the critical point towards the centre and integrate; the LP solution is the x* for which the centre is regular
# (v/x → 2/3 as x → 0). Partial derivatives by central differences.
D_(x, v, α) = (x - v)^2 - 1
G_(x, v, α) = (α * (x - v) - 2 / x) * (x - v)
H_(x, v, α) = (α - 2(x - v) / x) * (x - v) * α
function grad(f, x, v, α; e = 1e-6)
    ((f(x + e, v, α) - f(x - e, v, α)) / 2e, (f(x, v + e, α) - f(x, v - e, α)) / 2e, (f(x, v, α + e) - f(x, v, α - e)) / 2e)
end
function slopes(xs)
    vs = xs + 1; as = -2 / xs
    Dx, Dv, _ = grad(D_, xs, vs, as); Gx, Gv, Ga = grad(G_, xs, vs, as); Hx, Hv, Ha = grad(H_, xs, vs, as)
    sof(p) = (Hx + Hv * p) / (Dx + Dv * p - Ha)
    f(p) = p * (Dx + Dv * p) - (Gx + Gv * p + Ga * sof(p))
    # roots of f by scanning + bisection
    roots = Float64[]; ps = range(-20, 20, length = 40001)
    for i in 1:length(ps)-1
        a, b = ps[i], ps[i+1]; fa, fb = f(a), f(b)
        (isfinite(fa) && isfinite(fb) && fa * fb < 0) || continue
        for _ in 1:80; m = (a + b) / 2; f(m) * fa > 0 ? (a = m; fa = f(m)) : (b = m); end
        abs(f((a + b) / 2)) < 1e-6 && push!(roots, (a + b) / 2)
    end
    vs, as, [(p, sof(p)) for p in roots]
end
"Integrate from the critical point x* inward (towards x = 0) along slope (p, s); return v/x at x = −1e-3."
function inward(xs, p, s; dx = 1e-4)
    vs = xs + 1; as = -2 / xs
    x = xs + dx; v = vs + dx * p; α = as + dx * s; h = 1e-5
    while x < -1e-3
        k1 = rhs(x, v, α); k2 = rhs(x + h / 2, v + h / 2 * k1[1], α + h / 2 * k1[2])
        k3 = rhs(x + h / 2, v + h / 2 * k2[1], α + h / 2 * k2[2]); k4 = rhs(x + h, v + h * k3[1], α + h * k3[2])
        v += h / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]); α += h / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2]); x += h
        (isfinite(v) && isfinite(α) && abs(α) < 1e6) || return (NaN, NaN)
    end
    (v / x, α)
end
println("\ncritical-point method: x*, slope roots p, and v/x, alpha at x = -1e-3 (regular centre: v/x = 2/3)")
for xs in (-2.0, -2.2, -2.3, -2.3185, -2.33, -2.4, -2.6)
    vs, as, rts = slopes(xs)
    for (p, s) in rts
        vx, a0 = inward(xs, p, s)
        @printf("  x* = %8.4f  p = %9.5f  s = %9.5f   v/x = %10.5f  alpha = %10.5f\n", xs, p, s, vx, a0)
    end
end

# ---- the sonic point is a node: every α0 below the bisection edge reaches the sonic line with both numerators → 0.
# For each α0: where it meets the sonic line, the numerator there, and the asymptotes after stepping across.
println("\nalpha0 family through the node: x*, Nv/|x-v| at the stop, alpha x^2 and v at x = -1e4")
for a in (1.55, 1.60, 1.65, 1.667, 1.68, 1.70, 1.72)
    c, x1, v1, a1 = shoot(a)
    x, v, α = outer(a; xend = -1e4)
    @printf("  alpha0 = %.4f: %-6s x* = %8.5f  Nv = %+.2e  alpha x^2 = %8.5f  v = %9.5f\n", a, c, x1, a1 * (x1 - v1) - 2 / x1, α * x^2, v)
end

# ---- select LP by analyticity at the node: curves arrive with slope v' → p1 = 1/|x*| plus a non-analytic correction
# C|Δx|^(λ−1) in v' that dominates the analytic O(Δx) term near x*. The LP solution has C = 0: find the α0 where
# (v' − p1) at a fixed small distance from the sonic line changes sign. (p1 evaluated at the stopping x, where the
# curve is within D > −δ of the sonic line.)
function nodeprobe(a; δ = 1e-4)
    x = -1e-4; v = 2x / 3; α = a; h = -1e-6
    while true
        D = (x - v)^2 - 1
        D > -δ && break
        k1 = rhs(x, v, α); k2 = rhs(x + h / 2, v + h / 2 * k1[1], α + h / 2 * k1[2])
        k3 = rhs(x + h / 2, v + h / 2 * k2[1], α + h / 2 * k2[2]); k4 = rhs(x + h, v + h * k3[1], α + h * k3[2])
        v += h / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]); α += h / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2]); x += h
        h = D > -0.05 ? -1e-7 : -1e-6
    end
    xs = x - (1 + (x - v)) / (1 - rhs(x, v, α)[1])          # Newton estimate of where x − v = −1
    rhs(x, v, α)[1] - 1 / abs(xs), xs
end
println("\nanalyticity probe: v' − 1/|x*| near the sonic line (δ = 1e-4, 1e-5)")
for a in (1.60, 1.64, 1.66, 1.67, 1.68, 1.70)
    g4, xs = nodeprobe(a; δ = 1e-4); g5, _ = nodeprobe(a; δ = 1e-5)
    @printf("  alpha0 = %.3f  x* = %.5f  dv' (1e-4) = %+.4e  dv' (1e-5) = %+.4e\n", a, xs, g4, g5)
end
for δ in (1e-4, 1e-5, 1e-6)
    lo_, hi_ = 1.64, 1.70
    for _ in 1:40; m = (lo_ + hi_) / 2; nodeprobe(m; δ)[1] > 0 ? (lo_ = m) : (hi_ = m); end
    aLP = (lo_ + hi_) / 2; _, xs = nodeprobe(aLP; δ)
    xo, vo, αo = outer(aLP; xend = -1e5)
    @printf("LP (analytic at the node, δ = %.0e): alpha0 = %.5f, x* = %.5f, alpha x^2 -> %.4f, v -> %.4f (x = %.0f)\n",
            δ, aLP, xs, αo * xo^2, vo, xo)
end
