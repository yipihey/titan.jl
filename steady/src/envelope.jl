# Quasi-static quasi-star envelope: hydrostatic, constant luminosity L (from the central BH),
# local MLT convection (Böhm-Vitense / Kippenhahn & Weigert cubic), point-mass BH at the centre.
# Independent variable x = ln P, integrated inward from the outer boundary.

Base.@kwdef struct QSParams
    M::Float64                 # total mass (BH + envelope), g
    MBH::Float64               # BH (core) mass, g
    L::Float64                 # luminosity through the envelope, erg/s
    αMLT::Float64 = 0.8        # matches MESA-QUEST inlist
    ε_in::Float64 = 1e-3       # inner radius in units of R where m must equal MBH
    tab::KapTable
end

LEdd_es(M; X = COMP.X) = 4π * c * G * M / (0.2 * (1 + X))

"MLT temperature gradient and convective velocity."
function mlt(P, T, r, m, L, p::QSParams)
    th = thermo(P, T)
    th === nothing && return (∇ = NaN,)
    ρ = th.ρ
    κ = kap(p.tab, ρ, T)
    g = G * m / r^2
    HP = P / (ρ * g)
    ∇rad = 3κ * L * P / (16π * arad * c * G * m * T^4)
    Γloc = κ * L / (4π * c * G * m)
    if ∇rad <= th.∇ad
        return (∇ = ∇rad, ∇rad = ∇rad, ∇ad = th.∇ad, v = 0.0, cs = th.cs, ρ = ρ, κ = κ,
                HP = HP, Γ = Γloc, β = th.β, Γ1 = th.Γ1)
    end
    l = p.αMLT * HP
    U = 3arad * c * T^3 / (th.cP * ρ^2 * κ * l^2) * sqrt(8HP / (g * th.δ))
    W = ∇rad - th.∇ad
    f(ξ) = (ξ - U)^3 + (8U / 9) * (ξ^2 - U^2 - W)
    a, b = U, sqrt(W + U^2)
    for _ in 1:200
        mid = 0.5(a + b)
        f(mid) > 0 ? (b = mid) : (a = mid)
        (b - a) < 1e-14 * b && break
    end
    ξ = 0.5(a + b)
    ∇ = ξ^2 - U^2 + th.∇ad
    v = sqrt(max(g * th.δ * (ξ - U)^2 * l^2 / (8HP), 0.0))
    return (∇ = ∇, ∇rad = ∇rad, ∇ad = th.∇ad, v = v, cs = th.cs, ρ = ρ, κ = κ, HP = HP,
            Γ = Γloc, β = th.β, Γ1 = th.Γ1)
end

function rhs(x, y, p)
    r, m, lT = y
    P = exp(x); T = exp(lT)
    q = mlt(P, T, r, m, p.L, p)
    isnan(q.∇) && return (NaN, NaN, NaN), q
    drdx = -P * r^2 / (G * m * q.ρ)
    return (drdx, 4π * r^2 * q.ρ * drdx, q.∇), q
end

# Dormand–Prince 5(4) adaptive integrator for the 3-component state (tuples, allocation-free)
const DP_c = (0, 1/5, 3/10, 4/5, 8/9, 1, 1)
const DP_a = ((), (1/5,), (3/40, 9/40), (44/45, -56/15, 32/9),
              (19372/6561, -25360/2187, 64448/6561, -212/729),
              (9017/3168, -355/33, 46732/5247, 49/176, -5103/18656),
              (35/384, 0, 500/1113, 125/192, -2187/6784, 11/84))
const DP_b  = (35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0)
const DP_bs = (5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40)

tadd(y, ks, coeffs, h) = ntuple(i -> y[i] + h * sum(coeffs[j] * ks[j][i] for j in eachindex(coeffs); init = 0.0), 3)

"""
Integrate inward from (P_b, T_b) at r = R, m = M. Stops at r = ε_in R or when m ≤ MBH.
Returns (status, m_at_rin, profile) where status ∈ (:reached_rin, :mass_exhausted, :failed).
"""
function integrate_inward(R, Pb, Tb, p::QSParams; rtol = 1e-7, save = false, maxsteps = 200000)
    x = log(Pb); y = (R, p.M, log(Tb))
    hstep = 0.01
    rin = p.ε_in * R
    prof = save ? Dict{Symbol,Vector{Float64}}(k => Float64[] for k in
           (:r, :m, :P, :T, :ρ, :∇, :∇ad, :∇rad, :v, :cs, :κ, :Γ, :β, :Γ1, :HP)) : nothing
    for step in 1:maxsteps
        k1, q = rhs(x, y, p)
        isnan(k1[1]) && return (:failed, NaN, prof)
        if save
            for (kk, vv) in zip((:r, :m, :P, :T, :ρ, :∇, :∇ad, :∇rad, :v, :cs, :κ, :Γ, :β, :Γ1, :HP),
                                (y[1], y[2], exp(x), exp(y[3]), q.ρ, q.∇, q.∇ad, q.∇rad, q.v, q.cs, q.κ, q.Γ, q.β, q.Γ1, q.HP))
                push!(prof[kk], vv)
            end
        end
        while true
            ks = Vector{NTuple{3,Float64}}(undef, 7); ks[1] = k1
            bad = false
            for s in 2:7
                ys = tadd(y, ks[1:s-1], DP_a[s], hstep)
                ks[s] = rhs(x + DP_c[s] * hstep, ys, p)[1]
                if isnan(ks[s][1]) || ys[1] <= 0
                    bad = true; break
                end
            end
            if bad
                hstep /= 4; hstep < 1e-12 && return (:failed, NaN, prof); continue
            end
            y5 = tadd(y, ks, DP_b, hstep); y4 = tadd(y, ks, DP_bs, hstep)
            err = maximum(abs(y5[i] - y4[i]) / (rtol * (abs(y5[i]) + 1e-30) + (i == 3 ? rtol : 0.0)) for i in 1:3)
            if err <= 1 && y5[1] > 0
                # check crossing of r = rin or m = MBH within the step (linear interpolation)
                if y5[1] <= rin
                    w = (y[1] - rin) / (y[1] - y5[1])
                    return (:reached_rin, y[2] + w * (y5[2] - y[2]), prof)
                end
                if y5[2] <= p.MBH
                    return (:mass_exhausted, p.MBH, prof)
                end
                x += hstep; y = y5
                hstep *= min(4.0, 0.9 * err^(-0.2) + 1e-3)
                break
            else
                hstep *= max(0.1, 0.9 * err^(-0.2))
                hstep < 1e-12 && return (:failed, NaN, prof)
            end
        end
    end
    return (:failed, NaN, prof)
end
