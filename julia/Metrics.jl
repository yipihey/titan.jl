"""
    Metrics

Tier-(c) benchmark observables for SOD (planar Riemann problem) and the
SUBC/SUPC radiating shocks, computed from light per-step snapshots and
linearly interpolated in time to fixed target times (TITAN never lands on a
prescribed time, and two implementations take different time steps).
Definitions are deliberately simple and documented in PORTING_NOTES.md §P2.3.
"""
module Metrics

const CAR = 7.5659141e-15        # radiation constant a (titan.dat)

struct Snap
    jstep::Int; time::Float64
    r::Vector{Float64}           # interfaces ngrs..ngre+1
    rho::Vector{Float64}; T::Vector{Float64}; Er::Vector{Float64}; p::Vector{Float64}
    u::Vector{Float64}           # interfaces
    tetot::Float64; tee::Float64; tew::Float64
end

"Snapshot of the live model (needs the TitanFortran module `TF`)."
function snap(TF)
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    I = a:b+1; Z = a:b
    Snap(Int(TF.getscalar(:jstep)), TF.getscalar(:timen), TF.var(:rn)[I], TF.var(:dn)[Z],
         TF.var(:tn)[Z], TF.var(:ern)[Z], TF.var(:pgn)[Z], TF.var(:un)[I],
         TF.getscalar(:tetot), TF.getscalar(:tee), TF.getscalar(:tew))
end

centers(s) = 0.5 .* (s.r[1:end-1] .+ s.r[2:end])
Trad(s) = (max.(s.Er, 0.0) ./ CAR) .^ 0.25
"linear interpolation of y(x) at x0 (x increasing)"
function lin(x, y, x0)
    x0 <= x[1] && return y[1]; x0 >= x[end] && return y[end]
    i = searchsortedlast(x, x0)
    y[i] + (y[i+1] - y[i]) * (x0 - x[i]) / (x[i+1] - x[i])
end
"position where y crosses level going from index i to i+1 (linear in x)"
cross(x, y, i, level) = x[i] + (x[i+1] - x[i]) * (level - y[i]) / (y[i+1] - y[i])

"Position of the steepest outward density drop, refined to the half-jump crossing."
function front(rc, ρ; kmin = 1, kmax = length(ρ) - 1)
    # steepest drop per unit LENGTH (per-zone drops are dominated by the coarse
    # zones the adaptive grid leaves in smooth regions, e.g. a rarefaction fan)
    k = kmin - 1 + argmax([(ρ[i] - ρ[i+1]) / (rc[i+1] - rc[i]) for i in kmin:kmax])
    hi = maximum(ρ[max(1, k-2):k]); lo = minimum(ρ[k+1:min(end, k+3)])
    lev = 0.5 * (hi + lo)
    for i in max(1, k-3):min(length(ρ)-1, k+3)
        (ρ[i] - lev) * (ρ[i+1] - lev) <= 0 && ρ[i] != ρ[i+1] && return (cross(rc, ρ, i, lev), k, hi, lo)
    end
    (rc[k], k, hi, lo)
end

"Radiating-shock observables (SUBC/SUPC). T_thr: precursor-front threshold [K]."
function shock_metrics(s; T_thr = 170.0)
    rc = centers(s); Tr = Trad(s)
    rsh, k, ρpost, ρpre = front(rc, s.rho)
    rp = s.r[1]; W = rsh - rp                      # shocked-layer width (piston..shock)
    kspike = argmax(s.T)
    iprec = findlast(>(T_thr), s.T)
    (; time = s.time, r_shock = rsh, width = W, rho_post = ρpost, rho_pre = ρpre,
       T_spike = s.T[kspike], r_spike = rc[kspike], Trad_shock = lin(rc, Tr, rsh),
       T_post = lin(rc, s.T, rsh - 0.1W), Trad_post = lin(rc, Tr, rsh - 0.1W),
       T_pre = lin(rc, s.T, rsh + 0.05W), Trad_pre = lin(rc, Tr, rsh + 0.05W),
       precursor = iprec === nothing ? 0.0 : max(0.0, rc[iprec] - rsh),
       E_err = 0.0)
end

"SOD observables: shock, contact, rarefaction head/tail positions, plateau p*, u*, ρ*L, ρ*R."
function sod_metrics(s)
    rc = centers(s); ρ = s.rho; p = s.p; n = length(ρ)
    dx = [rc[i+1] - rc[i] for i in 1:n-1]
    # shock: steepest pressure drop per length (p is continuous at the contact)
    ksh = argmax([(p[i] - p[i+1]) / dx[i] for i in 1:n-1])
    plev = 0.5 * (maximum(p[max(1, ksh-3):ksh]) + minimum(p[ksh+1:min(n, ksh+4)]))
    rsh = rc[ksh]
    for i in max(1, ksh-4):min(n-1, ksh+4)
        (p[i] - plev) * (p[i+1] - plev) <= 0 && p[i] != p[i+1] && (rsh = cross(rc, p, i, plev); break)
    end
    # contact: density drops but pressure does not: maximise (Δlnρ − Δlnp)/Δx left of the shock
    score = [((log(ρ[i]) - log(ρ[i+1])) - (log(p[i]) - log(p[i+1]))) / dx[i] for i in 1:n-1]
    sel = [rc[i] < rsh - 1e-6 ? score[i] : -Inf for i in 1:n-1]
    kct = argmax(sel)
    ρL2 = lin(rc, ρ, rc[kct] - 0.02); ρR2 = lin(rc, ρ, 0.5 * (rc[kct] + rsh))
    rlev = 0.5 * (ρL2 + ρR2); rct = rc[kct]
    for i in max(1, kct-4):min(n-1, kct+4)
        (ρ[i] - rlev) * (ρ[i+1] - rlev) <= 0 && ρ[i] != ρ[i+1] && (rct = cross(rc, ρ, i, rlev); break)
    end
    # head/tail: interpolated crossings of 0.99 ρ_L and 1.01 ρ*_L (not zone-quantised)
    ihead = findfirst(<(0.99 * ρ[1]), ρ)
    ihead = ihead === nothing || ihead == 1 ? ihead : cross(rc, ρ, ihead - 1, 0.99 * ρ[1])
    itail = findfirst(<=(1.01 * ρL2), ρ)
    itail = itail === nothing || itail == 1 ? itail : cross(rc, ρ, itail - 1, 1.01 * ρL2)
    mid = 0.5 * (rct + rsh)
    (; time = s.time, r_shock = rsh, r_contact = rct,
       r_head = ihead === nothing ? NaN : ihead, r_tail = itail === nothing ? NaN : itail,
       p_star = lin(rc, p, mid), u_star = lin(s.r, s.u, mid), rho_2 = ρL2, rho_3 = ρR2, E_err = 0.0)
end

"Interpolate the metric NamedTuples (one per step) to time t0."
function at_time(ms::Vector, t0)
    ts = [m.time for m in ms]
    i = searchsortedlast(ts, t0)
    (i < 1 || i >= length(ms)) && return nothing
    a, b = ms[i], ms[i+1]; w = (t0 - a.time) / (b.time - a.time)
    NamedTuple{keys(a)}(Tuple(k == :time ? t0 : a[k] + w * (b[k] - a[k]) for k in keys(a)))
end

"Energy error relative to the larger of |internal+kinetic+rad| and |boundary work|."
energy_err(s, s0) = (s.tetot - s0.tetot) / max(abs(s.tee), abs(s.tew), abs(s0.tee))

"Exact Sod solution (γ = 1.4, ρ,p = (1,1)|(0.125,0.1), diaphragm x0) — reference positions."
function sod_exact(t; x0 = 1.5)
    γ = 1.4; ρL, pL, ρR, pR = 1.0, 1.0, 0.125, 0.1
    cL = sqrt(γ * pL / ρL); cR = sqrt(γ * pR / ρR)
    f(p, ρk, pk, ck) = p > pk ? (p - pk) * sqrt((2 / ((γ + 1) * ρk)) / (p + (γ - 1) / (γ + 1) * pk)) :
                                  2ck / (γ - 1) * ((p / pk)^((γ - 1) / 2γ) - 1)
    lo, hi = 1e-6, 1.0
    for _ in 1:200
        m = 0.5(lo + hi); (f(m, ρL, pL, cL) + f(m, ρR, pR, cR) > 0) ? (hi = m) : (lo = m)
    end
    ps = 0.5(lo + hi); us = 0.5 * (f(ps, ρR, pR, cR) - f(ps, ρL, pL, cL))  # uL = uR = 0
    S = sqrt((γ + 1) / 2γ * ps / pR + (γ - 1) / 2γ) * cR
    cs = cL * (ps / pL)^((γ - 1) / 2γ)
    (; p_star = ps, u_star = us, rho_2 = ρL * (ps / pL)^(1 / γ),
       rho_3 = ρR * ((ps / pR) + (γ - 1) / (γ + 1)) / ((γ - 1) / (γ + 1) * (ps / pR) + 1),
       r_shock = x0 + S * t, r_contact = x0 + us * t, r_head = x0 - cL * t, r_tail = x0 + (us - cs) * t)
end

end # module
