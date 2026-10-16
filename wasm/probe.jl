module Probe
using ..QSW
const Q = QSW
p01() = Q.params(0.1, 0.8, 1e-3)
thermo_nad(P::Float64, T::Float64)::Float64 = (th = Q.thermo(P, T); th === nothing ? NaN : th.∇ad)
mlt_nabla(P::Float64, T::Float64, r::Float64, m::Float64)::Float64 = (p = p01(); Q.mlt(P, T, r, m, p.L, p).∇)
bc_Pb(R::Float64, Mdot::Float64)::Float64 = (bc = Q.outer_bc(R, p01(); Mdot = Mdot * Q.Msun / Q.yr); bc === nothing ? NaN : bc.Pb)
bc_tau(R::Float64, Mdot::Float64)::Float64 = Q.infall_tau(Mdot * Q.Msun / Q.yr, 1e6Q.Msun, R, 1.4e44, p01().tab).τ
function integ(R::Float64, Pb::Float64, Tb::Float64)::Float64
    st, m, _ = Q.integrate_inward(R, Pb, Tb, p01())
    st == :failed ? NaN : m
end
end
module Probe2
using ..QSW
const Q = QSW
p01() = Q.params(0.1, 0.8, 1e-3)
rhs_dr(x::Float64, r::Float64, m::Float64, lT::Float64)::Float64 = Q.rhs(x, (r, m, lT), p01())[1][1]
rhs_rho(x::Float64, r::Float64, m::Float64, lT::Float64)::Float64 = (q = Q.rhs(x, (r, m, lT), p01())[2]; hasproperty(q, :ρ) ? q.ρ : NaN)
mlt_full(P::Float64, T::Float64, r::Float64, m::Float64)::Float64 = (p = p01(); q = Q.mlt(P, T, r, m, p.L, p); hasproperty(q, :v) ? q.v : -1.0)
end
module Probe3
using ..QSW
import ..QSW: rhs, tadd_n, DP_A, DP_C, DP_B, DP_BS, QSParams
const Q = QSW
const MARK = zeros(Int32, 4)        # [step, stage, phase, accepted]
function integ_dbg(R, Pb, Tb, p::QSParams; rtol = 1e-7, maxsteps = 200000)
    x = log(Pb); y = (R, p.M, log(Tb)); hstep = 0.01; rin = p.ε_in * R
    for step in 1:maxsteps
        MARK[1] = step; MARK[3] = 1
        k1, q = rhs(x, y, p)
        MARK[3] = 2
        isnan(k1[1]) && return (:failed, NaN)
        while true
            ks = Vector{NTuple{3,Float64}}(undef, 7); ks[1] = k1
            bad = false
            for s in 2:7
                MARK[2] = s; MARK[3] = 3
                ys = tadd_n(y, ks, DP_A[s], hstep, s - 1)
                MARK[3] = 4
                ks[s] = rhs(x + DP_C[s] * hstep, ys, p)[1]
                MARK[3] = 5
                if isnan(ks[s][1]) || ys[1] <= 0
                    bad = true; break
                end
            end
            MARK[3] = 6
            if bad
                hstep /= 4; hstep < 1e-12 && return (:failed, NaN); continue
            end
            y5 = tadd_n(y, ks, DP_B, hstep, 7); y4 = tadd_n(y, ks, DP_BS, hstep, 7)
            MARK[3] = 7
            err = maximum(abs(y5[i] - y4[i]) / (rtol * (abs(y5[i]) + 1e-30) + (i == 3 ? rtol : 0.0)) for i in 1:3)
            MARK[3] = 8
            if err <= 1 && y5[1] > 0
                y5[1] <= rin && return (:reached_rin, y[2] + (y[1] - rin) / (y[1] - y5[1]) * (y5[2] - y[2]))
                y5[2] <= p.MBH && return (:mass_exhausted, p.MBH)
                x += hstep; y = y5; MARK[4] += 1
                hstep *= min(4.0, 0.9 * err^(-0.2) + 1e-3)
                break
            else
                hstep *= max(0.1, 0.9 * err^(-0.2))
                hstep < 1e-12 && return (:failed, NaN)
            end
        end
    end
    return (:failed, NaN)
end
function integd(R::Float64, Pb::Float64, Tb::Float64)::Float64
    st, m = integ_dbg(R, Pb, Tb, Q.params(0.1, 0.8, 1e-3))
    st == :failed ? NaN : m
end
mark(i::Int32)::Int32 = MARK[i]
end
module Roots2
using ..QSW, ..Probe3
const Q = QSW
const ROOTS = zeros(16); const NROOT = Ref(0)
function sgn(R, p, Mdot)
    bc = Q.outer_bc(R, p; Mdot = Mdot)
    bc === nothing && return 0
    st, m = Probe3.integ_dbg(R, bc.Pb, bc.Tb, p)
    st == :failed && return 0
    st == :mass_exhausted && return -1
    m > p.MBH ? 1 : -1
end
function all_roots2(fBH::Float64, Mdot_msun_yr::Float64, α::Float64, ε::Float64)::Int32
    p = Q.params(fBH, α, ε); Mdot = Mdot_msun_yr * Q.Msun / Q.yr
    n = 105; lR0 = 14.0; dlR = 0.025; NROOT[] = 0
    sprev = sgn(10^lR0, p, Mdot)
    for i in 1:n-1
        a = lR0 + dlR * (i - 1); b = lR0 + dlR * i
        sb = sgn(10^b, p, Mdot)
        if sprev != 0 && sb != 0 && sprev != sb && NROOT[] < 16
            lo, hi, sa = a, b, sprev
            for _ in 1:40
                mid = 0.5(lo + hi); s = sgn(10^mid, p, Mdot)
                s == 0 && (hi = mid; continue)
                s == sa ? (lo = mid) : (hi = mid)
            end
            NROOT[] += 1; ROOTS[NROOT[]] = 0.5(lo + hi)
        end
        sprev = sb
    end
    Int32(NROOT[])
end
root2(k::Int32)::Float64 = 10^ROOTS[k]
end
