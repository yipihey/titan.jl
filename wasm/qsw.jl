# The steady quasi-star model (titan.jl/steady/src) for WasmTarget: the original sources, with the one EOS patch of eosw.jl.
module QSW
using ForwardDiff
const SRC = "/home/claude/titan.jl/steady/src"
include(joinpath(SRC, "constants.jl"))
let src = read(joinpath(SRC, "eos.jl"), String)
    ps = ["@inline function _fractions" => "function _fractions",
          # thermo(): derivatives from ForwardDiff.jacobian instead of hand-seeded EOSTag Duals (see eos_rhoT_jac)
          "    st = eos_rhoT(ρ, T)\n    return (ρ = ρ," => "    st = eos_rhoT_jac(ρ, T)\n    return (ρ = ρ,"]
    for p in ps; occursin(p.first, src) || error("patch target missing: " * p.first); src = replace(src, p); end
    include_string(QSW, src)
end
"Same quantities as eos_rhoT, from two ForwardDiff.gradient calls in (ln ρ, ln T) exactly as eos_rhoT_ref does
(log P and s), instead of hand-seeded EOSTag Duals. Agrees with eos_rhoT to rounding."
function eos_rhoT_jac(ρ, T)
    x = [log(ρ), log(T)]
    gP = ForwardDiff.gradient(v -> log(Ptot(exp(v[1]), exp(v[2]))), x)
    gS = ForwardDiff.gradient(v -> entropy(exp(v[1]), exp(v[2])), x)
    P = Ptot(ρ, T)
    χρ, χT = gP[1], gP[2]; Sρ, ST = gS[1], gS[2]
    ∇ad = 1 / (χT - χρ * ST / Sρ); Γ1 = χρ - χT * Sρ / ST; δ = χT / χρ; cP = ST - Sρ * δ
    return (P = P, ∇ad = ∇ad, Γ1 = Γ1, δ = δ, cP = cP, cs = sqrt(Γ1 * P / ρ), β = 1 - arad * T^4 / (3P))
end
include(joinpath(SRC, "opacity.jl"))
let src = read(joinpath(SRC, "envelope.jl"), String)
    # Type-stable Dormand–Prince stages (same coefficients as Float64, same summation order: bit-identical results).
    # The original indexes a heterogeneous tuple of tuples with a runtime stage number (dynamic dispatch) and
    # copies ks[1:s-1] (memmove); WasmTarget lowers neither.
    ps = ["tadd(y, ks[1:s-1], DP_a[s], hstep)" => "tadd_n(y, ks, DP_A[s], hstep, s - 1)",
          "ks[s] = rhs(x + DP_c[s] * hstep, ys, p)[1]" => "ks[s] = rhs(x + DP_C[s] * hstep, ys, p)[1]",
          # DP_b / DP_bs contain Int zeros: a Union{Int64,Float64} element in the generator traps in Wasm (illegal cast)
          "y5 = tadd(y, ks, DP_b, hstep); y4 = tadd(y, ks, DP_bs, hstep)" => "y5 = tadd_n(y, ks, DP_B, hstep, 7); y4 = tadd_n(y, ks, DP_BS, hstep, 7)"]
    for p in ps; occursin(p.first, src) || error("patch target missing: " * p.first); src = replace(src, p); end
    src *= """
    const DP_C = Float64.(DP_c)
    const DP_B = Float64.(DP_b); const DP_BS = Float64.(DP_bs)
    const DP_A = ntuple(s -> ntuple(j -> j < s ? Float64(DP_a[s][j]) : 0.0, 6), 7)
    tadd_n(y, ks, coeffs, h, n) = ntuple(i -> (acc = 0.0; for j in 1:n; acc += coeffs[j] * ks[j][i]; end; y[i] + h * acc), 3)
    """
    include_string(QSW, src)
end
include(joinpath(SRC, "boundary.jl"))

# ---- browser-side table: allocated and filled at run time from kaptable.bin
const TAB = Ref(KapTable(zeros(1), zeros(1), zeros(1, 1), zeros(1, 1)))
function table_init(nT::Int32, nR::Int32, lT0::Float64, dlT::Float64, lR0::Float64, dlR::Float64)::Int32
    lT = [lT0 + dlT * (i - 1) for i in 1:nT]; lR = [lR0 + dlR * (j - 1) for j in 1:nR]
    TAB[] = KapTable(lT, lR, zeros(Int(nT), Int(nR)), zeros(Int(nT), Int(nR)))
    Int32(0)
end
table_set(i::Int32, j::Int32, lk::Float64, lkP::Float64)::Int32 = (TAB[].lk[i, j] = lk; TAB[].lkP[i, j] = lkP; Int32(0))
kap_R(ρ::Float64, T::Float64)::Float64 = kap(TAB[], ρ, T)

params(fBH, α, ε) = (M = 1e6Msun; QSParams(M = M, MBH = fBH * M, L = LEdd_es(M), αMLT = α, ε_in = ε, tab = TAB[]))
"m(r_in) for trial R (NaN on failure; −1 if mass exhausted) — the shooting function."
function shoot(R::Float64, fBH::Float64, Mdot_msun_yr::Float64, α::Float64, ε::Float64)::Float64
    p = params(fBH, α, ε)
    bc = outer_bc(R, p; Mdot = Mdot_msun_yr * Msun / yr)
    bc === nothing && return NaN
    st, mrin, _ = integrate_inward(R, bc.Pb, bc.Tb, p)
    st == :failed ? NaN : st == :mass_exhausted ? -1.0 : mrin
end
"R of the first root (as solve_R does), or NaN."
function first_root(fBH::Float64, Mdot_msun_yr::Float64, α::Float64, ε::Float64)::Float64
    s = solve_R(params(fBH, α, ε); Mdot = Mdot_msun_yr * Msun / yr, lRmin = 15.0, lRmax = 16.8, nscan = 37)
    s === nothing ? NaN : s.R
end
end

module Roots
using ..QSW
const QS = QSW
const ROOTS = zeros(16)
const NROOT = Ref(0)
"All steady solutions on the scan log R = 14:0.025:16.6 with 40 bisections per sign change (as steady/scripts/branches.jl,
using the original shoot_sign). Stores the radii; returns their number."
function all_roots(fBH::Float64, Mdot_msun_yr::Float64, α::Float64, ε::Float64)::Int32
    p = QS.params(fBH, α, ε); Mdot = Mdot_msun_yr * QS.Msun / QS.yr
    n = 105; lR0 = 14.0; dlR = 0.025
    NROOT[] = 0
    sprev = QS.shoot_sign(10^lR0, p; Mdot = Mdot)[1]
    for i in 1:n-1
        a = lR0 + dlR * (i - 1); b = lR0 + dlR * i
        sb = QS.shoot_sign(10^b, p; Mdot = Mdot)[1]
        if sprev != 0 && sb != 0 && sprev != sb && NROOT[] < 16
            lo, hi, sa = a, b, sprev
            for _ in 1:40
                mid = 0.5(lo + hi); s = QS.shoot_sign(10^mid, p; Mdot = Mdot)[1]
                s == 0 && (hi = mid; continue)
                s == sa ? (lo = mid) : (hi = mid)
            end
            NROOT[] += 1; ROOTS[NROOT[]] = 0.5(lo + hi)
        end
        sprev = sb
    end
    Int32(NROOT[])
end
root(k::Int32)::Float64 = 10^ROOTS[k]
end
