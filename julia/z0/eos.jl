# Z = 0 equation of state: ideal H/He gas with LTE (Saha) ionization + blackbody radiation.
# Valid at the low densities of quasi-star envelopes (no degeneracy, no pressure ionization,
# no molecules — H2 is negligible above ~3500 K at these densities).
#
# v2 (AD-ready): every function is generic in its number type, so ForwardDiff Duals
# (including nested Duals for second derivatives) propagate through it. The two implicit
# solves (charge balance for n_e, and ρ(P,T)) are done in Float64 and then polished with two
# full Newton steps in the caller's number type. Starting from the exact root, k Newton steps
# make derivatives exact through order 2^k − 1 (here: through 3rd order).
# Primary API for a hydro code with (ρ, T) as variables: `eos_rhoT(ρ, T)`.

using ForwardDiff

struct Comp
    X::Float64
    Y::Float64
end
const COMP = Comp(0.76, 0.24)          # matches MESA-QUEST inlist (initial_y = 0.24)

"Strip all Dual layers → Float64."
pval(x::ForwardDiff.Dual) = pval(ForwardDiff.value(x))
pval(x::Real) = Float64(x)

"Two Newton steps on f(x) = 0 from a Float64 root x0, in the number type of f's parameters."
function newton_polish(f, x0; nsteps = 2)
    x = x0 + zero(f(x0))                # promote to the Dual type carried by f
    for _ in 1:nsteps
        fx = f(x)
        x = x - fx / ForwardDiff.derivative(f, x)
    end
    return x
end

λe3(T) = (2π * me * kB * T / h^2)^1.5  # inverse cube of electron thermal wavelength

function _saha_rates(T)
    kT = kB * T
    L = λe3(T)
    return L * exp(-χH / kT), 4L * exp(-χHe1 / kT), L * exp(-χHe2 / kT)
end

"Fractions for given ne: returns (xH, f0, y1, y2) and the charge sum."
@inline function _fractions(ne, nH, nHe, SH, S1, S2)
    xH = SH / (ne + SH)
    r1 = S1 / ne; r2 = S2 / ne
    f0 = 1 / (1 + r1 + r1 * r2)
    y1 = r1 * f0; y2 = r1 * r2 * f0
    return xH, f0, y1, y2, nH * xH + nHe * (y1 + 2y2)
end

"Float64 charge-balance solve (robust bisection in ln n_e)."
function _ne_float(ρ::Float64, T::Float64, comp)
    nH = comp.X * ρ / mH; nHe = comp.Y * ρ / mHe
    SH, S1, S2 = _saha_rates(T)
    nemax = nH + 2nHe
    lo = log(nemax * 1e-30 + 1e-300); hi = log(nemax)
    for _ in 1:200
        mid = 0.5(lo + hi)
        q = _fractions(exp(mid), nH, nHe, SH, S1, S2)[5]
        (exp(mid) > q) ? (hi = mid) : (lo = mid)
        hi - lo < 1e-14 && break
    end
    return exp(0.5(lo + hi))
end

"Number densities for given ρ, T (generic number types)."
function saha(ρ, T; comp = COMP)
    nH  = comp.X * ρ / mH
    nHe = comp.Y * ρ / mHe
    SH, S1, S2 = _saha_rates(T)
    ne0 = _ne_float(pval(ρ), pval(T), comp)
    # residual in ln ne: g(lne) = ln(charge sum) − lne
    g(lne) = log(_fractions(exp(lne), nH, nHe, SH, S1, S2)[5]) - lne
    lne = (ρ isa ForwardDiff.Dual || T isa ForwardDiff.Dual) ? newton_polish(g, log(ne0)) : log(ne0)
    ne = exp(lne)
    xH, f0, y1, y2, _ = _fractions(ne, nH, nHe, SH, S1, S2)
    return (ne = ne, nH0 = nH * (1 - xH), nHp = nH * xH, nHe0 = nHe * f0,
            nHe1 = nHe * y1, nHe2 = nHe * y2, nH = nH, nHe = nHe, xH = xH)
end

Pgas(ρ, T) = (s = saha(ρ, T); (s.nH + s.nHe + s.ne) * kB * T)
Prad(T) = arad * T^4 / 3
Ptot(ρ, T) = Pgas(ρ, T) + Prad(T)

"Specific internal energy (erg/g): translational + ionization + radiation."
function eint(ρ, T)
    s = saha(ρ, T)
    etr = 1.5 * (s.nH + s.nHe + s.ne) * kB * T
    eion = s.nHp * χH + s.nHe1 * χHe1 + s.nHe2 * (χHe1 + χHe2)
    return (etr + eion) / ρ + arad * T^4 / ρ
end

"Specific entropy (erg/g/K), Sackur–Tetrode sum over species + radiation."
function entropy(ρ, T)
    s = saha(ρ, T)
    kT = kB * T
    q(m) = (2π * m * kT / h^2)^1.5
    term(n, g, m) = pval(n) > 1e-280 ? n * (2.5 + log(g * q(m) / n)) : zero(n)
    S = term(s.nH0, 2, mH) + term(s.nHp, 1, mH) + term(s.nHe0, 1, mHe) +
        term(s.nHe1, 2, mHe) + term(s.nHe2, 1, mHe) + term(s.ne, 2, me)
    return kB * S / ρ + 4arad * T^3 / (3ρ)
end

"Solve P_total(ρ, T) = P for ρ (generic). Returns NaN if P ≤ P_rad(T)."
function rho_PT(P, T)
    Pf, Tf = pval(P), pval(T)
    Pg = Pf - Prad(Tf)
    Pg <= 0 && return oftype(P * T, NaN)
    lo = log(Pg * 0.5mH / (kB * Tf)); hi = log(Pg * 1.4mH / (kB * Tf))
    lr = 0.5(lo + hi)
    for _ in 1:200
        f = log(Pgas(exp(lr), Tf)) - log(Pg)
        abs(f) < 1e-14 && break
        f > 0 ? (hi = lr) : (lo = lr)
        dfdlr = ForwardDiff.derivative(x -> log(Pgas(exp(x), Tf)), lr)
        lrn = lr - f / dfdlr
        lr = (lo < lrn < hi) ? lrn : 0.5(lo + hi)
    end
    if P isa ForwardDiff.Dual || T isa ForwardDiff.Dual
        g(x) = log(Ptot(exp(x), T)) - log(P)
        lr = newton_polish(g, lr)
    end
    return exp(lr)
end

# --- quantities from a precomputed Saha state (so one Saha solve serves P, e and s) ---------
_Pgas(s, T) = (s.nH + s.nHe + s.ne) * kB * T
_eint(s, ρ, T) = (1.5 * (s.nH + s.nHe + s.ne) * kB * T + s.nHp * χH + s.nHe1 * χHe1 +
                  s.nHe2 * (χHe1 + χHe2)) / ρ + arad * T^4 / ρ
function _entropy(s, ρ, T)
    kT = kB * T
    q(m) = (2π * m * kT / h^2)^1.5
    term(n, g, m) = pval(n) > 1e-280 ? n * (2.5 + log(g * q(m) / n)) : zero(n)
    S = term(s.nH0, 2, mH) + term(s.nHp, 1, mH) + term(s.nHe0, 1, mHe) +
        term(s.nHe1, 2, mHe) + term(s.nHe2, 1, mHe) + term(s.ne, 2, me)
    return kB * S / ρ + 4arad * T^3 / (3ρ)
end

struct EOSTag end

"""
Full EOS state at (ρ, T) with exact first derivatives, for hydro codes using (ρ, T) as variables.
One Saha solve carries a 2-partial Dual in (ln ρ, ln T) through P, e and s.
Returns P, e, s, ∂P/∂ρ|T, ∂P/∂T|ρ, ∂e/∂ρ|T, c_V = ∂e/∂T|ρ, χρ, χT, Γ1, ∇ad, δ, c_P, c_s, β.
Generic: it can itself be differentiated by ForwardDiff (nested Duals, unique inner tag).
"""
function eos_rhoT(ρ, T)
    V = promote_type(typeof(ρ), typeof(T), Float64)
    tag = ForwardDiff.Tag(EOSTag(), V)  # registers the tag's ordering number before any inner tag
    D = ForwardDiff.Dual{typeof(tag)}
    ρd = D(V(ρ), V(ρ), zero(V))          # ∂/∂lnρ
    Td = D(V(T), zero(V), V(T))          # ∂/∂lnT
    sd = saha(ρd, Td)
    Pd = _Pgas(sd, Td) + arad * Td^4 / 3
    ed = _eint(sd, ρd, Td)
    Sd = _entropy(sd, ρd, Td)
    val(x) = ForwardDiff.value(x); par(x, i) = ForwardDiff.partials(x, i)
    P = val(Pd); e = val(ed); s = val(Sd)
    χρ = par(Pd, 1) / P; χT = par(Pd, 2) / P
    Sρ = par(Sd, 1); ST = par(Sd, 2)
    ∇ad = 1 / (χT - χρ * ST / Sρ)
    Γ1  = χρ - χT * Sρ / ST
    δ   = χT / χρ
    cP  = ST - Sρ * δ
    return (P = P, e = e, s = s, dPdρ = par(Pd, 1) / ρ, dPdT = par(Pd, 2) / T,
            dedρ = par(ed, 1) / ρ, cV = par(ed, 2) / T,
            χρ = χρ, χT = χT, ∇ad = ∇ad, Γ1 = Γ1, δ = δ, cP = cP, cs = sqrt(Γ1 * P / ρ),
            β = 1 - arad * T^4 / (3P))
end

"Reference implementation (three separate gradients); used only in tests."
function eos_rhoT_ref(ρ, T)
    x = [log(ρ), log(T)]
    gP = ForwardDiff.gradient(v -> log(Ptot(exp(v[1]), exp(v[2]))), x)
    gS = ForwardDiff.gradient(v -> entropy(exp(v[1]), exp(v[2])), x)
    gE = ForwardDiff.gradient(v -> eint(exp(v[1]), exp(v[2])), x)
    P = Ptot(ρ, T)
    χρ, χT = gP[1], gP[2]; Sρ, ST = gS[1], gS[2]
    return (P = P, e = eint(ρ, T), s = entropy(ρ, T), dPdρ = χρ * P / ρ, dPdT = χT * P / T,
            dedρ = gE[1] / ρ, cV = gE[2] / T, Γ1 = χρ - χT * Sρ / ST, ∇ad = 1 / (χT - χρ * ST / Sρ))
end

"Thermodynamic derivatives at (P, T) — same fields as v1, now via AD (exact)."
function thermo(P, T)
    ρ = rho_PT(P, T)
    isnan(pval(ρ)) && return nothing
    st = eos_rhoT(ρ, T)
    return (ρ = ρ, ∇ad = st.∇ad, δ = st.δ, Γ1 = st.Γ1, cP = st.cP, cs = st.cs, β = st.β)
end
