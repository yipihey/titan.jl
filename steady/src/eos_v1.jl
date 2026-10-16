# Z = 0 equation of state: ideal H/He gas with LTE (Saha) ionization + blackbody radiation.
# Valid at the low densities of quasi-star envelopes (no degeneracy, no pressure ionization,
# no molecules — H2 is negligible above ~3500 K at these densities, see notes).

struct Comp
    X::Float64
    Y::Float64
end
const COMP = Comp(0.76, 0.24)          # matches MESA-QUEST inlist (initial_y = 0.24)

λe3(T) = (2π * me * kB * T / h^2)^1.5  # inverse cube of electron thermal wavelength

"Number densities for given ρ, T. Returns NamedTuple."
function saha(ρ, T; comp = COMP)
    nH  = comp.X * ρ / mH
    nHe = comp.Y * ρ / mHe
    kT  = kB * T
    L   = λe3(T)
    SH  = L * exp(-χH / kT)            # g-factor ratio 1
    S1  = 4L * exp(-χHe1 / kT)
    S2  = L * exp(-χHe2 / kT)
    # charge balance f(ne) = ne - [nH xH + nHe (y1 + 2 y2)] = 0, monotone in ne; bisect in log
    function charges(ne)
        xH = SH / (ne + SH)
        r1 = S1 / ne; r2 = S2 / ne
        f0 = 1 / (1 + r1 + r1 * r2)
        y1 = r1 * f0; y2 = r1 * r2 * f0
        return nH * xH + nHe * (y1 + 2y2), xH, f0, y1, y2
    end
    nemax = nH + 2nHe
    lo = log(nemax * 1e-30 + 1e-300); hi = log(nemax)
    for _ in 1:200
        mid = 0.5(lo + hi)
        q = charges(exp(mid))[1]
        (exp(mid) > q) ? (hi = mid) : (lo = mid)
        hi - lo < 1e-13 && break
    end
    ne = exp(0.5(lo + hi))
    _, xH, f0, y1, y2 = charges(ne)
    return (ne = ne, nH0 = nH * (1 - xH), nHp = nH * xH, nHe0 = nHe * f0,
            nHe1 = nHe * y1, nHe2 = nHe * y2, nH = nH, nHe = nHe, xH = xH)
end

Pgas(ρ, T) = (s = saha(ρ, T); (s.nH + s.nHe + s.ne) * kB * T)
Prad(T) = arad * T^4 / 3

"Solve P_total(ρ, T) = P for ρ. Returns NaN if P < P_rad(T)."
function rho_PT(P, T)
    Pg = P - Prad(T)
    Pg <= 0 && return NaN
    # initial guess with μ = 0.6 / 1.25 bracket, Newton in log ρ with bisection safeguard
    lo = log(Pg * 0.5mH / (kB * T)); hi = log(Pg * 1.4mH / (kB * T))
    lr = 0.5(lo + hi)
    for _ in 1:100
        f = log(Pgas(exp(lr), T)) - log(Pg)
        abs(f) < 1e-13 && break
        f > 0 ? (hi = lr) : (lo = lr)
        d = 1e-6
        dfdlr = (log(Pgas(exp(lr + d), T)) - log(Pgas(exp(lr - d), T))) / 2d
        lrn = lr - f / dfdlr
        lr = (lo < lrn < hi) ? lrn : 0.5(lo + hi)
    end
    return exp(lr)
end

"Specific entropy (erg/g/K), Sackur–Tetrode sum over species + radiation."
function entropy(ρ, T)
    s = saha(ρ, T)
    kT = kB * T
    q(m) = (2π * m * kT / h^2)^1.5
    term(n, g, m) = n > 1e-280 ? n * (2.5 + log(g * q(m) / n)) : 0.0
    S = term(s.nH0, 2, mH) + term(s.nHp, 1, mH) + term(s.nHe0, 1, mHe) +
        term(s.nHe1, 2, mHe) + term(s.nHe2, 1, mHe) + term(s.ne, 2, me)
    return kB * S / ρ + 4arad * T^3 / (3ρ)
end

Ptot(ρ, T) = Pgas(ρ, T) + Prad(T)

"""
Thermodynamic derivatives at (P, T). Derivatives are taken in (lnρ, lnT), where the EOS is
smooth even when β = Pgas/P ≪ 1, then converted:
  χρ = ∂lnP/∂lnρ|T, χT = ∂lnP/∂lnT|ρ, Sρ = ∂S/∂lnρ|T, ST = ∂S/∂lnT|ρ
  ∇ad = 1/(χT − χρ ST/Sρ),  Γ1 = χρ − χT Sρ/ST,  δ = χT/χρ,  cP = ST − Sρ δ
"""
function thermo(P, T)
    ρ = rho_PT(P, T)
    isnan(ρ) && return nothing
    d = 1e-4
    lρ, lT = log(ρ), log(T)
    lnP(a, b) = log(Ptot(exp(a), exp(b)))
    S(a, b)   = entropy(exp(a), exp(b))
    χρ = (lnP(lρ + d, lT) - lnP(lρ - d, lT)) / 2d
    χT = (lnP(lρ, lT + d) - lnP(lρ, lT - d)) / 2d
    Sρ = (S(lρ + d, lT) - S(lρ - d, lT)) / 2d
    ST = (S(lρ, lT + d) - S(lρ, lT - d)) / 2d
    ∇ad = 1 / (χT - χρ * ST / Sρ)
    Γ1  = χρ - χT * Sρ / ST
    δ   = χT / χρ
    cP  = ST - Sρ * δ
    return (ρ = ρ, ∇ad = ∇ad, δ = δ, Γ1 = Γ1, cP = cP, cs = sqrt(Γ1 * P / ρ),
            β = 1 - Prad(T) / P)
end
