"""
    StellarEOS

Gas equation of state for Pop III stellar interiors and protostars (S2): ρ = 1e-12 … 1e4 g cm⁻³, T = 3e3 … 3e9 K (below ~2.5e3 K Z0Physics' Saha underflows; no molecules).
Composition by mass fractions X (H), Y (He), Z (C/N/O mix `ZMIX`), all of them at the caller's choice; X + Y + Z = 1.
Radiation is NOT in the gas EOS (TITAN evolves E_r): `eos` (unless `radiation = true`) and `gas_pe` return gas-only quantities; `Prad`, `erad`,
`srad` are separate. Everything is generic in the number type (ForwardDiff Duals pass through, nested as in Z0Physics).

Two thermodynamically complete branches, each giving p, e, s per unit mass at (ρ, T):
* S (Saha): Z0Physics' LTE H/He gas (ideal, no degeneracy), generalised to arbitrary X, Y (`Z0Physics.Comp`); metals are
  neutral ideal atoms. Valid for ρ ≲ 1e-2, T ≲ 3e5 K.
* E (fully ionised): ideal ions (H⁺, He²⁺, C, N, O bare nuclei, Sackur–Tetrode, g = 1 as in Z0Physics) + electrons (and,
  optionally, positrons) with arbitrary degeneracy and relativity from the generalised Fermi–Dirac integrals
  F_k(η, θ) = ∫₀^∞ x^k (1 + θx/2)^{1/2} / (e^{x−η} + 1) dx, θ = kT/mc², by composite Gauss–Legendre quadrature (see `fd`);
  plus the Yakovlev & Shalybkov (1989) uniform-background Coulomb correction (as in Timmes–Swesty). Energy zero:
  neutral atoms (as Z0Physics), i.e. E includes the total ionisation energy of every ion.
The two are joined by blending the FREE ENERGY, f = (1−w) f_S + w f_E, with w(ρ, T) = 1 − (1−w_T)(1−w_ρ) a C² polynomial
in log T (log T = 5.3 … 5.9) and in log ρ (pressure ionisation, log ρ = −1 … 0.7): p = ρ²∂f/∂ρ, s = −∂f/∂T, e = f + Ts.
Hence p, e, s are exactly consistent (Maxwell relations hold to round-off) and Γ1, ∇_ad are continuous (f ∈ C²).

API: `eos(ρ, T; X, Y, Z)` → NamedTuple (P, e, s, cV, χρ, χT, Γ1, ∇ad, δ, cP, cs, ne, η, …); `gas_pe(ρ, T; …)` → (p, e)
(the residual's hot path); `rho_PT(P, T; …)` inverts; `Prad`, `erad`, `srad`.
"""
module StellarEOS

import ..Z0Physics
using ForwardDiff
const Z0 = Z0Physics
const kB = Z0.kB; const h = Z0.h; const me = Z0.me; const cl = Z0.c; const mH = Z0.mH; const mHe = Z0.mHe
const arad = Z0.arad; const eV = Z0.eV; const χH = Z0.χH; const χHe1 = Z0.χHe1; const χHe2 = Z0.χHe2
const amu = 1.66053906660e-24
const qe = 4.80320471e-10                    # esu
const mc2 = me * cl^2
const A0 = 8π * sqrt(2.0) * (me * cl / h)^3  # n = A0 θ^{3/2} (F_{1/2} + θ F_{3/2}), P = A0 mc² (2/3) θ^{5/2} (…)

# =====================================================================================================================
# Generalised Fermi–Dirac integrals
# =====================================================================================================================
"Gauss–Legendre nodes and weights on [−1, 1] (Newton on P_n)."
function gauleg(n)
    x = zeros(n); w = zeros(n)
    for i in 1:n
        z = cos(π * (i - 0.25) / (n + 0.5))
        for _ in 1:100
            p1 = 1.0; p2 = 0.0
            for j in 1:n; p3 = p2; p2 = p1; p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j; end
            pp = n * (z * p1 - p2) / (z^2 - 1)
            dz = p1 / pp; z -= dz
            abs(dz) < 1e-16 && break
        end
        p1 = 1.0; p2 = 0.0
        for j in 1:n; p3 = p2; p2 = p1; p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j; end
        pp = n * (z * p1 - p2) / (z^2 - 1)
        x[i] = z; w[i] = 2 / ((1 - z^2) * pp^2)
    end
    x, w
end

# Master grid for η < ETAC: panels on x ∈ [0, 95]; the first one in s = √x (removes the x^{1/2} endpoint behaviour),
# the others in x with width ≤ 5 (the Fermi factor has poles at x = η ± iπ). 16 nodes per panel.
const NGL = 16
const ETAC = 48.0
const EDGES = [0.0, 1, 3, 6, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95]
struct LowGrid
    x::Vector{Float64}; c12::Vector{Float64}; c32::Vector{Float64}; c52::Vector{Float64}; last::Vector{Int}
end
function lowgrid()
    t, w = gauleg(NGL)
    xs = Float64[]; c12 = Float64[]; c32 = Float64[]; c52 = Float64[]; last = Int[]
    for p in 1:length(EDGES)-1
        a, b = EDGES[p], EDGES[p+1]
        for j in 1:NGL
            if p == 1                                   # x = s², s ∈ [0, 1]: dx = 2s ds, x^{1/2} = s
                s = 0.5 * (t[j] + 1); ws = 0.5 * w[j]; x = s^2; W = 2s * ws; sx = s
            else
                x = a + (b - a) * 0.5 * (t[j] + 1); W = (b - a) * 0.5 * w[j]; sx = sqrt(x)
            end
            push!(xs, x); push!(c12, W * sx); push!(c32, W * sx * x); push!(c52, W * sx * x^2)
        end
        push!(last, length(xs))
    end
    LowGrid(xs, c12, c32, c52, last)
end
const LOW = lowgrid()

# High-η grid: F_k = G_k(η) + ∫₀^U [g(η+u) − g(η−u)] /(e^u + 1) du, ∂F_k/∂η = ∫₀^U [g(η+u) + g(η−u)] w(u) du,
# w = e^u/(e^u+1)²; G_k = ∫₀^η g dx = 2η^{k+1} ∫₀¹ t^{2k+1} (1 + θη t²/2)^{1/2} dt.
const UEDGES = [0.0, 2, 5, 10, 18, 30, 46]
struct HighGrid
    u::Vector{Float64}; W::Vector{Float64}; fu::Vector{Float64}; wu::Vector{Float64}
    t::Vector{Float64}; wt::Vector{Float64}
end
function highgrid()
    t, w = gauleg(NGL); us = Float64[]; Ws = Float64[]
    for p in 1:length(UEDGES)-1
        a, b = UEDGES[p], UEDGES[p+1]
        for j in 1:NGL; push!(us, a + (b - a) * 0.5 * (t[j] + 1)); push!(Ws, (b - a) * 0.5 * w[j]); end
    end
    fu = [1 / (exp(u) + 1) for u in us]; wu = [exp(-u) / (1 + exp(-u))^2 for u in us]
    t2, w2 = gauleg(32)
    HighGrid(us, Ws, fu, wu, 0.5 .* (t2 .+ 1), 0.5 .* w2)
end
const HIGH = highgrid()

"""
    fd(η, θ) -> (F12, F32, F52, D12, D32, D52)

F_k(η, θ) for k = ½, 3/2, 5/2 and D_k = ∂F_k/∂η (analytic quadrature sums). Generic in the number type of η and θ (θ ≥ 0).
η < 48: composite Gauss–Legendre on the master grid truncated at x = max(η, 0) + 46; η ≥ 48: edge-form around the
Fermi energy (above). Accuracy ≲ 1e-13 relative (tested against an independent double-exponential quadrature).
"""
function fd(η::Number, θ::Number)
    T = promote_type(typeof(η), typeof(θ), Float64)
    ηv = Z0.pval(η)
    if ηv < ETAC
        xmax = max(ηv, 0.0) + 46.0
        np = 1
        while EDGES[np+1] < xmax; np += 1; end
        n = LOW.last[np]
        s12 = zero(T); s32 = zero(T); s52 = zero(T); d12 = zero(T); d32 = zero(T); d52 = zero(T)
        hθ = θ / 2
        @inbounds for j in 1:n
            x = LOW.x[j]
            f = inv(exp(x - η) + 1)
            q = sqrt(1 + hθ * x)
            ft = f * q; dt = f * (1 - f) * q
            a = LOW.c12[j]; b = LOW.c32[j]; c = LOW.c52[j]
            s12 += a * ft; s32 += b * ft; s52 += c * ft
            d12 += a * dt; d32 += b * dt; d52 += c * dt
        end
        return (s12, s32, s52, d12, d32, d52)
    else
        g12 = zero(T); g32 = zero(T); g52 = zero(T)
        a = θ * η / 2
        # ∫₀¹ t^{2k+1} (1 + a t²)^{1/2} dt: one panel for a ≤ 4, else panels [0, τ], [τ, 4τ], … with τ = 1/√a (ultra-relativistic)
        av = Z0.pval(a)
        lo = 0.0
        hi = av <= 4.0 ? 1.0 : min(1.0, 1 / sqrt(av))
        while true
            @inbounds for j in eachindex(HIGH.t)
                t = lo + (hi - lo) * HIGH.t[j]; q = sqrt(1 + a * t * t); t2 = t * t; w = (hi - lo) * HIGH.wt[j] * q * 2
                g12 += w * t2; g32 += w * t2 * t2; g52 += w * t2 * t2 * t2
            end
            hi >= 1.0 && break
            lo = hi; hi = min(1.0, 4hi)
        end
        sη = sqrt(η); η2 = η * η
        G12 = η * sη * g12; G32 = η2 * sη * g32; G52 = η2 * η * sη * g52
        s12 = G12; s32 = G32; s52 = G52; d12 = zero(T); d32 = zero(T); d52 = zero(T)
        hθ = θ / 2
        @inbounds for j in eachindex(HIGH.u)
            u = HIGH.u[j]; W = HIGH.W[j]; fu = HIGH.fu[j] * W; wu = HIGH.wu[j] * W
            xp = η + u; xm = η - u
            sp = sqrt(xp); sm = sqrt(xm); qp = sqrt(1 + hθ * xp); qm = sqrt(1 + hθ * xm)
            p12 = sp * qp; m12 = sm * qm
            p32 = xp * p12; m32 = xm * m12
            p52 = xp * p32; m52 = xm * m32
            s12 += fu * (p12 - m12); s32 += fu * (p32 - m32); s52 += fu * (p52 - m52)
            d12 += wu * (p12 + m12); d32 += wu * (p32 + m32); d52 += wu * (p52 + m52)
        end
        return (s12, s32, s52, d12, d32, d52)
    end
end

# =====================================================================================================================
# Composition
# =====================================================================================================================
const ZMIX = (C = 0.20, N = 0.06, O = 0.74)           # default mass fractions inside Z (CNO)
const MC = 12.011amu; const MN = 14.007amu; const MO = 15.999amu
# total ionisation energies (eV): C 489.99, N 667.05, O 871.41
const XC = 489.99eV; const XN = 667.05eV; const XO = 871.41eV

"Per-gram species data: mass fractions, atomic masses, charges, total ionisation energies."
struct Comp
    X::Float64; Y::Float64; Z::Float64
    x::NTuple{5,Float64}; m::NTuple{5,Float64}; zc::NTuple{5,Float64}; χ::NTuple{5,Float64}
    νion::Float64      # ions per gram
    νe::Float64        # electrons per gram when fully ionised
    eion::Float64      # total ionisation energy per gram [erg/g]
    νZ::Float64        # metal atoms per gram
end
function composition(X, Y, Z, zmix = ZMIX)
    s = zmix.C + zmix.N + zmix.O
    xs = (X, Y, Z * zmix.C / s, Z * zmix.N / s, Z * zmix.O / s)
    ms = (mH, mHe, MC, MN, MO); zs = (1.0, 2.0, 6.0, 7.0, 8.0); χs = (χH, χHe1 + χHe2, XC, XN, XO)
    νion = 0.0; νe = 0.0; eion = 0.0; νZ = 0.0
    for i in 1:5
        ν = xs[i] / ms[i]; νion += ν; νe += zs[i] * ν; eion += χs[i] * ν; i > 2 && (νZ += ν)
    end
    Comp(X, Y, Z, xs, ms, zs, χs, νion, νe, eion, νZ)
end
"Single-species ion composition (mean A, Z; χ = 0), for comparison with one-ion-species EOS (Helmholtz)."
function composition_single(abar, zbar)
    m = abar * amu
    Comp(0.0, 0.0, 0.0, (1.0, 0, 0, 0, 0), (m, 1.0, 1.0, 1.0, 1.0), (zbar, 0, 0, 0, 0), (0.0, 0, 0, 0, 0),
         1 / m, zbar / m, 0.0, 0.0)
end
"Mean atomic weight Ā = 1/(Σ x_i/A_i) (with masses in units of amu) and Z̄ per ion of a composition."
abar(c::Comp) = 1 / (c.νion * amu)
"Electrons per baryon-ish: Y_e = ν_e amu (fully ionised)."
Ye(c::Comp) = c.νe * amu

# =====================================================================================================================
# Fully ionised branch (E)
# =====================================================================================================================
"Order of derivatives the number type must carry exactly: depth of Dual nesting."
dualdepth(::Type{<:Real}) = 0
dualdepth(::Type{<:ForwardDiff.Dual{Tg,V,N}}) where {Tg,V,N} = 1 + dualdepth(V)
npolish(::Type{T}) where T = (d = dualdepth(T); d <= 1 ? 1 : 2)   # k Newton steps from the exact root: exact through 2^k − 1

"Dimensionless electron/positron gas moments: N (n/A0), P (P/(A0 mc²)), E (kinetic E/(A0 mc²)) and dN/dη."
@inline function fgas(η, θ)
    F12, F32, F52, D12, D32, D52 = fd(η, θ)
    θ15 = θ * sqrt(θ); θ25 = θ15 * θ
    N = θ15 * (F12 + θ * F32)
    P = (2 / 3) * θ25 * (F32 + 0.5θ * F52)
    E = θ25 * (F32 + θ * F52)
    dN = θ15 * (D12 + θ * D32)
    return N, P, E, dN
end
const PAIRCUT = -46.0       # include positrons if their kinetic η₊ = −η − 2/θ exceeds this (e^{η₊} ≳ 1e-20)

"Net (n₋ − n₊)/A0 and derivative, Float64, for the η solve."
function netN(η::Float64, θ::Float64, pairs::Bool)
    N, _, _, dN = fgas(η, θ)
    if pairs
        ηp = -η - 2 / θ
        if ηp > PAIRCUT
            Np, _, _, dNp = fgas(ηp, θ)
            return N - Np, dN + dNp
        end
    end
    return N, dN
end

"Solve n_net(η, θ) = target (Float64), safeguarded Newton on ln n."
function solve_eta(target::Float64, θ::Float64, pairs::Bool)
    r = target / (A0 * θ^1.5 * 0.5 * sqrt(π))               # n/n_Q (nonrelativistic scale)
    ηd = (1.5 * target / (A0 * θ^1.5))^(2 / 3)
    η = r < 0.3 ? log(r) + r / (2 * sqrt(2.0)) : (r > 3.0 ? ηd : 0.5 * (ηd + log1p(r)))
    pairs && θ > 0.05 && r < 1e-3 && (η = max(η, -1 / θ + 1e-3))
    lo = -Inf; hi = Inf
    for it in 1:200
        n, dn = netN(η, θ, pairs)
        R = n > 0 ? log(n / (target / A0)) : -Inf
        R > 0 ? (hi = η) : (lo = η)
        abs(R) < 1e-14 && return η
        ηn = n > 0 ? η - R * n / dn : η + 1.0 + abs(η)
        if !(lo < ηn < hi)
            ηn = isfinite(lo) && isfinite(hi) ? 0.5 * (lo + hi) : (isfinite(hi) ? η - 1 - abs(η) : η + 1 + abs(η))
        end
        abs(ηn - η) < 1e-15 * max(1.0, abs(η)) && return ηn
        η = ηn
    end
    return η
end

# Yakovlev & Shalybkov (1989) OCP Coulomb free energy per ion in kT, φ(Γ), u = Γ φ' (the fits used in the Timmes–Swesty
# EOS). The two fits are joined C² over Γ ∈ [0.6, 1.7] in ln Γ; the Γ ≥ 1 branch is shifted by the constant needed to make
# φ continuous at Γ = 1 (0.0074; fixes its entropy offset).
const CA1 = -0.898004; const CB1 = 0.96786; const CC1 = 0.220703; const CD1 = -0.86097; const CE1 = 2.5269
const CA2 = 0.29561; const CB2 = 1.9885; const CC2 = 0.288675
const CKAP = (-2CC2 + CA2 / CB2) - (CA1 + 4CB1 + (-4CC1) - CE1)
smoother(x) = (x = clamp(x, zero(x), one(x)); x^3 * (10 - 15x + 6x^2))
dsmoother(x) = (x <= 0 || x >= 1) ? zero(x) : 30 * x^2 * (1 - x)^2
"φ(Γ) and u(Γ) = Γ dφ/dΓ."
function coulomb_phi(Γ)
    Γ4 = Γ^0.25
    φlo = -2CC2 * Γ * sqrt(Γ) + CA2 / CB2 * Γ^CB2
    ulo = -3CC2 * Γ * sqrt(Γ) + CA2 * Γ^CB2
    φhi = CA1 * Γ + 4CB1 * Γ4 - 4CC1 / Γ4 + CD1 * log(Γ) - CE1 + CKAP
    uhi = CA1 * Γ + CB1 * Γ4 + CC1 / Γ4 + CD1
    l = log(Γ); l0 = log(0.6); l1 = log(1.7)
    xx = (l - l0) / (l1 - l0); σ = smoother(xx); dσ = dsmoother(xx) / (l1 - l0)      # dσ/dlnΓ
    φ = (1 - σ) * φlo + σ * φhi
    u = (1 - σ) * ulo + σ * uhi + dσ * (φhi - φlo)
    return φ, u
end

"""
Gas p, e, s, n_e, η of the fully ionised branch (per cm³ / per gram). `pairs`: e± pairs; `coulomb`: YS89 correction.
Energy includes the ionisation energy of all ions and the pair-creation energy 2mc² n₊ (electron rest mass excluded).
"""
function branch_E(ρ, T, c::Comp; pairs::Bool = true, coulomb::Bool = true)
    kT = kB * T
    θ = kT / mc2
    ntarget = ρ * c.νe
    nt = Z0.pval(ntarget); θv = Z0.pval(θ)
    pairs = pairs && θv > 0.015
    η0 = solve_eta(nt, θv, pairs)
    Tn = typeof(ntarget * θ)
    η = η0 + zero(Tn)
    ns = npolish(Tn)
    for _ in 1:ns                                          # Newton polish in the caller's number type
        N, _, _, dN = fgas(η, θ)
        if pairs && Z0.pval(-η - 2 / θ) > PAIRCUT
            Np, _, _, dNp = fgas(-η - 2 / θ, θ)
            N = N - Np; dN = dN + dNp
        end
        η = η - (N - ntarget / A0) / dN
    end
    N, P, E, _ = fgas(η, θ)
    ηp = -η - 2 / θ
    haspos = pairs && Z0.pval(ηp) > PAIRCUT
    Pe = A0 * mc2 * P; Ek = A0 * mc2 * E
    Sk = Pe + Ek - kT * A0 * N * η                          # T s_e (electrons)
    Ee = Ek
    if haspos
        Np, Pp, Ep, _ = fgas(ηp, θ)
        Pe += A0 * mc2 * Pp
        Sk += A0 * mc2 * (Pp + Ep) - kT * A0 * Np * ηp
        Ee += A0 * mc2 * Ep + 2 * mc2 * A0 * Np
    end
    # ions: ideal, Sackur–Tetrode with g = 1
    pion = ρ * c.νion * kT
    eion_th = 1.5 * c.νion * kT
    sion = zero(Tn)
    lT = log(2π * kT / h^2)
    for i in 1:5
        c.x[i] > 0 || continue
        ν = c.x[i] / c.m[i]
        sion += ν * kB * (2.5 + 1.5 * (lT + log(c.m[i])) - log(ρ * ν))
    end
    p = Pe + pion; e = Ee / ρ + eion_th + c.eion; s = Sk / (ρ * T) + sion
    if coulomb
        zbar = c.νe / c.νion
        ni = ρ * c.νion
        a = cbrt(3 / (4π * ni))
        Γ = zbar^2 * qe^2 / (kT * a)
        φ, u = coulomb_phi(Γ)
        ν = c.νion
        pc = ni * kT * u / 3; ec = ν * kT * u; sc = -ν * kB * (φ - u)
        p += pc; e += ec; s += sc
    end
    return (p = p, e = e, s = s, ne = ntarget, η = η)
end

# =====================================================================================================================
# Saha branch (S)
# =====================================================================================================================
function branch_S(ρ, T, c::Comp)
    # (X = Y = 0, pure metals: Z0Physics has no electron source; treat as neutral atoms)
    s = c.X + c.Y > 0 ? Z0.saha(ρ, T; comp = Z0.Comp(c.X, c.Y)) :
        (ne = zero(ρ * T), nH0 = zero(ρ * T), nHp = zero(ρ * T), nHe0 = zero(ρ * T), nHe1 = zero(ρ * T), nHe2 = zero(ρ * T), nH = zero(ρ * T), nHe = zero(ρ * T), xH = zero(ρ * T))
    kT = kB * T
    nZ = ρ * c.νZ
    ntot = s.nH + s.nHe + s.ne + nZ
    p = ntot * kT
    e = (1.5 * ntot * kT + s.nHp * χH + s.nHe1 * χHe1 + s.nHe2 * (χHe1 + χHe2)) / ρ
    q(m) = (2π * m * kT / h^2)^1.5
    term(n, g, m) = Z0.pval(n) > 1e-280 ? n * (2.5 + log(g * q(m) / n)) : zero(n)
    S = term(s.nH0, 2, mH) + term(s.nHp, 1, mH) + term(s.nHe0, 1, mHe) + term(s.nHe1, 2, mHe) + term(s.nHe2, 1, mHe) +
        term(s.ne, 2, me)
    for i in 3:5
        c.x[i] > 0 && (S += term(ρ * c.x[i] / c.m[i], 1, c.m[i]))
    end
    sg = kB * S / ρ
    nQ = 2 * (2π * me * kT / h^2)^1.5
    return (p = p, e = e, s = sg, ne = s.ne, η = log(max(s.ne, 1e-300) / nQ))
end

# =====================================================================================================================
# Blend and API
# =====================================================================================================================
"Blend window: log₁₀T ∈ [T0, T0+ΔT] (Saha → fully ionised); log₁₀ρ ∈ [R0, R0+ΔR] (pressure ionisation)."
const BLEND = (T0 = 5.3, ΔT = 0.6, R0 = -1.0, ΔR = 1.7)
const LN10 = log(10.0)

"w and its partials ∂w/∂lnρ, ∂w/∂lnT (analytic)."
function blend_w(ρ, T, b = BLEND)
    xT = (log10(T) - b.T0) / b.ΔT; xR = (log10(ρ) - b.R0) / b.ΔR
    wT = smoother(xT); wR = smoother(xR)
    w = 1 - (1 - wT) * (1 - wR)
    dwdlnT = (1 - wR) * dsmoother(xT) / (b.ΔT * LN10)
    dwdlnρ = (1 - wT) * dsmoother(xR) / (b.ΔR * LN10)
    return w, dwdlnρ, dwdlnT
end

"""
    gas_state(ρ, T, comp; pairs, coulomb, blend) -> (p, e, s, ne, η)

Blended gas state per unit volume / mass. The blend is of the free energy: with Δf = f_E − f_S,
p = (1−w)p_S + w p_E + ρ Δf ∂w/∂lnρ,  s = (1−w)s_S + w s_E − (Δf/T) ∂w/∂lnT,  e = (1−w)e_S + w e_E − Δf ∂w/∂lnT.
"""
function gas_state(ρ, T, c::Comp; pairs::Bool = true, coulomb::Bool = true, blend = BLEND)
    w, wρ, wT = blend_w(ρ, T, blend)
    wv = Z0.pval(w)
    if wv <= 0
        return branch_S(ρ, T, c)
    elseif wv >= 1
        return branch_E(ρ, T, c; pairs, coulomb)
    end
    S = branch_S(ρ, T, c); E = branch_E(ρ, T, c; pairs, coulomb)
    Δf = (E.e - E.s * T) - (S.e - S.s * T)
    p = (1 - w) * S.p + w * E.p + ρ * Δf * wρ
    s = (1 - w) * S.s + w * E.s - Δf / T * wT
    e = (1 - w) * S.e + w * E.e - Δf * wT
    return (p = p, e = e, s = s, ne = (1 - w) * S.ne + w * E.ne, η = (1 - w) * S.η + w * E.η)
end

"Gas pressure and specific internal energy only (what `TitanMaterials.material_gas` needs)."
function gas_pe(ρ, T; X = 0.76, Y = 0.24, Z = 0.0, zmix = ZMIX, pairs = true, coulomb = true)
    g = gas_state(ρ, T, composition(X, Y, Z, zmix); pairs, coulomb)
    return g.p, g.e
end

Prad(T) = arad * T^4 / 3
erad(ρ, T) = arad * T^4 / ρ
srad(ρ, T) = 4arad * T^3 / (3ρ)

struct EOSTag end

"""
    eos(ρ, T; X, Y, Z, zmix, pairs, coulomb) -> NamedTuple

Gas-only: P, e, s, cV = ∂e/∂T|ρ, χρ = ∂lnP/∂lnρ|T, χT = ∂lnP/∂lnT|ρ, Γ1, ∇ad, δ = χT/χρ, cP, cs, ne (net electrons cm⁻³),
η (electron degeneracy parameter μ_kin/kT, rest mass excluded), plus dPdρ, dPdT, dedρ, cVs = T ∂s/∂T|ρ (for checks).
One evaluation with a two-partial Dual in (ln ρ, ln T) (unique tag, so it nests inside an outer ForwardDiff Jacobian).
"""
function eos(ρ, T; X = 0.76, Y = 0.24, Z = 0.0, zmix = ZMIX, pairs = true, coulomb = true, radiation = false,
             comp = composition(X, Y, Z, zmix))
    radiation || return thermo((r, t) -> gas_state(r, t, comp; pairs, coulomb), ρ, T)
    return thermo(ρ, T) do r, t
        g = gas_state(r, t, comp; pairs, coulomb)
        (p = g.p + Prad(t), e = g.e + erad(r, t), s = g.s + srad(r, t), ne = g.ne, η = g.η)
    end
end

"Thermodynamic derivatives of any state function `fstate(ρ, T) -> (p, e, s, ne, η)` (as `eos`)."
function thermo(fstate, ρ, T)
    V = promote_type(typeof(ρ), typeof(T), Float64)
    tag = ForwardDiff.Tag(EOSTag(), V)
    D = ForwardDiff.Dual{typeof(tag)}
    ρd = D(V(ρ), V(ρ), zero(V)); Td = D(V(T), zero(V), V(T))
    g = fstate(ρd, Td)
    val(x) = ForwardDiff.value(x); par(x, i) = ForwardDiff.partials(x, i)
    P = val(g.p); e = val(g.e); s = val(g.s)
    χρ = par(g.p, 1) / P; χT = par(g.p, 2) / P
    Sρ = par(g.s, 1); ST = par(g.s, 2)
    ∇ad = 1 / (χT - χρ * ST / Sρ); Γ1 = χρ - χT * Sρ / ST; δ = χT / χρ; cP = ST - Sρ * δ
    return (P = P, e = e, s = s, cV = par(g.e, 2) / T, χρ = χρ, χT = χT, Γ1 = Γ1, ∇ad = ∇ad, δ = δ, cP = cP,
            cs = sqrt(Γ1 * P / ρ), ne = val(g.ne), η = val(g.η), dPdρ = par(g.p, 1) / ρ, dPdT = par(g.p, 2) / T,
            dedρ = par(g.e, 1) / ρ, cVs = ST)
end

"""
    rho_PT(P, T; X, Y, Z) -> ρ

Density at given GAS pressure and temperature (Newton in ln ρ with the exact χρ, bracketed; then polished in the caller's
number type so that derivatives w.r.t. P, T are exact). NaN if no root in [1e-14, 1e5].
"""
function rho_PT(P, T; X = 0.76, Y = 0.24, Z = 0.0, zmix = ZMIX, kw...)
    Pf, Tf = Z0.pval(P), Z0.pval(T)
    c = composition(X, Y, Z, zmix)
    lo = log(1e-14); hi = log(1e5)
    lr = log(Pf / (kB * Tf) * 1.2 * mH)                       # ideal-gas guess
    lr = clamp(lr, lo, hi)
    f(x) = log(gas_state(exp(x), Tf, c; kw...).p) - log(Pf)
    for _ in 1:100
        fx = f(lr)
        fx > 0 ? (hi = lr) : (lo = lr)
        abs(fx) < 1e-14 && break
        dfd = ForwardDiff.derivative(f, lr)
        lrn = lr - fx / dfd
        lr = (lo < lrn < hi) ? lrn : 0.5 * (lo + hi)
    end
    abs(f(lr)) < 1e-10 || return oftype(P * T, NaN)
    if P isa ForwardDiff.Dual || T isa ForwardDiff.Dual
        g(x) = log(gas_state(exp(x), T, c; kw...).p) - log(P)
        x = lr + zero(g(lr))
        for _ in 1:npolish(typeof(x))
            x = x - g(x) / ForwardDiff.derivative(g, x)
        end
        lr = x
    end
    return exp(lr)
end

end # module
