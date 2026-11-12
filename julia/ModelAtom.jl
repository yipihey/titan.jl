"""
    ModelAtom

Non-LTE model atoms for the one radiation physics shared by the run and the snapshot spectrum (R2 of
docs/POPIII_RADIATION_PLAN.md). Not used by any run yet.

An `Atom` is a set of bound levels of one ion stage plus its continuum (the next ion's ground state): energies, statistical
weights, ionisation energies, bound–bound radiative data (A, f), bound–free cross sections σ_i(ν) and electron-collision rate
coefficients. `rates` builds the transition-rate matrix of one point; `solve` gives the statistical-equilibrium level
populations at fixed T, n_e, total density of the element in this ion pair, a continuum mean intensity J(ν) and line escape
probabilities β (Sobolev or case A/B). `line_emissivity` gives the escaping 4πj of a line.

Detailed balance is built in, so that with J = B_ν(T), β arbitrary and any collision data the populations are exactly LTE
(Saha–Boltzmann); this is the first test (`modelatom_test.jl`). Radiative recombination is the Milne inverse of each level's
photo-ionisation, with stimulated recombination in the given J; collisional de-excitation and three-body recombination by
detailed balance.

Hydrogenic atoms (`hydrogenic(Z; nmax)`, H I and He II), l-collapsed levels n = 1 … nmax (g = 2n²): exact absorption
oscillator strengths (Z0Physics.hyd_f, Gordon's formula for n ≤ 12), bound–free σ = Kramers × g_bf (Stobbe's exact ground
state, Menzel–Pekeris for n ≥ 2: Z0Physics.gaunt_bf), two-photon 2s → 1s (A = 8.2291 Z⁶ s⁻¹, Nussbaumer & Schmutz 1984) with
the 2s share 2/8 of n = 2 (l-mixed limit). The l-collapsed atom is exact where l-states are collisionally mixed (n_e ≳ 1e8–1e10
for the levels that matter), which is the envelope and infall regime; at nebular densities it misses the l-resolved cascade
(Storey & Hummer 1995) by several per cent. Collision rates: PLACEHOLDER van Regemorter (excitation) and Seaton (ionisation)
unless `coll = :chianti` (CHIANTI n = 1 → n′ effective collision strengths: Anderson et al. 2002 for H I,
Ballance for He II; the CHIANTI files have no transitions from 2p, so n ≥ 2 → n′ and all ionisation use Johnson 1972, read from
the ADS scan and checked: 1–2 within 1–2 % of Anderson et al.). `coll = :johnson`: Johnson for everything.

Tested (`modelatom_test.jl`, 12 Nov): J = B gives LTE to 1e-8 (round-off); α_A, α_B within 0.2–0.6 % of Storey & Hummer 1995 with
nmax = 200 (−2 % / −4 % for H I / He II at nmax = 40); case-B Hβ within 0–10 % of SH95 for n_e = 1e8–1e14 (Hα up to +22 % at
1e13, collision-sensitive); He II 4686 and 1640 are 13–16 % low even at 1e14: their n = 3, 4 levels are not l-mixed (He II
l-mixing needs ~Z⁶ higher density), so the l-collapsed cascade misses the 4f → 3d channel. To do: l-resolved low levels with
proton l-mixing (Pengelly & Seaton 1964); the He I atom from the CHIANTI he_1 files.
"""
module ModelAtom

import ..Z0Physics
using LinearAlgebra
const Z0 = Z0Physics
const h = Z0.h; const kB = Z0.kB; const c = Z0.c; const eV = Z0.eV; const me = Z0.me
const ecgs = 4.80320471e-10
const ACOEF = 8π^2 * ecgs^2 / (me * c^3)         # A_ul = ACOEF ν² (g_l/g_u) f_lu  [7.422e-22 cgs]
const BCOEF = π * ecgs^2 / (me * c)               # πe²/(m c): ∫σ_ν dν = BCOEF f  [0.02654 cm² Hz]

struct Level
    E::Float64          # excitation energy above the ground state of this ion [erg]
    g::Float64
    χ::Float64          # ionisation energy from this level to the next ion's ground state [erg]
    label::String
end
struct Atom
    name::String
    levels::Vector{Level}
    gplus::Float64                          # statistical weight of the next ion's ground state
    f::Matrix{Float64}                      # absorption oscillator strengths f[l, u] (l < u), 0 = no radiative line
    A::Matrix{Float64}                      # spontaneous rates A[u, l]
    A2::Vector{Tuple{Int,Int,Float64}}      # two-photon decays (u, l, A)
    σ::Vector{Function}                     # σ_i(ν) [cm²], 0 below threshold
    qexc::Function                          # (l, u, T) -> upward collisional rate coefficient [cm³ s⁻¹]
    qion::Function                          # (i, T) -> collisional ionisation rate coefficient [cm³ s⁻¹]
    coll::Symbol                            # provenance tag of the collision data (:approx placeholder, :chianti)
end
nlev(a::Atom) = length(a.levels)
νline(a::Atom, l, u) = (a.levels[u].E - a.levels[l].E) / h

"n_i*/(n_e n_+): Saha–Boltzmann factor of level i [cm³]."
Φ(a::Atom, i, T) = a.levels[i].g / (2a.gplus) * exp(a.levels[i].χ / (kB * T)) / Z0.λe3(T)

# ---------------------------------------------------------------- hydrogenic atoms
"e^y E₁(y) (Abramowitz & Stegun 5.1.53 for y ≤ 1, 5.1.56 above)."
function expE1(y)
    if y <= 1
        E1 = -log(y) - 0.57721566 + 0.99999193y - 0.24991055y^2 + 0.05519968y^3 - 0.00976004y^4 + 0.00107857y^5
        return exp(y) * E1
    end
    (y^2 + 2.334733y + 0.250621) / (y^2 + 3.330657y + 1.681534) / y
end

"""
    hydrogenic(Z; nmax = 20, χ1 = Z == 1 ? Z0.χH : Z0.χHe2)

H I (Z = 1) or He II (Z = 2), l-collapsed n = 1 … nmax. Level energies E_n = χ1 (1 − 1/n²) with the ion's own (reduced-mass)
ground-state ionisation energy χ1.
"""
function hydrogenic(Z::Int; nmax = 20, χ1 = Z == 1 ? Z0.χH : Z0.χHe2, coll = :approx, chianti = nothing)
    lv = [Level(χ1 * (1 - 1 / n^2), 2.0n^2, χ1 / n^2, "n=$n") for n in 1:nmax]
    N = nmax; f = zeros(N, N); A = zeros(N, N)
    for l in 1:N-1, u in l+1:N
        f[l, u] = Z0.hyd_f(l, u)                                   # f is independent of Z (hydrogenic scaling)
        ν = (lv[u].E - lv[l].E) / h
        A[u, l] = ACOEF * ν^2 * (lv[l].g / lv[u].g) * f[l, u]
    end
    A2 = N >= 2 ? [(2, 1, 8.2291 * Z^6 * 2 / 8)] : Tuple{Int,Int,Float64}[]
    νth = [lv[n].χ / h for n in 1:N]
    σ = Function[(ν -> ν < νth[n] ? 0.0 : Z0.KRAMERS * Z^4 / (n^5 * ν^3) * Z0.gaunt_bf(n, ν / νth[n])) for n in 1:N]
    # PLACEHOLDER collision data (to be replaced by CHIANTI): van Regemorter with P = 0.276 e^y E₁(y) (neutral), floor 0.2 (ions);
    # Seaton ionisation with ḡ = 0.1, 0.2, 0.3 for n = 1, 2, ≥ 3.
    function qexc(l, u, T)
        f[l, u] > 0 || return 0.0
        ΔE = lv[u].E - lv[l].E; y = ΔE / (kB * T)
        P = Z == 1 ? 0.276 * expE1(y) : max(0.2, 0.276 * expE1(y))
        1.2521e-4 / sqrt(T) * (Z0.χH / ΔE) * f[l, u] * P * exp(-y)          # 8π/√3 · (2π a0² … ) form as in NLTEH.qexc
    end
    function qion(i, T)
        u0 = lv[i].χ / (kB * T); gbar = i == 1 ? 0.1 : i == 2 ? 0.2 : 0.3
        1.55e13 / sqrt(T) * gbar * σ[i](νth[i] * 1.000001) * exp(-u0) / u0
    end
    if coll === :johnson || coll === :chianti
        # Johnson 1972 for every pair and every ionisation; CHIANTI's n-collapsed Υ replaces the pairs it contains
        cup = coll === :chianti ? chianti_collapsed(chianti)[1] : nothing
        qexc2(l, u, T) = (Υ = cup === nothing ? nothing : cup(l, u, T);
                          Υ === nothing ? johnson_exc(l, u, T; Z) : q_from_ups(Υ, lv[l].g, lv[u].E - lv[l].E, T))
        qion2(i, T) = johnson_ion(i, T; Z)
        return Atom(Z == 1 ? "H I" : "He II", lv, 1.0, f, A, A2, σ, qexc2, qion2, coll)
    end
    Atom(Z == 1 ? "H I" : "He II", lv, 1.0, f, A, A2, σ, qexc, qion, :approx)
end

# ---------------------------------------------------------------- collision data
"""
Vriens & Smeets (1980, Phys. Rev. A 22, 940) electron-impact rate coefficients for hydrogen, excitation p → q and ionisation
from p [cm³ s⁻¹], kT and energies in eV. TRANSCRIBED FROM MEMORY of the widely reproduced fit formulas; checked here only
against the CHIANTI (Anderson et al. 2002) rates for n ≤ 5 (`modelatom_test.jl`), not against the paper. For hydrogenic ions
the rates are scaled as K_Z(T) = K_H(T/Z²)/Z³.
"""
function vs_exc(p::Int, q::Int, T; Z = 1, f = Z0.hyd_f(p, q))   # NOT USED: disagrees with Anderson et al. by 0.2–3×, transcription unverified
    Ry = 13.605693; kT = kB * (T / Z^2) / eV
    Ep = Ry / p^2; Epq = Ry * (1 / p^2 - 1 / q^2)
    A = 2Ry / Epq * f
    bp = 1.4 * log(p) / p - 0.7 / p - 0.51 / p^2 + 1.16 / p^3 - 0.55 / p^4
    B = 4Ry^2 / q^3 * (1 / Epq^2 + 4Ep / (3Epq^3) + bp * Ep^2 / Epq^4)
    s = abs(q - p); Δ = exp(-B / A) + 0.06 * s^2 / (p * q^2)
    Γ = Ry * log(1 + p^3 * kT / Ry) * (3 + 11 * (s / p)^2) / (6 + 1.6q * s + 0.3 / s^2 + 0.8 * q^1.5 / sqrt(s) * abs(s - 0.6))
    1.6e-7 * sqrt(kT) / (kT + Γ) * exp(-Epq / kT) * (A * log(0.3kT / Ry + Δ) + B) / Z^3
end
function vs_ion(p::Int, T; Z = 1)
    kT = kB * (T / Z^2) / eV; Ep = 13.605693 / p^2
    9.56e-6 * kT^-1.5 * exp(-Ep / kT) / (Ep^2.33 + 4.38 * Ep^1.72 * kT + 1.32 * Ep * kT^2) / Z^3
end

"Exponential integrals E₀ … E₂ at t > 0 (E₀ = e^{−t}/t; E₁ from Abramowitz & Stegun via expE1; E₂ by recursion, Johnson eq. 9)."
function expints(t)
    e = exp(-t); E1 = expE1(t) * e
    (e / t, E1, e - t * E1)
end
"""
Johnson (1972, ApJ 174, 227) electron-impact rate coefficients for hydrogen [cm³ s⁻¹], read from the ADS scan:
excitation n → n′ (eq. 36 with eqs. 11, 23, 25–26, 30–32, 37–38) and ionisation from n (eq. 39 with 20, 24, 40–42).
Oscillator strengths: the exact hydrogenic f (Johnson uses his eq. 3 fit). Hydrogenic ions scaled as K_Z(T) = K_H(T/Z²)/Z³
(approximate).
"""
function johnson_exc(n::Int, np::Int, T; Z = 1)
    Tz = T / Z^2; kT = kB * Tz; x = 1 - (n / np)^2; Enn = Z0.χH * x / n^2
    rn = n == 1 ? 0.45 : 1.94 * n^-1.57; bn = n == 1 ? -0.603 : (4.0 - 18.63 / n + 36.24 / n^2 - 28.09 / n^3) / n
    A = 2n^2 / x * Z0.hyd_f(n, np)
    B = 4n^4 / np^3 / x^2 * (1 + 4 / (3x) + bn / x^2)
    y = Enn / kT; z = rn * x + y
    _, E1y, E2y = expints(y); _, E1z, E2z = expints(z)
    pref = sqrt(8kT / (π * me)) * 2n^2 / x * π * 5.29177211e-9^2 * y^2
    pref * (A * ((1 / y + 0.5) * E1y - (1 / z + 0.5) * E1z) + (B - A * log(2n^2 / x)) * (E2y / y - E2z / z)) / Z^3
end
function johnson_ion(n::Int, T; Z = 1)
    Tz = T / Z^2; kT = kB * Tz; In = Z0.χH / n^2
    rn = n == 1 ? 0.45 : 1.94 * n^-1.57; bn = n == 1 ? -0.603 : (4.0 - 18.63 / n + 36.24 / n^2 - 28.09 / n^3) / n
    g = n == 1 ? (1.1330, -0.4059, 0.07014) : n == 2 ? (1.0785, -0.2319, 0.02947) :
        (0.9935 + 0.2328 / n - 0.1296 / n^2, -(0.6282 - 0.5598 / n + 0.5299 / n^2) / n, (0.3887 - 1.181 / n + 1.470 / n^2) / n^2)
    An = 32 / (3 * sqrt(3) * π) * n * (g[1] / 3 + g[2] / 4 + g[3] / 5)
    Bn = 2 / 3 * n^2 * (5 + bn)
    y = In / kT; z = rn + y
    E0y, E1y, E2y = expints(y); E0z, E1z, E2z = expints(z)
    ξ(E0, E1, E2) = E0 - 2E1 + E2
    pref = sqrt(8kT / (π * me)) * 2n^2 * π * 5.29177211e-9^2 * y^2
    pref * (An * (E1y / y - E1z / z) + (Bn - An * log(2n^2)) * (ξ(E0y, E1y, E2y) - ξ(E0z, E1z, E2z))) / Z^3
end

"""
CHIANTI level (.elvlc) and effective collision strength (.scups, Burgess & Tully 1992 scaled spline) files. `read_elvlc`
returns (index, configuration, 2S+1, L, J, E [cm⁻¹], observed if given else theoretical); `read_scups` the transitions
(i, j, ΔE [Ry], gf, high-T limit, ttype, C, x knots, y knots); `ups(t, T)` the effective collision strength Υ(T).
"""
function read_elvlc(path)
    out = []
    for ln in eachline(path)
        strip(ln) == "-1" && break
        p = split(ln); length(p) < 7 && continue
        i = parse(Int, p[1]); Eo = parse(Float64, p[end-1]); Et = parse(Float64, p[end])
        push!(out, (i = i, conf = p[2], mult = parse(Int, p[end-4]), L = p[end-3], J = parse(Float64, p[end-2]), E = Eo > 0 || i == 1 ? Eo : Et))
    end
    out
end
const NUMRE = r"-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+"
function read_scups(path)
    L = readlines(path); out = []; k = 1
    while k + 2 <= length(L)
        strip(L[k]) == "-1" && break
        h_ = [m.match for m in eachmatch(NUMRE, L[k])]
        i, j = parse(Int, h_[1]), parse(Int, h_[2]); dE, gf, lim = parse.(Float64, h_[3:5]); nspl, tt = parse(Int, h_[6]), parse(Int, h_[7]); C = parse(Float64, h_[8])
        x = parse.(Float64, [m.match for m in eachmatch(NUMRE, L[k+1])]); y = parse.(Float64, [m.match for m in eachmatch(NUMRE, L[k+2])])
        length(x) == nspl && length(y) == nspl || error("read_scups: bad record at line $k of $path")
        push!(out, (; i, j, dE, gf, lim, tt, C, x, y)); k += 3
    end
    out
end
"Natural cubic spline through (x_k, y_k), evaluated at t."
function nspline(x, y, t)
    n = length(x); h_ = diff(x); a = zeros(n); b = zeros(n); cc = zeros(n); d = zeros(n)
    # second derivatives M with M_1 = M_n = 0 (tridiagonal)
    for k in 2:n-1; a[k] = h_[k-1]; b[k] = 2(h_[k-1] + h_[k]); cc[k] = h_[k]; d[k] = 6((y[k+1] - y[k]) / h_[k] - (y[k] - y[k-1]) / h_[k-1]); end
    M = zeros(n)
    if n > 2
        T_ = Tridiagonal(a[3:n-1], b[2:n-1], cc[2:n-2]); M[2:n-1] = T_ \ d[2:n-1]
    end
    k = clamp(searchsortedlast(x, t), 1, n - 1); hk = h_[k]; u = (x[k+1] - t) / hk; w = (t - x[k]) / hk
    M[k] * (x[k+1] - t)^3 / (6hk) + M[k+1] * (t - x[k])^3 / (6hk) + (y[k] - M[k] * hk^2 / 6) * u + (y[k+1] - M[k+1] * hk^2 / 6) * w
end
function ups(t, T)
    kt = T / (157887.5 * t.dE); C = t.C
    if t.tt == 1
        x = 1 - log(C) / log(kt + C); return nspline(t.x, t.y, x) * log(kt + ℯ)
    elseif t.tt == 2
        x = kt / (kt + C); return nspline(t.x, t.y, x)
    elseif t.tt == 3
        x = kt / (kt + C); return nspline(t.x, t.y, x) / (kt + 1)
    elseif t.tt == 4
        x = 1 - log(C) / log(kt + C); return nspline(t.x, t.y, x) * log(kt + C)
    end
    error("ups: transition type $(t.tt) not implemented")
end
"Principal quantum number of a CHIANTI configuration string (the last shell: '3p' -> 3, '1s.2p' -> 2)."
nprin(conf) = parse(Int, match(r"(\d+)[spdfghi]\d*$", conf).captures[1])

"""
n-collapsed effective collision strengths Υ(1, n′; T) of a hydrogenic ion from CHIANTI files `stem`.elvlc/.scups: sums over the
fine-structure and l levels of the upper shell. Only n = 1 → n′ (the files have no transitions from 2p). Returns a function (p, q, T) -> Υ, or nothing where the pair is not in the file.
"""
function chianti_collapsed(stem)
    lv = read_elvlc(stem * ".elvlc"); nof = Dict(l.i => nprin(l.conf) for l in lv)
    tr = read_scups(stem * ".scups"); pairs = Dict{Tuple{Int,Int},Vector{Any}}()
    for t in tr
        p, q = nof[t.i], nof[t.j]; p == q && continue
        min(p, q) == 1 || continue          # the H I and He II files hold transitions from 1s and 2s only (none from 2p):
                                            # an n = 2 → n′ sum would be incomplete, so only n = 1 → n′ (complete) is used
        push!(get!(pairs, (min(p, q), max(p, q)), Any[]), t)
    end
    (p, q, T) -> (k = (min(p, q), max(p, q)); haskey(pairs, k) ? sum(ups(t, T) for t in pairs[k]) : nothing), maximum(values(nof))
end
"Upward rate coefficient from an effective collision strength: q_lu = 8.629e-6 Υ e^{−ΔE/kT}/(g_l √T) [cm³ s⁻¹]."
q_from_ups(Υ, gl, ΔE, T) = 8.629e-6 * Υ * exp(-ΔE / (kB * T)) / (gl * sqrt(T))

# ---------------------------------------------------------------- rates
"""
Frequency nodes (midpoints in ln ν) and weights dν for the continuum of level i at T: from the threshold ν0 to
ν0·max(50, 1 + 80 kT/χ), `nν` points. One grid serves the photo-ionisation rate and the recombination coefficient, so that
J = B_ν(T) gives detailed balance to round-off.
"""
function bf_grid(a::Atom, i, T; nν = 1500)
    χ = a.levels[i].χ; ν0 = χ / h; L = log(max(50.0, 1 + 80kB * T / χ)); dl = L / nν
    ν = [ν0 * exp((k - 0.5) * dl) for k in 1:nν]
    (ν, ν .* dl)
end

"""
Photo-ionisation rate R_ic = 4π∫σJ/(hν)dν and recombination coefficient α_i = (n_i*/n_e n_+) 4π∫σ/(hν)(2hν³/c² + J)e^{−hν/kT}dν
(spontaneous + stimulated, per n_e n_+) of level i in the mean intensity J(ν) (J = nothing: no radiation).
"""
function bf_rates(a::Atom, i, T; J = nothing)
    ν, w = bf_grid(a, i, T); kT = kB * T; χ = a.levels[i].χ
    pre = a.levels[i].g / (2a.gplus) / Z0.λe3(T)           # Φ_i e^{−χ/kT}: the exponentials combine without overflow
    α = 0.0; R = 0.0
    @inbounds for k in eachindex(ν)
        s = a.σ[i](ν[k]); s == 0 && continue
        Jv = J === nothing ? 0.0 : J(ν[k])
        α += w[k] * 4π * s / (h * ν[k]) * (2h * ν[k]^3 / c^2 + Jv) * exp(-(h * ν[k] - χ) / kT)
        R += w[k] * 4π * s / (h * ν[k]) * Jv
    end
    (R, pre * α)
end

"""
    rates(a, T, ne; J, Jline, β) -> P

Rate matrix P[i, j] (s⁻¹, i → j; index nlev+1 = continuum) at temperature T and electron density ne.
`J(ν)`: continuum mean intensity (nothing = 0). `β(l, u)`: line escape probability; the line rates are the escape-probability
(Sobolev) form with the local continuum J_c = J(ν_lu): net rate n_u A β − β (n_l B_lu − n_u B_ul) J_c, i.e.
P[u, l] = β (A + B_ul J_c), P[l, u] = β B_lu J_c. β = 0 makes the line exactly thick (case B for the Lyman series).
"""
function rates(a::Atom, T, ne; J = nothing, β = (l, u) -> 1.0)
    N = nlev(a); P = zeros(N + 1, N + 1)
    for l in 1:N-1, u in l+1:N
        A = a.A[u, l]
        if A > 0
            b = β(l, u); ν = νline(a, l, u)
            Jc = J === nothing ? 0.0 : J(ν)
            Bul = A * c^2 / (2h * ν^3); Blu = Bul * a.levels[u].g / a.levels[l].g
            P[u, l] += b * (A + Bul * Jc); P[l, u] += b * Blu * Jc
        end
        qu = a.qexc(l, u, T)
        if qu > 0
            P[l, u] += ne * qu
            P[u, l] += ne * qu * (a.levels[l].g / a.levels[u].g) * exp((a.levels[u].E - a.levels[l].E) / (kB * T))
        end
    end
    for (u, l, A2) in a.A2
        P[u, l] += A2            # (two-photon: no stimulated term; its inverse, two-photon absorption, is left out)
    end
    for i in 1:N
        R, α = bf_rates(a, i, T; J)
        qi = a.qion(i, T)
        P[i, N+1] += R + ne * qi
        P[N+1, i] += ne * α + ne * ne * qi * Φ(a, i, T)     # radiative (per n_+) and three-body recombination: rates per ion
    end
    P
end

"""
    solve(a, T, ne, ntot; J, β) -> (n = level populations, nplus = next-ion density)

Statistical equilibrium at fixed T and n_e: Σ_j n_j P_ji = n_i Σ_j P_ij for the bound levels and the continuum, with
Σ n_i + n_+ = ntot replacing the continuum equation.
"""
function solve(a::Atom, T, ne, ntot; kw...)
    P = rates(a, T, ne; kw...); N = nlev(a) + 1
    M = P' - Diagonal(vec(sum(P; dims = 2)))      # dn/dt = M n
    M[N, :] .= 1.0; rhs = zeros(N); rhs[N] = ntot
    x = M \ rhs
    (n = x[1:end-1], nplus = x[end], P)
end

"""
    solve_upper(a, T, ne, nplus; nfix = [0.0], kw...) -> n

Populations of the levels above the k = length(nfix) lowest ones, whose densities are given by `nfix`, for a given ion
density n_+ (the recombination-line
problem of a photo-ionised gas, where the ionisation balance is set elsewhere): Σ_j n_j P_ji + n_+ P_ci + n1 P_1i = n_i Σ_j P_ij.
"""
function solve_upper(a::Atom, T, ne, nplus; nfix = [0.0], kw...)
    P = rates(a, T, ne; kw...); N = nlev(a); k = length(nfix)
    idx = k+1:N; M = zeros(N - k, N - k); rhs = zeros(N - k)
    for (r, i) in enumerate(idx)
        M[r, r] -= sum(P[i, :])
        for (q, j) in enumerate(idx); j == i && continue; M[r, q] += P[j, i]; end
        rhs[r] = -(nplus * P[N+1, i] + sum(nfix[j] * P[j, i] for j in 1:k))
    end
    vcat(nfix, M \ rhs)
end

"LTE populations (Saha–Boltzmann) at T, ne for the same total."
function lte(a::Atom, T, ne, ntot)
    w = [Φ(a, i, T) * ne for i in 1:nlev(a)]      # n_i / n_+
    nplus = ntot / (1 + sum(w)); (n = w .* nplus, nplus)
end

"Escaping line emissivity 4πj = n_u A_ul β h ν_lu [erg cm⁻³ s⁻¹] of the line l–u."
line_emissivity(a::Atom, n, l, u; β = 1.0) = n[u] * a.A[u, l] * β * h * νline(a, l, u)

"Planck function B_ν(T)."
Bν(ν, T) = 2h * ν^3 / c^2 / expm1(h * ν / (kB * T))

end # module
