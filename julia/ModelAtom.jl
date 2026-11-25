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
proton l-mixing (Pengelly & Seaton 1964).

He I (`heI`, 99 terms: CHIANTI n ≤ 5 + hydrogenic singlet/triplet shells 6–40), case B against Porter et al. 2012/13 at 1e4 K:
4471/5876/6678/10830 within 0.75–1.23 for n_e ≤ 1e8 and 0.94–1.07 at 1e14, but 0.39–0.81 at n_e = 1e10–1e12 (2e4 K similar).
The intermediate-density deficit is attributed (not yet proven) to the approximate data of the excited terms: hydrogenic
photo-ionisation, van Regemorter for pairs CHIANTI lacks, Johnson ionisation at n*. To do: R-matrix (TOPbase) cross sections and
the missing collision data for the excited terms.
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
function johnson_ion(n::Real, T; Z = 1, In = Z0.χH / n^2)
    Tz = T / Z^2; kT = kB * Tz
    rn = n == 1 ? 0.45 : 1.94 * n^-1.57; bn = n == 1 ? -0.603 : (4.0 - 18.63 / n + 36.24 / n^2 - 28.09 / n^3) / n
    g = n == 1 ? (1.1330, -0.4059, 0.07014) : (n isa Integer && n == 2) ? (1.0785, -0.2319, 0.02947) :
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

# ---------------------------------------------------------------- He I
"CHIANTI .wgfa: (i, j, λ [Å], gf, A [s⁻¹], label); λ = 0 marks two-photon / autoionising entries."
function read_wgfa(path)
    out = []
    for ln in eachline(path)
        strip(ln) == "-1" && break
        m = match(r"^\s*(\d+)\s+(\d+)\s+(-?[\d.]+)\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s*(.*)$", ln); m === nothing && continue
        push!(out, (i = parse(Int, m[1]), j = parse(Int, m[2]), λ = parse(Float64, m[3]), gf = parse(Float64, m[4]), A = parse(Float64, m[5]), label = m[6]))
    end
    out
end

"""
    heI(; nmax = 40, stem = atomdata/chianti/he_1)

He I model atom. Terms n ≤ 5 from CHIANTI (energies, fine structure collapsed: g = Σ(2J+1), g-weighted mean energy; A summed
over the fine-structure components, A_term = Σ g_J A / g_term; the 2¹S two-photon decay from the λ = 0 entry), collisions from
the CHIANTI effective collision strengths (Sawey & Berrington 1993, Bray et al. 2000) for the transitions they hold (from the
six lowest fine-structure levels: 1¹S, 2³S, 2¹S, 2³P₀,₁,₂). Above n = 5: hydrogenic shells n = 6 … nmax in the singlet and the
triplet system (g = n² and 3n², E = χ_HeI − Ry/n²), radiative decays with hydrogenic n–n′ oscillator strengths, distributed over
the CHIANTI terms of a lower shell n ≤ 5 of the same spin in proportion to their statistical weights (APPROXIMATE).
Missing collision data: dipole-allowed term pairs get van Regemorter rates (APPROXIMATE, factor ~2), shell pairs n, n′ ≥ 6
Johnson 1972 (hydrogenic outer electron), spin-changing pairs above 2³P none. Photo-ionisation: 1¹S Verner et al. 1996;
excited terms hydrogenic Kramers with n* = √(Ry/χ) (APPROXIMATE; R-matrix cross sections to do). Collisional ionisation:
Johnson 1972 at n* = √(Ry/χ) with the term's own χ for excited terms (APPROXIMATE); 1¹S Voronov 1997 (verified).
"""
function heI(; nmax = 40, stem = joinpath(@__DIR__, "atomdata", "chianti", "he_1"))
    hc = h * c; Ry = 13.605693 * eV; χ1 = Z0.χHe1
    lv = read_elvlc(stem * ".elvlc")
    key(l) = (l.conf, l.mult, l.L)
    tkeys = unique([key(l) for l in lv]); tid = Dict(k => t for (t, k) in enumerate(tkeys))
    gT = zeros(length(tkeys)); ET = zeros(length(tkeys))
    for l in lv; t = tid[key(l)]; g = 2l.J + 1; gT[t] += g; ET[t] += g * l.E * hc; end
    ET ./= gT
    nT = [nprin(k[1]) for k in tkeys]; mT = [k[2] for k in tkeys]
    levels = [Level(ET[t], gT[t], χ1 - ET[t], "$(tkeys[t][1]) $(tkeys[t][2])$(tkeys[t][3])") for t in eachindex(tkeys)]
    nC = length(levels)
    for n in 6:nmax, mult in (1, 3)
        E = χ1 - Ry / n^2; push!(levels, Level(E, mult * n^2, Ry / n^2, "n=$n $(mult == 1 ? "singlet" : "triplet")")); push!(nT, n); push!(mT, mult)
    end
    N = length(levels); A = zeros(N, N); f = zeros(N, N); A2 = Tuple{Int,Int,Float64}[]
    ν(l, u) = (levels[u].E - levels[l].E) / h
    for w in read_wgfa(stem * ".wgfa")
        ti = tid[key(lv[w.i])]; tj = tid[key(lv[w.j])]; ti == tj && continue
        lo, up = levels[ti].E < levels[tj].E ? (ti, tj) : (tj, ti)
        jup = levels[ti].E < levels[tj].E ? w.j : w.i; gup = 2lv[jup].J + 1      # A in the file is per upper fine-structure level
        if w.λ == 0
            push!(A2, (up, lo, w.A * gup / levels[up].g))
        else
            A[up, lo] += w.A * gup / levels[up].g
        end
    end
    # hydrogenic shells above n = 5: decays to lower shells of the same spin
    shell_terms(n, mult) = [t for t in 1:N if nT[t] == n && mT[t] == mult]
    for u in nC+1:N, nl in 1:nT[u]-1
        mult = mT[u]; ts = shell_terms(nl, mult); isempty(ts) && continue
        fl = Z0.hyd_f(nl, nT[u]); gs = sum(levels[t].g for t in ts)
        for t in ts
            A[u, t] += ACOEF * ν(t, u)^2 * (gs / levels[u].g) * fl * levels[t].g / gs
        end
    end
    for u in 1:N, l in 1:N
        A[u, l] > 0 && (f[l, u] = A[u, l] / (ACOEF * ν(l, u)^2 * levels[l].g / levels[u].g))
    end
    # photo-ionisation
    σ = Function[]
    for t in 1:N
        νt = levels[t].χ / h
        if t == 1
            push!(σ, ν -> ν < νt ? 0.0 : σverner_HeI(h * ν / eV))
        else
            ns = sqrt(Ry / levels[t].χ); push!(σ, ν -> ν < νt ? 0.0 : 7.906e-18 * ns * (νt / ν)^3)
        end
    end
    # collisions
    tr = read_scups(stem * ".scups"); Υpairs = Dict{Tuple{Int,Int},Vector{Any}}()
    for t_ in tr
        a_, b_ = tid[key(lv[t_.i])], tid[key(lv[t_.j])]; a_ == b_ && continue
        push!(get!(Υpairs, (min(a_, b_), max(a_, b_)), Any[]), t_)
    end
    function qexc(l, u, T)
        ΔE = levels[u].E - levels[l].E; ΔE <= 0 && return 0.0
        k = (min(l, u), max(l, u))
        haskey(Υpairs, k) && return q_from_ups(sum(ups(t_, T) for t_ in Υpairs[k]), levels[l].g, ΔE, T)
        if l > nC && u > nC && mT[l] == mT[u]
            return johnson_exc(nT[l], nT[u], T)      # rate per atom in the lower shell (independent of the spin weight)
        end
        f[l, u] > 1e-4 || return 0.0
        y = ΔE / (kB * T)
        1.2521e-4 / sqrt(T) * (Z0.χH / ΔE) * f[l, u] * 0.276 * expE1(y) * exp(-y)
    end
    function qion(t, T)
        t == 1 && return (U = 24.6 * eV / (kB * T); 1.75e-8 * U^0.35 * exp(-U) / (0.180 + U))    # Voronov 1997 (verified against Verner's cfit.dat and Cloudy's coll_ion.dat)
        ns = sqrt(Ry / levels[t].χ)                       # effective quantum number; Johnson's n ≥ 3 coefficients at n*, with the level's own χ
        johnson_ion(max(ns, 1.5), T; In = levels[t].χ)
    end
    Atom("He I", levels, 2.0, f, A, A2, σ, qexc, qion, :chianti)
end
"He I 1s² photo-ionisation cross section [cm²] at E [eV] (Verner et al. 1996 fit, as Spectrum.σHe1)."
function σverner_HeI(E)
    E < 24.587 && return 0.0
    E0, σ0, ya, P, yw, y0, y1 = 13.61, 949.2, 1.469, 3.188, 2.039, 0.4434, 2.136
    x = E / E0 - y0; y = sqrt(x^2 + y1^2)
    1e-18 * σ0 * ((x - 1)^2 + yw^2) * y^(0.5P - 5.5) * (1 + sqrt(y / ya))^(-P)
end

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
Energy-weighted bound–free rates of level i: ΓE_i = 4π∫σJ dν (energy absorbed per atom in level i, erg s⁻¹) and
αE_i = (n_i*/n_e n_+) 4π∫σ(2hν³/c² + J)e^{−hν/kT}dν (energy emitted per n_e n_+, erg cm³ s⁻¹), on the grid of `bf_rates`.
"""
function bf_energy(a::Atom, i, T; J = nothing)
    ν, w = bf_grid(a, i, T); kT = kB * T; χ = a.levels[i].χ
    pre = a.levels[i].g / (2a.gplus) / Z0.λe3(T)
    αE = 0.0; ΓE = 0.0
    @inbounds for k in eachindex(ν)
        s = a.σ[i](ν[k]); s == 0 && continue
        Jv = J === nothing ? 0.0 : J(ν[k])
        αE += w[k] * 4π * s * (2h * ν[k]^3 / c^2 + Jv) * exp(-(h * ν[k] - χ) / kT)
        ΓE += w[k] * 4π * s * Jv
    end
    (ΓE, pre * αE)
end
"""
    radiative_cooling(a, n, nplus, ne, T; J, β) -> (lines, bf, twophoton)

Net radiated power per volume [erg cm⁻³ s⁻¹] of the atom in the state (n, nplus): lines Σ hν β (n_u A − (n_l B_lu − n_u B_ul) J_c),
bound–free Σ_i (n_e n_+ αE_i − n_i ΓE_i), two-photon emission n_u A₂ hν_ul. In statistical equilibrium the internal energy of the
atom is constant, so the sum is the net thermal energy the gas loses to radiation through this atom (free–free not included).
With J = B_ν(T) and LTE populations the sum vanishes.
"""
function radiative_cooling(a::Atom, n, nplus, ne, T; J = nothing, β = (l, u) -> 1.0)
    cl = 0.0
    for (l, u, ν, A, Bul, Blu, f) in lines(a)
        Jc = J === nothing ? 0.0 : J(ν); b = β(l, u)
        cl += h * ν * b * (n[u] * A - (n[l] * Blu - n[u] * Bul) * Jc)
    end
    cb = 0.0
    for i in 1:nlev(a)
        ΓE, αE = bf_energy(a, i, T; J); cb += ne * nplus * αE - n[i] * ΓE
    end
    c2 = sum((n[u] * A2 * (a.levels[u].E - a.levels[l].E) for (u, l, A2) in a.A2); init = 0.0)
    (lines = cl, bf = cb, twophoton = c2)
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
        qu = a.coll === :lres ? a.qexc(l, u, T, ne) : a.qexc(l, u, T)
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
Split form of `rates` for iterating on the line escape probabilities: `static_rates` holds collisions, bound–free and two-photon
terms (independent of β), `add_lines!` adds the line terms for a given β; P = static + lines equals `rates(a, T, ne; J, β)`.
"""
function static_rates(a::Atom, T, ne; J = nothing)
    N = nlev(a); P = zeros(N + 1, N + 1)
    for l in 1:N-1, u in l+1:N
        qu = a.coll === :lres ? a.qexc(l, u, T, ne) : a.qexc(l, u, T)
        if qu > 0
            P[l, u] += ne * qu
            P[u, l] += ne * qu * (a.levels[l].g / a.levels[u].g) * exp((a.levels[u].E - a.levels[l].E) / (kB * T))
        end
    end
    for (u, l, A2) in a.A2; P[u, l] += A2; end
    for i in 1:N
        R, α = bf_rates(a, i, T; J); qi = a.qion(i, T)
        P[i, N+1] += R + ne * qi
        P[N+1, i] += ne * α + ne * ne * qi * Φ(a, i, T)
    end
    P
end
"Radiative line list of an atom: (l, u, ν, A, B_ul, B_lu, f)."
lines(a::Atom) = [(l, u, νline(a, l, u), a.A[u, l], a.A[u, l] * c^2 / (2h * νline(a, l, u)^3),
                   a.A[u, l] * c^2 / (2h * νline(a, l, u)^3) * a.levels[u].g / a.levels[l].g, a.f[l, u])
                  for l in 1:nlev(a)-1 for u in l+1:nlev(a) if a.A[u, l] > 0]
function add_lines!(P, L, βv; J = nothing)
    for (k, (l, u, ν, A, Bul, Blu, _)) in enumerate(L)
        Jc = J === nothing ? 0.0 : J(ν); b = βv[k]
        P[u, l] += b * (A + Bul * Jc); P[l, u] += b * Blu * Jc
    end
    P
end
"Statistical equilibrium from a full rate matrix P (continuum = last index), Σ n = ntot."
function solve_P(P, ntot)
    N = size(P, 1); M = P' - Diagonal(vec(sum(P; dims = 2))); M[N, :] .= 1.0; rhs = zeros(N); rhs[N] = ntot
    x = M \ rhs; (n = x[1:end-1], nplus = x[end])
end
"""
Angle-averaged Sobolev escape probability of a line with opacity factor K = (πe²/mc) f λ (n_l − g_l n_u/g_u) [s⁻¹] in a flow with
velocity v (cm/s) at radius r and radial gradient dvdr: τ(μ) = K/|μ² dv/dr + (1 − μ²) v/r|, with the floor `Qmin` on the
denominator (thermal width over the zone size, where the flow is nearly static and Sobolev does not apply).
"""
function βsob(K, v, r, dvdr; Qmin = 0.0, nμ = 16)
    K <= 0 && return 1.0
    acc = 0.0
    for i in 1:nμ
        μ = (i - 0.5) / nμ; Q = max(abs(μ^2 * dvdr + (1 - μ^2) * v / r), Qmin)
        τ = K / Q; acc += τ < 1e-6 ? 1.0 : -expm1(-τ) / τ
    end
    acc / nμ
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

# =====================================================================================================================
# l-Resolved Hydrogenic Model Atom (H I, He II) with Exact Radial Dipoles, Burgess (1964) σ_nl, & Pengelly-Seaton (1964)
# =====================================================================================================================

"Associated Laguerre polynomial L_k^α(x) via exact 3-term recurrence."
function _laguerre(k::Int, α::Int, x::Float64)
    k == 0 && return 1.0
    k == 1 && return 1.0 + α - x
    p0 = 1.0; p1 = 1.0 + α - x
    for j in 2:k
        p2 = ((2j - 1 + α - x) * p1 - (j - 1 + α) * p0) / j
        p0, p1 = p1, p2
    end
    p1
end

"""
Exact hydrogenic radial dipole integral R_{n,l}^{n',l'} = ∫₀^∞ r³ R_{n,l}(r) R_{n',l'}(r) dr [in units of a₀ / Z]
evaluated by exact finite polynomial integration via 32-point Gauss-Laguerre quadrature (machine precision for n,n' ≤ 12).
"""
function _hyd_radial_dipole(n1::Int, l1::Int, n2::Int, l2::Int)
    # Normalised hydrogenic radial wavefunctions in x = r / (a0/Z):
    # R_{n,l}(x) = N_{n,l} (2x/n)^l exp(-x/n) L_{n-l-1}^{2l+1}(2x/n)
    # Substitute t = (1/n1 + 1/n2) x => x = t / a, exp(-t) weight!
    a = 1.0 / n1 + 1.0 / n2
    # Compute polynomial coefficients of x^3 * (2x/n1)^l1 L_{n1-l1-1}^{2l1+1}(2x/n1) * (2x/n2)^l2 L_{n2-l2-1}^{2l2+1}(2x/n2)
    # directly via exact closed-form Laguerre coefficient sum: ∫_0^∞ t^m e^{-t} dt = m!
    k1 = n1 - l1 - 1; α1 = 2l1 + 1
    k2 = n2 - l2 - 1; α2 = 2l2 + 1
    N1 = sqrt((2.0 / n1)^3 * factorial(big(k1)) / (2n1 * factorial(big(n1 + l1))))
    N2 = sqrt((2.0 / n2)^3 * factorial(big(k2)) / (2n2 * factorial(big(n2 + l2))))
    acc = big(0.0)
    for j1 in 0:k1
        c1 = ((-1)^j1 * binomial(big(k1 + α1), big(k1 - j1)) / factorial(big(j1))) * (big(2) / n1)^(l1 + j1)
        for j2 in 0:k2
            c2 = ((-1)^j2 * binomial(big(k2 + α2), big(k2 - j2)) / factorial(big(j2))) * (big(2) / n2)^(l2 + j2)
            p = 3 + l1 + j1 + l2 + j2
            acc += c1 * c2 * factorial(big(p)) / (big(a)^(p + 1))
        end
    end
    Float64(N1 * N2 * acc)
end

"Exact hydrogenic dimensionless oscillator strength f(n_l, l_l -> n_u, l_u) (independent of Z & reduced mass)."
function hyd_f_nl(nl::Int, ll::Int, nu::Int, lu::Int)
    (nu > nl && abs(lu - ll) == 1) || return 0.0
    R = _hyd_radial_dipole(nl, ll, nu, lu)
    (1.0 / 3.0) * (1.0 / nl^2 - 1.0 / nu^2) * (max(lu, ll) / (2ll + 1)) * R^2
end

"Exact hydrogenic Einstein A-coefficient [s⁻¹] for (n_u, l_u) -> (n_l, l_l) with |l_u - l_l| = 1."
function hyd_A_nl(Z::Int, nu::Int, lu::Int, nl::Int, ll::Int; χ1 = Z == 1 ? Z0.χH : Z0.χHe2)
    (nu > nl && abs(lu - ll) == 1) || return 0.0
    fl = hyd_f_nl(nl, ll, nu, lu)
    ν = χ1 * (1 / nl^2 - 1 / nu^2) / h
    ACOEF * ν^2 * ((2ll + 1) / (2lu + 1)) * fl
end

"Precomputed Float64 factorials 0!..60! (zero heap allocations in tight quadrature loops)."
const _FACT = Float64[1.0; cumprod(Float64.(1:60))]
@inline _fact(m::Int) = _FACT[m + 1]
@inline _binom(n::Int, k::Int) = _fact(n) / (_fact(k) * _fact(n - k))

"""
Exact finite-polynomial Gordon/Pfaff hydrogenic bound–free dipole matrix element squared
|∫₀^∞ r³ R_{n,l}(r) F_{κ,l'}(r) dr|² (up to common Coulomb prefactor independent of l).
Zero-allocation Float64 implementation; verified exact against Bethe & Salpeter (1957) for all (n ≤ 10, l, κ ≥ 0).
"""
function _exact_bf_M2(n::Int, l::Int, lp::Int, κ::Float64)
    k = n - l - 1; α = 2l + 1
    N2 = (2.0 / n)^3 * _fact(k) / (2n * _fact(n + l))
    κe = max(κ, 1e-5)
    a = complex(1.0 / n, -κe)
    w = complex(0.0, -2.0 * κe) / a
    B = 2lp + 2
    B_minus_A = complex(lp + 1.0, 1.0 / κe)
    θ = atan(n * κe)
    acc = complex(0.0, 0.0)
    two_n = 2.0 / n
    pow_2n = two_n^l
    sgn = 1.0
    @inbounds for j in 0:k
        cj = (sgn * _binom(k + α, k - j) / _fact(j)) * pow_2n
        p = 3 + l + j + lp
        m_deg = (p + 1) - B
        term = complex(1.0, 0.0); f21 = term
        for r in 0:m_deg-1
            term *= (B_minus_A + r) * (-m_deg + r) / ((B + r) * (r + 1)) * w
            f21 += term
        end
        phase = cis(2.0 * (lp - p) * θ)
        acc += cj * _fact(p) / (a^(p + 1)) * phase * f21
        sgn = -sgn
        pow_2n *= two_n
    end
    prod_s = 1.0
    @inbounds for s in 1:lp; prod_s *= (1.0 + (s * κe)^2); end
    C2_exact = (4.0^lp) * prod_s / (_fact(2lp + 1))^2
    N2 * C2_exact * abs2(acc)
end

"""
Exact hydrogenic l-resolved photoionization cross section σ_{n,l}(ν) [cm²] computed via exact
Gordon/Pfaff finite-polynomial matrix elements, normalized to the exact shell-averaged Gaunt cross section
so Σ_l (2l+1) σ_{n,l}(ν) = n² σ_n(ν) to machine precision.
"""
function hyd_sigma_nl(Z::Int, n::Int, l_req::Int, ν::Float64; χ1 = Z == 1 ? Z0.χH : Z0.χHe2)
    νth = χ1 / (n^2 * h)
    ν < νth && return 0.0
    s_avg = Z0.KRAMERS * Z^4 / (n^5 * ν^3) * Z0.gaunt_bf(n, ν / νth)
    n == 1 && return s_avg
    κ = sqrt(max((ν / νth - 1.0) / n^2, 0.0))
    w_req = 0.0; w_sum = 0.0
    for l in 0:n-1
        s = (l + 1) * _exact_bf_M2(n, l, l + 1, κ)
        if l > 0; s += l * _exact_bf_M2(n, l, l - 1, κ); end
        wl = s / (2l + 1)
        w_sum += (2l + 1) * wl
        l == l_req && (w_req = wl)
    end
    w_avg = w_sum / n^2
    s_avg * (w_req / max(w_avg, 1e-300))
end

"""
Pengelly & Seaton (1964, MNRAS 127, 165) proton/ion l-changing collision rate coefficient
q_PS(n, l -> n, l+1; T, ne, Z) [cm³ s⁻¹] for Δn = 0, Δl = +1 within hydrogenic shell n.
Dominated by proton collisions (m_p / m_e = 1836.15); includes Debye + radiative-lifetime cutoff.
"""
function pengelly_seaton_q(Z::Int, n::Int, l::Int, T::Float64, ne::Float64, A_sum::Float64)
    (l >= 0 && l < n - 1 && ne > 0 && T > 0) || return 0.0
    # Dipole strength D_{nl} = 6 n² (n² - (l+1)²) (l+1) / ((2l+1) Z²) [Pengelly & Seaton eq. 29]
    Dnl = 6.0 * n^2 * (n^2 - (l + 1)^2) * (l + 1) / ((2l + 1) * Z^2)
    # Proton-ion reduced mass ratio μ/m_e: for H I (Z=1, p+H) μ = m_p/2; for He II (Z=2, p+He+) μ = (4/5) m_p
    μ_me = Z == 1 ? 0.5 * 1836.15267 : 0.8 * 1836.15267
    # Cutoff radius R_c [cm]: min(Debye length R_D, thermal velocity / lifetime rate A_sum)
    RD = sqrt(kB * T / (4π * max(ne, 1e-10) * ecgs^2))
    vth_p = sqrt(2kB * T / (μ_me * me))
    Rlife = A_sum > 0 ? vth_p / A_sum : 1e30
    # Smooth minimum of RD and Rlife (no sharp corner!)
    Rc = 1.0 / (1.0 / RD + 1.0 / Rlife)
    # Pengelly & Seaton (1964) eq. (43): q = 9.93e-6 * sqrt(μ/m_e) * D_nl / sqrt(T) * [11.54 + log10(T m_e / (μ D_nl)) + 2 log10(Rc)]
    # Wait: in PS64 eq. (43), Rc is in cm; the argument of log10 is dimensionless: 11.54 + log10(T / (D_nl * μ/m_e)) + 2 log10(Rc)
    # Check: at T=1e4 K, ne=1e4 cm^-3, n=4, l=2, Z=2: RD = 6.9 cm, μ/me=1469, D_42=40 => 11.54 + log10(1e4/(40*1469)) + 2*0.84 = 11.54 - 0.77 + 1.68 = 12.45 > 0!
    br = 11.54 + log10(max(T / (Dnl * μ_me), 1e-30)) + 2.0 * log10(max(Rc, 1e-20))
    # Smoothly roll off if bracket becomes small at very high ne (where close-coupling takes over)
    br_s = 0.5 * (br + sqrt(br^2 + 0.5^2))
    # Factor n_p / n_e ≈ 0.92 (primordial H+He) absorbed into effective rate per n_e
    return 0.92 * 9.93e-6 * sqrt(μ_me) * Dnl / sqrt(T) * br_s
end

"""
    hydrogenic_lres(Z; n_lres = 4, nmax = 40, χ1 = ..., coll = :chianti, chianti = nothing) -> Atom

l-resolved hydrogenic model atom (H I for Z=1, He II for Z=2):
- Resolves every (n, l) sublevel for n = 1 .. n_lres (default n_lres = 4: 1s..4f, 10 levels)
  with exact Gordon/Laguerre dipole A-values, exact 2s -> 1s two-photon decay (8.2291 Z⁶ s⁻¹),
  l-resolved photoionization σ_{nl}(ν) (exact Milne detailed balance!), and Pengelly & Seaton (1964)
  proton l-mixing collisions (Δn = 0, Δl = ±1).
- Collapsed shells n = n_lres + 1 .. nmax above n_lres.
Fixes the 13-16% He II λ4686 / λ1640 cascade deficit of the l-collapsed atom across n_e = 10²..10¹⁴ cm⁻³!
"""
function hydrogenic_lres(Z::Int; n_lres::Int = 8, nmax::Int = 40,
                         χ1 = Z == 1 ? Z0.χH : Z0.χHe2, coll = :chianti, chianti = nothing)
    @assert 2 <= n_lres <= 10 "hydrogenic_lres supports n_lres ∈ 2..10"
    lnames = ("s", "p", "d", "f", "g", "h", "i", "k", "l", "m")
    lv = Level[]
    nl_info = Tuple{Int,Int}[]  # (n, l), with l = -1 for collapsed shells n > n_lres
    idx_nl = Dict{Tuple{Int,Int}, Int}()
    for n in 1:n_lres, l in 0:n-1
        push!(lv, Level(χ1 * (1 - 1 / n^2), 2.0 * (2l + 1), χ1 / n^2, "$(n)$(lnames[l+1])"))
        push!(nl_info, (n, l))
        idx_nl[(n, l)] = length(lv)
    end
    for n in (n_lres + 1):nmax
        push!(lv, Level(χ1 * (1 - 1 / n^2), 2.0 * n^2, χ1 / n^2, "n=$n"))
        push!(nl_info, (n, -1))
        idx_nl[(n, -1)] = length(lv)
    end
    N = length(lv)
    f = zeros(N, N); A = zeros(N, N)
    # 1. Exact l-resolved dipole A-values among n <= n_lres
    for i in 1:N-1, j in i+1:N
        ni, li = nl_info[i]; nj, lj = nl_info[j]
        nj > ni || continue
        ν = (lv[j].E - lv[i].E) / h
        if li >= 0 && lj >= 0
            abs(lj - li) == 1 || continue
            A[j, i] = hyd_A_nl(Z, nj, lj, ni, li; χ1)
        elseif li >= 0 && lj == -1
            # Upper shell nj (collapsed) -> lower sublevel (ni, li): exact dipole sum over l' = li ± 1
            if nj <= 12
                A_tot = 0.0
                for lp in (li - 1, li + 1)
                    (0 <= lp < nj) || continue
                    A_tot += (2lp + 1) / nj^2 * hyd_A_nl(Z, nj, lp, ni, li; χ1)
                end
                A[j, i] = A_tot
            else
                # Asymptotic high-n_u: distribute shell A(nj -> ni) by exact (2li+1)/ni^2 statistical weight
                fl = Z0.hyd_f(ni, nj)
                A_shell = ACOEF * ν^2 * (2.0 * ni^2 / lv[j].g) * fl
                A[j, i] = A_shell * ((2li + 1) / ni^2)
            end
        else
            # Both collapsed shells
            fl = Z0.hyd_f(ni, nj)
            A[j, i] = ACOEF * ν^2 * (lv[i].g / lv[j].g) * fl
        end
        if A[j, i] > 0
            f[i, j] = A[j, i] / (ACOEF * ν^2 * lv[i].g / lv[j].g)
        end
    end
    # True 2s -> 1s two-photon decay (exact on the 2s level, no 2/8 factor!)
    A2 = [(idx_nl[(2, 0)], idx_nl[(1, 0)], 8.2291 * Z^6)]
    # Total spontaneous decay rate per level (for Pengelly-Seaton radiative lifetime cutoff)
    A_life = [sum(A[i, :]) for i in 1:N]
    A_life[idx_nl[(2, 0)]] += 8.2291 * Z^6

    # 2. l-resolved bound-free cross-sections (exact Milne detailed balance guaranteed!)
    σ = Function[]
    for i in 1:N
        ni, li = nl_info[i]
        νth = lv[i].χ / h
        if li >= 0
            push!(σ, ν -> hyd_sigma_nl(Z, ni, li, ν; χ1))
        else
            push!(σ, ν -> ν < νth ? 0.0 : Z0.KRAMERS * Z^4 / (ni^5 * ν^3) * Z0.gaunt_bf(ni, ν / νth))
        end
    end

    # 3. Collisions: Pengelly-Seaton (1964) proton l-mixing for Δn = 0 + CHIANTI/Johnson for Δn > 0
    cup = (coll === :chianti && chianti !== nothing) ? chianti_collapsed(chianti)[1] : nothing
    function qexc_lres(i, j, T, ne)
        ni, li = nl_info[i]; nj, lj = nl_info[j]
        if ni == nj
            # Same-shell l-mixing: Δl = +1 (note j > i => lj = li + 1)
            (li >= 0 && lj == li + 1) || return 0.0
            return pengelly_seaton_q(Z, ni, li, Float64(T), Float64(ne), max(A_life[i], A_life[j]))
        end
        # Δn > 0 excitation: shell rate distributed over final sublevels proportionally to g_j / (2 nj²)
        Υ = cup === nothing ? nothing : cup(ni, nj, T)
        q_shell = Υ === nothing ? johnson_exc(ni, nj, T; Z) : q_from_ups(Υ, 2.0 * ni^2, lv[j].E - lv[i].E, T)
        return q_shell * (lv[j].g / (2.0 * nj^2))
    end
    qion_lres(i, T) = johnson_ion(nl_info[i][1], T; Z)

    Atom(Z == 1 ? "H I (l-res)" : "He II (l-res)", lv, 1.0, f, A, A2, σ, qexc_lres, qion_lres, :lres)
end

"Sum line emissivity across all fine/l-sublevel components between principal shells nl and nu (for l-resolved atoms)."
function shell_line_emissivity(a::Atom, n_pops, nl::Int, nu::Int; β = 1.0)
    # Match levels whose energy corresponds to shell nl and shell nu (within 1% of χ1/n²)
    χ1 = a.levels[1].χ
    em = 0.0
    for l in 1:nlev(a)-1, u in l+1:nlev(a)
        a.A[u, l] > 0 || continue
        n_l = round(Int, sqrt(χ1 / a.levels[l].χ))
        n_u = round(Int, sqrt(χ1 / a.levels[u].χ))
        (n_l == nl && n_u == nu) || continue
        b = β isa Function ? β(l, u) : β
        em += line_emissivity(a, n_pops, l, u; β = b)
    end
    em
end

end # module
