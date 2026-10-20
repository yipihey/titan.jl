"""
    NucNet

Population III project, S1: a 19-species thermonuclear network (pp chains, deuterium burning, CNO I/II, 3α, ¹²C(α,γ))
with JINA REACLIB rates, plasma screening, neutrino-loss-corrected energy release, and implicit integrators at fixed
(ρ, T) and at fixed (ρ, u) (T implicit), in the structure of `ChemNet.step`.

- **Rates**: REACLIB 2 format, read from `ENV["NUCNET_REACLIB"]` (default ~/Projects/popiii_lit/reaclib/
  reaclib_default2_20250330, the JINA REACLIB snapshot shipped with pynucastro). λ = Σ_sets exp(a0 + a1/T9 + a2 T9^(-1/3) +
  a3 T9^(1/3) + a4 T9 + a5 T9^(5/3) + a6 ln T9), all non-reverse sets of a reaction summed. Units N_A^(n−1)<σv>; the
  abundance rate is ρ^(n−1) λ Π Y_i / (identical-particle factorial). Reverse (photodisintegration) rates are not included:
  for T9 ≤ 0.3 the Q values of all included channels (≥ 1.94 MeV for p-captures, 7.16 MeV for 12C(α,γ), 7.275 for 3α) give
  factors exp(−Q/kT) ≤ 1e-4 relative to the forward rate and no channel is near equilibrium. Hence not needed here.
- **Electron capture / pep** (density dependent): REACLIB's 'ec' sets (⁷Be(e⁻,ν)⁷Li, p+p+e⁻→d ν) are N_A<σv> per unit electron density; λ_eff = λ ρ Y_e.
- **Energy**: ε = N_A Σ_r (Q_r − ⟨E_ν⟩_r) rate_r, ε_ν = N_A Σ ⟨E_ν⟩ rate_r.
- **Screening**: `:chugunov` (Chugunov, DeWitt & Yakovlev 2007, PRD 76, 025028; multi-component form as in
  pynucastro/MESA, Yakovlev et al. 2006), `:weak` (Salpeter 1954), `:none`.
"""
module NucNet

using ForwardDiff, LinearAlgebra, Printf

export SP, IDX, NET, burn, burn_coupled, evolve, steady, energy, ydot, composition, rates_table, screen_factor

const NA = 6.02214076e23
const MeV = 1.602176634e-6                 # erg
const kB = 1.380649e-16
const hbar = 1.054571817e-27
const qe = 4.803204712570263e-10           # esu
const mu_u = 1.66053906660e-24             # g
const REACLIB_FILE = get(ENV, "NUCNET_REACLIB", joinpath(homedir(), "Projects", "popiii_lit", "reaclib", "reaclib_default2_20250330"))

# --- species ---------------------------------------------------------------------------------------------------------
const SP = [:p, :d, :he3, :he4, :li7, :be7, :b8, :c12, :c13, :n13, :n14, :n15, :o15, :o16, :o17, :f17, :o18, :f18, :f19]
const NS = length(SP)
const IDX = Dict(s => i for (i, s) in enumerate(SP))
const ZZ = Float64[1, 1, 2, 2, 3, 4, 5, 6, 6, 7, 7, 7, 8, 8, 8, 9, 8, 9, 9]
const AA = Float64[1, 2, 3, 4, 7, 7, 8, 12, 13, 13, 14, 15, 15, 16, 17, 17, 18, 18, 19]
# Atomic mass excesses [keV], AME2020 (Wang et al. 2021, Chin. Phys. C 45, 030003), entered by hand for the Q cross-check.
const MEXC = Dict(:n => 8071.3171, :p => 7288.9705, :d => 13135.7216, :t => 14949.8109, :he3 => 14931.2181, :he4 => 2424.9156,
                  :li7 => 14907.1052, :be7 => 15768.99, :b8 => 22921.58, :c12 => 0.0, :c13 => 3125.0089, :n13 => 5345.48,
                  :n14 => 2863.4171, :n15 => 101.4382, :o15 => 2855.604, :o16 => -4737.0014, :o17 => -808.763,
                  :f17 => 1951.70, :o18 => -782.816, :f18 => 873.13, :f19 => -1487.44)
const MEX = [MEXC[s] for s in SP]

# --- REACLIB parser ----------------------------------------------------------------------------------------------------
struct RSet
    chapter::Int
    reac::Vector{String}; prod::Vector{String}
    label::String; resflag::Char; reverse::Bool
    Q::Float64
    a::NTuple{7,Float64}
end
const NREAC = (1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4)       # reactants per REACLIB chapter 1..11

"Parse the REACLIB 2 file; keep only sets whose nuclei are all in `keep`."
function parse_reaclib(path::AbstractString, keep)
    lines = readlines(path); out = RSet[]; i = 1
    while i <= length(lines) - 3
        ch = tryparse(Int, strip(lines[i]))
        if ch === nothing || !(1 <= ch <= 11); i += 1; continue; end
        l2 = rpad(lines[i+1], 80)
        names = [strip(l2[a:a+4]) for a in 6:5:31 if !isempty(strip(l2[a:a+4]))]
        if all(n -> n in keep, names)
            nr = NREAC[ch]
            lab = strip(l2[44:47]); fl = l2[48]; rev = l2[49] == 'v'
            Q = parse(Float64, strip(l2[53:64]))
            c = Float64[]
            for k in 1:13:40; push!(c, parse(Float64, strip(lines[i+2][k:k+12]))); end
            for k in 1:13:27; push!(c, parse(Float64, strip(lines[i+3][k:k+12]))); end
            push!(out, RSet(ch, String.(names[1:nr]), String.(names[nr+1:end]), String(lab), fl, rev, Q, NTuple{7,Float64}(c)))
        end
        i += 4
    end
    out
end

"""
    reaclib_lambda(sets, T9)

N_A^(n−1)<σv> [cm^(3(n−1)) mol^(1−n) s^-1] or decay rate [s^-1] summed over `sets`.
"""
@inline function reaclib_lambda(sets::Vector{NTuple{7,Float64}}, T9)
    t13 = cbrt(T9); ti = 1 / T9; lT = log(T9); s = zero(T9)
    for a in sets
        s += exp(a[1] + a[2] * ti + a[3] / t13 + a[4] * t13 + a[5] * T9 + a[6] * t13^5 + a[7] * lT)
    end
    s
end

# --- reaction table ------------------------------------------------------------------------------------------------------
struct Reaction
    name::String
    group::Symbol                          # :pp, :d, :cno, :he
    rin::NTuple{3,Int}                     # indices of rate-law reactants (0 = none)
    nin::Int
    fact::Float64                          # identical-particle factorial
    net::Vector{Tuple{Int,Float64}}        # net change of Y per event
    sets::Vector{NTuple{7,Float64}}
    Q::Float64                             # MeV (REACLIB or lumped)
    Enu::Float64                           # MeV mean neutrino energy
    kind::Symbol                           # :std, :ec (electron capture: × ρ Ye)
    screens::Vector{NTuple{4,Float64}}     # (Z1, A1, Z2, A2) pairs, sequential
    weak::Int                              # change of nuclear charge per event (expected: 0 strong, -1 weak)
    source::String
end

const SETS = Ref{Vector{RSet}}(RSet[])
function load_sets()
    isempty(SETS[]) && (SETS[] = parse_reaclib(REACLIB_FILE, Set(["p", "d", "t", "n", "he3", "he4", "li7", "be7", "b8", "c12", "c13",
        "n13", "n14", "n15", "o15", "o16", "o17", "f17", "o18", "f18", "f19"])))
    SETS[]
end

function pick(rin, rout, label)
    r = sort(String.(rin)); p = sort(String.(rout))
    s = filter(x -> !x.reverse && sort(x.reac) == r && sort(x.prod) == p && (label === nothing || x.label == label), load_sets())
    isempty(s) && error("NucNet: no REACLIB set for $rin -> $rout ($label)")
    s
end

function mkreaction(name, group, rin, rout; label = nothing, Enu = 0.0, kind = :std, net = nothing, Q = nothing, screen = :auto, wk = nothing)
    ss = pick(rin, rout, label)
    ri = Int[IDX[Symbol(s)] for s in rin]; ro = [Symbol(s) for s in rout]
    nin = length(ri); ri3 = (ri..., ntuple(_ -> 0, 3 - nin)...)
    cnt = Dict{Int,Int}(); for i in ri; cnt[i] = get(cnt, i, 0) + 1; end
    fact = prod(factorial(v) for v in values(cnt); init = 1)
    if net === nothing
        d = Dict{Int,Float64}()
        for i in ri; d[i] = get(d, i, 0.0) - 1; end
        for s in ro; d[IDX[s]] = get(d, IDX[s], 0.0) + 1; end
        netl = [(k, v) for (k, v) in d if v != 0]
    else
        netl = [(IDX[k], Float64(v)) for (k, v) in net if v != 0]
    end
    q = Q === nothing ? ss[1].Q : Q
    scr = NTuple{4,Float64}[]
    if screen === :auto
        if nin >= 2
            a, b = ri[1], ri[2]; push!(scr, (ZZ[a], AA[a], ZZ[b], AA[b]))
        end
    else
        scr = screen
    end
    dz = sum(ZZ[k] * v for (k, v) in netl)
    Reaction(name, group, ri3, nin, Float64(fact), netl, [s.a for s in ss], q, Enu, kind, scr, Int(round(dz)), join(unique(s.label * s.resflag for s in ss), ","))
end

# mean neutrino energies [MeV]: Bahcall (1989, Neutrino Astrophysics) Table 6.1 / Kippenhahn-Weigert (1990) §18.x:
# pp 0.265 (0.2668), pep 1.442 (line), 7Be EC 0.814 (90% 0.862 + 10% 0.384), 8B 6.71 (6.735), 13N 0.707, 15O 0.997, 17F 0.999,
# 18F 0.382 (mean of beta+ branch; 0.634 endpoint), hep 9.628 (Bahcall: 9.625)
function build_net()
    R = Reaction[]
    add(args...; kw...) = push!(R, mkreaction(args...; kw...))
    # pp chains
    add("p(p,e+nu)d", :pp, ["p", "p"], ["d"]; label = "bet+", Enu = 0.265)
    add("p(p e-,nu)d (pep)", :pp, ["p", "p"], ["d"]; label = "ec", Enu = 1.442, kind = :ec)
    add("d(p,g)he3", :pp, ["p", "d"], ["he3"])
    add("he3(he3,2p)he4", :pp, ["he3", "he3"], ["p", "p", "he4"])
    add("he3(a,g)be7", :pp, ["he4", "he3"], ["be7"])
    add("be7(e-,nu)li7", :pp, ["be7"], ["li7"]; label = "ec", Enu = 0.814, kind = :ec)
    add("li7(p,a)a", :pp, ["p", "li7"], ["he4", "he4"])
    add("be7(p,g)b8", :pp, ["p", "be7"], ["b8"])
    add("b8(e+nu)2a", :pp, ["b8"], ["he4", "he4"]; Enu = 6.71)
    add("he3(p,e+nu)a (hep)", :pp, ["p", "he3"], ["he4"]; label = "bet+", Enu = 9.628)
    # deuterium: d(d,n)he3 with the neutron captured at once by a proton (n(p,g)d: tau ~ 1e-7 s vs 880 s decay), net d + p -> he3;
    # d(d,p)t with t -> he3 (12.3 yr) taken as instantaneous
    ddn = pick(["d", "d"], ["n", "he3"], nothing)
    add("d(d,n)he3 [n+p->d]", :d, ["d", "d"], ["n", "he3"]; net = Dict(:d => -1, :p => -1, :he3 => 1), Q = ddn[1].Q + 2.22457)
    add("d(d,p)t [t->he3]", :d, ["d", "d"], ["p", "t"]; net = Dict(:d => -2, :p => 1, :he3 => 1), Q = pick(["d", "d"], ["p", "t"], nothing)[1].Q + 0.0186, Enu = 0.0129, wk = 1)
    add("d(he3,p)he4", :d, ["d", "he3"], ["p", "he4"])
    # CNO
    add("c12(p,g)n13", :cno, ["p", "c12"], ["n13"])
    add("n13(e+nu)c13", :cno, ["n13"], ["c13"]; Enu = 0.707)
    add("c13(p,g)n14", :cno, ["p", "c13"], ["n14"])
    add("n14(p,g)o15", :cno, ["p", "n14"], ["o15"])
    add("o15(e+nu)n15", :cno, ["o15"], ["n15"]; Enu = 0.997)
    add("n15(p,a)c12", :cno, ["p", "n15"], ["he4", "c12"])
    add("n15(p,g)o16", :cno, ["p", "n15"], ["o16"])
    add("o16(p,g)f17", :cno, ["p", "o16"], ["f17"])
    add("f17(e+nu)o17", :cno, ["f17"], ["o17"]; Enu = 0.999)
    add("o17(p,a)n14", :cno, ["p", "o17"], ["he4", "n14"])
    add("o17(p,g)f18", :cno, ["p", "o17"], ["f18"])
    add("f18(e+nu)o18", :cno, ["f18"], ["o18"]; Enu = 0.382)
    add("f18(p,a)o15", :cno, ["p", "f18"], ["he4", "o15"])
    add("o18(p,a)n15", :cno, ["p", "o18"], ["he4", "n15"])
    add("o18(p,g)f19", :cno, ["p", "o18"], ["f19"])
    add("f19(p,a)o16", :cno, ["p", "f19"], ["he4", "o16"])
    # helium burning (3a screened as a+a, then a+be8)
    add("a(2a,g)c12 (3a)", :he, ["he4", "he4", "he4"], ["c12"]; screen = [(2.0, 4.0, 2.0, 4.0), (2.0, 4.0, 4.0, 8.0)])
    add("c12(a,g)o16", :he, ["he4", "c12"], ["o16"])
    R
end
const NET = build_net()
const NR = length(NET)

# --- screening -------------------------------------------------------------------------------------------------------------
"Plasma state (composition-dependent part), computed from float abundances."
struct Plasma
    ρ::Float64; ytot::Float64; zbar::Float64; z2bar::Float64; ne::Float64; zeta::Float64
end
function Plasma(ρ, Y)
    yt = 0.0; zy = 0.0; z2y = 0.0
    for i in 1:NS; yt += Y[i]; zy += ZZ[i] * Y[i]; z2y += ZZ[i]^2 * Y[i]; end
    Plasma(ρ, yt, zy / yt, z2y / yt, ρ * NA * zy, sqrt(z2y + zy))
end

function smooth_clip(x, limit, start)
    lower, upper = limit < start ? (limit, x) : (x, limit)
    x < min(limit, start) && return lower
    x > max(limit, start) && return upper
    f = (1 - cos(π * (x - min(limit, start)) / (start - limit))) / 2
    (1 - f) * lower + f * upper
end

"""
    chugunov07(T, pl, z1, a1, z2, a2) -> h = ln(screening enhancement)

Chugunov, DeWitt & Yakovlev (2007, PRD 76, 025028) eq. 19 with Γ̃ from eq. 21, in the multi-component form of Yakovlev et al.
(2006) / Chugunov & DeWitt (2009) as implemented in pynucastro `screening/screen.py` (BSD-3; re-implemented here):
Z² → Z1Z2, n_i → n_e/z̃³ with z̃ = (Z1^(1/3)+Z2^(1/3))/2, m_i → 2μ12 m_u. T/T_p is clipped to ≥ 0.1 and Γ ≤ 600 smoothly.
"""
function chugunov07(T, pl::Plasma, z1, a1, z2, a2)
    mu12 = a1 * a2 / (a1 + a2); zz = z1 * z2; zt = 0.5 * (cbrt(z1) + cbrt(z2))
    ni = pl.ne / zt^3; mi = 2 * mu12 * mu_u
    Tp = hbar / kB * qe * sqrt(4π * zz * ni / mi)
    gfac = qe^2 / kB * cbrt(4π / 3) * cbrt(pl.ne)
    Tn = smooth_clip(T / Tp, 0.1, 0.2)
    Γ = smooth_clip(gfac * zz / (zt * Tn * Tp), 600.0, 590.0)
    ζ = cbrt(4 / (3π^2 * Tn^2))
    β = 0.41 - 0.6 / Γ; γc = 0.06 + 2.2 / Γ
    gt = Γ / cbrt(1 + ζ * (0.022 + ζ * (β + γc * ζ)))
    A1 = 2.7822; A2 = 98.34; A3 = sqrt(3.0) - A1 / sqrt(A2); B1 = -1.7476; B2 = 66.07; B3 = 1.12; B4 = 65.0
    h = gt^1.5 * (A1 / sqrt(A2 + gt) + A3 / (1 + gt)) + B1 * gt^2 / (B2 + gt) + B3 * gt^2 / (B4 + gt^2)
    max(h, zero(h))
end

"Weak screening, Salpeter (1954): ln f = 0.188 Z1 Z2 ζ ρ^(1/2) T6^(-3/2), ζ² = Σ (Z² + Z) Y."
weak_h(T, pl::Plasma, z1, a1, z2, a2) = 0.188 * z1 * z2 * pl.zeta * sqrt(pl.ρ) / (T * 1e-6)^1.5

"ln of the screening factor of reaction `r` (sum over its sequential pair screens)."
function screen_h(r::Reaction, T, pl::Plasma, mode::Symbol)
    h = zero(T)
    mode === :none && return h
    for (z1, a1, z2, a2) in r.screens
        h += mode === :chugunov ? chugunov07(T, pl, z1, a1, z2, a2) : weak_h(T, pl, z1, a1, z2, a2)
    end
    h
end
screen_factor(T, ρ, Y, z1, a1, z2, a2; mode = :chugunov) =
    exp(mode === :chugunov ? chugunov07(T, Plasma(ρ, Y), z1, a1, z2, a2) : mode === :weak ? weak_h(T, Plasma(ρ, Y), z1, a1, z2, a2) : 0.0)

# --- electron-capture type rates --------------------------------------------------------------------------------------------
# REACLIB 'ec' sets (7Be(e-,nu)7Li and p+p+e- -> d nu) are rates per unit electron density: N_A<sigma v> [cm^3 mol^-1 s^-1],
# so lambda_eff = lambda_REACLIB * rho * Ye  (n_e = rho Ye N_A).  For 7Be this is the Bahcall (1962)/CF88 expression
# 1.34e-10 T9^(-1/2) [1 - 0.537 T9^(1/3) + 3.86 T9^(2/3) + 0.0027/T9 exp(2.515e-3/T9)] (checked numerically in the test).
@inline be7ec_bahcall(T) = (T9 = T * 1e-9; 1.34e-10 / sqrt(T9) * (1 - 0.537 * cbrt(T9) + 3.86 * cbrt(T9)^2 + 0.0027 / T9 * exp(2.515e-3 / T9)))

# --- rates and ydot ----------------------------------------------------------------------------------------------------------
"""
    rates_k(T, ρ, pl, ye; screen = :chugunov) -> k

Effective event-rate coefficients k_r = λ_r f_screen ρ^(n_r−1) / fact_r [s^-1 per Π Y]; the event rate is k_r Π Y_in.
"""
function rates_k(T, ρ, pl::Plasma, ye; screen = :chugunov)
    T9 = T * 1e-9
    k = Vector{typeof(T9 * 1.0)}(undef, NR)
    for (ir, r) in enumerate(NET)
        λ = reaclib_lambda(r.sets, T9)
        r.kind === :ec && (λ *= ρ * ye)
        if !isempty(r.screens)
            λ *= exp(screen_h(r, T, pl, screen))
        end
        k[ir] = λ * ρ^(r.nin - 1) / r.fact
    end
    k
end

@inline event(r::Reaction, k, Y) = (v = k; r.rin[1] > 0 && (v *= Y[r.rin[1]]); r.rin[2] > 0 && (v *= Y[r.rin[2]]); r.rin[3] > 0 && (v *= Y[r.rin[3]]); v)

"dY/dt [mol/g/s] for state Y with coefficients k."
function ydot(Y, k, frozen = nothing)
    f = zeros(promote_type(eltype(Y), eltype(k)), NS)
    for (ir, r) in enumerate(NET)
        ev = event(r, k[ir], Y)
        for (i, c) in r.net; f[i] += c * ev; end
    end
    frozen === nothing || (for i in frozen; f[i] = 0; end)
    f
end

"Analytic Jacobian ∂(dY/dt)/∂Y at fixed coefficients k."
function jac!(J, Y, k, frozen = nothing)
    fill!(J, 0.0)
    for (ir, r) in enumerate(NET)
        for m in 1:r.nin
            d = k[ir]
            for l in 1:r.nin; l == m || (d *= Y[r.rin[l]]); end
            j = r.rin[m]
            for (i, c) in r.net; J[i, j] += c * d; end
        end
    end
    frozen === nothing || (for i in frozen; J[i, :] .= 0; end)
    J
end

"""
    energy(Y, ρ, T, pl; screen) -> (ε, ε_ν, εgroup)

ε_nuc (neutrino-loss corrected) and ε_ν in erg/g/s; `εgroup` a NamedTuple of ε by group (:pp, :d, :cno, :he).
"""
function energy(Y, k, ye = nothing)
    ε = zero(promote_type(eltype(Y), eltype(k))); εn = zero(ε); g = Dict(:pp => zero(ε), :d => zero(ε), :cno => zero(ε), :he => zero(ε))
    for (ir, r) in enumerate(NET)
        ev = event(r, k[ir], Y)
        e = NA * MeV * (r.Q - r.Enu) * ev
        ε += e; εn += NA * MeV * r.Enu * ev; g[r.group] += e
    end
    (ε, εn, (; pp = g[:pp], d = g[:d], cno = g[:cno], he = g[:he]))
end

"Abundance vector Y from a dict of mass fractions (e.g. `Dict(:p=>0.75, :he4=>0.25)`)."
function composition(X::AbstractDict)
    Y = zeros(NS); for (s, x) in X; Y[IDX[s]] += x / AA[IDX[s]]; end; Y
end

ye_of(Y) = sum(ZZ[i] * Y[i] for i in 1:NS)

# --- integrators ---------------------------------------------------------------------------------------------------------------
const YFLOOR = 1e-30

"Backward-Euler solve Y - Yo - h f(Y, k) = 0 at fixed coefficients k (damped Newton, analytic Jacobian). Returns (Y, converged)."
function solveY(Yo, h, k, fr, tol, maxit, afloor)
    Y = copy(Yo); J = zeros(NS, NS); M = zeros(NS, NS)
    for it in 1:maxit
        f = ydot(Y, k, fr); r = Y .- Yo .- h .* f
        jac!(J, Y, k, fr); M .= -h .* J; for i in 1:NS; M[i, i] += 1; end
        δ = try -(M \ r) catch; return (Y, false); end
        all(isfinite, δ) || return (Y, false)
        λ = 1.0
        for i in 1:NS; Y[i] > afloor && Y[i] + λ * δ[i] < 0 && (λ = min(λ, 0.9 * Y[i] / -δ[i])); end   # species below afloor: LU noise, clamp instead
        Y .+= λ .* δ; for i in 1:NS; Y[i] < 0 && (Y[i] = 0.0); end
        λ == 1.0 && all(i -> abs(δ[i]) <= tol * max(Y[i], afloor), 1:NS) && return (Y, true)
    end
    (Y, false)
end

"""
    burn(Y0, ρ, T, dt; screen = :chugunov, frozen = Int[], tol = 1e-10, maxit = 50) -> (Y, ε, εν, εg, nev, nsub)

Backward-Euler step of the network over `dt` at fixed ρ, T with the analytic Jacobian (verified against ForwardDiff in
the tests), Newton damping that keeps Y ≥ 0, recursive halving on Newton failure. Screening uses the composition at the
start of each substep (lagged). `nev[r]` = integrated number of events of reaction r (mol/g). ε, εν are time averages over dt (Σ h_k ε(Y_{k+1}) / dt), so that
Σ_i ΔY_i M_i c² = −(ε + εν) dt/N_A exactly up to the REACLIB-vs-mass-table Q differences. `frozen`: species held constant
(for steady-state/equilibrium studies).
"""
function burn(Y0, ρ, T, dt; screen = :chugunov, frozen = Int[], tol = 1e-10, maxit = 50, afloor = 1e-16)
    acc = Ref(0.0); accn = Ref(0.0); accg = Ref(zeros(4)); nsub = Ref(0); nev = zeros(NR); nfail = Ref(0)
    fr = isempty(frozen) ? nothing : collect(Int, frozen)
    function attempt(Yo, h)
        k = rates_k(T, ρ, Plasma(ρ, Yo), ye_of(Yo); screen)
        Y, ok = solveY(Yo, h, k, fr, tol, maxit, afloor)
        (Y, ok, k)
    end
    function adv(Yo, h, depth)
        Y, ok, k = attempt(Yo, h)
        if ok
            for i in 1:NS; Y[i] = max(Y[i], 0.0); end
            e, en, g = energy(Y, k)
            acc[] += h * e; accn[] += h * en; accg[] .+= h .* [g.pp, g.d, g.cno, g.he]; nsub[] += 1
            for (ir, r) in enumerate(NET); nev[ir] += h * event(r, k[ir], Y); end
            return Y
        end
        (depth >= 30 || (nfail[] += 1) > 400) && error("NucNet.burn: Newton failed (ρ = $ρ, T = $T, dt = $dt)")
        adv(adv(Yo, h / 2, depth + 1), h / 2, depth + 1)
    end
    Y = adv(collect(Float64, Y0), dt, 0)
    g = accg[] ./ dt
    (Y = Y, ε = acc[] / dt, εν = accn[] / dt, εg = (; pp = g[1], d = g[2], cno = g[3], he = g[4]), nev = nev, nsub = nsub[])
end

"""
    burn_coupled(Y0, ρ, u, T0, dt, eos_u; screen = :chugunov, tol = 1e-10, maxit = 50) -> (Y, T, ε, εν, nev, nsub)

Network and temperature together over `dt` at fixed ρ, as `ChemNet.step`: unknowns (Y, ln T), residuals
Y − Y0 − h f(Y, T) = 0 and eos_u(ρ, T, Y) − u − h ε_nuc(Y, T) = 0 (`u` the specific internal energy at the start; `eos_u(ρ, T, Y)`
returns the specific internal energy [erg/g], must accept ForwardDiff duals). Exact ForwardDiff Jacobian over (Y, ln T) with
the screening composition lagged within a substep. Recursive halving on failure.
"""
function burn_coupled(Y0, ρ, u, T0, dt, eos_u; screen = :chugunov, tol = 1e-10, maxit = 50)
    acc = Ref(0.0); accn = Ref(0.0); nsub = Ref(0); nev = zeros(NR); nfail = Ref(0)
    function attempt(Yo, To, uo, h)
        pl = Plasma(ρ, Yo)
        resid(z) = begin
            Y = z[1:NS]; T = exp(z[NS+1])
            k = rates_k(T, ρ, pl, ye_of(Y); screen)
            f = ydot(Y, k); e = energy(Y, k)[1]
            vcat(Y .- Yo .- h .* f, (eos_u(ρ, T, Y) - uo - h * e) / abs(uo))
        end
        k0 = rates_k(To, ρ, pl, ye_of(Yo); screen)
        Yg, okg = solveY(Yo, h, k0, nothing, 1e-8, 30, 1e-16)           # fixed-T solution as the initial guess
        z = vcat(okg ? Yg : Yo, log(To))
        for it in 1:maxit
            r = resid(z); J = ForwardDiff.jacobian(resid, z)
            δ = try -(J \ r) catch; return (z, false); end
            all(isfinite, δ) || return (z, false)
            λ = 1.0
            for i in 1:NS; z[i] > 1e-20 && z[i] + λ * δ[i] < 0 && (λ = min(λ, 0.9 * z[i] / -δ[i])); end
            abs(δ[NS+1]) * λ > 0.3 && (λ = 0.3 / abs(δ[NS+1]))
            if all(i -> abs(δ[i]) <= tol * max(z[i], 1e-16), 1:NS) && abs(δ[NS+1]) < tol
                z = z .+ δ; for i in 1:NS; z[i] < 0 && (z[i] = 0.0); end
                return (z, true)
            end
            z = z .+ λ .* δ; for i in 1:NS; z[i] < 0 && (z[i] = 0.0); end
        end
        (z, false)
    end
    function adv(Yo, To, uo, h, depth)
        z, ok = attempt(Yo, To, uo, h)
        if ok
            Y = max.(z[1:NS], 0.0); T = exp(z[NS+1])
            k = rates_k(T, ρ, Plasma(ρ, Yo), ye_of(Y); screen); e, en, _ = energy(Y, k)
            acc[] += h * e; accn[] += h * en; nsub[] += 1
            for (ir, r) in enumerate(NET); nev[ir] += h * event(r, k[ir], Y); end
            return (Y, T, uo + h * e)
        end
        (depth >= 30 || (nfail[] += 1) > 400) && error("NucNet.burn_coupled: Newton failed (ρ = $ρ, T0 = $T0, dt = $dt)")
        Y1, T1, u1 = adv(Yo, To, uo, h / 2, depth + 1)
        adv(Y1, T1, u1, h / 2, depth + 1)
    end
    Y, T, _ = adv(collect(Float64, Y0), T0, u, dt, 0)
    (Y = Y, T = T, ε = acc[] / dt, εν = accn[] / dt, nev = nev, nsub = nsub[])
end

"""
    evolve(Y0, ρ, T, tend; rtol = 0.05, dt0 = 1.0, screen, frozen, callback = nothing) -> (Y, steps)

Adaptive-step wrapper over `burn` (first-order backward Euler: each step is accepted if no abundance above `ythr`
changes by more than `rtol`, else halved; accepted steps grow by 1.3). `callback(t, Y, res)` is called after every step.
"""
function evolve(Y0, ρ, T, tend; rtol = 0.05, dt0 = 1.0, screen = :chugunov, frozen = Int[], callback = nothing, ythr = 1e-20, grow = 1.3)
    Y = collect(Float64, Y0); t = 0.0; dt = dt0; n = 0
    while t < tend * (1 - 1e-12)
        dt = min(dt, tend - t)
        res = burn(Y, ρ, T, dt; screen, frozen)
        ch = maximum(abs(res.Y[i] - Y[i]) / max(Y[i], ythr) for i in 1:NS if Y[i] > ythr; init = 0.0)
        if ch > rtol && dt > 1e-6
            dt /= 2; continue
        end
        Y = res.Y; t += dt; n += 1
        callback === nothing || callback(t, Y, res)
        dt *= grow
    end
    (Y = Y, steps = n)
end

"""
    steady(Y0, ρ, T; frozen, screen, tmax = 1e30) -> Y

Pseudo-transient continuation (dt × 3 per accepted burn step from 1 s) until Y no longer changes (rel. 1e-10), with
the species in `frozen` held fixed (e.g. p and he4 for the equilibrium of a catalytic cycle).
"""
function steady(Y0, ρ, T; frozen = Int[], screen = :chugunov, tmax = 1e30, ρT = nothing)
    Y = collect(Float64, Y0); dt = 1.0
    while dt < tmax
        res = burn(Y, ρ, T, dt; screen, frozen)
        ch = maximum(abs(res.Y[i] - Y[i]) / max(res.Y[i], 1e-25) for i in 1:NS)
        Y = res.Y; dt *= 3
        ch < 1e-10 && dt > 1e12 && break
    end
    Y
end

"Table of effective rates λ (N_A^(n−1)<σv> or s^-1, unscreened) at temperature T."
rates_table(T) = [(r.name, reaclib_lambda(r.sets, T * 1e-9)) for r in NET]

end # module
