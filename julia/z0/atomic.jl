# Atomic data and line physics of H I, He I, He II for the Z = 0 opacities (see docs/Z0_EOS_OPACITY.md).
#
#   gaunt_ff     thermally averaged, merged relativistic free-free Gaunt factor, van Hoof et al. 2014 (MNRAS 444, 420) tables
#   gaunt_bf     hydrogenic bound-free: n = 1 exact (Stobbe 1930), n ≥ 2 Seaton 1960 second-order expansion
#   hyd_f        n-averaged hydrogenic oscillator strengths, exact (Gordon 1929 radial integrals in BigFloat); Menzel–Pekeris
#                Kramers–Gaunt form beyond FTAB_N
#   occupation   Hummer & Mihalas 1988 occupation probabilities (charged and neutral perturbers), with the Hubeny, Hummer &
#                Lanz 1994 dissolved-level bound-free fraction
#   stark_pattern exact linear Stark components (parabolic states from the manifold's z matrix, dipole from Gordon integrals)
#   stark_profile quasi-static Holtsmark ions + electrons: impact (Griem) core, quasi-static wing, Doppler + radiative damping
#   HEI_LINES    NIST ASD He I multiplets (`hei_multiplets.tsv`)
using LinearAlgebra

const Ry∞   = 13.605693 * eV
const a0cm  = 5.29177211e-9
const ecgs  = 4.80320471e-10
const ħ     = h / 2π
const πe2mc = 0.0265400      # π e²/(m c) [cm² Hz]
const χ_Hn  = χH             # H I threshold used for level energies (consistent with Saha)

# ---- free-free Gaunt factors ------------------------------------------------------------------------------------------------
struct GauntTable; lg0::Float64; lu0::Float64; d::Float64; g::Matrix{Float64}; end   # g[iu, ig]
function _read_gaunt(path)
    nums = Float64[]; hdr = Float64[]
    for ln in eachline(path)
        startswith(lstrip(ln), "#") && continue
        s = strip(split(ln, "#")[1]); isempty(s) && continue
        v = parse.(Float64, split(s))
        if length(hdr) < 7; append!(hdr, v); else; append!(nums, v); end
    end
    ng, nu = Int(hdr[2]), Int(hdr[3])
    @assert length(nums) == ng * nu "gaunt table size"
    GauntTable(hdr[4], hdr[5], hdr[6], permutedims(reshape(nums, ng, nu)))
end
const GFF = Ref{Any}(nothing)
gff_tables() = (GFF[] === nothing && (GFF[] = (_read_gaunt(joinpath(@__DIR__, "gauntff_merged_Z01.dat")), _read_gaunt(joinpath(@__DIR__, "gauntff_merged_Z02.dat")))); GFF[])
@inline function _lag4(t)   # 4-point Lagrange weights at t ∈ [0,1], nodes −1, 0, 1, 2
    (-t * (t - 1) * (t - 2) / 6, (t + 1) * (t - 1) * (t - 2) / 2, -(t + 1) * t * (t - 2) / 2, (t + 1) * t * (t - 1) / 6)
end
"Thermally averaged free-free Gaunt factor for charge Z ∈ (1, 2) at temperature T and frequency ν."
function gaunt_ff(Z::Int, T, ν)
    t = gff_tables()[Z]; kT = kB * T
    x = (log10(Z^2 * Ry∞ / kT) - t.lg0) / t.d; y = (log10(h * ν / kT) - t.lu0) / t.d
    nu, ng = size(t.g)
    ix = clamp(floor(Int, x), 1, ng - 3); iy = clamp(floor(Int, y), 1, nu - 3)
    wx = _lag4(clamp(x - ix, 0.0, 1.0)); wy = _lag4(clamp(y - iy, 0.0, 1.0)); s = 0.0
    @inbounds for a in 1:4, b in 1:4; s += wy[a] * wx[b] * t.g[iy + a, ix + b]; end
    s
end

# ---- bound-free Gaunt factors ------------------------------------------------------------------------------------------------
"Hydrogenic bound-free Gaunt factor of level n at x = ν/ν_n (held at the threshold value for x < 1, dissolved levels)."
function gaunt_bf(n::Int, x)
    x = max(x, 1.0 + 1e-10)
    if n == 1
        ε = sqrt(x - 1)      # Stobbe's exact ground-state cross section relative to Kramers (edge value 0.7973)
        # σ = 6.304e-18 x⁻⁴ exp(4 − 4 atan(ε)/ε)/(1 − exp(−2π/ε)) relative to Kramers 7.906e-18 x⁻³
        return 0.79732 * exp(4 - 4 * atan(ε) / ε) / (x * (1 - exp(-2π / ε)))
    end
    u = x - 1; a = (n * (u + 1))^(-2 / 3)
    max(0.05, 1 + 0.1728 * a * (u - 1) - 0.0496 * a^2 * (u^2 + 4u / 3 + 1))
end
const KRAMERS = 2.8153e29      # σ_n(ν) = KRAMERS Z⁴/(n⁵ ν³) g

# ---- hydrogenic radial integrals and oscillator strengths ----------------------------------------------------------------------
function _hyp(a, b, c, x)   # terminating ₂F₁(−n, −m; c; x) with BigFloat
    s = big(1.0); t = big(1.0); k = 0
    while k < 400
        (a + k == 0 || b + k == 0) && break
        t *= (a + k) * (b + k) / ((c + k) * (k + 1)) * x; s += t; k += 1
    end
    s
end
"Gordon's ∫ R_{n,l} r R_{n',l−1} r² dr in Bohr radii (Bethe & Salpeter 63.2); n ≠ n', l ≥ 1."
function _gordon(n::Int, l::Int, np::Int)
    setprecision(BigFloat, 512) do
        x = big(-4) * n * np / big(n - np)^2
        F1 = _hyp(big(-n + l + 1), big(-np + l), big(2l), x)
        F2 = _hyp(big(-n + l - 1), big(-np + l), big(2l), x)
        pre = sqrt(factorial(big(n + l)) * factorial(big(np + l - 1)) / (factorial(big(n - l - 1)) * factorial(big(np - l))))
        p = big(4 * n * np)^(l + 1) * big(n - np)^(n + np - 2l - 2) / big(n + np)^(n + np)
        v = (iseven(np - l) ? 1 : -1) / (4 * factorial(big(2l - 1))) * pre * p * (F1 - (big(n - np) / big(n + np))^2 * F2)
        Float64(v)
    end
end
const RCACHE = Dict{NTuple{4,Int},Float64}(); const RLOCK = ReentrantLock()
function radial(n, l, np, lp)
    lock(RLOCK) do
        get!(RCACHE, (n, l, np, lp)) do
            lp == l - 1 ? _gordon(n, l, np) : _gordon(np, lp, n)
        end
    end
end
const FTAB_N = 40
const FTAB = Ref{Any}(nothing)
function _f_exact(n, np)
    s = 0.0
    for l in 0:n-1, lp in (l - 1, l + 1)
        (lp < 0 || lp >= np) && continue
        s += (1 / 3) * (1 / n^2 - 1 / np^2) * max(l, lp) * radial(n, l, np, lp)^2 / n^2
    end
    s
end
function ftable()
    FTAB[] === nothing || return FTAB[]
    lock(RLOCK) do
        FTAB[] === nothing || return
        F = zeros(FTAB_N, FTAB_N)
        for n in 1:12, np in n+1:FTAB_N; F[n, np] = _f_exact(n, np); end
        FTAB[] = F
    end
    FTAB[]
end
"n-averaged absorption oscillator strength n → n' > n, hydrogenic (independent of Z)."
function hyd_f(n::Int, np::Int)
    (n <= 12 && np <= FTAB_N) && return ftable()[n, np]
    r = (n / np)^2     # Kramers × Menzel–Pekeris first-order bound-bound Gaunt factor
    g = 1 - 0.1728 * (1 + r) / (1 - r)^(2 / 3) * n^(-2 / 3)
    (32 / (3 * sqrt(3) * π)) / n^5 / np^3 * (1 / n^2 - 1 / np^2)^-3 * max(g, 0.1)
end

# ---- occupation probabilities --------------------------------------------------------------------------------------------------
"""
Hummer & Mihalas 1988 occupation probability of hydrogenic level n of nuclear charge Z: charged perturbers (density Nion,
with charges folded in by the caller as Σ n_i) and neutral perturbers (density Nneu, radius r₁ = a0).
"""
function occupation(n::Real, Z::Int, Nion, Nneu)
    χn = Z^2 * Ry∞ / n^2
    K = n <= 3 ? 1.0 : (16 / 3) * n / (n + 1)^2
    q = sqrt(Float64(Z)) * ecgs^2 / (sqrt(K) * χn)      # H&M 4.71 (Z_p = 1)
    rn = n^2 * a0cm / Z
    exp(-(4π / 3) * (16 * Nion * q^3 + Nneu * (rn + a0cm)^3))
end

# ---- He I multiplets -------------------------------------------------------------------------------------------------------
struct HeILine; Ei::Float64; gi::Int; Ek::Float64; gk::Int; f::Float64; end
const HEIT = Ref{Any}(nothing)
function hei_data()
    HEIT[] === nothing || return HEIT[]
    Ls = HeILine[]; lev = Dict{Tuple{Float64,Int},Int}()
    for ln in eachline(joinpath(@__DIR__, "hei_multiplets.tsv"))
        startswith(ln, "#") && continue
        p = split(ln); Ei = parse(Float64, p[1]); gi = parse(Int, p[2]); Ek = parse(Float64, p[3]); gk = parse(Int, p[4]); gf = parse(Float64, p[5])
        lev[(Ei, gi)] = 1; lev[(Ek, gk)] = 1
        gf / gi >= 1e-4 && Ek > Ei && push!(Ls, HeILine(Ei, gi, Ek, gk, gf / gi))
    end
    levels = sort!([k for k in keys(lev)])
    HEIT[] = (; lines = Ls, levels)
end

# ---- Stark patterns -----------------------------------------------------------------------------------------------------------
"""
Eigenstates of z within the hydrogenic manifold n at fixed m ≥ 0 (|m| symmetric): q (= n₁ − n₂) values and coefficients on
|n, l, m⟩, l = |m|, …, n − 1. z in units of a0/Z: ⟨n l+1 m|z|n l m⟩ = −(3/2) n √(n² − (l+1)²) √(((l+1)² − m²)/((2l+1)(2l+3))).
"""
function _manifold(n, m)
    am = abs(m); ls = am:n-1; N = length(ls)
    N == 1 && return (qs = [0.0], C = ones(1, 1), ls = collect(ls))
    d = [-(1.5) * n * sqrt(n^2 - (l + 1)^2) * sqrt(((l + 1)^2 - m^2) / ((2l + 1) * (2l + 3))) for l in ls[1:end-1]]
    E = eigen(SymTridiagonal(zeros(N), d))
    (qs = E.values ./ (1.5 * n), C = E.vectors, ls = collect(ls))
end
_angz(l, m) = sqrt(((l + 1)^2 - m^2) / ((2l + 1) * (2l + 3)))                 # ⟨l+1 m|cosθ|l m⟩
_angp_up(l, m) = -sqrt((l + m + 1) * (l + m + 2) / (2 * (2l + 1) * (2l + 3)))  # ⟨l+1 m+1|(x+iy)/√2 · r̂|l m⟩
_angp_dn(l, m) = sqrt((l - m) * (l - m - 1) / (2 * (2l + 1) * (2l - 1)))       # ⟨l−1 m+1|(x+iy)/√2 · r̂|l m⟩
const PATCACHE = Dict{Tuple{Int,Int},Any}()
"""
Linear Stark pattern of n → n': vectors X (shift in units of (e a0 F/Z)/h, ≥ 0 after folding), I (relative strengths, ΣI = 1) and
the total line strength Σ|⟨r⟩|² (a0²/Z², for the check against hyd_f). Exact for n ≤ 8, n' ≤ 16; otherwise an equally spaced model.
"""
function stark_pattern(n::Int, np::Int)
    lock(RLOCK) do
        get!(PATCACHE, (n, np)) do
            if n > 8 || np > 16
                Xm = 1.5 * (np^2 - n^2) * 0.9; X = collect(range(0.1, 1.0, length = 6)) .* Xm
                return (; X, I = fill(1 / 6, 6), S = NaN)
            end
            d = Dict{Float64,Float64}()
            for m in -(n - 1):(n - 1)
                a = _manifold(n, m)
                for dm in (0, 1, -1)
                    mp = m + dm; abs(mp) > np - 1 && continue
                    b = _manifold(np, mp)
                    # dipole matrix in the |l⟩ bases
                    D = zeros(length(b.ls), length(a.ls))
                    for (ia, l) in enumerate(a.ls), (ib, lp) in enumerate(b.ls)
                        abs(lp - l) == 1 || continue
                        R = radial(n, l, np, lp); A = 0.0
                        if dm == 0
                            A = lp == l + 1 ? _angz(l, m) : _angz(lp, m)
                        elseif dm == 1
                            A = lp == l + 1 ? _angp_up(l, m) : _angp_dn(l, m)
                        else       # (x − iy)/√2 lowering: conjugate of raising from (lp, mp)
                            A = lp == l + 1 ? _angp_dn(lp, mp) : _angp_up(lp, mp)
                        end
                        D[ib, ia] = R * A
                    end
                    M = transpose(b.C) * D * a.C
                    for ia in eachindex(a.qs), ib in eachindex(b.qs)
                        s = M[ib, ia]^2; s < 1e-14 && continue
                        x = abs(1.5 * (np * b.qs[ib] - n * a.qs[ia])); x = round(x * 1e6) / 1e6
                        d[x] = get(d, x, 0.0) + s
                    end
                end
            end
            X = sort!(collect(keys(d))); Iv = [d[x] for x in X]; S = sum(Iv)
            (; X, I = Iv ./ S, S)
        end
    end
end

# ---- Holtsmark distribution -------------------------------------------------------------------------------------------------
const HOLT = Ref{Any}(nothing)
"W(β) of Holtsmark, ∫₀^∞ W dβ = 1, by quadrature of (2/π) β ∫ x sin(βx) exp(−x^{3/2}) dx; asymptotic series for β > 12."
function holtsmark_table()
    HOLT[] === nothing || return HOLT[]
    lock(RLOCK) do
        HOLT[] === nothing || return
        lb = range(log(1e-2), log(12.0), length = 240); W = zeros(length(lb))
        dx = 2e-3; xs = dx/2:dx:40.0; e = [x * exp(-x^1.5) for x in xs]
        for (i, l) in enumerate(lb)
            β = exp(l); s = 0.0
            @inbounds for j in eachindex(xs); s += e[j] * sin(β * xs[j]); end
            W[i] = 2 / π * β * s * dx
        end
        HOLT[] = (; lb0 = first(lb), dlb = step(lb), lW = log.(max.(W, 1e-300)), n = length(lb))
    end
    HOLT[]
end
"Holtsmark W(β), β ≥ 0."
function holtsmark(β)
    β <= 0 && return 0.0
    t = holtsmark_table()
    β > 12 && return 1.496 * β^-2.5 * (1 + 5.107 * β^-1.5 + 14.43 * β^-3)
    β < 1e-2 && return 4 / (3π) * β^2
    x = (log(β) - t.lb0) / t.dlb; i = clamp(floor(Int, x) + 1, 1, t.n - 1); f = x - i + 1
    exp((1 - f) * t.lW[i] + f * t.lW[i + 1])
end

# ---- line profile -------------------------------------------------------------------------------------------------------------
"Pseudo-Voigt (Thompson, Cox & Hastings 1987): Doppler 1/e half width D, Lorentz half width γ [Hz]; normalised over x."
function pvoigt(x, D, γ)
    fG = 2 * sqrt(log(2)) * D; fL = 2γ
    f = (fG^5 + 2.69269fG^4 * fL + 2.42843fG^3 * fL^2 + 4.47163fG^2 * fL^3 + 0.07842fG * fL^4 + fL^5)^0.2
    r = fL / f; η = clamp(1.36603r - 0.47719r^2 + 0.11116r^3, 0.0, 1.0); hw = f / 2
    η * (hw / π) / (x^2 + hw^2) + (1 - η) * exp(-log(2) * x^2 / hw^2) * sqrt(log(2) / π) / hw
end
# The Griem no-quenching, no-interference estimate overshoots: with the exact ion pattern this one factor reproduces the Gigosos et al.
# 2003 simulated FWHM of Hβ to 3 % and of Hα to 20 % over n_e = 1e16–1e18 cm⁻³, T = 1e4 K (z0 atomic checks in stellaropacity_test.jl).
const IMPACT_CAL = 0.2
"Griem-type electron-impact half width [Hz] of n → n' of charge Z; lnΛ between the Weisskopf radius and the Debye length, ≥ 1."
function impact_width(nl, nu, Z, ne, T)
    kT = kB * T; v = sqrt(8kT / (π * me))
    λD = sqrt(kT / (8π * ne * ecgs^2)); ρmin = ħ * nu^2 / (me * v * Z)
    lnΛ = max(1.0, log(λD / ρmin))
    C = (9 / 4) * (nl^2 * max(nl^2 - 1, 0) + nu^2 * (nu^2 - 1)) / Z^2
    IMPACT_CAL * (4π / 3) * (ħ / me)^2 * sqrt(2me / (π * kT)) * ne * C * lnΛ / (2π)
end
struct LineProf; ν0::Float64; S::Float64; lx0::Float64; dlx::Float64; lφ::Vector{Float64}; φ0::Float64; end
"φ(|Δν|) for one line on log-spaced offsets, S = (π e²/mc) f n_l w_u/w_l (per Hz φ)."
function line_profile(ν0, S, D, γrad, γe, Xs, Is, C_i, C_e)
    W = max(D, γrad + γe, isempty(Xs) ? 0.0 : maximum(Xs) * C_i)
    x0 = 1e-3 * W; x1 = 1.5 * ν0; nx = 90
    lxs = range(log(x0), log(x1), length = nx); lφ = zeros(nx)
    # group the pattern (≤ 10 bins by shift, log-spaced)
    gX = Float64[]; gI = Float64[]; I0 = 0.0
    if !isempty(Xs)
        nz = findall(>(1e-6), Xs); I0 = sum(Is) - sum(Is[nz]); 
        if !isempty(nz)
            a = minimum(Xs[nz]); b = maximum(Xs[nz]); nb = min(10, length(nz))
            edges = exp.(range(log(a) - 1e-9, log(b) + 1e-9, length = nb + 1))
            for j in 1:nb
                id = [k for k in nz if edges[j] <= Xs[k] < edges[j + 1]]
                isempty(id) && continue
                w = sum(Is[id]); push!(gI, w); push!(gX, (sum(Is[id] .* Xs[id] .^ 1.5) / w)^(2 / 3))
            end
        end
    else
        I0 = 1.0
    end
    βn = exp.(range(log(0.02), log(60.0), length = 80)); dlβ = log(βn[2] / βn[1]); Wβ = holtsmark.(βn)
    γcore = γrad + γe
    for (ix, lx) in enumerate(lxs)
        x = exp(lx); s = I0 * pvoigt(x, D, γcore)
        for (Ig, Xg) in zip(gI, gX)
            A = Xg * C_i; (A <= 0) && continue
            if max(D, γcore) < 0.03 * max(x, A) && x > 3A    # far wing: pure quasi-static, kernel negligible
                s += Ig * 0.5 * holtsmark(x / A) / A
            else
                t = 0.0
                for (b, wb) in zip(βn, Wβ)
                    sh = b * A; t += 0.5 * wb * b * dlβ * (pvoigt(x - sh, D, γcore) + pvoigt(x + sh, D, γcore))
                end
                s += Ig * t
            end
        end
        # electrons: impact Lorentzian (in γcore above) goes over to their quasi-static wing where that is smaller
        if γe > 0 && !isempty(gX) && C_e > 0
            L = (γe / π) / x^2; He = 0.0
            for (Ig, Xg) in zip(gI, gX); Ae = Xg * C_e; He += Ig * 0.5 * holtsmark(x / Ae) / Ae; end
            x > 3 * max(D, γcore) && He < L && (s += He - L)
        end
        lφ[ix] = log(max(s, 1e-300))
    end
    LineProf(ν0, S, lxs[1], step(lxs), lφ, exp(lφ[1]))
end
@inline function prof_eval(p::LineProf, x)
    ax = abs(x); ax <= exp(p.lx0) && return p.φ0
    t = (log(ax) - p.lx0) / p.dlx; i = floor(Int, t) + 1; n = length(p.lφ)
    i >= n && return 0.0
    f = t - i + 1; exp((1 - f) * p.lφ[i] + f * p.lφ[i + 1])
end
"Two-sided normalised profile value including the far-red (ν/ν0)⁴ Rayleigh limit."
@inline function prof_val(p::LineProf, ν)
    v = prof_eval(p, ν - p.ν0); ν < p.ν0 ? v * (ν / p.ν0)^4 : v
end
