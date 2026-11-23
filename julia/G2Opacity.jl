"""
    G2Opacity

G2 ("one non-LTE material everywhere", docs/POPIII_G2_MATERIAL.md): the continuum means of a network zone
(PopIIIMaterial.nlte_means) with the continuum ingredients that the equilibrium-table side (Z0Physics, used for the table's κ_P
above 10^3.8 K; OPAL for its κ_R) has and the network's code lacked. Same arguments and same conventions as nlte_means
(number densities from the zone's own composition; `hbfP` = H bound–free in κ_E/κ_P or not; κ_E weighted with B_ν(Tw)), plus the
helium ions. Switch `G2OPAC=1` (default off; off = no call, runs bit-identical): PopIIIMaterial.nlte_means / nlte_means_fast return
this instead.

Added against nlte_means (sources: docs/Z0_EOS_OPACITY.md and julia/z0/; each is Z0Physics' own function, evaluated on the network's
populations instead of Saha's):
- frequency grid to 25 kT/hc at 3e5 K (nlte_means' grid stops at 1e5 cm⁻¹ = 12.4 eV: above ~3e4 K the Rosseland and Planck weights
  lie beyond it and its means were the low-frequency free–free alone, too high by up to 10³);
- level bookkeeping: `n1` is the zone's neutral H (all levels, ChemNet's x_H), n ≥ 2 follow from `n2`, the ground state is the rest,
  n1 − Σ_{n≥2} p_n (nlte_means took all neutral H as ground state and added the excited levels on top: above ~3e4 K, where the
  partition function is ≫ 2, that overstated the bound–free opacity by the partition function, ×250 at 1e5 K); H⁻ free–free and
  Rayleigh on the ground state; `n2_lte` gives the LTE n = 2 for an equilibrium composition;
- H bound–free with Gaunt factors (n = 1 Stobbe exact, n ≥ 2 Seaton fit: Z0Physics.gaunt_bf) instead of Kramers g = 1, levels
  n = 2 … 16 from the zone's n = 2 by Boltzmann × Hummer & Mihalas (1988) occupation probabilities (Z0Physics.occupation), and the
  dissolved-level pseudo-continua below each edge (Hubeny, Hummer & Lanz 1994), as Z0Physics.kappa_nu;
- He I ground bound–free (Verner et al. 1996, Z0Physics.σHeI_ground) on the zone's neutral He; He II bound–free (hydrogenic Z = 2,
  levels in Boltzmann with the zone's He⁺, occupation probabilities, Gaunt factors);
- free–free of H⁺ and He⁺ (Z = 1) and He²⁺ (Z = 2) with the van Hoof et al. (2014) Gaunt factors (Z0Physics.gaunt_ff) instead of
  H⁺ alone with g = 1;
- Rayleigh cross sections of H and H₂ held at their 1300 Å values at shorter wavelengths (Z0Physics.σRay; the fits diverge towards Lyα).
Not added (see the doc): bound–bound lines of H/He (in OPAL, and Z0Physics' line lists are LTE and ~0.05–4 s per point), He⁻ and H₂⁻
free–free, H₂ lines, pressure effects on H⁻.
"""
module G2Opacity

import ..Z0Physics, ..ContinuumOpacity
import ForwardDiff
const Z0 = Z0Physics; const COp = ContinuumOpacity

const ON = get(ENV, "G2OPAC", "0") in ("1", "2")     # 2: the fast form (`means_fast`)
const NLEV = 16
const TTOP = 3e5            # the grid reaches 25 kT/hc at this T (the HOT branch's T1)

const FREQ = Ref{Any}(nothing)
"Frequency-only quantities (computed once): midpoints, widths, λ, ν³, σ(H⁻ bf), capped Rayleigh σ of H₂ and H, σ(He I ground)."
function freq()
    FREQ[] === nothing || return FREQ[]
    νs = COp.grid(; n = 1500, T = TTOP); n = length(νs) - 1
    ν = [0.5 * (νs[i] + νs[i+1]) for i in 1:n]; dν = [νs[i+1] - νs[i] for i in 1:n]; σ = ν ./ COp.c; λμ = 1e4 ./ σ
    λc = max.(λμ, 0.13)                                           # 1300 Å
    FREQ[] = (; ν, dν, σ, λμ, ν3 = ν .^ 3, σhm = COp.σ_Hm_bf.(λμ), ray2 = COp.σ_ray_H2.(λc .* 1e3), ray1 = COp.σ_ray_H.(λc .* 1e4),
                σHeI = Z0.σHeI_ground.(ν))
end

"""
    means(ρ, T; n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, hbfP, Tw = T) -> (κ_E, κ_P, κ_R) [cm² g⁻¹]

As PopIIIMaterial.nlte_means (n1: H atoms, taken as the ground state; n2: population of n = 2 [cm⁻³]; nHe: neutral He) with the
additions listed in the module docstring.
"""
function means(ρ, T; n1, n2, nH2, nHe, nHep = 0.0, nHepp = 0.0, nHp, ne, nHm, hbfP, Tw = T)
    ρ = Float64(ρ); T = Float64(T); F = freq(); t = COp.tables(); kT = COp.kB * T
    sE = 0.0; sP = 0.0; sB = 0.0; sR = 0.0; sD = 0.0; sW = 0.0; wJ = Tw != T
    Nion = nHp + nHep + nHepp; Nneu = n1 + nHe + nH2
    χH = Z0.χH; χ2 = Z0.χHe2; KR = Z0.KRAMERS
    # level populations: H n ≥ 2 from n2 (Boltzmann relative to n = 2, occupation probabilities w_n/w_2); He II from He⁺ in LTE
    wH = ntuple(n -> Z0.occupation(n, 1, Nion, Nneu), NLEV); wHe = ntuple(n -> Z0.occupation(n, 2, Nion, Nneu), NLEV)
    pH0 = ntuple(n -> n == 1 ? 0.0 : n2 * (n^2 / 4) * exp(-χH * (0.25 - 1 / n^2) / kT) * wH[n] / max(wH[2], 1e-300), NLEV)
    ng = max(n1 - sum(pH0), 0.0)                                  # ground state
    pH = ntuple(n -> n == 1 ? ng : pH0[n], NLEV)
    UHe = sum(n^2 * exp(-χ2 * (1 - 1 / n^2) / kT) * wHe[n] for n in 1:NLEV)
    pHe = ntuple(n -> nHep * n^2 * exp(-χ2 * (1 - 1 / n^2) / kT) * wHe[n] / max(UHe, 1e-300), NLEV)
    νH = χH / COp.h; νHe = χ2 / COp.h
    @inbounds for i in eachindex(F.ν)
        ν = F.ν[i]; dν = F.dν[i]; σ = F.σ[i]; λμ = F.λμ[i]; ν3 = F.ν3[i]
        stim = -expm1(-COp.h * ν / kT)
        acia = (COp.k_h2h2(t, σ, T) * nH2^2 + COp.kcia(t.H2He, σ, T) * nH2 * nHe + COp.kcia(t.H2H, σ, T) * nH2 * n1)      # CIA at every frequency point, always
        # free–free: H⁻ (John 1988), H⁺ + He⁺ (Z = 1) and He²⁺ (Z = 2) with van Hoof Gaunt factors
        aff = COp.k_Hm_ff(λμ, T) * ng * ne * kT
        aff += 3.692e8 * ne / (sqrt(T) * ν3) * stim * ((nHp + nHep) * Z0.gaunt_ff(1, T, ν) + (nHepp > 0 ? 4nHepp * Z0.gaunt_ff(2, T, ν) : 0.0))
        # H bound–free (Gaunt factors, occupation probabilities, dissolved levels)
        ab = 0.0
        for n in 1:NLEV
            pH[n] <= 0 && continue
            νn = νH / n^2
            if ν >= νn
                ab += pH[n] * KR / (n^5 * ν3) * Z0.gaunt_bf(n, ν / νn)
            elseif n < NLEV
                q = 1 / n^2 - ν / νH
                if q <= 1 / 1600
                    d = 1 - Z0.occupation(1 / sqrt(max(q, 1e-12)), 1, Nion, Nneu) / max(wH[n], 1e-300)
                    d > 0 && (ab += d * pH[n] * KR / (n^5 * ν3) * Z0.gaunt_bf(n, 1.0))
                end
            end
        end
        ab *= stim
        # helium bound–free: He I ground (Verner), He II hydrogenic Z = 2
        ahe = nHe * F.σHeI[i]
        if nHep > 0 && ν >= νHe / NLEV^2
            for n in 1:NLEV
                νn = νHe / n^2
                ν >= νn && (ahe += pHe[n] * KR * 16 / (n^5 * ν3) * Z0.gaunt_bf(n, ν / νn))
            end
        end
        ahe *= stim
        a0 = acia + aff + ahe
        σhm = F.σhm[i] * stim
        sca = F.ray2[i] * nH2 + F.ray1[i] * ng + COp.σT * ne
        B = COp.Bν(ν, T); D = COp.dBdT(ν, T); h = hbfP ? ab : 0.0
        BJ = wJ ? COp.Bν(ν, Tw) : B
        sE += (a0 + σhm * nHm + h) * BJ * dν; sW += BJ * dν; sP += (a0 + σhm * nHm + h) * B * dν; sB += B * dν
        sR += D / max(a0 + σhm * nHm + ab + sca, 1e-300) * dν; sD += D * dν
    end
    (sE / sW / ρ, sP / sB / ρ, sD / sR / ρ)
end

"""
LTE population of H n = 2 [cm⁻³] for `nH0` neutral H atoms (all levels): Boltzmann with occupation probabilities over n ≤ $(NLEV)
(Nion, Nneu as in `means`).
"""
function n2_lte(T, nH0, Nion, Nneu)
    kT = COp.kB * T; χH = Z0.χH
    w = ntuple(n -> Z0.occupation(n, 1, Nion, Nneu), NLEV)
    U = sum(n^2 * exp(-χH * (1 - 1 / n^2) / kT) * w[n] / w[1] for n in 1:NLEV)
    nH0 * 4 * exp(-χH * 0.75 / kT) * w[2] / w[1] / U
end

# ------------------------------------------------------------------ fast form (G2OPAC=2): ~280 frequency nodes, T tables
"""
`G2OPAC=2`: `means_fast` instead of `means` — the same physics on ~280 frequency nodes (log-uniform 10 cm⁻¹ … 25 kT/hc at 3e5 K,
plus node pairs at every H, He II (n ≤ 6), He I and H⁻ threshold), with the T-dependent kernels (van Hoof g_ff for Z = 1, 2, H⁻ ff,
the three CIA pairs) tabulated on log T 1.0 … 6.5 (Δ 0.01) and interpolated linearly in log T, and the T-independent cross sections
(H and He II bound–free with Gaunt factors, He I, H⁻ bf, Rayleigh) precomputed per node. Generic in the number type of T and of the
populations (ForwardDiff Duals pass: exact derivatives of the means as computed). Trapezoid rule on the nodes.
"""
const FAST = Ref{Any}(nothing)
const FLT0 = 1.0; const FDLT = 0.01; const FNT = 551; const NFAST = parse(Int, get(ENV, "G2NFAST", "400"))
function fastgrid()
    FAST[] === nothing || return FAST[]
    c = COp.c; h = COp.h
    σmax = 25 * COp.kB * TTOP / (h * c)
    σs = collect(10 .^ range(1.0, log10(σmax), length = NFAST))
    edges = Float64[]
    for n in 1:6; push!(edges, Z0.χH / n^2 / (h * c)); push!(edges, Z0.χHe2 / n^2 / (h * c)); end
    push!(edges, Z0.χHe1 / (h * c)); push!(edges, 1e4 / 1.6419)
    for e in edges, f in (1 - 3e-4, 1 + 3e-4, 1.003, 1.03, 0.997, 0.97); push!(σs, e * f); end
    sort!(σs); filter!(x -> 10.0 <= x <= σmax, σs); unique!(σs)
    ν = σs .* c; N = length(ν)
    w = zeros(N)                                                       # trapezoid weights
    for i in 1:N-1; d = ν[i+1] - ν[i]; w[i] += d / 2; w[i+1] += d / 2; end
    λμ = 1e4 ./ σs; λc = max.(λμ, 0.13)
    KR = Z0.KRAMERS
    SH = zeros(NLEV, N); SHe = zeros(NLEV, N)
    for i in 1:N, n in 1:NLEV
        νn = Z0.χH / h / n^2; ν[i] >= νn && (SH[n, i] = KR / (n^5 * ν[i]^3) * Z0.gaunt_bf(n, ν[i] / νn))
        νe = Z0.χHe2 / h / n^2; ν[i] >= νe && (SHe[n, i] = 16KR / (n^5 * ν[i]^3) * Z0.gaunt_bf(n, ν[i] / νe))
    end
    # dissolved-level windows below the H edges (Hubeny, Hummer & Lanz 1994, as in `means`)
    DIS = Tuple{Int,Int,Float64,Float64}[]
    for i in 1:N, n in 1:NLEV-1
        νn = Z0.χH / h / n^2; ν[i] < νn || continue
        q = 1 / n^2 - ν[i] / (Z0.χH / h); q <= 1 / 1600 && push!(DIS, (i, n, 1 / sqrt(max(q, 1e-12)), KR / (n^5 * ν[i]^3) * Z0.gaunt_bf(n, 1.0)))
    end
    t = COp.tables()
    G1 = zeros(FNT, N); G2 = zeros(FNT, N); HF = zeros(FNT, N); C1 = zeros(FNT, N); C2 = zeros(FNT, N); C3 = zeros(FNT, N)
    Threads.@threads for j in 1:FNT
        T = 10^(FLT0 + (j - 1) * FDLT)
        for i in 1:N
            G1[j, i] = Z0.gaunt_ff(1, T, ν[i]); G2[j, i] = Z0.gaunt_ff(2, T, ν[i]); HF[j, i] = COp.k_Hm_ff(λμ[i], T)
            C1[j, i] = COp.k_h2h2(t, σs[i], T); C2[j, i] = COp.kcia(t.H2He, σs[i], T); C3[j, i] = COp.kcia(t.H2H, σs[i], T)
        end
    end
    FAST[] = (; ν, w, ν3 = ν .^ 3, σhm = COp.σ_Hm_bf.(λμ), ray2 = COp.σ_ray_H2.(λc .* 1e3), ray1 = COp.σ_ray_H.(λc .* 1e4),
                σHeI = Z0.σHeI_ground.(ν), SH, SHe, DIS, G1, G2, HF, C1, C2, C3)
end
"B_ν(T_w) of the in-group shape: the plain form for Floats; a Dual T_w (T_c carrying its derivative) takes the form e^{-u}/(1 − e^{-u}), whose derivative does not overflow in the Wien tail (expm1² of u > 355)."
@inline _bnuw(ν, T) = COp.Bν(ν, T)
@inline _bnuw(ν, T::ForwardDiff.Dual) = (u = COp.h * ν / (COp.kB * T); u > 700 ? zero(u) : 2COp.h * ν^3 / COp.c^2 * exp(-u) / (-expm1(-u)))
@inline function _tint(A, j, a, i)                                     # linear in log T between rows j, j+1
    (1 - a) * A[j, i] + a * A[j+1, i]
end
struct MeansWS{T}
    wH::Vector{T}; wHe::Vector{T}; pH::Vector{T}; pHe::Vector{T}; dis::Vector{T}
    MeansWS{T}(nlev::Int, nnu::Int) where {T} = new{T}(Vector{T}(undef, nlev), Vector{T}(undef, nlev), Vector{T}(undef, nlev), Vector{T}(undef, nlev), zeros(T, nnu))
end
const MWS_F64 = Vector{MeansWS{Float64}}()
const MWS_POOL = Vector{Dict{DataType,Any}}()
const MWS_LOCK = ReentrantLock()
function _init_mws!(nnu::Int)
    nt = Threads.maxthreadid()
    length(MWS_F64) >= nt && return
    lock(MWS_LOCK) do
        while length(MWS_F64) < nt
            push!(MWS_F64, MeansWS{Float64}(NLEV, nnu))
            push!(MWS_POOL, Dict{DataType,Any}())
        end
    end
end
@inline function get_mws(::Type{Float64}, nnu::Int)
    tid = Threads.threadid()
    tid > length(MWS_F64) && _init_mws!(nnu)
    w = @inbounds MWS_F64[tid]
    length(w.dis) < nnu && resize!(w.dis, nnu)
    fill!(w.dis, 0.0)
    w
end
@inline function get_mws(::Type{T}, nnu::Int) where {T}
    tid = Threads.threadid()
    tid > length(MWS_POOL) && _init_mws!(nnu)
    d = @inbounds MWS_POOL[tid]
    w0 = get(d, T, nothing)
    w = w0 === nothing ? (w2 = MeansWS{T}(NLEV, nnu); d[T] = w2; w2) : (w0::MeansWS{T})
    length(w.dis) < nnu && resize!(w.dis, nnu)
    fill!(w.dis, zero(T))
    w
end

"`means` on the fast grid (see FAST); same arguments and conventions, generic in the number type."
means_fast(ρ, T; n1, n2, nH2, nHe, nHep = 0.0, nHepp = 0.0, nHp, ne, nHm, hbfP, Tw = T) =
    _means_fast(fastgrid(), ρ, T, n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, hbfP, Tw)
function _means_fast(F, ρ, T, n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, hbfP, Tw)     # (function barrier: F concrete)
    kB = COp.kB; h = COp.h; kT = kB * T
    Nion = nHp + nHep + nHepp; Nneu = n1 + nHe + nH2
    χH = Z0.χH; χ2 = Z0.χHe2
    Tp = typeof(one(n1) * one(n2) * one(nHep) * one(Nion) * one(Nneu) * one(T))
    N = length(F.ν)
    mws = get_mws(Tp, N)
    wH = mws.wH; wHe = mws.wHe; pH = mws.pH; pHe = mws.pHe; dis = mws.dis
    lwH = ntuple(n -> Z0.ln_occupation(n, 1, Nion, Nneu), NLEV)
    for n in 1:NLEV; wH[n] = exp(lwH[n]); wHe[n] = Z0.occupation(n, 2, Nion, Nneu); end
    sx = zero(Tp)
    for n in 2:NLEV; pH[n] = n2 * (n^2 / 4) * exp(-χH * (0.25 - 1 / n^2) / kT) * wH[n] / max(wH[2], 1e-300); sx += pH[n]; end
    ng = max(n1 - sx, zero(Tp)); pH[1] = ng
    UHe = zero(Tp)
    for n in 1:NLEV; pHe[n] = n^2 * exp(-χ2 * (1 - 1 / n^2) / kT) * wHe[n]; UHe += pHe[n]; end
    for n in 1:NLEV; pHe[n] *= nHep / UHe; end
    lT = log10(ForwardDiff_value(T)); x = clamp((lT - FLT0) / FDLT, 0.0, FNT - 1 - 1e-9); j = floor(Int, x) + 1
    a = clamp((log10(T) - FLT0) / FDLT - (j - 1), 0.0, 1.0)            # (generic weight; clamped at the table ends)
    wJ = Tw != T
    # dissolved-level contributions per node (sparse)
    for (i, n, ns, σe) in F.DIS
        d = 1 - exp(Z0.ln_occupation(ns, 1, Nion, Nneu) - lwH[n])         # (the ratio of two occupation probabilities, from their logarithms: both can be 1e-100)
        d > 0 && (dis[i] += d * pH[n] * σe)
    end
    z = zero(n1 * T * ρ); sE = z; sP = z; sB = z; sR = z; sD = z; sW = z
    @inbounds for i in 1:N
        ν = F.ν[i]; dν = F.w[i]; u = h * ν / kT
        em = exp(-u); stim = 1 - em
        B = 2h * ν^3 / COp.c^2 * em / stim; D = B * u / (T * stim)
        ab = dis[i]
        for n in 1:NLEV; s_ = F.SH[n, i]; s_ > 0 && (ab += pH[n] * s_); end
        ab *= stim
        ahe = nHe * F.σHeI[i]
        if nHep > 0
            for n in 1:NLEV; s_ = F.SHe[n, i]; s_ > 0 && (ahe += pHe[n] * s_); end
        end
        ahe *= stim
        aff = _tint(F.HF, j, a, i) * ng * ne * kT +
              3.692e8 * ne / (sqrt(T) * F.ν3[i]) * stim * ((nHp + nHep) * _tint(F.G1, j, a, i) + 4nHepp * _tint(F.G2, j, a, i))
        acia = _tint(F.C1, j, a, i) * nH2^2 + _tint(F.C2, j, a, i) * nH2 * nHe + _tint(F.C3, j, a, i) * nH2 * n1
        a0 = acia + aff + ahe
        σhm = F.σhm[i] * stim
        sca = F.ray2[i] * nH2 + F.ray1[i] * ng + COp.σT * ne
        hb = hbfP ? ab : zero(ab)
        BJ = wJ ? _bnuw(ν, Tw) : B
        sE += (a0 + σhm * nHm + hb) * BJ * dν; sW += BJ * dν; sP += (a0 + σhm * nHm + hb) * B * dν; sB += B * dν
        sR += D / max(a0 + σhm * nHm + ab + sca, 1e-300) * dν; sD += D * dν
    end
    (_qd(sE, sW) / ρ, _qd(sP, sB) / ρ, _qd(sD, sR) / ρ)
end
include(joinpath(@__DIR__, "HScatter.jl"))               # exact H(1s) Rayleigh + Raman cross sections (MGROUPS only; nothing computed unless used)
"True if the band [lo, hi] (Hz) lies entirely inside the range of the H(1s) Rayleigh/Raman table (0.80 … 0.9996 χ_H: between Lyα and the Lyman limit)."
lwband(lo, hi) = (νH = Z0.χH / COp.h; lo >= HScatter.EPS0 * νH && hi <= HScatter.EPS1 * νH)
const LWSUB = Ref{Any}(nothing)
"""
Sub-grid of a group lying between Lyα and the Lyman limit (`lwband`; docs §16): cells of LWSTRIDE cells of the HScatter table each, with
the H(1s) extinction σ_t = 1/⟨1/(σ_1s + σ_Ram)⟩ (harmonic over the cell: its transmission) and the Raman part σ_r = σ_t ⟨σ_Ram/(σ_1s + σ_Ram)⟩
of it, and for the smooth rest of the extinction the bracketing G2 nodes of the group (linear interpolation in ν). Against the full table
(7044 cells in 11.2 eV … χ_H(1 − 1/1600)) the group means of `group_fast` agree to 1e-3 where H(1s) dominates the extinction and to 2 %
(χ_g) anywhere (docs §16). `edges`: interior group edges [Hz]. Per group `nothing` or (ν, st, sr, nodes, ia, wa).
"""
const LWSTRIDE = 32
function lwsub(edges, gidx, G)
    s = LWSUB[]; (s !== nothing && s.edges == edges) && return s.sub
    F = fastgrid(); H = HScatter.table(); νH = Z0.χH / COp.h
    sub = Vector{Any}(nothing, G)
    for g in 2:G-1
        lwband(edges[g-1], edges[g]) || continue
        js = [j for j in eachindex(H.ε) if edges[g-1] <= H.ε[j] * νH < edges[g]]
        nodes = [i for i in eachindex(F.ν) if gidx[i] == g]
        (isempty(js) || length(nodes) < 2) && continue
        νn = F.ν[nodes]; ν = Float64[]; st = Float64[]; sr = Float64[]; ia = Int[]; wa = Float64[]
        for a in 1:LWSTRIDE:length(js)
            c = js[a:min(a + LWSTRIDE - 1, length(js))]; h = sum(1 / (H.s1[j] + H.sr[j]) for j in c); x = sum(H.ε[j] for j in c) / length(c) * νH
            push!(ν, x); push!(st, length(c) / h); push!(sr, sum(H.sr[j] / (H.s1[j] + H.sr[j]) for j in c) / h)
            q = clamp(searchsortedlast(νn, x), 1, length(nodes) - 1)
            push!(ia, q); push!(wa, clamp((x - νn[q]) / (νn[q+1] - νn[q]), 0.0, 1.0))
        end
        sub[g] = (; ν, st, sr, nodes, ia, wa)
    end
    LWSUB[] = (; edges = copy(edges), sub); sub
end
"""
MGROUPS (docs/POPIII_MULTIGROUP.md): per-group means of `means_fast`'s physics on the same nodes and with the same populations, plus the
grey means of the same call (Float64 only). `gidx[i]` is the group of node i. In-group weights: emission B_ν(T), absorption B_ν(Tw),
extinction ∂B_ν/∂T(T) — the weights of the grey means, so that Σ_g b_g(T) κP_g = κP and Σ_g b_g(Tw) κE_g = κE hold at every T (a group whose
weight underflows returns 0 and is taken as grey by the caller; its Planck fraction is negligible there). (A floor of the in-group shape
temperature broke the emission identity below it: 1000 K cloud gas then over-emitted into the optical group by orders of magnitude, t_mg27.)
Returns (κE_g, κP_g, χ_g, κ1s_g, κE, κP, χ): κ1s_g the H(1s) bound–free absorption per mass with the in-group shape B_ν(max(Tw, Tmin)) (the
star's light, also in unlit cold gas where B_ν(Tw) of the ionising groups underflows), whether or not hbfP puts H bound–free into κE_g.
With `edges` (interior group edges [Hz]) a group between Lyα and the Lyman limit (`lwsub`) takes the exact H(1s) Rayleigh + Raman cross
sections (HScatter; Kokubo 2024) on the table's cells instead of the low-energy Rayleigh expansion of the G2 nodes (invalid above ~8 eV,
docs §16): its χ_g is the harmonic mean of the total extinction (smooth extinction of the nodes + n_1s(σ_1s + σ_Raman)) with the weight
∂B_ν/∂T at the shape temperature max(Tw, Tmin) of the light the group carries, and the 9th return value κ_Ram,g is the Raman absorption
per mass in the in-group spectrum φ_ν ∝ B_ν/χ_ν of that harmonic mean, κ_Ram = ∫κ_Ram,ν φ_ν dν/∫φ_ν dν = χ_g ⟨σ_Ram n_1s/χ_ν⟩ (convergent
in the unresolved line cores, where the photons are trapped and do not travel; 0 for every other group).
8th return value: κ2_g, the bound–free absorption of H n ≥ 2 (and their dissolved levels) per mass with the shape of κ1s_g. 10th: σ̄_1s,g
[cm² per ground-state atom], the mean of σ_1s(ν)(1 − e^{−hν/kT}) in the group with the emission weight B_ν(T): the Milne (detailed
balance) emission of recombinations to the ground state into group g is 4π n_1* σ̄_1s,g b_g(T) σT⁴/π per volume, with the LTE (Saha)
ground-state density n_1* = n_e n_p λ_e³ e^{χ_H/kT} (the Milne relation, Milne 1924, 1924PMag...47..209M; recombination continua of
nebulae in this form: Brown & Mathews 1970, 1970ApJ...160..939B). Checks (docs §17): a Planck field in LTE is in exact balance with the
absorption of the same group; the photon-number integral of the same cross-section gives α_1 = α_A − α_B.
"""
group_fast(ρ, T, gidx, G; n1, n2, nH2, nHe, nHep = 0.0, nHepp = 0.0, nHp, ne, nHm, hbfP, Tw = T, Tmin = 3000.0, edges = nothing) =
    _group_fast(fastgrid(), ρ, T, gidx, G, n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, hbfP, Tw, Float64(Tmin), edges)
@inline _Bf(ν, t) = (u = COp.h * ν / (COp.kB * t); u > 700 ? zero(u) : 2COp.h * ν^3 / COp.c^2 * exp(-u) / (-expm1(-u)))       # (= B_ν; the form e^{-u}/(1 − e^{-u}) keeps the Duals' derivative finite in the Wien tail, where expm1(u)² overflows)
@inline _Df(ν, t) = (u = COp.h * ν / (COp.kB * t); u > 700 ? 0.0 : 2COp.h * ν^3 / COp.c^2 * u * exp(-u) / (t * (-expm1(-u))^2))
"(Σ ∂B/∂T, Σ ∂B/∂T/χ_ν, Σ B/χ_ν, Σ n σ_Ram B/χ_ν) over the cells of an `lwsub` group at the shape temperature T1s (xo: the nodes' smooth extinction)."
function _lwsums(ν::Vector{Float64}, st::Vector{Float64}, srm::Vector{Float64}, nodes::Vector{Int}, ia::Vector{Int}, wa::Vector{Float64}, xo::AbstractVector, ng, T1s)
    sd = zero(ng * T1s * xo[1]); sr = sd; sφ = sd; sφr = sd
    @inbounds for m in eachindex(ν)
        i0 = nodes[ia[m]]; i1 = nodes[ia[m]+1]
        χν = (1 - wa[m]) * xo[i0] + wa[m] * xo[i1] + ng * st[m]; χν > 0 || continue
        Dν = _Df(ν[m], T1s); Bν = _Bf(ν[m], T1s)
        sd += Dν; sr += Dν / χν; sφ += Bν / χν; sφr += ng * srm[m] * Bν / χν
    end
    sd, sr, sφ, sφr
end
"""
a/b with the derivative (a′ − (a/b) b′)/b: ForwardDiff's quotient divides by b², which underflows for the Planck weights of the Wien tail (1e-190: χ = ΣD/Σ(D/ext) of the ionising groups in cold gas).
Floats: a/b.
"""
@inline _qd(a, b) = a / b
@inline function _qd(a::ForwardDiff.Dual{T,V,N}, b::ForwardDiff.Dual{T,V,N}) where {T,V,N}
    bv = ForwardDiff.value(b); v = ForwardDiff.value(a) / bv
    ForwardDiff.Dual{T}(v, (ForwardDiff.partials(a) - v * ForwardDiff.partials(b)) / bv)
end
@inline _qd(a::ForwardDiff.Dual{T,V,N}, b::Real) where {T,V,N} = ForwardDiff.Dual{T}(ForwardDiff.value(a) / b, ForwardDiff.partials(a) / b)
@inline _qd(a::Real, b::ForwardDiff.Dual{T,V,N}) where {T,V,N} = (bv = ForwardDiff.value(b); v = a / bv; ForwardDiff.Dual{T}(v, -v * ForwardDiff.partials(b) / bv))
"Task-local work vectors of `_group_fast` (three of the grid's length): the evaluation runs in the threaded zone sweeps, where per-call vectors cost more in GC than the arithmetic."
_g2scratch(::Type{Float64}, N) = get!(() -> (zeros(N), zeros(N), zeros(N)), task_local_storage(), :G2O_group_scratch)::NTuple{3,Vector{Float64}}
_g2scratch(::Type{Tp}, N) where {Tp} = (zeros(Tp, N), zeros(Tp, N), zeros(Tp, N))
"""
`group_fast` is generic in the number type (ForwardDiff Duals in any of ρ, T, T_w and the populations): the means, the group means and their derivatives are
exact derivatives of the Float64 evaluation (the node table is Float64; the node loop, the group sums and the quotients carry the Duals).
"""
function _group_fast(F, ρ, T, gidx::Vector{Int}, G::Int, n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, hbfP::Bool, Tw, Tmin::Float64, edges)     # (function barrier: F concrete)
    Tp = typeof(one(ρ) * one(T) * one(n1) * one(n2) * one(nH2) * one(nHe) * one(nHep) * one(nHepp) * one(nHp) * one(ne) * one(nHm) * one(Tw))
    kB = COp.kB; h = COp.h; kT = kB * T
    Ts = T; Tws = Tw; T1s = max(Tw, Tmin)
    Nion = nHp + nHep + nHepp; Nneu = n1 + nHe + nH2
    χH = Z0.χH; χ2 = Z0.χHe2
    lwH = ntuple(n -> Z0.ln_occupation(n, 1, Nion, Nneu), NLEV)       # (as means_fast: the same populations, so Σ_g b_g κ_g is the grey mean)
    wH = [exp(lwH[n]) for n in 1:NLEV]; wHe = [Z0.occupation(n, 2, Nion, Nneu) for n in 1:NLEV]
    pH = zeros(Tp, NLEV); pHe = zeros(Tp, NLEV); sx = zero(Tp)
    for n in 2:NLEV; pH[n] = n2 * (n^2 / 4) * exp(-χH * (0.25 - 1 / n^2) / kT) * wH[n] / max(wH[2], 1e-300); sx += pH[n]; end
    ng = max(n1 - sx, zero(Tp)); pH[1] = ng
    UHe = zero(Tp)
    for n in 1:NLEV; pHe[n] = n^2 * exp(-χ2 * (1 - 1 / n^2) / kT) * wHe[n]; UHe += pHe[n]; end
    for n in 1:NLEV; pHe[n] *= nHep / UHe; end
    x = clamp((log10(T) - FLT0) / FDLT, 0.0, FNT - 1 - 1e-9); j = floor(Int, ForwardDiff_value(x)) + 1; a = clamp(x - (j - 1), 0.0, 1.0)
    wJ = Tw != T
    N = length(F.ν); dis, dis1, xo = _g2scratch(Tp, N); fill!(dis, zero(Tp)); fill!(dis1, zero(Tp))
    for (i, n, ns, σe) in F.DIS
        d = 1 - exp(Z0.ln_occupation(ns, 1, Nion, Nneu) - lwH[n])
        d > 0 && (dis[i] += d * pH[n] * σe; n == 1 && (dis1[i] += d * pH[n] * σe))
    end
    sub = edges === nothing ? nothing : lwsub(edges, gidx, G)                # (xo: every node is set in the loop below when sub is given)
    gE = zeros(Tp, G); gW = zeros(Tp, G); gP = zeros(Tp, G); gB = zeros(Tp, G); gR = zeros(Tp, G); gD = zeros(Tp, G); g1 = zeros(Tp, G); g1W = zeros(Tp, G); g2 = zeros(Tp, G)
    # the same with the floored shape temperature, used for a group whose actual weight underflows (cold gas, high groups): its Planck
    # fraction is then at the floor, so emission does not matter, but absorption and extinction of the light it carries (the star's) do —
    # a grey fallback gave the LW band the cold gas's H₂ CIA κ_E and absorbed it within a few zones (t_mg34)
    fE = zeros(Tp, G); fW = zeros(Tp, G); fP = zeros(Tp, G); fB = zeros(Tp, G); fR = zeros(Tp, G); fD = zeros(Tp, G); Tf = max(T, Tmin)
    g1P = zeros(Tp, G); f1P = zeros(Tp, G)                                   # H(1s) bound–free with the emission weight B_ν(T) (Milne emission)
    sE = zero(Tp); sP = zero(Tp); sB = zero(Tp); sR = zero(Tp); sD = zero(Tp); sW = zero(Tp)
    @inbounds for i in 1:N
        ν = F.ν[i]; dν = F.w[i]; u = h * ν / kT; em = exp(-u); stim = 1 - em
        B = 2h * ν^3 / COp.c^2 * em / stim; D = B * u / (T * stim)
        ab = dis[i]; a1 = pH[1] * F.SH[1, i]
        for n in 1:NLEV; s_ = F.SH[n, i]; s_ > 0 && (ab += pH[n] * s_); end
        ab *= stim; a1 *= stim; a2 = ab - a1 - dis1[i] * stim               # (a2: H n ≥ 2 bound–free and their dissolved levels; the ground
                                                                             #  state's dissolved Lyman pseudo-continuum, a sliver just below
                                                                             #  13.6 eV, would dominate an arithmetic LW-band mean)
        ahe = nHe * F.σHeI[i]
        if nHep > 0
            for n in 1:NLEV; s_ = F.SHe[n, i]; s_ > 0 && (ahe += pHe[n] * s_); end
        end
        ahe *= stim
        aff = _tint(F.HF, j, a, i) * ng * ne * kT +
              3.692e8 * ne / (sqrt(T) * F.ν3[i]) * stim * ((nHp + nHep) * _tint(F.G1, j, a, i) + 4nHepp * _tint(F.G2, j, a, i))
        acia = _tint(F.C1, j, a, i) * nH2^2 + _tint(F.C2, j, a, i) * nH2 * nHe + _tint(F.C3, j, a, i) * nH2 * n1
        a0 = acia + aff + ahe
        σhm = F.σhm[i] * stim
        sca = F.ray2[i] * nH2 + F.ray1[i] * ng + COp.σT * ne
        hb = hbfP ? ab : 0.0
        BJ = wJ ? _Bf(ν, Tw) : B
        aE = a0 + σhm * nHm + hb; aP = a0 + σhm * nHm + hb; ext = a0 + σhm * nHm + ab + sca
        sE += aE * BJ * dν; sW += BJ * dν; sP += aP * B * dν; sB += B * dν
        sR += D / max(ext, 1e-300) * dν; sD += D * dν
        g = gidx[i]; sub === nothing || (xo[i] = ext - F.ray1[i] * ng); Bs = _Bf(ν, Ts); Bw = _Bf(ν, Tws); Ds = _Df(ν, Ts)
        B1 = _Bf(ν, T1s); Bfl = _Bf(ν, Tf); Dfl = _Df(ν, Tf)
        fE[g] += aE * B1 * dν; fW[g] += B1 * dν; fP[g] += aP * Bfl * dν; fB[g] += Bfl * dν; fR[g] += Dfl / max(ext, 1e-300) * dν; fD[g] += Dfl * dν
        gE[g] += aE * Bw * dν; gW[g] += Bw * dν; g1[g] += a1 * B1 * dν; g1W[g] += B1 * dν; g2[g] += a2 * B1 * dν
        gP[g] += aP * Bs * dν; gB[g] += Bs * dν; s1 = F.SH[1, i] * stim; g1P[g] += s1 * Bs * dν; f1P[g] += s1 * Bfl * dν
        gR[g] += Ds / max(ext, 1e-300) * dν; gD[g] += Ds * dν
    end
    q = (a_, b_) -> b_ > 0 ? _qd(a_, b_) / ρ : zero(Tp)
    uf = w -> !(w > 1e-280)                                                  # (the group's actual weight underflowed)
    kram = zeros(Tp, G); χs = zeros(Tp, G)
    if sub !== nothing
        for g in 1:G
            S = sub[g]; S === nothing && continue
            sd, sr, sφ, sφr = _lwsums(S.ν, S.st, S.sr, S.nodes, S.ia, S.wa, xo, ng, T1s)
            sr > 0 && (χs[g] = _qd(sd, sr) / ρ); sφ > 0 && (kram[g] = _qd(sφr, sφ) / ρ)
        end
    end
    ([uf(gW[g]) ? q(fE[g], fW[g]) : q(gE[g], gW[g]) for g in 1:G], [uf(gB[g]) ? q(fP[g], fB[g]) : q(gP[g], gB[g]) for g in 1:G],
     [χs[g] > 0 ? χs[g] : uf(gD[g]) ? (fR[g] > 0 ? _qd(fD[g], fR[g]) / ρ : zero(Tp)) : (gR[g] > 0 ? _qd(gD[g], gR[g]) / ρ : zero(Tp)) for g in 1:G],
     [q(g1[g], g1W[g]) for g in 1:G], _qd(sE, sW) / ρ, _qd(sP, sB) / ρ, _qd(sD, sR) / ρ, [q(g2[g], g1W[g]) for g in 1:G], kram,
     [uf(gB[g]) ? (fB[g] > 0 ? _qd(f1P[g], fB[g]) : zero(Tp)) : _qd(g1P[g], gB[g]) for g in 1:G])
end
"Group index of every fast node for interior edges `edges` [Hz] (ascending)."
group_index(edges) = (F = fastgrid(); [searchsortedlast(edges, ν) + 1 for ν in F.ν])
ForwardDiff_value(x::ForwardDiff.Dual) = ForwardDiff_value(ForwardDiff.value(x))
ForwardDiff_value(x) = x
const FASTON = get(ENV, "G2OPAC", "0") == "2"
"`means_fast` with G2OPAC=2, else `means`."
means_any(ρ, T; kw...) = FASTON ? means_fast(ρ, T; kw...) : means(ρ, T; kw...)

end # module
