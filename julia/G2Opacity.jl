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
    means(ρ, T; n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, nHmS, hbfP, Tw = T) -> (κ_E, κ_P, κ_R) [cm² g⁻¹]

As PopIIIMaterial.nlte_means (n1: H atoms, taken as the ground state; n2: population of n = 2 [cm⁻³]; nHe: neutral He) with the
additions listed in the module docstring.
"""
function means(ρ, T; n1, n2, nH2, nHe, nHep = 0.0, nHepp = 0.0, nHp, ne, nHm, nHmS, hbfP, Tw = T)
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
    cia = nH2 > 1e-3 * n1; acia = 0.0
    νH = χH / COp.h; νHe = χ2 / COp.h
    @inbounds for i in eachindex(F.ν)
        ν = F.ν[i]; dν = F.dν[i]; σ = F.σ[i]; λμ = F.λμ[i]; ν3 = F.ν3[i]
        stim = -expm1(-COp.h * ν / kT)
        cia && i % 8 == 1 && (acia = COp.k_h2h2(t, σ, T) * nH2^2 + COp.kcia(t.H2He, σ, T) * nH2 * nHe + COp.kcia(t.H2H, σ, T) * nH2 * n1)
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
        sE += (a0 + σhm * nHm + h) * BJ * dν; sW += BJ * dν; sP += (a0 + σhm * nHmS + h) * B * dν; sB += B * dν
        sR += D / max(a0 + σhm * nHm + ab + sca, 1e-300) * dν; sD += D * dν
    end
    (sE / sW / ρ, sP / sB / ρ, sD / sR / ρ)
end

"""
LTE population of H n = 2 [cm⁻³] for `nH0` neutral H atoms (all levels): Boltzmann with occupation probabilities over n ≤ $(NLEV)
(Nion, Nneu as in `means`).
"""
function n2_lte(T, nH0, Nion, Nneu)
    kT = COp.kB * Float64(T); χH = Z0.χH
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
@inline function _tint(A, j, a, i)                                     # linear in log T between rows j, j+1
    (1 - a) * A[j, i] + a * A[j+1, i]
end
"`means` on the fast grid (see FAST); same arguments and conventions, generic in the number type."
means_fast(ρ, T; n1, n2, nH2, nHe, nHep = 0.0, nHepp = 0.0, nHp, ne, nHm, nHmS, hbfP, Tw = T) =
    _means_fast(fastgrid(), ρ, T, n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, nHmS, hbfP, Tw)
function _means_fast(F, ρ, T, n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, nHmS, hbfP, Tw)     # (function barrier: F concrete)
    kB = COp.kB; h = COp.h; kT = kB * T
    Nion = nHp + nHep + nHepp; Nneu = n1 + nHe + nH2
    χH = Z0.χH; χ2 = Z0.χHe2
    Tp = typeof(one(n1) * one(n2) * one(nHep) * one(Nion) * one(Nneu) * one(T))
    wH = Vector{Tp}(undef, NLEV); wHe = Vector{Tp}(undef, NLEV); pH = Vector{Tp}(undef, NLEV); pHe = Vector{Tp}(undef, NLEV)
    for n in 1:NLEV; wH[n] = Z0.occupation(n, 1, Nion, Nneu); wHe[n] = Z0.occupation(n, 2, Nion, Nneu); end
    sx = zero(Tp)
    for n in 2:NLEV; pH[n] = n2 * (n^2 / 4) * exp(-χH * (0.25 - 1 / n^2) / kT) * wH[n] / max(wH[2], 1e-300); sx += pH[n]; end
    ng = max(n1 - sx, zero(Tp)); pH[1] = ng
    UHe = zero(Tp)
    for n in 1:NLEV; pHe[n] = n^2 * exp(-χ2 * (1 - 1 / n^2) / kT) * wHe[n]; UHe += pHe[n]; end
    for n in 1:NLEV; pHe[n] *= nHep / UHe; end
    lT = log10(ForwardDiff_value(T)); x = clamp((lT - FLT0) / FDLT, 0.0, FNT - 1 - 1e-9); j = floor(Int, x) + 1
    a = clamp((log10(T) - FLT0) / FDLT - (j - 1), 0.0, 1.0)            # (generic weight; clamped at the table ends)
    cia = nH2 > 1e-3 * n1; wJ = Tw != T
    N = length(F.ν)
    # dissolved-level contributions per node (sparse)
    dis = zeros(Tp, N)
    for (i, n, ns, σe) in F.DIS
        d = 1 - Z0.occupation(ns, 1, Nion, Nneu) / max(wH[n], 1e-300)
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
        acia = cia ? (_tint(F.C1, j, a, i) * nH2^2 + _tint(F.C2, j, a, i) * nH2 * nHe + _tint(F.C3, j, a, i) * nH2 * n1) : zero(aff)
        a0 = acia + aff + ahe
        σhm = F.σhm[i] * stim
        sca = F.ray2[i] * nH2 + F.ray1[i] * ng + COp.σT * ne
        hb = hbfP ? ab : zero(ab)
        BJ = wJ ? COp.Bν(ν, Tw) : B
        sE += (a0 + σhm * nHm + hb) * BJ * dν; sW += BJ * dν; sP += (a0 + σhm * nHmS + hb) * B * dν; sB += B * dν
        sR += D / max(a0 + σhm * nHm + ab + sca, 1e-300) * dν; sD += D * dν
    end
    (sE / sW / ρ, sP / sB / ρ, sD / sR / ρ)
end
"""
MGROUPS (docs/POPIII_MULTIGROUP.md): per-group means of `means_fast`'s physics on the same nodes and with the same populations, plus the
grey means of the same call (Float64 only). `gidx[i]` is the group of node i. In-group weights: emission B_ν(T), absorption B_ν(Tw),
extinction ∂B_ν/∂T(T) — the weights of the grey means, so that Σ_g b_g(T) κP_g = κP and Σ_g b_g(Tw) κE_g = κE hold at every T (a group whose
weight underflows returns 0 and is taken as grey by the caller; its Planck fraction is negligible there). (A floor of the in-group shape
temperature broke the emission identity below it: 1000 K cloud gas then over-emitted into the optical group by orders of magnitude, t_mg27.)
Returns (κE_g, κP_g, χ_g, κ1s_g, κE, κP, χ): κ1s_g the H(1s) bound–free absorption per mass with the in-group shape B_ν(max(Tw, Tmin)) (the
star's light, also in unlit cold gas where B_ν(Tw) of the ionising groups underflows), whether or not hbfP puts H bound–free into κE_g.
"""
function group_fast(ρ, T, gidx, G; n1, n2, nH2, nHe, nHep = 0.0, nHepp = 0.0, nHp, ne, nHm, nHmS, hbfP, Tw = T, Tmin = 3000.0)
    F = fastgrid(); kB = COp.kB; h = COp.h; ρ = Float64(ρ); T = Float64(T); Tw = Float64(Tw); kT = kB * T
    Ts = T; Tws = Tw; T1s = max(Tw, Tmin)
    Nion = nHp + nHep + nHepp; Nneu = n1 + nHe + nH2
    χH = Z0.χH; χ2 = Z0.χHe2
    wH = [Z0.occupation(n, 1, Nion, Nneu) for n in 1:NLEV]; wHe = [Z0.occupation(n, 2, Nion, Nneu) for n in 1:NLEV]
    pH = zeros(NLEV); pHe = zeros(NLEV); sx = 0.0
    for n in 2:NLEV; pH[n] = n2 * (n^2 / 4) * exp(-χH * (0.25 - 1 / n^2) / kT) * wH[n] / max(wH[2], 1e-300); sx += pH[n]; end
    ng = max(n1 - sx, 0.0); pH[1] = ng
    UHe = 0.0
    for n in 1:NLEV; pHe[n] = n^2 * exp(-χ2 * (1 - 1 / n^2) / kT) * wHe[n]; UHe += pHe[n]; end
    for n in 1:NLEV; pHe[n] *= nHep / UHe; end
    x = clamp((log10(T) - FLT0) / FDLT, 0.0, FNT - 1 - 1e-9); j = floor(Int, x) + 1; a = clamp(x - (j - 1), 0.0, 1.0)
    cia = nH2 > 1e-3 * n1; wJ = Tw != T
    N = length(F.ν); dis = zeros(N); dis1 = zeros(N)
    for (i, n, ns, σe) in F.DIS
        d = 1 - Z0.occupation(ns, 1, Nion, Nneu) / max(wH[n], 1e-300)
        d > 0 && (dis[i] += d * pH[n] * σe; n == 1 && (dis1[i] += d * pH[n] * σe))
    end
    gE = zeros(G); gW = zeros(G); gP = zeros(G); gB = zeros(G); gR = zeros(G); gD = zeros(G); g1 = zeros(G); g1W = zeros(G); g2 = zeros(G)
    # the same with the floored shape temperature, used for a group whose actual weight underflows (cold gas, high groups): its Planck
    # fraction is then at the floor, so emission does not matter, but absorption and extinction of the light it carries (the star's) do —
    # a grey fallback gave the LW band the cold gas's H₂ CIA κ_E and absorbed it within a few zones (t_mg34)
    fE = zeros(G); fW = zeros(G); fP = zeros(G); fB = zeros(G); fR = zeros(G); fD = zeros(G); Tf = max(T, Tmin)
    sE = 0.0; sP = 0.0; sB = 0.0; sR = 0.0; sD = 0.0; sW = 0.0
    Bf(ν, t) = (u = h * ν / (kB * t); u > 700 ? 0.0 : 2h * ν^3 / COp.c^2 / expm1(u))
    Df(ν, t) = (u = h * ν / (kB * t); u > 700 ? 0.0 : 2h * ν^3 / COp.c^2 * u * exp(-u) / (t * (-expm1(-u))^2))
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
        acia = cia ? (_tint(F.C1, j, a, i) * nH2^2 + _tint(F.C2, j, a, i) * nH2 * nHe + _tint(F.C3, j, a, i) * nH2 * n1) : 0.0
        a0 = acia + aff + ahe
        σhm = F.σhm[i] * stim
        sca = F.ray2[i] * nH2 + F.ray1[i] * ng + COp.σT * ne
        hb = hbfP ? ab : 0.0
        BJ = wJ ? COp.Bν(ν, Tw) : B
        aE = a0 + σhm * nHm + hb; aP = a0 + σhm * nHmS + hb; ext = a0 + σhm * nHm + ab + sca
        sE += aE * BJ * dν; sW += BJ * dν; sP += aP * B * dν; sB += B * dν
        sR += D / max(ext, 1e-300) * dν; sD += D * dν
        g = gidx[i]; Bs = Bf(ν, Ts); Bw = Bf(ν, Tws); Ds = Df(ν, Ts)
        B1 = Bf(ν, T1s); Bfl = Bf(ν, Tf); Dfl = Df(ν, Tf)
        fE[g] += aE * B1 * dν; fW[g] += B1 * dν; fP[g] += aP * Bfl * dν; fB[g] += Bfl * dν; fR[g] += Dfl / max(ext, 1e-300) * dν; fD[g] += Dfl * dν
        gE[g] += aE * Bw * dν; gW[g] += Bw * dν; g1[g] += a1 * B1 * dν; g1W[g] += B1 * dν; g2[g] += a2 * B1 * dν
        gP[g] += aP * Bs * dν; gB[g] += Bs * dν
        gR[g] += Ds / max(ext, 1e-300) * dν; gD[g] += Ds * dν
    end
    q(a_, b_) = b_ > 0 ? a_ / b_ / ρ : 0.0
    uf(w) = !(w > 1e-280)                                                    # (the group's actual weight underflowed)
    ([uf(gW[g]) ? q(fE[g], fW[g]) : q(gE[g], gW[g]) for g in 1:G], [uf(gB[g]) ? q(fP[g], fB[g]) : q(gP[g], gB[g]) for g in 1:G],
     [uf(gD[g]) ? (fR[g] > 0 ? fD[g] / fR[g] / ρ : 0.0) : (gR[g] > 0 ? gD[g] / gR[g] / ρ : 0.0) for g in 1:G],
     [q(g1[g], g1W[g]) for g in 1:G], sE / sW / ρ, sP / sB / ρ, sD / sR / ρ, [q(g2[g], g1W[g]) for g in 1:G])
end
"Group index of every fast node for interior edges `edges` [Hz] (ascending)."
group_index(edges) = (F = fastgrid(); [searchsortedlast(edges, ν) + 1 for ν in F.ν])
ForwardDiff_value(x::ForwardDiff.Dual) = ForwardDiff_value(ForwardDiff.value(x))
ForwardDiff_value(x) = x
const FASTON = get(ENV, "G2OPAC", "0") == "2"
"`means_fast` with G2OPAC=2, else `means`."
means_any(ρ, T; kw...) = FASTON ? means_fast(ρ, T; kw...) : means(ρ, T; kw...)

end # module
