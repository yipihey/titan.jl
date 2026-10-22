"""
    ContinuumFF

Free–free (and the H₂⁺ bound–free) continuum absorption of metal-free gas that `ContinuumOpacity`
lacked, with their data tables:

- **Thermally averaged free–free Gaunt factors** ḡ_ff(γ², u), γ² = Z² Ry/kT, u = hν/kT, bilinear in
  (log γ², log u) on the 81 × 146 grid (log γ² from −6, log u from −16, step 0.2 dex): the merged
  relativistic tables of van Hoof, Ferland, Williams, Volk, Chatzikos, Lykins & Porter (2015, MNRAS 449,
  2112) for Z = 1, 2 (`data/gauntff_merged_Z01.dat`, `_Z02.dat`; non-relativistic below, relativistic
  above ~1e8 K, valid to 1e10 K) by default, and the non-relativistic table of van Hoof, Williams, Volk,
  Chatzikos, Ferland, Abel & Stancil (2014, MNRAS 444, 420) (`data/vanHoof2014-nr-gauntff.dat`) with
  `relativistic = false`.  Frequency-integrated ⟨g_ff⟩(γ²) (`data/gauntff_freqint.dat`, same authors) for
  the total emission.  All files from data.nublado.org/gauntff (fetched 21 Oct 2026; BSD-style licence
  inside) [read-src].  Replaces g = 1 in the hydrogenic free–free of H⁺, He⁺ (Z = 1) and He²⁺ (Z = 2).
  ⟨g_ff⟩ agrees with the Cen (1992) integrated fit used by the cooling coefficients to ≤ 1.7 % for Z = 1
  and 2 over 3e3–3e8 K [run], so that fit is kept there.
- **He⁻ free–free** (e⁻ + He + γ): John (1994, MNRAS 269, 871), K(λ, θ) [1e-26 cm⁴ dyn⁻¹] per
  neutral He atom per unit electron pressure, stimulated emission included; λ 5063–151 878 Å,
  θ = 5040/T = 0.5–3.6 (1400–10 080 K).  Transcription: Korg.jl `absorption_He.jl` (OCR of the paper)
  [read-src-code].  Outside the table: above 10 080 K extrapolated linearly in θ from the θ = 0.5, 0.6
  columns to θ = 0.25 (20 160 K) and held constant beyond (He⁻ ff is then negligible against the H⁺ and
  He⁺ free–free anyway); below 1400 K clamped; ∝ λ² beyond the red edge; 0 below 5063 Å.
- **H₂⁻ free–free** (e⁻ + H₂ + γ): John (1975, MNRAS 172, 305) as tabulated in Turbospectrum's
  `jonabs_vac_v07.3.dat` component 22 ("extr long wav BPz961002": the 210 μm point is Plez's
  extrapolation), K(λ, T) [cm⁴ dyn⁻¹] per H₂ molecule per unit electron pressure, λ 5001–2.1e6 Å,
  T 1000–15 000 K, log–log interpolation as the ILOGL/ILOGT flags say [read-src-code].  Turbospectrum
  applies the stimulated-emission factor to this component afterwards (it divides it out only for the
  He⁻/C⁻/N⁻/O⁻ tables), so the table EXCLUDES stimulated emission and `α_H2m_ff` multiplies by
  (1 − e^{−hν/kT}) [inferred from detabs.f].  Bell (1980, J. Phys. B 13, 1859) supersedes John 1975
  by ~10–20 %; its table was not available here — flagged.
- **H₂⁺ bound–free (photodissociation) and free–free** (H + H⁺ + γ): Stancil (1994, ApJ 430, 360),
  σ_bf [cm²] per H₂⁺ molecule (n_H₂⁺ = n_H n_H⁺/K in LTE, K tabulated) and the ff "cross section"
  per n_H n_H⁺ [cm⁵], λ 500/700–20 000 Å,
  T 3150–25 200 K.  Transcription: Korg.jl `Stancil1994.jl` [read-src-code].  With an H₂⁺ density from
  the network the bf term uses it directly; clamped in T; ff ∝ λ² beyond 2 μm; bf 0 beyond 2 μm.

All coefficients return absorption INCLUDING stimulated emission.
"""
module ContinuumFF

const h = 6.62607015e-27; const c = 2.99792458e10; const kB = 1.380649e-16
const Ry = 13.605693e0 * 1.602176634e-12          # [erg]
const DATA = joinpath(@__DIR__, "data")

# ------------------------------------------------------------------ van Hoof et al. 2014 Gaunt factors
struct GauntTable
    lg2_0::Float64; lu_0::Float64; d::Float64; ng::Int; nu::Int
    g::Matrix{Float64}            # [iu, ig]
end
"Read a van Hoof Gaunt-factor file: 2014 non-relativistic (magic 20140210) or 2015 merged (20140510, one extra Z line)."
function read_gaunt(path = joinpath(DATA, "vanHoof2014-nr-gauntff.dat"))
    lines = [l for l in readlines(path) if !startswith(strip(l), "#") && !isempty(strip(l))]
    hdr(i) = split(lines[i])
    magic = parse(Int, hdr(1)[1]); @assert magic in (20140210, 20140510)
    ng, nu = parse(Int, hdr(2)[1]), parse(Int, hdr(2)[2])
    lg2_0 = parse(Float64, hdr(3)[1]); lu_0 = parse(Float64, hdr(4)[1]); d = parse(Float64, hdr(5)[1])
    off = magic == 20140510 ? 6 : 5                    # merged files carry the atomic number Z
    g = zeros(nu, ng)
    for i in 1:nu
        v = parse.(Float64, split(lines[off + i])); @assert length(v) == ng
        g[i, :] = v
    end
    GauntTable(lg2_0, lu_0, d, ng, nu, g)
end
const GAUNT = Ref{Any}(nothing)                       # 2014 non-relativistic
const GAUNT_REL = Dict{Int,Any}()                     # 2015 merged relativistic, per Z
gaunt_table() = (GAUNT[] === nothing && (GAUNT[] = read_gaunt()); GAUNT[]::GauntTable)
gaunt_table_rel(Z) = get!(GAUNT_REL, Z) do
    read_gaunt(joinpath(DATA, "gauntff_merged_Z0$(Z).dat"))
end::GauntTable
"""
    gaunt_ff(ν, T, Z; relativistic = true) -> ḡ_ff

Thermally averaged free–free Gaunt factor, bilinear in (log γ², log u) with γ² = Z² Ry/kT and
u = hν/kT; clamped at the table edges (log γ² ∈ [−6, 10], log u ∈ [−16, 13]).  Default: the merged
relativistic tables of van Hoof et al. (2015) for Z = 1, 2; `relativistic = false` (or other Z): the
non-relativistic 2014 table.  The two agree to < 0.5 % below 1e7 K.
"""
function gaunt_ff(ν, T, Z; relativistic::Bool = true)
    t = (relativistic && Z in (1, 2)) ? gaunt_table_rel(Z) : gaunt_table()
    x = (log10(Z^2 * Ry / (kB * T)) - t.lg2_0) / t.d + 1
    y = (log10(h * ν / (kB * T)) - t.lu_0) / t.d + 1
    x = clamp(x, 1.0, t.ng - 1e-9); y = clamp(y, 1.0, t.nu - 1e-9)
    i = min(floor(Int, x), t.ng - 1); j = min(floor(Int, y), t.nu - 1); a = x - i; b = y - j
    (1 - a) * (1 - b) * t.g[j, i] + a * (1 - b) * t.g[j, i + 1] + (1 - a) * b * t.g[j + 1, i] + a * b * t.g[j + 1, i + 1]
end
"""
    α_ff_hydrogenic(ν, T, Z, nion, ne) [cm⁻¹]

Kramers free–free with the van Hoof Gaunt factor, stimulated emission included:
3.6919e8 Z² ḡ n_e n_ion T^{−1/2} ν^{−3} (1 − e^{−hν/kT}) (Rybicki & Lightman eq. 5.18b).
"""
α_ff_hydrogenic(ν, T, Z, nion, ne) = 3.6919e8 * Z^2 * gaunt_ff(ν, T, Z) * ne * nion / (sqrt(T) * ν^3) * (-expm1(-h * ν / (kB * T)))

# frequency-integrated ⟨g_ff⟩(γ²) (van Hoof et al. 2014, gauntff_freqint.dat): the Gaunt factor of the TOTAL free–free emission
const GAUNT_INT = Ref{Any}(nothing)
function gaunt_int_table()
    if GAUNT_INT[] === nothing
        lines = [l for l in readlines(joinpath(DATA, "gauntff_freqint.dat")) if !startswith(strip(l), "#") && !isempty(strip(l))]
        @assert parse(Int, split(lines[1])[1]) == 20140424
        rows = [parse.(Float64, split(l)) for l in lines[5:end]]
        GAUNT_INT[] = ([r[1] for r in rows], [r[2] for r in rows])
    end
    GAUNT_INT[]::Tuple{Vector{Float64},Vector{Float64}}
end
"""
    gaunt_ff_total(T, Z) -> ⟨g_ff⟩

Frequency-integrated thermally averaged free–free Gaunt factor (van Hoof et al. 2014), linear in log γ²
= log(Z² Ry/kT), clamped.  Total emission ε_ff = 1.426e-27 Z² √T ⟨g_ff⟩ n_e n_ion erg cm⁻³ s⁻¹; the
Cen (1992) fit 1.1 + 0.34 e^{−(5.5 − log T/Z²)²/3} of the cooling coefficients matches it to ≤ 1.7 %.
"""
function gaunt_ff_total(T, Z)
    lg2, g = gaunt_int_table(); x = clamp(log10(Z^2 * Ry / (kB * T)), lg2[1], lg2[end])
    i = clamp(searchsortedlast(lg2, x), 1, length(lg2) - 1); a = (x - lg2[i]) / (lg2[i + 1] - lg2[i])
    (1 - a) * g[i] + a * g[i + 1]
end

"Bilinear interpolation on an (x, y) table (vectors ascending); clamped in both; `logx`/`logy`: interpolate in log."
function interp2(xs, ys, M, x, y; logx = false, logy = false)
    fx = logx ? log : identity; fy = logy ? log : identity
    i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1); j = clamp(searchsortedlast(ys, y), 1, length(ys) - 1)
    a = clamp((fx(x) - fx(xs[i])) / (fx(xs[i + 1]) - fx(xs[i])), 0.0, 1.0)
    b = clamp((fy(y) - fy(ys[j])) / (fy(ys[j + 1]) - fy(ys[j])), 0.0, 1.0)
    (1 - a) * (1 - b) * M[i, j] + a * (1 - b) * M[i + 1, j] + (1 - a) * b * M[i, j + 1] + a * b * M[i + 1, j + 1]
end

# ------------------------------------------------------------------ He⁻ free–free (John 1994)
const HEM_THETA = [0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6]
const HEM_LAM = [5063.0, 5695.0, 6509.0, 7594.0, 9113.0, 11391.0, 15188.0, 18225.0, 22782.0, 30376.0, 36451.0, 45564.0, 60751.0, 91127.0, 113900.0, 151878.0]                 # [Å]
const HEM_K = [0.033 0.036 0.043 0.049 0.055 0.061 0.066 0.072 0.078 0.1 0.121;
    0.041 0.045 0.053 0.061 0.067 0.074 0.081 0.087 0.094 0.12 0.145;
    0.053 0.059 0.069 0.077 0.086 0.094 0.102 0.109 0.117 0.148 0.178;
    0.072 0.079 0.092 0.103 0.114 0.124 0.133 0.143 0.152 0.19 0.227;
    0.102 0.113 0.131 0.147 0.16 0.173 0.186 0.198 0.21 0.258 0.305;
    0.159 0.176 0.204 0.227 0.247 0.266 0.283 0.3 0.316 0.38 0.444;
    0.282 0.311 0.36 0.4 0.435 0.466 0.495 0.522 0.547 0.643 0.737;
    0.405 0.447 0.518 0.576 0.625 0.67 0.71 0.747 0.782 0.91 1.03;
    0.632 0.698 0.808 0.899 0.977 1.045 1.108 1.165 1.218 1.405 1.574;
    1.121 1.239 1.435 1.597 1.737 1.86 1.971 2.073 2.167 2.49 2.765;
    1.614 1.783 2.065 2.299 2.502 2.681 2.842 2.99 3.126 3.592 3.979;
    2.52 2.784 3.226 3.593 3.91 4.193 4.448 4.681 4.897 5.632 6.234;
    4.479 4.947 5.733 6.387 6.955 7.46 7.918 8.338 8.728 10.059 11.147;
    10.074 11.128 12.897 14.372 15.653 16.798 17.838 18.795 19.685 22.747 25.268;
    15.739 17.386 20.151 22.456 24.461 26.252 27.882 29.384 30.782 35.606 39.598;
    27.979 30.907 35.822 39.921 43.488 46.678 49.583 52.262 54.757 63.395 70.58]        # [1e-26 cm⁴ dyn⁻¹], rows λ, columns θ
"""
    k_Hem_ff(λÅ, T) [cm⁴ dyn⁻¹]

He⁻ free–free coefficient per neutral He atom per unit electron pressure (stimulated emission
included): α = k · n_He · n_e k T.  John (1994); θ > 3.6 clamped; θ < 0.5 (T > 10 080 K) extrapolated
linearly in θ from the θ = 0.5, 0.6 columns down to θ = 0.25 and held beyond (positive by construction);
0 below 5063 Å; ∝ λ² beyond 151 878 Å.
"""
function k_Hem_ff(λÅ, T)
    λÅ < HEM_LAM[1] && return 0.0
    θ = 5040 / T; λc = min(λÅ, HEM_LAM[end]); scale = (λÅ / λc)^2
    if θ >= HEM_THETA[1]
        return 1e-26 * interp2(HEM_LAM, HEM_THETA, HEM_K, λc, min(θ, HEM_THETA[end])) * scale
    end
    k1 = interp2(HEM_LAM, HEM_THETA, HEM_K, λc, HEM_THETA[1]); k2 = interp2(HEM_LAM, HEM_THETA, HEM_K, λc, HEM_THETA[2])
    θe = max(θ, 0.25)
    1e-26 * max(k1 + (k2 - k1) * (θe - HEM_THETA[1]) / (HEM_THETA[2] - HEM_THETA[1]), 0.0) * scale
end

# ------------------------------------------------------------------ H₂⁻ free–free (John 1975 via Turbospectrum)
const H2M_LAM = [5001.4, 7502.07, 10002.74, 15004.1, 25006.82, 50013.63, 100027.27, 150040.89, 200054.52, 2100572.5]   # [Å]
const H2M_T = [1000.0, 2500.0, 5000.0, 7500.0, 10000.0, 12500.0, 15000.0]
const H2M_K = [2.1e-27 1.23e-27 8.93e-28 7.48e-28 6.47e-28 5.58e-28 4.82e-28;
    4.19e-27 2.65e-27 2.07e-27 1.76e-27 1.52e-27 1.3e-27 1.12e-27;
    6.95e-27 4.7e-27 3.8e-27 3.22e-27 2.77e-27 2.37e-27 2.03e-27;
    1.46e-26 1.09e-26 8.94e-27 7.5e-27 6.38e-27 5.42e-27 4.63e-27;
    3.92e-26 3.19e-26 2.59e-26 2.14e-26 1.8e-26 1.52e-26 1.3e-26;
    1.62e-25 1.34e-25 1.06e-25 8.69e-26 7.29e-26 6.13e-26 5.21e-26;
    6.79e-25 5.46e-25 4.27e-25 3.48e-25 2.93e-25 2.46e-25 2.09e-25;
    1.528e-24 1.228e-24 9.61e-25 7.83e-25 6.59e-25 5.53e-25 4.7e-25;
    2.71e-24 2.18e-24 1.71e-24 1.39e-24 1.17e-24 9.84e-25 8.36e-25;
    2.994e-22 2.408e-22 1.883e-22 1.535e-22 1.292e-22 1.085e-22 9.217e-23]      # [cm⁴ dyn⁻¹] per H₂ per unit P_e, EXCLUDING stimulated emission; rows λ, columns T
"""
    k_H2m_ff(λÅ, T) [cm⁴ dyn⁻¹]

H₂⁻ free–free coefficient per H₂ molecule per unit electron pressure, without stimulated emission
(John 1975; log–log in λ and T; clamped in T to 1000–15 000 K; 0 below 5001 Å).
"""
function k_H2m_ff(λÅ, T)
    λÅ < H2M_LAM[1] && return 0.0
    interp2(H2M_LAM, H2M_T, H2M_K, min(λÅ, H2M_LAM[end]), clamp(T, H2M_T[1], H2M_T[end]); logx = true, logy = true) *
        (λÅ > H2M_LAM[end] ? (λÅ / H2M_LAM[end])^2 : 1.0)
end
"H₂⁻ free–free absorption [cm⁻¹], stimulated emission included."
α_H2m_ff(ν, T, nH2, ne) = k_H2m_ff(c / ν * 1e8, T) * nH2 * ne * kB * T * (-expm1(-h * ν / (kB * T)))

# ------------------------------------------------------------------ H₂⁺ bound–free and free–free (Stancil 1994)
const H2P_T = [3150.0, 4200.0, 5040.0, 6300.0, 8400.0, 12600.0, 16800.0, 25200.0]
const H2P_K = [9.6e+18, 9.7683e+19, 2.9997e+20, 8.9599e+20, 2.6532e+21, 8.4501e+21, 1.6852999999999999e+22, 4.2895e+22]       # n_H n_H⁺/n_H₂⁺ [cm⁻³] (LTE), at H2P_T
const H2P_LAM_FF = [700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 2950.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 20000.0, 30000.0, 40000.0, 50000.0, 110000.0, 150000.0, 200000.0]   # [Å]
const H2P_LAM_BF = [500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 20000.0, 30000.0, 45000.0, 50000.0, 100000.0, 150000.0, 200000.0]   # [Å]
const H2P_FF = 1e-39 .* [0.0174 0.0154 0.0142 0.013 0.0116 0.01 0.0091 0.00808;
    0.028 0.0246 0.0227 0.0207 0.0184 0.0158 0.0143 0.0126;
    0.0394 0.0346 0.0319 0.029 0.0257 0.022 0.0199 0.0175;
    0.0514 0.0451 0.0416 0.0378 0.0336 0.0287 0.0259 0.0227;
    0.064 0.0562 0.0519 0.0471 0.0418 0.0357 0.0322 0.0283;
    0.077 0.0676 0.0624 0.0567 0.0504 0.0431 0.0389 0.0341;
    0.0903 0.0794 0.0733 0.0666 0.0592 0.0506 0.0457 0.0401;
    0.104 0.0914 0.0843 0.0767 0.0682 0.0584 0.0527 0.0464;
    0.1177 0.1035 0.0956 0.0869 0.0773 0.0663 0.0599 0.0527;
    0.1317 0.1158 0.107 0.0973 0.0866 0.0743 0.0672 0.0592;
    0.1456 0.1281 0.1184 0.1078 0.096 0.0824 0.0746 0.0658;
    0.1597 0.1405 0.1299 0.1183 0.1054 0.0906 0.0821 0.0725;
    0.1737 0.153 0.1414 0.1288 0.1149 0.0988 0.0896 0.0793;
    0.1877 0.1654 0.1529 0.1394 0.1243 0.1071 0.0972 0.0861;
    0.2017 0.1777 0.1644 0.1499 0.1338 0.1154 0.1048 0.093;
    0.2156 0.1901 0.1759 0.1605 0.1433 0.1237 0.1125 0.0998;
    0.2294 0.2023 0.1873 0.1709 0.1527 0.1319 0.1201 0.1068;
    0.2431 0.2145 0.1987 0.1814 0.1622 0.1402 0.1277 0.1137;
    0.2568 0.2266 0.2099 0.1917 0.1716 0.1485 0.1354 0.1206;
    0.2703 0.2387 0.2211 0.202 0.1809 0.1567 0.143 0.1276;
    0.2836 0.2506 0.2322 0.2123 0.1902 0.1649 0.1506 0.1345;
    0.2969 0.2624 0.2433 0.2225 0.1994 0.1731 0.1582 0.1414;
    0.31 0.2741 0.2542 0.2325 0.2086 0.1812 0.1657 0.1483;
    0.3165 0.2799 0.2596 0.2376 0.2131 0.1853 0.1695 0.1518;
    0.323 0.2857 0.265 0.2425 0.2177 0.1893 0.1732 0.1552;
    0.3858 0.3419 0.3176 0.2913 0.2621 0.229 0.2103 0.1894;
    0.4451 0.3952 0.3677 0.3378 0.3048 0.2674 0.2463 0.2228;
    0.5011 0.4457 0.4152 0.3821 0.3456 0.3044 0.2812 0.2555;
    0.5539 0.4935 0.4603 0.4243 0.3848 0.3401 0.315 0.2873;
    0.6511 0.5821 0.5442 0.5032 0.4583 0.4078 0.3795 0.3484;
    0.7388 0.6625 0.6207 0.5756 0.5262 0.471 0.4402 0.4064;
    0.8186 0.7362 0.6911 0.6425 0.5895 0.5303 0.4975 0.4615;
    0.8918 0.8042 0.7563 0.7049 0.6488 0.5864 0.5518 0.5141;
    0.9596 0.8675 0.8173 0.7633 0.7047 0.6395 0.6036 0.5644;
    1.46 1.345 1.283 1.217 1.146 1.068 1.026 0.9815;
    1.805 1.682 1.617 1.547 1.474 1.394 1.351 1.306;
    2.2 2.075 2.009 1.939 1.865 1.786 1.745 1.7;
    2.313 2.188 2.122 2.052 1.979 1.901 1.859 1.816;
    3.18 3.062 3.0 2.936 2.869 2.798 2.761 2.723;
    3.811 3.7 3.642 3.582 3.52 3.456 3.422 3.387;
    4.323 4.218 4.164 4.108 4.05 3.99 3.958 3.926]        # [cm⁵] per n_H n_H⁺, rows λ, columns T
const H2P_BF = 1e-18 .* [7.34e-05 0.000143 0.000204 0.000287 0.000389 0.000497 0.000549 0.000598;
    0.01 0.015 0.0186 0.023 0.0276 0.0319 0.0337 0.0353;
    0.1676 0.1965 0.2105 0.2215 0.2266 0.2246 0.2211 0.2163;
    0.8477 0.8199 0.7797 0.7183 0.6376 0.5477 0.5037 0.4622;
    2.1113 1.8166 1.6135 1.3823 1.1403 0.9157 0.8176 0.7313;
    3.4427 2.8069 2.4213 2.0136 1.6137 1.2616 1.1129 0.9845;
    4.347 3.5155 3.0224 2.507 2.0062 1.5685 1.3846 1.2262;
    4.6981 3.8841 3.3806 2.8402 2.3019 1.8203 1.6147 1.4358;
    4.6169 3.9763 3.5361 3.0358 2.512 2.0239 1.8096 1.6202;
    4.2811 3.884 3.5476 3.1272 2.6535 2.1858 1.9729 1.7809;
    3.8331 3.685 3.466 3.1434 2.7388 2.3082 2.1033 1.9138;
    3.3624 3.4344 3.3301 3.1098 2.784 2.402 2.2106 2.0287;
    2.9167 3.167 3.1662 3.0438 2.7985 2.4707 2.2958 2.1244;
    2.5172 2.9031 2.9909 2.9573 2.7898 2.5175 2.361 2.2019;
    2.1697 2.6538 2.8149 2.86 2.7657 2.5495 2.4128 2.2682;
    1.8726 2.4243 2.6444 2.7573 2.73 2.5682 2.4516 2.3222;
    1.6208 2.216 2.4821 2.6522 2.6843 2.5725 2.4748 2.3597;
    1.4085 2.0289 2.3297 2.5473 2.6314 2.5653 2.4852 2.3839;
    1.2295 1.8617 2.1883 2.4453 2.5755 2.5517 2.4889 2.4012;
    1.0785 1.7129 2.0581 2.3479 2.5188 2.5347 2.4885 2.4144;
    0.9508 1.5804 1.9387 2.2553 2.4622 2.5145 2.484 2.423;
    0.8424 1.4623 1.8293 2.168 2.4063 2.492 2.4762 2.4277;
    0.7501 1.3568 1.7296 2.0865 2.353 2.4701 2.4687 2.4327;
    0.6712 1.2626 1.6392 2.0118 2.3045 2.4525 2.466 2.4433;
    0.6033 1.1783 1.5577 1.9443 2.262 2.4415 2.4707 2.4626;
    0.5448 1.1029 1.4843 1.8834 2.2256 2.4369 2.4825 2.4901;
    0.3477 0.8233 1.2037 1.6442 2.0832 2.4343 2.5567 2.64;
    0.2412 0.6469 1.0069 1.4533 1.9363 2.3659 2.5349 2.6658;
    0.1782 0.5283 0.8605 1.2951 1.792 2.263 2.4604 2.6222;
    0.1389 0.4487 0.7547 1.178 1.6859 2.1963 2.4228 2.6179;
    0.0927 0.3396 0.612 1.0109 1.5255 2.0874 2.3564 2.6031;
    0.0685 0.2749 0.5175 0.8894 1.3923 1.97 2.2587 2.5323;
    0.0542 0.2335 0.4534 0.8033 1.294 1.88 2.1824 2.476;
    0.0448 0.2041 0.4046 0.7334 1.2067 1.7876 2.0941 2.3965;
    0.0382 0.1813 0.3635 0.6701 1.1194 1.6816 1.9827 2.2828;
    0.0159 0.0901 0.1982 0.3951 0.7108 1.1448 1.3953 1.659;
    0.01 0.0596 0.1325 0.2699 0.4954 0.8132 1.0003 1.1999;
    0.00688 0.0425 0.0962 0.1994 0.3723 0.6216 0.771 0.9325;
    0.00641 0.04 0.0908 0.1889 0.354 0.5932 0.7371 0.8932;
    0.00356 0.0229 0.0526 0.111 0.2109 0.3582 0.4479 0.5463;
    0.0025 0.0161 0.0373 0.079 0.1506 0.2567 0.3216 0.3929;
    0.0019 0.0123 0.0286 0.0607 0.1161 0.1982 0.2487 0.3042]        # σ_bf [cm²] per H₂⁺ molecule, rows λ, columns T
"""
LTE equilibrium constant n_H n_H⁺/n_H₂⁺ [cm⁻³] (Stancil 1994): ln K linear in 1/T between the nodes
(Arrhenius; K ∝ e^{−D₀/kT} with D₀(H₂⁺) = 2.65 eV) and extrapolated the same way below 3150 K; clamped above 25 200 K.
"""
function K_H2p(T)
    T >= H2P_T[end] && return H2P_K[end]
    i = clamp(searchsortedlast(H2P_T, T), 1, length(H2P_T) - 1)
    x1 = 1 / H2P_T[i]; x2 = 1 / H2P_T[i + 1]; a = (1 / T - x1) / (x2 - x1)
    exp((1 - a) * log(H2P_K[i]) + a * log(H2P_K[i + 1]))
end
"""
    α_H2p(ν, T, nH, nHp; nH2p = nothing) -> (α_bf, α_ff) [cm⁻¹], stimulated emission included

H₂⁺ photodissociation (bf, σ_bf [cm²] per H₂⁺ molecule) and H + H⁺ free–free (per n_H n_H⁺, [cm⁵])
(Stancil 1994).  Without `nH2p` the bf term assumes LTE H₂⁺, n_H₂⁺ = n_H n_H⁺/K(T); with the network's
n_H₂⁺ it uses that.  Cross sections clamped in T to 3150–25 200 K; ff ∝ λ² beyond 2 μm and 0 below
700 Å; bf 0 outside 500–20 000 Å.
"""
function α_H2p(ν, T, nH, nHp; nH2p = nothing)
    λ = c / ν * 1e8; Tc = clamp(T, H2P_T[1], H2P_T[end]); stim = -expm1(-h * ν / (kB * T))
    bf = 0.0; ff = 0.0
    if H2P_LAM_BF[1] <= λ <= H2P_LAM_BF[end]
        σ = interp2(H2P_LAM_BF, H2P_T, H2P_BF, λ, Tc)              # σ_bf [cm²] per H₂⁺
        bf = σ * (nH2p === nothing ? nH * nHp / K_H2p(T) : nH2p)
    end
    if λ >= H2P_LAM_FF[1]
        λc = min(λ, H2P_LAM_FF[end])
        ff = interp2(H2P_LAM_FF, H2P_T, H2P_FF, λc, Tc) * (λ / λc)^2 * nH * nHp
    end
    (bf * stim, ff * stim)
end

end # module
