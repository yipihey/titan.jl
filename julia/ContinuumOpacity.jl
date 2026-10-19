"""
    ContinuumOpacity

Population III project, M3c: frequency-dependent continuum opacity of metal-free gas, κ_ν(ρ, T, composition), and its
Planck and Rosseland means. Sources (data outside the repository, `~/Projects/popiii_lit/`):

- **CIA** (collision-induced absorption), HITRAN (Karman et al. 2019 compilation; files `main/*_2011.cia`):
  H₂–H₂ (200–3000 K, 20–10 000 cm⁻¹; Abel et al. 2011/2012), H₂–He (200–9900 K, 20–20 000 cm⁻¹; Abel et al. 2012),
  H₂–H (1000–2500 K, 100–10 000 cm⁻¹; Gustafsson & Frommhold 2003). α_ν = k(ν, T) n₁ n₂ [cm⁻¹], k in cm⁵.
  Outside the tabulated T range k is held at the nearest table temperature (flagged: H₂–H₂ above 3000 K, where RA04 used
  Borysow et al. 2001 up to 7000 K — not in HITRAN main).
- **H⁻ bound–free**: John (1988, A&A 193, 189) fit to Wishart (1979): σ = 1e-18 λ³ (1/λ − 1/λ₀)^{3/2}
  Σ_{n=1}^{6} C_n (1/λ − 1/λ₀)^{(n−1)/2} cm², λ in μm, λ₀ = 1.6419 μm (coefficients as in the AHMIN routine of the
  "oracle" code and in RH); × n_H⁻ (1 − e^{−hν/kT}).
- **H⁻ free–free**: John (1988) Table 3a/b (coefficients cross-checked in RH `hydrogen.c` and NEMESIS `hmin_ff.f`):
  k = 1e-29 Σ θ^{(n+1)/2}(A_n λ² + B_n + C_n/λ + D_n/λ² + E_n/λ³ + F_n/λ⁴) cm⁴ dyn⁻¹, θ = 5040/T, λ in μm;
  α = k n_H P_e (stimulated emission included in the fit).
- **H (neutral) bound–free and H⁺ + e free–free**: hydrogenic (Kramers, Gaunt factors = 1) — adequate only as an
  order of magnitude; matters above ≈ 5000 K.
- **Rayleigh**: H₂ — RH (Victor & Dalgarno 1969; Tarafdar & Vardya 1973 fit) σ[Mb] = (87.79 + (1.323e6 +
  2.245e10/λ²)/λ²)/λ⁴, λ in nm; H — σ = 5.799e-13/λ⁴ + 1.422e-6/λ⁶ + 2.784/λ⁸ cm², λ in Å (Dalgarno 1962) [memory,
  unverified].
- **Thomson scattering** σ_T n_e.

Scattering enters the Rosseland mean (extinction), not the Planck mean (absorption only).
"""
module ContinuumOpacity

const h = 6.62607015e-27; const c = 2.99792458e10; const kB = 1.380649e-16; const me = 9.1093837e-28
const σT = 6.6524587e-25; const eV = 1.602176634e-12

# ------------------------------------------------------------------ CIA tables
struct CIA
    T::Vector{Float64}                 # temperatures of the blocks
    ν::Vector{Vector{Float64}}         # wavenumbers [cm⁻¹] per block
    k::Vector{Vector{Float64}}         # cm⁵ per block
end
function read_cia(path)
    lines = readlines(path); i = 1; Ts = Float64[]; νs = Vector{Float64}[]; ks = Vector{Float64}[]
    while i <= length(lines)
        hdr = split(lines[i]); i += 1
        isempty(hdr) && continue
        n = parse(Int, hdr[4]); T = parse(Float64, hdr[5])
        ν = zeros(n); k = zeros(n)
        for j in 1:n
            v = split(lines[i]); i += 1
            ν[j] = parse(Float64, v[1]); k[j] = parse(Float64, v[2])
        end
        push!(Ts, T); push!(νs, ν); push!(ks, k)
    end
    # several blocks can share a temperature (different spectral ranges or sources): keep them; sort by T
    p = sortperm(Ts); CIA(Ts[p], νs[p], ks[p])
end
"k(σ̃, T) [cm⁵]: on each block, linear in wavenumber (0 outside its range); in T, linear between the two nearest blocks
that cover σ̃ (held at the nearest one outside)."
function kcia(t::CIA, σ, T)
    lo = 0; hi = 0; klo = 0.0; khi = 0.0
    for b in eachindex(t.T)
        ν = t.ν[b]
        (σ < ν[1] || σ > ν[end]) && continue
        j = clamp(searchsortedfirst(ν, σ), 2, length(ν)); w = (σ - ν[j-1]) / (ν[j] - ν[j-1])
        kv = max((1 - w) * t.k[b][j-1] + w * t.k[b][j], 0.0)
        if t.T[b] <= T; lo = b; klo = kv; else; (hi == 0) && (hi = b; khi = kv); end
    end
    lo == 0 && hi == 0 && return 0.0
    lo == 0 && return khi
    hi == 0 && return klo
    w = (T - t.T[lo]) / (t.T[hi] - t.T[lo]); (1 - w) * klo + w * khi
end

const DIR = joinpath(homedir(), "Projects", "popiii_lit", "CIA")
const TABS = Ref{Any}(nothing)
tables() = (TABS[] === nothing && (TABS[] = (H2H2 = read_cia(joinpath(DIR, "H2-H2_2011.cia")),
    H2He = read_cia(joinpath(DIR, "H2-He_2011.cia")), H2H = read_cia(joinpath(DIR, "H2-H_2011.cia")))); TABS[])

# ------------------------------------------------------------------ H⁻ and hydrogen
const Cbf = (152.519, 49.534, -118.858, 92.536, -34.194, 4.982)
"H⁻ photodetachment cross section [cm²] (John 1988), λ in μm."
function σ_Hm_bf(λ)
    (λ >= 1.6419 || λ < 0.125) && return 0.0
    x = 1 / λ - 1 / 1.6419; s = 0.0
    for n in 6:-1:1; s = Cbf[n] + s * sqrt(x); end
    1e-18 * λ^3 * x^1.5 * s
end
const Aff1 = (0.0, 2483.346, -3449.889, 2200.040, -696.271, 88.283)
const Bff1 = (0.0, 285.827, -1158.382, 2427.719, -1841.400, 444.517)
const Cff1 = (0.0, -2054.291, 8746.523, -13651.105, 8624.970, -1863.864)
const Dff1 = (0.0, 2827.776, -11485.632, 16755.524, -10051.530, 2095.288)
const Eff1 = (0.0, -1341.537, 5303.609, -7510.494, 4400.067, -901.788)
const Fff1 = (0.0, 208.952, -812.939, 1132.738, -655.020, 132.985)
const Aff2 = (518.1021, 473.2636, -482.2089, 115.5291); const Bff2 = (-734.8666, 1443.4137, -737.1616, 169.6374)
const Cff2 = (1021.1775, -1977.3395, 1096.8827, -245.649); const Dff2 = (-479.0721, 922.3575, -521.1341, 114.243)
const Eff2 = (93.1373, -178.9275, 101.7963, -21.9972); const Fff2 = (-6.4285, 12.3600, -7.0571, 1.5097)
"H⁻ free–free coefficient [cm⁴ dyn⁻¹] per H atom per unit electron pressure (John 1988), λ in μm."
function k_Hm_ff(λ, T)
    θ = 5040 / T; s = 0.0
    if λ > 0.3645
        for n in 1:6
            s += θ^((n + 1) / 2) * (λ^2 * Aff1[n] + Bff1[n] + Cff1[n] / λ + Dff1[n] / λ^2 + Eff1[n] / λ^3 + Fff1[n] / λ^4)
        end
    elseif λ > 0.1823
        for n in 1:4
            s += θ^((n + 1) / 2) * (λ^2 * Aff2[n] + Bff2[n] + Cff2[n] / λ + Dff2[n] / λ^2 + Eff2[n] / λ^3 + Fff2[n] / λ^4)
        end
    end
    max(1e-29 * s, 0.0)
end
"Hydrogenic bound–free of neutral H from all levels n with hν ≥ χ/n² (Kramers, g = 1, LTE level populations) and
H⁺ + e free–free (Kramers): absorption coefficients [cm⁻¹], stimulated emission included."
function α_H_kramers(ν, T, nH, nHp, ne)
    χ = 13.598eV; σ0 = 7.906e-18                      # Kramers at threshold for n = 1, scales n (ν_n/ν)³
    α = 0.0; Z = 0.0
    for n in 1:30
        Z += 2n^2 * exp(-χ * (1 - 1 / n^2) / (kB * T))
    end
    for n in 1:30
        h * ν >= χ / n^2 || continue
        pn = 2n^2 * exp(-χ * (1 - 1 / n^2) / (kB * T)) / Z
        α += nH * pn * σ0 * n * (χ / n^2 / (h * ν))^3
    end
    αff = 3.692e8 * ne * nHp / (sqrt(T) * ν^3)          # Kramers free–free, Z = 1, g = 1 [cm⁻¹]
    (α + αff) * (-expm1(-h * ν / (kB * T)))
end
"Rayleigh cross sections [cm²]: H₂ (RH fit, λ in nm → σ in Mb) and H (Dalgarno 1962 [memory])."
σ_ray_H2(λnm) = (l2 = 1 / λnm^2; (87.79 + (1.323e6 + 2.245e10 * l2) * l2) * l2^2 * 1e-18)
σ_ray_H(λÅ) = 5.799e-13 / λÅ^4 + 1.422e-6 / λÅ^6 + 2.784 / λÅ^8

"""
    alpha(ν, T; nH, nH2, nHe, nHp, ne, nHm) -> (α_abs, α_sca, parts) [cm⁻¹]

Continuum absorption (true absorption, stimulated emission included) and scattering coefficients at frequency ν [Hz].
`parts` is a named tuple of the absorption contributions (for diagnostics).
"""
function alpha(ν, T; nH, nH2, nHe, nHp, ne, nHm)
    t = tables(); σ = ν / c; λμ = 1e4 / σ; stim = -expm1(-h * ν / (kB * T))
    a_h2h2 = kcia(t.H2H2, σ, T) * nH2^2
    a_h2he = kcia(t.H2He, σ, T) * nH2 * nHe
    a_h2h = kcia(t.H2H, σ, T) * nH2 * nH
    a_hmbf = σ_Hm_bf(λμ) * nHm * stim
    a_hmff = k_Hm_ff(λμ, T) * nH * ne * kB * T
    a_h = α_H_kramers(ν, T, nH, nHp, ne)
    sca = σ_ray_H2(λμ * 1e3) * nH2 + σ_ray_H(λμ * 1e4) * nH + σT * ne
    abs_ = a_h2h2 + a_h2he + a_h2h + a_hmbf + a_hmff + a_h
    (abs_, sca, (cia_h2h2 = a_h2h2, cia_h2he = a_h2he, cia_h2h = a_h2h, hm_bf = a_hmbf, hm_ff = a_hmff, h = a_h))
end

Bν(ν, T) = 2h * ν^3 / c^2 / expm1(h * ν / (kB * T))
dBdT(ν, T) = (x = h * ν / (kB * T); em = exp(-x); 2h * ν^3 / c^2 * x * em / (T * (1 - em)^2))   # overflow-safe for x ≫ 1

"Frequency grid [Hz] for means: log-uniform in wavenumber from 10 to 1e5 cm⁻¹ (resolves the CIA bands, H⁻ edge)."
grid(; n = 4000) = [c * 10^x for x in range(1.0, 5.0, length = n)]

"""
    means(ρ, T; nH, nH2, nHe, nHp, ne, nHm, νs = grid()) -> (; κP, κR, κR_abs) [cm² g⁻¹]

Planck mean of the absorption, Rosseland mean of the extinction (and of the absorption alone), trapezoid on νs.
"""
function means(ρ, T; nH, nH2, nHe, nHp, ne, nHm, νs = grid())
    sP = 0.0; sB = 0.0; sR = 0.0; sRa = 0.0; sD = 0.0
    for i in 1:length(νs)-1
        ν = 0.5 * (νs[i] + νs[i+1]); dν = νs[i+1] - νs[i]
        a, s, _ = alpha(ν, T; nH, nH2, nHe, nHp, ne, nHm)
        B = Bν(ν, T); D = dBdT(ν, T)
        sP += a * B * dν; sB += B * dν
        sR += D / max(a + s, 1e-300) * dν; sRa += D / max(a, 1e-300) * dν; sD += D * dν
    end
    (; κP = sP / sB / ρ, κR = sD / sR / ρ, κR_abs = sD / sRa / ρ)
end

end # module
