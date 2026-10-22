# Z = 0 Rosseland-mean opacity in LTE, computed from frequency-dependent continuum sources:
#   electron scattering (Thomson), H Rayleigh scattering,
#   H bound-free (levels n = 1..NLEV_H, hydrogenic Kramers cross sections, gaunt = 1),
#   He II bound-free (hydrogenic Z = 2, n = 1..NLEV_HE2), He I ground-state bound-free,
#   free-free on H+, He+ (Z=1) and He++ (Z=2) with gaunt = 1,
#   H− bound-free and free-free (John 1988 fits).
# NOT included: bound–bound lines (H, He), H2/H2+/CIA molecules, gaunt factors ≠ 1.
# Lines raise κ_R by up to a factor of a few around 1e4–1e5 K in OPAL-type tables, so treat
# the result as a lower bound there. Molecules matter only below ~3000–3500 K.

const NLEV_H   = 12
const NLEV_HE2 = 8
include(joinpath(@__DIR__, "..", "ContinuumFF.jl"))
"""
Thermally averaged free–free Gaunt factors (van Hoof et al. 2014, `ContinuumFF.gaunt_ff`) in `kappa_nu`
instead of g = 1.  Default `false`: the shipped Z0 tables (S3) were built with g = 1; switching on
changes κ by the Gaunt factor (1.1–1.8 at 1e4–1e6 K) and requires `build_table` and the S3
re-validation.  Float64 only (the table lookup is not AD-aware).
"""
const GAUNT_FF = Ref(false)

# --- H− cross sections (John 1988, A&A 193, 189) -----------------------------------------
"H− bound-free cross section (cm²), λ in μm."
function σHm_bf(λμ)
    λ0 = 1.6419
    (λμ >= λ0 || λμ < 0.125) && return 0.0
    C = (152.519, 49.534, -118.858, 92.536, -34.194, 4.982)
    x = 1 / λμ - 1 / λ0
    f = sum(C[n] * x^((n - 1) / 2) for n in 1:6)
    return 1e-18 * λμ^3 * x^1.5 * f
end

const _HmA1 = ((0.0,2483.346,-3449.889,2200.040,-696.271,88.283),
               (0.0,285.827,-1158.382,2427.719,-1841.400,444.517),
               (0.0,-2054.291,8746.523,-13651.105,8624.970,-1863.864),
               (0.0,2827.776,-11485.632,16755.524,-10051.530,2095.288),
               (0.0,-1341.537,5303.609,-7510.494,4400.067,-901.788),
               (0.0,208.952,-812.939,1132.738,-655.020,132.985))
const _HmA2 = ((518.1021,473.2636,-482.2089,115.5291,0.0,0.0),
               (-734.8666,1443.4137,-737.1616,169.6374,0.0,0.0),
               (1021.1775,-1977.3395,1096.8827,-245.6490,0.0,0.0),
               (-479.0721,922.3575,-521.1341,114.2430,0.0,0.0),
               (93.1373,-178.9275,101.7963,-21.9972,0.0,0.0),
               (-6.4285,12.3600,-7.0920,1.5097,0.0,0.0))
"H− free-free absorption per neutral H atom per unit electron pressure (cm⁴/dyn), incl. stimulated emission."
function kHm_ff(λμ, T)
    λμ < 0.1823 && return 0.0
    A = λμ > 0.3645 ? _HmA1 : _HmA2
    θ = 5040.0 / T
    s = 0.0
    for n in 1:6
        s += θ^((n + 1) / 2) * (λμ^2 * A[1][n] + A[2][n] + A[3][n] / λμ + A[4][n] / λμ^2 +
                                A[5][n] / λμ^3 + A[6][n] / λμ^4)
    end
    return 1e-29 * s
end

"Rayleigh scattering on H I (Dalgarno), λ in Å; capped below 1300 Å (near Lyα)."
function σRay(λÅ)
    λ = max(λÅ, 1300.0)
    return 5.799e-13 / λ^4 + 1.422e-6 / λ^6 + 2.784 / λ^8
end

# --- monochromatic opacity (cm⁻¹) --------------------------------------------------------
"Level populations (cm⁻³) of H I and He II, hydrogenic, LTE within each ion."
function level_pops(T, s)
    kT = kB * T
    wH  = [2n^2 * exp(-χH  * (1 - 1 / n^2) / kT) for n in 1:NLEV_H]
    wHe = [2n^2 * exp(-χHe2 * (1 - 1 / n^2) / kT) for n in 1:NLEV_HE2]
    return s.nH0 .* wH ./ sum(wH), s.nHe1 .* wHe ./ sum(wHe)
end

"""
Absorption (incl. stimulated emission) and scattering coefficients (cm⁻¹) at frequency ν.
`pH`, `pHe` are level populations from `level_pops`.
"""
function kappa_nu(ν, T, s, pH, pHe)
    kT  = kB * T
    u   = h * ν / kT
    stim = -expm1(-u)
    λcm = c / ν; λμ = λcm * 1e4; λÅ = λcm * 1e8
    ab = 0.0
    for n in 1:NLEV_H
        ν >= χH / (h * n^2) && (ab += pH[n] * 2.815e29 / (n^5 * ν^3))
    end
    for n in 1:NLEV_HE2
        ν >= χHe2 / (h * n^2) && (ab += pHe[n] * 2.815e29 * 16 / (n^5 * ν^3))
    end
    ν0 = χHe1 / h
    ν >= ν0 && (ab += s.nHe0 * 7.4e-18 * (ν0 / ν)^1.7)
    ab *= stim
    if GAUNT_FF[]
        ab += 3.69e8 * s.ne * ((s.nHp + s.nHe1) * ContinuumFF.gaunt_ff(ν, T, 1) + 4s.nHe2 * ContinuumFF.gaunt_ff(ν, T, 2)) / sqrt(T) / ν^3 * stim
    else
        ab += 3.69e8 * s.ne * (s.nHp + s.nHe1 + 4s.nHe2) / sqrt(T) / ν^3 * stim
    end
    nHm = s.nH0 * s.ne * exp(χHm / kT) / (4λe3(T))
    ab += nHm * σHm_bf(λμ) * stim
    ab += s.nH0 * kHm_ff(λμ, T) * s.ne * kT
    sc = s.ne * σT + s.nH0 * σRay(λÅ)
    return ab, sc
end

const UGRID = exp.(range(log(1e-3), log(60.0), length = 4000))

"""
Frequency grid (u = hν/kT) with every opacity discontinuity inserted as a breakpoint, so that
midpoint quadrature never straddles an edge: H I edges n = 1..NLEV_H, He II edges n = 1..NLEV_HE2,
He I ground edge, H− bound-free threshold (1.6419 μm) and cutoff (0.125 μm), the H− free-free fit
switches (0.3645, 0.1823 μm) and the Rayleigh cap (1300 Å). Without this, edges sliding through
the fixed grid as T changes gave ~1 % quadrature noise in κ_P (κ_R was barely affected).
"""
function ugrid_with_edges(T)
    kT = kB * T
    edges = Float64[]
    for n in 1:NLEV_H;   push!(edges, χH / (n^2 * kT)); end
    for n in 1:NLEV_HE2; push!(edges, χHe2 / (n^2 * kT)); end
    push!(edges, χHe1 / kT)
    for λμ in (1.6419, 0.125, 0.3645, 0.1823, 0.13)
        push!(edges, h * c / (λμ * 1e-4) / kT)
    end
    filter!(u -> UGRID[1] < u < UGRID[end], edges)
    return sort!(unique!(vcat(UGRID, edges)))
end

"Rosseland and Planck means (cm²/g) at (ρ, T). Planck mean uses absorption only."
function kappa_means(ρ, T)
    s = saha(ρ, T)
    pH, pHe = level_pops(T, s)
    numR = 0.0; denR = 0.0; numP = 0.0; denP = 0.0
    ug = ugrid_with_edges(T)
    for i in 1:length(ug)-1
        u = 0.5(ug[i] + ug[i+1]); du = ug[i+1] - ug[i]
        ν = u * kB * T / h
        wR = u^4 * exp(u) / expm1(u)^2          # ∝ dB_ν/dT
        wP = u^3 / expm1(u)                     # ∝ B_ν
        ab, sc = kappa_nu(ν, T, s, pH, pHe)
        numR += wR / (ab + sc) * du; denR += wR * du
        numP += wP * ab * du;        denP += wP * du
    end
    return denR / numR / ρ, numP / denP / ρ
end
kappa_ross(ρ, T) = kappa_means(ρ, T)[1]
kappa_planck(ρ, T) = kappa_means(ρ, T)[2]

# --- table on (log T, log R), R = ρ / T6³ -------------------------------------------------
struct KapTable
    lT::Vector{Float64}
    lR::Vector{Float64}
    lk::Matrix{Float64}        # log10 κ_Rosseland
    lkP::Matrix{Float64}       # log10 κ_Planck (absorption)
end
KapTable(lT, lR, lk) = KapTable(lT, lR, lk, fill(NaN, size(lk)))   # v1 tables (no Planck)

function build_table(; lT = 3.3:0.01:8.0, lR = -12.0:0.05:1.0)
    lTv = collect(lT); lRv = collect(lR)
    lk = zeros(length(lTv), length(lRv)); lkP = similar(lk)
    Threads.@threads for i in eachindex(lTv)
        T = 10^lTv[i]
        for j in eachindex(lRv)
            ρ = 10^lRv[j] * (T / 1e6)^3
            kR, kP = kappa_means(ρ, T)
            lk[i, j] = log10(kR); lkP[i, j] = log10(kP)
        end
    end
    return KapTable(lTv, lRv, lk, lkP)
end

# Catmull–Rom cubic weights: interpolating and C¹ across cell edges.
@inline function _crw(t)
    t2 = t * t; t3 = t2 * t
    return ((-t3 + 2t2 - t) / 2, (3t3 - 5t2 + 2) / 2, (-3t3 + 4t2 + t) / 2, (t3 - t2) / 2)
end

"""
C¹ bicubic (Catmull–Rom) interpolation of a log10 table in (log T, log R); generic in the
number type, so ForwardDiff derivatives are continuous across cell edges. Outside the table
the arguments are clamped (derivative then zero) — `in_table` reports whether that happened.
"""
function _interp(tab::KapTable, M::Matrix{Float64}, ρ, T)
    lT = log10(T); lR = log10(ρ) - 3(lT - 6)
    nT, nR = length(tab.lT), length(tab.lR)
    x = (lT - tab.lT[1]) / (tab.lT[2] - tab.lT[1])
    y = (lR - tab.lR[1]) / (tab.lR[2] - tab.lR[1])
    x = clamp(x, 0, nT - 1 - 1e-9)
    if EXTRAP_LOWR[] && pval(y) < 0
        # below the lowest log R: linear extrapolation in log R with the slope between the first two columns
        # (fully ionized, thin gas: free-free/bound-free ∝ ρ, slope 1 for κ_P; electron scattering, slope 0 for κ_R)
        v0 = _cr2(M, x, zero(y), nT, nR); v1 = _cr2(M, x, one(y), nT, nR)
        return 10^(v0 + y * (v1 - v0))
    end
    y = clamp(y, 0, nR - 1 - 1e-9)
    return 10^_cr2(M, x, y, nT, nR)
end
"Catmull–Rom value of table M at fractional node coordinates (x, y) (both inside the table)."
function _cr2(M, x, y, nT, nR)
    i = floor(Int, pval(x)) + 1; j = floor(Int, pval(y)) + 1
    j = min(j, nR - 1); i = min(i, nT - 1)
    wx = _crw(x - (i - 1)); wy = _crw(y - (j - 1))
    v = zero(x * y)
    for a in 1:4, b in 1:4
        ii = clamp(i - 2 + a, 1, nT); jj = clamp(j - 2 + b, 1, nR)
        v += wx[a] * wy[b] * M[ii, jj]
    end
    v
end
"""
EXTRAP_LOWR[]: below the table's lowest log R (= log ρ − 3 log T₆) extrapolate log κ linearly in log R instead of
clamping. Clamping evaluates the opacity of gas up to 1e3 times denser: in a 1e7–1e8 K, 1e-17 g cm⁻³ accretion
corona it overstates the Planck mean, hence free-free cooling, by that factor. Default true.
"""
const EXTRAP_LOWR = Ref(true)
kap(tab::KapTable, ρ, T) = _interp(tab, tab.lk, ρ, T)
kapP(tab::KapTable, ρ, T) = _interp(tab, tab.lkP, ρ, T)

"Old bilinear lookup (v1), kept for regression comparison."
function kap_bilinear(tab::KapTable, ρ, T)
    lT = log10(T); lR = log10(ρ) - 3(lT - 6)
    x = clamp((lT - tab.lT[1]) / (tab.lT[2] - tab.lT[1]), 0, length(tab.lT) - 1.000001)
    y = clamp((lR - tab.lR[1]) / (tab.lR[2] - tab.lR[1]), 0, length(tab.lR) - 1.000001)
    i = floor(Int, pval(x)) + 1; j = floor(Int, pval(y)) + 1
    fx = x - (i - 1); fy = y - (j - 1)
    v = (1 - fx) * (1 - fy) * tab.lk[i, j] + fx * (1 - fy) * tab.lk[i+1, j] +
        (1 - fx) * fy * tab.lk[i, j+1] + fx * fy * tab.lk[i+1, j+1]
    return 10^v
end

in_table(tab::KapTable, ρ, T) = (lT = log10(pval(T)); lR = log10(pval(ρ)) - 3(lT - 6);
    tab.lT[1] <= lT <= tab.lT[end] && tab.lR[1] <= lR <= tab.lR[end])
