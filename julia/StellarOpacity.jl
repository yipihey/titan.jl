"""
    StellarOpacity

Population III project, S3: Rosseland and Planck mean opacities κ_R(ρ, T, X), κ_P(ρ, T, X) of metal-free gas
(Y = 1 − X, Z_CNO = 0) from the molecular regime (log T = 1.5) to log T = 9.5, AD-ready (ForwardDiff Duals pass through
every function), C¹ in log T, log R (R = ρ/T₆³) and X.

Runtime = one precomputed table `stellar/kaptable_stellar.bin.gz` (built by `stellaropacity_build.jl`, which documents the
sources and the blends): log10 κ_R and log10 κ_P on a uniform grid in (log T, log R) for the ten X nodes of the OPAL
GN93hz Z = 0 tables, Catmull–Rom bicubic in (log T, log R) and a shape-preserving cubic Hermite (PCHIP) in X. Region
flags (`coverage`) tell which physics a point rests on. Outside the grid the arguments are clamped (derivative 0);
`in_table` reports that.

Rosseland mean, by log T:
- 1.5–3.8   own continuum means (ContinuumOpacity: CIA, H⁻, Rayleigh, Thomson) in chemical equilibrium (H₂ ⇌ 2H,
            H ⇌ H⁺ + e⁻, He ⇌ He⁺ + e⁻, Li electrons, H⁻ by Saha); no atomic or H₂ lines;
- 3.8–4.1   log-κ blend (smoothstep) of the own means and OPAL;
- 4.1–8.7   OPAL (Iglesias & Rogers 1996, GN93hz, Z = 0), where it exists; beyond its edges (log R > 1, R < −8,
            the hot high-R corner, T > 8.7) OPAL extended by the analytic high-T model, anchored to OPAL's last node;
- all T ≥ 1e4: harmonic combination with the degenerate-electron conduction (Cassisi et al. 2007).
Planck mean: own means below, continuum Planck mean of H/He (Saha, Gaunt factors 1, no lines) above (see S3 notes).

The analytic high-T model (`kappa_high`) lives here: it is also used standalone.
"""
module StellarOpacity

using ForwardDiff
using Serialization

export kappa_R, kappa_R_rad, kappa_cond, kappa_P, kappas, log_kappa_R, log_kappa_P, kappa_high, kappa_es, in_table, coverage, XNODES

const XNODES = [0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0]
const LT0 = 1.5;  const DLT = 0.025; const NT = 321      # log T 1.5 … 9.5
const LR0 = -15.0; const DLR = 0.25;  const NR = 81       # log R −15 … +5
const TABLEFILE = joinpath(@__DIR__, "stellar", "kaptable_stellar.bin.gz")

val(x::Real) = Float64(x)
val(x::ForwardDiff.Dual) = val(ForwardDiff.value(x))
smoothstep(x) = (x = clamp(x, 0, 1); x * x * (3 - 2x))

# ------------------------------------------------------------------ analytic high-T model
"""
Thermal (Compton) reduction of electron scattering, f(T) = 1/(1 + (T/4.5e8)^0.95): the Buchler & Yueh (1976)/Paczyński
(1983) form; exponent 0.95 (instead of 0.86) refitted to OPAL's R = 10⁻⁸ column (pure scattering, 1e6–3e8 K, X = 0–1,
agrees to 1 %).
"""
thermal_factor(T) = 1 / (1 + (T / 4.5e8)^0.95)
"Electron degeneracy reduction of the scattering opacity, Paczyński (1983)/Iben (1975): 1/(1 + 2.7e11 ρ/T²)."
degen_factor(ρ, T) = 1 / (1 + 2.7e11 * ρ / T^2)
"Electron-scattering opacity [cm² g⁻¹]: 0.2(1+X) × thermal × degeneracy factors."
kappa_es(ρ, T, X) = 0.2 * (1 + X) * thermal_factor(T) * degen_factor(ρ, T)
# Empirical free–free parameters (fit to OPAL GN93hz Z = 0, X = 0, 0.35, 0.7, 1; log T 6.8–8.5, log R ≥ −3; rms 0.063 dex):
# κ_ff = g(T) 3.68e22 (1+X) ρ T^−3.5 / (1 + b 2.7e11 ρ/T²)^p,  g(T) = g0 + g1 (log T − 7.5) + g2 (log T − 7.5)²
const FF = (g0 = 1.904472940665797, g1 = 0.28529763042300693, g2 = -0.17981939175863548, b = 22.23520786829821, p = 0.36585663768350685)
"Free–free (Kramers form with an OPAL-fitted effective Gaunt factor and degeneracy suppression) [cm² g⁻¹]."
function kappa_ff(ρ, T, X)
    lT = clamp(log10(T), 5.0, 9.5)
    g = FF.g0 + FF.g1 * (lT - 7.5) + FF.g2 * (lT - 7.5)^2
    g * 3.68e22 * (1 + X) * ρ * T^-3.5 / (1 + FF.b * 2.7e11 * ρ / T^2)^FF.p
end
"""
High-T radiative Rosseland opacity of fully ionised H/He: electron scattering (Compton + degeneracy) + free–free.
Empirical (fitted to OPAL where it exists, rms 0.06 dex over log T 6.8–8.5, log R ≥ −3, worst 0.4 dex at the degenerate
edge); within the table it is only used as the T-shape beyond OPAL's edges, anchored to OPAL there.
"""
kappa_high(ρ, T, X) = kappa_es(ρ, T, X) + kappa_ff(ρ, T, X)

# ------------------------------------------------------------------ table
struct Table
    lkR::Array{Float64,3}; lkP::Array{Float64,3}     # [iT, iR, iX]: radiative log10 κ_R (no conduction), log10 κ_P
    flag::Array{UInt8,3}
    lKH::Matrix{Float64}; lKHe::Matrix{Float64}      # Cassisi et al. 2007 log10 K_cond [erg/cm/s/K], [log ρ −6…11.5, log T 3…10] step 0.25
end
const TAB = Ref{Union{Table,Nothing}}(nothing)

function load_table(path = TABLEFILE)
    io = IOBuffer(read(`gzip -dc $path`))
    nT, nR, nX = read(io, Int64), read(io, Int64), read(io, Int64)
    (nT == NT && nR == NR && nX == length(XNODES)) || error("StellarOpacity: table shape $nT×$nR×$nX does not match the module")
    rd(n) = collect(reinterpret(Float64, read(io, 8n)))
    lkR = reshape(rd(nT * nR * nX), nT, nR, nX); lkP = reshape(rd(nT * nR * nX), nT, nR, nX)
    flag = reshape(read(io, nT * nR * nX), nT, nR, nX)
    lKH = reshape(rd(71 * 29), 71, 29); lKHe = reshape(rd(71 * 29), 71, 29)
    Table(lkR, lkP, flag, lKH, lKHe)
end
table() = (TAB[] === nothing && (TAB[] = load_table()); TAB[]::Table)

@inline cr(p0, p1, p2, p3, t) = p1 + 0.5t * (p2 - p0 + t * (2p0 - 5p1 + 4p2 - p3 + t * (3(p1 - p2) + p3 - p0)))
"Catmull–Rom bicubic value of A[:, :, k] at fractional grid coordinates (x, y) (1-based, edge-clamped), generic in x, y."
@inline function bicubic(A, k, x, y)
    n, m = size(A, 1), size(A, 2)
    xv = val(x); yv = val(y)
    i = clamp(floor(Int, xv), 1, n - 1); j = clamp(floor(Int, yv), 1, m - 1)
    tx = x - i; ty = y - j
    i0 = max(i - 1, 1); i3 = min(i + 2, n)
    @inbounds col(jj) = (jc = clamp(jj, 1, m); cr(A[i0, jc, k], A[i, jc, k], A[i+1, jc, k], A[i3, jc, k], tx))
    cr(col(j - 1), col(j), col(j + 1), col(j + 2), ty)
end

# PCHIP slope at node m (1-based) of the X nodes with values ys (function m -> value)
@inline function pchip_slope(m, y)
    N = length(XNODES)
    if m == 1
        h1 = XNODES[2] - XNODES[1]; h2 = XNODES[3] - XNODES[2]
        d1 = (y(2) - y(1)) / h1; d2 = (y(3) - y(2)) / h2
        s = ((2h1 + h2) * d1 - h1 * d2) / (h1 + h2)
        val(s) * val(d1) <= 0 ? zero(s) * 1 : (val(d1) * val(d2) <= 0 && abs(val(s)) > 3abs(val(d1)) ? 3d1 : s)
    elseif m == N
        h1 = XNODES[N] - XNODES[N-1]; h2 = XNODES[N-1] - XNODES[N-2]
        d1 = (y(N) - y(N - 1)) / h1; d2 = (y(N - 1) - y(N - 2)) / h2
        s = ((2h1 + h2) * d1 - h1 * d2) / (h1 + h2)
        val(s) * val(d1) <= 0 ? zero(s) * 1 : (val(d1) * val(d2) <= 0 && abs(val(s)) > 3abs(val(d1)) ? 3d1 : s)
    else
        hl = XNODES[m] - XNODES[m-1]; hr = XNODES[m+1] - XNODES[m]
        dl = (y(m) - y(m - 1)) / hl; dr = (y(m + 1) - y(m)) / hr
        val(dl) * val(dr) <= 0 && return zero(dl) * 1
        w1 = 2hr + hl; w2 = hr + 2hl
        (w1 + w2) / (w1 / dl + w2 / dr)
    end
end
"Shape-preserving (PCHIP, C¹) interpolation in X of the node values `y(m)`."
@inline function xinterp(y, X)
    N = length(XNODES)
    Xv = clamp(val(X), 0.0, 1.0)
    k = clamp(searchsortedlast(XNODES, Xv), 1, N - 1)
    h = XNODES[k+1] - XNODES[k]
    Xc = X < 0 ? zero(X) : (X > 1 ? one(X) : X)
    t = (Xc - XNODES[k]) / h
    y0 = y(k); y1 = y(k + 1); d0 = pchip_slope(k, y); d1 = pchip_slope(k + 1, y)
    t2 = t * t; t3 = t2 * t
    (2t3 - 3t2 + 1) * y0 + (t3 - 2t2 + t) * h * d0 + (-2t3 + 3t2) * y1 + (t3 - t2) * h * d1
end

@inline coords(lT, lR) = (clamp((lT - LT0) / DLT, 0, NT - 1 - 1e-9) + 1, clamp((lR - LR0) / DLR, 0, NR - 1 - 1e-9) + 1)
"log10 κ_R at (log T, log R, X) (arguments clamped to the table)."
function log_kappa_R(lT, lR, X)
    A = table().lkR; x, y = coords(lT, lR)
    xinterp(m -> bicubic(A, m, x, y), X)
end
"log10 κ_P at (log T, log R, X)."
function log_kappa_P(lT, lR, X)
    A = table().lkP; x, y = coords(lT, lR)
    xinterp(m -> bicubic(A, m, x, y), X)
end
@inline lTlR(ρ, T) = (lT = log10(T); (lT, log10(ρ) - 3 * (lT - 6)))
"log10 κ_R of the radiative opacity alone (OPAL/own means/high-T model; no conduction)."
log_kappa_rad(lT, lR, X) = log_kappa_R(lT, lR, X)

# --- conduction (Cassisi et al. 2007 tables, Z = 1 and 2, as the mixture 1/K = X/K_H + Y/K_He)
const CRHO0 = -6.0; const CT0 = 3.0; const CD = 0.25
@inline function cr_cond(M, x, y)
    n, m = size(M)
    j = clamp(floor(Int, val(x)), 1, m - 1); i = clamp(floor(Int, val(y)), 1, n - 1)
    tx = x - j; ty = y - i
    @inbounds row(ii) = (ic = clamp(ii, 1, n); cr(M[ic, clamp(j - 1, 1, m)], M[ic, j], M[ic, j+1], M[ic, clamp(j + 2, 1, m)], tx))
    cr(row(i - 1), row(i), row(i + 1), row(i + 2), ty)
end
"log10 K of one table at (log ρ, log T): Catmull–Rom, T clamped to 10³–10¹⁰, log ρ extended linearly beyond the table."
function log_K(M, lρ, lT)
    n, m = size(M)
    x = clamp((lT - CT0) / CD + 1, 1, m - 1e-9)
    y = (lρ - CRHO0) / CD + 1
    top = n - 1e-9
    if y > top
        return cr_cond(M, x, top) + (cr_cond(M, x, top) - cr_cond(M, x, top - 1)) * (y - top)
    end
    cr_cond(M, x, max(y, 1))
end
"Conduction switch-on (smoothsteps in log ρ: −6…−5, log T: 4…5): the tables assume ionised matter."
wcond(lρ, lT) = smoothstep(lρ + 6) * smoothstep(lT - 4)
"Electron-conduction opacity κ_cond = 16σT³/(3ρK) [cm² g⁻¹] (unswitched: valid for ionised matter)."
function kappa_cond(ρ, T, X)
    t = table(); lT = log10(T); lρ = log10(ρ)
    iK = X * 10^(-log_K(t.lKH, lρ, lT)) + (1 - X) * 10^(-log_K(t.lKHe, lρ, lT))
    16 * 5.670374e-5 * T^3 * iK / (3ρ)
end
"""
Rosseland mean including conduction: 1/κ = 1/κ_rad + w_c/κ_cond (the heat flux is carried by radiation and electrons in
parallel; κ_cond is defined by the same form as κ_rad, F = −16σT³/(3κρ) ∇T). Equals `kappa_R_rad` where w_c = 0.
"""
function kappa_R(ρ, T, X)
    lT, lR = lTlR(ρ, T)
    κr = 10^log_kappa_R(lT, lR, X)
    w = wcond(log10(ρ), lT)
    val(w) == 0 ? κr : 1 / (1 / κr + w / kappa_cond(ρ, T, X))
end
"Radiative Rosseland mean (without conduction) [cm² g⁻¹]."
kappa_R_rad(ρ, T, X) = (t = lTlR(ρ, T); 10^log_kappa_R(t[1], t[2], X))
"Planck mean opacity (absorption only) [cm² g⁻¹]."
kappa_P(ρ, T, X) = (t = lTlR(ρ, T); 10^log_kappa_P(t[1], t[2], X))
"(κ_R with conduction, κ_P)."
kappas(ρ, T, X) = (kappa_R(ρ, T, X), kappa_P(ρ, T, X))
"Whether (ρ, T) is inside the table (else the arguments were clamped and the derivative is zero in that direction)."
function in_table(ρ, T)
    lT = log10(val(T)); lR = log10(val(ρ)) - 3 * (lT - 6)
    LT0 <= lT <= LT0 + DLT * (NT - 1) && LR0 <= lR <= LR0 + DLR * (NR - 1)
end
"""
Region flag at (log T, log R, X), nearest node: 1 = OPAL data, 2 = beyond OPAL's edges (extended/anchored OPAL, or high-T
model), 3 = blend zone (own means ↔ OPAL, or near OPAL's cold-corner holes), 4 = own continuum means only, 0 = none.
"""
function coverage(lT, lR, X)
    t = table(); x, y = coords(lT, lR); k = argmin(abs.(collect(XNODES) .- X))
    t.flag[round(Int, x), round(Int, y), k]
end

end # module
