"""
    LymanWings

S7c-2: scattering cross section of ground-state atomic hydrogen from the infrared through the Lyman series (Lyα … Lyman limit) and above, and the
trapped-Lyα radiation pressure of McKee & Tan (2008). Self-contained (no dependencies).

**Cross section.** Kramers–Heisenberg formula for 1s → (np) → 1s, exact in the non-relativistic dipole limit, in Hartree atomic units
(ω in Ha, the Lyman limit at ω = 1/2):

    σ(ω) = σ_T ω⁴ |S(ω)|²,        S(ω) = Σ_{n≥2} f_n / (ω_n² − ω² − i ω γ_n)  +  ∫_{1/2}^∞ (df/dε) dε / (ε² − ω²)       (P.V. above the edge)

f_n = 2⁸ n⁵ (n−1)^{2n−4} / (3 (n+1)^{2n+4}) are the exact 1s–np oscillator strengths (f_2 = 0.41620), ω_n = (1 − 1/n²)/2, df/dε the exact bound–free
oscillator-strength density of the 1s level (Stobbe's formula: σ_bf = 6.304e-18 cm² at threshold), γ_n = A_n/(2.4189e-17 s) the total decay rate of np
(Lorentz damping; only the cores and the near wings of the lines depend on it). Limits: S(0) = 4.5 (the static polarizability 4.5 a₀³), the sum rule Σ f = 1 (σ → σ_T for ω → ∞,
a bound electron scatters like a free one), σ = σ_T ω⁴ S(0)² × (1 + O(ω²)) far below Lyα — Lee & Kim (2005, MNRAS 358, 1472) expansion σ/σ_T = (ω/ω_l)⁴ Σ c_p (ω/ω_l)^{2p},
c_0 = 1.2654 (checked against Gavrila 1967 there). On resonance the damped sum gives the textbook peak 3λ²/2π and the line integral π e² f_n/(m c) (checked in `s7c2_test.jl`).

**Thermal cores.** The scattering is coherent in the rest frame of the atom: the cross section of the gas is the convolution of the atom's σ(ω) with the Maxwellian line-of-sight velocity
distribution. For the first `NV` lines the resonant term |S_n|² (a Lorentzian times the slowly varying Rayleigh factor) is replaced by the same times the ratio of the Voigt to the Lorentz profile
(H(a,x)/(√π Δν_D)) / φ_L(ν) (same line integral), exact in the Doppler core and the wing alike, with Δν_D = ν_n b/c, b² = 2kT/m_H + v_turb². The interference with the rest of the sum and the
dispersive wing are kept as in the atom. Not included: hyperfine/fine structure, Lyα excitation of n = 2 atoms, the velocity field inside the line (Sobolev): all matter only within |x| ≲ 10.

**Frequency redistribution (how it is used).** In the transfer the scattering is *coherent in the observer's frame* (no frequency change; the direct beam is removed per extinction, the scattered
power feeds the diffuse field at the same frequency when it is on). In the wings (|x| ≫ 1, where the force on a thick infall comes from) this is the correct limit (Adams 1972; Neufeld 1990: the
redistribution in the observer frame by the thermal motion is small compared with the scale of the wing); in the Doppler cores it is wrong, where photons would diffuse in frequency: the
trapped-Lyα force of the core photons is added separately (`lya_pressure`, MT08 eq. 22). The momentum removed from a beam by an atom is the extinction momentum, wing or core.

**Tabulation.** `table(edges, Ts)`: bin-averaged σ̄(bin, T) with the Voigt cores (fine composite quadrature in the bins that contain or neighbour a line centre; 3-point Gauss–Legendre
elsewhere); `sigma_bin(tab, i, T)` interpolates in ln T. `line_edges` returns bin edges that resolve the lines (geometric steps in |ν − ν_n|).
"""
module LymanWings

# ------------------------------------------------------------------ constants (cgs and Hartree units)
const c = 2.99792458e10; const hP = 6.62607015e-27; const eV = 1.602176634e-12; const kB = 1.380649e-16
const σT = 6.6524587e-25; const a0sq = 2.8002852e-17                      # Thomson cross section, a₀² [cm²]
const νHa = 6.579683920502e15                                             # Hartree frequency E_h/h [Hz]
const νRy = 3.2880513e15                                                  # Rydberg frequency of hydrogen (reduced mass) [Hz]
const tau_au = 2.4188843265857e-17                                        # atomic unit of time [s]
const αfs = 1 / 137.035999084
const mH = 1.6735575e-24
const πe2mc = 0.02654008                                                  # π e²/(m c) [cm² Hz]
const E0 = νRy / νHa                                                       # Lyman limit in Ha (reduced mass): 0.49973
const NLINES = 300                                                        # explicit lines of the sum (the rest is added as a smooth tail)
const NV = 12                                                             # lines with a thermal (Voigt) core

# ------------------------------------------------------------------ atomic data
"Oscillator strength 1s → np (exact hydrogen)."
f_np(n) = exp(8log(2) + 5log(n) + (2n - 4) * log(n - 1) - log(3) - (2n + 4) * log(n + 1))
"Transition frequency of Ly-n [Hz] (reduced-mass hydrogen)."
ν_line(n) = νRy * (1 - 1 / n^2)
"Line frequency in Ha (the electron-mass unit; ν/ν_Ha), i.e. including the reduced-mass shift."
ω_line(n) = ν_line(n) / νHa
"A(np → 1s) [s⁻¹] from f: (2/3) ω² f α³ in atomic units."
A_np1s(n) = (2 / 3) * ω_line(n)^2 * f_np(n) * αfs^3 / tau_au
# branching ratio of np to 1s (the rest goes to ns, nd): 2p 1, 3p 0.882, 4p 0.84, ≥ 5p 0.81 (NIST A_ki, to a few %; they only set the damping width)
branch(n) = n == 2 ? 1.0 : n == 3 ? 0.882 : n == 4 ? 0.84 : 0.81
"Total decay rate of np [s⁻¹]."
A_tot(n) = A_np1s(n) / branch(n)

const FN = [n >= 2 ? f_np(n) : 0.0 for n in 1:NLINES]
const WN = [n >= 2 ? ω_line(n) : 0.0 for n in 1:NLINES]
const GN = [n >= 2 ? A_tot(n) * tau_au : 0.0 for n in 1:NLINES]            # γ_n in Ha

# ------------------------------------------------------------------ continuum: df/dε and its polarizability integral
"Stobbe's 1s photoionisation cross section [cm²] at photon energy ε [Ha] (threshold 1/2 Ha, reduced-mass shift ignored)."
function σ_bf(ε)
    e0 = E0; ε <= e0 && return 0.0
    η = sqrt(ε / e0 - 1)
    6.3043e-18 * (e0 / ε)^4 * exp(4 - 4 * atan(η) / η) / (-expm1(-2π / η))
end
"df/dε [1/Ha] of the 1s continuum (σ = 2π² α (df/dε) a₀² Ha)."
dfde(ε) = σ_bf(ε) / a0sq / (2π^2 * αfs)

# Gauss–Legendre nodes (n points) on [-1, 1] by Newton iteration
function gl(n)
    x = zeros(n); w = zeros(n)
    for i in 1:n
        z = cos(π * (i - 0.25) / (n + 0.5))
        for _ in 1:100
            p1 = 1.0; p2 = 0.0
            for j in 1:n; p3 = p2; p2 = p1; p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j; end
            pp = n * (z * p1 - p2) / (z^2 - 1); dz = p1 / pp; z -= dz
            abs(dz) < 1e-15 && break
        end
        p1 = 1.0; p2 = 0.0
        for j in 1:n; p3 = p2; p2 = p1; p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j; end
        pp = n * (z * p1 - p2) / (z^2 - 1)
        x[i] = -z; w[i] = 2 / ((1 - z^2) * pp^2)
    end
    (x, w)
end
const GLX, GLW = gl(64)
const G5X, G5W = gl(5)

"∫_{1/2}^∞ (df/dε)/(ε² − ω²) dε for ω < 1/2, and its principal value for ω ≥ 1/2 [1/Ha²]."
function S_cont(ω)
    e0 = E0
    if ω < e0 * 0.999999
        s = 0.0
        for i in eachindex(GLX)
            y = 0.5 * (GLX[i] + 1); ε = e0 / y
            s += 0.5 * GLW[i] * dfde(ε) / (ε^2 - ω^2) * e0 / y^2
        end
        return s
    end
    h(ε) = dfde(ε) / (ε + ω)
    d = ω - e0; s = 0.0
    if d > 0        # symmetric part around ω: ∫_0^d [h(ω+u) − h(ω−u)]/u du
        for i in eachindex(GLX)
            u = 0.5 * d * (GLX[i] + 1)
            s += 0.5 * d * GLW[i] * (h(ω + u) - h(ω - u)) / u
        end
    end
    # remainder from ω + d to ∞ (regular): ε = (ω+d)/y
    L = ω + d
    for i in eachindex(GLX)
        y = 0.5 * (GLX[i] + 1); ε = L / y
        s += 0.5 * GLW[i] * h(ε) / (ε - ω) * L / y^2
    end
    s
end

"Bound–bound part of S without line `skip` (0 = none): Σ_n f_n/(ω_n² − ω² − iωγ_n), plus the analytic tail beyond NLINES [complex, Ha⁻²]."
function S_lines(ω; skip = 0)
    s = 0.0 + 0.0im
    for n in 2:NLINES
        n == skip && continue
        s += FN[n] / (WN[n]^2 - ω^2 - im * ω * GN[n])
    end
    # tail: f_n ≈ f_∞/n³ (f_∞ = 1.5625.. from df/dε at threshold), ω_n → ω_∞ ; Σ_{n>N} ≈ ∫ dn f_∞/n³ → f_∞/(2 (N + 1/2)²)
    ftail = dfde(E0) / (2 * (NLINES + 0.5)^2)
    s + ftail / (E0^2 - ω^2)
end
S_total(ω; skip = 0) = S_lines(ω; skip) + S_cont(ω)
"Resonant term of line n [complex]."
S_res(ω, n) = FN[n] / (WN[n]^2 - ω^2 - im * ω * GN[n])

"Coherent Rayleigh/Kramers–Heisenberg cross section [cm²] at frequency ν [Hz] (damped, no thermal motion)."
σ_coh(ν) = (ω = ν / νHa; σT * ω^4 * abs2(S_total(ω)))

# ------------------------------------------------------------------ Voigt function (Humlíček 1982, W4; relative error < 1e-4)
function voigt_H(a, x)
    t = complex(a, -abs(x)); s = abs(x) + a
    if s >= 15
        w = t * 0.5641896 / (0.5 + t * t)
    elseif s >= 5.5
        u = t * t; w = t * (1.410474 + u * 0.5641896) / (0.75 + u * (3 + u))
    elseif a >= 0.195 * abs(x) - 0.176
        w = (16.4955 + t * (20.20933 + t * (11.96482 + t * (3.778987 + t * 0.5642236)))) /
            (16.4955 + t * (38.82363 + t * (39.27121 + t * (21.69274 + t * (6.699398 + t)))))
    else
        u = t * t
        w = exp(u) - t * (36183.31 - u * (3321.9905 - u * (1540.787 - u * (219.0313 - u * (35.76683 - u * (1.320522 - u * 0.56419)))))) /
                     (32066.6 - u * (24322.84 - u * (9022.228 - u * (2186.181 - u * (364.2191 - u * (61.57037 - u * (1.841439 - u)))))))
    end
    real(w)
end

# ------------------------------------------------------------------ thermal line profile
"Doppler width Δν_D [Hz] of line n for gas temperature T [K] and extra turbulent speed vt [cm/s]."
Δν_D(n, T; vt = 0.0) = ν_line(n) / c * sqrt(2kB * T / mH + vt^2)
"Resonant atom cross section |σ_T ω⁴ S_n|² of line n at ν [Hz]."
σ_res(ν, n) = (ω = ν / νHa; σT * ω^4 * abs2(S_res(ω, n)))
"Lorentz profile (per Hz) of line n, same width and centre as the damped term: HWHM A/(4π)."
φ_L(ν, n) = (Γ = A_tot(n); (Γ / (4π^2)) / ((ν - ν_line(n))^2 + (Γ / (4π))^2))
"Voigt profile (per Hz) of line n."
function φ_V(ν, n, T; vt = 0.0)
    D = Δν_D(n, T; vt); a = A_tot(n) / (4π * D); x = (ν - ν_line(n)) / D
    voigt_H(a, x) / (sqrt(π) * D)
end
"Correction σ_res,n(ν) (φ_V/φ_L − 1) of the thermal motion to the cross section of the gas near line n [cm²]."
function σ_thermal_corr(ν, n, T; vt = 0.0)
    pl = φ_L(ν, n)
    σ_res(ν, n) * (φ_V(ν, n, T; vt) / pl - 1)
end
"σ(ν, T) of the gas [cm²]: damped KH cross section of the atom + the thermal-core corrections of the first NV lines (those within 300 Doppler widths)."
function σ_gas(ν, T; vt = 0.0, nv = NV)
    s = σ_coh(ν)
    for n in 2:nv
        D = Δν_D(n, T; vt)
        abs(ν - ν_line(n)) < 300D || continue
        s += σ_thermal_corr(ν, n, T; vt)
    end
    s
end

# ------------------------------------------------------------------ bin averages
"Composite Gauss–Legendre integral of f over [a, b] with break points geometric around `ν0` from `d0` (the part of the interval on each side of ν0)."
function geo_integral(f, a, b, ν0, d0)
    pts = Float64[a, b]
    ν0 > a && ν0 < b && push!(pts, ν0)
    d = d0
    while d < 2 * max(abs(a - ν0), abs(b - ν0))
        for p in (ν0 - d, ν0 + d); p > a && p < b && push!(pts, p); end
        d *= 2
    end
    sort!(unique!(pts))
    s = 0.0
    for i in 1:length(pts)-1
        lo, hi = pts[i], pts[i+1]; hw = 0.5 * (hi - lo); mid = 0.5 * (hi + lo)
        for k in eachindex(G5X); s += G5W[k] * hw * f(mid + hw * G5X[k]); end
    end
    s
end
const GL3X = (-sqrt(3 / 5), 0.0, sqrt(3 / 5)); const GL3W = (5 / 9, 8 / 9, 5 / 9)
function gl3(f, a, b)
    hw = 0.5 * (b - a); m = 0.5 * (a + b); s = 0.0
    for k in 1:3; s += GL3W[k] * f(m + hw * GL3X[k]); end
    0.5 * s      # mean over the bin (weights sum to 2)
end

"""
    bin_average(a, b, T; nclose = 3, vt = 0, nv = NV) -> σ̄ [cm²]

Mean of σ_gas over the frequency bin [a, b] (Hz). A bin that contains a line centre of the first `nv` lines, or lies within `nclose` bin widths of one, is integrated in two parts: the
resonant term with its Voigt profile (composite quadrature, geometric about the line centre) and the remainder (interference, other lines, continuum) by 3-point Gauss–Legendre; other bins by
3 points, plus the Voigt-core corrections of the lines within 300 Doppler widths.
"""
function bin_average(a, b, T; nclose = 3, vt = 0.0, nv = NV)
    W = b - a; m = 0.5 * (a + b)
    nearest = 0; dmin = Inf
    for n in 2:nv
        d = abs(m - ν_line(n)); d < dmin && (dmin = d; nearest = n)
    end
    if nearest > 0 && dmin <= nclose * W + 0.5 * W
        n = nearest; ν0 = ν_line(n); D = Δν_D(n, T; vt)
        # resonant part with the Voigt profile: σ_res,n(ν) φ_V(ν)/φ_L(ν) = σ_res(ν) × ratio
        fres = ν -> σ_res(ν, n) * (φ_V(ν, n, T; vt) / φ_L(ν, n))
        res = geo_integral(fres, a, b, ν0, 0.25 * D) / W
        # remainder: |S_n + S_r|² − |S_n|² = 2 Re(S_n* S_r) + |S_r|²
        frem = ν -> (ω = ν / νHa; Sr = S_total(ω; skip = n); Sn = S_res(ω, n); σT * ω^4 * (2 * real(conj(Sn) * Sr) + abs2(Sr)))
        rem = gl3(frem, a, b)
        # other lines' thermal corrections within the window (rare): add by 3 points
        extra = 0.0
        for k in 2:nv
            k == n && continue
            abs(m - ν_line(k)) < 300 * Δν_D(k, T; vt) + W && (extra += gl3(ν -> σ_thermal_corr(ν, k, T; vt), a, b))
        end
        return res + rem + extra
    end
    gl3(ν -> σ_gas(ν, T; vt, nv), a, b)
end

"""
    line_edges(n1, n2; dmin = 1e-5, per_decade = 6, dmax = 0.1) -> Vector{Float64}

Frequency edges [Hz] that resolve the lines Ly-n, n = n1 … n2: a pair of edges at ν_n(1 ± δ) for δ geometric from `dmin` to min(`dmax`, half the distance to the neighbours)
with `per_decade` steps per decade. The central bin is [ν_n(1−dmin), ν_n(1+dmin)] (the Doppler core of the cold gas, Δν_D/ν = 1.1e-5 at 1000 K).
"""
function line_edges(n1, n2; dmin = 1e-5, per_decade = 6, dmax = 0.1)
    e = Float64[]
    for n in n1:n2
        ν0 = ν_line(n)
        dl = n > 2 ? (ν0 - ν_line(n - 1)) / ν0 / 2 : dmax; dr = n < NLINES ? (ν_line(n + 1) - ν0) / ν0 / 2 : dmax
        for (d, sgn) in ((min(dl, dmax), -1), (min(dr, dmax), +1))
            δ = dmin
            while δ < d
                push!(e, ν0 * (1 + sgn * δ)); δ *= 10^(1 / per_decade)
            end
        end
        push!(e, ν0 * (1 - dmin), ν0 * (1 + dmin))
    end
    sort!(unique!(e))
end

# ------------------------------------------------------------------ the table
"""
σ̄ per bin and temperature (`sig[bin, iT]`, `lnT`) and, for the *structured* bins (within 3 % of a line centre of the first NV lines), a sub-sampling of σ(ν, T) (`sub_ν`, `sub_w` = fraction of
the bin width, `sub_sig[bin][q, iT]`) for the trapped-radiation integrals. Other bins have a single sample (the bin average).
"""
struct Table
    edges::Vector{Float64}; lnT::Vector{Float64}; sig::Matrix{Float64}
    structured::BitVector
    sub_ν::Vector{Vector{Float64}}; sub_w::Vector{Vector{Float64}}; sub_sig::Vector{Matrix{Float64}}
end

"Nodes and weights (fractions of the bin width) of the composite quadrature used in `bin_average` for the resonant term: geometric break points about ν0."
function geo_nodes(a, b, ν0, d0)
    pts = Float64[a, b]
    ν0 > a && ν0 < b && push!(pts, ν0)
    d = d0
    while d < 2 * max(abs(a - ν0), abs(b - ν0))
        for p in (ν0 - d, ν0 + d); p > a && p < b && push!(pts, p); end
        d *= 2
    end
    sort!(unique!(pts))
    νq = Float64[]; wq = Float64[]
    for i in 1:length(pts)-1
        lo, hi = pts[i], pts[i+1]; hw = 0.5 * (hi - lo); mid = 0.5 * (hi + lo)
        for k in eachindex(G5X); push!(νq, mid + hw * G5X[k]); push!(wq, G5W[k] * hw / (b - a)); end
    end
    (νq, wq)
end

"""
    table(edges; Ts = 10 .^ range(1, 5.5, length = 37), vt = 0.0, near = 0.03) -> Table

Bin-averaged σ(ν, T) of the gas for bins with edges `edges` [Hz]. Bins away from the lines (more than 300 Doppler widths of the hottest tabulated T) are T-independent and are computed once.
"""
function table(edges; Ts = 10 .^ range(1, 5.5, length = 37), vt = 0.0, near = 0.03, threads = true)
    nb = length(edges) - 1; nT = length(Ts); sig = zeros(nb, nT)
    Tmax = maximum(Ts)
    nearb = falses(nb); structured = falses(nb); nearest = zeros(Int, nb)
    for i in 1:nb
        m = 0.5 * (edges[i] + edges[i+1]); W = edges[i+1] - edges[i]
        dmin = Inf
        for n in 2:NV
            d = abs(m - ν_line(n))
            d < 300 * Δν_D(n, Tmax; vt) + 4W && (nearb[i] = true)
            d < near * ν_line(n) && (structured[i] = true)
            d < dmin && (dmin = d; nearest[i] = n)
        end
    end
    sub_ν = [Float64[] for _ in 1:nb]; sub_w = [Float64[] for _ in 1:nb]; sub_sig = [zeros(0, 0) for _ in 1:nb]
    body(i) = begin
        a, b = edges[i], edges[i+1]
        if nearb[i]
            for (jt, T) in enumerate(Ts); sig[i, jt] = bin_average(a, b, T; vt); end
        else
            s = bin_average(a, b, Ts[1]; vt); sig[i, :] .= s
        end
        if structured[i]
            n = nearest[i]
            νq, wq = geo_nodes(a, b, ν_line(n), 0.25 * Δν_D(n, Ts[1]; vt))
            sub_ν[i] = νq; sub_w[i] = wq
            S = zeros(length(νq), nT)
            for (jt, T) in enumerate(Ts), q in eachindex(νq); S[q, jt] = σ_gas(νq[q], T; vt); end
            sub_sig[i] = S
        end
    end
    if threads && Threads.nthreads() > 1
        Threads.@threads for i in 1:nb; body(i); end
    else
        for i in 1:nb; body(i); end
    end
    Table(collect(Float64, edges), log.(Ts), sig, structured, sub_ν, sub_w, sub_sig)
end
function _tw(tab::Table, T)
    lt = clamp(log(T), tab.lnT[1], tab.lnT[end]); nT = length(tab.lnT)
    nT == 1 && return (1, 0.0)
    j = clamp(searchsortedlast(tab.lnT, lt), 1, nT - 1)
    (j, (lt - tab.lnT[j]) / (tab.lnT[j+1] - tab.lnT[j]))
end
_interp(s1, s2, w) = (s1 > 0 && s2 > 0) ? exp((1 - w) * log(s1) + w * log(s2)) : (1 - w) * s1 + w * s2
"σ̄ of bin i at temperature T (log–log interpolation in T)."
function sigma_bin(tab::Table, i, T)
    nT = length(tab.lnT)
    if nT == 1 || tab.sig[i, 1] == tab.sig[i, nT]; return tab.sig[i, 1]; end
    j, w = _tw(tab, T)
    _interp(tab.sig[i, j], tab.sig[i, j+1], w)
end
"σ at the sub-samples of the structured bin i at temperature T (`sub_ν`, `sub_w`, σ)."
function sigma_sub(tab::Table, i, T)
    S = tab.sub_sig[i]; j, w = _tw(tab, T)
    σ = [_interp(S[q, j], S[q, min(j + 1, size(S, 2))], w) for q in 1:size(S, 1)]
    (tab.sub_ν[i], tab.sub_w[i], σ)
end

# ------------------------------------------------------------------ trapped radiation pressure (McKee & Tan 2008, eq. 22 and Appendix B), generalised to the Lyman series
"""
    line_plateau(n, Neff, T; vt = 0, eps = 0) -> (Jplat, xe, aτ)

MT08 eqs. B7–B8, B13–B14: the isotropic mean intensity J/F (in units of the injected flux) of photons injected in the Doppler core and the diffusion regime of line Ly-n of a slab of effective H I
column Neff, `0.411 (aτ̄)^{1/3} / (min[1, (aτ̄/450)^{1/3}] + Γ)`, Γ = 1.04 (ε τ̄)^{1/3} the two-photon destruction (Lyα only; ε the probability per scattering). `xe = (aτ̄)^{1/3}` is
the half-width of the escape region in Doppler units. Returns J = 0 for aτ̄ < 1 (the line is not trapped).
"""
function line_plateau(n, Neff, T; vt = 0.0, eps = 0.0)
    D = Δν_D(n, T; vt); a = A_tot(n) / (4π * D); s̄ = πe2mc * FN[n] / D; τ̄ = s̄ * Neff; aτ = a * τ̄
    aτ < 1 && return (0.0, 0.0, aτ)
    xe = aτ^(1 / 3); sat = min(1.0, (aτ / 450)^(1 / 3))
    Γ = n == 2 ? 1.04 * (eps * τ̄)^(1 / 3) : 0.0
    (0.411 * xe / (sat + Γ), xe, aτ)
end
"Damping parameter a and σ̄ = π e² f/(m c Δν_D) of Lyα (kept for the tests)."
function lya_params(T; vt = 0.0)
    D = Δν_D(2, T; vt)
    (A_tot(2) / (4π * D), πe2mc * FN[2] / D)
end
"Rayleigh factor f(ν) of MT08 eq. A4 (= 1 at line centre; 0.5 at 0.8 ν₀)."
f_ray(ν) = (ν0 = ν_line(2); 4ν^4 / (ν0^2 * (ν0 + ν)^2))

"""
    jratio(ν, σ, Neff, T; vt = 0, eps = 0) -> J/F

Isotropic intensity J (frequency-integrated) per unit flux F of photons injected at frequency ν into the neutral gas of effective column `Neff` with scattering cross section `σ` [cm²]
at ν: the harmonic blend of the coherent-scattering limit J_coh = (3/4π) τ_ν F, τ_ν = Neff σ (MT08 B10: the force of the beam τ F/c, no saturation because the scattered photons stay in the
net flux) and the diffusion plateau of the nearest line (`line_plateau`), `1/J = 1/J_coh + 1/J_plat` — MT08's B11 for Lyα (their 5.41 x̂²/f term is 1/J_coh with the Lorentz wing), extended to the
KH cross section (all lines and the Rayleigh continuum) and to the lines Ly-β … with the plateau of the line nearest in Doppler widths.
"""
function jratio(ν, σ, Neff, T; vt = 0.0, eps = 0.0, nv = NV)
    Jc = 3 / (4π) * Neff * σ
    Jc <= 0 && return 0.0
    nbest = 2; xbest = Inf
    for n in 2:nv
        x = abs(ν - ν_line(n)) / Δν_D(n, T; vt); x < xbest && (xbest = x; nbest = n)
    end
    Jp = line_plateau(nbest, Neff, T; vt, eps)[1]
    Jp == 0 && return Jc
    Jc * Jp / (Jc + Jp)
end

"""
    lya_pressure(νs, Fν, Neff, T; Bν = nothing, vt = 0.0, eps = 0.0) -> P_iso [dyn cm⁻²]

MT08 eq. 22 / B22 for the Lyα line alone (the original formula, B11): P_iso = (4π/3c) ∫ dν_i min[B_ν(T*), (J/F)(ν_i) F_ν(ν_i)] with `νs` ascending (Hz), `Fν` the incident flux per Hz at the
position and (J/F) = 0.411 (aτ̄)^{1/3}/(min[1, (aτ̄/450)^{1/3}] + 5.41 x̂²/f(ν) + Γ(x̂)).
"""
function lya_pressure(νs, Fν, Neff, T; Bν = nothing, vt = 0.0, eps = 0.0, xhmax = 40.0)
    a, s̄ = lya_params(T; vt); τ̄ = s̄ * Neff; aτ = a * τ̄
    aτ < 1 && return 0.0
    D = Δν_D(2, T; vt); xe = aτ^(1 / 3); sat = min(1.0, (aτ / 450)^(1 / 3))
    s = 0.0
    for i in 1:length(νs)-1
        νm = 0.5 * (νs[i] + νs[i+1]); dν = νs[i+1] - νs[i]
        x̂ = (νm - ν_line(2)) / D / xe
        abs(x̂) > xhmax && continue
        Γ = abs(x̂) < 1 ? 1.04 * (eps * τ̄)^(1 / 3) : 0.0
        r = 0.411 * xe / (sat + 5.41 * x̂^2 / f_ray(νm) + Γ)
        j = r * 0.5 * (Fν[i] + Fν[i+1])
        Bν === nothing || (j = min(j, 0.5 * (Bν[i] + Bν[i+1])))
        s += j * dν
    end
    4π / (3c) * s
end

end # module
