# S8b: the 21 cm hyperfine line of atomic hydrogen as an (opt-in) diagnostic of the Pop III collapse / feedback runs.
#   include("HI21.jl"); using .HI21      (module HI21; no dependencies; nothing in the hydro or chemistry reads it)
#
# (1) Spin temperature per zone,   T_s^-1 = (T_γ^-1 + x_c T_k^-1 + x_α T_c^-1)/(1 + x_c + x_α)       [Field 1958; Furlanetto, Oh & Briggs 2006 (FOB06) eq. 20]
#       x_c = (T*/A10 T_γ) Σ_i n_i κ10^i(T_k)   i = H (Zygelman 2005 as tabulated in FOB06 Table 3), e (Furlanetto & Furlanetto 2007a Table 1), p (FF07b Table 1),
#             H₂ (no data: κ = `KH2` × κ_HH, default 0), He, He⁺ (not calculated / closed shell: 0)
#       x_α = 4 P_α T*/(27 A10 T_γ) S_α     P_α = 4π ∫ σ_ν J_ν/hν dν  = Lyα scattering rate per H atom (Field 1958; Wouthuysen 1952)
#       S_α, T_c: Hirata (2006) fits (eqs 40–42) when the input Lyα is the *continuum* intensity near the line (spectral distortion by recoil, fine structure, spin
#       diffusion, incomplete relaxation of the colour temperature); S_α = 1 and the supplied T_c when the input is the actual (e.g. Monte Carlo) line intensity.
#     `ts_rates` is the same thing without the T* ≪ T_γ linearisation (explicit 0↔1 rates incl. stimulated terms); `relax` is the time-dependent level population
#     (rate n₁/n₀ with the same rates), with the equilibration time 1/(R_up + R_down).
# (2) Line transfer: 21 cm opacity α_ν = (3 λ³ A10/8π) n₀ (1 − e^{−T*/T_s}) φ_v (thermal Gaussian profile in velocity), Sobolev τ per zone, rays through a spherical cloud
#     (impact parameters, velocity-resolved, source function T_s, background T_γ): T_b(v) = T_γ e^{−τ} + ∫ T_s e^{−(τ−τ')} dτ'; observed δT_b = (T_b − T_γ)/(1+z).
module HI21

using Printf

# ------------------------------------------------------------------ constants (cgs)
const kB = 1.380649e-16; const hP = 6.62607015e-27; const c = 2.99792458e10; const mH = 1.6735575e-24
const NU10 = 1.420405751768e9                     # 21 cm line frequency [Hz]
const TSTAR = hP * NU10 / kB                      # 0.0682 K
const A10 = 2.85e-15                              # s⁻¹
const LAM21 = c / NU10                            # 21.106 cm
const LAM_A = 1215.668e-8; const NU_A = c / LAM_A # Lyα
const A_A = 6.2649e8                              # Lyα Einstein A [s⁻¹]
const GAM_A = A_A / (4π)                          # Hirata's γ [Hz]
const SIGINT_A = 3 * LAM_A^2 * A_A / (8π)         # ∫σ_ν dν of Lyα = π e² f/(m_e c) = 1.105e-2 cm² Hz (g_u/g_l = 3)
const G = 6.67430e-8

# ------------------------------------------------------------------ κ10 tables [cm³ s⁻¹]
# H–H: Zygelman (2005) for T_K ≤ 300 K, Sigurdson (priv. comm., FOB06 method) above; FOB06 Table 3.
const T_HH = [1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000.0]
const K_HH = 1e-13 .* [1.38, 1.43, 2.71, 6.60, 14.7, 28.8, 91.0, 178, 273, 367, 538, 686, 814, 925, 1020, 1110, 1190, 1750, 2090, 2560, 2910, 3310, 4270, 4970, 6030, 6870, 7870.0]
# H–e⁻ (all partial waves below the n = 2 threshold; Furlanetto & Furlanetto 2007a Table 1 = FOB06 Table 4)
const T_EH = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000.0]
const K_EH = 1e-9 .* [0.239, 0.337, 0.530, 0.746, 1.05, 1.63, 2.26, 3.11, 4.59, 5.92, 7.15, 7.71, 8.17, 8.32, 8.37, 8.29, 8.11]
# H–p (Furlanetto & Furlanetto 2007b Table 1; Krstić et al. 2004 cross sections)
const T_PH = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000.0]
const K_PH = 1e-9 .* [0.4028, 0.4517, 0.4301, 0.3699, 0.3172, 0.3047, 0.3379, 0.4043, 0.5471, 0.7051, 0.9167, 1.070, 1.301, 1.480, 1.695, 1.975, 2.201]

"Log–log interpolation, clamped to the table ends (above the last row the rate is held constant: excitation to n ≥ 2 is not in the tables; below, the rate is flat or slowly varying)."
function loglog(Tt, Kt, T)
    T <= Tt[1] && return Kt[1]; T >= Tt[end] && return Kt[end]
    j = searchsortedlast(Tt, T); w = log(T / Tt[j]) / log(Tt[j+1] / Tt[j])
    exp((1 - w) * log(Kt[j]) + w * log(Kt[j+1]))
end
κ_HH(T) = loglog(T_HH, K_HH, T)
κ_eH(T) = loglog(T_EH, K_EH, T)
κ_pH(T) = loglog(T_PH, K_PH, T)

# ------------------------------------------------------------------ background temperature at 21 cm
"""
    Tgamma(Tcmb, z; Ar = 0, β = -2.6) -> T_γ at the source (z) [K]
CMB plus an optional power-law radio excess, T_rad(ν_obs) = Ar (ν_obs/78.1 MHz)^β [K] at the observed frequency ν_obs = ν21/(1+z) (Feng & Holder 2018 form);
its source-frame value is (1+z) T_rad. Default: CMB only.
"""
Tgamma(Tcmb, z = 0.0; Ar = 0.0, β = -2.6) = Ar == 0 ? Tcmb : Tcmb + (1 + z) * Ar * (NU10 / (1 + z) / 78.1e6)^β

# ------------------------------------------------------------------ coupling coefficients
"""
    xcoll(Tk, Tg; nH, ne, np, nH2 = 0, KH2 = 0.0) -> (xc, parts)
Collisional coupling x_c = (T*/A10 T_γ) Σ n_i κ10^i. `KH2`: H–H₂ coefficient as a multiple of κ_HH (no calculation exists; ortho/para H₂ and H spin exchange needs
an open-shell partner, so the default is 0).
"""
function xcoll(Tk, Tg; nH, ne = 0.0, np = 0.0, nH2 = 0.0, KH2 = 0.0)
    f = TSTAR / (A10 * Tg)
    xh = f * nH * κ_HH(Tk); xe = f * ne * κ_eH(Tk); xp = f * np * κ_pH(Tk); x2 = f * nH2 * KH2 * κ_HH(Tk)
    (xh + xe + xp + x2, (; H = xh, e = xe, p = xp, H2 = x2))
end

"x_α for a Lyα scattering rate per atom `P` [s⁻¹] and scattering factor `S` (1 if J is the actual line intensity)."
xalpha(P, Tg, S = 1.0) = 4 * P * TSTAR / (27 * A10 * Tg) * S
"Lyα scattering rate per H atom from the mean intensity of photons near the line, `Jν` [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹]: P = 4π (Jν/hν_α) ∫σ dν  (flat J across the line)."
P_from_J(Jν) = 4π * Jν / (hP * NU_A) * SIGINT_A
"Lyα Gunn–Peterson-like depth τ_GP = 3 A_α n_HI λ_α³/(8π g) for a velocity gradient `g` = |dv/dr| (Hubble: H) [s⁻¹]."
tauGP(nHI, g) = 3 * A_A * nHI * LAM_A^3 / (8π * max(g, 1e-30))

"Hirata (2006) eq. 40: S̃_α(T_k, T_s, τ_GP) (valid T_k, T_s ≥ 2 K, 1e5 ≤ τ_GP ≤ 1e7; arguments outside are clamped)."
function S_tilde(Tk, Ts, τ)
    Tk = max(Tk, 2.0); Ts = max(Ts, 2.0); τ = clamp(τ, 1e5, 1e7)
    ξ = (1e-7 * τ)^(1 / 3) * Tk^(-2 / 3)
    (1 - 0.0631789 / Tk + 0.115995 / Tk^2 - 0.401403 / (Ts * Tk) + 0.336463 / (Ts * Tk^2)) / (1 + 2.98394ξ + 1.53583ξ^2 + 3.85289ξ^3)
end
"Hirata (2006) eq. 42: effective colour temperature."
Tc_eff(Tk, Ts) = (Tk = max(Tk, 2.0); 1 / (1 / Tk + 0.405535 / Tk * (1 / max(Ts, 2.0) - 1 / Tk)))

"Spin temperature, linearised form (FOB06 eq. 20)."
ts_lin(Tg, Tk, xc, xa, Tc) = (1 + xc + xa) / (1 / Tg + xc / Tk + xa / Tc)

"""
    spin_temperature(Tg, Tk, xc; P = 0, Tc = NaN, kind = :continuum, τGP = 2e6, iter = 12) -> (Ts, xα, Tc, S)
`kind = :continuum`: P from the continuum intensity at Lyα, Hirata's S̃_α and T_c^eff (implicit in T_s: fixed-point iteration). `kind = :actual`: P is the true scattering rate (Monte Carlo
line transfer), S = 1, T_c as given (NaN → T_k). With P = 0 the result is the Field (1958) collisional + CMB formula.
"""
function spin_temperature(Tg, Tk, xc; P = 0.0, Tc = NaN, kind = :continuum, τGP = 2e6, iter = 12)
    P <= 0 && return (ts_lin(Tg, Tk, xc, 0.0, Tk), 0.0, Tk, 1.0)
    if kind === :actual
        Tcv = isnan(Tc) ? Tk : Tc; xa = xalpha(P, Tg)
        return (ts_lin(Tg, Tk, xc, xa, Tcv), xa, Tcv, 1.0)
    end
    Ts = Tg; xa = 0.0; Tcv = Tk; S = 1.0
    for _ in 1:iter
        Tcv = Tc_eff(Tk, Ts); S = S_tilde(Tk, Ts, τGP); xa = xalpha(P, Tg, S)
        Tn = ts_lin(Tg, Tk, xc, xa, Tcv)
        abs(Tn - Ts) < 1e-9 * Ts && (Ts = Tn; break); Ts = Tn
    end
    (Ts, xa, Tcv, S)
end

# ------------------------------------------------------------------ explicit rates; time-dependent level population
"""
    rates(Tg, Tk, xc, xa, Tc) -> (up, down)
Per-atom rates [s⁻¹] 0→1 (per n₀) and 1→0 (per n₁): up = 3A10 n̄ + 3(C10 e^{-T*/Tk} + P10 e^{-T*/Tc}), down = A10(1 + n̄) + C10 + P10, with n̄ = 1/(e^{T*/Tγ} − 1),
C10 = x_c A10 T_γ/T*, P10 = x_α A10 T_γ/T* (= 4P_α/27). The upper-level fraction y = n₁/n_HI relaxes to up/(up + down) at the rate up + down.
"""
function rates(Tg, Tk, xc, xa, Tc)
    nbar = 1 / expm1(TSTAR / Tg)
    C = xc * A10 * Tg / TSTAR; Pf = xa * A10 * Tg / TSTAR
    (3A10 * nbar + 3C * exp(-TSTAR / Tk) + 3Pf * exp(-TSTAR / Tc), A10 * (1 + nbar) + C + Pf)
end
ts_from_y(y) = y <= 0 ? 0.0 : -TSTAR / log(y / (3 * (1 - y)))
y_from_ts(Ts) = (r = 3 * exp(-TSTAR / Ts); r / (1 + r))
"Equilibrium spin temperature from the explicit rates (no T* ≪ T_γ linearisation)."
function ts_rates(Tg, Tk, xc, xa, Tc)
    up, dn = rates(Tg, Tk, xc, xa, Tc); ts_from_y(up / (up + dn))
end
"Time-dependent level population: y₀ (upper fraction) after `dt` seconds with constant coefficients → (y, T_s, t_eq)."
function relax(y0, dt, Tg, Tk, xc, xa, Tc)
    up, dn = rates(Tg, Tk, xc, xa, Tc); λ = up + dn; ye = up / λ
    y = ye + (y0 - ye) * exp(-λ * dt)
    (y, ts_from_y(y), 1 / λ)
end

# ------------------------------------------------------------------ line opacity
"Line-centre-normalised opacity prefactor: α_v = κ0(n_HI, T_s) φ_v, φ_v = exp(-(v/b)²)/(√π b) per cm/s; κ0 = (3 λ³ A10/8π) n₀ (1 − e^{-T*/T_s})."
κ0(nHI, Ts) = 3 * LAM21^3 * A10 / (8π) * nHI / (1 + 3exp(-TSTAR / Ts)) * (-expm1(-TSTAR / Ts))
bth(Tk) = sqrt(2 * kB * Tk / mH)
"Sobolev optical depth for a velocity gradient `g` [s⁻¹] (τ = κ0/g); thermal depth over a path `L` at line centre (τ = κ0 L/(√π b))."
tau_sobolev(nHI, Ts, g) = κ0(nHI, Ts) / max(g, 1e-30)
tau_thermal(nHI, Ts, Tk, L) = κ0(nHI, Ts) * L / (sqrt(π) * bth(Tk))

# ------------------------------------------------------------------ cosmology (for angular sizes)
"Angular-diameter distance [cm] in flat ΛCDM (default: ABN02: Ω_m = 0.3, h = 0.7)."
function DA(z; Ωm = 0.3, h = 0.7)
    H0 = h * 100 * 1e5 / 3.0857e24; n = 2000; s = 0.0
    for i in 1:n; zz = (i - 0.5) * z / n; s += 1 / sqrt(Ωm * (1 + zz)^3 + 1 - Ωm); end
    c / H0 * s * z / n / (1 + z)
end
Hz(z; Ωm = 0.3, h = 0.7) = h * 100 * 1e5 / 3.0857e24 * sqrt(Ωm * (1 + z)^3 + 1 - Ωm)

# ------------------------------------------------------------------ Lyα input (pluggable)
"""
    LyaIn(P, Tc, kind)
Per-zone Lyα scattering rate per H atom `P` [s⁻¹], colour temperature `Tc` [K] (NaN → T_k) and the meaning of `P`: `:continuum` (intensity of the photons near the line without the
line's own scattering: Hirata S̃_α and T_c^eff are applied) or `:actual` (the true line intensity, e.g. from Monte Carlo transfer: S = 1, T_c as given).
"""
struct LyaIn
    P::Vector{Float64}; Tc::Vector{Float64}; kind::Symbol
end
lya_none(n) = LyaIn(zeros(n), fill(NaN, n), :continuum)
"From RadTransfer's zone mean intensity at the Lyα bin, `Jν` [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹] per zone (continuum, T_c = T_k relaxed through Hirata's fit)."
lya_from_Jnu(Jν) = LyaIn(P_from_J.(Jν), fill(NaN, length(Jν)), :continuum)
"From the Monte Carlo module (S8a): J_α and T_c per zone (the actual line intensity)."
lya_actual(Jν, Tc) = LyaIn(P_from_J.(Jν), collect(float.(Tc)), :actual)
"Hook: a function (state NamedTuple with r, T, ρ, nHI, …) -> LyaIn; set it to LyaMC's wrapper when available. Default: no Lyα."
const LYA_PROVIDER = Ref{Any}(nothing)
"""
    use_lyamc!(f)
Switch the Lyα input to the Monte Carlo module (S8a, `LyaMC.jl`): `f(s)` takes a snapshot NamedTuple (r, u, ρ, T, y, m) and returns `(; J, Tc)` per zone — `J` the mean intensity of the
*actual* line radiation as a number intensity × hν_α (erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹ at the line, the same units as RadTransfer's J) and `Tc` its
colour temperature [K]. With it S_α = 1 and no Hirata correction is applied (the Monte Carlo already contains recoil and the line profile). The hook is not wired to LyaMC yet: the API of that module was not
available when this was written.
"""
use_lyamc!(f) = (LYA_PROVIDER[] = s -> (r = f(s); lya_actual(r.J, r.Tc)); nothing)

# ------------------------------------------------------------------ per-zone analysis
const yHe = 1 / 12
"Number densities per zone from a snapshot (ρ, y = (H⁺, H⁻, H₂, H₂⁺) per H nucleus as `ycomp` of popiii_m2c.jl)."
function densities(ρ, y)
    μH = mH * (1 + 4yHe); nuc = ρ / μH
    xH = 1 - y[1] - y[2] - 2y[3] - 2y[4]
    (; nHI = max(xH, 0.0) * nuc, np = y[1] * nuc, ne = max(y[1] + y[4] - y[2], 0.0) * nuc, nH2 = y[3] * nuc, nuc)
end

"""
    zones(s; Tcmb, z, lya = nothing, KH2 = 0.0, neq = nothing, dt = 0.0, radio = (0.0, -2.6)) -> NamedTuple of per-zone vectors
`s`: snapshot NamedTuple (r [n+1], u [n+1], ρ, T, y, m [n+1]) of popiii_m2c.jl. `lya`: `LyaIn` (default none). Returns T_s (equilibrium, linearised), T_s from the explicit rates, x_c parts, x_α,
Sobolev / thermal / combined τ21 (velocity gradient: max(|du/dr|, |u/r|, H(z))), the equilibration time t_eq and t_ff = √(3π/32Gρ). `neq`: upper-level fractions carried from
the previous snapshot (time-dependent mode; advanced by `dt`).
"""
function zones(s; Tcmb, z, lya = nothing, KH2 = 0.0, neq = nothing, dt = 0.0, radio = (0.0, -2.6))
    n = length(s.T); Tg = Tgamma(Tcmb, z; Ar = radio[1], β = radio[2])
    L = lya === nothing ? lya_none(n) : lya
    o = (; nHI = zeros(n), ne = zeros(n), np = zeros(n), nH2 = zeros(n), Tk = collect(float.(s.T)), Tg = fill(Tg, n), xc = zeros(n), xcH = zeros(n), xce = zeros(n), xcp = zeros(n),
         xa = zeros(n), Tc = zeros(n), S = ones(n), Ts = zeros(n), Ts_rates = zeros(n), Ts_neq = fill(NaN, n), y = zeros(n), teq = zeros(n), tff = zeros(n), grad = zeros(n),
         tau_sob = zeros(n), tau_th = zeros(n), tau = zeros(n), Ts_use = zeros(n), tauGP = zeros(n), P = copy(L.P), dr = zeros(n))
    Hh = Hz(z)
    for j in 1:n
        d = densities(s.ρ[j], s.y[j]); o.nHI[j] = d.nHI; o.ne[j] = d.ne; o.np[j] = d.np; o.nH2[j] = d.nH2
        Tk = s.T[j]; xc, parts = xcoll(Tk, Tg; nH = d.nHI, ne = d.ne, np = d.np, nH2 = d.nH2, KH2 = KH2)
        o.xc[j] = xc; o.xcH[j] = parts.H; o.xce[j] = parts.e; o.xcp[j] = parts.p
        rc = 0.5 * (s.r[j] + s.r[j+1]); dr = s.r[j+1] - s.r[j]; o.dr[j] = dr
        g = max(abs((s.u[j+1] - s.u[j]) / dr), abs(0.5 * (s.u[j] + s.u[j+1]) / max(rc, 1e-30)), Hh); o.grad[j] = g
        τg = tauGP(d.nHI, g); o.tauGP[j] = τg
        Ts, xa, Tc, S = spin_temperature(Tg, Tk, xc; P = L.P[j], Tc = L.Tc[j], kind = L.kind, τGP = τg)
        o.Ts[j] = Ts; o.xa[j] = xa; o.Tc[j] = Tc; o.S[j] = S
        up, dn = rates(Tg, Tk, xc, xa, Tc); o.Ts_rates[j] = ts_from_y(up / (up + dn)); o.teq[j] = 1 / (up + dn); o.y[j] = up / (up + dn)
        if neq !== nothing
            yn, Tn, _ = relax(neq[j], dt, Tg, Tk, xc, xa, Tc); o.y[j] = yn; o.Ts_neq[j] = Tn
        end
        Tsu = isnan(o.Ts_neq[j]) ? Ts : o.Ts_neq[j]; o.Ts_use[j] = Tsu
        o.tau_sob[j] = tau_sobolev(d.nHI, Tsu, g); o.tau_th[j] = tau_thermal(d.nHI, Tsu, Tk, dr)
        # a proper profile: the Sobolev depth holds where the resonance region (b/g) is thinner than the zone, the thermal depth over the zone otherwise
        o.tau[j] = 1 / (1 / max(o.tau_sob[j], 1e-300) + 1 / max(o.tau_th[j], 1e-300))
        o.tff[j] = sqrt(3π / (32 * G * s.ρ[j]))
    end
    o
end

# ------------------------------------------------------------------ rays through the spherical cloud
"""
    spectrum(r, u, Tk, Ts, nHI, p, vgrid; Tbg, hubble = 0.0) -> Tb(v)
Brightness temperature (source frame) of the ray with impact parameter `p` [cm] through the shells `r` (n+1 interfaces), radial velocity `u` at the interfaces [cm/s], zone T_k, T_s, n_HI
(constant per zone), as a function of the velocity v = −u·μ (positive: receding from the observer) on `vgrid` (uniform, [cm/s]). The ray is integrated from the far side to the observer
(+s) in sub-steps of ≤ b/2 in velocity: T_b ← T_b e^{−Δτ} + T_s (1 − e^{−Δτ}) starting from `Tbg` (T_γ). `hubble` adds the Hubble flow H r to u.
"""
function spectrum(r, u, Tk, Ts, nHI, p, vgrid; Tbg, hubble = 0.0)
    n = length(Tk); nv = length(vgrid); v0 = vgrid[1]; dv = vgrid[2] - vgrid[1]
    Tb = fill(float(Tbg), nv)
    p >= r[end] && return Tb
    kfirst = max(searchsortedlast(r, p), 1)           # shell containing p (or the first shell)
    uu(k, rr) = (w = (rr - r[k]) / (r[k+1] - r[k]); (u[k] + w * (u[k+1] - u[k])) + hubble * rr)
    function segment!(k, sa, sb)                      # integrate s from sa to sb (sa < sb: toward the observer; the caller orders far → near)
        (nHI[k] <= 0 || sb == sa) && return
        b = bth(Tk[k]); k0 = κ0(nHI[k], Ts[k]); T_s = Ts[k]
        # sub-step count from the velocity change along the chord
        tv = 0.0; vprev = 0.0; m = 8
        for i in 0:m
            s = sa + (sb - sa) * i / m; rr = sqrt(p^2 + s^2); v = -uu(k, clamp(rr, r[k], r[k+1])) * s / max(rr, 1e-300)
            i > 0 && (tv += abs(v - vprev)); vprev = v
        end
        ns = clamp(ceil(Int, 2tv / b) + 2, 4, 4000); ds = (sb - sa) / ns
        for i in 1:ns
            s = sa + (i - 0.5) * ds; rr = clamp(sqrt(p^2 + s^2), r[k], r[k+1]); vc = -uu(k, rr) * s / max(rr, 1e-300)
            i1 = max(1, floor(Int, (vc - 5b - v0) / dv) + 1); i2 = min(nv, ceil(Int, (vc + 5b - v0) / dv) + 1)
            f = k0 * abs(ds) / (sqrt(π) * b)
            for iv in i1:i2
                x = (vgrid[iv] - vc) / b; dτ = f * exp(-x * x)
                dτ > 1e-30 && (Tb[iv] = Tb[iv] * exp(-dτ) + T_s * (-expm1(-dτ)))
            end
        end
    end
    # far side (s < 0): shells from the outside in
    for k in n:-1:kfirst
        r[k+1] > p || continue
        sb = sqrt(max(r[k+1]^2 - p^2, 0.0)); sa = r[k] > p ? sqrt(r[k]^2 - p^2) : 0.0
        if r[k] > p; segment!(k, -sb, -sa); else; segment!(k, -sb, sb); end     # the shell containing p is crossed in one go
    end
    # near side
    for k in (r[kfirst] > p ? kfirst : kfirst + 1):n
        r[k] > p || continue
        segment!(k, sqrt(r[k]^2 - p^2), sqrt(r[k+1]^2 - p^2))
    end
    Tb
end

"""
    cloud_spectra(s, o; Tcmb… ) -> NamedTuple
Velocity-resolved ray calculation of the whole cloud from the per-zone result `o` of `zones`: impact parameters p = 0, r_1 … r_n, area-weighted mean spectrum over the disc of radius R_out,
the central-ray spectrum, δT_b(v) in the source frame (T_b − T_γ) and observed (÷(1+z)). `rmax_frac`: only the rays with p ≤ rmax_frac R are used (the outer envelope is a very large area of
nearly CMB-temperature gas; the disc averages of a cloud and of its envelope differ).
"""
function cloud_spectra(s, o; z, hubble = 0.0, dvfrac = 4.0, nvmax = 4000, pmax = Inf)
    n = length(o.Tk); Tg = o.Tg[1]
    u = collect(float.(s.u)); r = collect(float.(s.r))
    bmin = minimum(bth(o.Tk[j]) for j in 1:n if o.nHI[j] > 0; init = bth(100.0))
    bmax = maximum(bth(o.Tk[j]) for j in 1:n if o.nHI[j] > 0; init = bth(100.0))
    vlim = maximum(abs, u) + (hubble > 0 ? hubble * r[end] : 0.0) + 6bmax
    dv = max(bmin / dvfrac, 2vlim / nvmax)
    vg = collect(-vlim:dv:vlim)
    ps = unique(sort([0.0; filter(<(min(pmax, r[end])), r[1:n])]))
    # annulus areas
    Rout = min(pmax, r[end]); pm = [0.0; 0.5 .* (ps[1:end-1] .+ ps[2:end]); Rout]
    w = [π * (pm[i+1]^2 - pm[i]^2) for i in eachindex(ps)]
    dT = zeros(length(ps), length(vg))
    for (i, p) in enumerate(ps)
        Tb = spectrum(r, u, o.Tk, o.Ts_use, o.nHI, p, vg; Tbg = Tg, hubble = hubble)
        dT[i, :] = Tb .- Tg
    end
    mean = vec(sum(dT .* w, dims = 1)) ./ sum(w)
    (; v = vg, p = ps, w, dT, mean, Rout, Tg, z)
end

"Equivalent width ∫ δT dv [K km/s] of a spectrum on a uniform grid."
eqwidth(v, dT) = sum(dT) * (v[2] - v[1]) / 1e5

# ------------------------------------------------------------------ one snapshot: zones + rays + integrated quantities
const Msun = 1.98847e33; const AU = 1.495978707e13; const YR = 3.15576e7
const PCORE = 1e4 * AU                    # aperture radius of the 'core' (collapsing cloud) mean: 0.048 pc

"""
    analyze(s; Tcmb, z, lya = nothing, neq = nothing, dt = 0.0, rays = true, hubble = 0.0, KH2 = 0.0, radio = (0.0, -2.6), pmax = Inf) -> (o, sp, sm)
`o`: per-zone result of `zones`; `sp`: ray result of `cloud_spectra` (or nothing); `sm`: scalar summary (NamedTuple; see `SUMMARY_FIELDS`). Observed brightness temperatures are
(T_b − T_γ)/(1 + z) in mK; the flux density uses the cloud's projected disc πR²/D_A² (R = outer radius of the snapshot).
"""
function analyze(s; Tcmb, z, lya = nothing, neq = nothing, dt = 0.0, rays = true, hubble = 0.0, KH2 = 0.0, radio = (0.0, -2.6), pmax = Inf)
    o = zones(s; Tcmb = Tcmb, z = z, lya = lya, neq = neq, dt = dt, KH2 = KH2, radio = radio)
    n = length(o.Tk); Tg = o.Tg[1]
    dm = [s.m[j+1] - s.m[j] for j in 1:n]
    xHm = [densities(s.ρ[j], s.y[j]).nHI * mH / s.ρ[j] for j in 1:n]       # neutral-H mass fraction of the zone (H atoms only)
    mHI = max.(dm .* xHm, 0.0); MHI = sum(mHI)
    w = MHI > 0 ? mHI ./ MHI : zeros(n)
    Vz = [4π / 3 * (s.r[j+1]^3 - s.r[j]^3) for j in 1:n]
    L21 = sum(0.75 * A10 * hP * NU10 * o.nHI[j] * Vz[j] for j in 1:n)                  # optically thin spontaneous emission rate [erg/s] (T_s ≫ T*)
    jmax = argmax(o.tau); abs_ = o.Ts .< Tg; Mabs = sum(mHI[abs_]); Mem = MHI - Mabs
    sp = nothing; dmin = dmax = vmin = vmax = cmin = cmax = ew = Sp = Sint = Kmn = Kmx = Sk = NaN
    D = DA(z); R = s.r[end]; θ = R / D * 206264.806247
    if rays
        sp = cloud_spectra(s, o; z = z, hubble = hubble, pmax = pmax)
        obs = sp.mean ./ (1 + z) .* 1e3; cen = sp.dT[1, :] ./ (1 + z) .* 1e3        # mK
        imn = argmin(obs); imx = argmax(obs)
        dmin = obs[imn]; vmin = sp.v[imn] / 1e5; dmax = obs[imx]; vmax = sp.v[imx] / 1e5
        cmin = minimum(cen); cmax = maximum(cen)
        ew = eqwidth(sp.v, obs)                                                          # mK km/s (disc mean)
        # flux density of the disc-mean spectrum at the observer: S_ν = (2 k ν_obs² δT_obs/c²) π R²/D_A²  [erg s⁻¹ cm⁻² Hz⁻¹] → nJy (1 Jy = 1e-23)
        νo = NU10 / (1 + z); Ω = π * sp.Rout^2 / D^2
        Sfun(Tm) = 2kB * νo^2 * Tm / c^2 * Ω / 1e-23 * 1e9
        ic = findall(p -> p < PCORE, sp.p); wc = sp.w[ic]; kd = vec(sum(sp.dT[ic, :] .* wc, dims = 1)) ./ sum(wc) ./ (1 + z) .* 1e3     # core aperture p < PCORE (mean over the disc)
        Kmn = minimum(kd); Kmx = maximum(kd); Sk = Sfun((abs(Kmn) > abs(Kmx) ? Kmn : Kmx) * 1e-3 * sum(wc) / (π * sp.Rout^2))
        Sp = Sfun((abs(obs[imn]) > abs(obs[imx]) ? obs[imn] : obs[imx]) * 1e-3); Sint = Sfun(1.0) * ew * 1e-3 * 1e5 * νo / c   # peak (signed, larger |δT|) [nJy]; ∫S dν [nJy Hz] = S(1 K) · (∫δT dv [K cm/s]) ν_o/c
    end
    teqratio = o.teq ./ o.tff
    sm = (; ncen = log10(s.ρ[1] / mH), Tcen = s.T[1], Rout_AU = R / AU, MHI = MHI / Msun, Tk_mw = sum(w .* o.Tk), Ts_mw = sum(w .* o.Ts_use), Ts_min = minimum(o.Ts_use[o.nHI .> 0]; init = Inf),
           Ts_max = maximum(o.Ts_use[o.nHI .> 0]; init = 0.0), tau_max = o.tau[jmax], r_taumax_AU = 0.5 * (s.r[jmax] + s.r[jmax+1]) / AU, Mabs = Mabs / Msun, Mem = Mem / Msun, L21 = L21,
           teq_max_over_tff = maximum(teqratio[o.nHI .> 0]; init = 0.0), f_noneq = sum(w[teqratio .> 1]; init = 0.0), dTb_disc_min = dmin, v_dmin = vmin, dTb_disc_max = dmax, v_dmax = vmax,
           dTb_cen_min = cmin, dTb_cen_max = cmax, dTb_core_min = Kmn, dTb_core_max = Kmx, S_core_nJy = Sk, EW = ew, S_peak_nJy = Sp, S_int_nJyHz = Sint, theta_arcsec = θ, DA_Mpc = D / 3.0857e24, Tgamma = Tg)
    (o, sp, sm)
end

const SUMMARY_FIELDS = (:ncen, :Tcen, :Rout_AU, :MHI, :Tk_mw, :Ts_mw, :Ts_min, :Ts_max, :tau_max, :r_taumax_AU, :Mabs, :Mem, :L21, :teq_max_over_tff, :f_noneq, :dTb_disc_min, :v_dmin, :dTb_disc_max,
                        :v_dmax, :dTb_cen_min, :dTb_cen_max, :dTb_core_min, :dTb_core_max, :S_core_nJy, :EW, :S_peak_nJy, :S_int_nJyHz, :theta_arcsec, :Tgamma)
summary_header() = "label,step,t_yr," * join(string.(SUMMARY_FIELDS), ",")
summary_row(label, step, t, sm) = @sprintf("%s,%d,%.8e,", label, step, t / YR) * join([@sprintf("%.6e", getfield(sm, k)) for k in SUMMARY_FIELDS], ",")

const PROFILE_HEADER = "label,j,r_AU,m_Msun,n_HI,T_k,T_s,T_gamma,x_c,x_cH,x_ce,x_cp,x_a,T_c,tau_sob,tau_th,tau21,t_eq_yr,t_eq_over_tff,dTb_ray_min_mK,dTb_ray_max_mK,v_ray_min_kms"
"Per-zone profile rows; the ray columns refer to the ray whose impact parameter is the zone's inner interface (rays are only available if `sp` is not nothing)."
function profile_rows(label, s, o, sp)
    n = length(o.Tk); rows = String[]
    for j in 1:n
        dmn = dmx = vmn = NaN
        if sp !== nothing
            i = findfirst(p -> p >= s.r[j] * (1 - 1e-12), sp.p)
            if i !== nothing
                d = sp.dT[i, :] ./ (1 + sp.z) .* 1e3; im = argmin(d); dmn = d[im]; dmx = maximum(d); vmn = sp.v[im] / 1e5
            end
        end
        push!(rows, @sprintf("%s,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e",
            label, j, 0.5 * (s.r[j] + s.r[j+1]) / AU, 0.5 * (s.m[j] + s.m[j+1]) / Msun, o.nHI[j], o.Tk[j], o.Ts_use[j], o.Tg[j], o.xc[j], o.xcH[j], o.xce[j], o.xcp[j], o.xa[j], o.Tc[j],
            o.tau_sob[j], o.tau_th[j], o.tau[j], o.teq[j] / YR, o.teq[j] / o.tff[j], dmn, dmx, vmn))
    end
    rows
end
"Spectrum file rows: v [km/s], disc-mean δT_b (observed, mK), the central ray and rays at p ≈ 1e-3, 1e-2, 0.1 R."
function spectrum_rows(sp)
    ps = [sp.p[argmin(abs.(sp.p .- f * sp.Rout))] for f in (1e-3, 1e-2, 0.1)]; ii = [findfirst(==(p), sp.p) for p in ps]
    f = 1e3 / (1 + sp.z); rows = String[]
    hdr = "v_kms,dTb_mean_mK,dTb_p0_mK," * join([@sprintf("dTb_p%.3gAU_mK", p / AU) for p in ps], ",")
    for k in eachindex(sp.v)
        push!(rows, @sprintf("%.5f,%.6e,%.6e,", sp.v[k] / 1e5, sp.mean[k] * f, sp.dT[1, k] * f) * join([@sprintf("%.6e", sp.dT[i, k] * f) for i in ii], ","))
    end
    hdr, rows
end

# ------------------------------------------------------------------ inline use from a running driver (popiii_m2c.jl, HI21=1)
mutable struct Inline
    fsum::IO; fprof::IO; run::String; Tcmb::Float64; z::Float64; rays::Bool; yprev::Any; tprev::Float64; neq::Bool; n::Int
end
const INLINE = Ref{Any}(nothing)
"""
    open_inline!(dir, run; Tcmb, z = Tcmb/2.7255 − 1, rays = true, neq = false)
Opens `<dir>/hi21_<run>.csv` and `…_prof.csv` for `inline!`. The Lyα input is `LYA_PROVIDER[]` (a function (snapshot) → `LyaIn`; default none).
"""
function open_inline!(dir, run; Tcmb, z = Tcmb / 2.7255 - 1, rays = true, neq = false)
    mkpath(dir)
    fs = open(joinpath(dir, "hi21_$(run).csv"), "w"); println(fs, summary_header())
    fp = open(joinpath(dir, "hi21_$(run)_prof.csv"), "w"); println(fp, PROFILE_HEADER)
    INLINE[] = Inline(fs, fp, run, Tcmb, z, rays, nothing, 0.0, neq, 0)
end
"One record from the running driver: step `k`, time `t` [s], interface radii/velocities/masses (n+1), zone ρ, T, y (as the SNAPS of popiii_m2c.jl)."
function inline!(k, t, r, u, ρ, T, y, m)
    I = INLINE[]; I === nothing && return nothing
    t0 = time(); s = (; step = k, t, r, u, ρ, T, y, m)
    lya = LYA_PROVIDER[] === nothing ? nothing : LYA_PROVIDER[](s)
    o, sp, sm = analyze(s; Tcmb = I.Tcmb, z = I.z, lya = lya, neq = I.neq ? I.yprev : nothing, dt = I.yprev === nothing ? 0.0 : t - I.tprev, rays = I.rays)
    I.yprev = o.y; I.tprev = t; I.n += 1; label = "$(I.run)_$(I.n)"
    println(I.fsum, summary_row(label, k, t, sm)); foreach(l -> println(I.fprof, l), profile_rows(label, s, o, sp)); flush(I.fsum); flush(I.fprof)
    @printf("HI21 step %d: <Ts>_m = %.1f K, max tau21 = %.2e, dTb(disc) %.2f .. %.2f mK, %.1f s\n", k, sm.Ts_mw, sm.tau_max, sm.dTb_disc_min, sm.dTb_disc_max, time() - t0)
    nothing
end

end # module
