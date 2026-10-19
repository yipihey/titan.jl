"""
    PrimordialChem

Population III project, M2: non-equilibrium chemistry of metal-free gas (H, H⁺, H⁻, H₂, H₂⁺, e⁻; He neutral with
n_He/n_H = 1/12 as RHFC02 App. A1.2), the equation of state of RHFC02 App. A1.3, and a one-zone free-fall driver.

Abundances are x_i = n_i/n_H (n_H = number density of H nuclei). H nuclei and charge are conserved by construction:
the unknowns are y = (x_H⁺, x_H⁻, x_H₂, x_H₂⁺), and
  x_H = 1 − x_H⁺ − x_H⁻ − 2x_H₂ − 2x_H₂⁺,   x_e = x_H⁺ + x_H₂⁺ − x_H⁻.
Rates come from a `RateSet` (a function k(T) → named tuple of rate coefficients); the network itself is written once.
All functions are generic in the number type (ForwardDiff Duals pass through).
"""
module PrimordialChem

using Printf

const kB = 1.380649e-16; const mH = 1.6735575e-24; const eV = 1.602176634e-12; const G = 6.67430e-8
const yHe = 1 / 12                                  # n_He / n_H (RHFC02 App. A1.2)
const χH2 = 4.48eV; const χH = 13.598eV; const χHm = 0.754eV; const χH2p = 2.65eV   # binding energies (see notes)

"Energy offset: u + UOFF_TH = energy relative to cold H₂ (+ neutral He), positive (SCVH95's zero)."
const UOFF_TH = (χH + χH2 / 2) / (mH * (1 + 4yHe))
"Mass per H nucleus [g] (H + He, metal free)."
const μH = mH * (1 + 4yHe)

xH(y) = 1 - y[1] - y[2] - 2y[3] - 2y[4]
xe(y) = y[1] + y[4] - y[2]
"Particles per H nucleus (H, H⁺, H⁻, H₂, H₂⁺, e⁻, He)."
ntot(y) = xH(y) + y[1] + y[2] + y[3] + y[4] + xe(y) + yHe

"""
Rate equations dy/dt for y = (x_H⁺, x_H⁻, x_H₂, x_H₂⁺) at H-nucleus density nH and rate coefficients k (named tuple,
cgs; missing reactions are zero). Reaction labels follow the rate sets (see `rates_*`).
"""
function dydt(y, nH, k)
    H = xH(y); Hp, Hm, H2, H2p = y; e = xe(y)
    g(s) = hasproperty(k, s) ? getproperty(k, s) : zero(eltype(y))
    # rates per H nucleus per unit time (x_i x_j n_H k for two-body, x_i x_j x_k n_H² k for three-body)
    r_ion   = g(:ion_e) * H * e * nH            # H + e → H⁺ + 2e
    r_rec   = g(:rec) * Hp * e * nH             # H⁺ + e → H + γ
    r_Hm    = g(:Hm_form) * H * e * nH          # H + e → H⁻ + γ
    r_H2Hm  = g(:H2_Hm) * H * Hm * nH           # H + H⁻ → H₂ + e
    r_Hmn   = g(:Hm_Hp) * Hm * Hp * nH          # H⁻ + H⁺ → 2H
    r_Hmd   = g(:Hm_e) * Hm * e * nH            # H⁻ + e → H + 2e
    r_H2p   = g(:H2p_form) * H * Hp * nH        # H + H⁺ → H₂⁺ + γ
    r_H2H2p = g(:H2_H2p) * H2p * H * nH         # H₂⁺ + H → H₂ + H⁺
    r_H2pH2 = g(:H2p_H2) * H2 * Hp * nH         # H₂ + H⁺ → H₂⁺ + H
    r_3b    = g(:tb_H) * H^3 * nH^2             # 3H → H₂ + H
    r_3b2   = g(:tb_H2) * H^2 * H2 * nH^2       # 2H + H₂ → 2H₂
    r_dH    = g(:cd_H) * H2 * H * nH            # H₂ + H → 3H
    r_dH2   = g(:cd_H2) * H2 * H2 * nH          # H₂ + H₂ → 2H + H₂
    r_HH    = g(:ion_HH) * H * H * nH           # H + H → H⁺ + H + e
    dHp  = r_ion + r_HH - r_rec - r_Hmn - r_H2p + r_H2H2p - r_H2pH2
    dHm  = r_Hm - r_H2Hm - r_Hmn - r_Hmd
    dH2  = r_H2Hm + r_H2H2p - r_H2pH2 + r_3b + r_3b2 - r_dH - r_dH2
    dH2p = r_H2p - r_H2H2p + r_H2pH2
    # (chemical heating needs no term here: the binding energies are part of the internal energy, eos_rhfc02)
    (dHp, dHm, dH2, dH2p)
end

"""
Backward-Euler step of the network over dt at fixed (nH, T), Newton on the 4 unknowns with a damped update that keeps
every abundance non-negative. Returns the new y. Recursive interval halving on failure.
"""
function step(y0, nH, T, dt, rates; tol = 1e-10, maxit = 50)
    k = rates(T)
    function attempt(y0, h)
        y = collect(y0)
        for it in 1:maxit
            f = collect(dydt(y, nH, k))
            R = y .- collect(y0) .- h .* f
            J = zeros(eltype(y), 4, 4); ε = 1e-8
            for j in 1:4
                yp = copy(y); d = ε * max(abs(y[j]), 1e-20); yp[j] += d
                J[:, j] = (yp .- collect(y0) .- h .* collect(dydt(yp, nH, k)) .- R) ./ d
            end
            δ = -(J \ R)
            # damping: no abundance below zero, including the eliminated x_H and x_e (linear in y)
            λ = 1.0
            for j in 1:4; y[j] + λ * δ[j] < 0 && (λ = min(λ, 0.9 * y[j] / -δ[j])); end
            for (f, df) in ((xH(y), -δ[1] - δ[2] - 2δ[3] - 2δ[4]), (xe(y), δ[1] + δ[4] - δ[2]))
                f + λ * df < 0 && (λ = min(λ, 0.9 * f / -df))
            end
            y .+= λ .* δ
            maximum(abs.(λ .* δ) ./ max.(abs.(y), 1e-20)) < tol && λ == 1.0 && return (y, true)
        end
        (y, false)
    end
    # adaptive: on failure split the interval in halves, recursively (backward Euler is L-stable; the halves only help
    # Newton, by keeping the start point closer to the solution)
    function adv(y0, h, depth)
        y, ok = attempt(y0, h)
        ok && return Tuple(y)
        depth >= 30 && error("PrimordialChem.step: Newton failed (nH = $nH, T = $T, dt = $dt)")
        adv(adv(y0, h / 2, depth + 1), h / 2, depth + 1)
    end
    adv(y0, dt, 0)
end

"""
Chemistry and temperature together over dt at fixed ρ and fixed internal energy u (thermal + chemical): unknowns y and
ln T; residuals y − y0 − dt·dydt(y, n_H, k(T)) and (u(ρ, T, y) − u + dt Λ(ρ, T, y)/ρ)/u_s, u_s = k T/μ_H (`cooling`
= Λ [erg cm⁻³ s⁻¹], optional, backward Euler). Damped Newton (abundances, x_H,
x_e ≥ 0), recursive interval halving on failure. The binding energy released or absorbed by the reactions changes T
inside the same solve: no splitting error between chemistry and heating (dissociation at T ≳ 2000 K is limited by
the thermal energy available, which an operator split with rates at a fixed T overshoots).
Returns (y, T).
"""
function step_coupled(y0, ρ, u, T0, dt, rates; eos = eos_rhfc02, cooling = nothing, tol = 1e-10, maxit = 60)
    nH = ρ / μH
    function resid(z, y0, h)
        y = (z[1], z[2], z[3], z[4]); T = exp(z[5])
        f = dydt(y, nH, rates(T))
        [y[1] - y0[1] - h * f[1], y[2] - y0[2] - h * f[2], y[3] - y0[3] - h * f[3], y[4] - y0[4] - h * f[4],
         (eos(ρ, T, y)[1] - u + (cooling === nothing ? 0.0 : h * cooling(ρ, T, y) / ρ)) / (kB * T / μH)]
    end
    function attempt(y0, T0, h)
        z = [y0[1], y0[2], y0[3], y0[4], log(T0)]
        for it in 1:maxit
            R = resid(z, y0, h); J = zeros(5, 5)
            for j in 1:5
                zp = copy(z); d = j == 5 ? 1e-7 : 1e-7 * max(abs(z[j]), 1e-20); zp[j] += d
                J[:, j] = (resid(zp, y0, h) .- R) ./ d
            end
            δ = -(J \ R); all(isfinite, δ) || return (z, false)
            λ = 1.0
            for j in 1:4; z[j] + λ * δ[j] < 0 && (λ = min(λ, 0.9 * z[j] / -δ[j])); end
            y = (z[1], z[2], z[3], z[4])
            for (f, df) in ((xH(y), -δ[1] - δ[2] - 2δ[3] - 2δ[4]), (xe(y), δ[1] + δ[4] - δ[2]))
                f + λ * df < 0 && (λ = min(λ, 0.9 * f / -df))
            end
            abs(δ[5]) * λ > 0.5 && (λ = 0.5 / abs(δ[5]))          # at most a factor e^0.5 in T per iteration
            z .+= λ .* δ
            λ == 1.0 && maximum(abs.(δ[1:4]) ./ max.(abs.(z[1:4]), 1e-20)) < tol && abs(δ[5]) < tol && return (z, true)
        end
        (z, false)
    end
    function adv(y0, T0, h, depth)
        z, ok = attempt(y0, T0, h)
        ok && return ((z[1], z[2], z[3], z[4]), exp(z[5]))
        depth >= 30 && error("PrimordialChem.step_coupled: Newton failed (nH = $nH, T0 = $T0, dt = $dt)")
        y1, T1 = adv(y0, T0, h / 2, depth + 1)
        adv(y1, T1, h / 2, depth + 1)
    end
    adv(y0, T0, dt, 0)
end

"""
RHFC02 (App. A1.3, Eqs. A6–A9) internal energy per unit mass: u = N_tot k T/(γ_ad − 1) − N_tot Σ y_i χ_i, with
γ_ad the number-weighted mean of the species' γ (5/3 monatomic, 7/5 for H₂⁺, H₂ from Eq. A9:
1/(γ_H₂ − 1) = 5/2 + x² eˣ/(eˣ − 1)², x = 6100 K/T). Chemical energy relative to free H atoms, protons and electrons
at rest: −χ per bound particle (H₂: 4.48 eV; H: 13.6 eV relative to p + e; H⁻: 0.754 eV relative to H + e; H₂⁺).
**Interpretation [design]:** RHFC02 write −N_tot Σ y_i χ_i with y_i = N_i/N_tot, i.e. −Σ N_i χ_i; the zero of χ
(here: fully ionised H, free electrons) is our choice; only differences matter.
Returns (u [erg/g], p [dyn cm⁻²]) at ρ, T, y.
"""
function eos_rhfc02(ρ, T, y)
    nH = ρ / μH; N = ntot(y)                         # particles per H nucleus
    H = xH(y); Hp, Hm, H2, H2p = y; e = xe(y)
    xx = 6100 / T
    inv_g_H2 = 5 / 2 + xx^2 * exp(xx) / expm1(xx)^2
    # number-weighted mean of 1/(γ−1) — "weighted average of the adiabatic indexes" (A1.3); averaging 1/(γ−1) is
    # the one that makes u_th = Σ N_i kT/(γ_i − 1) [interpretation, see notes]
    sum_inv = (H + Hp + Hm + e + yHe) * 1.5 + H2 * inv_g_H2 + H2p * 2.5
    uth = sum_inv * kB * T / μH
    uchem = -(H * χH + Hm * (χH + χHm) + H2 * (2χH + χH2) + H2p * (2χH - χH + χH2p)) / μH   # see docstring
    p = N * nH * kB * T
    (uth + uchem, p)
end

"T from u at fixed ρ, y (bisection in ln T, as RHFC02 A2.1 (bisection); u is monotone in T at fixed y)."
function T_of_u(u, ρ, y; T0 = 300.0, eos = eos_rhfc02)
    lo, hi = T0 / 2, T0 * 2
    while eos(ρ, lo, y)[1] > u; lo /= 4; lo < 1e-6 && error("T_of_u: u below the T → 0 limit"); end
    while eos(ρ, hi, y)[1] < u; hi *= 4; hi > 1e12 && error("T_of_u: no bracket"); end
    for _ in 1:200
        m = sqrt(lo * hi); eos(ρ, m, y)[1] > u ? (hi = m) : (lo = m)
        hi / lo - 1 < 1e-14 && break
    end
    sqrt(lo * hi)
end

# ------------------------------------------------------------------ rate sets
# Sources and checks: ~/Projects/popiii_lit/rates_extraction.md (AAZN97 = Abel, Anninos, Zhang & Norman 1997, New Astron.
# 2, 181, appendix tables; PSS83 = Palla, Salpeter & Stahler 1983 as quoted by RA04 eq. 31 and RA05 Table 1.1, agreeing
# with the ADS OCR of PSS83 Table 1). T in K unless noted; AAZN97 fits marked "eV" take T_eV = T/11605 (AAZN97's factor).
const KeV = 11605.0
_lnpoly(lT, c) = (s = zero(lT); for (i, a) in enumerate(c); s += a * lT^(i - 1); end; s)
"AAZN97 k1, H + e → H⁺ + 2e (Janev et al. 1987), T in eV; zero below 1000 K as AAZN97."
k1_aazn(T) = T < 1000 ? zero(T) : exp(_lnpoly(log(T / KeV), (-32.71396786, 13.536556, -5.73932875, 1.56315498, -0.2877056,
    3.48255977e-2, -2.63197617e-3, 1.11954395e-4, -2.03914985e-6)))
"AAZN97 k2, H⁺ + e → H + γ, case A (fit to Ferland et al. 1992), T in eV."
k2_aazn(T) = exp(_lnpoly(log(T / KeV), (-28.6130338, -0.72411256, -2.02604473e-2, -2.38086188e-3, -3.21260521e-4,
    -1.42150291e-5, 4.98910892e-6, 5.75561414e-7, -1.85676704e-8, -3.07113524e-9)))
"AAZN97 k7, H + e → H⁻ + γ (Wishart 1979), T in K."
function k7_aazn(T)
    l = log10(T)
    T <= 6000 ? 1.429e-18 * T^0.7620 * T^(0.1523l) * T^(-3.274e-2 * l^2) :
                3.802e-17 * T^(0.1998l) * 10^(4.0415e-5 * l^6 - 5.447e-3 * l^4)
end
"AAZN97 k8, H + H⁻ → H₂ + e (Janev et al. 1987), T in eV above 0.1 eV, 1.428e-9 below."
k8_aazn(T) = T / KeV < 0.1 ? 1.428e-9 + zero(T) : exp(_lnpoly(log(T / KeV), (-20.06913897, 0.22898, 3.5998377e-2, -4.55512e-3,
    -3.10511544e-4, 1.0732940e-4, -8.36671960e-6, 2.23830623e-7)))
"PSS83 k4, 3H → H₂ + H [cm⁶ s⁻¹] (RA04 eq. 31, RA05 k8)."
k4_pss(T) = 5.5e-29 / T
"PSS83 k5, H₂ + H → 3H [cm³ s⁻¹] (RA05 k9, printed cm⁶ there: typo; PSS83 OCR cm³). Detailed-balance partner of k4."
k5_pss(T) = 6.5e-7 / sqrt(T) * exp(-52000 / T) * (1 - exp(-6000 / T))
"RA04's modification of k5 (RA04 eq. 32), to match Martin, Schwarz & Mandy 1996 at n ≳ 1e9."
k5_ra04(T) = 2.2e-9 * T^0.2 * exp(-51800 / T) * (1 - exp(-6000 / T))
"PSS83 k9, H + H → H⁺ + e + H = 1.7e-4 × PSS83 k8 (Lotz 1967: 5.8e-11 T^½ e^(−158000/T)) (RA05 k12, k13)."
k9_pss(T) = 1.7e-4 * 5.8e-11 * sqrt(T) * exp(-158000 / T)

const hP = 6.62607015e-27; const D0_H2 = 4.478eV; const θrot_H2 = 85.3; const θvib_H2 = 6332.0
"""
Equilibrium constant of H₂ ⇌ 2H, K = n_H²/n_H₂ [cm⁻³] = 16/(Z_rot Z_vib) (π m_H k T/h²)^{3/2} e^{−D0/kT}: rigid rotor
(θ_rot = 85.3 K, ortho:para nuclear weights 3:1), harmonic vibration (θ_vib = 6332 K), D0 = 4.478 eV, H electron spin
2 per atom. (rates_extraction.md §C.3.) The PSS83 pair k5/k4 equals 0.25 K to 2 % over 1000–5000 K [run]: PSS83's
detailed balance lacks a factor 4 (probably the H spin degeneracy [interp]).
"""
function K_H2(T)
    Zr = zero(T); for J in 0:60; Zr += (iseven(J) ? 1 : 3) * (2J + 1) * exp(-J * (J + 1) * θrot_H2 / T); end
    Zv = 1 / (1 - exp(-θvib_H2 / T))
    16 / (Zr * Zv) * (π * mH * kB * T / hP^2)^1.5 * exp(-D0_H2 / (kB * T))
end
"Dissociation rate in detailed balance with a formation rate `tb` (default PSS83 k4): k5 = k4 K_H₂ [cm³ s⁻¹]."
k5_db(T; tb = k4_pss) = tb(T) * K_H2(T)

const me = 9.1093837e-28
"Saha constant of H ⇌ H⁺ + e⁻: n_H⁺ n_e/n_H = (2π m_e k T/h²)^{3/2} e^{−13.6 eV/kT} (statistical weights 1·2/2) [cm⁻³]."
K_ion(T) = (2π * me * kB * T / hP^2)^1.5 * exp(-χH / (kB * T))
"⁷Li/H by number, SBBN (Coc et al. 2004), as Mayer & Duschl 2005 Table 1."
const yLi = 4.15e-10
"Saha constant of Li ⇌ Li⁺ + e⁻ (χ = 5.3917 eV; g(Li⁺) g(e)/g(Li) = 1·2/2) [cm⁻³]."
K_Li(T) = (2π * me * kB * T / hP^2)^1.5 * exp(-5.3917eV / (kB * T))
"""
Total free-electron density [cm⁻³] with lithium: n_e = n_e,H + n_Li⁺, n_Li⁺ = n_Li K_Li/(K_Li + n_e) (Saha; collisional
ionisation equilibrium is fast at the densities where it matters). Li supplies the electrons in metal-free gas below
≈ 3500 K at high density (Mayer & Duschl 2005 §2.3; M3c: 15× H's at 3000 K, 1e18 cm⁻³); negligible for p and u
(Li/H = 4e-10) and, at the network's temperatures (≲ 1600 K), below the network's own x_e.
"""
function ne_with_Li(nH, T, neH)
    nLi = yLi * nH; K = K_Li(T); ne = neH
    for _ in 1:60
        ne_new = neH + nLi * K / (K + ne + 1e-300)
        abs(ne_new - ne) <= 1e-12 * ne_new && return ne_new
        ne = ne_new
    end
    ne
end

"""
Equilibrium composition at (ρ, T) (RHFC02 A1.2 for ρ ≥ ρ_eq, p ≤ 1e4 dyn cm⁻²: Saha equations): H₂ ⇌ 2H (K_H2) and
H ⇌ H⁺ + e⁻ (K_ion); H⁻ and H₂⁺ set to 0; He neutral. x_H² n_H = K_H2 x_H₂, x_H⁺² n_H = K_ion x_H, x_H + x_H⁺ + 2x_H₂ = 1.
Solved for x_H by bisection (monotone). Returns y = (x_H⁺, 0, x_H₂, 0).
**Not yet:** SCVH95 tables for p ≥ 1e6 dyn cm⁻² (pressure ionisation, non-ideal effects).
"""
function y_saha(ρ, T)
    nH = ρ / μH; A = K_H2(T) / nH; B = K_ion(T) / nH
    A == 0 && return (0.0, 0.0, 0.5, 0.0)                 # cold: fully molecular
    g(x) = x + sqrt(B * x) + 2x^2 / A - 1                 # x = x_H
    lo, hi = 0.0, 1.0
    for _ in 1:200; m = (lo + hi) / 2; g(m) > 0 ? (hi = m) : (lo = m); end
    x = (lo + hi) / 2
    xp = sqrt(B * x)
    (xp, 0.0, max(1 - x - xp, 0.0) / 2, 0.0)        # x_H₂ from conservation (x²/A underflows when K_H2 → 0 at low T)
end
"""
Equilibrium-chemistry step at fixed ρ and u over dt (RHFC02 above ρ_eq): find T with u(ρ, T, y_saha(ρ, T)) =
u − dt Λ(ρ, T, y)/ρ (bisection in ln T). Returns (y, T).
"""
function step_equilibrium(ρ, u, T0, dt; eos = eos_rhfc02, cooling = nothing)
    F(T) = (y = y_saha(ρ, T); eos(ρ, T, y)[1] - u + (cooling === nothing ? 0.0 : dt * cooling(ρ, T, y) / ρ))
    lo, hi = T0 / 1.5, T0 * 1.5
    while F(lo) > 0; lo /= 1.5; lo < 1.0 && error("step_equilibrium: no bracket below"); end
    while F(hi) < 0; hi *= 1.5; hi > 1e9 && error("step_equilibrium: no bracket above"); end
    for _ in 1:200; m = sqrt(lo * hi); F(m) > 0 ? (hi = m) : (lo = m); hi / lo - 1 < 1e-13 && break; end
    T = sqrt(lo * hi); (y_saha(ρ, T), T)
end

"""
RHFC02 Table A1 ("faithful"): AAZN97 k1, k2 (case A), k7, k8; PSS83 k4, k5, k6 = k4/8, k7 = k5/8, k9. No H₂⁺ and no
He reactions (RHFC02 lists none). `k5` selects the dissociation rate (default PSS83; `k5_ra04` for RA04's).
"""
rates_rhfc02(; k5 = k5_pss, tb = k4_pss) = T -> (ion_e = k1_aazn(T), rec = k2_aazn(T), Hm_form = k7_aazn(T), H2_Hm = k8_aazn(T),
    tb_H = tb(T), tb_H2 = tb(T) / 8, cd_H = k5(T), cd_H2 = k5(T) / 8, ion_HH = k9_pss(T))

# ------------------------------------------------------------------ H₂ line cooling
"Hollenbach & McKee 1979 LTE H₂ cooling per molecule [erg s⁻¹] (RA05 eq. 1.23 = RA04 eq. 7)."
function Λ_H2_LTE_HM(T)
    T3 = T / 1000
    9.5e-22 * T3^3.76 / (1 + 0.12 * T3^2.1) * exp(-(0.13 / T3)^3) + 3e-24 * exp(-0.51 / T3) +
        6.7e-19 * exp(-5.86 / T3) + 1.6e-18 * exp(-11.7 / T3)
end
"Galli & Palla 1998 eq. A.7, H₂–H low-density limit per molecule per H atom [erg cm³ s⁻¹], 10 ≤ T ≤ 1e4 K (clamped)."
function Λ_H2_0_GP(T)
    l = log10(clamp(T, 10.0, 1e4))
    10^(-103.0 + 97.59l - 48.05l^2 + 10.80l^3 - 0.9032l^4)
end
"""
Thin H₂ cooling [erg cm⁻³ s⁻¹], RA05 eq. 1.22 (GP98 A.5–A.6 read consistently): per molecule
Λ_LTE/(1 + Λ_LTE/(n_H Λ_0)) with n_H the atomic-H density, times n_H₂.
"""
function cool_H2_thin(ρ, T, y)
    nH = ρ / μH; nHat = xH(y) * nH; nH2 = y[3] * nH
    L = Λ_H2_LTE_HM(T); L0 = Λ_H2_0_GP(T)
    nH2 * L / (1 + L / max(nHat * L0, 1e-300))
end
"""
H₂ cooling with a radiation-bath floor: Λ(T) − Λ(T_bath) at the same density and composition, so the gas cannot cool
below T_bath (here the CMB, T_bath = 2.725 (1 + z) K; ABN02's start is z = 19: 54.5 K). [design] Enzo/GA08-style
convention; no source in the RHFC02 set (RHFC02 starts at T ≫ T_CMB).
"""
cool_floor(f, Tb) = (ρ, T, y) -> f(ρ, T, y) - f(ρ, Tb, y)
"RA04 eq. 22: optically thick H₂ lines ≈ thin × min(1, (n/n0)^−0.45), n = ρ/m_H, n0 = 8e9 cm⁻³."
cool_H2_ra04(ρ, T, y) = cool_H2_thin(ρ, T, y) * min(1.0, (ρ / mH / 8e9)^-0.45)
"""
RA04 CIE cooling (eq. 5 with RA04's escape factor, §3.3 / eq. A1): ρ L_CIE (1 − e^−τ)/τ, L_CIE = A ρ T⁴ X f_H2
[erg g⁻¹ s⁻¹], τ = (ρ/3.3e-8 g cm⁻³)^2.8; A = 0.072 (eq. 5), 0.054 in the escape-factor fit to RHFC02 (A1).
f_H2 = fraction of H nuclei in H₂ = 2 x_H₂; X = 1/(1 + 4 y_He).
"""
function cool_CIE_ra04(ρ, T, y; A = 0.054)
    X = 1 / (1 + 4yHe); τ = (ρ / 3.3e-8)^2.8
    esc = τ < 1e-6 ? 1.0 - τ / 2 : -expm1(-τ) / τ
    ρ * A * ρ * T^4 * X * 2y[3] * esc
end

"""
    onezone(; nH0, T0, y0, rates, cooling, nHmax, dlnρ = 0.01, eos = eos_rhfc02, out = nothing)

One-zone free-fall collapse (Omukai 2000 style): dρ/dt = ρ/t_ff, t_ff = (3π/32Gρ)^½; du/dt = (p/ρ²) dρ/dt − Λ/ρ with
u the full internal energy (thermal + chemical; chemical heating/cooling enters through u_chem at fixed u). Steps of
Δln ρ = dlnρ; each step: (1) compression + cooling over dt (trapezoid in p, backward Euler in Λ, solved for T), (2) chemistry and T solved together over dt at fixed u (`step_coupled`). `cooling(ρ, T, y)` returns Λ [erg cm⁻³ s⁻¹]. `out(ρ, T, y, t)` is called every step.
Returns the vector of (nH, T, y).
"""
function onezone(; nH0, T0, y0, rates, cooling, nHmax, dlnρ = 0.01, eos = eos_rhfc02, out = nothing)
    ρ = nH0 * μH; T = T0; y = y0; t = 0.0
    u = eos(ρ, T, y)[1]
    track = [(nH0, T, y)]
    while ρ / μH < nHmax
        tff = sqrt(3π / (32G * ρ)); dt = dlnρ * tff
        ρn = ρ * exp(dlnρ)
        # (1) compression (p dV, trapezoid in p) and cooling, implicit in T (backward Euler for Λ: the cooling time can be
        #     ≪ the step): solve F(T) = u(ρn, T, y) − u − p̄ (1/ρ − 1/ρn) + Λ(ρn, T, y)/ρn dt = 0 in ln T, bracketing
        p0 = eos(ρ, T, y)[2]
        F(Tx) = eos(ρn, Tx, y)[1] - u - 0.5 * (p0 + eos(ρn, Tx, y)[2]) * (1 / ρ - 1 / ρn) + cooling(ρn, Tx, y) / ρn * dt
        lo, hi = T / 10, T * 10
        while F(lo) > 0; lo /= 10; lo < 1e-3 && error("onezone: no temperature (cooling)"); end
        while F(hi) < 0; hi *= 10; end
        for _ in 1:200                                  # bisection in ln T (robust; F is monotone in T in practice)
            m = sqrt(lo * hi); F(m) > 0 ? (hi = m) : (lo = m)
            hi / lo - 1 < 1e-13 && break
        end
        Tn = sqrt(lo * hi)
        un = eos(ρn, Tn, y)[1]
        # (2) chemistry and T together at fixed u (step_coupled)
        yn, Tn = step_coupled(y, ρn, un, Tn, dt, rates; eos)
        ρ, T, y, u, t = ρn, Tn, yn, un, t + dt
        push!(track, (ρ / μH, T, y))
        out === nothing || out(ρ, T, y, t)
    end
    track
end

end # module
