"""
    FreeEnergyEOS

Self-consistent analytical Helmholtz free-energy model for astrophysical H/He/CNO gas across
ρ = 10⁻²⁴ … 10⁴ g cm⁻³, T = 10 … 3×10⁹ K, and arbitrary composition (X, Y, Z), using the chemical picture:

    f_tot(T, {n_i}) = f_id,heavy + f_int+ub + f_id,e± + f_Coul + f_HS + f_rep     [erg cm⁻³]

Components (all C^∞/C² analytical, generic in number type, nested ForwardDiff-safe):
1. Ideal translational free energy of heavy species (Boltzmann, exact atomic/molecular masses).
2. Internal partition functions with self-consistent Potekhin (1996, Phys. Plasmas 3, 4156) quantum
   unbounding occupation probabilities w_n = exp(-n_c v_n) for H I (n ≤ 50), He I (NIST n ≤ 2 + n ≤ 35),
   He II (n ≤ 35), and ExoMol RACPPK 301-level H₂ (C² blended to exact rigid-rotor × harmonic below 800 K).
   Because Σ N_n ln w_n + F_ub/kT = 0 identically at level equilibrium, F_int+ub = Σ n_s (E0_s - kT ln Z_{w,s})
   satisfies all Maxwell relations automatically and produces smooth pressure ionization/dissociation.
3. Ideal Fermi–Dirac relativistic electron + e± pair gas (StellarEOS.electron_f): exact degeneracy chemical
   potential shift (Pauli blocking of recombination at high density) and pair creation at T ≳ 10⁸ K.
4. Non-ideal Coulomb + electron exchange-correlation + ion-electron screening:
   - f_ee(r_s, t): Karasiev, Sjostrom, Dufty & Trickey (2014, PRL 112, 076403) RPIMC/QMC fit (exact DH limit).
   - f_ii(Γ_i): Potekhin & Chabrier (2000, PRE 62, 8554) Caillol MC + DH + Abe fit (C^∞ for all Γ).
   - f_ie(Γ_e, r_s, Z_i): Potekhin & Chabrier (2000, PRE 62, 8554) finite-T HNC + relativistic TF fit.
   Combined via the linear mixing rule (Chabrier & Potekhin 1998, §VI): exact two-component Debye–Hückel limit
   at Γ ≪ 1 and exact dense screened OCP at Γ ≫ 1.
5. Dense neutral fluid interactions (Ross 1979 / Saumon & Chabrier 1991 / Potekhin 1996):
   - Entropic excluded volume: exact Boublík–Mansoori–Carnahan–Starling–Leland (BMCSL, Mansoori et al. 1971)
     multi-component hard-sphere free energy for H₂, H, He with Barker–Henderson soft-sphere diameters d_B(T)
     smoothly density-compressed via s(η_0) = (1 + (η_0/η_*)²)^(-1/6) so η < η_* < 1 strictly for all ρ.
   - Short-range soft-core repulsive energy f_rep = ε_rep ξ_0 η_c² (T-independent quantum core repulsion),
     reproducing the dense cold molecular H₂ compression curve and smooth pressure metallization at ρ ~ 1–2 g cm⁻³
     while remaining strictly convex (χ_ρ > 0) at all densities.

All thermodynamic quantities derive strictly from f_tot via ForwardDiff + exact Implicit Function Theorem:
  μ_i = ∂f_tot/∂n_i,   p = Σ n_i μ_i - f_tot,   s = -(∂f_tot/∂T)/ρ,   u = (f_tot + T ρ s)/ρ.
Non-ideal chemical potential shifts Δμ_i^ex = ∂f_ex/∂n_i feed directly into ChemNet's detailed balance
ln K_eff = ln K_ideal - Σ Δμ_i^ex / (k_B T), so non-LTE kinetics relax to the exact minimum of f_tot.
"""
module FreeEnergyEOS

import ..Z0Physics, ..StellarEOS
using ForwardDiff, LinearAlgebra, Printf, Serialization

const Z0 = Z0Physics
const SE = StellarEOS
const kB = SE.kB; const hP = SE.h; const ħ = hP / (2π); const c_light = SE.cl
const eV = SE.eV; const qe = SE.qe; const me = SE.me; const mH = SE.mH; const mHe = SE.mHe
const amu = SE.amu; const a0 = 5.29177211e-9; const Ha = qe^2 / a0      # 1 Hartree = 27.211386 eV
const Ry = 0.5 * Ha
const mp = 1.67262192e-24

# Ionization / dissociation energies [erg] (matching ChemNet energy zero: H⁺ = 0, He = 0)
const χH   = 13.598434 * eV
const EA_H = 0.754 * eV
const D0_H2 = 4.478 * eV
const D0_H2p = 2.65 * eV
const χHe1 = 24.587387 * eV
const χHe2 = 54.417763 * eV
const E0_H    = -χH
const E0_Hp   = 0.0
const E0_H2   = -(2χH + D0_H2)
const E0_Hm   = -(χH + EA_H)
const E0_H2p  = -(χH + D0_H2p)
const E0_He   = 0.0
const E0_Hep  = χHe1
const E0_Hepp = χHe1 + χHe2

# ---------------------------------------------------------------------------------------------------------------------
# 1. Exact Potekhin & Chabrier (2000, PRE 62, 8554) ion-ion OCP & ion-electron screening free energies
# ---------------------------------------------------------------------------------------------------------------------
const PC_A1 = -0.907347
const PC_A2 = 0.62849
const PC_A3 = -sqrt(3.0) / 2 - PC_A1 / sqrt(PC_A2)
const PC_B1 = 4.50e-3
const PC_B2 = 170.0
const PC_B3 = -8.4e-5
const PC_B4 = 3.70e-3
const PC_CINF = (18.0 / 175.0) * (12.0 / π)^(2.0 / 3.0)

"Potekhin & Chabrier (2000) eq. (16): liquid OCP reduced excess free energy f_ii(Γ) = F_ii / (N_i k_B T). C^∞ for Γ > 0."
@inline function f_ii_pc00(Γ)
    sG = sqrt(Γ)
    if Z0.pval(Γ) < 1e-6
        t = Γ / PC_A2
        term1 = PC_A1 * (sG * Γ / sqrt(PC_A2)) * (2.0 / 3.0 - t * (0.2 - t * (3.0 / 28.0)))
        term2 = 2 * PC_A3 * (sG * Γ) * (1.0 / 3.0 - Γ * (0.2 - Γ / 7.0))
        term3 = PC_B1 * (Γ * Γ / (2 * PC_B2))
        term4 = 0.5 * PC_B3 * log1p(Γ * Γ / PC_B4)
        return term1 + term2 + term3 + term4
    end
    q = sG / sqrt(PC_A2)
    t1 = PC_A1 * (sqrt(Γ * (PC_A2 + Γ)) - PC_A2 * asinh(q))
    t2 = 2 * PC_A3 * (sG - atan(sG))
    t3 = PC_B1 * (Γ - PC_B2 * log1p(Γ / PC_B2))
    t4 = 0.5 * PC_B3 * log1p(Γ * Γ / PC_B4)
    return t1 + t2 + t3 + t4
end

"Potekhin & Chabrier (2000) eqs. (22)-(24): ion-electron screening reduced free energy f_ie(Γ_e, r_s, x_r, Z)."
@inline function f_ie_pc00(Γe, rs, xr, Z::Float64)
    cDH = (Z / sqrt(3.0)) * ((1 + Z)^1.5 - 1 - Z^1.5)
    Z13 = cbrt(Z)
    cTF = PC_CINF * Z * Z * Z13 * (1 - 1 / Z13 + 0.2 / sqrt(Z))
    a = 1.11 * Z^0.475
    lnZ = log(Z)
    b = 0.2 + 0.078 * lnZ * lnZ
    ν = 1.16 + 0.08 * lnZ
    sGe = sqrt(Γe)
    Geν = Γe^ν
    zrs = Z / rs
    g1 = 1 + 0.78 / (21 + Γe * zrs * zrs * zrs) * sqrt(Γe / Z)
    rs2_over_den = Z0.pval(rs) > 1.0 ? 1.0 / (6.0 + 1.0 / (rs * rs)) : (rs * rs) / (1.0 + 6.0 * rs * rs)
    g2_over_rs = 1.0 / rs + ((Z - 1.0) / 9.0) * (1.0 + 1.0 / (0.001 * Z * Z + 2 * Γe)) * rs2_over_den
    xr2 = xr * xr
    if Z0.pval(xr) > 10.0
        ix = 1.0 / xr
        h1 = (0.2 + ix * ix) / (0.2 + 0.37 / sqrt(Z) + 0.18 * ix / sqrt(sqrt(Z)) + ix * ix)
    else
        h1 = (1.0 + 0.2 * xr2) / (1.0 + 0.18 * xr / sqrt(sqrt(Z)) + 0.37 * xr2 / sqrt(Z) + 0.2 * xr2)
    end
    h2 = 1.0 / sqrt(1.0 + xr2)
    num = cDH * sGe + cTF * a * Geν * g1 * h1
    den = 1.0 + (b * sGe + a * Geν * g2_over_rs) * h2
    return -Γe * num / den
end

# ---------------------------------------------------------------------------------------------------------------------
# 2. Exact Karasiev, Sjostrom, Dufty & Trickey (2014, PRL 112, 076403) electron exchange-correlation f_ee(r_s, t)
# ---------------------------------------------------------------------------------------------------------------------
const K_LAM = cbrt(4.0 / (9π))
const K_B1 = 0.283997; const K_B2 = 48.932154; const K_B3 = 0.370919; const K_B4 = 61.095357
const K_B5 = sqrt(1.5) / K_LAM * K_B3
const K_C1 = 0.870089; const K_C2 = 0.193077;  const K_C3 = 2.414644
const K_D1 = 0.579824; const K_D2 = 94.537454; const K_D3 = 97.839603; const K_D4 = 59.939999; const K_D5 = 24.388037
const K_E1 = 0.212036; const K_E2 = 16.731249; const K_E3 = 28.485792; const K_E4 = 34.028876; const K_E5 = 17.235515

@inline function f_ee_ksdt(rs, t, Γe)
    tv = Z0.pval(t)
    if tv > 1.0
        u = 1.0 / t; u2 = u * u; u4 = u2 * u2
        su = sqrt(u)
        tu = tanh(u); tsu = tanh(su)
        a = 0.610887 * tu * (1.7035 - 0.09227 * u + 3.04363 * u2 + 0.75 * u4) / (5.1105 + 8.31051 * u2 + u4)
        b = tsu * (K_B3 + K_B2 * u2 + K_B1 * u4) / (K_B5 + K_B4 * u2 + u4)
        d = tsu * (K_D3 + K_D2 * u2 + K_D1 * u4) / (K_D5 + K_D4 * u2 + u4)
        e = tu * (K_E3 + K_E2 * u2 + K_E1 * u4) / (K_E5 + K_E4 * u2 + u4)
        c = (K_C1 + K_C2 * exp(-K_C3 * u)) * e
    else
        t2 = t * t; t4 = t2 * t2
        invt = 1.0 / t
        tu = tv < 1e-2 ? one(t) : tanh(invt)
        tsu = tv < 1e-4 ? one(t) : tanh(sqrt(invt))
        a = 0.610887 * tu * (0.75 + 3.04363 * t2 - 0.09227 * t * t2 + 1.7035 * t4) / (1.0 + 8.31051 * t2 + 5.1105 * t4)
        b = tsu * (K_B1 + K_B2 * t2 + K_B3 * t4) / (1.0 + K_B4 * t2 + K_B5 * t4)
        d = tsu * (K_D1 + K_D2 * t2 + K_D3 * t4) / (1.0 + K_D4 * t2 + K_D5 * t4)
        e = tu * (K_E1 + K_E2 * t2 + K_E3 * t4) / (1.0 + K_E4 * t2 + K_E5 * t4)
        c = (K_C1 + K_C2 * (tv < 1e-3 ? zero(t) : exp(-K_C3 * invt))) * e
    end
    srs = sqrt(rs)
    return -Γe * (a + b * srs + c * rs) / (1.0 + d * srs + e * rs)
end

"Combined non-ideal Coulomb free energy density f_Coul = (F_ee + F_ii + F_ie) / V [erg cm⁻³] via linear mixing rule."
function f_coulomb(ne, ions, T)
    Z0.pval(ne) <= 1e-30 && return zero(ne * T)
    kT = kB * T
    ae = cbrt(3.0 / (4π * ne))
    rs = ae / a0
    Γe = qe * qe / (ae * kT)
    pF = ħ * cbrt(3π * π * ne)
    xr = pF / (me * c_light)
    TF_nr = (pF * pF) / (2 * me * kB)
    t = T / TF_nr
    f = ne * f_ee_ksdt(rs, t, Γe)
    for (ni, Zi) in ions
        # Note: even if ni == 0, ForwardDiff tangent may be nonzero! For ni == 0, ni * finite -> exact derivative!
        Z0.pval(ni) < 0 && continue
        Γi = Γe * Zi * cbrt(Zi * Zi)
        f += ni * (f_ii_pc00(Γi) + f_ie_pc00(Γe, rs, xr, Zi))
    end
    return kT * f
end

# ---------------------------------------------------------------------------------------------------------------------
# 3. Dense neutral fluid: C^∞ density-softened BMCSL + soft-core repulsive energy (Ross 1979 / SC91 / Potekhin 1996)
# ---------------------------------------------------------------------------------------------------------------------
mutable struct HSPar
    # Barker-Henderson T-dependent soft-sphere diameters (in Bohr radii a0), Potekhin (1996) eq. (4.18)
    d0_H::Float64;  T_H::Float64;  c1_H::Float64;  c2_H::Float64;  dmin_H::Float64
    d0_H2::Float64; T_H2::Float64; c1_H2::Float64; c2_H2::Float64; dmin_H2::Float64
    d0_He::Float64; T_He::Float64; c1_He::Float64; c2_He::Float64; dmin_He::Float64
    # Smooth BMCSL packing-fraction ceiling η_* (s(η_0) = (1 + (η_0/η_*)²)^(-1/6) -> η < η_* strictly for all ρ)
    η_star::Float64
    # Short-range soft-core repulsive wall f_rep = ε_rep * ξ0 * η_c² (with core diameters d_c in a0)
    ε_rep_eV::Float64
    dc_H2::Float64; dc_H::Float64; dc_He::Float64
    # Unbounding orbital-radius multipliers (1.0 = exact Potekhin 1996 eq. 4.21)
    α_ub_H::Float64
    α_ub_He::Float64
    α_ub_H2::Float64
end
const HSPAR = Ref(HSPar(
    4.6, 3.25e4, 4.0, 0.5, 0.9,    # H  (a0)
    5.6, 1.8e4,  1.6, 0.5, 1.3,    # H₂ (a0, Silvera-Goldman / SC91 σ ≈ 5.5 a0 at 300 K)
    4.5, 4.0e4,  2.2, 0.5, 0.9,    # He (a0, Aziz He-He core)
    0.50,                          # η_* < 1 (fluid packing ceiling)
    18.0,                          # ε_rep [eV]
    2.45, 1.85, 1.95,              # soft-core diameters d_c [a0] for H₂, H, He
    1.0, 1.0, 1.0                  # α_ub (exact Potekhin 1996)
))

@inline function _d_soft(T, d0, T0, c1, c2, dmin)
    t = T / T0
    d = d0 / (1.0 + log1p(c1 * sqrt(t)) + t / (1.0 + log1p(c2 * t)))
    a0 * (dmin + 0.5 * (d - dmin) + 0.5 * sqrt((d - dmin)^2 + (0.15 * dmin)^2))
end
@inline d_H(T)  = (p = HSPAR[]; _d_soft(T, p.d0_H,  p.T_H,  p.c1_H,  p.c2_H,  p.dmin_H))
@inline d_H2(T) = (p = HSPAR[]; _d_soft(T, p.d0_H2, p.T_H2, p.c1_H2, p.c2_H2, p.dmin_H2))
@inline d_He(T) = (p = HSPAR[]; _d_soft(T, p.d0_He, p.T_He, p.c1_He, p.c2_He, p.dmin_He))

"BMCSL helper h(x) = (x + (1-x)² ln(1-x)) / (x² (1-x)²), C^∞ at x = 0 via 7-term Taylor series for x < 1e-3."
@inline function _bmcsl_h(x)
    omx = 1.0 - x
    if Z0.pval(x) < 1e-3
        num = 1.5 - x * (1.0/3.0 + x * (1.0/12.0 + x * (1.0/30.0 + x * (1.0/60.0 + x * (1.0/105.0 + x / 168.0)))))
        return num / (omx * omx)
    end
    return (x + omx * omx * log1p(-x)) / (x * x * omx * omx)
end

"""
Exact C^∞ neutral fluid excess free energy density f_HS + f_rep [erg cm⁻³] for (nH2, nH, nHe) at T.
Well-defined and differentiable even when all n_i → 0 (returns 0 with exact zero gradient at ξ0 = 0).
"""
function f_hs_bmcsl(nH2, nH, nHe, nHp, nHep, nHepp, T)
    ξ0 = nH2 + nH + nHe
    Z0.pval(ξ0) <= 1e-30 && return zero(ξ0 * T)
    p = HSPAR[]
    d2 = d_H2(T); d1 = d_H(T); dhe = d_He(T)
    ξ1_0 = 0.5 * (nH2 * d2 + nH * d1 + nHe * dhe)
    ξ2_0 = π * (nH2 * d2 * d2 + nH * d1 * d1 + nHe * dhe * dhe)
    η0   = (π / 6.0) * (nH2 * d2^3 + nH * d1^3 + nHe * dhe^3)
    # Smooth C^∞ density compression: η = η0 / sqrt(1 + (η0/η_*)²) < η_* strictly for all ρ ∈ (0, ∞)!
    q = η0 / p.η_star
    inv_f = 1.0 / sqrt(1.0 + q * q)     # s³ = inv_f; s = cbrt(inv_f)
    s = cbrt(inv_f)
    ξ1 = s * ξ1_0
    ξ2 = s * s * ξ2_0
    η  = inv_f * η0
    Φ = -ξ0 * log1p(-η) + ξ1 * ξ2 / (1.0 - η) + (ξ2 * ξ2 * ξ2 / (36π)) * _bmcsl_h(η)
    # Soft-core repulsive energy (strictly convex ∝ n³, T-independent quantum overlap repulsion):
    # Note: ηc includes surrounding ions with 2*v_ion_H = v_H2 and v_ion_He = v_He so compression does not
    # vanish as neutrals ionize, driving smooth pressure metallization at ρ ~ 1–3 g cm⁻³!
    v2 = (π / 6.0) * (p.dc_H2 * a0)^3
    v1 = (π / 6.0) * (p.dc_H  * a0)^3
    vh = (π / 6.0) * (p.dc_He * a0)^3
    ηc = nH2 * v2 + nH * v1 + nHe * vh + nHp * (0.5 * v2) + (nHep + nHepp) * vh
    f_rep = (p.ε_rep_eV * eV) * ξ0 * ηc * ηc
    return kB * T * Φ + f_rep
end

# ---------------------------------------------------------------------------------------------------------------------
# 4. Internal partition functions with self-consistent Potekhin (1996) quantum unbounding Z_w(T, n_c)
# ---------------------------------------------------------------------------------------------------------------------
const NH_LV = 50
const NHE_LV = 35
const HI_DATA = [(χH * (1.0 - 1.0 / n^2), 4.0 * n^2, (4π / 3.0) * (a0 * n * sqrt((7n^2 + 5) / 4.0))^3) for n in 1:NH_LV]
const HEI_DATA = let
    v1 = (4π / 3.0) * (0.78 * a0)^3
    v2 = (4π / 3.0) * (a0 * 2.0 * sqrt((7*4 + 5) / 4.0))^3
    lv = Tuple{Float64,Float64,Float64}[
        (0.0, 1.0, v1),
        (19.8196 * eV, 3.0, v2), (20.6158 * eV, 1.0, v2),
        (20.9641 * eV, 9.0, v2), (21.2180 * eV, 3.0, v2),
    ]
    for n in 3:NHE_LV
        push!(lv, (χHe1 - Ry / n^2, 4.0 * n^2, (4π / 3.0) * (a0 * n * sqrt((7n^2 + 5) / 4.0))^3))
    end
    lv
end
const HEII_DATA = [(χHe2 * (1.0 - 1.0 / n^2), 2.0 * n^2, (4π / 3.0) * (0.5 * a0 * n * sqrt((7n^2 + 5) / 4.0))^3) for n in 1:NHE_LV]

const H2_LEV = let
    f = joinpath(get(ENV, "POPIII_LIT", joinpath(homedir(), "Projects", "popiii_lit")), "H2lines", "1H2__RACPPK.states")
    if isfile(f)
        rows = [split(l) for l in readlines(f) if !isempty(strip(l))]
        (E = [parse(Float64, r[2]) * hP * c_light for r in rows], g = [parse(Float64, r[3]) for r in rows])
    else
        (E = Float64[], g = Float64[])
    end
end
const V_H2_UNIT = (4π / 3.0) * (1.48 * a0)^3

function lnZ_atom(data::Vector{Tuple{Float64,Float64,Float64}}, T, nc, α3::Float64)
    kT = kB * T
    v1 = data[1][3]
    shift = -nc * α3 * v1
    s = zero(T + nc)
    @inbounds for i in eachindex(data)
        Ei, gi, vi = data[i]
        arg = -Ei / kT - nc * α3 * (vi - v1)
        Z0.pval(arg) > -700.0 && (s += gi * exp(arg))
    end
    return shift + log(s)
end

function lnZ_H2(T, nc)
    α3 = HSPAR[].α_ub_H2^3
    kT = kB * T
    ub = -nc * α3 * V_H2_UNIT
    Tv = Z0.pval(T)
    if !isempty(H2_LEV.E) && Tv >= 1500.0
        s = zero(T)
        @inbounds for i in eachindex(H2_LEV.E)
            arg = -H2_LEV.E[i] / kT
            Z0.pval(arg) > -700.0 && (s += H2_LEV.g[i] * exp(arg))
        end
        return ub + log(s)
    end
    θrot = 85.3; θvib = 6332.0
    Zr = zero(T)
    for J in 0:60
        arg = -J * (J + 1) * θrot / T
        Z0.pval(arg) > -700.0 && (Zr += (iseven(J) ? 1.0 : 3.0) * (2J + 1) * exp(arg))
    end
    Zv = 1.0 / (-expm1(-θvib / T))
    lnZ_rr = log(Zr * Zv)
    if isempty(H2_LEV.E) || Tv <= 800.0
        return ub + lnZ_rr
    end
    x = (log(T) - log(800.0)) / (log(1500.0) - log(800.0))
    w = x * x * x * (10.0 - 15.0 * x + 6.0 * x * x)
    s = zero(T)
    @inbounds for i in eachindex(H2_LEV.E)
        arg = -H2_LEV.E[i] / kT
        Z0.pval(arg) > -700.0 && (s += H2_LEV.g[i] * exp(arg))
    end
    return ub + (1.0 - w) * lnZ_rr + w * log(s)
end

@inline _smoother5(x) = (x = clamp(x, zero(x), one(x)); x^3 * (10 - 15x + 6x^2))
@inline _metal_w(T) = _smoother5((log10(T) - 5.0) / 0.6)

# ---------------------------------------------------------------------------------------------------------------------
# 5. Master Helmholtz free energy density f_tot(T, nvec, ...) [erg cm⁻³]
# ---------------------------------------------------------------------------------------------------------------------
@inline _lam3(m, T) = (hP * hP / (2π * m * kB * T))^1.5
# Guard at n > 1e-100 so nested 2nd derivatives (1/n²) never underflow Float64 (~1e-308);
# species with n < 1e-100 cm⁻³ contribute < 1e-100 to f, P, e, s (84+ digits below machine ε).
@inline _fid_term(n, m, T, gnuc::Float64) = Z0.pval(n) > 1e-100 ? n * (log(n) + log(_lam3(m, T) / gnuc) - 1.0) : zero(n * T)

const DEF_MZ = (12.011amu, 14.007amu, 15.999amu)
const DEF_ZZ = (6.0, 7.0, 8.0)
const DEF_χZ = (489.99eV, 667.05eV, 871.41eV)

"""
Non-translational, non-ideal-electron Helmholtz free energy density:
    f_nt = f_int+ub + f_Coul + f_HS + f_rep
Smooth and finite with exact partial derivatives ∂f_nt/∂n_i even at n_i = 0!
"""
function f_nontrans_no_e(T, nH2, nH, nHp, nHe, nHep, nHepp;
                         nZ = (0.0, 0.0, 0.0), zZ = DEF_ZZ, χZ = DEF_χZ,
                         coulomb::Bool = true, hs::Bool = true)
    kT = kB * T
    ne = nHp + nHep + 2.0 * nHepp
    hasZ = (Z0.pval(nZ[1]) > 0 || Z0.pval(nZ[2]) > 0 || Z0.pval(nZ[3]) > 0)
    wZ = hasZ ? _metal_w(T) : zero(T)
    f_Z_int = zero(T * nH2)
    if hasZ
        for i in 1:3
            Z0.pval(nZ[i]) > 0 || continue
            f_Z_int += nZ[i] * wZ * χZ[i]
            ne += wZ * zZ[i] * nZ[i]
        end
    end
    nc = ne + nH + nH2 + nHe + nHep
    p = HSPAR[]
    # Note: do NOT skip species with n == 0 so ∂f_int/∂n_i is ALWAYS exact even for vanished species!
    f_int = f_Z_int +
        nH2  * (E0_H2  - kT * lnZ_H2(T, nc)) +
        nH   * (E0_H   - kT * lnZ_atom(HI_DATA, T, nc, p.α_ub_H^3)) +
        nHe  * (E0_He  - kT * lnZ_atom(HEI_DATA, T, nc, p.α_ub_He^3)) +
        nHep * (E0_Hep - kT * lnZ_atom(HEII_DATA, T, nc, p.α_ub_He^3)) +
        nHepp * E0_Hepp
    f_c = if coulomb
        if hasZ
            f_coulomb(ne, ((nHp, 1.0), (nHep, 1.0), (nHepp, 2.0), (wZ * nZ[1], zZ[1]), (wZ * nZ[2], zZ[2]), (wZ * nZ[3], zZ[3])), T)
        else
            f_coulomb(ne, ((nHp, 1.0), (nHep, 1.0), (nHepp, 2.0)), T)
        end
    else
        zero(f_int)
    end
    f_h = hs ? f_hs_bmcsl(nH2, nH, nHe, nHp, nHep, nHepp, T) : zero(f_int)
    return f_int + f_c + f_h
end

"""
    f_total(T, nH2, nH, nHp, nHe, nHep, nHepp; nZ, pairs, coulomb, hs) -> f [erg cm⁻³]

Total Helmholtz free energy density of H/He/(CNO) gas at temperature T and species number densities [cm⁻³].
"""
function f_total(T, nH2, nH, nHp, nHe, nHep, nHepp;
                 nZ = (0.0, 0.0, 0.0), mZ = DEF_MZ, zZ = DEF_ZZ, χZ = DEF_χZ,
                 pairs::Bool = true, coulomb::Bool = true, hs::Bool = true)
    kT = kB * T
    f_id = kT * (
        _fid_term(nH2,   2mH, T, 1.0) +
        _fid_term(nH,    mH,  T, 1.0) +
        _fid_term(nHp,   mp,  T, 2.0) +
        _fid_term(nHe,   mHe, T, 1.0) +
        _fid_term(nHep,  mHe - me, T, 1.0) +
        _fid_term(nHepp, mHe - 2me, T, 1.0)
    )
    ne = nHp + nHep + 2.0 * nHepp
    hasZ = (Z0.pval(nZ[1]) > 0 || Z0.pval(nZ[2]) > 0 || Z0.pval(nZ[3]) > 0)
    wZ = hasZ ? _metal_w(T) : zero(T)
    if hasZ
        for i in 1:3
            Z0.pval(nZ[i]) > 0 || continue
            f_id += kT * _fid_term(nZ[i], mZ[i], T, 1.0)
            ne += wZ * zZ[i] * nZ[i]
        end
    end
    f_e = (Z0.pval(ne) > 1e-100 || (pairs && Z0.pval(T) > 5e7)) ? SE.electron_f(ne, T; pairs)[1] : zero(f_id)
    f_nt = f_nontrans_no_e(T, nH2, nH, nHp, nHe, nHep, nHepp; nZ, zZ, χZ, coulomb, hs)
    return f_id + f_e + f_nt
end

"""
    f_excess(T, nH2, nH, nHp, nHe, nHep, nHepp; ...) -> f_ex [erg cm⁻³]

Excess (non-ideal) Helmholtz free energy density = f_tot - f_ideal_Boltzmann_unperturbed.
Its partial derivatives ∂f_ex/∂n_i are the exact non-ideal chemical potential shifts Δμ_i^ex
entering ChemNet detailed balance!
"""
function f_excess(T, nH2, nH, nHp, nHe, nHep, nHepp;
                  nZ = (0.0, 0.0, 0.0), zZ = DEF_ZZ,
                  pairs::Bool = true, coulomb::Bool = true, hs::Bool = true)
    kT = kB * T
    ne = nHp + nHep + 2.0 * nHepp
    hasZ = (Z0.pval(nZ[1]) > 0 || Z0.pval(nZ[2]) > 0 || Z0.pval(nZ[3]) > 0)
    wZ = hasZ ? _metal_w(T) : zero(T)
    if hasZ
        for i in 1:3; Z0.pval(nZ[i]) > 0 && (ne += wZ * zZ[i] * nZ[i]); end
    end
    nc = ne + nH + nH2 + nHe + nHep
    p = HSPAR[]
    z0 = zero(nc)
    f_ub_shift =
        nH2  * kT * (lnZ_H2(T, z0) - lnZ_H2(T, nc)) +
        nH   * kT * (lnZ_atom(HI_DATA, T, z0, p.α_ub_H^3) - lnZ_atom(HI_DATA, T, nc, p.α_ub_H^3)) +
        nHe  * kT * (lnZ_atom(HEI_DATA, T, z0, p.α_ub_He^3) - lnZ_atom(HEI_DATA, T, nc, p.α_ub_He^3)) +
        nHep * kT * (lnZ_atom(HEII_DATA, T, z0, p.α_ub_He^3) - lnZ_atom(HEII_DATA, T, nc, p.α_ub_He^3))
    fe_full, _, _, _ = SE.electron_f(ne, T; pairs)
    nQe = 2.0 / _lam3(me, T)
    fe_boltz = Z0.pval(ne) > 1e-290 ? ne * kT * (log(ne / nQe) - 1.0) : zero(fe_full)
    f_e_ex = fe_full - fe_boltz
    f_c = if coulomb
        if hasZ
            f_coulomb(ne, ((nHp, 1.0), (nHep, 1.0), (nHepp, 2.0), (wZ * nZ[1], zZ[1]), (wZ * nZ[2], zZ[2]), (wZ * nZ[3], zZ[3])), T)
        else
            f_coulomb(ne, ((nHp, 1.0), (nHep, 1.0), (nHepp, 2.0)), T)
        end
    else
        zero(f_e_ex)
    end
    f_h = hs ? f_hs_bmcsl(nH2, nH, nHe, nHp, nHep, nHepp, T) : zero(f_e_ex)
    return f_ub_shift + f_e_ex + f_c + f_h
end

# ---------------------------------------------------------------------------------------------------------------------
# 6. Equilibrium Solver: Unclamped Chemical-Potential Newton on R(z) = Δμ/kT = 0 + Exact IFT Derivatives
# ---------------------------------------------------------------------------------------------------------------------
@inline function _lse2(a, b)
    m = max(a, b)
    m == -Inf && return m
    return m + log(exp(a - m) + exp(b - m))
end

"Unpack log-ratios z = (z1, z2, z3, z4) into species densities AND exact non-underflowing log(nH), log(ne)."
function _unpack_z_full(z, nH_tot, nHe_tot, neZ)
    T0 = typeof(z[1] + nH_tot + nHe_tot + neZ)
    # H family
    mH_ = max(zero(z[1]), z[1], z[2])
    a0 = exp(-mH_); a2 = exp(z[1] - mH_); ap = exp(z[2] - mH_)
    sH = a0 + a2 + ap
    nH  = nH_tot * (a0 / sH)
    nH2 = 0.5 * nH_tot * (a2 / sH)
    nHp = nH_tot * (ap / sH)
    ln_nH  = Z0.pval(nH_tot) > 0 ? log(nH_tot) - mH_ - log(sH) : convert(T0, -1e6)
    ln_nHp = ln_nH + z[2]
    # He family
    z34 = z[3] + z[4]
    mHe_ = max(zero(z[3]), z[3], z34)
    b0 = exp(-mHe_); b1 = exp(z[3] - mHe_); b2 = exp(z34 - mHe_)
    sHe = b0 + b1 + b2
    nHe   = nHe_tot * (b0 / sHe)
    nHep  = nHe_tot * (b1 / sHe)
    nHepp = nHe_tot * (b2 / sHe)
    ln_nHe   = Z0.pval(nHe_tot) > 0 ? log(nHe_tot) - mHe_ - log(sHe) : convert(T0, -1e6)
    ln_nHep  = ln_nHe + z[3]
    ln_nHepp = ln_nHe + z34
    # Exact log(ne) via log-sum-exp (never underflows even at T = 100 K!)
    ln_ne = convert(T0, -1e6)
    if Z0.pval(nH_tot) > 0
        ln_ne = ln_nHp
    end
    if Z0.pval(nHe_tot) > 0
        ln_ne = Z0.pval(nH_tot) > 0 ? _lse2(ln_ne, _lse2(ln_nHep, ln_nHepp + log(2.0))) : _lse2(ln_nHep, ln_nHepp + log(2.0))
    end
    if Z0.pval(neZ) > 0
        ln_ne = _lse2(ln_ne, log(neZ))
    end
    return nH2, nH, nHp, nHe, nHep, nHepp, ln_nH, ln_ne
end

"Exact stoichiometric chemical-potential balance R(z) = Δμ/(k_B T). Smooth, unclamped, O(1) Jacobian everywhere!"
function _chem_resid(z, ρ, T, X::Float64, Y::Float64, Z::Float64;
                     pairs::Bool = true, coulomb::Bool = true, hs::Bool = true)
    nH_tot = ρ * X / mH
    nHe_tot = ρ * Y / mHe
    sZ = SE.ZMIX.C + SE.ZMIX.N + SE.ZMIX.O
    nZ = Z > 0 ? (ρ * Z * SE.ZMIX.C / sZ / SE.MC, ρ * Z * SE.ZMIX.N / sZ / SE.MN, ρ * Z * SE.ZMIX.O / sZ / SE.MO) : (zero(ρ), zero(ρ), zero(ρ))
    wZ = Z > 0 ? _metal_w(T) : zero(T)
    neZ = Z > 0 ? wZ * (6nZ[1] + 7nZ[2] + 8nZ[3]) : zero(ρ * T)
    nH2, nH, nHp, nHe, nHep, nHepp, ln_nH, ln_ne = _unpack_z_full(z, nH_tot, nHe_tot, neZ)
    kT = kB * T
    # Non-translational chemical potentials ∂f_nt/∂n_i (finite and exact even at n_i = 0):
    nvec = [nH2, nH, nHp, nHe, nHep, nHepp]
    g_nt = ForwardDiff.gradient(v -> f_nontrans_no_e(T, v[1], v[2], v[3], v[4], v[5], v[6]; nZ, coulomb, hs), nvec)
    # Electron chemical potential μ_e / (k_B T) without underflow:
    if Z0.pval(ln_ne) > -400.0 || Z0.pval(T) > 5e7
        ne_use = Z0.pval(ln_ne) > -700.0 ? exp(ln_ne) : zero(ln_ne)
        _, μe, _, _ = SE.electron_f(ne_use, T; pairs)
        μe_kT = μe / kT
    else
        μe_kT = ln_ne - log(2.0 / _lam3(me, T))
    end
    # Translational thermal wavelength ratios (constant per reaction):
    lH = _lam3(mH, T); lH2 = _lam3(2mH, T); lp = _lam3(mp, T)
    lHe = _lam3(mHe, T); lHep = _lam3(mHe - me, T); lHepp = _lam3(mHe - 2me, T)
    r1 = X > 0 ? ((z[1] - log(2.0) - ln_nH) + log(lH2 / (lH * lH)) + (g_nt[1] - 2g_nt[2]) / kT) : z[1]
    r2 = X > 0 ? (z[2] - log(2.0) + log(lp / lH) + μe_kT + (g_nt[3] - g_nt[2]) / kT) : z[2]
    r3 = Y > 0 ? (z[3] + log(lHep / lHe) + μe_kT + (g_nt[5] - g_nt[4]) / kT) : z[3]
    r4 = Y > 0 ? (z[4] + log(lHepp / lHep) + μe_kT + (g_nt[6] - g_nt[5]) / kT) : z[4]
    return [r1, r2, r3, r4]
end

"Ideal-gas Saha initial guess for z = (z1, z2, z3, z4) at (ρ, T)."
function _init_z(ρ::Float64, T::Float64, X::Float64, Y::Float64, Z::Float64)
    nH_tot = ρ * X / mH; nHe_tot = ρ * Y / mHe
    kT = kB * T
    lH = _lam3(mH, T); lH2 = _lam3(2mH, T); le = _lam3(me, T)
    lnZH = lnZ_atom(HI_DATA, T, 0.0, 1.0); lnZH2 = lnZ_H2(T, 0.0)
    lnZHe = lnZ_atom(HEI_DATA, T, 0.0, 1.0); lnZHep = lnZ_atom(HEII_DATA, T, 0.0, 1.0)
    lnKH2 = log(lH2 / (lH * lH)) + 2lnZH - lnZH2 - D0_H2 / kT
    lnSH = log(4.0 / le) - lnZH - χH / kT
    lnS1 = log(2.0 / le) + lnZHep - lnZHe - χHe1 / kT
    lnS2 = log(2.0 / le) - lnZHep - χHe2 / kT
    sZ = SE.ZMIX.C + SE.ZMIX.N + SE.ZMIX.O
    neZ = Z > 0 ? _metal_w(T) * ρ * Z * (6SE.ZMIX.C/SE.MC + 7SE.ZMIX.N/SE.MN + 8SE.ZMIX.O/SE.MO) / sZ : 0.0
    nemax = max(nH_tot + 2nHe_tot + neZ, 1e-30)
    lo = log(nemax) - 1000.0; hi = log(nemax) + 2.0
    function q_diff(lne)
        ne = exp(clamp(lne, -700.0, 700.0))
        r_p = exp(clamp(lnSH - lne, -300.0, 300.0))
        b = 1.0 + r_p
        invK = exp(clamp(-lnKH2, -300.0, 300.0))
        nH_ = nH_tot > 0 ? (2 * nH_tot / (b + sqrt(b * b + 8 * invK * nH_tot))) : 0.0
        nHp_ = nH_ * r_p
        r1 = exp(clamp(lnS1 - lne, -300.0, 300.0))
        r2 = exp(clamp(lnS2 - lne, -300.0, 300.0))
        f0 = 1.0 / (1.0 + r1 * (1.0 + r2))
        nHep_ = nHe_tot * f0 * r1; nHepp_ = nHep_ * r2
        return (nHp_ + nHep_ + 2nHepp_ + neZ) - ne
    end
    if q_diff(lo) <= 0
        lne = lo
    else
        for _ in 1:100
            m = 0.5 * (lo + hi)
            q_diff(m) > 0 ? (lo = m) : (hi = m)
            hi - lo < 1e-11 && break
        end
        lne = 0.5 * (lo + hi)
    end
    r_p = exp(clamp(lnSH - lne, -300.0, 300.0))
    b = 1.0 + r_p; invK = exp(clamp(-lnKH2, -300.0, 300.0))
    nH_ = nH_tot > 0 ? (2 * nH_tot / (b + sqrt(b * b + 8 * invK * nH_tot))) : 1e-30
    # Do not artificially clamp z to ±400: at cold T = 300 K, z4 ~ -2000 is the exact analytical root!
    z1 = X > 0 ? (log(2.0) - lnKH2 + log(max(nH_, 1e-300))) : 0.0
    z2 = X > 0 ? (lnSH - lne) : 0.0
    z3 = Y > 0 ? (lnS1 - lne) : 0.0
    z4 = Y > 0 ? (lnS2 - lne) : 0.0
    return [z1, z2, z3, z4]
end

"Core Newton solver for R(z) = 0 starting from z_start. Returns (z, ok, rnorm, f)."
function _solve_newton(z_start::Vector{Float64}, ρ::Float64, T::Float64, X::Float64, Y::Float64, Z::Float64;
                       pairs::Bool = true, coulomb::Bool = true, hs::Bool = true)
    nH_tot = ρ * X / mH; nHe_tot = ρ * Y / mHe
    sZ = SE.ZMIX.C + SE.ZMIX.N + SE.ZMIX.O
    nZ = Z > 0 ? (ρ * Z * SE.ZMIX.C / sZ / SE.MC, ρ * Z * SE.ZMIX.N / sZ / SE.MN, ρ * Z * SE.ZMIX.O / sZ / SE.MO) : (0.0, 0.0, 0.0)
    wZ = Z > 0 ? _metal_w(T) : 0.0
    neZ = Z > 0 ? wZ * (6nZ[1] + 7nZ[2] + 8nZ[3]) : 0.0
    rf(zv) = _chem_resid(zv, ρ, T, X, Y, Z; pairs, coulomb, hs)
    z = copy(z_start)
    ok = false
    rnorm = Inf
    res = ForwardDiff.DiffResults.JacobianResult(z)
    for it in 1:80
        ForwardDiff.jacobian!(res, rf, z)
        R = ForwardDiff.DiffResults.value(res)
        J = ForwardDiff.DiffResults.jacobian(res)
        rnorm = maximum(abs, R)
        if rnorm < 1e-11
            ok = true; break
        end
        ϕ = 0.5 * sum(abs2, R)
        δ = try
            -(J \ R)
        catch
            fill(NaN, 4)
        end
        gϕ = J' * R
        if !all(isfinite, δ) || dot(gϕ, δ) >= 0.0
            A = J' * J
            λ = max(1e-6 * tr(A), 1e-10)
            δ = -((A + λ * I) \ gϕ)
        end
        # Allow large steps along linear tail directions (where |z_k| > 40), keep ≤ 10 in active transition zone
        for k in 1:4
            cap = abs(z[k]) > 40.0 ? 200.0 : 10.0
            abs(δ[k]) > cap && (δ[k] = copysign(cap, δ[k]))
        end
        # Backtracking line search on ϕ = 0.5 ‖R‖²
        α = 1.0
        slope = dot(gϕ, δ)
        moved = false
        for _ in 1:25
            zt = z .+ α .* δ
            Rt = rf(zt)
            if all(isfinite, Rt)
                ϕt = 0.5 * sum(abs2, Rt)
                if ϕt < ϕ + 1e-4 * α * slope || ϕt < (1.0 - 1e-8) * ϕ
                    z = zt; moved = true; break
                end
            end
            α *= 0.5
        end
        if !moved
            # Try pure LM steepest descent step with stronger damping before giving up
            A = J' * J
            for λ in (1e-3, 1e-1, 10.0)
                δ2 = -((A + λ * max(tr(A), 1.0) * I) \ gϕ)
                m2 = maximum(abs, δ2); m2 > 4.0 && (δ2 .*= 4.0 / m2)
                zt = z .+ δ2
                Rt = rf(zt)
                if all(isfinite, Rt) && 0.5 * sum(abs2, Rt) < ϕ
                    z = zt; moved = true; break
                end
            end
            !moved && break
        end
    end
    nH2, nH, nHp, nHe, nHep, nHepp, _, _ = _unpack_z_full(z, nH_tot, nHe_tot, neZ)
    fval = f_total(T, nH2, nH, nHp, nHe, nHep, nHepp; nZ, pairs, coulomb, hs)
    return z, ok, rnorm, fval
end

"""
    equilibrate(ρ, T; X = 0.75, Y = 0.25, Z = 0.0, pairs = true, coulomb = true, hs = true, z0 = nothing)

Solve exact chemical equilibrium R(z) = Δμ/kT = 0 for `f_total` at `(ρ, T)` (Float64).
At dense cool conditions (ρ > 0.2 g cm⁻³, T < 2e5 K) checks both molecular and ionized basins to guarantee
the global minimum of `f_total`.
"""
function equilibrate(ρ::Float64, T::Float64; X::Float64 = 0.75, Y::Float64 = 0.25, Z::Float64 = 0.0,
                     pairs::Bool = true, coulomb::Bool = true, hs::Bool = true, z0 = nothing)
    z_init = z0 === nothing ? _init_z(ρ, T, X, Y, Z) : copy(z0)
    z, ok, rnorm, fval = _solve_newton(z_init, ρ, T, X, Y, Z; pairs, coulomb, hs)
    # At dense cool/warm conditions (pressure-ionization zone), also test half-ionized & ionized seeds to guarantee global minimum!
    if z0 === nothing && ρ > 0.2 && T < 2e5
        for z_seed in (
            [X > 0 ? -1.0 : 0.0, X > 0 ? 0.5 : 0.0, Y > 0 ? -5.0 : 0.0, Y > 0 ? -10.0 : 0.0],
            [X > 0 ? -15.0 : 0.0, X > 0 ? 8.0 : 0.0, Y > 0 ? 4.0 : 0.0, Y > 0 ? 2.0 : 0.0],
        )
            zi, oki, ri, fi = _solve_newton(z_seed, ρ, T, X, Y, Z; pairs, coulomb, hs)
            if oki && (!ok || fi < fval)
                z, ok, rnorm, fval = zi, oki, ri, fi
            end
        end
    end
    # If first attempt failed, try a neutral-atomic seed too:
    if !ok
        z_alt = [X > 0 ? -5.0 : 0.0, X > 0 ? -2.0 : 0.0, Y > 0 ? -5.0 : 0.0, Y > 0 ? -10.0 : 0.0]
        za, oka, ra, fa = _solve_newton(z_alt, ρ, T, X, Y, Z; pairs, coulomb, hs)
        if oka && (!ok || fa < fval)
            z, ok, rnorm, fval = za, oka, ra, fa
        end
    end
    nH_tot = ρ * X / mH; nHe_tot = ρ * Y / mHe
    sZ = SE.ZMIX.C + SE.ZMIX.N + SE.ZMIX.O
    nZ = Z > 0 ? (ρ * Z * SE.ZMIX.C / sZ / SE.MC, ρ * Z * SE.ZMIX.N / sZ / SE.MN, ρ * Z * SE.ZMIX.O / sZ / SE.MO) : (0.0, 0.0, 0.0)
    wZ = Z > 0 ? _metal_w(T) : 0.0
    neZ = Z > 0 ? wZ * (6nZ[1] + 7nZ[2] + 8nZ[3]) : 0.0
    nH2, nH, nHp, nHe, nHep, nHepp, _, _ = _unpack_z_full(z, nH_tot, nHe_tot, neZ)
    ne = nHp + nHep + 2nHepp + neZ
    return (; nH2, nH, nHp, nHe, nHep, nHepp, ne, nZ, z, ok, rnorm, f = fval)
end

"Exact Euler (p_gas, e_gas, s_gas) from f_total at arbitrary (ρ, T, z) — generic in Dual numbers for IFT!"
function _pes_from_z(ρ, T, z, X::Float64, Y::Float64, Z::Float64;
                     pairs::Bool = true, coulomb::Bool = true, hs::Bool = true, radiation::Bool = false)
    nH_tot = ρ * X / mH; nHe_tot = ρ * Y / mHe
    sZ = SE.ZMIX.C + SE.ZMIX.N + SE.ZMIX.O
    nZ = Z > 0 ? (ρ * Z * SE.ZMIX.C / sZ / SE.MC, ρ * Z * SE.ZMIX.N / sZ / SE.MN, ρ * Z * SE.ZMIX.O / sZ / SE.MO) : (zero(ρ), zero(ρ), zero(ρ))
    wZ = Z > 0 ? _metal_w(T) : zero(T)
    neZ = Z > 0 ? wZ * (6nZ[1] + 7nZ[2] + 8nZ[3]) : zero(ρ * T)
    nH2, nH, nHp, nHe, nHep, nHepp, _, _ = _unpack_z_full(z, nH_tot, nHe_tot, neZ)
    # Exact Euler pressure via V-scaling identity: f(λ n) = F(V/λ)/V => P = λ ∂f(λ n)/∂λ|_{λ=1} - f!
    # Using single scalar scaling λ instead of 9D gradient is 6× faster AND handles n_i=0 & nZ automatically!
    f0 = f_total(T, nH2, nH, nHp, nHe, nHep, nHepp; nZ, pairs, coulomb, hs)
    df_dλ = ForwardDiff.derivative(λ -> f_total(T, λ*nH2, λ*nH, λ*nHp, λ*nHe, λ*nHep, λ*nHepp;
                                                nZ = (λ*nZ[1], λ*nZ[2], λ*nZ[3]), pairs, coulomb, hs), 1.0)
    p_gas = df_dλ - f0
    df_dT = ForwardDiff.derivative(t -> f_total(t, nH2, nH, nHp, nHe, nHep, nHepp; nZ, pairs, coulomb, hs), T)
    s_gas = -df_dT / ρ
    e_gas = (f0 - T * df_dT) / ρ
    if radiation
        return p_gas + SE.Prad(T), e_gas + SE.erad(ρ, T), s_gas + SE.srad(ρ, T)
    end
    return p_gas, e_gas, s_gas
end

"""
    eos_lte(ρ, T; X = 0.75, Y = 0.25, Z = 0.0, pairs = true, coulomb = true, hs = true, radiation = false)

Full thermodynamically consistent LTE EOS at `(ρ, T, X, Y, Z)` via `f_total`.
Second derivatives (`χρ, χT, cV, dedρ, Γ1, ∇ad`) are computed analytically to machine precision
via the Implicit Function Theorem on `R(z*, ρ, T) = 0` + nested ForwardDiff (zero finite differences!).
"""
function eos_lte(ρ::Float64, T::Float64; X::Float64 = 0.75, Y::Float64 = 0.25, Z::Float64 = 0.0,
                 pairs::Bool = true, coulomb::Bool = true, hs::Bool = true, radiation::Bool = false)
    eq = equilibrate(ρ, T; X, Y, Z, pairs, coulomb, hs)
    z = eq.z
    # Implicit Function Theorem: J * dz = -(∂R/∂lnρ) dlnρ - (∂R/∂lnT) dlnT
    rf_z(zv) = _chem_resid(zv, ρ, T, X, Y, Z; pairs, coulomb, hs)
    J = ForwardDiff.jacobian(rf_z, z)
    dR_dlr = ForwardDiff.derivative(lr -> _chem_resid(z, exp(lr), T, X, Y, Z; pairs, coulomb, hs), log(ρ))
    dR_dlt = ForwardDiff.derivative(lt -> _chem_resid(z, ρ, exp(lt), X, Y, Z; pairs, coulomb, hs), log(T))
    dz_dlr = -(J \ dR_dlr)
    dz_dlt = -(J \ dR_dlt)
    # Propagate 2-component Dual for (lnρ, lnT) through _pes_from_z in a SINGLE exact pass!
    Ttag = typeof(ForwardDiff.Tag(_pes_from_z, Float64))
    lr_d = ForwardDiff.Dual{Ttag}(log(ρ), 1.0, 0.0)
    lt_d = ForwardDiff.Dual{Ttag}(log(T), 0.0, 1.0)
    z_d = [ForwardDiff.Dual{Ttag}(z[i], dz_dlr[i], dz_dlt[i]) for i in 1:4]
    Pd, ed, sd = _pes_from_z(exp(lr_d), exp(lt_d), z_d, X, Y, Z; pairs, coulomb, hs, radiation)
    P_ = ForwardDiff.value(Pd); e_ = ForwardDiff.value(ed); s_ = ForwardDiff.value(sd)
    dP_dlr, dP_dlt = ForwardDiff.partials(Pd).values
    de_dlr, de_dlt = ForwardDiff.partials(ed).values
    ds_dlr, ds_dlt = ForwardDiff.partials(sd).values
    χρ = dP_dlr / P_; χT = dP_dlt / P_
    cV = de_dlt / T
    dedρ = de_dlr / ρ
    Sρ = ds_dlr; ST = ds_dlt
    ∇ad = 1.0 / (χT - χρ * ST / Sρ)
    Γ1 = χρ - χT * Sρ / ST
    δ = χT / χρ
    cP = ST - Sρ * δ
    nH_tot = ρ * X / mH; nHe_tot = ρ * Y / mHe
    xHp = nH_tot > 0 ? eq.nHp / nH_tot : 0.0
    fH2 = nH_tot > 0 ? 2 * eq.nH2 / nH_tot : 0.0
    xHep = nHe_tot > 0 ? eq.nHep / nHe_tot : 0.0
    xHepp = nHe_tot > 0 ? eq.nHepp / nHe_tot : 0.0
    return (; P = P_, e = e_, s = s_, cV, cVs = ST, dedρ, χρ, χT, Γ1, ∇ad, δ, cP,
              cs = sqrt(max(Γ1 * P_ / ρ, 0.0)), ne = eq.ne, xHp, fH2, xHep, xHepp, ok = eq.ok)
end

# ---------------------------------------------------------------------------------------------------------------------
# 7. Exact Non-Ideal Chemical Potential Shifts Δμ_i^ex / (k_B T) for ChemNet Detailed Balance
# ---------------------------------------------------------------------------------------------------------------------
# Match ChemNet's ideal reference (n ≤ 15 electronic partition functions, χH = 13.598 eV) so that
# ln K_ideal(ChemNet) - Σ ν_i Δμ_i^ex/kT equals the exact FreeEnergyEOS equilibrium constant to 15 digits!
const _CN_NH = 15
const _CN_χH = 13.598 * eV
const _CN_E0_H  = -_CN_χH
const _CN_E0_H2 = -(2_CN_χH + D0_H2)
const _CN_HI_DATA = HI_DATA[1:_CN_NH]
const _CN_HEI_DATA = HEI_DATA[1:min(_CN_NH + 3, length(HEI_DATA))]
const _CN_HEII_DATA = HEII_DATA[1:_CN_NH]

@inline function _lnZ_ideal_cn(data::Vector{Tuple{Float64,Float64,Float64}}, T)
    kT = kB * T; s = zero(T)
    @inbounds for i in eachindex(data)
        arg = -data[i][1] / kT
        Z0.pval(arg) > -700.0 && (s += data[i][2] * exp(arg))
    end
    return log(s)
end

"""
Excess free energy density f_ex(T, nH2, nH, nHp, nHe, nHep, nHepp, ne) with `ne` as an independent 7th variable,
relative to ChemNet's ideal Boltzmann gas (n ≤ 15, unperturbed partition functions).
Its 7 partial derivatives divided by k_B T are the exact dimensionless chemical-potential shifts
(Δμ_H2, Δμ_H, Δμ_Hp, Δμ_He, Δμ_Hep, Δμ_Hepp, Δμ_e) / (k_B T) entering ChemNet detailed balance!
"""
function f_excess_7(T, nH2, nH, nHp, nHe, nHep, nHepp, ne;
                    nZ = (0.0, 0.0, 0.0), zZ = DEF_ZZ,
                    pairs::Bool = true, coulomb::Bool = true, hs::Bool = true)
    kT = kB * T
    hasZ = (Z0.pval(nZ[1]) > 0 || Z0.pval(nZ[2]) > 0 || Z0.pval(nZ[3]) > 0)
    wZ = hasZ ? _metal_w(T) : zero(T)
    nc = ne + nH + nH2 + nHe + nHep
    p = HSPAR[]
    z0 = zero(nc)
    # Unbounding + level-truncation + E0-reference shift relative to ChemNet's ideal gas:
    f_ub_shift =
        nH2  * ((E0_H2 - _CN_E0_H2) + kT * (lnZ_H2(T, z0) - lnZ_H2(T, nc))) +
        nH   * ((E0_H  - _CN_E0_H)  + kT * (_lnZ_ideal_cn(_CN_HI_DATA, T) - lnZ_atom(HI_DATA, T, nc, p.α_ub_H^3))) +
        nHe  * (kT * (_lnZ_ideal_cn(_CN_HEI_DATA, T)  - lnZ_atom(HEI_DATA, T, nc, p.α_ub_He^3))) +
        nHep * (kT * (_lnZ_ideal_cn(_CN_HEII_DATA, T) - lnZ_atom(HEII_DATA, T, nc, p.α_ub_He^3)))
    # Electron degeneracy + relativity + pairs relative to ideal non-rel Boltzmann electrons:
    f_e_ex = if Z0.pval(ne) > 1e-100 || (pairs && Z0.pval(T) > 5e7)
        fe_full, _, _, _ = SE.electron_f(ne, T; pairs)
        nQe = 2.0 / _lam3(me, T)
        fe_boltz = Z0.pval(ne) > 1e-100 ? ne * kT * (log(ne) - log(nQe) - 1.0) : zero(fe_full)
        fe_full - fe_boltz
    else
        zero(f_ub_shift)
    end
    # Coulomb + exchange-correlation + ion-electron screening:
    f_c = if coulomb
        if hasZ
            f_coulomb(ne, ((nHp, 1.0), (nHep, 1.0), (nHepp, 2.0), (wZ * nZ[1], zZ[1]), (wZ * nZ[2], zZ[2]), (wZ * nZ[3], zZ[3])), T)
        else
            f_coulomb(ne, ((nHp, 1.0), (nHep, 1.0), (nHepp, 2.0)), T)
        end
    else
        zero(f_ub_shift)
    end
    # Neutral fluid excluded volume + soft-core repulsion:
    f_h = hs ? f_hs_bmcsl(nH2, nH, nHe, nHp, nHep, nHepp, T) : zero(f_ub_shift)
    return f_ub_shift + f_e_ex + f_c + f_h
end

"""
    delta_mu_ex_7(T, nH2, nH, nHp, nHe, nHep, nHepp, ne; ...) -> NTuple{7}

Dimensionless non-ideal chemical potential shifts `(Δμ_H2, Δμ_H, Δμ_Hp, Δμ_He, Δμ_Hep, Δμ_Hepp, Δμ_e) / (k_B T)`
for ChemNet detailed balance (`ln K_eff = ln K_ideal - Σ ν_i Δμ_i^ex/kT`). Generic in Dual numbers!
"""
function delta_mu_ex_7(T, nH2, nH, nHp, nHe, nHep, nHepp, ne;
                       nZ = (0.0, 0.0, 0.0), zZ = DEF_ZZ,
                       pairs::Bool = true, coulomb::Bool = true, hs::Bool = true)
    # Promote cleanly so ForwardDiff gradient works seamlessly even when inputs are already Duals:
    ET = promote_type(typeof(T), typeof(nH2), typeof(nH), typeof(nHp), typeof(nHe), typeof(nHep), typeof(nHepp), typeof(ne), Float64)
    v0 = ET[nH2, nH, nHp, nHe, nHep, nHepp, ne]
    g = ForwardDiff.gradient(v -> f_excess_7(T, v[1], v[2], v[3], v[4], v[5], v[6], v[7]; nZ, zZ, pairs, coulomb, hs), v0)
    ikT = 1.0 / (kB * T)
    return (g[1] * ikT, g[2] * ikT, g[3] * ikT, g[4] * ikT, g[5] * ikT, g[6] * ikT, g[7] * ikT)
end

# ---------------------------------------------------------------------------------------------------------------------
# 8. Fast, Exact Non-LTE & Fully-Ionised EOS Directly from f_total (Nested-ForwardDiff Safe!)
# ---------------------------------------------------------------------------------------------------------------------
"Exact internal + unbounding (e_int_per_particle, p_ub_factor = -nc ∂lnZ/∂nc) for atomic level list."
function _atom_ep(data::Vector{Tuple{Float64,Float64,Float64}}, T, nc, α3::Float64)
    kT = kB * T; v1 = data[1][3]
    z = zero(T + nc); e = zero(T + nc); dv = zero(T + nc)
    @inbounds for i in eachindex(data)
        Ei, gi, vi = data[i]
        arg = -Ei / kT - nc * α3 * (vi - v1)
        if Z0.pval(arg) > -700.0
            t = gi * exp(arg)
            z += t
            e += t * Ei
            dv += t * (vi - v1)
        end
    end
    # E_int = ⟨E_i⟩;  -∂(ln Z_w)/∂nc = α3 * (v1 + ⟨vi - v1⟩) => p_ub_per_particle = kT * nc * α3 * (v1 + dv/z)
    return e / z, nc * α3 * (v1 + dv / z)
end

"Exact internal + unbounding (e_int_per_H2, p_ub_factor) from lnZ_H2(T, nc)."
function _h2_ep(T, nc)
    α3 = HSPAR[].α_ub_H2^3
    kT = kB * T
    # Exact e_int = kT * T * ∂(lnZ_H2(T, 0))/∂T via direct analytic sum (matches lnZ_H2 to 16 digits!)
    Tv = Z0.pval(T)
    e_rr = begin
        θrot = 85.3; θvib = 6332.0
        zr = zero(T); er = zero(T)
        for J in 0:60
            ej = J * (J + 1) * θrot
            arg = -ej / T
            if Z0.pval(arg) > -700.0
                t = (iseven(J) ? 1.0 : 3.0) * (2J + 1) * exp(arg)
                zr += t; er += t * ej
            end
        end
        xv = θvib / T
        ev = Z0.pval(xv) < 700.0 ? θvib / expm1(xv) : zero(T)
        kB * (er / zr + ev)
    end
    pub = nc * α3 * V_H2_UNIT
    if isempty(H2_LEV.E) || Tv <= 800.0
        return e_rr, pub
    end
    zp = zero(T); ep = zero(T)
    @inbounds for i in eachindex(H2_LEV.E)
        arg = -H2_LEV.E[i] / kT
        if Z0.pval(arg) > -700.0
            t = H2_LEV.g[i] * exp(arg)
            zp += t; ep += t * H2_LEV.E[i]
        end
    end
    e_pf = ep / zp
    Tv >= 1500.0 && return e_pf, pub
    # Inside [800, 1500] K blend: exact T * ∂lnZ/∂T = (1-w) e_rr + w e_pf + kT * (T dw/dT) * (ln Z_pf - ln Z_rr)!
    θrot = 85.3; θvib = 6332.0
    zr2 = zero(T)
    for J in 0:60
        arg = -J * (J + 1) * θrot / T
        Z0.pval(arg) > -700.0 && (zr2 += (iseven(J) ? 1.0 : 3.0) * (2J + 1) * exp(arg))
    end
    lnZ_rr = log(zr2 / (-expm1(-θvib / T)))
    lnZ_pf = log(zp)
    L = log(1500.0) - log(800.0)
    x = (log(T) - log(800.0)) / L
    w = x * x * x * (10.0 - 15.0 * x + 6.0 * x * x)
    T_dwdT = (30.0 * x * x * (1.0 - x) * (1.0 - x)) / L
    return (1.0 - w) * e_rr + w * e_pf + kT * T_dwdT * (lnZ_pf - lnZ_rr), pub
end

"""
    eos_nonlte_dens(ρ, T, nH2, nH, nHp, nHe, nHep, nHepp; nHm=0, nH2p=0, nLi=0, nLip=0, nZ=(0,0,0), pairs=true, coulomb=true, hs=true) -> (u, p)

Exact non-LTE (u [erg/g], p [dyn/cm²]) derived analytically from `f_total` (matches `_pes_from_z` to 15 digits
at equilibrium, 10–20× faster than end-to-end AD, and nested-ForwardDiff safe for Newton + TDC LResp!).
"""
function eos_nonlte_dens(ρ, T, nH2, nH, nHp, nHe, nHep, nHepp;
                         nHm = zero(ρ), nH2p = zero(ρ), nLi = zero(ρ), nLip = zero(ρ),
                         nZ = (zero(ρ), zero(ρ), zero(ρ)), zZ = DEF_ZZ, χZ = DEF_χZ,
                         pairs::Bool = true, coulomb::Bool = true, hs::Bool = true)
    kT = kB * T
    # 1. Heavy ideal translational pressure & energy + E0 + trace species (H⁻, H₂⁺, Li, Li⁺)
    n_heavy = nH2 + nH + nHp + nHe + nHep + nHepp + nHm + nH2p + nLi + nLip
    E_heavy = (1.5 * (n_heavy - nH2p) + 2.5 * nH2p) * kT +
              nH2 * E0_H2 + nH * E0_H + nHp * E0_Hp + nHm * E0_Hm + nH2p * E0_H2p +
              nHe * E0_He + nHep * E0_Hep + nHepp * E0_Hepp - nLi * (5.3917 * eV)
    ne = nHp + nHep + 2.0 * nHepp + nH2p + nLip - nHm
    hasZ = (Z0.pval(nZ[1]) > 0 || Z0.pval(nZ[2]) > 0 || Z0.pval(nZ[3]) > 0)
    wZ = hasZ ? _metal_w(T) : zero(T)
    if hasZ
        # Exact T-derivative of nZ * wZ(T) * χZ contributes -T ∂(wZ)/∂T to internal energy E = f - T ∂f/∂T!
        dwZ_dlnT = if 5.0 < Z0.pval(log10(T)) < 5.6
            x = (log10(T) - 5.0) / 0.6
            (30.0 * x * x * (1.0 - x) * (1.0 - x)) / (0.6 * log(10.0))
        else
            zero(T)
        end
        for i in 1:3
            Z0.pval(nZ[i]) > 0 || continue
            n_heavy += nZ[i]
            E_heavy += nZ[i] * (1.5 * kT + (wZ - dwZ_dlnT) * χZ[i])
            ne += wZ * zZ[i] * nZ[i]
        end
    end
    ne_pos = Z0.pval(ne) > 0 ? ne : zero(ne)
    # 2. Internal excitation + Potekhin quantum-unbounding pressure (exact derivatives of -n_s kT ln Z_w(T, nc)):
    nc = ne_pos + nH + nH2 + nHe + nHep
    hp = HSPAR[]
    p_ub = zero(n_heavy * kT)
    if Z0.pval(nH2) > 0
        e2, u2 = _h2_ep(T, nc)
        E_heavy += nH2 * e2; p_ub += nH2 * u2 * kT
    end
    if Z0.pval(nH) > 0
        e1, u1 = _atom_ep(HI_DATA, T, nc, hp.α_ub_H^3)
        E_heavy += nH * e1; p_ub += nH * u1 * kT
    end
    if Z0.pval(nHe) > 0
        eh0, uh0 = _atom_ep(HEI_DATA, T, nc, hp.α_ub_He^3)
        E_heavy += nHe * eh0; p_ub += nHe * uh0 * kT
    end
    if Z0.pval(nHep) > 0
        eh1, uh1 = _atom_ep(HEII_DATA, T, nc, hp.α_ub_He^3)
        E_heavy += nHep * eh1; p_ub += nHep * uh1 * kT
    end
    p = n_heavy * kT + p_ub
    E = E_heavy
    # 3. Ideal Fermi–Dirac relativistic e± gas (StellarEOS.electron_pe — exact match to electron_f!)
    Pe, Ee = SE.electron_pe(ne_pos, T; pairs)
    p += Pe; E += Ee
    # 4. Non-ideal Coulomb (PC00+KSDT14) + Neutral BMCSL/soft-core via exact 2-scalar Euler AD (skipped only when truly negligible!)
    if (coulomb && Z0.pval(ne_pos) > 1e6) || (hs && Z0.pval(nH2 + nH + nHe) > 1e12)
        f_ch(λ, t) = begin
            fc = if coulomb && Z0.pval(ne_pos) > 1e6
                if hasZ
                    f_coulomb(λ * ne_pos, ((λ * nHp, 1.0), (λ * nHep, 1.0), (λ * nHepp, 2.0),
                              (λ * wZ * nZ[1], zZ[1]), (λ * wZ * nZ[2], zZ[2]), (λ * wZ * nZ[3], zZ[3])), t)
                else
                    f_coulomb(λ * ne_pos, ((λ * nHp, 1.0), (λ * nHep, 1.0), (λ * nHepp, 2.0)), t)
                end
            else
                zero(λ * t * ne_pos)
            end
            fh = (hs && Z0.pval(nH2 + nH + nHe) > 1e12) ?
                f_hs_bmcsl(λ * nH2, λ * nH, λ * nHe, λ * nHp, λ * nHep, λ * nHepp, t) : zero(fc)
            fc + fh
        end
        f0 = f_ch(1.0, T)
        df_dλ = ForwardDiff.derivative(λ -> f_ch(λ, T), 1.0)
        df_dT = ForwardDiff.derivative(t -> f_ch(1.0, t), T)
        p += df_dλ - f0
        E += f0 - T * df_dT
    end
    return E / ρ, p
end

"Exact fully-ionised limit of `f_total` for arbitrary NucNet composition `comp` (H⁺, He²⁺, bare CNO, Fermi–Dirac e±, PC00+KSDT14 Coulomb)."
function eos_full_ion(ρ, T, comp; pairs::Bool = true, coulomb::Bool = true)
    nHp   = ρ * comp.x[1] / comp.m[1]
    nHepp = ρ * comp.x[2] / comp.m[2]
    nZ = (comp.x[3] > 0 ? ρ * comp.x[3] / comp.m[3] : zero(ρ),
          comp.x[4] > 0 ? ρ * comp.x[4] / comp.m[4] : zero(ρ),
          comp.x[5] > 0 ? ρ * comp.x[5] / comp.m[5] : zero(ρ))
    z = zero(ρ)
    return eos_nonlte_dens(ρ, T, z, z, nHp, z, z, nHepp; nZ, pairs, coulomb, hs = false)
end

# ---------------------------------------------------------------------------------------------------------------------
# 9. Self-Generated C² Bicubic B-Spline LTE Table (Built Directly from Our Own FreeEnergyEOS.eos_lte!)
# ---------------------------------------------------------------------------------------------------------------------
# Exact energy zero offset: cold H₂ + neutral He at 0 K (strictly positive u + UOFF_FE everywhere!)
const UOFF_FE = (0.75 / mH) * (χH + 0.5 * D0_H2)

"1D uniform cubic B-spline prefilter with natural boundary conditions (C² continuous value, 1st & 2nd derivatives!)."
function _bspline1d_prefilter(y::AbstractVector{Float64})
    n = length(y)
    # Solve tridiagonal system for interior c[2..n-1] with c[1] = y[1], c[n] = y[n], ghost c[0]=2c[1]-c[2], c[n+1]=2c[n]-c[n-1]
    # Interior equation: (1/6) c[i-1] + (4/6) c[i] + (1/6) c[i+1] = y[i]  =>  c[i-1] + 4 c[i] + c[i+1] = 6 y[i]
    c = zeros(Float64, n + 2)     # 1-based shifted: c[i+1] corresponds to knot i (i = 0..n+1)
    if n == 1
        c .= y[1]; return c
    elseif n == 2
        c[2] = y[1]; c[3] = y[2]; c[1] = 2c[2] - c[3]; c[4] = 2c[3] - c[2]; return c
    end
    m = n - 2
    dl = ones(Float64, m - 1); d = fill(4.0, m); du = ones(Float64, m - 1)
    rhs = [6.0 * y[i] for i in 2:n-1]
    rhs[1] -= y[1]; rhs[end] -= y[n]
    sol = Tridiagonal(dl, d, du) \ rhs
    c[2] = y[1]; c[n+1] = y[n]
    for i in 1:m; c[i+2] = sol[i]; end
    c[1] = 2c[2] - c[3]
    c[n+2] = 2c[n+1] - c[n]
    return c
end

"2D uniform cubic B-spline prefilter: returns (n+2) × (m+2) coefficient matrix C (C² everywhere!)."
function _bspline2d_prefilter(A::Matrix{Float64})
    n, m = size(A)
    tmp = zeros(Float64, n + 2, m)
    for j in 1:m
        tmp[:, j] .= _bspline1d_prefilter(view(A, :, j))
    end
    C = zeros(Float64, n + 2, m + 2)
    for i in 1:n+2
        C[i, :] .= _bspline1d_prefilter(view(tmp, i, :))
    end
    return C
end

@inline function _bs3(c0, c1, c2, c3, t)
    om = 1.0 - t; t2 = t * t; t3 = t2 * t
    return ((om * om * om) * c0 + (3t3 - 6t2 + 4.0) * c1 + (-3t3 + 3t2 + 3t + 1.0) * c2 + t3 * c3) / 6.0
end

"Evaluate 2D C² cubic B-spline from prefiltered (n+2)×(m+2) matrix C at grid coords (x, y) ∈ [1, n] × [1, m]."
@inline function bspline2d(C::Matrix{Float64}, x, y)
    n = size(C, 1) - 2; m = size(C, 2) - 2
    xc = clamp(x, 1.0, Float64(n)); yc = clamp(y, 1.0, Float64(m))
    i = clamp(floor(Int, Z0.pval(xc)), 1, n - 1); tx = xc - i
    j = clamp(floor(Int, Z0.pval(yc)), 1, m - 1); ty = yc - j
    # Knot i corresponds to C[i+1, :]; stencil uses i-1..i+2 -> C[i .. i+3]
    @inbounds begin
        r0 = _bs3(C[i,   j], C[i,   j+1], C[i,   j+2], C[i,   j+3], ty)
        r1 = _bs3(C[i+1, j], C[i+1, j+1], C[i+1, j+2], C[i+1, j+3], ty)
        r2 = _bs3(C[i+2, j], C[i+2, j+1], C[i+2, j+2], C[i+2, j+3], ty)
        r3 = _bs3(C[i+3, j], C[i+3, j+1], C[i+3, j+2], C[i+3, j+3], ty)
    end
    return _bs3(r0, r1, r2, r3, tx)
end

struct LTETable
    X::Float64; Y::Float64; Z::Float64
    lρ0::Float64; dlρ::Float64; nρ::Int
    lT0::Float64; dlT::Float64; nT::Int
    # Raw nodal arrays AND C² B-spline coefficient arrays:
    lp::Matrix{Float64};  lus::Matrix{Float64}
    fH2::Matrix{Float64}; xHp::Matrix{Float64}; zHe1::Matrix{Float64}; zHe2::Matrix{Float64}
    c_lp::Matrix{Float64}; c_lus::Matrix{Float64}
    c_fH2::Matrix{Float64}; c_xHp::Matrix{Float64}; c_zHe1::Matrix{Float64}; c_zHe2::Matrix{Float64}
end

"""
    build_lte_table(; X=0.75, Y=0.25, Z=0.0, lρ0=-20.0, lρ1=1.5, dlρ=0.05, lT0=1.7, lT1=6.2, dlT=0.02)

Build self-consistent C² B-spline LTE table directly from `FreeEnergyEOS.equilibrate` & `_pes_from_z`
(warm-started along isotherms, multithreaded).
"""
function build_lte_table(; X::Float64 = 0.75, Y::Float64 = 0.25, Z::Float64 = 0.0,
                         lρ0::Float64 = -20.0, lρ1::Float64 = 1.5, dlρ::Float64 = 0.05,
                         lT0::Float64 = 1.7, lT1::Float64 = 6.2, dlT::Float64 = 0.02)
    nρ = round(Int, (lρ1 - lρ0) / dlρ) + 1
    nT = round(Int, (lT1 - lT0) / dlT) + 1
    lp   = zeros(nρ, nT); lus  = zeros(nρ, nT)
    fH2  = zeros(nρ, nT); xHp  = zeros(nρ, nT)
    zHe1 = zeros(nρ, nT); zHe2 = zeros(nρ, nT)
    Threads.@threads for j in 1:nT
        T = 10.0^(lT0 + (j - 1) * dlT)
        z_prev = nothing
        for i in 1:nρ
            ρ = 10.0^(lρ0 + (i - 1) * dlρ)
            # Warm-start from previous density along isotherm (except across pressure-ionization threshold where global check runs)
            eq = (ρ > 0.2 && T < 2e5) ? equilibrate(ρ, T; X, Y, Z) : equilibrate(ρ, T; X, Y, Z, z0 = z_prev)
            if !eq.ok
                eq = equilibrate(ρ, T; X, Y, Z)
            end
            z_prev = eq.ok ? eq.z : nothing
            p_, e_, _ = _pes_from_z(ρ, T, eq.z, X, Y, Z)
            nH_tot = ρ * X / mH
            lp[i, j]   = log10(p_)
            lus[i, j]  = log10(max(e_ + UOFF_FE, 1.0))
            fH2[i, j]  = nH_tot > 0 ? clamp(2 * eq.nH2 / nH_tot, 0.0, 1.0) : 0.0
            xHp[i, j]  = nH_tot > 0 ? clamp(eq.nHp / nH_tot, 0.0, 1.0) : 0.0
            zHe1[i, j] = nH_tot > 0 ? max(eq.nHep / nH_tot, 0.0) : 0.0
            zHe2[i, j] = nH_tot > 0 ? max(eq.nHepp / nH_tot, 0.0) : 0.0
        end
    end
    return LTETable(X, Y, Z, lρ0, dlρ, nρ, lT0, dlT, nT,
                    lp, lus, fH2, xHp, zHe1, zHe2,
                    _bspline2d_prefilter(lp), _bspline2d_prefilter(lus),
                    _bspline2d_prefilter(fH2), _bspline2d_prefilter(xHp),
                    _bspline2d_prefilter(zHe1), _bspline2d_prefilter(zHe2))
end

const LTE_TAB = Ref{Any}(nothing)
const LTE_TAB_LOCK = ReentrantLock()
lte_cache_path() = get(ENV, "FE_EOS_CACHE", joinpath(@__DIR__, "..", "lib", "popiii", "fe_lte_table_X75_Y25_v1.jls"))

function lte_table()
    t = LTE_TAB[]; t === nothing || return t::LTETable
    lock(LTE_TAB_LOCK) do
        t2 = LTE_TAB[]; t2 === nothing || return t2::LTETable
        path = lte_cache_path()
        tb = if isfile(path)
            raw = deserialize(path)
            raw isa LTETable ? raw : LTETable(
                raw.X, raw.Y, raw.Z, raw.lρ0, raw.dlρ, raw.nρ, raw.lT0, raw.dlT, raw.nT,
                raw.lp, raw.lus, raw.fH2, raw.xHp, raw.zHe1, raw.zHe2,
                raw.c_lp, raw.c_lus, raw.c_fH2, raw.c_xHp, raw.c_zHe1, raw.c_zHe2
            )
        else
            tb_ = build_lte_table()
            try
                mkpath(dirname(path))
                # Serialize plain NamedTuple (module/package-UUID independent!)
                serialize(path, (; tb_.X, tb_.Y, tb_.Z, tb_.lρ0, tb_.dlρ, tb_.nρ, tb_.lT0, tb_.dlT, tb_.nT,
                                   tb_.lp, tb_.lus, tb_.fH2, tb_.xHp, tb_.zHe1, tb_.zHe2,
                                   tb_.c_lp, tb_.c_lus, tb_.c_fH2, tb_.c_xHp, tb_.c_zHe1, tb_.c_zHe2))
            catch
            end
            tb_
        end
        LTE_TAB[] = tb
        tb
    end
end

@inline _gx(t::LTETable, lρ) = (lρ - t.lρ0) / t.dlρ + 1.0
@inline _gy(t::LTETable, lT) = (lT - t.lT0) / t.dlT + 1.0

"C² B-spline LTE table EOS (u [erg/g], p [dyn/cm²]) matching FreeEnergyEOS.eos_lte; nested-ForwardDiff safe!"
@inline function lte_eos(ρ, T)
    t = lte_table()
    x = _gx(t, log10(ρ)); y = _gy(t, log10(T))
    return (10.0^bspline2d(t.c_lus, x, y) - UOFF_FE, 10.0^bspline2d(t.c_lp, x, y))
end

"C² B-spline LTE table H composition (xHp, 0, xH2, 0) per H nucleus (matches PopIIIEqEOS.composition interface)."
@inline function lte_composition(ρ, T)
    t = lte_table()
    x = _gx(t, log10(ρ)); y = _gy(t, log10(T))
    xp = clamp(bspline2d(t.c_xHp, x, y), zero(ρ * T), one(ρ * T))
    f2 = clamp(bspline2d(t.c_fH2, x, y), zero(ρ * T), one(ρ * T))
    return (xp, zero(xp), 0.5 * f2, zero(xp))
end

"C² B-spline LTE table He composition (n_Hep/n_H, n_Hepp/n_H) (matches PopIIIEqEOS.composition_he interface)."
@inline function lte_composition_he(ρ, T)
    t = lte_table()
    x = _gx(t, log10(ρ)); y = _gy(t, log10(T))
    yHe_max = (t.Y / mHe) / (t.X / mH)
    z1 = clamp(bspline2d(t.c_zHe1, x, y), zero(ρ * T), yHe_max)
    z2 = clamp(bspline2d(t.c_zHe2, x, y), zero(ρ * T), yHe_max)
    return (z1, z2)
end

"Solve u(ρ, T) + dt * cooling(ρ, T, y)/ρ = u0 for T using our self-generated LTE table (matches PopIIIEqEOS.step)."
function lte_step(ρ::Float64, u0::Float64, T0::Float64, dt::Float64; cooling = nothing)
    F(T) = lte_eos(ρ, T)[1] - u0 + (cooling === nothing ? 0.0 : dt * cooling(ρ, T, lte_composition(ρ, T)) / ρ)
    lT = log(clamp(T0, 50.0, 1e6))
    ok = false
    for _ in 1:40
        T = exp(lT); r = F(T)
        dr = ForwardDiff.derivative(lt -> F(exp(lt)), lT)
        (!isfinite(dr) || abs(dr) < 1e-30) && break
        δ = clamp(-r / dr, -0.5, 0.5)
        lT = clamp(lT + δ, log(50.0), log(1.5e6))
        if abs(δ) < 1e-11
            ok = true; break
        end
    end
    T1 = exp(lT)
    if !ok
        # Robust bisection fallback
        lo, hi = 50.0, 1.5e6
        if F(lo) <= 0 && F(hi) >= 0
            for _ in 1:80
                m = sqrt(lo * hi); F(m) > 0 ? (hi = m) : (lo = m)
                hi / lo - 1 < 1e-10 && break
            end
            T1 = sqrt(lo * hi)
        end
    end
    return lte_composition(ρ, T1), T1
end

end # module
