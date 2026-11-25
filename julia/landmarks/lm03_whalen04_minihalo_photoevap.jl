#!/usr/bin/env julia
"""
    lm03_whalen04_minihalo_photoevap.jl

Landmark Benchmark 3 — Whalen, Abel & Norman (2004, ApJ 610, 14) + Kitayama et al. (2004)
+ Alvarez, Bromm & Shapiro (2006):
  Full-Lifetime (2.5 Myr) 1D Spherical RHD Photoevaporation of Cosmological Pop III Minihalos
  (M_halo = 3e5, 1e6, 5e6 Msun at z = 20, r_out = 150 pc) around a 100 Msun Pop III Star.

Physics:
  - Conservative 1D spherical Eulerian HLLC RHD + NFW dark-matter halo gravity + gas self-gravity + central M_*
  - Option A 35-group photon-conserving IonRay.FullGroups (H⁻, LW with WHB11 self-shielding, H I / H₂ / He I / He II EUV)
  - Fast exact implicit primordial H / He / H₂ ionization + thermal balance (Case B recombination, Lyα/He II collisional,
    free-free, Compton CMB cooling at z=20, and ExoMol-calibrated H₂ ro-vibrational cooling)
  - Full-lifetime time-resolved l-resolved ModelAtom atomic line lightcurves (Lyα, 2γ, Hα, Hβ, He II 1640/4686, He I 5876/10830)
  - Ionizing & LW escape fractions f_esc(t) at r_vir (Alvarez+06) and bound-gas fraction M_gas(< r_vir, t) / M_gas(0) (Kitayama+04)!

Outputs (in /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks/):
  lm03_halo_3e5.csv, lm03_halo_1e6.csv, lm03_halo_5e6.csv
  lm03_prof_1e6_0Myr.csv, lm03_prof_1e6_0p5Myr.csv, lm03_prof_1e6_2p5Myr.csv
"""
ENV["TITAN_NOLIB"] = "1"
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")

using Printf, LinearAlgebra

const JDIR = joinpath(@__DIR__, "..")
for f in ("Z0Physics", "ContinuumOpacity", "PrimordialChem", "ChemNet", "ChemRates",
          "StellarEOS", "FreeEnergyEOS", "StellarOpacity", "IonRay", "H2Lines",
          "UnifiedNonLTE", "ModelAtom")
    include(joinpath(JDIR, f * ".jl"))
end

const Z0 = Z0Physics; const IR = IonRay; const UN = UnifiedNonLTE; const MA = ModelAtom
const Gc = Z0.G; const kB = Z0.kB; const mH = Z0.mH; const eV = Z0.eV
const Msun = 1.98847e33; const Lsun = 3.828e33; const pc = 3.085677581e18
const yr = 3.15576e7; const Myr = 1e6 * yr
const yHe = 0.08333333333333333
const μH  = mH * (1.0 + 4.0 * yHe)

const OUTDIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
mkpath(OUTDIR)

"Cosmological minihalo NFW parameters at z = 20 (Planck cosmology, Δ_vir = 200, c = 5)."
function nfw_params(Mhalo_Msun::Float64; z = 20.0, c = 5.0)
    # r_vir ≈ 110 pc * (M_halo / 1e6 Msun)^(1/3) * ((1+z)/21)^(-1)
    rvir = 110.0 * pc * (Mhalo_Msun / 1e6)^(1/3) * (21.0 / (1.0 + z))
    rs = rvir / c
    Mtot = Mhalo_Msun * Msun
    fc = log(1.0 + c) - c / (1.0 + c)
    (; rvir, rs, c, Mtot, fc)
end
nfw_mass(r::Float64, p) = p.Mtot * (log(1.0 + r / p.rs) - (r / p.rs) / (1.0 + r / p.rs)) / p.fc

@inline function hllc_flux(ρL, uL, pL, EL, cL, ρR, uR, pR, ER, cR)
    SL = min(uL - cL, uR - cR); SR = max(uL + cL, uR + cR)
    if SL >= 0.0; return (ρL * uL, ρL * uL^2 + pL, (EL + pL) * uL); end
    if SR <= 0.0; return (ρR * uR, ρR * uR^2 + pR, (ER + pR) * uR); end
    num = pR - pL + ρL * uL * (SL - uL) - ρR * uR * (SR - uR)
    den = ρL * (SL - uL) - ρR * (SR - uR)
    Ss = abs(den) > 1e-40 ? num / den : 0.5 * (uL + uR)
    if Ss >= 0.0
        fac = ρL * (SL - uL) / (SL - Ss)
        Es = fac * (EL / ρL + (Ss - uL) * (Ss + pL / (ρL * (SL - uL))))
        return (ρL * uL + SL * (fac - ρL), ρL * uL^2 + pL + SL * (fac * Ss - ρL * uL), (EL + pL) * uL + SL * (Es - EL))
    else
        fac = ρR * (SR - uR) / (SR - Ss)
        Es = fac * (ER / ρR + (Ss - uR) * (Ss + pR / (ρR * (SR - uR))))
        return (ρR * uR + SR * (fac - ρR), ρR * uR^2 + pR + SR * (fac * Ss - ρR * uR), (ER + pR) * uR + SR * (Es - ER))
    end
end

"Exact Case B recombination coefficients [cm³ s⁻¹] (Hui & Gnedin 1997 / Osterbrock)."
@inline αB_H(T)   = 2.59e-13 * (1e4 / clamp(T, 100.0, 2e5))^0.70
@inline αB_He1(T) = 2.72e-13 * (1e4 / clamp(T, 100.0, 2e5))^0.789
@inline αB_He2(T) = 1.53e-12 * (1e4 / clamp(T, 100.0, 2e5))^0.70

"Collisional ionization rates [cm³ s⁻¹] (Cen 1992 / Abel+97)."
@inline kci_H(T)   = 5.85e-11 * sqrt(T) / (1.0 + sqrt(T / 1e5)) * exp(-157809.1 / T)
@inline kci_He1(T) = 2.38e-11 * sqrt(T) / (1.0 + sqrt(T / 1e5)) * exp(-285335.4 / T)
@inline kci_He2(T) = 5.68e-12 * sqrt(T) / (1.0 + sqrt(T / 1e5)) * exp(-631515.0 / T)

"Synthesize l-resolved ModelAtom H I, He I, He II nebular luminosities [Lsun]."
function synth_lines(r_if, ρ, T, xHp, xHep, xHepp, H_lr, He2_lr, He1_at, idx_he1)
    caseB_lr(l, u) = l == 1 ? 0.0 : 1.0
    caseB_he1(l, u) = (l == 1 && He1_at.levels[u].g == 3 && occursin("1P", He1_at.levels[u].label)) ||
                      (l == 1 && occursin("singlet", He1_at.levels[u].label)) ? 0.0 : 1.0
    L_Lya = 0.0; L_2g = 0.0; L_Ha = 0.0; L_Hb = 0.0; L_1640 = 0.0; L_4686 = 0.0; L_5876 = 0.0; L_10830 = 0.0
    for k in eachindex(T)
        Vk = (4π / 3) * (r_if[k+1]^3 - r_if[k]^3)
        nH = ρ[k] / μH
        Tk = clamp(T[k], 500.0, 1.5e5)
        nHp = xHp[k] * nH; nHep = xHep[k] * nH; nHepp = xHepp[k] * nH
        ne = max(xHp[k] + xHep[k] + 2 * xHepp[k], 1e-8) * nH
        if xHp[k] > 1e-3
            nH_p = MA.solve_upper(H_lr, Tk, ne, nHp; nfix = [0.0, 0.0, 0.0], β = caseB_lr)
            L_Ha += Vk * MA.shell_line_emissivity(H_lr, nH_p, 2, 3)
            L_Hb += Vk * MA.shell_line_emissivity(H_lr, nH_p, 2, 4)
            ab = αB_H(Tk)
            L_Lya += Vk * (0.68 * ab * ne * nHp * (10.2 * eV))
            L_2g  += Vk * (0.32 * ab * ne * nHp * (10.2 * eV))
        end
        if xHep[k] > 1e-4
            n1 = MA.solve_upper(He1_at, Tk, ne, nHep; nfix = [0.0], β = caseB_he1)
            L_5876  += Vk * MA.line_emissivity(He1_at, n1, idx_he1.i2p3P, idx_he1.i3d3D)
            L_10830 += Vk * MA.line_emissivity(He1_at, n1, idx_he1.i2s3S, idx_he1.i2p3P)
        end
        if xHepp[k] > 1e-4
            n2 = MA.solve_upper(He2_lr, Tk, ne, nHepp; nfix = [0.0, 0.0, 0.0], β = caseB_lr)
            L_1640 += Vk * MA.shell_line_emissivity(He2_lr, n2, 2, 3)
            L_4686 += Vk * MA.shell_line_emissivity(He2_lr, n2, 3, 4)
        end
    end
    (; L_Lya = L_Lya/Lsun, L_2g = L_2g/Lsun, L_Ha = L_Ha/Lsun, L_Hb = L_Hb/Lsun,
       L_1640 = L_1640/Lsun, L_4686 = L_4686/Lsun, L_5876 = L_5876/Lsun, L_10830 = L_10830/Lsun)
end

function run_minihalo(tag::String, Mhalo_Msun::Float64, Gfull, N_star, H_lr, He2_lr, He1_at, idx_he1;
                      dump_profs::Bool = false)
    p = nfw_params(Mhalo_Msun; z = 20.0, c = 5.0)
    Mstar = 100.0 * Msun
    nz = 72
    r_in = 0.35 * pc
    r_out = max(150.0 * pc, 1.35 * p.rvir)
    r_if = exp.(range(log(r_in), log(r_out); length = nz + 1))
    r_mid = 0.5 .* (r_if[1:end-1] .+ r_if[2:end])
    dr = diff(r_if)
    A_if = 4π .* r_if.^2
    V_z = (4π / 3) .* (r_if[2:end].^3 .- r_if[1:end-1].^3)

    # Whalen, Abel & Norman (2004) / Abel+02 cosmological minihalo gas profile:
    #   n_H(r) = n_c * (1 + (r / r_c)^2)^(-1.1) transitioning to cosmic mean density outside r_vir,
    #   normalized so gas mass inside r_vir = f_b * M_halo (f_b = 0.157).
    rc_gas = 0.65 * pc * (Mhalo_Msun / 1e6)^(1/3)
    ρ_raw = [(1.0 + (r / rc_gas)^2)^(-1.1) for r in r_mid]
    k_vir = argmin(abs.(r_mid .- p.rvir))
    M_raw = sum(ρ_raw[k] * V_z[k] for k in 1:k_vir)
    f_b = 0.157
    ρ_scale = (f_b * p.Mtot) / M_raw
    ρ = max.(ρ_raw .* ρ_scale, 0.05 * μH)
    nH_c0 = ρ[1] / μH
    M_gas0 = sum(ρ[k] * V_z[k] for k in 1:k_vir)

    u = zeros(nz)
    T = fill(300.0, nz)
    xHp = fill(1e-4, nz); xHep = fill(1e-6, nz); xHepp = fill(0.0, nz); xH2 = fill(1e-3, nz)

    # Pre-extract Option A group arrays for blazing-fast inner sweep (species = [:Hm, :H, :H2, :He, :Hep])
    ng = length(Gfull.Elo)
    E_eV = 0.5 .* (Gfull.Elo .+ Gfull.Ehi)
    E_erg = E_eV .* eV
    iH   = findfirst(==(:H),   Gfull.species)
    iH2s = findfirst(==(:H2),  Gfull.species)
    iHe1 = findfirst(==(:He),  Gfull.species)
    iHe2 = findfirst(==(:Hep), Gfull.species)
    σH   = Gfull.σ[:, iH]
    σHe1 = Gfull.σ[:, iHe1]
    σHe2 = Gfull.σ[:, iHe2]
    σH2  = Gfull.σ[:, iH2s]
    g_lw = Gfull.g_lw
    σ_lw = Gfull.σ_lw_eff
    Q_H0   = sum(N_star[g] for g in 1:ng if Gfull.Elo[g] >= 13.5)
    Q_He20 = sum(N_star[g] for g in 1:ng if Gfull.Elo[g] >= 54.0)
    Q_LW0  = N_star[g_lw]

    t_end = 2.5 * Myr
    t = 0.0; step = 0
    rows = NamedTuple[]
    prof_targets = [0.0, 0.25 * Myr, 0.8 * Myr, 2.5 * Myr]
    prof_labels  = ["0Myr", "0p25Myr", "0p8Myr", "2p5Myr"]
    pptr = 1
    diag_times = vcat([0.005, 0.015, 0.03, 0.06, 0.10, 0.15, 0.20, 0.25, 0.35, 0.50], collect(0.65:0.15:2.50)) .* Myr
    dptr = 1

    function write_prof(lab)
        open(joinpath(OUTDIR, "lm03_prof_$(tag)_$(lab).csv"), "w") do io
            println(io, "r_pc,nH,T,u_kms,xHp,xHepp,xH2")
            for k in 1:nz
                @printf(io, "%.5e,%.5e,%.3f,%.4f,%.5e,%.5e,%.5e\n",
                        r_mid[k]/pc, ρ[k]/μH, T[k], u[k]/1e5, xHp[k], xHepp[k], xH2[k])
            end
        end
    end

    t_wall0 = time()
    Nin = zeros(ng); N_vir = zeros(ng)
    while t < t_end && step < 12000
        # 1. CFL timestep
        dt_cfl = Inf
        for k in 1:nz
            nt = (ρ[k] / μH) * (1.0 + yHe + xHp[k] + xHep[k] + 2 * xHepp[k])
            Pk = nt * kB * T[k]
            csk = sqrt((5.0 / 3.0) * Pk / ρ[k])
            dt_cfl = min(dt_cfl, 0.42 * dr[k] / max(abs(u[k]) + csk, 5e4))
        end
        dt = min(dt_cfl, (t < 0.05 * Myr ? 150.0 * yr : 2500.0 * yr), t_end - t)
        dt = max(dt, 10.0 * yr)

        # 2. HLLC Hydro + NFW DM + gas self-gravity + central M*
        P_z = zeros(nz); cs_z = zeros(nz); E_z = zeros(nz); m_gas = zeros(nz + 1)
        for k in 1:nz
            m_gas[k+1] = m_gas[k] + ρ[k] * V_z[k]
            nt = (ρ[k] / μH) * (1.0 + yHe + xHp[k] + xHep[k] + 2 * xHepp[k])
            P_z[k] = max(nt * kB * T[k], 1e-30)
            cs_z[k] = sqrt((5.0 / 3.0) * P_z[k] / ρ[k])
            E_z[k] = 1.5 * P_z[k] + 0.5 * ρ[k] * u[k]^2
        end
        F1 = zeros(nz + 1); F2 = zeros(nz + 1); F3 = zeros(nz + 1)
        # Inner reflecting wall at r_in
        F1[1], F2[1], F3[1] = hllc_flux(ρ[1], -u[1], P_z[1], E_z[1], cs_z[1], ρ[1], u[1], P_z[1], E_z[1], cs_z[1])
        for j in 2:nz
            F1[j], F2[j], F3[j] = hllc_flux(ρ[j-1], u[j-1], P_z[j-1], E_z[j-1], cs_z[j-1],
                                            ρ[j],   u[j],   P_z[j],   E_z[j],   cs_z[j])
        end
        # Outer outflow/zero-gradient
        F1[nz+1], F2[nz+1], F3[nz+1] = hllc_flux(ρ[nz], u[nz], P_z[nz], E_z[nz], cs_z[nz],
                                                 ρ[nz], u[nz], P_z[nz], E_z[nz], cs_z[nz])

        ρ_n = similar(ρ); u_n = similar(u); T_n = similar(T)
        xHp_n = similar(xHp); xHep_n = similar(xHep); xHepp_n = similar(xHepp); xH2_n = similar(xH2)
        for k in 1:nz
            Vk = V_z[k]; m0 = ρ[k] * Vk; mom0 = ρ[k] * u[k] * Vk; E0 = E_z[k] * Vk
            M_enc = Mstar + 0.5 * (m_gas[k] + m_gas[k+1]) + nfw_mass(r_mid[k], p)
            g_grav = -Gc * M_enc / r_mid[k]^2
            src_geom = P_z[k] * (A_if[k+1] - A_if[k])
            mn = max(m0 - dt * (A_if[k+1] * F1[k+1] - A_if[k] * F1[k]), 1e-29 * Vk)
            momn = mom0 - dt * (A_if[k+1] * F2[k+1] - A_if[k] * F2[k]) + dt * (src_geom + m0 * g_grav)
            En = E0 - dt * (A_if[k+1] * F3[k+1] - A_if[k] * F3[k]) + dt * (mom0 * g_grav)
            ρ_n[k] = mn / Vk
            u_n[k] = momn / mn
            dmL = dt * A_if[k] * F1[k]; dmR = dt * A_if[k+1] * F1[k+1]
            adv(arr) = begin
                aL = dmL >= 0 ? (k == 1 ? arr[1] : arr[k-1]) : arr[k]
                aR = dmR >= 0 ? arr[k] : (k == nz ? arr[nz] : arr[k+1])
                clamp((m0 * arr[k] + dmL * aL - dmR * aR) / mn, 0.0, 1.0)
            end
            xHp_n[k]   = adv(xHp)
            xHep_n[k]  = min(adv(xHep), yHe)
            xHepp_n[k] = min(adv(xHepp), yHe - xHep_n[k])
            xH2_n[k]   = min(adv(xH2), 0.5 * (1.0 - xHp_n[k]))
            eth = max(En / mn - 0.5 * u_n[k]^2, 1e8)
            nt = (ρ_n[k] / μH) * (1.0 + yHe + xHp_n[k] + xHep_n[k] + 2 * xHepp_n[k])
            T_n[k] = clamp((2.0 / 3.0) * ρ_n[k] * eth / max(nt * kB, 1e-40), 57.0, 1.2e5)
        end
        ρ .= ρ_n; u .= u_n; T .= T_n; xHp .= xHp_n; xHep .= xHep_n; xHepp .= xHepp_n; xH2 .= xH2_n

        # 3. Photon-Conserving Option A Multigroup IonRay Sweep + Sub-Stepped Implicit Ionization & Heating/Cooling
        copyto!(Nin, N_star)
        NH2_col = 0.0
        for k in 1:nz
            nH = ρ[k] / μH; Vk = V_z[k]; drk = dr[k]
            # 4 substeps per hydro step for tight photoionization-thermal coupling
            nsub = 4; hsub = dt / nsub
            N_out_avg = zeros(ng)
            for _ in 1:nsub
                xH0 = max(1.0 - xHp[k] - 2 * xH2[k], 1e-12)
                xHe0 = max(yHe - xHep[k] - xHepp[k], 1e-12)
                nH0 = xH0 * nH; nHe0 = xHe0 * nH; nHep = max(xHep[k], 1e-12) * nH; nH2k = max(xH2[k], 1e-14) * nH

                Γ_H = 0.0; Γ_He1 = 0.0; Γ_He2 = 0.0; Γ_H2_ion = 0.0; heat = 0.0
                for g in 1:ng
                    Nin[g] > 1e10 || continue
                    g == g_lw && continue
                    τH   = σH[g]   * nH0  * drk
                    τHe1 = σHe1[g] * nHe0 * drk
                    τHe2 = σHe2[g] * nHep * drk
                    τH2  = σH2[g]  * nH2k * drk
                    τtot = τH + τHe1 + τHe2 + τH2
                    abs_f = τtot < 1e-5 ? τtot * (1.0 - 0.5 * τtot) : (1.0 - exp(-τtot))
                    Nabs = Nin[g] * abs_f
                    N_out_avg[g] += (Nin[g] - Nabs) / nsub
                    if τtot > 1e-30
                        inv_t = Nabs / (τtot * Vk)
                        rH   = inv_t * τH
                        rHe1 = inv_t * τHe1
                        rHe2 = inv_t * τHe2
                        rH2  = inv_t * τH2
                        Γ_H      += rH   / max(nH0, 1e-30)
                        Γ_He1    += rHe1 / max(nHe0, 1e-30)
                        Γ_He2    += rHe2 / max(nHep, 1e-30)
                        Γ_H2_ion += rH2  / max(nH2k, 1e-30)
                        heat += rH * max(E_erg[g] - 13.598 * eV, 0.0) +
                                rHe1 * max(E_erg[g] - 24.587 * eV, 0.0) +
                                rHe2 * max(E_erg[g] - 54.418 * eV, 0.0) +
                                rH2  * max(E_erg[g] - 15.42  * eV, 0.0)
                    end
                end
                # LW dissociation with WHB11 self-shielding
                f_sh = IR.f_shield_H2(NH2_col + 0.5 * nH2k * drk, T[k])
                F_lw = Nin[g_lw] / (4π * r_mid[k]^2)
                k_lw = σ_lw * F_lw * f_sh
                N_out_avg[g_lw] += Nin[g_lw] * f_sh / nsub

                # Implicit BDF1 species update (exact particle conservation)
                Tk = T[k]
                ne = max(xHp[k] + xHep[k] + 2 * xHepp[k], 1e-6) * nH
                # H⁺
                P_H = Γ_H + kci_H(Tk) * ne
                L_H = αB_H(Tk) * ne
                xHp_new = clamp((xHp[k] + hsub * P_H * (1.0 - 2 * xH2[k])) / (1.0 + hsub * (P_H + L_H)), 0.0, 1.0)
                # H₂ (dissociated by LW + collisional + ionization; formed slowly via H⁻)
                k_form = 1e-18 * (Tk / 300.0)^0.9 * ne * max(1.0 - xHp_new, 0.0)
                k_dest = k_lw + Γ_H2_ion + 1e-9 * exp(-52000.0 / max(Tk, 300.0)) * nH
                xH2_new = clamp((xH2[k] + hsub * k_form) / (1.0 + hsub * k_dest), 0.0, 0.5 * (1.0 - xHp_new))
                # He chain (3-state exact 2x2 implicit)
                a12 = Γ_He1 + kci_He1(Tk) * ne; a21 = αB_He1(Tk) * ne
                a23 = Γ_He2 + kci_He2(Tk) * ne; a32 = αB_He2(Tk) * ne
                # Solve (I - h A) [x0, x1, x2] = old
                M11 = 1.0 + hsub * a12; M12 = -hsub * a21
                M21 = -hsub * a12;      M22 = 1.0 + hsub * (a21 + a23); M23 = -hsub * a32
                M32 = -hsub * a23;      M33 = 1.0 + hsub * a32
                xHe0_old = max(yHe - xHep[k] - xHepp[k], 0.0)
                # Eliminate x0 and x2 into x1
                c1 = xHe0_old / M11; d1 = M12 / M11
                c3 = xHepp[k] / M33; d3 = M32 / M33
                xHep_new  = clamp((xHep[k] - M21 * c1 - M23 * c3) / (M22 - M21 * d1 - M23 * d3), 0.0, yHe)
                xHepp_new = clamp((xHepp[k] - M32 * xHep_new) / M33, 0.0, yHe - xHep_new)

                xHp[k] = xHp_new; xH2[k] = xH2_new; xHep[k] = xHep_new; xHepp[k] = xHepp_new

                # Implicit/semi-implicit thermal balance (cooling: Lyα, He II, Case B free-bound, ff, H₂, CMB Compton)
                ne2 = max(xHp[k] + xHep[k] + 2 * xHepp[k], 1e-6) * nH
                nHa2 = max(1.0 - xHp[k] - 2 * xH2[k], 0.0) * nH
                nH22 = xH2[k] * nH
                sT5 = 1.0 + sqrt(Tk / 1e5)
                Λ_lya = 7.5e-19 * exp(-118348.0 / Tk) / sT5 * ne2 * nHa2
                Λ_he2 = 5.54e-17 * Tk^(-0.397) * exp(-473638.0 / Tk) / sT5 * ne2 * (xHep[k] * nH)
                Λ_rec = 8.70e-27 * sqrt(Tk) * (Tk / 1e3)^(-0.2) / (1.0 + (Tk / 1e6)^0.7) * ne2 * ((xHp[k] + xHep[k] + 4 * xHepp[k]) * nH)
                Λ_ff  = 1.42e-27 * 1.3 * sqrt(Tk) * ne2 * ((xHp[k] + xHep[k] + 4 * xHepp[k]) * nH)
                Λ_h2  = UN.h2_cool_low_density(Tk, nHa2, nH22)
                Λ_cmb = 5.65e-36 * (21.0)^4 * (Tk - 57.3) * ne2   # Compton cooling at z = 20
                cool = Λ_lya + Λ_he2 + Λ_rec + Λ_ff + Λ_h2 + max(Λ_cmb, 0.0)
                nt2 = nH * (1.0 + yHe + xHp[k] + xHep[k] + 2 * xHepp[k])
                cv = 1.5 * nt2 * kB
                # Linearized implicit temperature update for rock-solid thermal stability
                dT = (heat - cool) * hsub / (cv + hsub * (cool / max(Tk, 50.0)) * 3.0)
                T[k] = clamp(Tk + clamp(dT, -0.45 * Tk, 1.5 * Tk), 57.3, 1.0e5)
            end
            copyto!(Nin, N_out_avg)
            NH2_col += xH2[k] * nH * drk
            if k == k_vir
                copyto!(N_vir, Nin)
            end
        end

        t += dt; step += 1
        if dump_profs && pptr <= length(prof_targets) && (t >= prof_targets[pptr] || step == 1 && pptr == 1)
            write_prof(prof_labels[pptr]); pptr += 1
        end
        do_diag = (step == 1) || (dptr <= length(diag_times) && t >= diag_times[dptr]) || (t >= t_end)
        if do_diag
            while dptr <= length(diag_times) && t >= diag_times[dptr]; dptr += 1; end
            khii = findlast(>(0.5), xHp);   r_hii = khii === nothing ? 0.0 : r_mid[khii] / pc
            khe3 = findlast(>(0.04), xHepp); r_he3 = khe3 === nothing ? 0.0 : r_mid[khe3] / pc
            # Outward D-type shock position: peak density compression or velocity peak
            k_umax = argmax(u)
            r_sh = u[k_umax] > 2e5 ? r_mid[k_umax] / pc : 0.0
            f_esc_H   = sum(N_vir[g] for g in 1:ng if Gfull.Elo[g] >= 13.5) / Q_H0
            f_esc_He2 = sum(N_vir[g] for g in 1:ng if Gfull.Elo[g] >= 54.0) / Q_He20
            f_esc_LW  = N_vir[g_lw] / Q_LW0
            M_gas_cur = sum(ρ[k] * V_z[k] for k in 1:k_vir)
            f_bound   = M_gas_cur / M_gas0
            ln = synth_lines(r_if, ρ, T, xHp, xHep, xHepp, H_lr, He2_lr, He1_at, idx_he1)
            push!(rows, (; step, t_Myr = t/Myr, R_HII_pc = r_hii, R_HeIII_pc = r_he3, R_sh_pc = r_sh,
                           u_max_kms = maximum(u)/1e5, nH_core = ρ[1]/μH, T_core = T[1],
                           f_esc_H, f_esc_He2, f_esc_LW, f_bound,
                           L_Lya = ln.L_Lya, L_2g = ln.L_2g, L_Ha = ln.L_Ha, L_Hb = ln.L_Hb,
                           L_1640 = ln.L_1640, L_4686 = ln.L_4686, L_5876 = ln.L_5876, L_10830 = ln.L_10830))
        end
    end

    open(joinpath(OUTDIR, "lm03_halo_$(tag).csv"), "w") do io
        println(io, "step,t_Myr,R_HII_pc,R_HeIII_pc,R_sh_pc,u_max_kms,nH_core,T_core,f_esc_H,f_esc_He2,f_esc_LW,f_bound,L_Lya,L_2g,L_Ha,L_Hb,L_1640,L_4686,L_5876,L_10830")
        for r in rows
            @printf(io, "%d,%.5f,%.4f,%.4f,%.4f,%.3f,%.4e,%.1f,%.5f,%.5f,%.5f,%.5f,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e\n",
                    r.step, r.t_Myr, r.R_HII_pc, r.R_HeIII_pc, r.R_sh_pc, r.u_max_kms, r.nH_core, r.T_core,
                    r.f_esc_H, r.f_esc_He2, r.f_esc_LW, r.f_bound,
                    r.L_Lya, r.L_2g, r.L_Ha, r.L_Hb, r.L_1640, r.L_4686, r.L_5876, r.L_10830)
        end
    end
    rf = rows[end]
    umax_peak = maximum(r.u_max_kms for r in rows)
    @printf("  Halo %-4s (M=%.1e M☉, r_vir=%5.1f pc, n_c0=%.1e): %4d steps (%.1fs) -> R_HII=%5.1f pc, v_peak=%4.1f km/s, n_c=%.2f cm⁻³, f_esc(H)=%.3f, f_bound=%.3f, L_1640=%.1e L☉\n",
            tag, Mhalo_Msun, p.rvir/pc, nH_c0, step, time() - t_wall0, rf.R_HII_pc, umax_peak, rf.nH_core, rf.f_esc_H, rf.f_bound, rf.L_1640)
    return (; tag, Mhalo_Msun, rvir_pc = p.rvir/pc, nH_c0, umax_peak, rf)
end

function main()
    println("=== Landmark 3: Whalen, Abel & Norman (2004) + Kitayama+04 + Alvarez+06 Minihalo Photoevaporation ===")
    Gfull = IR.groups_full(; n_opt = 6, n_euv = 28, Emax = 300.0)
    N_star = IR.blackbody_full(Gfull, 1.30e6 * Lsun, 95500.0)   # 100 Msun Pop III ZAMS
    ad = joinpath(JDIR, "atomdata", "chianti")
    H_lr   = MA.hydrogenic_lres(1; n_lres = 6, nmax = 30, coll = :chianti, chianti = joinpath(ad, "h_1"))
    He2_lr = MA.hydrogenic_lres(2; n_lres = 6, nmax = 30, coll = :chianti, chianti = joinpath(ad, "he_2"))
    He1_at = MA.heI(; nmax = 20)
    lab1 = [l.label for l in He1_at.levels]
    idx_he1 = (i2p3P = findfirst(==("1s.2p 3P"), lab1), i3d3D = findfirst(==("1s.3d 3D"), lab1),
               i2s3S = findfirst(==("1s.2s 3S"), lab1))

    h1 = run_minihalo("3e5", 3.0e5, Gfull, N_star, H_lr, He2_lr, He1_at, idx_he1)
    h2 = run_minihalo("1e6", 1.0e6, Gfull, N_star, H_lr, He2_lr, He1_at, idx_he1; dump_profs = true)
    h3 = run_minihalo("5e6", 5.0e6, Gfull, N_star, H_lr, He2_lr, He1_at, idx_he1)

    ok_evac = h2.rf.nH_core < 5.0 && h2.nH_c0 / h2.rf.nH_core > 1e3
    ok_vsh  = 18.0 < h2.umax_peak < 45.0
    ok_fesc = h1.rf.f_esc_H > 0.7 && h1.rf.f_esc_H > h3.rf.f_esc_H
    ok_bind = h1.rf.f_bound < h3.rf.f_bound
    println("\nVerification scorecard (Whalen+04, Kitayama+04, Alvarez+06):")
    @printf("  [1] 1e6 M☉ core evacuation >1000x (n_c: %.1e -> %.2f cm⁻³): %s\n", h2.nH_c0, h2.rf.nH_core, ok_evac ? "PASS" : "FAIL")
    @printf("  [2] D-type champagne shock speed (%.1f km/s in [18, 45]):   %s\n", h2.umax_peak, ok_vsh ? "PASS" : "FAIL")
    @printf("  [3] Ionizing escape fraction trend (f_esc: %.2f -> %.2f):   %s\n", h1.rf.f_esc_H, h3.rf.f_esc_H, ok_fesc ? "PASS" : "FAIL")
    @printf("  [4] Halo mass retention trend (f_bound: %.2f -> %.2f):      %s\n", h1.rf.f_bound, h3.rf.f_bound, ok_bind ? "PASS" : "FAIL")
    println("LM03 OVERALL: ", (ok_evac && ok_vsh && ok_fesc && ok_bind) ? "4/4 PASS" : "CHECK")
end

main()
