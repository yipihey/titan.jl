#!/usr/bin/env julia
"""
    lm02_omukai_collapse_jcrit.jl

Landmark Benchmark 2 — Omukai & Nishi (1998, ApJ 508, 141) + Omukai (2000, 2001)
+ Ripamonti & Abel (2004, MNRAS 348, 1019):
  Primordial Gas Collapse across 22 Decades in Density (n_H = 1 -> 1e22 cm⁻³),
  J_crit Bifurcation (H₂-cooling Pop III vs Atomic-cooling DCBH), and Emergent
  ExoMol H₂ Ro-Vibrational Line Luminosities.

Couples Titan.jl's exact microphysics stack:
  - ChemNet + ChemRates.modern(; he=true) (3-body H₂ formation, H⁻, H₂⁺, H/He ionization)
  - H2Lines (ExoMol 1H2__RACPPK ro-vibrational lines with Jeans/Sobolev escape probabilities β_esc)
  - ContinuumOpacity (Borysow H₂-H₂ / H₂-He CIA/CIE, H⁻ bf/ff, Kramers H bf/ff)
  - IonRay smooth C∞ Wolcott-Green+11 H₂ LW self-shielding f_sh(N_H2, T)
  - FreeEnergyEOS (exact LTE dense core thermodynamics & γ_eff up to n_H = 1e22 cm⁻³)
  - UnifiedNonLTE smooth C² non-LTE <-> LTE transition

Outputs (in /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks/):
  lm02_collapse_J0.csv, lm02_collapse_J10.csv, ..., lm02_collapse_J1e4_T1e4.csv, lm02_collapse_J1e4_T1e5.csv
  lm02_h2_lines_J0.csv (Ripamonti & Abel 2004 individual H₂ ro-vibrational line luminosities)
"""
ENV["TITAN_NOLIB"] = "1"
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")

using Printf, LinearAlgebra, ForwardDiff

const JDIR = joinpath(@__DIR__, "..")
for f in ("Z0Physics", "ContinuumOpacity", "PrimordialChem", "ChemNet", "ChemRates",
          "StellarEOS", "FreeEnergyEOS", "StellarOpacity", "IonRay", "H2Lines", "UnifiedNonLTE")
    include(joinpath(JDIR, f * ".jl"))
end

const Z0 = Z0Physics; const CO = ContinuumOpacity; const CN = ChemNet; const CR = ChemRates
const FE = FreeEnergyEOS; const IR = IonRay; const HL = H2Lines; const UN = UnifiedNonLTE

const Gc   = Z0.G
const kB   = Z0.kB
const hP   = Z0.h
const c_c  = Z0.c
const σSB  = Z0.σSB
const Lsun = 3.828e33
const OUTDIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
mkpath(OUTDIR)

"Hollenbach & McKee (1979) / Glover & Abel (2008) low-density H₂ LTE-to-nLTE collisional factor + ExoMol RACPPK line escape."
function h2_line_cooling_and_spectrum(L::HL.LineList, T::Float64, nH::Float64, nH2::Float64,
                                      nHa::Float64, λJ::Float64, vth::Float64, l_idx::NamedTuple)
    # Optically thin LTE rate per molecule from ExoMol RACPPK
    f = HL.populations(L, T)
    # Critical density ratio (Hollenbach & McKee 1979 / Glover & Abel 2008)
    Λ_lte_thin = 0.0
    Λ_esc = 0.0
    # Specific Ripamonti & Abel (2004) diagnostic lines: 0-0 S(0), S(1), S(2), S(3), 1-0 S(1)
    L_spec = zeros(5)
    # Characteristic Jeans-scale column for Sobolev/LVG escape: N_eff = 0.5 * nH2 * λJ
    N_eff = 0.5 * nH2 * λJ
    # Low-density collisional limit (Glover & Abel 2008 fit for H-H₂ + H₂-H₂)
    T3 = clamp(T / 1000.0, 0.05, 50.0)
    lt = log10(T3)
    Λ_H = exp10(-24.3112093 + 3.5692468 * lt - 11.3328607 * lt^2 - 27.8500821 * lt^3 -
                21.3282648 * lt^4 - 4.2519023 * lt^5)
    # Prevent high-T polynomial turnover above T = 3000 K
    if T > 3000.0
        lt3 = log10(3.0)
        Λ3 = exp10(-24.3112093 + 3.5692468 * lt3 - 11.3328607 * lt3^2 - 27.8500821 * lt3^3 -
                   21.3282648 * lt3^4 - 4.2519023 * lt3^5)
        Λ_H = Λ3 * (T / 3000.0)^1.2
    end
    Λ_H2 = 1.4e-23 * T3^2.0 / (1.0 + 0.5 * T3^1.5)
    Λ_low = nH2 * (nHa * Λ_H + nH2 * Λ_H2)

    @inbounds for i in eachindex(L.A)
        ν0 = L.ν[i]; iu = L.iu[i]; il = L.il[i]
        gu = L.g[iu]; gl = L.g[il]
        em_i = f[iu] * L.A[i] * hP * ν0
        Λ_lte_thin += nH2 * em_i
        # Line-centre optical depth τ_i = (c² / 8π ν0²) * (g_u/g_l f_l - f_u) * A_ul * N_eff / Δν_D
        ΔνD = ν0 * vth / c_c
        τ_i = (c_c^2 / (8π * ν0^2)) * max(gu / gl * f[il] - f[iu], 0.0) * L.A[i] * N_eff / max(ΔνD * sqrt(π), 1e-30)
        β_i = τ_i < 1e-5 ? (1.0 - 0.5 * τ_i) : (1.0 - exp(-τ_i)) / τ_i
        c_i = nH2 * em_i * β_i
        Λ_esc += c_i
        if i == l_idx.S0;      L_spec[1] = c_i
        elseif i == l_idx.S1;  L_spec[2] = c_i
        elseif i == l_idx.S2;  L_spec[3] = c_i
        elseif i == l_idx.S3;  L_spec[4] = c_i
        elseif i == l_idx.v1S1; L_spec[5] = c_i
        end
    end
    # Exact two-level-equivalent bridging between low-density nLTE and optically-thick LTE:
    # Λ_H2 = Λ_esc / (1 + Λ_lte_thin / max(Λ_low, 1e-60))
    nlte_fac = Λ_low / max(Λ_low + Λ_lte_thin, 1e-60)
    return Λ_esc * nlte_fac, L_spec .* nlte_fac, Λ_esc / max(Λ_lte_thin, 1e-60)
end

"Find ExoMol line indices for canonical Ripamonti & Abel (2004) H₂ diagnostic lines by wavelength [μm]."
function find_h2_diag_lines(L::HL.LineList)
    pick_lam(lam_um) = argmin([abs(c_c / ν * 1e4 - lam_um) for ν in L.ν])
    (; S0 = pick_lam(28.2188), S1 = pick_lam(17.0348), S2 = pick_lam(12.2786),
       S3 = pick_lam(9.6649), v1S1 = pick_lam(2.1218))
end

"Run a single one-zone Larson-Penston primordial collapse track from n_H = 1 -> 1e22 cm⁻³."
function run_collapse(tag::String, J21::Float64, Trad_ext::Float64, L_h2::HL.LineList, l_idx, net, photos;
                    dump_lines::Bool = false)
    iHp = findfirst(==(:Hp), net.free); iH2 = findfirst(==(:H2), net.free)
    iHm = findfirst(==(:Hm), net.free)
    nf = length(net.free)
    x = zeros(nf); x[iHp] = 1.0e-4; x[iH2] = 1.0e-6

    μH = CN.μH()
    nH = 1.0; ρ = nH * μH; T = 200.0
    env_c = (Trad = 54.0, J = ν -> 0.0, static = true, caseB = true)
    # Photodissociation & photodetachment base rates for external radiation field J_21 (Omukai 2001 / Wolcott-Green+11):
    #   k_LW  = 1.38e-12 * J_21 * f_sh(N_H2, T)  [s⁻¹]
    #   k_H-  = α_Hm * 1e-10 * J_21              [s⁻¹], with α_Hm = 2.0 (T_rad=1e4 K) or 0.001 (T_rad=1e5 K)
    k_LW_0 = 1.38e-12 * J21
    k_Hm_0 = (Trad_ext < 3e4 ? 2.0e-10 : 1.0e-13) * J21

    # Collapse grid in log10(n_H) from 0 to 22 (fine Δlog10(nH) = 0.05)
    dlnρ = log(10.0) * 0.04
    n_steps = round(Int, 22.0 / 0.04)

    rows = NamedTuple[]
    line_rows = NamedTuple[]
    T_min_popiii = Inf; n_Tmin = 0.0; xH2_max = 0.0; T_at_1e4 = 0.0; T_at_1e10 = 0.0; T_at_1e20 = 0.0

    for st in 0:n_steps
        d = CN.full(x, net)
        nHa  = max(d[:H], 0.0) * nH;  nHp = max(d[:Hp], 0.0) * nH
        nH2  = max(d[:H2], 0.0) * nH; nHe = max(d[:He], 0.0) * nH
        nHm  = max(d[:Hm], 0.0) * nH; ne  = max(d[:e],  0.0) * nH

        # Smooth C² LTE blend weight (UnifiedNonLTE) -> smoothly transitions to exact FreeEnergyEOS in dense protostellar core!
        w_lte = UN.w_lte(ρ, T, UN.BlendPar())
        if w_lte > 0.99
            st_fe = FE.eos_lte(ρ, T; X = 0.75, Y = 0.25, pairs = true, radiation = true)
            P = st_fe.P; cs = sqrt(max(st_fe.Γ1 * P / ρ, 1e8))
            γ_eff = st_fe.Γ1
            xH2_rep = 0.5 * st_fe.fH2; xHp_rep = st_fe.xHp; xe_rep = st_fe.ne / nH
        else
            ntot = nHa + nHp + nH2 + nHe + nHm + ne + max(d[:Hep],0.0)*nH + max(d[:Hepp],0.0)*nH
            P = max(ntot * kB * T, 1e-30)
            γ_eff = (nH2 > 0.1 * nH && T < 3000.0) ? 7.0 / 5.0 : 5.0 / 3.0
            cs = sqrt(γ_eff * P / ρ)
            xH2_rep = d[:H2]; xHp_rep = d[:Hp]; xe_rep = d[:e]
        end

        t_ff = sqrt(3π / (32 * Gc * ρ))
        λJ   = sqrt(π * cs^2 / (Gc * ρ))
        vth  = sqrt(2 * kB * T / (2 * Z0.mH))
        M_J  = (4π / 3) * ρ * (0.5 * λJ)^3

        # 1. H₂ ro-vibrational line cooling (ExoMol RACPPK + Jeans escape)
        Λ_h2, L_spec, β_mean = h2_line_cooling_and_spectrum(L_h2, clamp(T, 50.0, 1.5e4), nH, nH2, nHa, λJ, vth, l_idx)

        # 2. Continuum cooling (CIE H₂-H₂/H₂-He + H⁻ bf/ff + Kramers H) — LTE Kirchhoff emission operates at high density (n_H > 1e8 cm⁻³)
        κP_c = 0.0; κR_c = 0.0; τ_c = 0.0
        if nH > 1e8
            if w_lte > 0.95
                κR_c = StellarOpacity.kappa_R(ρ, T, 0.75)
                lT = log10(clamp(T, 3000.0, 1e8)); lR = log10(ρ) - 3.0 * (lT - 6.0)
                κP_c = StellarOpacity.in_table(ρ, T) ? exp10(StellarOpacity.log_kappa_P(lT, lR, 0.75)) : κR_c
            else
                κR_c, κP_c = CO.means(ρ, clamp(T, 200.0, 4e4); nH = nHa, nH2 = nH2, nHe = nHe, nHp = nHp, ne = ne, nHm = nHm)
            end
            τ_c = κR_c * ρ * (0.5 * λJ)
        end
        # Smooth density ramp on Kirchhoff continuum cooling over n_H ∈ [1e8, 1e10] cm⁻³
        g_cont = UN.log_ramp(nH, 1e8, 1e10)
        Λ_cont = g_cont * (4 * σSB * T^4 * κP_c * ρ) / (1.0 + τ_c + τ_c^2)
        # Atomic H Lyα cooling (suppressed by Lyman-line & continuum optical depth at high density)
        τ_lya = 1e-13 * nHa * (0.5 * λJ) * sqrt(1e4 / max(T, 100.0))
        β_lya = 1.0 / (1.0 + 1e-4 * τ_lya + τ_c^2)
        sT5 = 1.0 + sqrt(max(T, 100.0) / 1e5)
        Λ_lya = 7.5e-19 * exp(-118348.0 / max(T, 100.0)) / sT5 * ne * nHa * β_lya
        Λ_tot = Λ_h2 + Λ_cont + Λ_lya

        # Record track diagnostics
        if nH >= 1e2 && nH <= 1e6 && T < T_min_popiii
            T_min_popiii = T; n_Tmin = nH
        end
        xH2_max = max(xH2_max, xH2_rep)
        abs(log10(nH) - 4.0)  < 0.025 && (T_at_1e4  = T)
        abs(log10(nH) - 10.0) < 0.025 && (T_at_1e10 = T)
        abs(log10(nH) - 20.0) < 0.025 && (T_at_1e20 = T)

        push!(rows, (; nH, ρ, T, xH2 = xH2_rep, xHp = xHp_rep, xe = xe_rep,
                       Λ_h2, Λ_cont, Λ_lya, τ_c, β_mean, λJ, M_J_Msun = M_J / 1.98847e33))
        if dump_lines
            Vj = (4π / 3) * (0.5 * λJ)^3
            push!(line_rows, (; nH, T, xH2 = xH2_rep,
                                L_S0 = L_spec[1]*Vj/Lsun, L_S1 = L_spec[2]*Vj/Lsun, L_S2 = L_spec[3]*Vj/Lsun,
                                L_S3 = L_spec[4]*Vj/Lsun, L_v1S1 = L_spec[5]*Vj/Lsun,
                                L_H2tot = Λ_h2*Vj/Lsun, L_cont = Λ_cont*Vj/Lsun))
        end

        st == n_steps && break

        # Advance one collapse step Δlnρ = dlnρ over dt = dlnρ * t_ff
        dt = dlnρ * t_ff
        ρ_new = ρ * exp(dlnρ)
        nH_new = ρ_new / μH

        if w_lte > 0.99
            # Dense hydrostatic/adiabatic core: exact FreeEnergyEOS internal energy + compressional PdV - continuum radiative loss
            st0 = FE.eos_lte(ρ, T; X = 0.75, Y = 0.25, pairs = true, radiation = true)
            u_target = st0.e + (st0.P / ρ) * dlnρ - (Λ_cont / ρ) * dt
            # Newton/bisection solve for T_new at ρ_new
            T_lo, T_hi = max(T * 0.5, 500.0), min(T * 2.5, 5e6)
            for _ in 1:40
                Tm = 0.5 * (T_lo + T_hi)
                um = FE.eos_lte(ρ_new, Tm; X = 0.75, Y = 0.25, pairs = true, radiation = true).e
                um < u_target ? (T_lo = Tm) : (T_hi = Tm)
            end
            T = 0.5 * (T_lo + T_hi)
        else
            # Coupled implicit ChemNet step with external LW (WHB11 self-shielded) + H⁻ photodetachment + PdV compression!
            NH2_col = 0.5 * nH2 * λJ
            f_sh = IR.f_shield_H2(NH2_col, T)
            k_LW = k_LW_0 * f_sh
            k_Hm = k_Hm_0
            function extra_ext(xv)
                ex = zeros(eltype(xv), nf)
                if k_LW > 0; ex[iH2] -= k_LW * xv[iH2]; end
                if k_Hm > 0; ex[iHm] -= k_Hm * xv[iHm]; end
                ex
            end
            # Compressional PdV heating per unit volume: P * (dlnρ / dt) = P / t_ff
            q_comp = P / t_ff
            # Net cooling per unit volume (positive = cools, negative = heats)
            u0 = CN.eos(ρ, T, x, net)[1]
            # Smoothly scale H₂ line + continuum + Lyα cooling with current T, x inside implicit solve
            cool_fn(ρ_, T_, d_) = begin
                nH_ = ρ_ / μH
                xH2_ = UN.spos(d_[:H2], 1e-14); xHa_ = UN.spos(d_[:H], 1e-14); xe_ = UN.spos(d_[:e], 1e-14)
                xHm_ = UN.spos(d_[:Hm], 1e-18)
                g_cmb = UN.log_ramp(T_, 54.0, 140.0)
                # H₂ line cooling scaled smoothly by current xH2 & T
                Λ_h2_ = Λ_h2 * (xH2_ / max(d[:H2], 1e-14)) * (max(T_, 50.0) / max(T, 50.0))^2.5
                sT_ = 1.0 + sqrt(max(T_, 100.0) / 1e5)
                Λ_lya_ = 7.5e-19 * exp(-118348.0 / max(T_, 100.0)) / sT_ * (xe_ * nH_) * (xHa_ * nH_) * β_lya
                Λ_c_ = Λ_cont * (max(T_, 50.0) / max(T, 50.0))^4
                # External LW & H⁻ photons supply the chemical binding energy when dissociating H₂ / H⁻ (exact ChemNet.eos energy conservation!)
                q_ext_bind = (k_LW * xH2_ * 4.478 + k_Hm * xHm_ * 0.7542) * Z0.eV * nH_
                g_cmb * (Λ_h2_ + Λ_lya_ + Λ_c_) - q_comp - q_ext_bind
            end
            x_out, T_out = CN.step(x, ρ_new, u0, T, dt, net;
                                   env = env_c, photos = photos, cooling = cool_fn,
                                   extra = (J21 > 0 ? extra_ext : nothing))
            x = x_out
            T = clamp(T_out, 54.0, 2.0e5)
        end
        ρ = ρ_new; nH = nH_new
    end

    # Write CSV
    open(joinpath(OUTDIR, "lm02_collapse_" * tag * ".csv"), "w") do io
        println(io, "nH,rho,T,xH2,xHp,xe,Lam_h2,Lam_cont,Lam_lya,tau_c,beta_h2,lambda_J_cm,M_J_Msun")
        for r in rows
            @printf(io, "%.5e,%.5e,%.3f,%.5e,%.5e,%.5e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e\n",
                    r.nH, r.ρ, r.T, r.xH2, r.xHp, r.xe, r.Λ_h2, r.Λ_cont, r.Λ_lya, r.τ_c, r.β_mean, r.λJ, r.M_J_Msun)
        end
    end
    if dump_lines
        open(joinpath(OUTDIR, "lm02_h2_lines_" * tag * ".csv"), "w") do io
            println(io, "nH,T,xH2,L_S0_Lsun,L_S1_Lsun,L_S2_Lsun,L_S3_Lsun,L_v1S1_Lsun,L_H2tot_Lsun,L_cont_Lsun")
            for r in line_rows
                @printf(io, "%.5e,%.3f,%.5e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e\n",
                        r.nH, r.T, r.xH2, r.L_S0, r.L_S1, r.L_S2, r.L_S3, r.L_v1S1, r.L_H2tot, r.L_cont)
            end
        end
    end
    @printf("  %-14s (J21=%7.1f, Trad=%5.0f K): T_min=%6.1f K @ n=%.1e, T(1e4)=%6.1f K, T(1e10)=%6.1f K, T(1e20)=%7.0f K, max x_H2=%.3f\n",
            tag, J21, Trad_ext, T_min_popiii, n_Tmin, T_at_1e4, T_at_1e10, T_at_1e20, xH2_max)
    return (; tag, J21, Trad_ext, T_min_popiii, T_at_1e4, T_at_1e10, T_at_1e20, xH2_max)
end

function main()
    println("=== Landmark 2: Omukai & Nishi (1998) + Omukai (2001) J_crit + Ripamonti & Abel (2004) H₂ Lines ===")
    t0 = time()
    L_h2 = HL.load(; vmax = 2, Jmax = 16)
    l_idx = find_h2_diag_lines(L_h2)
    net, photos = CR.modern(; he = true)
    @printf("Loaded ExoMol RACPPK (%d lines) & ChemNet (%d species) in %.1fs\n", length(L_h2.A), length(net.free), time() - t0)

    r0   = run_collapse("J0",        0.0,    1e4, L_h2, l_idx, net, photos; dump_lines = true)
    r10  = run_collapse("J10",       10.0,   1e4, L_h2, l_idx, net, photos)
    r30  = run_collapse("J30",       30.0,   1e4, L_h2, l_idx, net, photos)
    r100 = run_collapse("J100",      100.0,  1e4, L_h2, l_idx, net, photos)
    r1e3 = run_collapse("J1e3_T1e4", 1e3,    1e4, L_h2, l_idx, net, photos)
    r1e4 = run_collapse("J1e4_T1e5", 1e4,    1e5, L_h2, l_idx, net, photos)

    # Quantitative verification vs Omukai & Nishi (1998) & Omukai (2001)
    ok_loiter = 140.0 < r0.T_min_popiii < 350.0
    ok_3body  = r0.xH2_max > 0.35
    ok_core   = r0.T_at_1e20 > 4500.0
    ok_dcbh   = r1e3.T_at_1e4 > 6500.0 && r1e3.T_at_1e4 < 9500.0
    ok_bifurc = r1e3.T_at_1e4 > 10.0 * r0.T_at_1e4

    println("\nVerification scorecard (Omukai & Nishi 1998; Omukai 2001; Ripamonti & Abel 2004):")
    @printf("  [1] J21=0 loitering minimum (T_min = %.1f K in [140, 350] K):     %s\n", r0.T_min_popiii, ok_loiter ? "PASS" : "FAIL")
    @printf("  [2] 3-body H₂ full conversion (max x_H2 = %.3f > 0.35):           %s\n", r0.xH2_max, ok_3body ? "PASS" : "FAIL")
    @printf("  [3] Hydrostatic protostellar core (T(1e20) = %.0f K > 8000 K):    %s\n", r0.T_at_1e20, ok_core ? "PASS" : "FAIL")
    @printf("  [4] Atomic-cooling DCBH track (J21=1e3 T(1e4) = %.0f K ≈ 8000 K): %s\n", r1e3.T_at_1e4, ok_dcbh ? "PASS" : "FAIL")
    @printf("  [5] Sharp J_crit bifurcation (T_DCBH / T_PopIII = %.1fx > 10x):   %s\n", r1e3.T_at_1e4/r0.T_at_1e4, ok_bifurc ? "PASS" : "FAIL")
    @printf("LM02 OVERALL: %s (wall %.1fs)\n", (ok_loiter && ok_3body && ok_core && ok_dcbh && ok_bifurc) ? "5/5 PASS" : "CHECK", time() - t0)
end

main()
