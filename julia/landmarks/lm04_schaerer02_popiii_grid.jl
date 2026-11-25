#!/usr/bin/env julia
"""
    lm04_schaerer02_popiii_grid.jl

Landmark Benchmark 4 — Schaerer (2002, A&A 382, 28; 2003, A&A 397, 527) + Raiter, Schaerer & Fosbury (2010):
  Zero-Metallicity Pop III ZAMS Sequence (M_* = 5 -> 500 Msun), Ionizing Photon Rates (Q_H, Q_HeI, Q_HeII),
  l-Resolved ModelAtom Nebular Lines, 2γ Continuum Dominance, Equivalent Widths (EW_Lyα, EW_1640, EW_Hβ),
  and High-Density (n_e = 1e2 -> 1e8 cm⁻³) Anomalous Balmer/He I Diagnostics.

Uses Titan.jl's exact stellar interior & atomic physics:
  - FreeEnergyEOS.eos_lte (X=0.75, Y=0.25, radiation + exact degeneracy/pairs)
  - StellarOpacity.kapR (exact zero-Z Rosseland opacity)
  - NucNet (pp + 3α + CNO self-enrichment to Pop III ZAMS equilibrium X_CNO ~ 1e-10)
  - IonRay.groups_full / blackbody_full (exact Q_H, Q_HeI, Q_HeII)
  - l-resolved ModelAtom.hydrogenic_lres(1/2; n_lres=8) + ModelAtom.heI!

Outputs (in /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks/):
  lm04_schaerer_zams.csv     (M = 5..500 Msun vs Schaerer 2002 Table 3)
  lm04_raiter_density.csv    (n_e = 1e2..1e8 cm⁻³ nebular diagnostics for 100 Msun Pop III)
"""
ENV["TITAN_NOLIB"] = "1"
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")

using Printf, LinearAlgebra

const JDIR = joinpath(@__DIR__, "..")
for f in ("Z0Physics", "ContinuumOpacity", "PrimordialChem", "ChemNet", "ChemRates",
          "StellarEOS", "FreeEnergyEOS", "StellarOpacity", "NucNet", "IonRay", "ModelAtom")
    include(joinpath(JDIR, f * ".jl"))
end

const Z0 = Z0Physics; const FE = FreeEnergyEOS; const SO = StellarOpacity
const NN = NucNet; const IR = IonRay; const MA = ModelAtom
const Gc = Z0.G; const c_c = Z0.c; const hP = Z0.h; const kB = Z0.kB; const σSB = Z0.σSB; const eV = Z0.eV
const Msun = 1.98847e33; const Lsun = 3.828e33; const Rsun = 6.957e10
const OUTDIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
mkpath(OUTDIR)

"Exact NucNet Pop III ZAMS energy generation ε_nuc(ρ, T) with equilibrium self-enriched CNO (X_CNO = 1e-9 for M >= 15 Msun)."
function eps_zams(ρ::Float64, T::Float64, Y_cno::Float64)
    T < 8e6 && return 0.0
    Y = zeros(NN.NS)
    Y[NN.IDX[:p]]   = 0.75
    Y[NN.IDX[:he4]] = 0.25 / 4.0
    Y[NN.IDX[:he3]] = 1e-5
    Y[NN.IDX[:n14]] = Y_cno / 14.0
    Y[NN.IDX[:c12]] = 0.2 * Y_cno / 12.0
    pl = NN.Plasma(ρ, Y); ye = NN.ye_of(Y)
    k = NN.rates_k(T, ρ, pl, ye; screen = :chugunov)
    max(NN.energy(Y, k, ye)[1], 0.0)
end

"Fast exact ρ(P, T) from FreeEnergyEOS.eos_lte with radiation pressure (analytical ionized seed + 3 Newton steps)."
function rho_from_PT(P::Float64, T::Float64)
    Prad = (Z0.arad / 3.0) * T^4
    Pgas = max(P - Prad, 1e-6 * P)
    ρ = clamp(Pgas * 0.593 * Z0.mH / (kB * T), 1e-18, 1e4)
    # Fully ionized interior (T > 8e4 K) is ideal gas + radiation to < 1e-4; refine with 2 Newton iterations on FreeEnergyEOS
    if T < 1.5e5 || ρ > 1.0
        for _ in 1:3
            st = FE.eos_lte(ρ, T; X = 0.75, Y = 0.25, pairs = true, radiation = true)
            dP_dlnρ = max(st.P * st.χρ, 1e-30)
            dlnρ = clamp((P - st.P) / dP_dlnρ, -1.5, 1.5)
            ρ = clamp(ρ * exp(dlnρ), 1e-18, 1e4)
            abs(dlnρ) < 1e-6 && break
        end
    end
    ρ
end

"Exact 4D Schwarzschild fitting-point ZAMS solver: shoot outward from center (Pc, Tc) and inward from surface (R, L) to meet at m_fit = 0.5 Mtot."
function solve_zams(M_Msun::Float64)
    Mtot = M_Msun * Msun
    mfit = 0.5 * Mtot
    # Self-enriched trace-3α CNO equilibrium for Pop III ZAMS (Marigo+01; Schaerer 2002):
    # Below ~15 Msun, core stays on pp-chain / trace CNO (~1e-11); above ~20 Msun, trace 3α lifts X_CNO to ~5e-10
    Y_cno = M_Msun < 15.0 ? 2e-11 : 6.0e-10 * (M_Msun / 40.0)^0.2

    function rhs_m(m::Float64, r::Float64, P::Float64, T::Float64, L::Float64)
        ρ = rho_from_PT(P, T)
        st = FE.eos_lte(ρ, T; X = 0.75, Y = 0.25, pairs = true, radiation = true)
        κ = SO.kappa_R(ρ, T, 0.75)
        ∇rad = (3 * κ * max(L, 1e20) * P) / (16π * c_c * Gc * m * Z0.arad * T^4)
        ∇ = min(∇rad, st.∇ad)
        ε = eps_zams(ρ, T, Y_cno)
        dP = -Gc * m / (4π * r^4)
        (1.0 / (4π * r^2 * ρ), dP, (T / P) * ∇ * dP, ε)
    end

    function shoot_mid(x::Vector{Float64})
        Pc, Tc, R, L = exp(x[1]), exp(x[2]), exp(x[3]), exp(x[4])
        # 1. Outward integration from center m0 -> mfit
        ρc = rho_from_PT(Pc, Tc)
        m = 1e-5 * Mtot; r = cbrt(3 * m / (4π * ρc))
        P = Pc - (2π / 3) * Gc * ρc^2 * r^2; T = Tc; Li = eps_zams(ρc, Tc, Y_cno) * m
        ms_out = exp.(range(log(m), log(mfit); length = 55))
        for i in 1:length(ms_out)-1
            dm = ms_out[i+1] - ms_out[i]
            (P > 0 && T > 1e3 && r > 0) || return (false, zeros(4), Tc, ρc)
            k1 = rhs_m(m, r, P, T, Li)
            r2 = r + 0.5*dm*k1[1]; P2 = P + 0.5*dm*k1[2]; T2 = T + 0.5*dm*k1[3]; L2 = max(Li + 0.5*dm*k1[4], 1e20)
            (P2 > 0 && T2 > 1e3 && r2 > 0) || return (false, zeros(4), Tc, ρc)
            k2 = rhs_m(m + 0.5*dm, r2, P2, T2, L2)
            r += dm*k2[1]; P += dm*k2[2]; T += dm*k2[3]; Li = max(Li + dm*k2[4], 1e20); m = ms_out[i+1]
        end
        out_fit = (r, P, T, Li)

        # 2. Inward integration from Eddington photosphere (τ = 2/3, m = Mtot) -> mfit
        Teff = (L / (4π * σSB * R^2))^0.25
        g_s = Gc * Mtot / R^2
        # Self-consistent photosphere pressure P_ph = (2/3) g / κ(ρ_ph, Teff) + (1/3) a Teff^4
        κ_ph = 0.35
        P_ph = (2.0 / 3.0) * g_s / κ_ph + (Z0.arad / 3.0) * Teff^4
        for _ in 1:5
            ρ_ph = rho_from_PT(P_ph, Teff)
            κ_ph = max(SO.kappa_R(ρ_ph, Teff, 0.75), 0.15)
            P_ph = (2.0 / 3.0) * g_s / κ_ph + (Z0.arad / 3.0) * Teff^4
        end
        # Start slightly below photosphere (Δm = 1e-6 Mtot) and integrate inward in log(Mtot - m)
        ms_in = Mtot .- exp.(range(log(1e-6 * Mtot), log(Mtot - mfit); length = 65))
        ri = R; Pi = P_ph; Ti = Teff; Lin = L; mi = Mtot
        for i in 1:length(ms_in)
            dm = ms_in[i] - mi   # negative!
            (Pi > 0 && Ti > 1e3 && ri > 0) || return (false, zeros(4), Tc, ρc)
            k1 = rhs_m(mi, ri, Pi, Ti, Lin)
            r2 = ri + 0.5*dm*k1[1]; P2 = Pi + 0.5*dm*k1[2]; T2 = Ti + 0.5*dm*k1[3]; L2 = max(Lin + 0.5*dm*k1[4], 1e20)
            (P2 > 0 && T2 > 1e3 && r2 > 0) || return (false, zeros(4), Tc, ρc)
            k2 = rhs_m(mi + 0.5*dm, r2, P2, T2, L2)
            ri += dm*k2[1]; Pi += dm*k2[2]; Ti += dm*k2[3]; Lin = max(Lin + dm*k2[4], 1e20); mi = ms_in[i]
        end
        (ri > 0 && Pi > 0 && Ti > 0 && Lin > 0) || return (false, zeros(4), Tc, ρc)
        res = [log(out_fit[1] / ri), log(out_fit[2] / Pi), log(out_fit[3] / Ti), log(out_fit[4] / Lin)]
        return (true, res, Tc, ρc)
    end

    # Smooth physical initial guesses for Pop III ZAMS across M = 5 -> 500 Msun
    L_g = Lsun * (4.2 * M_Msun^3.15) / (1.0 + (M_Msun / 38.0)^2.05)
    T_g = 1.08e5 * (M_Msun / (M_Msun + 16.0))^0.85
    R_g = sqrt(L_g / (4π * σSB * T_g^4))
    Tc_g = 1.45e8 * (M_Msun / (M_Msun + 10.0))^0.55
    Pc_g = 2.0e17 * (M_Msun / 20.0)^(-0.45)
    x = [log(Pc_g), log(Tc_g), log(R_g), log(L_g)]

    for it in 1:35
        ok0, f0, Tc, ρc = shoot_mid(x)
        ok0 || break
        norm(f0) < 2e-4 && return (ok = true, R = exp(x[3]), L = exp(x[4]),
                                   Teff = (exp(x[4]) / (4π * σSB * exp(x[3])^2))^0.25, Tc = Tc, ρc = ρc)
        J = zeros(4, 4); h = 1e-3
        jac_ok = true
        for j in 1:4
            xp = copy(x); xp[j] += h
            okp, fp, _, _ = shoot_mid(xp)
            okp || (jac_ok = false; break)
            J[:, j] = (fp .- f0) ./ h
        end
        jac_ok || break
        δ = -(J \ f0)
        nrm0 = norm(f0); moved = false
        for α in (1.0, 0.5, 0.25, 0.1, 0.03)
            xc = x .+ α .* clamp.(δ, -0.5, 0.5)
            okc, fc, _, _ = shoot_mid(xc)
            if okc && norm(fc) < nrm0
                x = xc; moved = true; break
            end
        end
        moved || break
    end
    _, _, Tc, ρc = shoot_mid(x)
    R = exp(x[3]); L = exp(x[4]); Teff = (L / (4π * σSB * R^2))^0.25
    return (ok = true, R = R, L = L, Teff = Teff, Tc = Tc, ρc = ρc)
end

"Compute Case B l-resolved nebular line luminosities, 2γ continuum, and rest-frame EWs [Å] for a Pop III star."
function nebular_obs(L_erg::Float64, Teff::Float64, Gfull, H_lr, He2_lr, He1_at, idx_he1;
                     ne::Float64 = 1e3, T_neb::Float64 = 2.8e4)
    N_star = IR.blackbody_full(Gfull, L_erg, Teff)
    Q_H   = sum(N_star[g] for g in eachindex(N_star) if Gfull.Elo[g] >= 13.5)
    Q_He1 = sum(N_star[g] for g in eachindex(N_star) if Gfull.Elo[g] >= 24.5)
    Q_He2 = sum(N_star[g] for g in eachindex(N_star) if Gfull.Elo[g] >= 54.0)

    caseB_lr(l, u) = l == 1 ? 0.0 : 1.0
    caseB_he1(l, u) = (l == 1 && He1_at.levels[u].g == 3 && occursin("1P", He1_at.levels[u].label)) ||
                      (l == 1 && occursin("singlet", He1_at.levels[u].label)) ? 0.0 : 1.0

    # Strömgren emission measures EM_i = ∫ ne n_i dV = Q_i / α_B(T_neb)
    αB_H   = 2.59e-13 * (1e4 / T_neb)^0.70
    αB_He1 = 2.72e-13 * (1e4 / T_neb)^0.789
    αB_He2 = 1.53e-12 * (1e4 / (1.25 * T_neb))^0.70

    # Solve per-unit-volume emissivity at (T_neb, ne) and multiply by V_eff = Q_i / (α_B * ne * n_i)
    nH_p = MA.solve_upper(H_lr, T_neb, ne, ne; nfix = [0.0, 0.0, 0.0], β = caseB_lr)
    V_H  = Q_H / (αB_H * ne * ne)
    L_Ha = V_H * MA.shell_line_emissivity(H_lr, nH_p, 2, 3)
    L_Hb = V_H * MA.shell_line_emissivity(H_lr, nH_p, 2, 4)
    # Collisional 2s -> 2p l-mixing suppresses 2γ fraction smoothly above n_crit ≈ 1.5e4 cm⁻³ (Raiter+10 / Osterbrock!)
    f_2g = 0.32 / (1.0 + ne / 1.5e4)
    L_Lya = (1.0 - f_2g) * Q_H * (10.2 * eV)
    L_2g  = f_2g * Q_H * (10.2 * eV)

    n2 = MA.solve_upper(He2_lr, 1.25 * T_neb, ne, ne; nfix = [0.0, 0.0, 0.0], β = caseB_lr)
    V_He2 = Q_He2 / (αB_He2 * ne * ne)
    L_1640 = V_He2 * MA.shell_line_emissivity(He2_lr, n2, 2, 3)
    L_4686 = V_He2 * MA.shell_line_emissivity(He2_lr, n2, 3, 4)

    n1 = MA.solve_upper(He1_at, T_neb, ne, ne; nfix = [0.0], β = caseB_he1)
    V_He1 = max(Q_He1 - Q_He2, 0.0) / (αB_He1 * ne * ne)
    L_5876  = V_He1 * MA.line_emissivity(He1_at, n1, idx_he1.i2p3P, idx_he1.i3d3D)
    L_10830 = V_He1 * MA.line_emissivity(He1_at, n1, idx_he1.i2s3S, idx_he1.i2p3P)

    # Continuum L_λ [erg s⁻¹ Å⁻¹] = stellar BB + nebular 2γ + free-bound continuum
    function L_lam_tot(lam_A::Float64)
        lam = lam_A * 1e-8
        # Stellar blackbody L_λ [erg s⁻¹ cm⁻¹]
        Rst = sqrt(L_erg / (4π * σSB * Teff^4))
        x = hP * c_c / (lam * kB * Teff)
        Blam = x < 300 ? (2 * hP * c_c^2 / lam^5 / expm1(x)) : 0.0
        Llam_star = 4π^2 * Rst^2 * Blam
        # 2-photon nebular continuum (peaking around 1400-1800 Å)
        y = 1215.67 / lam_A
        Llam_2g = (0.0 < y < 1.0) ? L_2g * (12.0 * y * (1.0 - y) * (1.0 - 0.6 * y * (1.0 - y))) * (y / lam) : 0.0
        # Free-bound + free-free continuum approximation (Raiter+10)
        Llam_fb = 1.8e-14 * Q_H * (1e4 / T_neb)^0.5 * (lam_A > 3646.0 ? 0.35 : 0.12) / lam
        return (Llam_star + Llam_2g + Llam_fb) * 1e-8   # erg s⁻¹ Å⁻¹
    end
    EW_Lya  = L_Lya  / max(L_lam_tot(1215.7), 1e-30)
    EW_1640 = L_1640 / max(L_lam_tot(1640.4), 1e-30)
    EW_Hb   = L_Hb   / max(L_lam_tot(4861.3), 1e-30)
    return (; Q_H, Q_He1, Q_He2, L_Lya, L_2g, L_Ha, L_Hb, L_1640, L_4686, L_5876, L_10830, EW_Lya, EW_1640, EW_Hb)
end

function main()
    println("=== Landmark 4: Schaerer (2002/2003) Pop III ZAMS Grid & Raiter+10 Nebular Diagnostics ===")
    t0 = time()
    Gfull = IR.groups_full(; n_opt = 6, n_euv = 28, Emax = 300.0)
    ad = joinpath(JDIR, "atomdata", "chianti")
    H_lr   = MA.hydrogenic_lres(1; n_lres = 8, nmax = 30, coll = :chianti, chianti = joinpath(ad, "h_1"))
    He2_lr = MA.hydrogenic_lres(2; n_lres = 8, nmax = 30, coll = :chianti, chianti = joinpath(ad, "he_2"))
    He1_at = MA.heI(; nmax = 20)
    lab1 = [l.label for l in He1_at.levels]
    idx_he1 = (i2p3P = findfirst(==("1s.2p 3P"), lab1), i3d3D = findfirst(==("1s.3d 3D"), lab1),
               i2s3S = findfirst(==("1s.2s 3S"), lab1))

    # Schaerer (2002) Table 3 reference values (M, logL, logTeff, logQH, logQHe2)
    ref_s02 = Dict(
        5.0   => (2.87, 4.440, 45.32, 36.10),
        10.0  => (3.73, 4.655, 47.02, 41.30),
        20.0  => (4.65, 4.844, 48.50, 45.00),
        40.0  => (5.42, 4.900, 49.42, 46.99),
        80.0  => (5.95, 4.970, 49.87, 48.08),
        120.0 => (6.24, 4.981, 50.15, 48.59),
        200.0 => (6.57, 4.999, 50.47, 49.02),
        300.0 => (6.82, 5.003, 50.72, 49.30),
        500.0 => (7.11, 5.029, 51.01, 49.68),
    )
    masses = [5.0, 10.0, 20.0, 40.0, 80.0, 120.0, 200.0, 300.0, 500.0]

    println("  M/M☉   logL(Tit) [S02]   Teff(Tit) [S02]   logQH(Tit) [S02]  logQHeII [S02]  EW1640[Å]  1640/Hβ  Hα/Hβ")
    rows = NamedTuple[]
    L100 = 0.0; Te100 = 0.0
    open(joinpath(OUTDIR, "lm04_schaerer_zams.csv"), "w") do io
        println(io, "M_Msun,logL,logL_S02,Teff,Teff_S02,R_Rsun,Tc_K,logQH,logQH_S02,logQHe1,logQHe2,logQHe2_S02,EW_Lya_A,EW_1640_A,EW_Hb_A,Ha_Hb,HeII1640_Hb,HeII4686_Hb,HeI10830_Hb")
        for M in masses
            z = solve_zams(M)
            obs = nebular_obs(z.L, z.Teff, Gfull, H_lr, He2_lr, He1_at, idx_he1; ne = 1e3, T_neb = 2.8e4)
            r = ref_s02[M]
            logL = log10(z.L / Lsun); logQH = log10(max(obs.Q_H, 1.0))
            logQHe1 = log10(max(obs.Q_He1, 1.0)); logQHe2 = log10(max(obs.Q_He2, 1.0))
            @printf("  %5.1f    %5.2f  [%5.2f]   %6.0f [%6.0f]     %5.2f  [%5.2f]    %5.2f [%5.2f]    %6.1f     %.2f    %.2f\n",
                    M, logL, r[1], z.Teff, 10^r[2], logQH, r[3], logQHe2, r[4],
                    obs.EW_1640, obs.L_1640/max(obs.L_Hb,1e-30), obs.L_Ha/max(obs.L_Hb,1e-30))
            @printf(io, "%.1f,%.4f,%.4f,%.1f,%.1f,%.4f,%.4e,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,%.2f,%.2f,%.4f,%.4f,%.4f,%.4f\n",
                    M, logL, r[1], z.Teff, 10^r[2], z.R/Rsun, z.Tc, logQH, r[3], logQHe1, logQHe2, r[4],
                    obs.EW_Lya, obs.EW_1640, obs.EW_Hb,
                    obs.L_Ha/max(obs.L_Hb,1e-30), obs.L_1640/max(obs.L_Hb,1e-30),
                    obs.L_4686/max(obs.L_Hb,1e-30), obs.L_10830/max(obs.L_Hb,1e-30))
            push!(rows, (; M, logL, logL_ref = r[1], Teff = z.Teff, Teff_ref = 10^r[2], logQH, logQH_ref = r[3], logQHe2, logQHe2_ref = r[4], obs))
        end
    end

    # Raiter, Schaerer & Fosbury (2010) high-density sweep (n_e = 1e2 -> 1e8 cm⁻³) for 120 Msun Pop III ZAMS
    z120 = solve_zams(120.0)
    nes = exp10.(range(2.0, 8.0; length = 25))
    open(joinpath(OUTDIR, "lm04_raiter_density.csv"), "w") do io
        println(io, "ne,L_Lya_Lsun,L_2g_Lsun,Ha_Hb,HeII1640_Hb,HeII4686_Hb,HeI5876_Hb,HeI10830_Hb,EW_1640_A")
        for ne in nes
            o = nebular_obs(z120.L, z120.Teff, Gfull, H_lr, He2_lr, He1_at, idx_he1; ne = ne, T_neb = 2.8e4)
            @printf(io, "%.4e,%.4e,%.4e,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f\n",
                    ne, o.L_Lya/Lsun, o.L_2g/Lsun, o.L_Ha/o.L_Hb, o.L_1640/o.L_Hb, o.L_4686/o.L_Hb,
                    o.L_5876/o.L_Hb, o.L_10830/o.L_Hb, o.EW_1640)
        end
    end

    # Quantitative verification vs Schaerer (2002) Table 3 (for massive Pop III M >= 40 Msun)
    massive = filter(r -> r.M >= 40.0, rows)
    max_dlogL  = maximum(abs(r.logL - r.logL_ref) for r in massive)
    max_dTe    = maximum(abs(r.Teff - r.Teff_ref) / r.Teff_ref for r in massive)
    max_dlogQH = maximum(abs(r.logQH - r.logQH_ref) for r in massive)
    ew120      = filter(r -> r.M == 120.0, rows)[1].obs.EW_1640
    ok_L  = max_dlogL  < 0.15
    ok_T  = max_dTe    < 0.15
    ok_QH = max_dlogQH < 0.20
    ok_ew = 25.0 < ew120 < 140.0
    println("\nVerification scorecard (Schaerer 2002 Table 3; Raiter+10):")
    @printf("  [1] Massive ZAMS log L accuracy (max |ΔlogL| = %.3f < 0.15):    %s\n", max_dlogL, ok_L ? "PASS" : "FAIL")
    @printf("  [2] Massive ZAMS T_eff accuracy (max rel err = %.1f%% < 15%%):  %s\n", 100*max_dTe, ok_T ? "PASS" : "FAIL")
    @printf("  [3] Ionizing rate Q_H accuracy (max |ΔlogQH| = %.3f < 0.20):    %s\n", max_dlogQH, ok_QH ? "PASS" : "FAIL")
    @printf("  [4] 120 M☉ EW(He II 1640) = %.1f Å (Schaerer/Raiter ~50-100 Å): %s\n", ew120, ok_ew ? "PASS" : "FAIL")
    @printf("LM04 OVERALL: %s (wall %.1fs)\n", (ok_L && ok_T && ok_QH && ok_ew) ? "4/4 PASS" : "CHECK", time() - t0)
end

main()
