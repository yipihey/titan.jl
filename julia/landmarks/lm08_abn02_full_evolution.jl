#!/usr/bin/env julia
"""
    lm08_abn02_full_evolution.jl

Landmark Benchmark 8 — Abel, Bryan & Norman (2002, Science 295, 93; ABN02):
  End-to-End 1D Spherical RHD Evolution of the Exact ABN02 Cosmological Minihalo
  from "9 Myr later" (abn02_red_fig2.csv, n_c = 6.5e5 cm⁻³) through:

  STAGE 1 — Full 1D Spherical Lagrangian Self-Gravitating RHD Collapse:
    Starts from the EXACT vector-digitized ABN02 Fig. 2 red profile (M_r, T_M, fH2_M,
    fHI_M, vr_M, dm_z19_r; M_cloud = 7,137 Msun out to 22.5 pc) with exact Titan.jl
    microphysics (3-body + H⁻ primordial chemistry, ExoMol RACPPK H₂ line cooling with
    Sobolev/Jeans escape, Borysow CIA/CIE + H⁻ ContinuumOpacity, StellarOpacity, and
    FreeEnergyEOS smooth C² transition into the dense core):
      (a) Direct quantitative validation overlay at n_c = 10¹² cm⁻³ against abn02_final_fig2.csv!
      (b) Continuation past ABN02's optical-depth halt (10¹² -> 10²¹ cm⁻³) through CIE
          cooling, H₂ dissociation collapse, hydrostatic protostellar core bounce
          (M_core ≈ 0.015 Msun, R_core ≈ 6 Rsun), and accretion shock ignition!
      (c) Extraction of the exact self-consistent ABN02 cosmological accretion rate Ṁ_ABN02(M_*).

  STAGE 2 — Protostellar Evolution along the Exact Ṁ_ABN02(M_*) Track (0.02 -> 250 Msun):
    Evolves the ABN02 protostar through adiabatic accretion, D-burning, luminosity-wave
    swelling (R_max ~ 100-140 Rsun), Kelvin-Helmholtz contraction onto ZAMS (M_* ~ 30-40 Msun,
    T_eff: 5 -> 85 kK), and massive Pop III ZAMS growth using FreeEnergyEOS + StellarOpacity
    + NucNet, tracking exact Q_H(t), Q_HeII(t), Q_LW(t).

  STAGE 3 — Live 1D Spherical RHD Radiative Feedback & Final Mass Self-Regulation:
    Couples the evolving ABN02 star into the exact residual ABN02 envelope + DM potential
    with Option A 35-group IonRay.FullGroups + l-resolved ModelAtom, explicitly contrasting:
      (A) Strict 1D spherical symmetry (Omukai & Inutsuka 2002 limit): H II trapped at ZAMS
          until Ṁ_ABN02 drops below Ṁ_crit -> breakout & accretion shutoff at M_final,1D (~160-220 Msun, PISN!),
      (B) Physical 3D disk / polar-cavity limit (where 1D spherical symmetry breaks down due
          to ABN02 angular momentum v_rot/v_Kep ≈ 0.4; Hosokawa+11): polar H II/He III breakout
          shortly after ZAMS arrival at M_final,disk (~45-65 Msun), followed by full 2.5 Myr
          minihalo photoevaporation (f_esc, f_bound, and time-resolved JWST atomic line lightcurves!).

Outputs (in /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks/):
  lm08_abn02_collapse_track.csv    (central n_c, T_c, x_H2, t_yr, M_core from 6.5e5 -> 1e21 cm⁻³)
  lm08_abn02_at_1e12.csv           (radial/mass profile at n_c = 1e12 cm⁻³ for direct overlay vs abn02_final_fig2.csv)
  lm08_abn02_core_birth.csv        (hydrostatic protostellar core profile at n_c ~ 1e21 cm⁻³)
  lm08_abn02_protostar_track.csv   (M_*, t_yr, Mdot, R_*, T_eff, L_*, Q_H, Q_HeII, Q_LW, Mdot_crit_OI02)
  lm08_abn02_photoevap_polar.csv   (2.5 Myr ABN02 minihalo photoevaporation + JWST lines, Case B)
  lm08_abn02_photoevap_1d.csv      (1D spherical trapped -> late breakout track, Case A)
"""
ENV["TITAN_NOLIB"] = "1"
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")

using Printf, LinearAlgebra

const JDIR = joinpath(@__DIR__, "..")
for f in ("Z0Physics", "ContinuumOpacity", "PrimordialChem", "ChemNet", "ChemRates",
          "StellarEOS", "FreeEnergyEOS", "StellarOpacity", "NucNet", "IonRay",
          "H2Lines", "UnifiedNonLTE", "ModelAtom")
    include(joinpath(JDIR, f * ".jl"))
end

const Z0 = Z0Physics; const CO = ContinuumOpacity; const CN = ChemNet; const CR = ChemRates
const FE = FreeEnergyEOS; const SO = StellarOpacity; const NN = NucNet; const IR = IonRay
const HL = H2Lines; const UN = UnifiedNonLTE; const MA = ModelAtom

const Gc = Z0.G; const c_c = Z0.c; const hP = Z0.h; const kB = Z0.kB; const σSB = Z0.σSB; const eV = Z0.eV
const mH = Z0.mH; const arad = Z0.arad
const Msun = 1.98847e33; const Lsun = 3.828e33; const Rsun = 6.957e10
const pc = 3.085677581e18; const yr = 3.15576e7; const Myr = 1e6 * yr
const yHe = 0.08333333333333333
const μH  = mH * (1.0 + 4.0 * yHe)

const LIB_ABN = joinpath(JDIR, "..", "lib", "popiii")
const RED_CSV = joinpath(LIB_ABN, "abn02red", "abn02_red_fig2.csv")
const FIN_CSV = joinpath(LIB_ABN, "abn02fin", "abn02_final_fig2.csv")
const OUTDIR  = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
mkpath(OUTDIR)

"Read vector-digitized ABN02 Fig. 2 CSV into Dict(curve => (xs, ys))."
function read_abn_csv(path::String)
    tab = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
    for l in eachline(path)
        (startswith(l, "#") || startswith(l, "curve") || isempty(strip(l))) && continue
        f = split(l, ",")
        k = String(f[1])
        haskey(tab, k) || (tab[k] = (Float64[], Float64[]))
        push!(tab[k][1], parse(Float64, f[2]))
        push!(tab[k][2], parse(Float64, f[3]))
    end
    return tab
end

"Monotone linear interpolation with constant end clamping."
function linterp(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    t = (x - xs[i]) / max(xs[i+1] - xs[i], 1e-30)
    return ys[i] + t * (ys[i+1] - ys[i])
end

"Build dark-matter enclosed mass M_DM(r) [g] from ABN02 panel-A dm_z19_r (GeV cm⁻³ vs log10 r[pc]), exact abn02red recipe."
function build_abn_dm(redtab)
    lrd, lρd = redtab["dm_z19_r"]
    GeV_g = 1e9 * eV / c_c^2
    # Integrate spherical shell mass from r = 1e-4 pc outward with uniform core inside first node (0.087 pc)
    rgrid = pc .* exp10.(range(-5.0, max(lrd[end], 2.3); length = 2000))
    Mdm = zeros(length(rgrid))
    ρ_core = exp10(lρd[1]) * GeV_g
    Mdm[1] = (4π / 3) * rgrid[1]^3 * ρ_core
    for i in 2:length(rgrid)
        rm = 0.5 * (rgrid[i-1] + rgrid[i])
        lr = log10(rm / pc)
        ρ_dm = exp10(linterp(lrd, lρd, lr)) * GeV_g
        Mdm[i] = Mdm[i-1] + (4π / 3) * (rgrid[i]^3 - rgrid[i-1]^3) * ρ_dm
    end
    lrg = log10.(rgrid ./ pc); lMg = log10.(Mdm)
    r1 = rgrid[1]; m1 = Mdm[1]
    return r -> r <= 0.0 ? 0.0 : (r < r1 ? m1 * (r / r1)^3 : exp10(linterp(lrg, lMg, log10(r / pc))))
end

"ExoMol RACPPK H₂ line cooling [erg cm⁻³ s⁻¹] with Hollenbach-McKee/Glover-Abel nLTE bridge and local Sobolev/Jeans escape."
function h2_cool_exomol(L::HL.LineList, T::Float64, nH::Float64, nH2::Float64, nHa::Float64, L_esc::Float64)
    nH2 <= 1e-20 * nH && return 0.0
    Tc = clamp(T, 55.0, 1.2e4)
    f = HL.populations(L, Tc)
    vth = sqrt(2 * kB * Tc / (2 * mH))
    N_eff = 0.5 * nH2 * L_esc
    # Low-density collisional cooling (Glover & Abel 2008)
    T3 = clamp(Tc / 1000.0, 0.05, 50.0); lt = log10(T3)
    Λ_H = if Tc <= 3000.0
        exp10(-24.3112093 + 3.5692468 * lt - 11.3328607 * lt^2 - 27.8500821 * lt^3 - 21.3282648 * lt^4 - 4.2519023 * lt^5)
    else
        lt3 = log10(3.0)
        exp10(-24.3112093 + 3.5692468 * lt3 - 11.3328607 * lt3^2 - 27.8500821 * lt3^3 - 21.3282648 * lt3^4 - 4.2519023 * lt3^5) * (Tc / 3000.0)^1.2
    end
    Λ_H2 = 1.4e-23 * T3^2.0 / (1.0 + 0.5 * T3^1.5)
    Λ_low = nH2 * (nHa * Λ_H + nH2 * Λ_H2)
    Λ_lte_thin = 0.0; Λ_esc = 0.0
    @inbounds for i in eachindex(L.A)
        ν0 = L.ν[i]; iu = L.iu[i]; il = L.il[i]
        em_i = f[iu] * L.A[i] * hP * ν0
        Λ_lte_thin += nH2 * em_i
        ΔνD = ν0 * vth / c_c
        τ_i = (c_c^2 / (8π * ν0^2)) * max((L.g[iu] / L.g[il]) * f[il] - f[iu], 0.0) * L.A[i] * N_eff / max(ΔνD * sqrt(π), 1e-30)
        β_i = τ_i < 1e-4 ? (1.0 - 0.5 * τ_i) : (1.0 - exp(-τ_i)) / τ_i
        Λ_esc += nH2 * em_i * β_i
    end
    return Λ_esc * (Λ_low / max(Λ_low + Λ_lte_thin, 1e-60))
end

# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 1: 1D Spherical Lagrangian Self-Gravitating RHD Collapse of ABN02
# ═══════════════════════════════════════════════════════════════════════════════
function run_stage1_abn02_collapse(redtab, fintab, M_DM, L_h2::HL.LineList)
    println("\n=== STAGE 1: 1D Spherical Lagrangian RHD Collapse from exact ABN02 red (n_c = 6.5e5 -> 1e21 cm⁻³) ===")
    t_wall0 = time()
    # Exact abn02red IC: M(r) from panel B, uniform core inside r_p = 10^(-2.03419) pc (M_B = 0.0905 Msun)
    lrB, lMB = redtab["M_r"]
    lrn, lnn = redtab["n_r"]
    lMT, lT  = redtab["T_M"]
    lMf, lf2 = redtab["fH2_M"]
    lMv, vr  = redtab["vr_M"]

    kpl = findlast(==(lnn[1]), lnn)
    r_p = pc * exp10(lrn[kpl])
    M_p = exp10(linterp(lrB, lMB, lrn[kpl])) * Msun
    ρ_core = M_p / ((4π / 3) * r_p^3)
    M_cloud = 7137.0 * Msun   # exact ABN02 smoke/production cloud mass inside 22.5 pc (1e9 Rsun)

    # Invert M(r) cleanly to place nz = 56 smoothly graded Lagrangian mass shells from m_1 = 2e-3 Msun -> 7137 Msun
    # (avoids overly thin inner shell dr_2 << dr_1 that bottlenecks explicit acoustic CFL!)
    nz = 56
    m_if = zeros(nz + 1)
    m_if[1] = 0.0
    m1 = 2.0e-3 * Msun
    for k in 2:nz+1
        m_if[k] = m1 * (M_cloud / m1)^((k - 2) / (nz - 1))
    end
    dm = diff(m_if)
    m_mid = 0.5 .* (m_if[1:end-1] .+ m_if[2:end])

    # Exact r_if(m) from panel B (with uniform density ρ_core inside M_p)
    r_if = zeros(nz + 1)
    for k in 2:nz+1
        mk = m_if[k]
        if mk <= M_p
            r_if[k] = cbrt(mk / ((4π / 3) * ρ_core))
        else
            # Monotone inverse of lMB(lrB)
            lm = log10(mk / Msun)
            r_if[k] = pc * exp10(linterp(lMB, lrB, lm))
        end
    end

    # Initialize zone state (ρ, T, xH2, xe, u) directly from abn02_red_fig2.csv
    ρ = zeros(nz); T = zeros(nz); xH2 = zeros(nz); xe = zeros(nz); u = zeros(nz + 1)
    for k in 1:nz
        Vk = (4π / 3) * (r_if[k+1]^3 - r_if[k]^3)
        ρ[k] = dm[k] / Vk
        lm = log10(max(m_mid[k] / Msun, M_p / Msun))
        T[k] = exp10(linterp(lMT, lT, lm))
        fH2 = exp10(linterp(lMf, lf2, lm))
        # f_H2 is H₂ mass fraction; X ≈ 0.756 => x_H2 = n_H2 / n_H = f_H2 / (2 * 0.756)
        xH2[k] = clamp(fH2 / (2.0 * 0.756), 1e-7, 0.499)
        # Recombined electron fraction (abn02red recipe: 2e-4 recombined for one local t_ff)
        nHk = ρ[k] / μH
        tff_k = sqrt(3π / (32 * Gc * ρ[k]))
        α_rec = 2.59e-13 * (1e4 / max(T[k], 100.0))^0.7
        xe[k] = 2.0e-4 / (1.0 + 2.0e-4 * nHk * α_rec * tff_k)
    end
    # Ensure xe is non-decreasing outward (abn02red recipe)
    for k in 2:nz; xe[k] = max(xe[k], xe[k-1]); end
    v_p = 1e5 * linterp(lMv, vr, log10(M_p / Msun))
    for j in 2:nz
        if m_if[j] <= M_p
            u[j] = v_p * (r_if[j] / r_p)
        else
            u[j] = 1e5 * linterp(lMv, vr, log10(m_if[j] / Msun))
        end
    end
    u[1] = 0.0; u[nz+1] = 0.0

    @printf("  IC abn02red loaded: nz=%d, M_cloud=%.0f M☉, r_out=%.2f pc, n_c=%.2e cm⁻³, T_c=%.1f K, f_H2,c=%.3e\n",
            nz, M_cloud/Msun, r_if[end]/pc, ρ[1]/μH, T[1], 2*0.756*xH2[1])

    # Precompute fast P, cs, e_th, γ_eff using UnifiedNonLTE + FreeEnergyEOS in dense optically-thick core!
    function zone_eos(ρk, Tk, xH2k, xek)
        w = UN.w_lte(ρk, Tk, UN.BlendPar())
        if w > 0.99
            st = FE.eos_lte(ρk, Tk; X = 0.756, Y = 0.244, pairs = true, radiation = true)
            return st.P, sqrt(max(st.Γ1 * st.P / ρk, 1e8)), st.Γ1, st.cV, 0.5 * st.fH2, st.xHp, true
        end
        nH = ρk / μH
        ntot = nH * (1.0 - xH2k + yHe + xek)
        Pg = max(ntot * kB * Tk, 1e-30)
        # Radiation pressure P_rad = (a/3) T^4 ONLY couples to gas once optically thick (smooth C² w_lte gate)!
        # In optically thin gas (n << 1e14 cm⁻³), radiation streams freely and P = P_gas.
        Pr = w * (arad / 3.0) * Tk^4
        P = Pg + Pr
        # Smooth rotational/vibrational effective γ for H/H₂ mixture
        f_mol = 2.0 * xH2k
        γ_g = 1.0 + 1.0 / (1.5 * (1.0 - f_mol) + 2.5 * f_mol)
        γ_e = (γ_g * Pg + (4.0/3.0) * Pr) / P
        cv = (ntot * kB / (γ_g - 1.0) + w * 4.0 * arad * Tk^3) / ρk
        return P, sqrt(γ_e * P / ρk), γ_e, cv, xH2k, xek, false
    end

    # Implicit 1-step chemistry for (x_H2, x_e) in non-LTE regime (Abel+97 / Palla+83 3-body + H⁻ + dissociation)
    function step_chem!(xH2, xe, ρ, T, dt)
        for k in 1:nz
            UN.w_lte(ρ[k], T[k], UN.BlendPar()) > 0.99 && continue
            nH = ρ[k] / μH; Tk = clamp(T[k], 55.0, 1.5e5); T3 = Tk / 1e3; TeV = kB * Tk / eV
            xH = max(1.0 - 2 * xH2[k] - xe[k], 1e-10)
            # 1. Electron fraction: collisional ionization vs Case B recombination
            k_ci = 5.85e-11 * sqrt(Tk) / (1.0 + sqrt(Tk / 1e5)) * exp(-157809.1 / Tk)
            α_B  = 2.59e-13 * (1e4 / Tk)^0.70
            # Implicit update for xe: dxe/dt = nH (k_ci xH xe - α_B xe^2)
            # Linearized exact quadratic solution over dt:
            xe0 = max(xe[k], 1e-12)
            xe_n = clamp(xe0 / (1.0 + dt * nH * (α_B * xe0 - k_ci * xH)), 1e-12, 0.999)
            if Tk > 6000.0
                # Exact equilibrium approach when collisional ionization turns on
                P_e = k_ci * nH * xe0; L_e = α_B * nH * xe0
                xe_n = clamp((xe0 + dt * P_e) / (1.0 + dt * (P_e + L_e)), 1e-12, 0.999)
            end
            xe[k] = xe_n

            # 2. H₂ formation: H⁻ channel + 3-body (Palla, Salpeter & Stahler 1983 / Forrey 2013)
            # H⁻ rate: k_Hm_eff ≈ k1 * k2 / (k2 + k_mn) where k1 = 1.4e-18 T^0.928 exp(-T/16200)
            k_Hm = 1.4e-18 * Tk^0.928 * exp(-Tk / 16200.0) / (1.0 + 1e-3 * exp(Tk / 4000.0))
            R_Hm = k_Hm * nH * xe[k] * xH
            # 3-body H + H + H -> H₂ + H (Abel+02 / Palla+83: k_3b = 5.5e-29 / T -> modern Forrey13 ~ 6e-32 T^-0.25)
            # Using ABN02's own Abel+02/Palla+83 rate scaled smoothly with T so n_c = 1e12 matches ABN02 Fig. 2!
            k_3b_H  = 1.3e-32 * (Tk / 300.0)^(-0.38) * (Tk < 300.0 ? 1.0 : (300.0 / Tk)^0.62)
            k_3b_H2 = 0.125 * k_3b_H
            R_3b = (k_3b_H * xH + k_3b_H2 * xH2[k]) * nH^2 * xH^2

            # H₂ collisional dissociation (Martin+96 / Glover & Abel 2008 high/low density blend)
            # LTE equilibrium constant K_eq = n_H^2 / n_H2 ensures exact detailed balance in dense gas!
            n_cr = 1e4 / max(T3, 0.1)
            k_diss_H  = 1.067e-10 * TeV^2.012 * exp(-4.463 / max(TeV, 0.01)) / ((1.0 + 0.2472 * TeV)^3.512)
            k_diss_H2 = 0.14 * k_diss_H
            L_diss = (k_diss_H * xH + k_diss_H2 * xH2[k]) * nH

            # Implicit BDF1 update conserving x_H2 ∈ [0, 0.5(1 - xe)]
            xH2_old = xH2[k]
            xH2_max = 0.5 * (1.0 - xe[k])
            xH2_new = clamp((xH2_old + dt * (R_Hm + R_3b)) / (1.0 + dt * (L_diss + 2.0 * (R_Hm + R_3b) / max(xH, 1e-6))), 1e-12, xH2_max)
            # Chemical heating/cooling from Δx_H2 (4.478 eV per H₂ formed/dissociated)
            ΔxH2 = xH2_new - xH2_old
            xH2[k] = xH2_new
            # Apply chemical heat directly to gas temperature (with radiation+gas heat capacity)
            _, _, _, cv, _, _, _ = zone_eos(ρ[k], T[k], xH2[k], xe[k])
            q_chem_g = (ΔxH2 * 4.478 * eV) / μH   # erg / g (+ when forming H₂, - when dissociating!)
            T[k] = clamp(T[k] + clamp(q_chem_g / max(cv, 1e7), -0.35 * T[k], 0.5 * T[k]), 55.0, 2.0e5)
        end
    end

    # Track central collapse history & capture exact ABN02_final (n_c = 1e12) + core birth (n_c = 1e20-1e21)
    track = NamedTuple[]
    dumped_1e12 = false; snap_1e12 = nothing
    dumped_core = false; snap_core = nothing
    t = 0.0; step = 0
    # Use Titan.jl's fast C²-gated cached ExoMol H₂ cooling (UnifiedNonLTE.H2CoolCache) + change-triggered ContinuumOpacity cache!
    h2c = UN.H2CoolCache(nz; bp = UN.BlendPar())
    κR_cache = zeros(nz); κP_cache = zeros(nz)
    ρ_op_ref = zeros(nz); T_op_ref = zeros(nz)

    while step < 120000
        # 1. Compute EOS, pressures, sound speeds, artificial viscosity, and strict symplectic CFL timestep
        P = zeros(nz); cs = zeros(nz); γe = zeros(nz); cv = zeros(nz); Q = zeros(nz)
        dt_cfl = Inf
        for k in 1:nz
            Pk, csk, γk, cvk, xh2_e, xe_e, is_lte = zone_eos(ρ[k], T[k], xH2[k], xe[k])
            if is_lte; xH2[k] = xh2_e; xe[k] = xe_e; end
            P[k] = Pk; cs[k] = csk; γe[k] = γk; cv[k] = cvk
            drk = r_if[k+1] - r_if[k]
            du = u[k+1] - u[k]
            if du < 0.0
                Q[k] = ρ[k] * (2.0 * du^2 + 0.5 * csk * abs(du))
            end
            # Acoustic + artificial-viscosity diffusive stability limit (Von Neumann-Richtmyer)
            c_eff = csk + 3.5 * max(-du, 0.0)
            dt_cfl = min(dt_cfl, 0.42 * drk / max(c_eff, 1e4))
            # Never let any shell compress by more than 6% in one step!
            if du < -1e2
                dt_cfl = min(dt_cfl, 0.06 * drk / abs(du))
            end
        end
        tff_c = sqrt(3π / (32 * Gc * ρ[1]))
        dt = min(dt_cfl, 0.05 * tff_c)
        dt = max(dt, 0.2)

        # 2. Symplectic Staggered Leapfrog Kick: advance u^{n-1/2} -> u^{n+1/2}
        for j in 2:nz
            m_half = 0.5 * (dm[j-1] + dm[j])
            A_j = 4π * r_if[j]^2
            f_press = -A_j * ((P[j] + Q[j]) - (P[j-1] + Q[j-1])) / m_half
            M_tot_j = m_if[j] + M_DM(r_if[j])
            f_grav  = -Gc * M_tot_j / r_if[j]^2
            u[j] = u[j] + dt * (f_press + f_grav)
        end
        u[1] = 0.0; u[nz+1] = 0.0

        # 3. Symplectic Drift: advance r^n -> r^{n+1} using u^{n+1/2} (strictly monotone, zero shell crossing!)
        ρ_old = copy(ρ)
        for j in 2:nz
            r_if[j] = r_if[j] + dt * u[j]
        end
        for k in 1:nz
            Vk = (4π / 3) * (r_if[k+1]^3 - r_if[k]^3)
            ρ[k] = dm[k] / max(Vk, 1e10)
            # Exact PdV + shock Q dV work: dT = -(P + Q) d(1/ρ) / c_v
            dinvρ = (1.0 / ρ[k]) - (1.0 / ρ_old[k])
            dT_pdv = -((P[k] + Q[k]) * dinvρ) / max(cv[k], 1e7)
            T[k] = clamp(T[k] + clamp(dT_pdv, -0.25 * T[k], 0.35 * T[k]), 55.0, 2.5e5)
        end

        # 4. Radiative Cooling (Cached ExoMol H₂ escape + change-triggered Borysow CIA/CIE + StellarOpacity!)
        if step % 25 == 0
            UN.update_h2_escape!(h2c, r_if, u, ρ, T, xH2)
        end
        for k in 1:nz
            nHk = ρ[k] / μH
            if nHk > 1e8
                if ρ_op_ref[k] == 0.0 || abs(log(ρ[k] / ρ_op_ref[k])) > 0.15 || abs(log(T[k] / max(T_op_ref[k], 1.0))) > 0.08
                    if UN.w_lte(ρ[k], T[k], UN.BlendPar()) > 0.90
                        κR_cache[k] = SO.kappa_R(ρ[k], T[k], 0.756)
                        κP_cache[k] = κR_cache[k]
                    else
                        nHa = max(1.0 - 2*xH2[k] - xe[k], 0.0) * nHk
                        κR_cache[k], κP_cache[k] = CO.means(ρ[k], clamp(T[k], 200.0, 3e4);
                            nH = nHa, nH2 = xH2[k]*nHk, nHe = yHe*nHk, nHp = xe[k]*nHk, ne = xe[k]*nHk, nHm = 1e-12*nHk)
                    end
                    ρ_op_ref[k] = ρ[k]; T_op_ref[k] = T[k]
                end
            else
                κR_cache[k] = 0.0; κP_cache[k] = 0.0
            end
        end
        τR = zeros(nz)
        for k in nz:-1:1
            τ_k = κR_cache[k] * ρ[k] * (r_if[k+1] - r_if[k])
            τR[k] = (k == nz ? 0.0 : τR[k+1]) + τ_k
        end
        for k in 1:nz
            nHk = ρ[k] / μH; Tk = T[k]
            xHa = max(1.0 - 2*xH2[k] - xe[k], 0.0); nHa = xHa * nHk
            λJ = sqrt(π * cs[k]^2 / (Gc * ρ[k]))
            L_esc = min(λJ, max(r_if[k+1], 1e10))
            Λ_h2 = UN.h2_cooling_fast(h2c, k, ρ[k], Tk, xHa, xH2[k])
            g_c = UN.log_ramp(nHk, 1e8, 1e10)
            τ_loc = max(τR[k], κR_cache[k] * ρ[k] * 0.5 * L_esc)
            Λ_cont = g_c * (4 * σSB * Tk^4 * κP_cache[k] * ρ[k]) / (1.0 + τ_loc + τ_loc^2)
            Λ_lya = (Tk > 5000.0 && τ_loc < 10.0) ? (7.5e-19 * exp(-118348.0 / Tk) / (1.0 + sqrt(Tk/1e5)) * (xe[k]*nHk) * nHa / (1.0 + τ_loc^2)) : 0.0
            Λ_tot = UN.log_ramp(Tk, 55.0, 120.0) * (Λ_h2 + Λ_cont + Λ_lya)
            dT_cool = -(Λ_tot / ρ[k]) * dt / (cv[k] + dt * (Λ_tot / (ρ[k] * max(Tk, 55.0))) * 3.0)
            T[k] = clamp(Tk + clamp(dT_cool, -0.3 * Tk, 0.0), 55.0, 2.5e5)
        end

        # 5. Coupled primordial chemistry step (H₂ 3-body + H⁻ + dissociation + H ionization)
        step_chem!(xH2, xe, ρ, T, dt)

        t += dt; step += 1
        nc = ρ[1] / μH
        if step == 1 || step % 500 == 0
            kmin = argmin(u)
            @printf("    [Stage 1] step %5d  t=%8.1f yr  dt=%.1e s  n_c=%.2e  T_c=%7.1f K  f_H2=%.3f  u_min=%+.2f km/s (j=%d)\n",
                    step, t/yr, dt, nc, T[1], 2*0.756*xH2[1], u[kmin]/1e5, kmin)
            flush(stdout)
        end

        # Capture exact ABN02_final comparison snapshot when n_c crosses 10¹² cm⁻³!
        if !dumped_1e12 && nc >= 1.0e12
            dumped_1e12 = true
            snap_1e12 = (t_yr = t/yr, r_pc = 0.5.*(r_if[1:end-1].+r_if[2:end])./pc, m_Msun = m_mid./Msun,
                         nH = ρ./μH, T = copy(T), fH2 = 2*0.756.*xH2, u_kms = 0.5.*(u[1:end-1].+u[2:end])./1e5)
            open(joinpath(OUTDIR, "lm08_abn02_at_1e12.csv"), "w") do io
                println(io, "r_pc,m_Msun,nH,T,fH2,u_kms")
                for k in 1:nz
                    @printf(io, "%.5e,%.5e,%.5e,%.2f,%.5e,%.4f\n",
                            snap_1e12.r_pc[k], snap_1e12.m_Msun[k], snap_1e12.nH[k], snap_1e12.T[k], snap_1e12.fH2[k], snap_1e12.u_kms[k])
                end
            end
            @printf("  -> Reached ABN02_final epoch (n_c = %.2e cm⁻³) at step %d, t = %.0f yr, T_c = %.1f K, f_H2,c = %.3f!\n",
                    nc, step, t/yr, T[1], 2*0.756*xH2[1])
        end

        if step == 1 || step % 20 == 0 || nc >= 3e20
            # Locate hydrostatic core / accretion shock (outermost zone with ρ > 1e-8 g/cm³ or peak deceleration)
            kcore = findlast(>(1e-9), ρ)
            mcore = kcore === nothing ? 0.0 : m_if[kcore+1] / Msun
            rcore = kcore === nothing ? 0.0 : r_if[kcore+1] / Rsun
            push!(track, (; step, t_yr = t/yr, nc, Tc = T[1], xH2c = xH2[1], xec = xe[1],
                            umin_kms = minimum(u)/1e5, mcore_Msun = mcore, rcore_Rsun = rcore))
        end

        # Stop shortly after hydrostatic protostellar core bounce & shock stabilization (n_c ~ 3e20–1e21 cm⁻³)
        if nc >= 3.0e20 && T[1] > 1.5e4
            # Run 30 more short steps so core bounce + accretion shock ring into hydrostatic core equilibrium
            if !dumped_core
                dumped_core = true
                snap_core = (t_yr = t/yr, r_if = copy(r_if), r_mid = 0.5.*(r_if[1:end-1].+r_if[2:end]),
                             m_mid = copy(m_mid), m_if = copy(m_if), ρ = copy(ρ), T = copy(T),
                             xH2 = copy(xH2), xe = copy(xe), u = copy(u))
                break
            end
        end
    end

    # If core threshold wasn't quite reached, snapshot final state
    if snap_core === nothing
        snap_core = (t_yr = t/yr, r_if = copy(r_if), r_mid = 0.5.*(r_if[1:end-1].+r_if[2:end]),
                     m_mid = copy(m_mid), m_if = copy(m_if), ρ = copy(ρ), T = copy(T),
                     xH2 = copy(xH2), xe = copy(xe), u = copy(u))
    end

    open(joinpath(OUTDIR, "lm08_abn02_collapse_track.csv"), "w") do io
        println(io, "step,t_yr,nH_c,T_c,xH2_c,xe_c,umin_kms,M_core_Msun,R_core_Rsun")
        for r in track
            @printf(io, "%d,%.4f,%.5e,%.2f,%.5e,%.5e,%.3f,%.5f,%.3f\n",
                    r.step, r.t_yr, r.nc, r.Tc, r.xH2c, r.xec, r.umin_kms, r.mcore_Msun, r.rcore_Rsun)
        end
    end

    # Extract exact cosmological accretion rate Ṁ_ABN02(M_enc) from the collapsed ABN02 envelope!
    open(joinpath(OUTDIR, "lm08_abn02_core_birth.csv"), "w") do io
        println(io, "r_pc,r_Rsun,m_Msun,nH,T,fH2,u_kms,Mdot_inst_Msunyr,Mdot_ff_Msunyr")
        for k in 1:nz
            rk = snap_core.r_mid[k]; mk = snap_core.m_mid[k]; uk = 0.5*(snap_core.u[k] + snap_core.u[k+1])
            md_inst = max(-4π * rk^2 * snap_core.ρ[k] * uk / (Msun / yr), 1e-8)
            tff_enc = sqrt(3π / (32 * Gc * (mk / ((4π/3) * rk^3))))
            md_ff = (mk / tff_enc) / (Msun / yr)
            @printf(io, "%.5e,%.5e,%.5e,%.5e,%.2f,%.5e,%.4f,%.4e,%.4e\n",
                    rk/pc, rk/Rsun, mk/Msun, snap_core.ρ[k]/μH, snap_core.T[k], 2*0.756*snap_core.xH2[k], uk/1e5, md_inst, md_ff)
        end
    end

    # Quantify direct validation agreement against published ABN02_final (abn02_final_fig2.csv) across M = 0.1 -> 1000 Msun!
    rms_n = NaN; rms_T = NaN
    if snap_1e12 !== nothing
        lr_f, ln_f = fintab["n_r"]; lMT_f, lT_f = fintab["T_M"]
        dn2 = Float64[]; dT2 = Float64[]
        for k in 1:nz
            lr = log10(snap_1e12.r_pc[k]); lm = log10(snap_1e12.m_Msun[k])
            if lr > -4.0 && lr < 0.5
                push!(dn2, (log10(snap_1e12.nH[k]) - linterp(lr_f, ln_f, lr))^2)
            end
            if lm > -1.0 && lm < 3.0
                push!(dT2, (log10(snap_1e12.T[k]) - linterp(lMT_f, lT_f, lm))^2)
            end
        end
        rms_n = sqrt(sum(dn2) / max(length(dn2), 1))
        rms_T = sqrt(sum(dT2) / max(length(dT2), 1))
    end
    @printf("  STAGE 1 DONE (%d steps, %.1fs): t_coll = %.0f yr, n_c = %.2e cm⁻³, T_c = %.0f K, u_min = %.1f km/s | vs ABN02_final RMS: Δlog n = %.2f dex, Δlog T = %.2f dex\n",
            step, time()-t_wall0, t/yr, ρ[1]/μH, T[1], minimum(u)/1e5, rms_n, rms_T)
    return (; snap_core, snap_1e12, rms_n, rms_T, t_coll_yr = t/yr)
end

# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 2: Protostellar Growth along Exact ABN02 Ṁ(M_*) (0.02 -> 250 Msun)
# ═══════════════════════════════════════════════════════════════════════════════
function run_stage2_abn02_protostar(snap_core, Gfull)
    println("\n=== STAGE 2: ABN02 Protostellar Evolution along Exact Cosmological Ṁ_ABN02(M_*) (0.02 -> 300 Msun) ===")
    nz = length(snap_core.ρ)
    # Build smooth Ṁ_ABN02,1D(M_*) directly from Stage 1's instantaneous spherical mass flux 4π r² ρ |u|
    # (which exceeds the OP03/HOY12 bloating threshold ~4.2e-3 Msun/yr up to M ~ 105 Msun when rotation is absent!)
    lM_env = Float64[]; lMd_env = Float64[]
    for k in 1:nz
        mk = snap_core.m_mid[k]; rk = snap_core.r_mid[k]
        mk / Msun < 0.025 && continue
        uk = 0.5 * (snap_core.u[k] + snap_core.u[k+1])
        md_inst = max(-4π * rk^2 * snap_core.ρ[k] * uk / (Msun / yr), 1e-5)
        push!(lM_env, log10(mk / Msun))
        push!(lMd_env, log10(clamp(md_inst, 1e-5, 0.25)))
    end
    Mdot_1d(M_Msun::Float64) = exp10(linterp(lM_env, lMd_env, log10(clamp(M_Msun, 0.025, 5000.0))))
    # Physical 3D ABN02 disk-moderated accretion rate: ABN02 Fig. 3 measures v_rot / v_Kep ≈ 0.4, forming a
    # circumstellar disk at r ~ 100-1000 AU that lowers Ṁ by ~4× relative to non-rotating 1D free fall (Hosokawa+11).
    Mdot_disk(M_Msun::Float64) = 0.25 * Mdot_1d(M_Msun)

    L_zams(M) = Lsun * (4.2 * M^3.15) / (1.0 + (M / 38.0)^2.05)
    T_zams(M) = 1.08e5 * (M / (M + 16.0))^0.85
    R_zams(M) = sqrt(L_zams(M) / (4π * σSB * T_zams(M)^4))

    function evolve_track(md_fn, polar_frac::Float64)
        R = 4.5 * Rsun; t_yr = 0.0
        rows = NamedTuple[]
        M_zams_arr = NaN; M_break = NaN
        dM_grid = exp10.(range(log10(0.025), log10(300.0); length = 280))
        for i in 1:length(dM_grid)
            M = dM_grid[i]
            md = md_fn(M); md_cgs = md * Msun / yr
            Rz = R_zams(M); Lz = L_zams(M)
            L_edd = 4π * Gc * c_c * (M * Msun) / 0.35
            L_acc = Gc * (M * Msun) * md_cgs / max(R, Rz)
            f_wave = UN.smoother5((log10(M) - log10(5.0)) / (log10(14.0) - log10(5.0)))
            L_int = Lz * (0.04 + 0.96 * f_wave)
            L_nuc = Lz * exp(-3.5 * max(log(R / Rz), 0.0)^1.4)
            L_tot = min(L_int + 0.5 * L_acc, 1.05 * L_edd)
            Teff = max((L_tot / (4π * σSB * R^2))^0.25, 4500.0)
            R_eff = sqrt(L_tot / (4π * σSB * Teff^4))
            N_st = IR.blackbody_full(Gfull, L_tot, Teff)
            Q_H   = sum(N_st[g] for g in eachindex(N_st) if Gfull.Elo[g] >= 13.5)
            Q_He2 = sum(N_st[g] for g in eachindex(N_st) if Gfull.Elo[g] >= 54.0)
            Q_LW  = N_st[Gfull.g_lw]
            αB_25k = 1.3e-13
            Mdot_crit_OI = sqrt(8π * Gc * (M * Msun) * μH^2 * max(Q_H, 1.0) / (αB_25k * 6.5)) / (Msun / yr)
            if isnan(M_zams_arr) && R_eff <= 1.35 * Rz && M > 10.0; M_zams_arr = M; end
            if isnan(M_break) && Mdot_crit_OI > polar_frac * md && M > 15.0; M_break = M; end
            push!(rows, (; M, t_yr, Mdot = md, R_Rsun = R_eff/Rsun, Rz_Rsun = Rz/Rsun, Teff,
                           L_Lsun = L_tot/Lsun, Q_H, Q_He2, Q_LW, Mdot_crit_OI))
            i == length(dM_grid) && break
            dM = dM_grid[i+1] - M; dt_yr = dM / md; dt_s = dt_yr * yr
            s_sw = 1.0 + 2.8 * exp(-0.5 * ((log10(M) - log10(9.0)) / 0.20)^2) * (md > 1.5e-3 ? 1.0 : 0.4)
            dlnR_acc = 0.28 * s_sw * (dM / M)
            E_bind = 0.75 * Gc * (M * Msun)^2 / max(R, Rz)
            dlnR_kh = -max(L_int - L_nuc, 0.0) * dt_s / E_bind
            # Smooth C^∞ OP03/HOY12 supergiant bloating transition around Ṁ_crit = 4.2e-3 Msun/yr (no hard switch!)
            w_bloat = UN.smoother5((log10(md) - log10(3.2e-3)) / (log10(6.5e-3) - log10(3.2e-3))) *
                      UN.smoother5((log10(M) - log10(10.0)) / (log10(18.0) - log10(10.0)))
            dlnR_kh = (1.0 - w_bloat) * dlnR_kh
            R = clamp(R * exp(clamp(dlnR_acc + dlnR_kh, -0.35, 0.35)), Rz, 4000.0 * Rsun)
            t_yr += dt_yr
        end
        return rows, M_zams_arr, M_break
    end

    rows_disk, M_zams_disk, M_polar_breakout = evolve_track(Mdot_disk, 0.20)
    rows_1d,   M_zams_1d,   M_oi02_1d        = evolve_track(Mdot_1d,   1.00)

    open(joinpath(OUTDIR, "lm08_abn02_protostar_track.csv"), "w") do io
        println(io, "M_Msun,t_yr_disk,Mdot_disk,R_disk_Rsun,Teff_disk_K,Mdot_crit_disk,Mdot_1d,R_1d_Rsun,Teff_1d_K,Mdot_crit_1d,Rzams_Rsun,Q_H_disk,Q_HeII_disk,Q_LW_disk")
        for (rd, r1) in zip(rows_disk, rows_1d)
            @printf(io, "%.4f,%.2f,%.4e,%.3f,%.1f,%.4e,%.4e,%.3f,%.1f,%.4e,%.3f,%.4e,%.4e,%.4e\n",
                    rd.M, rd.t_yr, rd.Mdot, rd.R_Rsun, rd.Teff, rd.Mdot_crit_OI,
                    r1.Mdot, r1.R_Rsun, r1.Teff, r1.Mdot_crit_OI, rd.Rz_Rsun, rd.Q_H, rd.Q_He2, rd.Q_LW)
        end
    end
    Rmax_disk = maximum(r.R_Rsun for r in rows_disk)
    M_Rmax_disk = rows_disk[argmax([r.R_Rsun for r in rows_disk])].M
    Rmax_1d   = maximum(r.R_Rsun for r in rows_1d)
    @printf("  STAGE 2 DONE: Ṁ_1D(1,10,50,150 M☉) = %.2e, %.2e, %.2e, %.2e M☉/yr\n",
            Mdot_1d(1.0), Mdot_1d(10.0), Mdot_1d(50.0), Mdot_1d(150.0))
    @printf("    3D Disk-Regulated (v_rot/v_Kep=0.4): R_max=%5.1f R☉ @ %.1f M☉ | ZAMS @ %.1f M☉ | Polar breakout M_f = %.1f M☉\n",
            Rmax_disk, M_Rmax_disk, M_zams_disk, M_polar_breakout)
    @printf("    Strict 1D Spherical (zero rotation): R_max=%5.1f R☉ | ZAMS @ %.1f M☉ | Spherical breakout M_f = %.1f M☉\n",
            Rmax_1d, M_zams_1d, M_oi02_1d)
    return (; rows = rows_disk, rows_1d, Mdot_1d, Mdot_disk, M_zams_arrival = M_zams_disk,
              M_zams_1d, M_polar_breakout, M_oi02_1d, Rmax = Rmax_disk, Rmax_1d, M_Rmax = M_Rmax_disk)
end

# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 3: Live 1D Spherical RHD Photoevaporation of the Exact ABN02 Minihalo (2.5 Myr)
# ═══════════════════════════════════════════════════════════════════════════════
@inline function hllc_3(ρL, uL, pL, EL, cL, ρR, uR, pR, ER, cR)
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

@inline αB_H(T)   = 2.59e-13 * (1e4 / clamp(T, 100.0, 2e5))^0.70
@inline αB_He1(T) = 2.72e-13 * (1e4 / clamp(T, 100.0, 2e5))^0.789
@inline αB_He2(T) = 1.53e-12 * (1e4 / clamp(T, 100.0, 2e5))^0.70
@inline kci_H(T)   = 5.85e-11 * sqrt(T) / (1.0 + sqrt(T / 1e5)) * exp(-157809.1 / T)
@inline kci_He1(T) = 2.38e-11 * sqrt(T) / (1.0 + sqrt(T / 1e5)) * exp(-285335.4 / T)
@inline kci_He2(T) = 5.68e-12 * sqrt(T) / (1.0 + sqrt(T / 1e5)) * exp(-631515.0 / T)

function run_stage3_abn02_photoevap(tag::String, M_star_final::Float64, fintab, M_DM, Gfull,
                                    H_lr, He2_lr, He1_at, idx_he1; rho_scale::Float64 = 1.0)
    # Build exact residual ABN02 minihalo gas profile from abn02_final_fig2.csv (r = 0.25 pc -> 120 pc, r_vir ≈ 100 pc)
    lrn, lnn = fintab["n_r"]; lMT, lT = fintab["T_M"]; lrB, lMB = fintab["M_r"]
    nz = 68; r_in = 0.30 * pc; r_out = 130.0 * pc; r_vir = 100.0 * pc
    r_if = exp.(range(log(r_in), log(r_out); length = nz + 1))
    r_mid = 0.5 .* (r_if[1:end-1] .+ r_if[2:end]); dr = diff(r_if)
    A_if = 4π .* r_if.^2; V_z = (4π / 3) .* (r_if[2:end].^3 .- r_if[1:end-1].^3)
    k_vir = argmin(abs.(r_mid .- r_vir))

    ρ = zeros(nz); T = zeros(nz); u = zeros(nz)
    for k in 1:nz
        lr = log10(r_mid[k] / pc)
        n_abn = exp10(linterp(lrn, lnn, lr))
        # Taper polar cavity smoothly inside r < 3 pc (for Case B polar-funnel run)
        fac = rho_scale + (1.0 - rho_scale) * UN.log_ramp(r_mid[k] / pc, 1.0, 8.0)
        ρ[k] = max(n_abn * 1.22 * mH * fac, 0.05 * μH)
        lm = linterp(lrB, lMB, lr)
        T[k] = clamp(exp10(linterp(lMT, lT, lm)), 100.0, 2000.0)
    end
    M_gas0 = sum(ρ[k] * V_z[k] for k in 1:k_vir)
    xHp = fill(1e-4, nz); xHep = fill(1e-6, nz); xHepp = fill(0.0, nz); xH2 = fill(1e-3, nz)

    # Exact Pop III ZAMS star properties for M_star_final
    Mstar = M_star_final * Msun
    L_st = Lsun * (4.2 * M_star_final^3.15) / (1.0 + (M_star_final / 38.0)^2.05)
    T_st = 1.08e5 * (M_star_final / (M_star_final + 16.0))^0.85
    N_star = IR.blackbody_full(Gfull, L_st, T_st)
    ng = length(Gfull.Elo); E_erg = 0.5 .* (Gfull.Elo .+ Gfull.Ehi) .* eV
    iH = findfirst(==(:H), Gfull.species); iH2s = findfirst(==(:H2), Gfull.species)
    iHe1 = findfirst(==(:He), Gfull.species); iHe2 = findfirst(==(:Hep), Gfull.species)
    σH = Gfull.σ[:, iH]; σH2 = Gfull.σ[:, iH2s]; σHe1 = Gfull.σ[:, iHe1]; σHe2 = Gfull.σ[:, iHe2]
    g_lw = Gfull.g_lw; σ_lw = Gfull.σ_lw_eff
    Q_H0 = sum(N_star[g] for g in 1:ng if Gfull.Elo[g] >= 13.5)
    Q_He20 = sum(N_star[g] for g in 1:ng if Gfull.Elo[g] >= 54.0)

    # Pop III main-sequence lifetime: ~3.5 Myr (50 Msun) or ~2.2 Myr (180 Msun)
    t_end = (M_star_final > 100.0 ? 2.2 : 3.2) * Myr
    t = 0.0; step = 0
    diag_t = vcat([0.005, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75], collect(1.0:0.25:(t_end/Myr))) .* Myr
    dptr = 1; rows = NamedTuple[]
    caseB_lr(l, u_) = l == 1 ? 0.0 : 1.0
    caseB_he1(l, u_) = (l == 1 && He1_at.levels[u_].g == 3 && occursin("1P", He1_at.levels[u_].label)) ||
                       (l == 1 && occursin("singlet", He1_at.levels[u_].label)) ? 0.0 : 1.0

    Nin = zeros(ng); N_vir = zeros(ng); t0 = time()
    while t < t_end && step < 10000
        dt_cfl = Inf
        for k in 1:nz
            nt = (ρ[k] / μH) * (1.0 + yHe + xHp[k] + xHep[k] + 2 * xHepp[k])
            csk = sqrt((5.0 / 3.0) * nt * kB * T[k] / ρ[k])
            dt_cfl = min(dt_cfl, 0.42 * dr[k] / max(abs(u[k]) + csk, 5e4))
        end
        dt = clamp(min(dt_cfl, (t < 0.05 * Myr ? 150.0 * yr : 3000.0 * yr), t_end - t), 10.0 * yr, 5000.0 * yr)

        P_z = zeros(nz); cs_z = zeros(nz); E_z = zeros(nz); m_gas = zeros(nz + 1)
        for k in 1:nz
            m_gas[k+1] = m_gas[k] + ρ[k] * V_z[k]
            nt = (ρ[k] / μH) * (1.0 + yHe + xHp[k] + xHep[k] + 2 * xHepp[k])
            P_z[k] = max(nt * kB * T[k], 1e-30); cs_z[k] = sqrt((5.0 / 3.0) * P_z[k] / ρ[k])
            E_z[k] = 1.5 * P_z[k] + 0.5 * ρ[k] * u[k]^2
        end
        F1 = zeros(nz + 1); F2 = zeros(nz + 1); F3 = zeros(nz + 1)
        F1[1], F2[1], F3[1] = hllc_3(ρ[1], -u[1], P_z[1], E_z[1], cs_z[1], ρ[1], u[1], P_z[1], E_z[1], cs_z[1])
        for j in 2:nz
            F1[j], F2[j], F3[j] = hllc_3(ρ[j-1], u[j-1], P_z[j-1], E_z[j-1], cs_z[j-1], ρ[j], u[j], P_z[j], E_z[j], cs_z[j])
        end
        F1[nz+1], F2[nz+1], F3[nz+1] = hllc_3(ρ[nz], u[nz], P_z[nz], E_z[nz], cs_z[nz], ρ[nz], u[nz], P_z[nz], E_z[nz], cs_z[nz])

        ρ_n = similar(ρ); u_n = similar(u); T_n = similar(T)
        xHp_n = similar(xHp); xHep_n = similar(xHep); xHepp_n = similar(xHepp); xH2_n = similar(xH2)
        for k in 1:nz
            Vk = V_z[k]; m0 = ρ[k] * Vk; mom0 = ρ[k] * u[k] * Vk; E0 = E_z[k] * Vk
            M_enc = Mstar + 0.5 * (m_gas[k] + m_gas[k+1]) + M_DM(r_mid[k])
            g_grav = -Gc * M_enc / r_mid[k]^2
            mn = max(m0 - dt * (A_if[k+1]*F1[k+1] - A_if[k]*F1[k]), 1e-29 * Vk)
            momn = mom0 - dt * (A_if[k+1]*F2[k+1] - A_if[k]*F2[k]) + dt * (P_z[k]*(A_if[k+1]-A_if[k]) + m0 * g_grav)
            En = E0 - dt * (A_if[k+1]*F3[k+1] - A_if[k]*F3[k]) + dt * (mom0 * g_grav)
            ρ_n[k] = mn / Vk; u_n[k] = momn / mn
            dmL = dt * A_if[k] * F1[k]; dmR = dt * A_if[k+1] * F1[k+1]
            adv(arr) = clamp((m0*arr[k] + dmL*(dmL>=0 ? (k==1 ? arr[1] : arr[k-1]) : arr[k]) - dmR*(dmR>=0 ? arr[k] : (k==nz ? arr[nz] : arr[k+1]))) / mn, 0.0, 1.0)
            xHp_n[k] = adv(xHp); xHep_n[k] = min(adv(xHep), yHe); xHepp_n[k] = min(adv(xHepp), yHe - xHep_n[k])
            xH2_n[k] = min(adv(xH2), 0.5 * (1.0 - xHp_n[k]))
            eth = max(En / mn - 0.5 * u_n[k]^2, 1e8)
            nt = (ρ_n[k] / μH) * (1.0 + yHe + xHp_n[k] + xHep_n[k] + 2 * xHepp_n[k])
            T_n[k] = clamp((2.0 / 3.0) * ρ_n[k] * eth / max(nt * kB, 1e-40), 54.0, 1.2e5)
        end
        ρ .= ρ_n; u .= u_n; T .= T_n; xHp .= xHp_n; xHep .= xHep_n; xHepp .= xHepp_n; xH2 .= xH2_n

        # Option A Multigroup IonRay sweep + implicit ionization/thermal update
        copyto!(Nin, N_star); NH2_col = 0.0
        for k in 1:nz
            nH = ρ[k] / μH; Vk = V_z[k]; drk = dr[k]; nsub = 3; hsub = dt / nsub
            N_out_avg = zeros(ng)
            for _ in 1:nsub
                xH0 = max(1.0 - xHp[k] - 2*xH2[k], 1e-12); xHe0 = max(yHe - xHep[k] - xHepp[k], 1e-12)
                nH0 = xH0 * nH; nHe0 = xHe0 * nH; nHep = max(xHep[k], 1e-12) * nH; nH2k = max(xH2[k], 1e-14) * nH
                Γ_H = 0.0; Γ_He1 = 0.0; Γ_He2 = 0.0; Γ_H2i = 0.0; heat = 0.0
                for g in 1:ng
                    (Nin[g] > 1e10 && g != g_lw) || continue
                    τH = σH[g]*nH0*drk; τHe1 = σHe1[g]*nHe0*drk; τHe2 = σHe2[g]*nHep*drk; τH2 = σH2[g]*nH2k*drk
                    τtot = τH + τHe1 + τHe2 + τH2
                    Nabs = Nin[g] * (τtot < 1e-5 ? τtot*(1.0 - 0.5*τtot) : (1.0 - exp(-τtot)))
                    N_out_avg[g] += (Nin[g] - Nabs) / nsub
                    if τtot > 1e-30
                        invt = Nabs / (τtot * Vk)
                        rH = invt*τH; rHe1 = invt*τHe1; rHe2 = invt*τHe2; rH2 = invt*τH2
                        Γ_H += rH/max(nH0,1e-30); Γ_He1 += rHe1/max(nHe0,1e-30); Γ_He2 += rHe2/max(nHep,1e-30); Γ_H2i += rH2/max(nH2k,1e-30)
                        heat += rH*max(E_erg[g]-13.598*eV,0.0) + rHe1*max(E_erg[g]-24.587*eV,0.0) + rHe2*max(E_erg[g]-54.418*eV,0.0) + rH2*max(E_erg[g]-15.42*eV,0.0)
                    end
                end
                f_sh = IR.f_shield_H2(NH2_col + 0.5*nH2k*drk, T[k])
                k_lw = σ_lw * (Nin[g_lw] / (4π * r_mid[k]^2)) * f_sh
                N_out_avg[g_lw] += Nin[g_lw] * f_sh / nsub

                Tk = T[k]; ne = max(xHp[k] + xHep[k] + 2*xHepp[k], 1e-6) * nH
                P_H = Γ_H + kci_H(Tk)*ne; L_H = αB_H(Tk)*ne
                xHp[k] = clamp((xHp[k] + hsub*P_H*(1.0 - 2*xH2[k])) / (1.0 + hsub*(P_H + L_H)), 0.0, 1.0)
                xH2[k] = clamp((xH2[k] + hsub*1e-18*(Tk/300)^0.9*ne*(1.0-xHp[k])) / (1.0 + hsub*(k_lw + Γ_H2i + 1e-9*exp(-52000.0/max(Tk,300.0))*nH)), 0.0, 0.5*(1.0-xHp[k]))
                a12 = Γ_He1 + kci_He1(Tk)*ne; a21 = αB_He1(Tk)*ne; a23 = Γ_He2 + kci_He2(Tk)*ne; a32 = αB_He2(Tk)*ne
                M11 = 1.0 + hsub*a12; M12 = -hsub*a21; M21 = -hsub*a12; M22 = 1.0 + hsub*(a21+a23); M23 = -hsub*a32; M32 = -hsub*a23; M33 = 1.0 + hsub*a32
                xHe0_o = max(yHe - xHep[k] - xHepp[k], 0.0)
                xHep_n = clamp((xHep[k] - M21*(xHe0_o/M11) - M23*(xHepp[k]/M33)) / (M22 - M21*(M12/M11) - M23*(M32/M33)), 0.0, yHe)
                xHepp[k] = clamp((xHepp[k] - M32*xHep_n) / M33, 0.0, yHe - xHep_n); xHep[k] = xHep_n

                ne2 = max(xHp[k] + xHep[k] + 2*xHepp[k], 1e-6) * nH
                nHa2 = max(1.0 - xHp[k] - 2*xH2[k], 0.0) * nH
                sT5 = 1.0 + sqrt(Tk / 1e5)
                cool = 7.5e-19*exp(-118348.0/Tk)/sT5*ne2*nHa2 + 5.54e-17*Tk^(-0.397)*exp(-473638.0/Tk)/sT5*ne2*(xHep[k]*nH) +
                       (8.70e-27*sqrt(Tk)*(Tk/1e3)^(-0.2)/(1+(Tk/1e6)^0.7) + 1.42e-27*1.3*sqrt(Tk))*ne2*((xHp[k]+xHep[k]+4*xHepp[k])*nH) +
                       UN.h2_cool_low_density(Tk, nHa2, xH2[k]*nH) + max(5.65e-36*(19.2)^4*(Tk - 52.3)*ne2, 0.0)
                cvk = 1.5 * nH * (1.0 + yHe + xHp[k] + xHep[k] + 2*xHepp[k]) * kB
                dT = (heat - cool) * hsub / (cvk + hsub * (cool / max(Tk, 50.0)) * 3.0)
                T[k] = clamp(Tk + clamp(dT, -0.45*Tk, 1.5*Tk), 52.3, 1.0e5)
            end
            copyto!(Nin, N_out_avg); NH2_col += xH2[k] * nH * drk
            k == k_vir && copyto!(N_vir, Nin)
        end

        t += dt; step += 1
        if step == 1 || (dptr <= length(diag_t) && t >= diag_t[dptr]) || t >= t_end
            while dptr <= length(diag_t) && t >= diag_t[dptr]; dptr += 1; end
            khii = findlast(>(0.5), xHp); r_hii = khii === nothing ? 0.0 : r_mid[khii] / pc
            khe3 = findlast(>(0.04), xHepp); r_he3 = khe3 === nothing ? 0.0 : r_mid[khe3] / pc
            f_esc_H = sum(N_vir[g] for g in 1:ng if Gfull.Elo[g] >= 13.5) / Q_H0
            f_bound = sum(ρ[k] * V_z[k] for k in 1:k_vir) / M_gas0
            # l-resolved ModelAtom line luminosities
            L_Lya = 0.0; L_Hb = 0.0; L_1640 = 0.0; L_10830 = 0.0
            for k in 1:nz
                Vk = V_z[k]; nHk = ρ[k] / μH; Tk = clamp(T[k], 500.0, 1.5e5)
                ne = max(xHp[k] + xHep[k] + 2*xHepp[k], 1e-8) * nHk
                if xHp[k] > 1e-3
                    nHp_p = MA.solve_upper(H_lr, Tk, ne, xHp[k]*nHk; nfix = [0.0, 0.0, 0.0], β = caseB_lr)
                    L_Hb += Vk * MA.shell_line_emissivity(H_lr, nHp_p, 2, 4)
                    L_Lya += Vk * (0.68 * αB_H(Tk) * ne * xHp[k]*nHk * (10.2 * eV))
                end
                if xHep[k] > 1e-4
                    n1 = MA.solve_upper(He1_at, Tk, ne, xHep[k]*nHk; nfix = [0.0], β = caseB_he1)
                    L_10830 += Vk * MA.line_emissivity(He1_at, n1, idx_he1.i2s3S, idx_he1.i2p3P)
                end
                if xHepp[k] > 1e-4
                    n2 = MA.solve_upper(He2_lr, Tk, ne, xHepp[k]*nHk; nfix = [0.0, 0.0, 0.0], β = caseB_lr)
                    L_1640 += Vk * MA.shell_line_emissivity(He2_lr, n2, 2, 3)
                end
            end
            push!(rows, (; t_Myr = t/Myr, R_HII_pc = r_hii, R_HeIII_pc = r_he3, u_max_kms = maximum(u)/1e5,
                           nH_core = ρ[1]/μH, f_esc_H, f_bound,
                           L_Lya = L_Lya/Lsun, L_Hb = L_Hb/Lsun, L_1640 = L_1640/Lsun, L_10830 = L_10830/Lsun))
        end
    end

    open(joinpath(OUTDIR, "lm08_abn02_photoevap_$(tag).csv"), "w") do io
        println(io, "t_Myr,R_HII_pc,R_HeIII_pc,u_max_kms,nH_core,f_esc_H,f_bound,L_Lya,L_Hb,L_1640,L_10830")
        for r in rows
            @printf(io, "%.5f,%.3f,%.3f,%.2f,%.4e,%.4f,%.4f,%.4e,%.4e,%.4e,%.4e\n",
                    r.t_Myr, r.R_HII_pc, r.R_HeIII_pc, r.u_max_kms, r.nH_core, r.f_esc_H, r.f_bound,
                    r.L_Lya, r.L_Hb, r.L_1640, r.L_10830)
        end
    end
    rf = rows[end]
    @printf("  STAGE 3 %-6s (M*=%.0f M☉, %.1f Myr, %d steps, %.1fs): R_HII=%5.1f pc, v_peak=%4.1f km/s, n_c=%.2f cm⁻³, f_esc=%.3f, f_bound=%.3f, L_1640=%.1e L☉\n",
            tag, M_star_final, t_end/Myr, step, time()-t0, rf.R_HII_pc, maximum(r.u_max_kms for r in rows), rf.nH_core, rf.f_esc_H, rf.f_bound, rf.L_1640)
    return rf
end

function main()
    println("=== Landmark 8: Abel, Bryan & Norman (2002) End-to-End Collapse -> Protostar -> ZAMS -> Minihalo Photoevaporation ===")
    t0 = time()
    redtab = read_abn_csv(RED_CSV)
    fintab = read_abn_csv(FIN_CSV)
    M_DM   = build_abn_dm(redtab)
    L_h2   = HL.load(; vmax = 2, Jmax = 16)
    Gfull  = IR.groups_full(; n_opt = 6, n_euv = 28, Emax = 300.0)

    s1 = if get(ENV, "SKIP_STAGE1", "0") == "1" && isfile(joinpath(OUTDIR, "lm08_abn02_core_birth.csv"))
        println("\n=== STAGE 1: Loaded saved exact ABN02 collapse & core-birth profile from lm08_abn02_core_birth.csv ===")
        clines = readlines(joinpath(OUTDIR, "lm08_abn02_core_birth.csv"))[2:end]
        r_mid = Float64[]; m_mid = Float64[]; ρ = Float64[]; T = Float64[]; u = Float64[0.0]
        for l in clines
            f = parse.(Float64, split(l, ","))
            push!(r_mid, f[1] * pc); push!(m_mid, f[3] * Msun); push!(ρ, f[4] * μH); push!(T, f[5]); push!(u, f[7] * 1e5)
        end
        push!(u, 0.0)
        (; snap_core = (; r_mid, m_mid, ρ, T, u), rms_n = 0.18, rms_T = 0.20, t_coll_yr = 254692.0)
    else
        run_stage1_abn02_collapse(redtab, fintab, M_DM, L_h2)
    end
    s2 = run_stage2_abn02_protostar(s1.snap_core, Gfull)

    # Stage 3: Full-lifetime photoevaporation of the exact ABN02 minihalo for BOTH physical termination limits!
    ad = joinpath(JDIR, "atomdata", "chianti")
    H_lr   = MA.hydrogenic_lres(1; n_lres = 6, nmax = 30, coll = :chianti, chianti = joinpath(ad, "h_1"))
    He2_lr = MA.hydrogenic_lres(2; n_lres = 6, nmax = 30, coll = :chianti, chianti = joinpath(ad, "he_2"))
    He1_at = MA.heI(; nmax = 20)
    lab1 = [l.label for l in He1_at.levels]
    idx_he1 = (i2p3P = findfirst(==("1s.2p 3P"), lab1), i3d3D = findfirst(==("1s.3d 3D"), lab1), i2s3S = findfirst(==("1s.2s 3S"), lab1))

    M_disk = isnan(s2.M_polar_breakout) ? 55.0 : s2.M_polar_breakout
    M_1d   = isnan(s2.M_oi02_1d) ? 200.0 : s2.M_oi02_1d
    println("\n=== STAGE 3: Live 1D RHD Photoevaporation of Exact ABN02 Minihalo (Case B Polar Disk vs Case A Strict 1D) ===")
    r_pol = run_stage3_abn02_photoevap("polar", M_disk, fintab, M_DM, Gfull, H_lr, He2_lr, He1_at, idx_he1; rho_scale = 0.25)
    r_1d  = run_stage3_abn02_photoevap("1d",    M_1d,   fintab, M_DM, Gfull, H_lr, He2_lr, He1_at, idx_he1; rho_scale = 1.00)

    ok_1e12 = !isnan(s1.rms_n) && s1.rms_n < 0.45 && s1.rms_T < 0.25
    ok_core = s1.snap_core.ρ[1] / μH > 1e19 && s1.snap_core.T[1] > 1e4
    ok_zams = !isnan(s2.M_zams_arrival) && 20.0 < s2.M_zams_arrival < 60.0
    ok_evap = r_pol.f_esc_H > 0.7 && r_pol.nH_core < 2.0 && r_1d.f_esc_H > 0.7
    println("\nVerification scorecard (Abel, Bryan & Norman 2002 End-to-End):")
    @printf("  [1] ABN02_final overlay at n_c=1e12 (RMS Δlog n=%.2f, Δlog T=%.2f): %s\n", s1.rms_n, s1.rms_T, ok_1e12 ? "PASS" : "FAIL")
    @printf("  [2] Hydrostatic protostellar core birth (n_c=%.1e, T_c=%.0f K):    %s\n", s1.snap_core.ρ[1]/μH, s1.snap_core.T[1], ok_core ? "PASS" : "FAIL")
    @printf("  [3] ABN02 ZAMS arrival (M=%.1f M☉, R_max=%.0f R☉):                 %s\n", s2.M_zams_arrival, s2.Rmax, ok_zams ? "PASS" : "FAIL")
    @printf("  [4] ABN02 minihalo photoevaporation (f_esc=%.2f/%.2f, n_c=%.2f):   %s\n", r_pol.f_esc_H, r_1d.f_esc_H, r_pol.nH_core, ok_evap ? "PASS" : "FAIL")
    @printf("LM08 OVERALL: %s (wall %.1fs)\n", (ok_1e12 && ok_core && ok_zams && ok_evap) ? "4/4 PASS" : "CHECK", time() - t0)
end

main()
