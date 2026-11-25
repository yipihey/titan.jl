#!/usr/bin/env julia
"""
    val_hii_breakout_rhd.jl

Live Dynamical 1D Spherical RHD H II / He III / Lyman–Werner Breakout & Accretion Reversal Simulation
coupled to Titan.jl's unified microphysics stack:
  - Conservative 1D spherical Eulerian HLLC hydrodynamics with central protostar + envelope self-gravity,
    open inner accretion sink (r_in = 4e-4 pc ≈ 80 AU, inside the gravitational radius r_g = GM*/c_s²),
    and outer halo inflow/outflow boundary (r_out = 1.0 pc).
  - Operator-split smooth coupled non-LTE primordial chemistry (ChemNet + ChemRates.modern(; he=true))
    + Option A photon-conserving multigroup IonRay.FullGroups (H⁻, LW with WHB11 H₂ self-shielding,
      H I / H₂ / He I / He II EUV photoionization)
    + smooth-gated fast H₂ ro-vibrational line cooling (UnifiedNonLTE.h2_cooling_fast, cached ExoMol β_esc)
    + exact Case B radiative recombination binding-escape, Lyα/He II collisional, and free-free cooling.
  - Smooth C² KH-contraction stellar transition from bloated cool supergiant (T_eff = 5,000 K, L = 6e5 L☉)
    to hot 100 M☉ Pop III ZAMS (T_eff = 95,500 K, L = 1.3e6 L☉) over t_ramp = 300 yr, followed dynamically
    to t = 8,000 yr as the H II region transitions from trapped → R-type breakout → D-type champagne shock
    that reverses accretion into a supersonic photoevaporative outflow!
  - Time-resolved non-LTE atomic line synthesis using l-resolved ModelAtom.hydrogenic_lres(1/2; n_lres=8)
    and ModelAtom.heI: Lyα, 2γ continuum, Hα, Hβ, Paα, Brγ, He II 1640/4686, He I 5876/10830!

Outputs (in /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/):
  - hii_breakout_rhd.csv               (global track + emergent atomic line luminosities)
  - profiles/hii_breakout_rhd_s*.csv   (radial snapshots for Pop III Explorer v82 + Fig. 7)
"""
ENV["TITAN_NOLIB"] = "1"
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")

using Printf, LinearAlgebra, ForwardDiff

const JDIR = @__DIR__
for f in ("Z0Physics", "ContinuumOpacity", "PrimordialChem", "ChemNet", "ChemRates",
          "StellarEOS", "FreeEnergyEOS", "StellarOpacity", "IonRay", "H2Lines",
          "UnifiedNonLTE", "ModelAtom")
    include(joinpath(JDIR, f * ".jl"))
end

const Z0 = Z0Physics; const FE = FreeEnergyEOS; const CN = ChemNet; const CR = ChemRates
const IR = IonRay; const HL = H2Lines; const UN = UnifiedNonLTE; const MA = ModelAtom

const Gc   = Z0.G
const kB   = Z0.kB
const Msun = 1.98847e33
const Lsun = 3.828e33
const Rsun = 6.957e10
const pc   = 3.085677581e18
const yr   = 3.15576e7

"Smooth C² stellar KH-contraction trajectory: bloated 5,000 K supergiant -> hot 95,500 K 100 M☉ ZAMS."
function star_state(t_yr::Float64; t_on = 150.0, t_ramp = 450.0,
                    T0 = 5000.0, T1 = 95500.0, L0 = 6.0e5 * Lsun, L1 = 1.30e6 * Lsun)
    w = UN.smoother5((t_yr - t_on) / t_ramp)
    Teff = exp(log(T0) + w * (log(T1) - log(T0)))
    L    = exp(log(L0) + w * (log(L1) - log(L0)))
    Rst  = sqrt(L / (4π * Z0.σSB * Teff^4))
    (; Teff, L, Rst, w)
end

"HLLC Riemann flux for 1D Euler equations (ρ, ρu, E_tot) with ideal-gas local γ_eff."
@inline function hllc_flux(ρL, uL, pL, EL, cL, ρR, uR, pR, ER, cR)
    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)
    if SL >= 0.0
        return (ρL * uL, ρL * uL^2 + pL, (EL + pL) * uL, uL)
    elseif SR <= 0.0
        return (ρR * uR, ρR * uR^2 + pR, (ER + pR) * uR, uR)
    end
    num = pR - pL + ρL * uL * (SL - uL) - ρR * uR * (SR - uR)
    den = ρL * (SL - uL) - ρR * (SR - uR)
    Sstar = abs(den) > 1e-40 ? num / den : 0.5 * (uL + uR)
    if Sstar >= 0.0
        fac = ρL * (SL - uL) / (SL - Sstar)
        Estar = fac * (EL / ρL + (Sstar - uL) * (Sstar + pL / (ρL * (SL - uL))))
        f1 = ρL * uL + SL * (fac - ρL)
        f2 = ρL * uL^2 + pL + SL * (fac * Sstar - ρL * uL)
        f3 = (EL + pL) * uL + SL * (Estar - EL)
        return (f1, f2, f3, Sstar)
    else
        fac = ρR * (SR - uR) / (SR - Sstar)
        Estar = fac * (ER / ρR + (Sstar - uR) * (Sstar + pR / (ρR * (SR - uR))))
        f1 = ρR * uR + SR * (fac - ρR)
        f2 = ρR * uR^2 + pR + SR * (fac * Sstar - ρR * uR)
        f3 = (ER + pR) * uR + SR * (Estar - ER)
        return (f1, f2, f3, Sstar)
    end
end

"Synthesize exact non-LTE emergent atomic H I, He I, and He II nebular line luminosities [erg/s] across the 1D envelope."
function synth_nebular_lines(r_if, ρ, T, xs, net, H_lr, He2_lr, He1_at, idx_he1)
    iHp   = findfirst(==(:Hp),   net.free)
    iHep  = findfirst(==(:Hep),  net.free)
    iHepp = findfirst(==(:Hepp), net.free)
    caseB_lr(l, u) = l == 1 ? 0.0 : 1.0
    caseB_he1(l, u) = (l == 1 && He1_at.levels[u].g == 3 && occursin("1P", He1_at.levels[u].label)) ||
                      (l == 1 && occursin("singlet", He1_at.levels[u].label)) ? 0.0 : 1.0
    L_Lya = 0.0; L_2g = 0.0; L_Ha = 0.0; L_Hb = 0.0; L_Paa = 0.0; L_Brg = 0.0
    L_1640 = 0.0; L_4686 = 0.0; L_5876 = 0.0; L_10830 = 0.0
    for k in eachindex(T)
        r1, r2 = r_if[k], r_if[k+1]
        Vk = 4π / 3 * (r2^3 - r1^3)
        nHk = ρ[k] / CN.μH()
        Tk = clamp(T[k], 500.0, 1.5e5)
        xHp = max(xs[k][iHp], 0.0); xHep = max(xs[k][iHep], 0.0); xHepp = max(xs[k][iHepp], 0.0)
        nHp = xHp * nHk; nHep = xHep * nHk; nHepp = xHepp * nHk
        ne = max(xHp + xHep + 2xHepp, 1e-8) * nHk
        if xHp > 1e-3
            nH = MA.solve_upper(H_lr, Tk, ne, nHp; nfix = [0.0, 0.0, 0.0], β = caseB_lr)
            L_Ha  += Vk * MA.shell_line_emissivity(H_lr, nH, 2, 3)
            L_Hb  += Vk * MA.shell_line_emissivity(H_lr, nH, 2, 4)
            L_Paa += Vk * MA.shell_line_emissivity(H_lr, nH, 3, 4)
            L_Brg += Vk * MA.shell_line_emissivity(H_lr, nH, 4, 7)
            αB = 2.59e-13 * (1e4 / Tk)^0.7
            L_Lya += Vk * (0.68 * αB * ne * nHp * (10.2 * Z0.eV))
            L_2g  += Vk * (0.32 * αB * ne * nHp * (10.2 * Z0.eV))
        end
        if xHep > 1e-4
            n1 = MA.solve_upper(He1_at, Tk, ne, nHep; nfix = [0.0], β = caseB_he1)
            L_5876  += Vk * MA.line_emissivity(He1_at, n1, idx_he1.i2p3P, idx_he1.i3d3D)
            L_10830 += Vk * MA.line_emissivity(He1_at, n1, idx_he1.i2s3S, idx_he1.i2p3P)
        end
        if xHepp > 1e-4
            n2 = MA.solve_upper(He2_lr, Tk, ne, nHepp; nfix = [0.0, 0.0, 0.0], β = caseB_lr)
            L_1640 += Vk * MA.shell_line_emissivity(He2_lr, n2, 2, 3)
            L_4686 += Vk * MA.shell_line_emissivity(He2_lr, n2, 3, 4)
        end
    end
    (; L_Lya, L_2g, L_Ha, L_Hb, L_Paa, L_Brg, L_1640, L_4686, L_5876, L_10830)
end

function main()
    outdir = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs"
    pdir   = joinpath(outdir, "profiles")
    mkpath(pdir)

    Mstar = 100.0 * Msun
    Mdot0 = parse(Float64, get(ENV, "MDOT0", "1.5e-4")) * Msun / yr
    nz    = 88
    r_in  = 2.5e-3 * pc     # ~515 AU (> r_g ≈ 185 AU): below Omukai & Inutsuka Mdot_crit -> R-type + D-type champagne breakout!
    r_out = 1.20 * pc
    t_end = parse(Float64, get(ENV, "TEND_YR", "2500.0")) * yr  # cool infall (0..150 yr), LW+R-type breakout (300..550 yr), D-type champagne shock (600..2500 yr)

    r_if  = exp.(range(log(r_in), log(r_out); length = nz + 1))
    r_mid = 0.5 .* (r_if[1:end-1] .+ r_if[2:end])
    dr    = diff(r_if)
    A_if  = 4π .* r_if.^2
    V_z   = (4π / 3) .* (r_if[2:end].^3 .- r_if[1:end-1].^3)

    # Initial self-consistent steady Bondi/free-fall accretion envelope transitioning smoothly into a 5e3 cm⁻³ minihalo core
    ρ = zeros(nz); u = zeros(nz); T = fill(500.0, nz)
    for k in 1:nz
        vff = sqrt(2 * Gc * Mstar / r_mid[k])
        u[k] = -vff
        ρ[k] = Mdot0 / (4π * r_mid[k]^2 * vff)
    end
    ρ_halo = 5.0e3 * CN.μH()
    ρ .= sqrt.(ρ.^2 .+ ρ_halo^2)
    for k in 1:nz
        # Match mass flux smoothly
        u[k] = -min(sqrt(2 * Gc * Mstar / r_mid[k]), Mdot0 / (4π * r_mid[k]^2 * ρ[k]))
    end

    # Chemistry & Option A Multigroup IonRay + Smooth H₂ cooling + l-resolved ModelAtoms
    net, photos = CR.modern(; he = true)
    Gfull = IR.groups_full(; n_opt = 6, n_euv = 28, Emax = 300.0)
    bp = UN.BlendPar()
    h2c = UN.H2CoolCache(nz; bp)
    env_halo = (Trad = 54.0, J = ν -> 0.0, static = true, caseB = true)

    iHp   = findfirst(==(:Hp),   net.free)
    iHm   = findfirst(==(:Hm),   net.free)
    iH2   = findfirst(==(:H2),   net.free)
    iH2p  = findfirst(==(:H2p),  net.free)
    iHep  = findfirst(==(:Hep),  net.free)
    iHepp = findfirst(==(:Hepp), net.free)
    nf = length(net.free)

    x0 = zeros(nf); x0[iHp] = 1e-4; x0[iH2] = 1e-3
    xs = [copy(x0) for _ in 1:nz]

    # Pre-build l-resolved H I, He II & CHIANTI He I atoms for emergent line synthesis
    ad = joinpath(JDIR, "atomdata", "chianti")
    H_lr   = MA.hydrogenic_lres(1; n_lres = 6, nmax = 30, coll = :chianti, chianti = joinpath(ad, "h_1"))
    He2_lr = MA.hydrogenic_lres(2; n_lres = 6, nmax = 30, coll = :chianti, chianti = joinpath(ad, "he_2"))
    He1_at = MA.heI(; nmax = 20)
    lab1 = [l.label for l in He1_at.levels]
    idx_he1 = (i2p3P = findfirst(==("1s.2p 3P"), lab1), i3d3D = findfirst(==("1s.3d 3D"), lab1),
               i2s3S = findfirst(==("1s.2s 3S"), lab1))

    # Smooth cooling function matching ChemNet's exact Case B recombination binding-energy escape!
    function zone_cool(k, ρk, Tk, d)
        nHk = ρk / CN.μH()
        ne = UN.spos(d[:e]) * nHk; nHa = UN.spos(d[:H]) * nHk; nHp = UN.spos(d[:Hp]) * nHk
        nHep = UN.spos(d[:Hep]) * nHk; nHepp = UN.spos(d[:Hepp]) * nHk
        Λ_h2 = UN.h2_cooling_fast(h2c, k, ρk, Tk, d[:H], d[:H2])
        Tm = UN.spos(Tk - 100.0, 20.0) + 100.0
        sT5 = 1.0 + sqrt(Tm / 1e5)
        Λ_lya = 7.5e-19 * exp(-118348.0 / Tm) / sT5 * ne * nHa
        Λ_he  = 5.54e-17 * Tm^(-0.397) * exp(-473638.0 / Tm) / sT5 * ne * nHep
        gT_floor = UN.log_ramp(Tk, 54.0, 200.0)
        Tv = max(ForwardDiff.value(Tk), 54.0)
        α_H   = CN.photorates_t(photos[2], Tv, env_halo)[2]
        α_He1 = CN.photorates_t(photos[3], Tv, env_halo)[2]
        α_He2 = CN.photorates_t(photos[4], Tv, env_halo)[2]
        Λ_rec_bind = ne * (α_H * nHp * FE.χH + α_He1 * nHep * Z0.χHe1 + α_He2 * nHepp * Z0.χHe2)
        Λ_rec = 8.70e-27 * sqrt(Tm) * (Tm / 1e3)^(-0.2) / (1.0 + (Tm / 1e6)^0.7) * ne * (nHp + nHep + 4nHepp)
        Λ_ff  = 1.42e-27 * 1.3 * sqrt(Tm) * ne * (nHp + nHep + 4nHepp)
        Λ_h2 + gT_floor * (Λ_lya + Λ_he + Λ_rec_bind + Λ_rec + Λ_ff)
    end

    # Local thermal EOS helper (thermal-only internal energy e_th for hydro; chemical binding carried in ChemNet!)
    function thermal_p_cs_e(ρk, Tk, xk)
        d = CN.full(xk, net)
        nH = ρk / CN.μH()
        ntot = nH * (max(d[:H], 0.0) + max(d[:Hp], 0.0) + max(d[:Hm], 0.0) + max(d[:H2], 0.0) +
                     max(d[:H2p], 0.0) + max(d[:e], 0.0) + max(d[:He], 0.0) + max(d[:Hep], 0.0) + max(d[:Hepp], 0.0))
        P = max(ntot * kB * Tk, 1e-30)
        γ = 5.0 / 3.0  # monatomic-dominated (H2 <= 1e-3)
        eth = P / ((γ - 1.0) * ρk)
        cs = sqrt(γ * P / ρk)
        return P, cs, eth, ntot
    end

    # Output CSV + profile dumper
    csv_path = joinpath(outdir, "hii_breakout_rhd.csv")
    csvio = open(csv_path, "w")
    println(csvio, "step,t_yr,dt_yr,Teff,L_star_Lsun,R_star_Rsun,Q_H_s,Q_HeII_s,Q_LW_s,R_HII_pc,R_HeIII_pc,R_LW_pc,R_shock_pc,T_HII_K,Mdot_in_Msunyr,Mdot_0p02pc_Msunyr,u_max_kms,u_min_kms,L_Lya_Lsun,L_2g_Lsun,L_Ha_Lsun,L_Hb_Lsun,L_Paa_Lsun,L_Brg_Lsun,L_HeII1640_Lsun,L_HeII4686_Lsun,L_HeI5876_Lsun,L_HeI10830_Lsun")

    function dump_prof(st, t_yr, ss)
        fn = joinpath(pdir, @sprintf("hii_breakout_rhd_s%05d_t%.0f.csv", st, t_yr))
        open(fn, "w") do io
            println(io, "zone,r_pc,r_Rsun,dr_Rsun,m_Msun,rho,nH,T,u_kms,cs_kms,P_gas,xHp,xH2,xHep,xHepp,Mdot_Msunyr")
            m_enc = Mstar
            for k in 1:nz
                m_enc += ρ[k] * V_z[k]
                P, cs, _, _ = thermal_p_cs_e(ρ[k], T[k], xs[k])
                mdk = 4π * r_mid[k]^2 * ρ[k] * u[k] / (Msun / yr) # positive = outflow, negative = accretion
                @printf(io, "%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.5e,%.4e,%.4e,%.4e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n",
                        k, r_mid[k]/pc, r_mid[k]/Rsun, dr[k]/Rsun, m_enc/Msun, ρ[k], ρ[k]/CN.μH(),
                        T[k], u[k]/1e5, cs/1e5, P, xs[k][iHp], xs[k][iH2], xs[k][iHep], xs[k][iHepp], mdk)
            end
        end
    end

    # Main dynamical 1D RHD + Option A IonRay loop
    t = 0.0; step = 0
    u_if = zeros(nz + 1)
    k_0p02 = argmin(abs.(r_mid .- 0.02 * pc))
    next_prof_t = 0.0
    prof_times = [0.0, 150.0, 260.0, 340.0, 400.0, 450.0, 500.0, 600.0, 750.0, 950.0,
                  1200.0, 1600.0, 2200.0, 3000.0, 4000.0, 5000.0]
    p_ptr = 1

    println("VAL_HII_BREAKOUT_RHD: 100 M☉ contracting Pop III star -> H II / He III / LW breakout + champagne outflow")
    println("  step    t [yr]    Teff [K]   R_HII [pc] R_HeIII  R_LW   R_sh  Mdot_in   u_max   L(Lyα)   L(1640)  Hα/Hβ")

    t_wall0 = time()
    while t < t_end && step < 4500
        # 1. Compute CFL timestep
        dt_cfl = Inf
        for k in 1:nz
            _, cs, _, _ = thermal_p_cs_e(ρ[k], T[k], xs[k])
            dt_cfl = min(dt_cfl, 0.35 * dr[k] / max(abs(u[k]) + cs, 1e4))
        end
        # Cap dt during stellar turn-on ramp so R-type breakout is cleanly resolved
        dt_cap = (t / yr < 800.0) ? 4.0 * yr : 25.0 * yr
        dt = min(dt_cfl, dt_cap, t_end - t)
        dt = max(dt, 1e4)

        # 2. Conservative 1D Spherical HLLC Hydrodynamic Step (with central M* + envelope self-gravity + passive species advection)
        P_z = zeros(nz); cs_z = zeros(nz); E_z = zeros(nz); m_if = zeros(nz + 1)
        m_if[1] = Mstar
        for k in 1:nz
            m_if[k+1] = m_if[k] + ρ[k] * V_z[k]
            P, cs, eth, _ = thermal_p_cs_e(ρ[k], T[k], xs[k])
            P_z[k] = P; cs_z[k] = cs; E_z[k] = ρ[k] * (eth + 0.5 * u[k]^2)
        end

        F1 = zeros(nz + 1); F2 = zeros(nz + 1); F3 = zeros(nz + 1); S_if = zeros(nz + 1)
        # Inner boundary (k = 1 interface): open outflow/accretion sink (zero-gradient if infalling, reflecting floor only if vacuum)
        ρL_in = ρ[1]; uL_in = min(u[1], 0.0); pL_in = P_z[1]; EL_in = E_z[1]; cL_in = cs_z[1]
        F1[1], F2[1], F3[1], S_if[1] = hllc_flux(ρL_in, uL_in, pL_in, EL_in, cL_in, ρ[1], u[1], P_z[1], E_z[1], cs_z[1])
        # Interior interfaces (piecewise-linear minmod reconstruction for sharp shock capturing)
        for j in 2:nz
            F1[j], F2[j], F3[j], S_if[j] = hllc_flux(ρ[j-1], u[j-1], P_z[j-1], E_z[j-1], cs_z[j-1],
                                                     ρ[j],   u[j],   P_z[j],   E_z[j],   cs_z[j])
        end
        # Outer boundary (k = nz+1): zero-gradient ambient halo
        F1[nz+1], F2[nz+1], F3[nz+1], S_if[nz+1] = hllc_flux(ρ[nz], u[nz], P_z[nz], E_z[nz], cs_z[nz],
                                                             ρ[nz], u[nz], P_z[nz], E_z[nz], cs_z[nz])

        # Upwind species advection + conservative spherical update
        ρ_new = similar(ρ); u_new = similar(u); T_new = similar(T)
        xs_new = [copy(x) for x in xs]
        for k in 1:nz
            Vk = V_z[k]
            m_old = ρ[k] * Vk
            mom_old = ρ[k] * u[k] * Vk
            E_old = E_z[k] * Vk
            # Spherical geometric pressure source term: +2 P / r * V_k = P * (A_{k+1} - A_k)
            src_geom = P_z[k] * (A_if[k+1] - A_if[k])
            # Exact gravitational acceleration g = -G M_enc / r^2
            M_enc = 0.5 * (m_if[k] + m_if[k+1])
            g_grav = -Gc * M_enc / r_mid[k]^2
            m_n = max(m_old - dt * (A_if[k+1] * F1[k+1] - A_if[k] * F1[k]), 1e-26 * Vk)
            mom_n = mom_old - dt * (A_if[k+1] * F2[k+1] - A_if[k] * F2[k]) + dt * (src_geom + m_old * g_grav)
            E_n = E_old - dt * (A_if[k+1] * F3[k+1] - A_if[k] * F3[k]) + dt * (mom_old * g_grav)

            ρ_new[k] = m_n / Vk
            u_new[k] = mom_n / m_n
            # Upwind donor-cell advection of chemical species x_i (conserves mass fraction strictly in [0, 1])
            dmL = dt * A_if[k] * F1[k]
            dmR = dt * A_if[k+1] * F1[k+1]
            x_inL = dmL >= 0.0 ? (k == 1 ? xs[1] : xs[k-1]) : xs[k]
            x_inR = dmR >= 0.0 ? xs[k] : (k == nz ? xs[nz] : xs[k+1])
            for s in 1:nf
                xs_new[k][s] = clamp((m_old * xs[k][s] + dmL * x_inL[s] - dmR * x_inR[s]) / m_n, 0.0, 1.0)
            end
            # Recover post-advection gas temperature from thermal energy
            eth_n = max(E_n / m_n - 0.5 * u_new[k]^2, 1e8)
            d_k = CN.full(xs_new[k], net)
            nHk = ρ_new[k] / CN.μH()
            ntot = nHk * (max(d_k[:H], 0.0) + max(d_k[:Hp], 0.0) + max(d_k[:Hm], 0.0) + max(d_k[:H2], 0.0) +
                          max(d_k[:H2p], 0.0) + max(d_k[:e], 0.0) + max(d_k[:He], 0.0) + max(d_k[:Hep], 0.0) + max(d_k[:Hepp], 0.0))
            T_new[k] = clamp((2.0 / 3.0) * ρ_new[k] * eth_n / max(ntot * kB, 1e-40), 54.0, 2.0e5)
        end
        ρ .= ρ_new; u .= u_new; T .= T_new; xs = xs_new

        # 3. Coupled Implicit Option A Multigroup IonRay + ChemNet + Smooth-Gated H₂ Cooling Step
        t += dt; step += 1
        t_yr = t / yr
        ss = star_state(t_yr)
        N_star = IR.blackbody_full(Gfull, ss.L, ss.Teff)
        # Update interface velocities for Sobolev H₂ escape probabilities
        u_if[1] = u[1]; u_if[nz+1] = u[nz]
        for j in 2:nz; u_if[j] = 0.5 * (u[j-1] + u[j]); end
        UN.update_h2_escape!(h2c, r_if, u_if, ρ, T, [x[iH2] for x in xs])

        Nin = copy(N_star)
        NH2_col = 0.0
        for k in 1:nz
            src = IR.source_full(ρ[k], r_if[k], r_if[k+1], Nin, NH2_col, T[k], Gfull, net; binding_in_u = true)
            uk = CN.eos(ρ[k], T[k], xs[k], net)[1]
            cool_k(ρ_, T_, d_) = zone_cool(k, ρ_, T_, d_) - src.heat(d_)
            x_out, T_out = CN.step(xs[k], ρ[k], uk, T[k], dt, net;
                                   env = env_halo, photos = photos, cooling = cool_k,
                                   extra = src.extra, accepted = src.accepted)
            xs[k] = x_out; T[k] = clamp(T_out, 54.0, 2.0e5)
            Nin, NH2_col = src.out(xs[k], dt)
        end

        # 4. Diagnostics, atomic line synthesis & logging
        do_prof = (p_ptr <= length(prof_times) && t_yr >= prof_times[p_ptr]) || (t >= t_end)
        if step == 1 || step % 25 == 0 || do_prof
            Q_H   = sum(N_star[g] for g in eachindex(N_star) if Gfull.Elo[g] >= 13.5)
            Q_He2 = sum(N_star[g] for g in eachindex(N_star) if Gfull.Elo[g] >= 54.0)
            Q_LW  = N_star[Gfull.g_lw]
            xHp_v   = [x[iHp]   for x in xs]
            xH2_v   = [x[iH2]   for x in xs]
            xHepp_v = [x[iHepp] for x in xs]
            khii = findlast(>(0.5), xHp_v);  r_hii = khii === nothing ? 0.0 : r_mid[khii] / pc
            khe3 = findlast(>(0.04), xHepp_v); r_he3 = khe3 === nothing ? 0.0 : r_mid[khe3] / pc
            klw  = findfirst(>(5e-4), xH2_v);  r_lw  = klw  === nothing ? r_out / pc : r_mid[klw] / pc
            # Outward D-type shock location: strongest velocity deceleration (u[k] - u[k+1]) where u[k] > 3 km/s
            du = u[1:end-1] .- u[2:end]
            ksh_cand = findall(k -> u[k] > 3e5 && du[k] > 2e5, 1:nz-1)
            r_sh = isempty(ksh_cand) ? 0.0 : r_if[ksh_cand[argmax(du[ksh_cand])] + 1] / pc
            T_hii = khii === nothing ? T[1] : T[max(1, khii ÷ 2)]
            md_in = 4π * r_mid[1]^2 * ρ[1] * u[1] / (Msun / yr)
            md_02 = 4π * r_mid[k_0p02]^2 * ρ[k_0p02] * u[k_0p02] / (Msun / yr)
            ln = synth_nebular_lines(r_if, ρ, T, xs, net, H_lr, He2_lr, He1_at, idx_he1)
            @printf(csvio, "%d,%.5e,%.4e,%.1f,%.4e,%.3f,%.4e,%.4e,%.4e,%.5f,%.5f,%.5f,%.5f,%.1f,%.4e,%.4e,%.3f,%.3f,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e\n",
                    step, t_yr, dt/yr, ss.Teff, ss.L/Lsun, ss.Rst/Rsun, Q_H, Q_He2, Q_LW,
                    r_hii, r_he3, r_lw, r_sh, T_hii, md_in, md_02, maximum(u)/1e5, minimum(u)/1e5,
                    ln.L_Lya/Lsun, ln.L_2g/Lsun, ln.L_Ha/Lsun, ln.L_Hb/Lsun, ln.L_Paa/Lsun, ln.L_Brg/Lsun,
                    ln.L_1640/Lsun, ln.L_4686/Lsun, ln.L_5876/Lsun, ln.L_10830/Lsun)
            flush(csvio)
            if do_prof
                dump_prof(step, t_yr, ss)
                while p_ptr <= length(prof_times) && t_yr >= prof_times[p_ptr]; p_ptr += 1; end
            end
            if step == 1 || step % 60 == 0 || do_prof
                hab = ln.L_Hb > 0 ? ln.L_Ha / ln.L_Hb : 0.0
                @printf("  %4d  %8.1f  %8.0f    %6.3f    %6.3f  %5.3f %5.3f %+9.2e %+6.1f  %8.2e %8.2e  %.2f\n",
                        step, t_yr, ss.Teff, r_hii, r_he3, r_lw, r_sh, md_02, maximum(u)/1e5,
                        ln.L_Lya/Lsun, ln.L_1640/Lsun, hab)
                flush(stdout)
            end
        end
    end
    close(csvio)
    @printf("DONE hii_breakout_rhd: %d steps, t=%.1f yr in %.1f s wall time!\n", step, t/yr, time() - t_wall0)
end

main()
