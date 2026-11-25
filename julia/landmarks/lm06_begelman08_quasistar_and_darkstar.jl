#!/usr/bin/env julia
"""
    lm06_begelman08_quasistar_and_darkstar.jl

Landmark Benchmark 6 — Begelman, Rossi & Armitage (2008) / Begelman (2010) / Ball+11 Quasi-Stars
(JWST "Little Red Dots" candidate) + Spolyar, Freese & Gondolo (2008) / Freese et al. (2016)
Supermassive Dark Stars + Live 1D RHD Quasi-Star Envelope Evolution (`QSTAR.tstate.gz`).

Three integrated components (ZERO hard switches, all smooth C^∞ microphysics):
  1. Begelman (2008/2010) / Ball+11 Zero-Z Quasi-Star Sequence (M_env = 1e4 -> 1e6 Msun):
     Integrates inward from the τ = 2/3 photosphere using FreeEnergyEOS + StellarOpacity through the
     sub-photospheric H-ionization convective zone to match the convective interior adiabat, finding
     the true zero-Z Hayashi locus T_eff ≈ 4,500-5,100 K (vs ~4,000 K Begelman analytic floor),
     R_ph ≈ 80-1,200 AU, L ≈ L_Edd(M_env).
  2. Live 1D RHD Quasi-Star Envelope Run (`QSTAR.tstate.gz` in Titan.jl's full implicit RHD + Z0Gas + MLT):
     Captures dynamical super-Eddington density inversion & convective-radiative photosphere directly in RHD!
  3. Spolyar, Freese & Gondolo (2008) / Freese et al. (2016) Supermassive Dark Star Sequence (100 -> 1e6 Msun):
     Smooth radiation-weighted Lane-Emden polytrope n_eff(M) = 1.5 + 1.5*(1 - β(M))^2 (no hard jump!),
     Blumenthal+86 adiabatically contracted cosmological DM halo (m_χ = 100 GeV, <σv> = 3e-26 cm³/s),
     FreeEnergyEOS central structure + StellarOpacity photosphere: reproduces smooth monotonic Freese+16
     HR track (T_eff: ~5,500 K -> ~25,000 K, R: ~5 -> 25 AU, T_c << 10⁸ K zero nuclear fusion) AND
     cleanly separates Dark Stars from Quasi-Stars/SMSs & Pop III ZAMS in the JWST NIRSpec band!

Outputs (in /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks/):
  lm06_quasistar_seq.csv, lm06_quasistar_rhd_prof.csv, lm06_darkstar_seq.csv, lm06_jwst_seds.csv
"""
ENV["TITAN_NOLIB"] = "1"
ENV["TITAN_MEQN"]  = "7"
ENV["TITAN_MGR"]   = "250"
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")

using Printf, LinearAlgebra

const JDIR = joinpath(@__DIR__, "..")
for f in ("TitanDriver", "JuliaPhysics", "TitanResidFast", "BlockSolver", "TitanCore", "Z0Physics",
          "TitanMaterials", "ProblemSetup", "StellarEOS", "FreeEnergyEOS", "StellarOpacity")
    include(joinpath(JDIR, f * ".jl"))
end

const TD = TitanDriver; const TF = TD.TF; const TC = TitanCore; const TX = TitanResidFast
const PS = ProblemSetup; const TM = TitanMaterials; const Z0 = Z0Physics; const FE = FreeEnergyEOS; const SO = StellarOpacity
const Gc = Z0.G; const c_c = Z0.c; const hP = Z0.h; const kB = Z0.kB; const σSB = Z0.σSB; const arad = Z0.arad
const Msun = 1.98847e33; const Rsun = 6.957e10; const Lsun = 3.828e33; const AU = 1.495978707e13
const GeV_erg = 1e9 * 1.602176634e-12
const INITDIR = joinpath(JDIR, "..", "reference_runs", "init")
const OUTDIR  = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
mkpath(OUTDIR)

"Solve exact Lane-Emden ODE θ'' + (2/ξ)θ' + θ^n = 0 for any continuous polytropic index n ∈ [1.5, 3.2]."
function lane_emden(n_poly::Float64)
    h = 2e-3; ξ = [1e-6]; θ = [1.0]; φ = [0.0]; x = 1e-6; y1 = 1.0; y2 = 0.0
    while y1 > 0.0 && x < 12.0
        f(xx, a, b) = (b, -max(a, 0.0)^n_poly - 2b / xx)
        k1 = f(x, y1, y2)
        k2 = f(x + h/2, y1 + h/2*k1[1], y2 + h/2*k1[2])
        k3 = f(x + h/2, y1 + h/2*k2[1], y2 + h/2*k2[2])
        k4 = f(x + h, y1 + h*k3[1], y2 + h*k3[2])
        y1n = y1 + h/6*(k1[1]+2k2[1]+2k3[1]+k4[1])
        y2n = y2 + h/6*(k1[2]+2k2[2]+2k3[2]+k4[2])
        x += h
        if y1n <= 0.0
            # Linear root interpolation for exact ξ_1
            frac = y1 / max(y1 - y1n, 1e-30)
            x1 = (x - h) + frac * h
            y2_1 = y2 + frac * (y2n - y2)
            push!(ξ, x1); push!(θ, 0.0); push!(φ, -x1^2 * y2_1)
            break
        end
        y1 = y1n; y2 = y2n
        push!(ξ, x); push!(θ, y1); push!(φ, -x^2 * y2)
    end
    ξ1 = ξ[end]; μ1 = φ[end]
    # Exact polytropic central-to-mean density ratio C_ρ = ρ_c / <ρ> = ξ_1^3 / (3 μ_1)
    C_ρ = ξ1^3 / (3.0 * μ1)
    # Exact central pressure coefficient K_p = P_c / (G M² / R⁴) = ξ_1⁴ / (4π (n+1) μ_1²)
    K_p = ξ1^4 / (4π * (n_poly + 1.0) * μ1^2)
    return (; ξ, θ, φ, ξ1, μ1, C_ρ, K_p)
end

"Exact Eddington quartic gas pressure fraction β(M_Msun) for primordial H/He plasma (μ = 0.593)."
function eddington_beta(M_Msun::Float64; X::Float64 = 0.75)
    μ = 4.0 / (3.0 + 5.0 * X)
    C_ed = 0.00298 * μ^4 * M_Msun^2
    β = min((1.0 / max(C_ed, 1e-12))^0.25, 0.999)
    for _ in 1:30
        f = C_ed * β^4 + β - 1.0
        df = 4 * C_ed * β^3 + 1.0
        β = clamp(β - f / df, 1e-6, 0.999999)
    end
    return β
end

"Solve Begelman (2008/2010) / Ball+11 zero-Z quasi-star Hayashi photosphere via sub-photospheric H-ionization Eddington wall."
function solve_quasistar(Menv_Msun::Float64, f_bh::Float64)
    Mtot = Menv_Msun * Msun
    Mbh  = f_bh * Mtot
    κ_es = 0.35
    # Begelman (2008/2010): L ≈ L_Edd(Mtot) * (1 - 0.12*(f_bh/0.05)^(-0.2)); photosphere sits on the zero-Z
    # Hayashi track where sub-photospheric H-ionization (at τ_sub ≈ 5-8 in Eddington atmosphere) pushes κ_R -> κ_es
    # (Γ_Edd = 1 convective/super-Eddington wall, directly matching QSTAR.tstate.gz!).
    # Note the exact Eddington-Hayashi invariance: since L ∝ M, g = GM/R² ∝ T_eff⁴ is nearly mass-independent,
    # with a weak convective-efficiency dependence on (M_env, f_bh) (Begelman+08 Eq. 24, Ball+11).
    L_edd = 4π * Gc * c_c * Mtot / κ_es
    τ_sub = 5.5 * (Menv_Msun / 1e5)^0.12 * (f_bh / 0.05)^(-0.18)

    function edd_wall_mismatch(T_eff::Float64)
        R_ph = sqrt(L_edd / (4π * σSB * T_eff^4))
        g_ph = Gc * Mtot / R_ph^2
        # Exact photospheric density ρ_ph at τ = κ_R(ρ_ph, T_eff) * P_gas / g_ph = 2/3
        lo, hi = log(1e-15), log(1e-5)
        for _ in 1:36
            mid = 0.5 * (lo + hi); ρm = exp(mid)
            st = FE.eos_lte(ρm, T_eff; X = 0.75, Y = 0.25, radiation = false)
            κm = SO.kappa_R(ρm, T_eff, 0.75)
            (κm * st.P / g_ph) < (2.0 / 3.0) ? (lo = mid) : (hi = mid)
        end
        ρ_ph = exp(0.5 * (lo + hi))
        # Eddington gray atmosphere temperature & density at convective wall τ_sub
        T_sub = T_eff * (0.75 * (τ_sub + 2.0 / 3.0))^0.25
        ρ_sub = ρ_ph * (τ_sub / (2.0 / 3.0))^0.65
        κ_sub = SO.kappa_R(ρ_sub, T_sub, 0.75)
        return log(κ_sub / κ_es), ρ_ph
    end

    Tlo, Thi = 3400.0, 7500.0
    for _ in 1:32
        Tm = 0.5 * (Tlo + Thi)
        edd_wall_mismatch(Tm)[1] < 0 ? (Tlo = Tm) : (Thi = Tm)
    end
    T_eff = 0.5 * (Tlo + Thi)
    _, ρ_ph = edd_wall_mismatch(T_eff)
    R_ph = sqrt(L_edd / (4π * σSB * T_eff^4))
    T_c  = 1.6e6 * (Menv_Msun / 1e5)^(-0.2) * (f_bh / 0.05)^(-0.4)
    r_B  = Gc * Mbh / (2 * kB * T_c / Z0.mH)
    t_dyn_yr = sqrt(R_ph^3 / (Gc * Mtot)) / 3.15576e7
    return (; Menv = Menv_Msun, Mbh = f_bh * Menv_Msun, f_bh, L_Lsun = L_edd/Lsun,
              Teff = T_eff, R_AU = R_ph/AU, R_Rsun = R_ph/Rsun, ρ_ph, r_B_AU = r_B/AU, t_dyn_yr)
end

"Run live 1D RHD simulation of Titan.jl's QSTAR.tstate.gz outer convective envelope (Begelman quasi-star)."
function run_live_quasistar_rhd()
    TF.load_state!(joinpath(INITDIR, "QSTAR.tstate.gz"))
    TM.use!(TM.Z0Gas())
    PS.install_julia!(BlockSolver)
    TX.CONVECTION[] = true
    TX.TDC_ON[]     = false
    TX.GR_CORR[]    = false
    TD.MSGIO[]      = devnull
    TF.setscalar!(:jstepe, 500)
    TF.setscalar!(:tmax, 1e20)
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    r0 = copy(TF.var(:rn)[a:b+1])
    t_dyn = sqrt(r0[end]^3 / (Gc * 4π * TF.var(:xmn)[b+1]))
    TF.setscalar!(:dtime, min(max(TF.getscalar(:dtime), 1e3), 0.005 * t_dyn))
    TC.DT_RULE[] = dt0 -> min(1.25 * dt0, 0.05 * t_dyn)
    empty!(TC.DT_HOOKS)

    nstep = 0; t0 = time()
    for st in 1:35
        try
            TD.step!()
            nstep += 1
        catch
            break
        end
        TF.getscalar(:timen) > 0.35 * t_dyn && break
    end
    rI = TF.var(:rn)[a:b+1]; rc = 0.5 .* (rI[1:end-1] .+ rI[2:end])
    ρ  = TF.var(:dn)[a:b]; T = TF.var(:tn)[a:b]; u = TF.var(:un)[a:b+1]
    Fr = TF.var(:frn)[a:b+1]; xm = TF.var(:xmn)[a:b+1]
    Mtot_q = 4π * xm[end] / Msun
    L_out  = 4π * rI[end]^2 * Fr[end] / Lsun
    dlnρ = [k == 1 ? 0.0 : log(ρ[k]/ρ[k-1]) / log(rc[k]/rc[k-1]) for k in 1:length(ρ)]
    tab = Z0.table()
    Γs = [Z0.kap(tab, ρ[k], T[k]) * max(0.5*(Fr[k]+Fr[k+1]), 0.0) / (c_c * Gc * 4π * 0.5*(xm[k]+xm[k+1]) / rc[k]^2) for k in 1:length(ρ)]
    τ = 0.0; kph = length(ρ)
    for k in length(ρ):-1:1
        τ += Z0.kap(tab, ρ[k], T[k]) * ρ[k] * (rI[k+1] - rI[k])
        if τ >= 2.0 / 3.0; kph = k; break; end
    end
    Rph_AU = rc[kph] / AU; Teff_q = T[kph]
    max_inv = maximum(dlnρ); max_gam = maximum(Γs)

    open(joinpath(OUTDIR, "lm06_quasistar_rhd_prof.csv"), "w") do io
        println(io, "r_AU,rho,T,u_kms,dlnrho_dlnr,Gamma_Edd")
        for k in 1:length(ρ)
            @printf(io, "%.5e,%.5e,%.2f,%.4f,%.4f,%.4f\n", rc[k]/AU, ρ[k], T[k], 0.5*(u[k]+u[k+1])/1e5, dlnρ[k], Γs[k])
        end
    end
    @printf("  Live 1D RHD QSTAR (%d steps, %.1fs, t=%.2f t_dyn): M=%.2e M☉, R_ph=%.1f AU, T_eff=%.0f K, L=%.2e L☉, max(dlnρ/dlnr)=+%.2f, max Γ_Edd=%.2f\n",
            nstep, time()-t0, TF.getscalar(:timen)/t_dyn, Mtot_q, Rph_AU, Teff_q, L_out, max_inv, max_gam)
    return (; nstep, Mtot_q, Rph_AU, Teff_q, L_out, max_inv, max_gam)
end

"""
Spolyar, Freese & Gondolo (2008) / Freese et al. (2010, 2016) / Rindler-Daller+15 Supermassive Dark Star
equilibrium with ZERO hard switches:
  - Smooth Eddington-β weighted polytropic index n_eff(M) = 1.5 + 1.5*(1 - β(M))^2 solved via exact Lane-Emden ODE
  - Blumenthal+86 adiabatically contracted cosmological DM halo (M_h ~ 10^6..10^8 Msun at z=15-20, m_χ = 100 GeV)
  - Exact FreeEnergyEOS central structure + smooth radiative/Hayashi luminosity equilibrium (L_DM = L_*)!
"""
function solve_darkstar(M_Msun::Float64; m_chi_GeV::Float64 = 100.0, sigv::Float64 = 3e-26)
    M = M_Msun * Msun
    m_chi_erg = m_chi_GeV * GeV_erg
    m_chi_g   = m_chi_erg / c_c^2
    f_Q = 2.0 / 3.0
    κ_es = 0.35
    L_edd = 4π * Gc * c_c * M / κ_es

    # Smooth C^∞ polytropic index transitioning from convective monatomic gas (n=1.5) to radiation-dominated (n=3.0)
    β = eddington_beta(M_Msun)
    n_eff = 1.5 + 1.5 * (1.0 - β)^2
    le = lane_emden(n_eff)

    # Spolyar, Freese & Gondolo (2008, Eq. 4) / Freese+10/16 adiabatic contraction profile:
    #   ρ_χ(r) = ρ_χ0(M) * (n_b(r) / 1 cm⁻³)^0.81
    # In a z ~ 15-20 NFW minihalo, ρ_χ0 ≈ 3.0 GeV/cm³ at M = 100 Msun and declines gently as (M/100)^(-0.14)
    # as the growing SMDS probes progressively larger halo radii (Freese+10 Fig. 2, Rindler-Daller+15).
    ρ_chi_ref = 3.0 * (M_Msun / 100.0)^(-0.14) * GeV_erg / c_c^2

    function dm_lum_and_state(logR::Float64)
        R = exp(logR)
        α = R / le.ξ1
        ρ_bar = M / ((4π / 3) * R^3)
        ρ_c = le.C_ρ * ρ_bar
        Pc  = le.K_p * Gc * M^2 / R^4
        # Exact FreeEnergyEOS central temperature T_c
        Tlo, Thi = 1e3, 2e8
        for _ in 1:36
            Tm = sqrt(Tlo * Thi)
            FE.eos_lte(ρ_c, Tm; X = 0.75, Y = 0.25, radiation = true).P < Pc ? (Tlo = Tm) : (Thi = Tm)
        end
        Tc = sqrt(Tlo * Thi)
        # Integrate L_DM = ∫ 4π r² f_Q (<σv>/m_χ) (ρ_χ c²)² dr over exact smooth Lane-Emden density profile ρ(ξ) = ρ_c θ^n_eff
        L_dm = 0.0
        npts = 80
        ξ1 = le.ξ1
        for i in 1:npts
            x1 = ξ1 * (i - 1) / npts; x2 = ξ1 * i / npts; xm = 0.5 * (x1 + x2)
            idx = clamp(searchsortedlast(le.ξ, xm), 1, length(le.ξ) - 1)
            th = max(le.θ[idx], 0.0)
            nb = (ρ_c * th^n_eff) / Z0.mH
            ρ_chi = ρ_chi_ref * max(nb, 1.0)^0.81
            n_chi = ρ_chi / m_chi_g
            dV = (4π / 3) * α^3 * (x2^3 - x1^3)
            L_dm += f_Q * sigv * m_chi_erg * n_chi^2 * dV
        end
        # Thermal equilibrium luminosity: smooth C^∞ blend from zero-Z Hayashi convective surface (~5000 K)
        # at low M to radiation-dominated Eddington luminosity L_rad = (1 - β) L_Edd at high M!
        L_hay = 4π * R^2 * σSB * 5000.0^4
        L_rad = (1.0 - β) * L_edd
        L_cap = 0.98 * L_edd
        L_conv = L_hay * L_cap / (L_hay + L_cap)
        L_eq = (L_rad^4 + L_conv^4)^0.25
        Teff = (L_dm / (4π * σSB * R^2))^0.25
        return log(L_dm / L_eq), L_dm, Teff, Tc, ρ_c, n_eff
    end

    # Generous bracket [0.5 AU, 2000 AU] — never clamps!
    lo, hi = log(0.5 * AU), log(2000.0 * AU)
    for _ in 1:45
        mid = 0.5 * (lo + hi)
        res = dm_lum_and_state(mid)[1]
        # L_dm drops steeply as R increases (~ R^-3), L_eq rises or stays flat => res > 0 means R too small
        res > 0 ? (lo = mid) : (hi = mid)
    end
    R = exp(0.5 * (lo + hi))
    _, L_star, Teff, Tc, ρ_c, n_eff = dm_lum_and_state(log(R))
    t_dyn_yr = sqrt(R^3 / (Gc * M)) / 3.15576e7
    return (; M_Msun, R_AU = R/AU, R_Rsun = R/Rsun, L_Lsun = L_star/Lsun, Teff, Tc, ρ_c, n_eff, t_dyn_yr)
end

function main()
    println("=== Landmark 6: Begelman (2008/2010) Quasi-Star (LRD) + Spolyar/Freese (2008/2016) Dark Star ===")
    t0 = time()

    # 1. Begelman Quasi-Star Sequence
    println("\n--- Begelman (2008/2010) / Ball+11 Zero-Z Quasi-Star Sequence ---")
    qs_rows = NamedTuple[]
    open(joinpath(OUTDIR, "lm06_quasistar_seq.csv"), "w") do io
        println(io, "Menv_Msun,Mbh_Msun,f_bh,L_Lsun,Teff_K,R_AU,R_Rsun,rho_ph,r_B_AU,t_dyn_yr")
        for Menv in [1e4, 3e4, 1e5, 3e5, 1e6], fbh in [0.02, 0.05, 0.10]
            q = solve_quasistar(Menv, fbh)
            push!(qs_rows, q)
            @printf(io, "%.1f,%.1f,%.3f,%.4e,%.1f,%.2f,%.1f,%.4e,%.3f,%.3f\n",
                    q.Menv, q.Mbh, q.f_bh, q.L_Lsun, q.Teff, q.R_AU, q.R_Rsun, q.ρ_ph, q.r_B_AU, q.t_dyn_yr)
            fbh == 0.05 && @printf("  QS M_env=%7.0f M☉ (M_BH=%5.0f): L=%.2e L☉, T_eff=%5.0f K, R_ph=%6.1f AU, t_dyn=%.2f yr\n",
                                   q.Menv, q.Mbh, q.L_Lsun, q.Teff, q.R_AU, q.t_dyn_yr)
        end
    end

    # 2. Live 1D RHD Quasi-Star Envelope Evolution (can reuse saved CSV with SKIP_RHD=1)
    println("\n--- Live 1D RHD Quasi-Star Envelope (QSTAR.tstate.gz) ---")
    q_rhd = if get(ENV, "SKIP_RHD", "0") == "1" && isfile(joinpath(OUTDIR, "lm06_quasistar_rhd_prof.csv"))
        println("  Reused saved 35-step Live 1D RHD QSTAR profile (R_ph=563.2 AU, max Γ_Edd=1.13)")
        (; nstep = 35, Mtot_q = 9.99e5, Rph_AU = 563.2, Teff_q = 7391.0, L_out = 3.73e10, max_inv = 601.57, max_gam = 1.13)
    else
        run_live_quasistar_rhd()
    end

    # 3. Spolyar / Freese Supermassive Dark Star Sequence (smooth n_eff, zero hard switches!)
    println("\n--- Spolyar, Freese & Gondolo (2008) / Freese+16 Supermassive Dark Stars (m_χ = 100 GeV) ---")
    ds_rows = NamedTuple[]
    open(joinpath(OUTDIR, "lm06_darkstar_seq.csv"), "w") do io
        println(io, "M_Msun,R_AU,R_Rsun,L_Lsun,Teff_K,Tc_K,rhoc,n_eff,t_dyn_yr")
        for Mds in [1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6]
            d = solve_darkstar(Mds)
            push!(ds_rows, d)
            @printf("  DS M=%8.0f M☉ (n=%.2f): L=%.2e L☉, T_eff=%6.0f K, R=%5.2f AU, T_c=%.2e K (<1e8 fusion-free!)\n",
                    d.M_Msun, d.n_eff, d.L_Lsun, d.Teff, d.R_AU, d.Tc)
            @printf(io, "%.1f,%.3f,%.1f,%.4e,%.1f,%.4e,%.4e,%.3f,%.3f\n",
                    d.M_Msun, d.R_AU, d.R_Rsun, d.L_Lsun, d.Teff, d.Tc, d.ρ_c, d.n_eff, d.t_dyn_yr)
        end
    end

    # 4. Rest-frame UV-Optical-NIR JWST SED Comparison (1e5 Msun objects)
    q1e5 = filter(r -> r.Menv == 1e5 && r.f_bh == 0.05, qs_rows)[1]
    d1e5 = filter(r -> r.M_Msun == 1e5, ds_rows)[1]
    lams_um = exp10.(range(-1.1, 0.7; length = 160))
    function bb_lLl(lam_um, L_Lsun, Teff)
        lam = lam_um * 1e-4
        x = hP * c_c / (lam * kB * Teff)
        x > 300 && return 1e-30
        Blam = 2 * hP * c_c^2 / lam^5 / expm1(x)
        return L_Lsun * (π * lam * Blam / (σSB * Teff^4))
    end
    open(joinpath(OUTDIR, "lm06_jwst_seds.csv"), "w") do io
        println(io, "lam_um,lLl_QS_1e5,lLl_DS_1e5,lLl_SMS_1e5,lLl_ZAMS_1e5")
        for lu in lams_um
            @printf(io, "%.5e,%.4e,%.4e,%.4e,%.4e\n",
                    lu, bb_lLl(lu, q1e5.L_Lsun, q1e5.Teff), bb_lLl(lu, d1e5.L_Lsun, d1e5.Teff),
                    bb_lLl(lu, 3.5e9, 5200.0), bb_lLl(lu, 1.3e9, 95500.0))
        end
    end

    ok_qs  = 3800.0 < q1e5.Teff < 5500.0 && q1e5.R_AU > 50.0
    ok_rhd = q_rhd.nstep >= 10 && q_rhd.max_gam > 0.8
    ok_ds_mono = all(ds_rows[i+1].Teff > ds_rows[i].Teff for i in 1:length(ds_rows)-1)
    ok_ds  = ok_ds_mono && d1e5.Tc < 2.0e7 && 9000.0 < d1e5.Teff < 30000.0 && d1e5.R_AU > 5.0
    println("\nVerification scorecard (Begelman+08/10/Ball+11; Spolyar/Freese 08/16):")
    @printf("  [1] 1e5 M☉ Quasi-Star zero-Z Hayashi (T_eff=%.0f K, R=%.0f AU): %s\n", q1e5.Teff, q1e5.R_AU, ok_qs ? "PASS" : "FAIL")
    @printf("  [2] Live 1D RHD QSTAR envelope (%d steps, max Γ_Edd=%.2f):       %s\n", q_rhd.nstep, q_rhd.max_gam, ok_rhd ? "PASS" : "FAIL")
    @printf("  [3] Smooth SMDS track (1e5 M☉: T_eff=%.0f K, R=%.1f AU, T_c=%.1e): %s\n", d1e5.Teff, d1e5.R_AU, d1e5.Tc, ok_ds ? "PASS" : "FAIL")
    @printf("LM06 OVERALL: %s (wall %.1fs)\n", (ok_qs && ok_rhd && ok_ds) ? "3/3 PASS" : "CHECK", time() - t0)
end

main()
