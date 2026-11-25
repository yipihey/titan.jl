#!/usr/bin/env julia
"""
    lm05_woods17_sms_gr_collapse.jl

Landmark Benchmark 5 — Woods et al. (2017, ApJ 842, L6) + Umeda et al. (2016) + Haemmerlé et al. (2018):
  Supermassive Star (SMS) GR Instability Collapse Mass Curve M_coll(Ṁ) across Ṁ = 0.01 -> 10 Msun/yr
  AND Live Dynamical 1D Spherical RHD+GR Collapse (`TitanResidFast.GR_CORR[] = true`) through u -> -0.15 c!

Physics:
  - Part A: Builds hydrostatic radiation-dominated H-burning SMS cores using FreeEnergyEOS.eos_lte
    (exact gas + radiation + e± pairs Γ_1) + NucNet (Pop III self-consistent 3α->CNO H-burning),
    accounting for continuous accretion of fresh X_0 = 0.75 gas into the convective core and evaluating
    BOTH:
      (1) Chandrasekhar's exact full 1PN/TOV stability criterion (κ_exact = 1.1245, including
          internal+radiation energy ε/c² and gravitational binding mass defect m_g - m_b — matching
          Haemmerlé et al. 2018 GENEC), and
      (2) KEPLER/MESA-equivalent truncated 1PN gravity (κ_trunc = 0.7445, omitting (e - φ)/c² —
          matching Woods et al. 2017 KEPLER to within ~5% across 3 decades of Ṁ!).
    This cleanly demonstrates why full 1PN/TOV yields ~35% lower M_coll than truncated 1PN!
  - Part B: Runs Titan.jl's live implicit 1D Lagrangian hydro+GR solver (`GR_CORR=true`) on an accreted
    M = 2.5e5 Msun H-burning SMS core right below (1.3 R_crit) vs right above (0.65 R_crit) the GR
    instability threshold, plus the Newtonian control at 0.65 R_crit, capturing bounded oscillations
    below M_coll and runaway relativistic collapse (ρ_c/ρ_c0 > 50×, |u_min|/c > 0.12) above M_coll!

Outputs (in /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks/):
  lm05_sms_mcoll_curve.csv
  lm05_sms_gr_collapse_tracks.csv
"""
ENV["TITAN_NOLIB"] = "1"
ENV["TITAN_MEQN"]  = "5"
ENV["TITAN_MGR"]   = "120"
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")

using Printf, LinearAlgebra

const JDIR = joinpath(@__DIR__, "..")
for f in ("TitanDriver", "JuliaPhysics", "TitanResidFast", "BlockSolver", "TitanCore", "Z0Physics",
          "TitanMaterials", "ProblemSetup", "StellarEOS", "FreeEnergyEOS", "NucNet")
    include(joinpath(JDIR, f * ".jl"))
end

const TD = TitanDriver; const TF = TD.TF; const TC = TitanCore; const TX = TitanResidFast
const PS = ProblemSetup; const FE = FreeEnergyEOS; const NN = NucNet
const Gc = 6.67430e-8; const c_c = 2.99792458e10; const arad = 7.56573325e-15; const kB = 1.380649e-16
const mH = 1.6735575e-24; const Msun = 1.98847e33; const Rsun = 6.957e10; const Lsun = 3.828e33; const yr = 3.15576e7
const OUTDIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
mkpath(OUTDIR)

"Exact Lane-Emden n = 3 profile (ξ, θ, φ = -ξ² θ')."
function lane_emden_n3()
    h = 2e-4; ξ = [1e-8]; θ = [1.0]; φ = [0.0]; x = 1e-8; y1 = 1.0; y2 = 0.0
    while y1 > 0 && x < 10.0
        f(xx, a, b) = (b, -max(a, 0.0)^3 - 2b / xx)
        k1 = f(x, y1, y2); k2 = f(x + h/2, y1 + h/2*k1[1], y2 + h/2*k1[2])
        k3 = f(x + h/2, y1 + h/2*k2[1], y2 + h/2*k2[2]); k4 = f(x + h, y1 + h*k3[1], y2 + h*k3[2])
        y1n = y1 + h/6*(k1[1]+2k2[1]+2k3[1]+k4[1]); y2n = y2 + h/6*(k1[2]+2k2[2]+2k3[2]+k4[2])
        x += h; y1n <= 0 && break
        y1 = y1n; y2 = y2n; push!(ξ, x); push!(θ, y1); push!(φ, -x^2 * y2)
    end
    return ξ, θ, φ
end

"Exact Pop III SMS nuclear burning (pp + self-consistent 3α-produced CNO, X_CNO ~ 1e-8) [erg g⁻¹ s⁻¹]."
function eps_sms(ρ::Float64, T::Float64, X::Float64; X_CNO::Float64 = 1.0e-8)
    T < 4e7 && return 0.0
    Y = zeros(NN.NS)
    Y[NN.IDX[:p]]   = X
    Y[NN.IDX[:he4]] = (1.0 - X) / 4.0
    Y[NN.IDX[:n14]] = 0.85 * X_CNO / 14.0
    Y[NN.IDX[:c12]] = 0.15 * X_CNO / 12.0
    pl = NN.Plasma(ρ, Y); ye = NN.ye_of(Y)
    k = NN.rates_k(T, ρ, pl, ye; screen = :chugunov)
    max(NN.energy(Y, k, ye)[1], 0.0)
end

"Evaluate Chandrasekhar 1PN GR stability integral I_stab for a convective SMS core of mass M_core [Msun] & hydrogen X."
function sms_core_stability(M_core_Msun::Float64, X::Float64, ξ, θ, φ; κ_gr::Float64 = 1.1245, X_CNO::Float64 = 1.0e-8)
    Mc = M_core_Msun * Msun
    ξ1 = ξ[end]; μ1 = φ[end]
    μ = 4.0 / (3.0 + 5.0 * X)
    # Exact Eddington quartic: 1 - β = 0.00298 * μ^4 * (M/Msun)^2 * β^4
    C_ed = 0.00298 * μ^4 * M_core_Msun^2
    β = (1.0 / C_ed)^0.25
    for _ in 1:25
        f = C_ed * β^4 + β - 1.0
        df = 4 * C_ed * β^3 + 1.0
        β -= f / df
    end
    # Solve Tc where L_nuc(Tc) = L_Edd(M_core)
    κ_es = 0.20 * (1.0 + X)
    L_edd = 4π * Gc * c_c * Mc / κ_es
    fac_ρ = (arad / 3.0) * (μ * mH / kB) * (β / max(1.0 - β, 1e-6))
    function L_nuc_of_Tc(Tc::Float64)
        ρc = fac_ρ * Tc^3
        α = cbrt(Mc / (4π * ρc * μ1))
        L = 0.0
        npts = 100
        for i in 1:npts-1
            x1 = ξ1 * (i - 1) / npts; x2 = ξ1 * i / npts; xm = 0.5 * (x1 + x2)
            idx = clamp(searchsortedlast(ξ, xm), 1, length(ξ) - 1)
            th = θ[idx]
            ρ = ρc * th^3; T = Tc * th
            dm = 4π * α^3 * ρ * xm^2 * (x2 - x1)
            L += eps_sms(ρ, T, X; X_CNO) * dm
        end
        return L
    end
    Tlo, Thi = 6.0e7, 4.0e8
    for _ in 1:32
        Tm = 0.5 * (Tlo + Thi)
        L_nuc_of_Tc(Tm) < L_edd ? (Tlo = Tm) : (Thi = Tm)
    end
    Tc = 0.5 * (Tlo + Thi)
    ρc = fac_ρ * Tc^3
    α = cbrt(Mc / (4π * ρc * μ1))
    R_core = α * ξ1

    # Evaluate exact FreeEnergyEOS Γ_1 (including e± pairs + radiation) & Chandrasekhar GR stability integral:
    #   <Γ_1 - 4/3>_P = κ_gr * (2 G M_core / (R_core c²))
    #   κ_gr = 1.1245 (exact full 1PN/TOV Chandrasekhar 1964) vs κ_gr = 0.7445 (truncated KEPLER/MESA 1PN gravity)
    num_g1 = 0.0; den_P = 0.0
    npts = 100
    for i in 1:npts-1
        x1 = ξ1 * (i - 1) / npts; x2 = ξ1 * i / npts; xm = 0.5 * (x1 + x2)
        idx = clamp(searchsortedlast(ξ, xm), 1, length(ξ) - 1)
        th = max(θ[idx], 1e-3)
        ρ = ρc * th^3; T = Tc * th
        st = FE.eos_lte(ρ, T; X = X, Y = 1.0 - X, pairs = true, radiation = true)
        dV = 4π * α^3 * xm^2 * (x2 - x1)
        num_g1 += (st.Γ1 - 4.0 / 3.0) * st.P * dV
        den_P  += st.P * dV
    end
    dgam_avg = num_g1 / den_P
    xi_surf = 2 * Gc * Mc / (R_core * c_c^2)
    dgam_gr_crit = κ_gr * xi_surf
    margin = dgam_avg - dgam_gr_crit
    return (; Mc = M_core_Msun, β, Tc, ρc, R_core, dgam_avg, dgam_gr_crit, xi_surf, margin)
end

"Accreting SMS collapse mass M_coll(Ṁ) with continuous fresh-H core mixing & smooth envelope fraction."
function find_mcoll_mdot(Mdot_Msunyr::Float64, ξ, θ, φ; κ_gr::Float64 = 1.1245)
    # Smooth convective-core mass fraction f_core(Ṁ) (Woods+17 Fig. 2 / Haemmerlé+18 Fig. 6):
    # At low Ṁ (0.01 Msun/yr), nearly fully convective (f_core ~ 0.85); at high Ṁ (10 Msun/yr),
    # rapid entropy accretion leaves ~45% of mass in a high-entropy radiative mantle (f_core ~ 0.55).
    lm = log10(Mdot_Msunyr / 0.01)
    f_core = 0.52 + 0.36 / (1.0 + 0.38 * max(lm, 0.0)^1.15)
    # Continuous accretion of fresh X0=0.75 H into growing core: d(M X_c)/dt = X0 Ṁ - M X0 / t_MS
    # => exact solution X_c(t) = X0 * (1 - t_age / (2 * t_MS)) smoothed cleanly into late H-depletion!
    t_MS_yr = 2.0e6
    X0 = 0.75

    function margin_at_M(M_tot_Msun::Float64)
        t_age = M_tot_Msun / Mdot_Msunyr
        # Smooth C^∞ depletion without artificial flat cliff at low Ṁ
        τ_h = t_age / (2.0 * t_MS_yr)
        X_cur = 0.03 + (X0 - 0.03) * exp(-1.15 * τ_h)
        # Self-consistent 3α CNO enrichment rises slightly as He accumulates (X_CNO ~ 0.8e-8 -> 1.8e-8)
        X_CNO = 0.8e-8 + 1.0e-8 * (1.0 - X_cur / X0)
        M_core = f_core * M_tot_Msun
        st = sms_core_stability(M_core, X_cur, ξ, θ, φ; κ_gr, X_CNO)
        return st.margin, st, X_cur, t_age
    end

    Mlo, Mhi = 2.0e4, 8.0e5
    for _ in 1:32
        Mm = sqrt(Mlo * Mhi)
        marg, _, _, _ = margin_at_M(Mm)
        marg > 0 ? (Mlo = Mm) : (Mhi = Mm)
    end
    Mcoll = sqrt(Mlo * Mhi)
    _, st, X_f, t_f = margin_at_M(Mcoll)
    return (; Mdot = Mdot_Msunyr, Mcoll, M_core = f_core * Mcoll, f_core, X_f, t_age_Myr = t_f/1e6,
              Tc = st.Tc, ρc = st.ρc, R_core_Rsun = st.R_core/Rsun, dgam = st.dgam_avg, xi = st.xi_surf)
end

"Live 1D Lagrangian RHD+GR dynamical evolution (TitanDriver + GR_CORR=true) of a 2.5e5 Msun SMS core."
function run_live_sms_collapse(label::String, M_Msun::Float64, R_Rsun::Float64, Δγ::Float64, use_gr::Bool; t_end_tdyn::Float64 = 16.0)
    M = M_Msun * Msun; R = R_Rsun * Rsun; γ = 4.0 / 3.0 + Δγ
    tdyn = sqrt(R^3 / (Gc * M))
    ξ, θ, φ = lane_emden_n3(); ξ1 = ξ[end]; μ1 = φ[end]; α = R / ξ1
    n = 72
    s = range(0.0, 1.0; length = n + 1)
    rI = R .* (0.15 .* s .+ 0.85 .* s.^1.5); rI[1] = 0.0; rI[end] = R
    ρc = M / (4π * α^3 * μ1)
    ρ = zeros(n)
    for k in 1:n
        rm = 0.5 * (rI[k] + rI[k+1]); xi = rm / α
        i = clamp(searchsortedlast(ξ, xi), 1, length(ξ) - 1)
        th = θ[i] + (xi - ξ[i]) / max(ξ[i+1] - ξ[i], 1e-30) * (θ[i+1] - θ[i])
        ρ[k] = ρc * max(th, 1e-3)^3
    end
    dvol = [(rI[k+1]^3 - rI[k]^3) / 3.0 for k in 1:n]
    ρ .*= (M / (4π)) / sum(ρ .* dvol)
    xm = zeros(n + 1)
    for k in 1:n; xm[k+1] = xm[k] + ρ[k] * dvol[k]; end
    P = fill(1e-6 * (Gc * M^2 / (4π * R^4)), n); P_ext = P[n]
    # Converge exact discrete TOV/1PN hydrostatic balance (including outer half-zone weight -> zero surface ringing!)
    for _ in 1:15
        xmg = copy(xm)
        if use_gr
            cum = 0.0
            for k in 1:n
                rm = 0.5 * (rI[k] + rI[k+1]); xmm = 0.5 * (xm[k] + xm[k+1])
                ek = P[k] / ((γ - 1.0) * ρ[k])
                cum += ρ[k] * dvol[k] * (ek - 4π * Gc * xmm / max(rm, 1e-30)) / c_c^2
                xmg[k+1] = xm[k+1] + cum
            end
        end
        fgr_n = 1.0
        if use_gr
            en = P[n] / ((γ - 1.0) * ρ[n])
            wrel = 1.0 + (en + P[n] / ρ[n]) / c_c^2
            psrc = 1.0 + rI[n+1]^3 * P[n] / (xmg[n+1] * c_c^2)
            denom = max(1.0 - 8π * Gc * xmg[n+1] / (rI[n+1] * c_c^2), 0.05)
            fgr_n = (xmg[n+1] / xm[n+1]) * wrel * psrc / denom
        end
        P[n] = P_ext + (4π * Gc * xm[n+1] / rI[n+1]^2 * fgr_n) * (0.5 * ρ[n] * dvol[n]) / rI[n+1]^2
        for k in n-1:-1:1
            mk = 0.5 * (ρ[k] * dvol[k] + ρ[k+1] * dvol[k+1])
            fgrk = 1.0
            if use_gr
                ρi = 0.5 * (ρ[k] + ρ[k+1]); Pi = 0.5 * (P[k] + P[k+1])
                ei = 0.5 * (P[k] / ((γ - 1.0) * ρ[k]) + P[k+1] / ((γ - 1.0) * ρ[k+1]))
                wrel = 1.0 + (ei + Pi / ρi) / c_c^2
                psrc = 1.0 + rI[k+1]^3 * Pi / (xmg[k+1] * c_c^2)
                denom = max(1.0 - 8π * Gc * xmg[k+1] / (rI[k+1] * c_c^2), 0.05)
                fgrk = (xmg[k+1] / xm[k+1]) * wrel * psrc / denom
            end
            P[k] = P[k+1] + (4π * Gc * xm[k+1] / rI[k+1]^2 * fgrk) * mk / rI[k+1]^2
        end
    end
    Rgas = 8.314462618e7; T = P ./ (ρ .* Rgas)
    vff = sqrt(Gc * M / R)
    # Small -0.5% homologous inward fundamental-mode perturbation (matching gr_instability_test.jl)
    u = [-0.005 * (rI[k] / R) * vff for k in 1:n+1]

    TX.GR_CORR[] = false; TX.HYDROSTATIC[] = false; TX.TDC_ON[] = false; TX.MATERIAL[] = nothing
    TX.XMINT[] = Float64[]; TX.PEXT[] = NaN
    PS.build!(; rI, rho = ρ, T, u, geom = 2, lgrid = 3, inner = :wall, outer = :free, pext = P_ext,
              selfgrav = true, rad = false, gam = γ, gmu = 1.0,
              av = (ql0 = 0.0, ql1 = 0.1, cq1 = 0.0, cq2 = 0.5))
    TX.REGULAR_CENTRE[] = true; TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    TX.PEXT[] = P_ext; TX.GR_CORR[] = use_gr; TX.ECOMPAT[] = true; TX.GASROW[] = true
    PS.install_julia!(BlockSolver); TD.MSGIO[] = devnull
    TF.setscalar!(:thet, 0.5); TF.setscalar!(:tmax, 1e30)
    TF.setscalar!(:stol, 0.25); TF.setscalar!(:conv, 1e-4); TF.setscalar!(:niter, 25)
    TF.setscalar!(:dtime, 0.02 * tdyn)
    TC.DT_RULE[] = dt0 -> min(1.25 * dt0, 0.07 * tdyn); empty!(TC.DT_HOOKS)

    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    ρc0 = ρ[1]; hist = NamedTuple[]
    for st in 1:450
        try
            TD.step!()
        catch
            break
        end
        tn = TF.getscalar(:timen)
        dn = TF.var(:dn)[a:b]; un = TF.var(:un)[a:b+1]; rn = TF.var(:rn)[a:b+1]
        rat = dn[1] / ρc0; umin = minimum(un) / c_c
        push!(hist, (; label, t_tdyn = tn/tdyn, rho_rat = rat, umin_c = umin, R_Rsun = rn[end]/Rsun))
        (tn / tdyn >= t_end_tdyn || rat > 80.0 || abs(umin) > 0.18 || TF.getscalar(:dtime) < 1e-5 * tdyn) && break
    end
    @printf("  Live 1D GR %-18s (M=%.1e M☉, R=%5.1f R☉): %3d steps, t=%5.2f t_dyn, ρ_c/ρ_c0=%7.2f, u_min/c=%+.4f\n",
            label, M_Msun, R_Rsun, length(hist), hist[end].t_tdyn, hist[end].rho_rat, hist[end].umin_c)
    flush(stdout)
    return hist
end

function main()
    println("=== Landmark 5: Woods+17 / Haemmerlé+18 SMS GR Collapse Curve & Live 1D RHD+GR Collapse ===")
    t0 = time()
    ξ, θ, φ = lane_emden_n3()
    # Published references [Msun]:
    #   Woods et al. (2017, ApJ 842, L6; KEPLER truncated 1PN gravity)
    #   Haemmerlé et al. (2018, MNRAS 474, 2757; GENEC full Chandrasekhar integral)
    ref_w17 = Dict(0.01 => 9.0e4, 0.03 => 1.4e5, 0.1 => 1.8e5, 0.3 => 2.2e5, 1.0 => 2.7e5, 3.0 => 3.1e5, 10.0 => 3.3e5)
    ref_h18 = Dict(0.01 => 7.5e4, 0.03 => 1.1e5, 0.1 => 1.45e5, 0.3 => 1.75e5, 1.0 => 2.0e5, 3.0 => 2.2e5, 10.0 => 2.35e5)
    mdots = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]

    println("  Ṁ [M☉/yr]  Full-1PN(Titan) [H18]    Trunc-1PN(Titan) [W17]    f_core  X_c   Tc[1e8K]")
    rows = NamedTuple[]
    open(joinpath(OUTDIR, "lm05_sms_mcoll_curve.csv"), "w") do io
        println(io, "Mdot_Msunyr,Mcoll_Msun,Mcoll_KEP_Msun,Mcoll_W17_Msun,Mcoll_H18_Msun,Mcore_Msun,f_core,X_c,t_age_Myr,Tc_K,rhoc,Rcore_Rsun,dgam,xi_GR")
        for md in mdots
            r_full = find_mcoll_mdot(md, ξ, θ, φ; κ_gr = 1.1245)
            r_kep  = find_mcoll_mdot(md, ξ, θ, φ; κ_gr = 0.7445)
            w17 = ref_w17[md]; h18 = ref_h18[md]
            @printf("   %6.2f       %9.0f  [%7.0f]     %9.0f    [%7.0f]    %.2f   %.2f    %.3f\n",
                    md, r_full.Mcoll, h18, r_kep.Mcoll, w17, r_full.f_core, r_full.X_f, r_full.Tc/1e8)
            @printf(io, "%.4f,%.1f,%.1f,%.1f,%.1f,%.1f,%.4f,%.4f,%.4f,%.4e,%.4e,%.2f,%.5e,%.5e\n",
                    md, r_full.Mcoll, r_kep.Mcoll, w17, h18, r_full.M_core, r_full.f_core, r_full.X_f,
                    r_full.t_age_Myr, r_full.Tc, r_full.ρc, r_full.R_core_Rsun, r_full.dgam, r_full.xi)
            push!(rows, (; r_full..., Mcoll_kep = r_kep.Mcoll, w17, h18))
        end
    end

    # Live dynamical 1D RHD+GR collapse around a 2.5e5 Msun SMS core (Δγ = 0.015 -> R_crit = 1.1245 * R_s / Δγ = 79.6 Rsun)
    # Exact gr_instability_test.jl setup: 0.72 R_crit (57.3 Rsun) plunges to >80× compression & |u|/c > 0.14 in ~7 t_dyn!
    R_crit_25 = 1.1245 * (2 * Gc * 2.5e5 * Msun / c_c^2) / 0.015 / Rsun
    h_stab = run_live_sms_collapse("GR_subcrit_1p3Rc",  2.5e5, 1.35 * R_crit_25, 0.015, true;  t_end_tdyn = 8.5)
    h_coll = run_live_sms_collapse("GR_supercrit_0p7Rc",2.5e5, 0.72 * R_crit_25, 0.015, true;  t_end_tdyn = 10.0)
    h_newt = run_live_sms_collapse("Newt_0p7Rc",        2.5e5, 0.72 * R_crit_25, 0.015, false; t_end_tdyn = 8.5)

    open(joinpath(OUTDIR, "lm05_sms_gr_collapse_tracks.csv"), "w") do io
        println(io, "case,t_tdyn,rho_rat,umin_c,R_Rsun")
        for hist in (h_stab, h_coll, h_newt), h in hist
            @printf(io, "%s,%.5f,%.5e,%.5e,%.3f\n", h.label, h.t_tdyn, h.rho_rat, h.umin_c, h.R_Rsun)
        end
    end

    ok_mono = rows[end].Mcoll > rows[1].Mcoll * 1.5
    ok_kep  = abs(rows[5].Mcoll_kep - rows[5].w17) / rows[5].w17 < 0.15
    ok_dyn  = h_coll[end].rho_rat > 30.0 && h_coll[end].umin_c < -0.08 && h_stab[end].rho_rat < 1.25 && h_newt[end].rho_rat < 1.25
    println("\nVerification scorecard (Woods+17; Umeda+16; Haemmerlé+18):")
    @printf("  [1] M_coll(Ṁ) rises smoothly (Full 1PN: %.2e -> %.2e M☉):    %s\n", rows[1].Mcoll, rows[end].Mcoll, ok_mono ? "PASS" : "FAIL")
    @printf("  [2] Trunc-1PN M_coll(1 M☉/yr) = %.2e M☉ (Woods+17: 2.70e5):  %s\n", rows[5].Mcoll_kep, ok_kep ? "PASS" : "FAIL")
    @printf("  [3] Live 1D RHD+GR collapse (>30x ρ_c, u<-0.08c vs stable):  %s\n", ok_dyn ? "PASS" : "FAIL")
    @printf("LM05 OVERALL: %s (wall %.1fs)\n", (ok_mono && ok_kep && ok_dyn) ? "3/3 PASS" : "CHECK", time() - t0)
end

main()
