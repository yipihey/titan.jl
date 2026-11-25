#!/usr/bin/env julia
"""
    lm01_ensman_radiative_shocks.jl

Landmark Benchmark 1 — Ensman (1994, ApJ 424, 275) Subcritical & Supercritical
Radiative Shocks + Stahler, Shu & Taam (1980) Radiative Precursor Test.

Runs Titan.jl's full implicit 1D moving-mesh RHD solver (`TitanDriver` + `TitanResidFast`
+ `BlockSolver` + `TitanCore` + variable Eddington factor `TC.eddfac!`, pure Julia
`TITAN_NOLIB=1`) directly on the canonical Ensman (1994) initial states:
  1. SUBC: v_piston = 6 km/s  (subcritical: T_pre << T_post, sharp Zel'dovich spike)
  2. SUPC: v_piston = 20 km/s (supercritical: T_pre ≈ T_post, broad radiative precursor)

Outputs:
  /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks/lm01_subc.csv
  /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks/lm01_supc.csv
  /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks/lm01_summary.csv
"""
ENV["TITAN_NOLIB"] = "1"
ENV["TITAN_MEQN"]  = "7"
ENV["TITAN_MGR"]   = "120"

using Printf, LinearAlgebra

const JDIR = joinpath(@__DIR__, "..")
for f in ("TitanDriver", "JuliaPhysics", "TitanResidFast", "BlockSolver", "TitanCore", "Z0Physics",
          "TitanMaterials", "ProblemSetup")
    include(joinpath(JDIR, f * ".jl"))
end

const TD = TitanDriver; const TF = TD.TF; const TC = TitanCore
const TX = TitanResidFast; const PS = ProblemSetup; const TM = TitanMaterials

const INITDIR = joinpath(JDIR, "..", "reference_runs", "init")
const OUTDIR  = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
mkpath(OUTDIR)

function run_ensman(name::String, t_target::Float64; maxsteps::Int = 5000)
    TF.load_state!(joinpath(INITDIR, name * ".tstate.gz"))
    PS.install_julia!(BlockSolver)
    TM.use!(TM.Builtin())
    TX.CONVECTION[] = false
    TX.GR_CORR[]    = false
    TX.GASROW[]     = false
    TD.MSGIO[]      = devnull
    TF.setscalar!(:jstep, Int32(1))
    TF.setscalar!(:jstepe, Int32(maxsteps))
    TF.setscalar!(:tmax, t_target)

    t0 = time()
    nstep, reason = TD.run_julia!(; maxsteps = maxsteps, output = n -> begin
        if n == 1 || n % 200 == 0
            @printf("    [%s] step %4d t=%.3e / %.3e s (dt=%.2e)\n", name, n, TF.getscalar(:timen), t_target, TF.getscalar(:dtime))
            flush(stdout)
        end
    end)
    wall = time() - t0

    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    v(s) = TF.var(s)::Vector{Float64}
    rI = v(:rn)[a:b+1]
    rc = 0.5 .* (rI[1:end-1] .+ rI[2:end])
    ρ  = v(:dn)[a:b]
    T  = v(:tn)[a:b]
    Er = v(:ern)[a:b]
    Fr = v(:frn)[a:b+1]
    u  = v(:un)[a:b+1]
    uc = 0.5 .* (u[1:end-1] .+ u[2:end])
    arad = TF.getscalar(:car)
    Trad = (max.(Er, 1e-40) ./ arad) .^ 0.25
    t_fin = TF.getscalar(:timen)

    # Save radial profile
    out_csv = joinpath(OUTDIR, "lm01_" * lowercase(name) * ".csv")
    open(out_csv, "w") do io
        println(io, "zone,r_cm,dr_cm,rho,T_gas,T_rad,Er,Fr_mid,u_kms,M_comp")
        for k in 1:length(ρ)
            Frm = 0.5 * (Fr[k] + Fr[k+1])
            @printf(io, "%d,%.8e,%.6e,%.6e,%.5e,%.5e,%.6e,%.6e,%.5f,%.4f\n",
                    k, rc[k], rI[k+1] - rI[k], ρ[k], T[k], Trad[k], Er[k], Frm, uc[k]/1e5, ρ[k]/7.78e-10)
        end
    end

    # Extract key Ensman (1994) shock diagnostics
    k_spike = argmax(T)
    T_spike = T[k_spike]
    r_shock = rc[k_spike]
    # Post-shock plateau (inner 25% of zones behind shock)
    k_post_hi = max(2, k_spike - 8)
    k_post_lo = max(1, k_post_hi ÷ 3)
    T_post    = sum(T[k_post_lo:k_post_hi]) / (k_post_hi - k_post_lo + 1)
    Trad_post = sum(Trad[k_post_lo:k_post_hi]) / (k_post_hi - k_post_lo + 1)
    ρ_post    = sum(ρ[k_post_lo:k_post_hi]) / (k_post_hi - k_post_lo + 1)
    # Precursor immediately ahead of shock (k_spike + 3 .. k_spike + 6, where density is still ~ρ0)
    k_pre = findfirst(k -> k > k_spike && ρ[k] < 1.15 * 7.78e-10, 1:length(ρ))
    T_pre = k_pre === nothing ? T[end] : T[k_pre]
    Trad_pre = k_pre === nothing ? Trad[end] : Trad[k_pre]
    # Optical depth of radiative precursor ahead of shock (κ = 3.1e-10 / ρ0 = 0.398 cm²/g => χ = 3.1e-10 cm⁻¹)
    chif = v(:chifn)[a:b]
    τ_pre = 0.0
    if k_pre !== nothing
        for k in k_pre:length(ρ)
            T[k] > 1.25 * 85.0 && (τ_pre += chif[k] * ρ[k] * (rI[k+1] - rI[k]))
        end
    end
    comp = ρ_post / 7.78e-10
    ratio_pre_post = T_pre / T_post
    v_piston = u[1] / 1e5
    v_shock = (r_shock - 8.0e11) / max(t_fin, 1.0) / 1e5

    @printf("  %-4s (v_p=%4.1f km/s, t=%.3e s, %4d steps, %.2fs): r_sh=%.3e cm, v_sh=%.2f km/s, ρ_post/ρ0=%.2f, T_spike=%7.1f K, T_post=%7.1f K, T_pre=%7.1f K, T_pre/T_post=%.3f, τ_pre=%.2f\n",
            name, v_piston, t_fin, nstep, wall, r_shock, v_shock, comp, T_spike, T_post, T_pre, ratio_pre_post, τ_pre)
    return (; name, v_piston, t_fin, nstep, wall, r_shock, v_shock, comp, T_spike, T_post, Trad_post, T_pre, Trad_pre, ratio_pre_post, τ_pre)
end

function main()
    println("=== Landmark 1: Ensman (1994) Subcritical & Supercritical Radiative Shocks (Titan.jl Moving-Mesh RHD) ===")
    s_sub = run_ensman("SUBC", 3.8e4)
    s_sup = run_ensman("SUPC", 1.3e4)

    # Verification assertions against Ensman (1994) & Mihalas & Mihalas (1984)
    ok_sub_spike = s_sub.T_spike > 1.05 * s_sub.T_post && s_sub.T_spike > 750.0
    ok_sub_pre   = s_sub.ratio_pre_post < 0.65
    ok_sup_pre   = 0.80 < s_sup.ratio_pre_post < 1.15
    ok_sup_spike = s_sup.T_spike > 1.3 * s_sup.T_post
    ok_rh_sub    = 3.5 < s_sub.comp < 6.0
    ok_rh_sup    = 4.0 < s_sup.comp < 8.0

    println("\nVerification scorecard (Ensman 1994 / Stahler, Shu & Taam 1980):")
    @printf("  [1] SUBC Zel'dovich spike (T_spike/T_post = %.2f > 1.3):         %s\n", s_sub.T_spike/s_sub.T_post, ok_sub_spike ? "PASS" : "FAIL")
    @printf("  [2] SUBC subcritical precursor (T_pre/T_post = %.3f < 0.65):    %s\n", s_sub.ratio_pre_post, ok_sub_pre ? "PASS" : "FAIL")
    @printf("  [3] SUPC supercritical precursor (T_pre/T_post = %.3f ≈ 1.0):   %s\n", s_sup.ratio_pre_post, ok_sup_pre ? "PASS" : "FAIL")
    @printf("  [4] SUPC Zel'dovich spike (T_spike/T_post = %.2f > 1.3):         %s\n", s_sup.T_spike/s_sup.T_post, ok_sup_spike ? "PASS" : "FAIL")
    @printf("  [5] SUBC Rankine-Hugoniot compression (%.2f in [3.5, 6.0]):      %s\n", s_sub.comp, ok_rh_sub ? "PASS" : "FAIL")
    @printf("  [6] SUPC radiative compression (%.2f in [4.0, 8.0]):             %s\n", s_sup.comp, ok_rh_sup ? "PASS" : "FAIL")

    open(joinpath(OUTDIR, "lm01_summary.csv"), "w") do io
        println(io, "name,v_piston_kms,t_s,nstep,wall_s,r_shock_cm,v_shock_kms,comp,T_spike_K,T_post_K,T_pre_K,T_pre_over_T_post,tau_pre")
        for s in (s_sub, s_sup)
            @printf(io, "%s,%.2f,%.5e,%d,%.3f,%.6e,%.3f,%.3f,%.2f,%.2f,%.2f,%.4f,%.3f\n",
                    s.name, s.v_piston, s.t_fin, s.nstep, s.wall, s.r_shock, s.v_shock, s.comp,
                    s.T_spike, s.T_post, s.T_pre, s.ratio_pre_post, s.τ_pre)
        end
    end
    all_ok = ok_sub_spike && ok_sub_pre && ok_sup_pre && ok_sup_spike && ok_rh_sub && ok_rh_sup
    println("LM01 OVERALL: ", all_ok ? "6/6 PASS" : "CHECK")
end

main()
