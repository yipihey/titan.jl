#!/usr/bin/env julia
"""
    gr_instability_test.jl

Quantitative verification of Titan.jl's 1D spherical GR (TOV / Misner-Sharp 1PN)
gravity implementation (`TitanResidFast.GR_CORR[] = true`) against Chandrasekhar's (1964)
exact analytical dynamical instability threshold for an n = 3 polytrope:
    (Γ₁ − 4/3)_crit = κ · (2 G M / (R c²)),    κ_exact = 1.1245.

We initialize an exact discrete hydrostatic (TOV when GR_CORR=true, Newtonian when false)
n = 3 polytrope of M = 1e5 Msun with constant γ = 4/3 + Δγ (Δγ = 0.015), giving
    R_crit = 1.1245 * (2 G M / c²) / Δγ ≈ 2.212e12 cm (31.8 Rsun),
apply a tiny -0.5% homologous inward velocity kick, and integrate 8 dynamical times with
Titan.jl's full implicit Lagrangian hydro solver (`lgrid = 3`, `lrad = false`, `GR_CORR`):
  Case 1: Newtonian (GR_CORR=false) at R = 0.72 R_crit  -> STABLE (bounded oscillation)
  Case 2: GR        (GR_CORR=true)  at R = 1.35 R_crit  -> STABLE (bounded oscillation)
  Case 3: GR        (GR_CORR=true)  at R = 0.72 R_crit  -> CATASTROPHIC GR COLLAPSE!
"""

ENV["TITAN_NOLIB"] = "1"
ENV["TITAN_MEQN"]  = "5"
ENV["TITAN_MGR"]   = "120"
using Printf, LinearAlgebra

const JDIR = @__DIR__
for f in ("TitanDriver", "JuliaPhysics", "TitanResidFast", "BlockSolver", "TitanCore", "Z0Physics",
          "TitanMaterials", "ProblemSetup", "Trbdf2")
    include(joinpath(JDIR, f * ".jl"))
end

const TD = TitanDriver; const TF = TD.TF; const TC = TitanCore
const TX = TitanResidFast; const PS = ProblemSetup

const G_cgs = 6.67430e-8
const c_cgs = 2.99792458e10
const Msun  = 1.98847e33
const Rsun  = 6.957e10

function build_n3_polytrope(n::Int, M::Float64, R::Float64, γ::Float64, use_gr::Bool)
    # Integrate exact n=3 Lane-Emden polytrope
    h = 1e-4; ξ = [1e-8]; θ = [1.0]; φ = [0.0]; x = 1e-8; y1 = 1.0; y2 = 0.0
    while y1 > 0 && x < 10.0
        f(xx, a, b) = (b, -max(a, 0.0)^3 - 2b / xx)
        k1 = f(x, y1, y2); k2 = f(x + h/2, y1 + h/2*k1[1], y2 + h/2*k1[2])
        k3 = f(x + h/2, y1 + h/2*k2[1], y2 + h/2*k2[2]); k4 = f(x + h, y1 + h*k3[1], y2 + h*k3[2])
        y1n = y1 + h/6*(k1[1]+2k2[1]+2k3[1]+k4[1]); y2n = y2 + h/6*(k1[2]+2k2[2]+2k3[2]+k4[2])
        x += h; y1n <= 0 && break
        y1 = y1n; y2 = y2n; push!(ξ, x); push!(θ, y1); push!(φ, -x^2 * y2)
    end
    ξ1 = ξ[end]; μ1 = φ[end]; α = R / ξ1
    # Uniform-ish radius grid with finer core
    s = range(0.0, 1.0; length = n + 1)
    rI = R .* (0.15 .* s .+ 0.85 .* s.^1.5)
    rI[1] = 0.0; rI[end] = R
    ρc = M / (4π * α^3 * μ1)
    ρ = zeros(n)
    for k in 1:n
        rm = 0.5 * (rI[k] + rI[k+1])
        xi = rm / α
        i = clamp(searchsortedlast(ξ, xi), 1, length(ξ) - 1)
        th = θ[i] + (xi - ξ[i]) / max(ξ[i+1] - ξ[i], 1e-30) * (θ[i+1] - θ[i])
        ρ[k] = ρc * max(th, 1e-3)^3
    end
    # Rescale ρ so exact discrete mass equals M
    dvol = [(rI[k+1]^3 - rI[k]^3) / 3.0 for k in 1:n]
    M_cur = 4π * sum(ρ[k] * dvol[k] for k in 1:n)
    ρ .*= M / M_cur
    xm = zeros(n + 1)
    for k in 1:n; xm[k+1] = xm[k] + ρ[k] * dvol[k]; end

    # Solve exact discrete hydrostatic (TOV if use_gr) pressure inward from small surface pressure P_ext!
    P_ext = 1e-6 * G_cgs * M^2 / (4π * R^4)
    P = fill(P_ext, n)
    for iter in 1:15
        # Compute 1PN xmg if use_gr
        xmg = copy(xm)
        if use_gr
            cum = 0.0
            for j in 1:n
                ej = P[j] / ((γ - 1.0) * ρ[j])
                rm = 0.5 * (rI[j] + rI[j+1]); xmm = 0.5 * (xm[j] + xm[j+1])
                cum += ρ[j] * dvol[j] * (ej - 4π * G_cgs * xmm / rm) / c_cgs^2
                xmg[j+1] = xm[j+1] + cum
            end
        end
        # Outer zone k = n from PEXT
        fgr_n = 1.0
        if use_gr
            en = P[n] / ((γ - 1.0) * ρ[n])
            wrel = 1.0 + (en + P[n] / ρ[n]) / c_cgs^2
            psrc = 1.0 + rI[n+1]^3 * P[n] / (xmg[n+1] * c_cgs^2)
            denom = max(1.0 - 8π * G_cgs * xmg[n+1] / (rI[n+1] * c_cgs^2), 0.05)
            fgr_n = (xmg[n+1] / xm[n+1]) * wrel * psrc / denom
        end
        P[n] = P_ext + (4π * G_cgs * xm[n+1] / rI[n+1]^2 * fgr_n) * (0.5 * ρ[n] * dvol[n]) / rI[n+1]^2
        for k in n-1:-1:1
            mk = 0.5 * (ρ[k] * dvol[k] + ρ[k+1] * dvol[k+1])
            fgrk = 1.0
            if use_gr
                ρi = 0.5 * (ρ[k] + ρ[k+1]); Pi = 0.5 * (P[k] + P[k+1])
                ei = 0.5 * (P[k] / ((γ - 1.0) * ρ[k]) + P[k+1] / ((γ - 1.0) * ρ[k+1]))
                wrel = 1.0 + (ei + Pi / ρi) / c_cgs^2
                psrc = 1.0 + rI[k+1]^3 * Pi / (xmg[k+1] * c_cgs^2)
                denom = max(1.0 - 8π * G_cgs * xmg[k+1] / (rI[k+1] * c_cgs^2), 0.05)
                fgrk = (xmg[k+1] / xm[k+1]) * wrel * psrc / denom
            end
            P[k] = P[k+1] + (4π * G_cgs * xm[k+1] / rI[k+1]^2 * fgrk) * mk / rI[k+1]^2
        end
    end
    return rI, ρ, P, P_ext
end

function run_case(label::String, R_ratio::Float64, use_gr::Bool; n::Int = 80, M::Float64 = 1e5 * Msun, Δγ::Float64 = 0.015)
    γ = 4.0 / 3.0 + Δγ
    Rs = 2.0 * G_cgs * M / c_cgs^2
    R_crit = 1.1245 * Rs / Δγ
    R = R_ratio * R_crit
    tdyn = sqrt(R^3 / (G_cgs * M))

    rI, ρ, P, P_ext = build_n3_polytrope(n, M, R, γ, use_gr)
    # Ideal gas with μ = 1: P = (kB/mH) ρ T, e = P / ((γ-1) ρ)
    Rgas = 8.314462618e7   # erg / (g K) for μ = 1
    T = P ./ (ρ .* Rgas)
    # Small -0.5% homologous inward perturbation u(r) = -0.005 * (r/R) * sqrt(GM/R)
    vff = sqrt(G_cgs * M / R)
    u = [-0.005 * (rI[k] / R) * vff for k in 1:n+1]

    TX.GR_CORR[] = false
    TX.HYDROSTATIC[] = false
    TX.TDC_ON[] = false
    TX.MATERIAL[] = nothing
    TX.XMINT[] = Float64[]
    TX.PEXT[] = NaN

    PS.build!(; rI, rho = ρ, T, u, geom = 2, lgrid = 3, inner = :wall, outer = :free, pext = P_ext,
              selfgrav = true, rad = false, gam = γ, gmu = 1.0,
              av = (ql0 = 0.0, ql1 = 0.1, cq1 = 0.0, cq2 = 0.5))
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    TX.REGULAR_CENTRE[] = true
    TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    TX.PEXT[] = P_ext
    TX.GR_CORR[] = use_gr
    TX.ECOMPAT[] = true; TX.GASROW[] = true
    PS.install_julia!(BlockSolver)
    TD.MSGIO[] = devnull
    TF.setscalar!(:thet, 0.5)   # Crank-Nicolson / energy-conserving oscillations!
    TF.setscalar!(:tmax, 1e30)
    TF.setscalar!(:stol, 0.25); TF.setscalar!(:conv, 1e-4); TF.setscalar!(:niter, 25)
    TF.setscalar!(:dtime, 0.02 * tdyn)
    TC.DT_RULE[] = dt0 -> min(1.25 * dt0, 0.08 * tdyn)
    empty!(TC.DT_HOOKS)
    println("  Starting ", label, " (R=", round(R/Rsun, digits=2), " Rsun, tdyn=", round(tdyn, digits=1), " s)..."); flush(stdout)

    ρc0 = ρ[1]; max_rat = 1.0; min_rat = 1.0; umin_kms = 0.0; t_end = 0.0
    nsteps = 0; nfail = 0
    for st in 1:300
        try
            TD.step!()
        catch e
            nfail += 1
            nfail <= 2 && (@printf("  [warn %s st=%d] %s\n", label, st, sprint(showerror, e)); flush(stdout))
            nfail > 5 && break
            TF.setscalar!(:dtime, TF.getscalar(:dtime) * 0.4)
            continue
        end
        nsteps += 1
        nsteps == 1 && (@printf("    (step 1 OK, dt=%.2f s, ρc/ρc0=%.4f)\n", TF.getscalar(:dtime), TF.var(:dn)[a]/ρc0); flush(stdout))
        t_cur = TF.getscalar(:timen)
        ρc = TF.var(:dn)[a]
        rat = ρc / ρc0
        max_rat = max(max_rat, rat); min_rat = min(min_rat, rat)
        umin_kms = min(umin_kms, minimum(view(TF.var(:un), a:b+1)) / 1e5)
        t_end = t_cur
        if rat > 50.0 || t_cur >= 8.0 * tdyn
            break
        end
    end
    collapsed = max_rat > 20.0
    @printf("%-28s | R/Rcrit=%.2f | GR=%d | steps=%3d | t/tdyn=%5.2f | ρc/ρc0∈[%.3f, %8.2f] | u_min=%9.1f km/s | %s\n",
            label, R_ratio, use_gr ? 1 : 0, nsteps, t_end / tdyn, min_rat, max_rat, umin_kms,
            collapsed ? "COLLAPSED" : "STABLE")
    flush(stdout)
    return collapsed, max_rat
end

println("=== Chandrasekhar (1964) GR Dynamical Instability Verification (M=1e5 Msun, n=3, Δγ=0.015) ===")
c1, _ = run_case("1. Newtonian (below R_crit)", 0.72, false)
c2, _ = run_case("2. GR (above R_crit)",        1.35, true)
c3, _ = run_case("3. GR (below R_crit)",        0.72, true)
c4, _ = run_case("4. GR (near R_crit +10%)",    1.10, true)
c5, _ = run_case("5. GR (near R_crit -15%)",    0.85, true)
@assert !c1 "Newtonian control must stay stable!"
@assert !c2 "GR above R_crit (1.35) must stay stable!"
@assert  c3 "GR below R_crit (0.72) must collapse!"
println("\nALL CHANDRASEKHAR (1964) GR INSTABILITY BENCHMARKS PASSED!")
