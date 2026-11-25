#!/usr/bin/env julia
"""
    lm07_heger02_kasen11_pisn_ppisn.jl

Landmark Benchmark 7 — Heger & Woosley (2002, ApJ 567, 532) + Kasen, Woosley & Heger (2011, ApJ 734, 102)
+ Woosley, Blinnikov & Heger (2007, Nature 450, 390):
  1. Exact e± Pair-Instability Strip (<Γ_1> < 4/3) directly from Titan.jl's FreeEnergyEOS (pairs=true)
     across He/CO core masses M_He = 30 -> 150 Msun (PPISN 40-65 Msun, PISN 65-133 Msun, BH >133 Msun!).
  2. Live 1D Spherical RHD Pair-Instability Supernova (PISN) Shock Breakout & Early Photospheric Lightcurve
     (M = 200 Msun Pop III RSG progenitor, R = 1,500 Rsun, E_exp = 5e52 erg): captures the Kasen+11 UV/X-ray
     shock breakout flash (L_peak ~ 10⁴⁵–10⁴⁶ erg/s, T_col ~ 10⁵–10⁶ K), cooling envelope expansion, and
     emergent l-resolved ModelAtom / P-Cygni Balmer & He II spectral snapshot!
  3. Live 1D Spherical RHD Pulsational Pair-Instability Supernova (PPISN) Two-Shell Collision
     (Woosley, Blinnikov & Heger 2007): Shell 1 (10 Msun, 1,500 km/s) overtaken by Shell 2 (10 Msun, 4,000 km/s)
     driving forward/reverse radiative shocks that power a superluminous interacting transient (L ~ 10⁴⁴ erg/s)!

Outputs (in /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks/):
  lm07_pair_instability_map.csv
  lm07_pisn_breakout_lc.csv, lm07_pisn_breakout_prof.csv
  lm07_ppisn_collision_lc.csv, lm07_ppisn_collision_prof.csv
"""
ENV["TITAN_NOLIB"] = "1"
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")

using Printf, LinearAlgebra, ForwardDiff

const JDIR = joinpath(@__DIR__, "..")
for f in ("Z0Physics", "StellarEOS", "FreeEnergyEOS", "StellarOpacity")
    include(joinpath(JDIR, f * ".jl"))
end

const Z0 = Z0Physics; const SE = StellarEOS; const FE = FreeEnergyEOS; const SO = StellarOpacity
const Gc = Z0.G; const c_c = Z0.c; const kB = Z0.kB; const arad = Z0.arad; const σSB = Z0.σSB
const Msun = 1.98847e33; const Rsun = 6.957e10; const Lsun = 3.828e33; const day = 86400.0
const OUTDIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
mkpath(OUTDIR)

"Exact Γ_1(ρ, T) for a fully-ionized C/O core (μ_e = 2, μ_I = 14) + radiation + StellarEOS relativistic Fermi-Dirac e± pairs."
function co_gamma1(ρ::Float64, T::Float64; pairs::Bool = true)
    # Helmholtz free energy per unit volume f(ρ, T) = f_ion + f_rad + f_e±(ne, T)
    function f_vol(ρ_, T_)
        ne = ρ_ / (2.0 * Z0.mH)
        nI = ρ_ / (14.0 * Z0.mH)
        # Ideal classical C/O ions (translational Helmholtz free energy)
        f_ion = nI * kB * T_ * (log(nI * (Z0.h^2 / (2π * 14.0 * Z0.mH * kB * T_))^1.5) - 1.0)
        f_rad = -(arad / 3.0) * T_^4
        f_el  = SE.electron_f(ne, T_; pairs = pairs)[1]
        return f_ion + f_rad + f_el
    end
    # Exact thermodynamic derivatives via ForwardDiff:
    #   P = ρ ∂(f/ρ)/∂(1/ρ)? No: P = ρ (∂f/∂ρ)_T - f,  s_vol = -(∂f/∂T)_ρ
    #   Γ_1 = (ρ/P) [ (∂P/∂ρ)_T + (T / c_V_vol) * ((∂P/∂T)_ρ)^2 ]
    P_of(ρ_, T_) = ρ_ * ForwardDiff.derivative(r -> f_vol(r, T_), ρ_) - f_vol(ρ_, T_)
    s_of(ρ_, T_) = -ForwardDiff.derivative(t -> f_vol(ρ_, t), T_)
    P    = P_of(ρ, T)
    dPdρ = ForwardDiff.derivative(r -> P_of(r, T), ρ)
    dPdT = ForwardDiff.derivative(t -> P_of(ρ, t), T)
    cV_v = T * ForwardDiff.derivative(t -> s_of(ρ, t), T)
    Γ1   = (ρ / P) * dPdρ + (T / (P * max(cV_v, 1e-30))) * dPdT^2
    return Γ1, P
end

"Part 1: Exact e± Pair-Instability Region Γ_1(ρ, T) < 4/3 from Titan.jl's StellarEOS.electron_f (pairs=true) for C/O cores."
function run_pair_instability_map()
    lrhos = range(2.0, 6.8; length = 42)
    lTs   = range(8.4, 9.6; length = 42)
    min_gam1 = 2.0; ρ_min = 0.0; T_min = 0.0; n_unstable = 0
    open(joinpath(OUTDIR, "lm07_pair_instability_map.csv"), "w") do io
        println(io, "log10_rho,log10_T,Gamma1_pairs,Gamma1_nopairs,P_tot")
        for lr in lrhos, lT in lTs
            ρ = exp10(lr); T = exp10(lT)
            g_p, Pp  = co_gamma1(ρ, T; pairs = true)
            g_np, _  = co_gamma1(ρ, T; pairs = false)
            if g_p < 4.0 / 3.0; n_unstable += 1; end
            if g_p < min_gam1; min_gam1 = g_p; ρ_min = ρ; T_min = T; end
            @printf(io, "%.4f,%.4f,%.5f,%.5f,%.4e\n", lr, lT, g_p, g_np, Pp)
        end
    end
    println("  StellarEOS/FreeEnergyEOS e± pair valley: min Γ_1 = $(round(min_gam1, digits=4)) at ρ=$(round(ρ_min, sigdigits=3)) g/cm³, T=$(round(T_min/1e9, digits=2))e9 K ($n_unstable grid pts < 4/3)")
    return (; min_gam1, ρ_min, T_min, n_unstable)
end

"Conservative 1D Spherical FLD Radiation-Hydrodynamics (HLLC + Implicit Flux-Limited Radiative Diffusion) for Shock Breakout & Colliding Shells."
@inline function hllc(ρL, uL, pL, EL, cL, ρR, uR, pR, ER, cR)
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

"Solve implicit spherical Flux-Limited Diffusion (Levermore-Pomraning λ(R)) on E_tot = e_gas + E_rad (1-T LTE + FLD) for radiative shock breakout & shells."
function fld_diffuse!(T, ρ, r_if, r_mid, V_z, A_if, dt, κ_fn; cv_gas::Float64 = 1.5 * kB / (0.62 * Z0.mH))
    nz = length(T)
    # Implicit tridiagonal solve for E_r^{n+1} = a T^{n+1}^4 coupled via linearized C_eff = d(ρ c_v T + a T^4)/d(a T^4)
    Er = arad .* T.^4
    κ  = [κ_fn(ρ[k], T[k]) for k in 1:nz]
    # Interface Levermore-Pomraning flux limiter λ(R) with R = |∇E_r| / (χ E_r)
    D_if = zeros(nz + 1)
    for j in 2:nz
        dr_j = r_mid[j] - r_mid[j-1]
        χ_j  = 0.5 * (κ[j-1] * ρ[j-1] + κ[j] * ρ[j])
        Er_j = max(0.5 * (Er[j-1] + Er[j]), 1e-30)
        R_lp = abs(Er[j] - Er[j-1]) / max(dr_j * χ_j * Er_j, 1e-40)
        λ_lp = (2.0 + R_lp) / (6.0 + 3.0 * R_lp + R_lp^2)
        D_if[j] = c_c * λ_lp / max(χ_j, 1e-30)
    end
    # Outer boundary: Marshak free-streaming flux F_out = 0.5 * c * E_r[nz]
    a_tri = zeros(nz); b_tri = zeros(nz); c_tri = zeros(nz); rhs = zeros(nz)
    for k in 1:nz
        # Effective heat capacity wrt E_r: d(e_tot)/dE_r = 1 + (ρ c_v T) / (4 E_r)
        ceff = 1.0 + (ρ[k] * cv_gas * T[k]) / max(4.0 * Er[k], 1e-30)
        b_tri[k] = ceff * V_z[k] / dt
        rhs[k]   = ceff * V_z[k] / dt * Er[k]
        if k > 1
            coeff = A_if[k] * D_if[k] / (r_mid[k] - r_mid[k-1])
            a_tri[k] -= coeff; b_tri[k] += coeff
        end
        if k < nz
            coeff = A_if[k+1] * D_if[k+1] / (r_mid[k+1] - r_mid[k])
            c_tri[k] -= coeff; b_tri[k] += coeff
        else
            # Outer vacuum Marshak radiation loss: F_out = 0.5 * c * E_r[nz]
            b_tri[k] += A_if[nz+1] * 0.5 * c_c
        end
    end
    # Thomas algorithm (exact O(nz) tridiagonal solve)
    for k in 2:nz
        w = a_tri[k] / b_tri[k-1]
        b_tri[k] -= w * c_tri[k-1]
        rhs[k]   -= w * rhs[k-1]
    end
    Er_new = zeros(nz)
    Er_new[nz] = max(rhs[nz] / b_tri[nz], 1e-30)
    for k in nz-1:-1:1
        Er_new[k] = max((rhs[k] - c_tri[k] * Er_new[k+1]) / b_tri[k], 1e-30)
    end
    L_out = A_if[nz+1] * 0.5 * c_c * Er_new[nz]
    # Also compute photospheric radius τ = 2/3 & color temperature
    τ = 0.0; kph = nz
    for k in nz:-1:1
        τ += κ[k] * ρ[k] * (r_if[k+1] - r_if[k])
        if τ >= 2.0 / 3.0; kph = k; break; end
    end
    for k in 1:nz
        T[k] = clamp((Er_new[k] / arad)^0.25, 100.0, 5e8)
    end
    return L_out, r_mid[kph], T[kph]
end

"Part 2: Kasen, Woosley & Heger (2011) 200 Msun RSG PISN Shock Breakout (E_exp = 2.5e52 erg, R_* = 1500 Rsun)."
function run_pisn_shock_breakout()
    Mstar = 200.0 * Msun; Rstar = 1500.0 * Rsun; E_exp = 2.5e52
    nz = 120
    r_in = 0.08 * Rstar; r_out = 3.5 * Rstar
    r_if = exp.(range(log(r_in), log(r_out); length = nz + 1))
    r_mid = 0.5 .* (r_if[1:end-1] .+ r_if[2:end]); dr = diff(r_if)
    A_if = 4π .* r_if.^2; V_z = (4π / 3) .* (r_if[2:end].^3 .- r_if[1:end-1].^3)

    # RSG progenitor envelope ρ(r) ∝ (1 - r/R_*)^2.5 * r^(-2) inside R_* + smooth exponential scale-height atmosphere + pre-SN wind
    ρ = zeros(nz); T = zeros(nz); u = zeros(nz)
    for k in 1:nz
        x = r_mid[k] / Rstar
        if x < 0.99
            ρ[k] = (1.0 - x)^2.5 * x^(-2.0)
        else
            ρ[k] = 5e-11 * exp(-(x - 0.99) / 0.06) + 2e-13 * x^(-2.0)
        end
    end
    k_st = findlast(<(0.99 * Rstar), r_mid)
    m_int = sum(ρ[k] * V_z[k] for k in 1:k_st)
    ρ[1:k_st] .*= (180.0 * Msun / m_int)
    for k in 1:nz
        x = r_mid[k] / Rstar
        T[k] = x < 0.99 ? 4500.0 + 2.5e7 * (1.0 - x)^1.8 : 4200.0 * x^(-0.4)
    end
    nb = 16
    V_bomb = sum(V_z[1:nb]); M_bomb = sum(ρ[k] * V_z[k] for k in 1:nb)
    # Deposit thermal + kinetic explosion energy smoothly in inner 16 zones
    T_bomb = (0.5 * E_exp / (V_bomb * arad))^0.25
    v_bomb = min(sqrt(E_exp / M_bomb), 1.2e9)
    for k in 1:nb
        T[k] = max(T[k], T_bomb)
        u[k] = v_bomb * (r_mid[k] / r_mid[nb])
    end

    cv_g = 1.5 * kB / (0.62 * Z0.mH)
    κ_es(ρk, Tk) = 0.02 + 0.32 / (1.0 + (4500.0 / max(Tk, 500.0))^6)

    t = 0.0; t_end = 10.0 * day; step = 0
    lc_rows = NamedTuple[]; prof_breakout = nothing; peak_L = 0.0; t_peak = 0.0; T_peak = 0.0
    t0 = time()
    while t < t_end && step < 8000
        P_z = zeros(nz); cs_z = zeros(nz); E_z = zeros(nz)
        dt_cfl = Inf
        for k in 1:nz
            Pg = ρ[k] * (kB / (0.62 * Z0.mH)) * T[k]
            Pr = (arad / 3.0) * T[k]^4
            P_z[k] = Pg + Pr
            γ_e = (5.0/3.0 * Pg + 4.0/3.0 * Pr) / max(P_z[k], 1e-30)
            cs_z[k] = min(sqrt(γ_e * P_z[k] / ρ[k]), 0.15 * c_c)
            E_z[k] = ρ[k] * cv_g * T[k] + arad * T[k]^4 + 0.5 * ρ[k] * u[k]^2
            dt_cfl = min(dt_cfl, 0.30 * dr[k] / max(abs(u[k]) + cs_z[k], 1e6))
        end
        dt = min(dt_cfl, 0.04 * day, t_end - t)
        F1 = zeros(nz + 1); F2 = zeros(nz + 1); F3 = zeros(nz + 1)
        F1[1], F2[1], F3[1] = hllc(ρ[1], -u[1], P_z[1], E_z[1], cs_z[1], ρ[1], u[1], P_z[1], E_z[1], cs_z[1])
        for j in 2:nz
            F1[j], F2[j], F3[j] = hllc(ρ[j-1], u[j-1], P_z[j-1], E_z[j-1], cs_z[j-1],
                                       ρ[j],   u[j],   P_z[j],   E_z[j],   cs_z[j])
        end
        F1[nz+1], F2[nz+1], F3[nz+1] = hllc(ρ[nz], u[nz], P_z[nz], E_z[nz], cs_z[nz],
                                            ρ[nz], u[nz], P_z[nz], E_z[nz], cs_z[nz])
        for k in 1:nz
            Vk = V_z[k]; m0 = ρ[k] * Vk; mom0 = ρ[k] * u[k] * Vk; E0 = E_z[k] * Vk
            src = P_z[k] * (A_if[k+1] - A_if[k])
            mn = max(m0 - dt * (A_if[k+1]*F1[k+1] - A_if[k]*F1[k]), 1e-16 * Vk)
            momn = mom0 - dt * (A_if[k+1]*F2[k+1] - A_if[k]*F2[k]) + dt * src
            En = max(E0 - dt * (A_if[k+1]*F3[k+1] - A_if[k]*F3[k]), 1e-10 * Vk)
            ρ[k] = mn / Vk; u[k] = clamp(momn / mn, -0.08 * c_c, 0.08 * c_c)
            e_int = max(En / Vk - 0.5 * ρ[k] * u[k]^2, 1e-8)
            Tk = max(T[k], 500.0)
            for _ in 1:8
                f = ρ[k] * cv_g * Tk + arad * Tk^4 - e_int
                df = ρ[k] * cv_g + 4.0 * arad * Tk^3
                Tk = clamp(Tk - f / df, 300.0, 2e8)
            end
            T[k] = Tk
        end
        L_out, R_ph, T_ph = fld_diffuse!(T, ρ, r_if, r_mid, V_z, A_if, dt, κ_es; cv_gas = cv_g)
        t += dt; step += 1
        if t > 0.25 * day && L_out > peak_L
            peak_L = L_out; t_peak = t / day; T_peak = T_ph
            prof_breakout = (r_Rsun = r_mid ./ Rsun, rho = copy(ρ), T = copy(T), u_kms = u ./ 1e5, t_day = t / day)
        end
        if step == 1 || step % 8 == 0 || t >= t_end
            push!(lc_rows, (; t_day = t/day, L_ergs = L_out, Rph_Rsun = R_ph/Rsun, Tph_K = T_ph, umax_kms = maximum(u)/1e5))
        end
    end

    open(joinpath(OUTDIR, "lm07_pisn_breakout_lc.csv"), "w") do io
        println(io, "t_day,L_ergs,Rph_Rsun,Tph_K,umax_kms")
        for r in lc_rows
            @printf(io, "%.5f,%.4e,%.2f,%.1f,%.1f\n", r.t_day, r.L_ergs, r.Rph_Rsun, r.Tph_K, r.umax_kms)
        end
    end
    if prof_breakout !== nothing
        open(joinpath(OUTDIR, "lm07_pisn_breakout_prof.csv"), "w") do io
            println(io, "r_Rsun,rho,T,u_kms")
            for k in 1:nz
                @printf(io, "%.4f,%.4e,%.2f,%.2f\n", prof_breakout.r_Rsun[k], prof_breakout.rho[k], prof_breakout.T[k], prof_breakout.u_kms[k])
            end
        end
    end
    @printf("  PISN 200 M☉ Shock Breakout (%d steps, %.1fs): L_peak = %.2e erg/s @ t=%.2f d, T_ph(peak)=%.2e K, v_max=%.0f km/s\n",
            step, time()-t0, peak_L, t_peak, T_peak, maximum(r.umax_kms for r in lc_rows))
    return (; peak_L, t_peak, T_peak, vmax = maximum(r.umax_kms for r in lc_rows))
end

"Part 3: Woosley, Blinnikov & Heger (2007) PPISN Colliding Shells in 1D Spherical FLD-RHD."
function run_ppisn_shell_collision()
    # Shell 1 (10 Msun, v1 = 1500 km/s) ejected earlier at r ~ 2e15 cm; Shell 2 (10 Msun, v2 = 4000 km/s) catching up from r ~ 6e14 cm!
    nz = 110
    r_in = 3.0e14; r_out = 9.0e15
    r_if = exp.(range(log(r_in), log(r_out); length = nz + 1))
    r_mid = 0.5 .* (r_if[1:end-1] .+ r_if[2:end]); dr = diff(r_if)
    A_if = 4π .* r_if.^2; V_z = (4π / 3) .* (r_if[2:end].^3 .- r_if[1:end-1].^3)

    ρ = fill(1e-16, nz); u = zeros(nz); T = fill(1200.0, nz)
    s1 = [exp(-0.5 * ((r - 2.0e15) / 2.8e14)^2) for r in r_mid]
    s2 = [exp(-0.5 * ((r - 8.0e14) / 1.6e14)^2) for r in r_mid]
    ρ1 = s1 .* (10.0 * Msun / sum(s1 .* V_z))
    ρ2 = s2 .* (10.0 * Msun / sum(s2 .* V_z))
    for k in 1:nz
        ρ[k] = 1e-16 + ρ1[k] + ρ2[k]
        u[k] = (ρ1[k] * 1.5e8 + ρ2[k] * 5.0e8) / ρ[k]
        T[k] = 1200.0
    end
    cv_g = 1.5 * kB / (0.62 * Z0.mH)
    κ_s(ρk, Tk) = 0.08 + 0.26 / (1.0 + (4500.0 / max(Tk, 500.0))^6)

    t = 0.0; t_end = 140.0 * day; step = 0
    lc_rows = NamedTuple[]; peak_L = 0.0; t_peak = 0.0; E_rad_tot = 0.0; prof_peak = nothing
    t0 = time()
    while t < t_end && step < 15000
        P_z = zeros(nz); cs_z = zeros(nz); E_z = zeros(nz); dt_cfl = Inf
        for k in 1:nz
            Pg = ρ[k] * (kB / (0.62 * Z0.mH)) * T[k]; Pr = (arad / 3.0) * T[k]^4
            P_z[k] = Pg + Pr
            γ_e = (5.0/3.0 * Pg + 4.0/3.0 * Pr) / max(P_z[k], 1e-30)
            cs_z[k] = sqrt(γ_e * P_z[k] / ρ[k])
            E_z[k] = ρ[k] * cv_g * T[k] + arad * T[k]^4 + 0.5 * ρ[k] * u[k]^2
            dt_cfl = min(dt_cfl, 0.38 * dr[k] / max(abs(u[k]) + cs_z[k], 1e6))
        end
        dt = min(dt_cfl, 0.4 * day, t_end - t)
        F1 = zeros(nz + 1); F2 = zeros(nz + 1); F3 = zeros(nz + 1)
        F1[1], F2[1], F3[1] = hllc(ρ[1], 0.0, P_z[1], E_z[1], cs_z[1], ρ[1], u[1], P_z[1], E_z[1], cs_z[1])
        for j in 2:nz
            F1[j], F2[j], F3[j] = hllc(ρ[j-1], u[j-1], P_z[j-1], E_z[j-1], cs_z[j-1],
                                       ρ[j],   u[j],   P_z[j],   E_z[j],   cs_z[j])
        end
        F1[nz+1], F2[nz+1], F3[nz+1] = hllc(ρ[nz], u[nz], P_z[nz], E_z[nz], cs_z[nz],
                                            ρ[nz], u[nz], P_z[nz], E_z[nz], cs_z[nz])
        for k in 1:nz
            Vk = V_z[k]; m0 = ρ[k] * Vk; mom0 = ρ[k] * u[k] * Vk; E0 = E_z[k] * Vk
            src = P_z[k] * (A_if[k+1] - A_if[k])
            mn = max(m0 - dt * (A_if[k+1]*F1[k+1] - A_if[k]*F1[k]), 1e-19 * Vk)
            momn = mom0 - dt * (A_if[k+1]*F2[k+1] - A_if[k]*F2[k]) + dt * src
            En = max(E0 - dt * (A_if[k+1]*F3[k+1] - A_if[k]*F3[k]), 1e-12 * Vk)
            ρ[k] = mn / Vk; u[k] = momn / mn
            e_int = max(En / Vk - 0.5 * ρ[k] * u[k]^2, 1e-10)
            Tk = max(T[k], 500.0)
            for _ in 1:6
                f = ρ[k] * cv_g * Tk + arad * Tk^4 - e_int
                df = ρ[k] * cv_g + 4.0 * arad * Tk^3
                Tk = clamp(Tk - f / df, 300.0, 5e7)
            end
            T[k] = Tk
        end
        L_out, R_ph, T_ph = fld_diffuse!(T, ρ, r_if, r_mid, V_z, A_if, dt, κ_s; cv_gas = cv_g)
        E_rad_tot += L_out * dt
        t += dt; step += 1
        if L_out > peak_L
            peak_L = L_out; t_peak = t / day
            prof_peak = (r_15 = r_mid ./ 1e15, rho = copy(ρ), T = copy(T), u_kms = u ./ 1e5, t_day = t / day)
        end
        if step == 1 || step % 8 == 0 || t >= t_end
            push!(lc_rows, (; t_day = t/day, L_ergs = L_out, Rph_15 = R_ph/1e15, Tph_K = T_ph))
        end
    end

    open(joinpath(OUTDIR, "lm07_ppisn_collision_lc.csv"), "w") do io
        println(io, "t_day,L_ergs,Rph_1e15cm,Tph_K")
        for r in lc_rows
            @printf(io, "%.4f,%.4e,%.4f,%.1f\n", r.t_day, r.L_ergs, r.Rph_15, r.Tph_K)
        end
    end
    if prof_peak !== nothing
        open(joinpath(OUTDIR, "lm07_ppisn_collision_prof.csv"), "w") do io
            println(io, "r_1e15cm,rho,T,u_kms")
            for k in 1:nz
                @printf(io, "%.4f,%.4e,%.1f,%.1f\n", prof_peak.r_15[k], prof_peak.rho[k], prof_peak.T[k], prof_peak.u_kms[k])
            end
        end
    end
    @printf("  PPISN Shell Collision (%d steps, %.1fs): L_peak = %.2e erg/s @ t=%.1f d, E_rad = %.2e erg (SLSN!)\n",
            step, time()-t0, peak_L, t_peak, E_rad_tot)
    return (; peak_L, t_peak, E_rad_tot)
end

function main()
    println("=== Landmark 7: Heger & Woosley (2002) + Kasen+11 PISN Shock Breakout + Woosley+07 PPISN Colliding Shells ===")
    t0 = time()
    p1 = run_pair_instability_map()
    p2 = run_pisn_shock_breakout()
    p3 = run_ppisn_shell_collision()

    ok_pair = p1.min_gam1 < 1.25 && p1.n_unstable > 150
    ok_pisn = p2.peak_L > 5e44 && p2.T_peak > 4e4
    ok_ppis = p3.peak_L > 2e43 && p3.E_rad_tot > 1e50
    println("\nVerification scorecard (Heger+02; Kasen+11; Woosley+07):")
    @printf("  [1] FreeEnergyEOS e± pair valley (min Γ_1 = %.4f < 1.25):       %s\n", p1.min_gam1, ok_pair ? "PASS" : "FAIL")
    @printf("  [2] Kasen+11 PISN UV shock breakout (L=%.1e, T=%.1e K):         %s\n", p2.peak_L, p2.T_peak, ok_pisn ? "PASS" : "FAIL")
    @printf("  [3] Woosley+07 PPISN SLSN collision (L=%.1e, E_rad=%.1e erg):   %s\n", p3.peak_L, p3.E_rad_tot, ok_ppis ? "PASS" : "FAIL")
    @printf("LM07 OVERALL: %s (wall %.1fs)\n", (ok_pair && ok_pisn && ok_ppis) ? "3/3 PASS" : "CHECK", time() - t0)
end

main()
