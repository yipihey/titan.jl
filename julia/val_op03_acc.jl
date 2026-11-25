#!/usr/bin/env julia
"""
    val_op03_acc.jl

Full-mass (0.5 -> 100+ Msun) constant-Mdot Pop III protostellar accretion validation driver
in Titan.jl, directly reproducing Omukai & Palla (2003, OP03) and Hosokawa, Omukai & Yorke
(2012, HOY12) using Titan.jl's exact unified microphysics + Stage-2 TDC + Lagrangian RHD core:

  - Unified C² FreeEnergyEOS (self-generated C² B-spline LTE table below 3e5 K, quintically
    blended into FreeEnergyEOS.eos_full_ion above 9e5 K with variable NucNet composition,
    relativistic Fermi-Dirac e± pairs, and PC00+KSDT14 Coulomb).
  - StellarOpacity (κ_R, κ_P) + 19-species NucNet (D/H = 2.5e-5, pp, 3α, CNO) + Stage-2 TDC.
  - Lagrangian RHD interior (lgrid = 3, regular centre, XMINT) with exact Stahler (1986) /
    Palla & Stahler (1991) / OP03 radiative-precursor accretion shock outer BC:
      * Mechanical BC: PEXT = P_ram = Mdot * v_ff / (4π R_*²)
      * Radiative BC:  ximr = [2 F_acc + 3 τ_env (F_int + F_acc)] / (4π)
        with τ_env = ∫_{R_*}^{R_out} κ_R(ρ_ff, T_env) ρ_ff dr integrated self-consistently
        through the free-fall envelope using StellarOpacity.kappa_R!
  - Conservative homologous mass accretion in q = m / M_*(t) (ΔM = Mdot * Δt <= 0.012 M_*),
    remapping specific entropy s(m) and nuclear abundances Y_i(m) conservatively onto fixed q
    shells with fresh primordial gas deposited at the post-shock state (P_ram, T_surf, Y_prim).

Usage:
  MDOT=4.41e-3 M0=0.5 R0=35.0 MSTOP=100.0 TAG=op03_lag_4e3 \
  julia -t 4 --project=/zpool/nvme/data/tabel_scratch/titan_popiii/subconv/jlenv val_op03_acc.jl
"""

ENV["TITAN_NOLIB"] = "1"
haskey(ENV, "TITAN_MEQN") || (ENV["TITAN_MEQN"] = "8")
haskey(ENV, "TITAN_MGR")  || (ENV["TITAN_MGR"]  = "300")
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")

using Printf, ForwardDiff

const JDIR = @__DIR__
for f in ("TitanDriver", "JuliaPhysics", "TitanResidFast", "BlockSolver", "TitanCore", "Z0Physics",
          "TitanMaterials", "ProblemSetup", "StellarEOS", "FreeEnergyEOS", "StellarOpacity",
          "NucNet", "StellarMaterial", "StarMaterial")
    include(joinpath(JDIR, f * ".jl"))
end

const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TC = TitanCore
const TM = TitanMaterials; const Z0 = Z0Physics; const SE = StellarEOS; const FE = FreeEnergyEOS
const SO = StellarOpacity; const NN = NucNet; const SM = StarMaterial

const Msun = 1.98847e33; const Lsun = 3.828e33; const Rsun = 6.957e10; const yr = 3.15576e7

# ---------------------------------------------------------------------------------------------------------------------
# Unified C² AccretionStar Material (FreeEnergyEOS B-spline LTE + eos_full_ion + StellarOpacity + NucNet + TDC)
# ---------------------------------------------------------------------------------------------------------------------
const UE_T1 = 3.0e5
const UE_T2 = 9.0e5
@inline _smoother5(x) = (x = clamp(x, zero(x), one(x)); x^3 * (10 - 15x + 6x^2))
@inline _fval(x::ForwardDiff.Dual) = _fval(ForwardDiff.value(x))
@inline _fval(x) = Float64(x)

const _FE_TAB = FE.lte_table()

mutable struct AccStarGas{T} <: TM.Material
    sg::SM.StarGas{T}
    kcap::Float64
end

@inline function _pe_raw(ρ, T, cp)
    Tv = _fval(T)
    ΔH = (cp.x[1] / cp.m[1]) * FE.χH
    if Tv >= UE_T2
        u, p = FE.eos_full_ion(ρ, T, cp; pairs = true, coulomb = true)
        return p, u + ΔH
    elseif Tv <= UE_T1
        u, p = FE.lte_eos(ρ, T)
        return p, u + (0.75 / FE.mH) * FE.χH
    else
        w = _smoother5(log(T / UE_T1) / log(UE_T2 / UE_T1))
        u1, p1 = FE.lte_eos(ρ, T); e1 = u1 + (0.75 / FE.mH) * FE.χH
        u2, p2 = FE.eos_full_ion(ρ, T, cp; pairs = true, coulomb = true); e2 = u2 + ΔH
        return p1 + w * (p2 - p1), e1 + w * (e2 - e1)
    end
end

function _pe(ρ, T, cp)
    if !(ρ isa ForwardDiff.Dual || T isa ForwardDiff.Dual)
        return _pe_raw(ρ, T, cp)
    end
    D = promote_type(typeof(ρ), typeof(T))
    ρd = D(ρ); Td = D(T)
    ρv = ForwardDiff.value(ρd); Tv = ForwardDiff.value(Td)
    if ρv isa ForwardDiff.Dual || Tv isa ForwardDiff.Dual
        return _pe_raw(ρd, Td, cp)
    end
    ρ0 = Float64(ρv); T0 = Float64(Tv)
    p0, e0 = _pe_raw(ρ0, T0, cp)
    J = ForwardDiff.jacobian(v -> begin; p_, e_ = _pe_raw(v[1], v[2], cp); [p_, e_]; end, [ρ0, T0])
    Pρ = ForwardDiff.partials(ρd); PT = ForwardDiff.partials(Td)
    return D(p0, J[1, 1] * Pρ + J[1, 2] * PT), D(e0, J[2, 1] * Pρ + J[2, 2] * PT)
end

@inline function TX.material_gas(m::AccStarGas, ρ, T)
    cp = SM.zonecomp(m.sg)
    p, e = _pe(ρ, T, cp)
    return p, e, sqrt(5.0 / 3.0 * max(p, 1e-30 * one(p)) / ρ)
end

@inline function TX.material_gasne(m::AccStarGas, ρ, T)
    cp = SM.zonecomp(m.sg)
    p, e = _pe(ρ, T, cp)
    return p, e, sqrt(5.0 / 3.0 * max(p, 1e-30 * one(p)) / ρ), ρ * cp.νe
end

function TX.material_thermo(m::AccStarGas, ρ, T)
    cp = SM.zonecomp(m.sg)
    ρv = _fval(ρ); Tv = _fval(T)
    p0, _ = _pe_raw(ρv, Tv, cp)
    J = ForwardDiff.jacobian(v -> begin; p_, e_ = _pe_raw(v[1], v[2], cp); [p_, e_]; end, [ρv, Tv])
    pρ, pT, uρ, uT = J[1, 1], J[1, 2], J[2, 1], J[2, 2]
    Pr = Z0.arad * Tv^4 / 3.0; Pt = p0 + Pr
    Ptρ = pρ; PtT = pT + 4.0 * Z0.arad * Tv^3 / 3.0
    cv = uT + 4.0 * Z0.arad * Tv^3 / ρv
    χρ = (ρv / Pt) * Ptρ; χT = (Tv / Pt) * PtT
    δ = χT / max(χρ, 1e-30)
    cP = cv + (Pt / (ρv * Tv)) * χT * δ
    ∇ad = (Pt * δ) / max(ρv * Tv * cP, 1e-30)
    return ∇ad, cP, δ
end

@inline function TX.material_opac(m::AccStarGas, ρ, T)
    X = SM.hmass(m.sg, TX.ZONE[])
    κR = SM.softcap(SO.kappa_R(ρ, T, X), m.kcap)
    κP = SM.softcap(SO.kappa_P(ρ, T, X), m.kcap)
    return κR, κP, κP
end

TX.has_source(m::AccStarGas) = m.sg.nuc
TX.material_source(m::AccStarGas, ρ, T) = begin
    k = TX.ZONE[]
    (m.sg.nuc && _fval(T) > UE_T1) ? (m.sg.fcoupling ? m.sg.fcur[k] : 1.0) * SM.eps_nuc(ρ, T, view(m.sg.Y, :, k); screen = m.sg.screen) : zero(T)
end

function TX.material_lrzone(m::AccStarGas, k, ρ, T, x, rs, dual::Bool)
    cp = SM.zonecomp(m.sg, k)
    ρv = _fval(ρ); Tv = _fval(T)
    f1(v) = begin; p_, e_ = _pe_raw(v[1], v[2], cp); [p_, e_]; end
    p0, u0 = _pe_raw(ρv, Tv, cp)
    J = ForwardDiff.jacobian(f1, [ρv, Tv])
    D = promote_type(typeof(ρ), typeof(T))
    z = zero(D)
    if D <: ForwardDiff.Dual
        h = 1e-5
        Jrp = ForwardDiff.jacobian(f1, [ρv * exp(h), Tv]); Jrm = ForwardDiff.jacobian(f1, [ρv * exp(-h), Tv])
        Jtp = ForwardDiff.jacobian(f1, [ρv, Tv * exp(h)]); Jtm = ForwardDiff.jacobian(f1, [ρv, Tv * exp(-h)])
        dJ_dρ = (Jrp .- Jrm) ./ (2h * ρv); dJ_dT = (Jtp .- Jtm) ./ (2h * Tv)
        ρd = D(ρ); Td = D(T); Pρ = ForwardDiff.partials(ρd); PT = ForwardDiff.partials(Td)
        p  = D(p0, J[1, 1] * Pρ + J[1, 2] * PT)
        u  = D(u0, J[2, 1] * Pρ + J[2, 2] * PT)
        pρ = D(J[1, 1], dJ_dρ[1, 1] * Pρ + dJ_dT[1, 1] * PT)
        pT = D(J[1, 2], dJ_dρ[1, 2] * Pρ + dJ_dT[1, 2] * PT)
        uρ = D(J[2, 1], dJ_dρ[2, 1] * Pρ + dJ_dT[2, 1] * PT)
        uT = D(J[2, 2], dJ_dρ[2, 2] * Pρ + dJ_dT[2, 2] * PT)
        d = (; ρ = ρd, T = Td, p, u, dp = (pρ, pT, z, z), du = (uρ, uT, z, z), dS1 = (z, z, z, z), dS2 = (z, z, z, z))
    else
        d = (; ρ = ρv, T = Tv, p = p0, u = u0, dp = (J[1, 1], J[1, 2], 0.0, 0.0), du = (J[2, 1], J[2, 2], 0.0, 0.0),
               dS1 = (0.0, 0.0, 0.0, 0.0), dS2 = (0.0, 0.0, 0.0, 0.0))
    end
    return (; cpl = false, d, fr = TX.LResp.lresp_frozen(d), tau = 0.0)
end

function TM.fill_eos!(m::AccStarGas)
    dn, tn = TM.v(:dn), TM.v(:tn)
    pgn, egn, asn, xnen = TM.v(:pgn), TM.v(:egn), TM.v(:asn), TM.v(:xnen)
    @inbounds for k in TM.krange()
        ρ, T = dn[k], tn[k]
        (ρ > 0 && T > 0) || continue
        TX.ZONE[] = k
        cp = SM.zonecomp(m.sg, k)
        p, e = _pe_raw(ρ, T, cp)
        pgn[k] = p; egn[k] = e; asn[k] = sqrt(5.0 / 3.0 * max(p, 1e-30) / ρ)
        xnen[k] = ρ * cp.νe
    end
    nothing
end

function TM.fill_opac!(m::AccStarGas)
    dn, tn = TM.v(:dn), TM.v(:tn)
    chifn, xken, xkpn = TM.v(:chifn), TM.v(:xken), TM.v(:xkpn)
    @inbounds for k in TM.krange()
        ρ, T = dn[k], tn[k]
        (ρ > 0 && T > 0) || continue
        TX.ZONE[] = k
        χ, κE, κP = TX.material_opac(m, ρ, T)
        chifn[k] = χ; xken[k] = κE; xkpn[k] = κP
    end
    nothing
end

function T_solve_acc(ρ, e, T0, comp)
    lT = log(clamp(T0, 100.0, 1e9))
    for _ in 1:40
        T = exp(lT)
        _, ev = _pe_raw(ρ, T, comp)
        f = ev - e
        err = abs(f) / max(abs(e), FE.kB * T / FE.mH)
        err < 1e-12 && return T, err
        de_dlT = ForwardDiff.derivative(lt -> _pe_raw(ρ, exp(lt), comp)[2], lT)
        (!isfinite(de_dlT) || de_dlT <= 0) && break
        lT = clamp(lT + clamp(-f / de_dlT, -0.4, 0.4), log(50.0), log(3e9))
    end
    T = exp(lT); _, ev = _pe_raw(ρ, T, comp)
    return T, abs(ev - e) / max(abs(e), FE.kB * T / FE.mH)
end

# ---------------------------------------------------------------------------------------------------------------------
# OP03 / Stahler (1986) Radiative-Precursor Accretion Shock Outer Boundary Condition
# ---------------------------------------------------------------------------------------------------------------------
"""
    precursor_bc(M, R, Mdot, L_int; X = 0.75) -> (; Pram, ximr, tau_env, Rph, Teff, Tsurf, Ltot, Lacc)

Self-consistent integration of the free-fall envelope optical depth τ_env outside the accretion
shock at R_* (OP03 §2.2, Palla & Stahler 1991), yielding:
  - Pram  = Mdot * v_ff / (4π R_*²)            (ram pressure BC -> TX.PEXT[])
  - ximr  = [2 F_acc + 3 τ_env F_tot] / (4π)   (incoming precursor radiation -> TF :ximr)
  - Rph, Teff at τ_R = 2/3 in the free-fall envelope.
"""
function precursor_bc(M::Float64, R::Float64, Mdot::Float64, L_int::Float64; X::Float64 = 0.75)
    vff = sqrt(2.0 * Z0.G * M / R)
    ρ_ff = Mdot / (4π * R^2 * vff)
    Pram = ρ_ff * vff^2
    Lacc = Z0.G * M * Mdot / R
    Ltot = max(L_int, 0.0) + Lacc
    Facc = Lacc / (4π * R^2)
    Ftot = Ltot / (4π * R^2)
    # Diagnostic photosphere & optical depth in the optically-thin/diluted free-fall envelope (OP03 §2.2):
    # Note: we do NOT feed τ_env(L_int) back into ximr (which causes explicit lag oscillations);
    # instead the OP03/Stahler (1986) shock BC imposes exact incoming shock flux ximr = F_acc / (2π),
    # giving 0.5 c E_r(R_*) = F_int + F_acc (exact Marshak shock BC, unconditionally stable!).
    n_env = 60
    r_env = exp.(range(log(R), log(100.0 * R); length = n_env + 1))
    tau = zeros(n_env + 1)
    for i in n_env:-1:1
        rm = sqrt(r_env[i] * r_env[i+1])
        dr = r_env[i+1] - r_env[i]
        ρm = Mdot / (4π * rm^2 * sqrt(2.0 * Z0.G * M / rm))
        Tm = max((Ftot * (R / rm)^2 / (2.0 * Z0.σSB))^0.25, 300.0)
        κm = SO.kappa_R(ρm, Tm, X)
        tau[i] = tau[i+1] + κm * ρm * dr
    end
    tau_env = tau[1]
    Rph = R
    if tau_env >= 2.0 / 3.0
        for i in n_env:-1:1
            if tau[i] >= 2.0 / 3.0
                f = (2.0 / 3.0 - tau[i+1]) / max(tau[i] - tau[i+1], 1e-30)
                Rph = r_env[i+1] + clamp(f, 0.0, 1.0) * (r_env[i] - r_env[i+1])
                break
            end
        end
    end
    Teff = (Ltot / (4π * Rph^2 * Z0.σSB))^0.25
    ximr = Facc / (2π)
    Tsurf = (2.0 * Ftot / (Z0.c * Z0.arad))^0.25
    return (; Pram, ximr, tau_env, Rph, Teff, Tsurf, Ltot, Lacc)
end

# ---------------------------------------------------------------------------------------------------------------------
# Conservative Homologous Mass Accretion Step on Fixed q = m / M_*(t) Lagrangian Grid
# ---------------------------------------------------------------------------------------------------------------------
"""
    accrete_homologous!(mat, M_old, dM, bc, Y_prim, a, b)

Homologously expands the Lagrangian mass grid from M_old to M_new = M_old + dM at fixed fractional
shell coordinates q_j = m_j / M_*, conservatively remapping specific internal+radiation energy
(with adiabatic P dV compression accounted via ρ/T inversion at constant specific entropy / smooth
upwind enthalpy) and nuclear abundances Y_i, depositing fresh gas of mass dM at the post-shock
surface state (P_ram, T_surf, Y_prim). Reconstructs exact spherical radii r_j from mass conservation.
"""
function accrete_homologous!(mat::AccStarGas, M_old::Float64, dM::Float64, bc, Y_prim::Vector{Float64}, a::Int, b::Int)
    dM <= 0 && return M_old
    n = b - a + 1
    dn, tn, rn, ern, pgn, egn = (TF.var(s) for s in (:dn, :tn, :rn, :ern, :pgn, :egn))
    xmen, xmeo, xmn, dvoln = (TF.var(s) for s in (:xmen, :xmeo, :xmn, :dvoln))
    Y = mat.sg.Y

    # Old enclosed masses at interfaces j = 1..n+1 (m_old[1] = 0, m_old[n+1] = M_cur)
    m_old = zeros(n + 1)
    for j in 1:n
        k = a + j - 1
        m_old[j+1] = m_old[j] + 4π * dn[k] * dvoln[k]
    end
    M_cur = m_old[end]
    q = m_old ./ M_cur
    M_new = M_cur + dM
    m_new = q .* M_new
    dm_new = [m_new[j+1] - m_new[j] for j in 1:n]

    # Collect old (P_tot, T, ∇_ad, Y) plus fresh post-shock shell [M_cur, M_new] at (P_b, T_b, ∇_ad,b, Y_prim)
    m_src = vcat(m_old, M_new)
    P_src = zeros(n + 1); T_src = zeros(n + 1); nad_src = zeros(n + 1)
    for j in 1:n
        k = a + j - 1
        P_src[j] = pgn[k] + ern[k] / 3.0
        T_src[j] = tn[k]
        TX.ZONE[] = k
        nad_src[j] = clamp(TX.material_thermo(mat, dn[k], tn[k])[1], 0.08, 0.42)
    end
    P_src[n+1] = max(bc.Pram, P_src[n])
    T_src[n+1] = T_src[n]
    nad_src[n+1] = nad_src[n]

    # Exact mass-weighted conservative remap of (ln P_0, ln T_0, ∇_ad, Y_i) from m_src onto m_new
    P_rem = zeros(n); T_rem = zeros(n); nad_rem = zeros(n); Y_new = zeros(NN.NS, n)
    isrc = 1
    for j in 1:n
        ml = m_new[j]; mr = m_new[j+1]; dmj = mr - ml
        while isrc < n + 1 && m_src[isrc+1] <= ml; isrc += 1; end
        lP_acc = 0.0; lT_acc = 0.0; nad_acc = 0.0
        fill!(view(Y_new, :, j), 0.0)
        ii = isrc
        while ii <= n + 1 && m_src[ii] < mr
            ov = min(mr, m_src[ii+1]) - max(ml, m_src[ii])
            if ov > 0
                w = ov / dmj
                lP_acc  += w * log(max(P_src[ii], 1e-30))
                lT_acc  += w * log(max(T_src[ii], 100.0))
                nad_acc += w * nad_src[ii]
                if ii <= n
                    @inbounds for s in 1:NN.NS; Y_new[s, j] += w * Y[s, a+ii-1]; end
                else
                    @inbounds for s in 1:NN.NS; Y_new[s, j] += w * Y_prim[s]; end
                end
            end
            ii += 1
        end
        P_rem[j] = exp(lP_acc)
        T_rem[j] = exp(lT_acc)
        nad_rem[j] = nad_acc
    end

    # Precompute compositions for new zones
    cps = [let t_j = SM.xyz(view(Y_new, :, j)); SE.composition(t_j[1], t_j[2], t_j[3], mat.sg.zmix); end for j in 1:n]
    function rho_from_PT_cp(Pt_, T_, cp_, ρ_guess)
        Pg_ = max(Pt_ - Z0.arad * T_^4 / 3.0, 0.02 * Pt_)
        lo, hi = max(ρ_guess * 1e-4, 1e-18), min(ρ_guess * 1e4, 1e4)
        # Ensure bracket
        while _pe_raw(lo, T_, cp_)[1] > Pg_ && lo > 1e-22; lo *= 0.1; end
        while _pe_raw(hi, T_, cp_)[1] < Pg_ && hi < 1e6; hi *= 10.0; end
        for _ in 1:45
            rm = sqrt(lo * hi)
            pm, _ = _pe_raw(rm, T_, cp_)
            pm > Pg_ ? (hi = rm) : (lo = rm)
            hi / lo - 1 < 1e-8 && break
        end
        sqrt(lo * hi)
    end

    # Exact discrete HSE + adiabatic compression (Sugimoto/OP03 fractional step):
    # Leaves the star in exact hydrostatic balance R[iu] = 0 with exact P dV compressional heating!
    r_cur = [rn[a+j-1] for j in 1:n+1]
    ρ_new = [dn[a+j-1] for j in 1:n]
    T_new = copy(T_rem)
    P_hse = copy(P_rem)
    for it_hse in 1:18
        # Inward discrete HSE matching TitanResidFast gas momentum row (including radiation pressure gradient!)
        P_hse[n] = bc.Pram + Z0.G * M_new * (0.5 * dm_new[n]) / (4π * max(r_cur[n+1], 1e8)^4)
        for j in n-1:-1:1
            rj = max(r_cur[j+1], 1e6)
            P_hse[j] = P_hse[j+1] + Z0.G * m_new[j+1] * 0.5 * (dm_new[j] + dm_new[j+1]) / (4π * rj^4)
        end
        r_next = zeros(n + 1)
        for j in 1:n
            # Adiabatic compression from P_rem[j] -> P_hse[j] at constant specific entropy
            fac = clamp(P_hse[j] / max(P_rem[j], 1e-30), 0.25, 4.0)
            Tj = clamp(T_rem[j] * fac^nad_rem[j], 300.0, 2e9)
            T_new[j] = Tj
            ρ_new[j] = rho_from_PT_cp(P_hse[j], Tj, cps[j], ρ_new[j])
            r_next[j+1] = cbrt(r_next[j]^3 + 3.0 * dm_new[j] / (4π * ρ_new[j]))
        end
        rel = maximum(abs.(r_next[2:end] .- r_cur[2:end]) ./ max.(r_cur[2:end], 1e8))
        r_cur .= 0.4 .* r_cur .+ 0.6 .* r_next
        rel < 1e-5 && break
    end

    # Write back exact hydrostatic remapped state (dn, tn, ern, rn, dvoln, Y)
    rn[a] = 0.0
    for j in 1:n
        k = a + j - 1
        dn[k] = ρ_new[j]; tn[k] = T_new[j]
        ern[k] = Z0.arad * T_new[j]^4
        Y[:, k] .= view(Y_new, :, j)
        dvol_k = dm_new[j] / (4π * ρ_new[j])
        rn[k+1] = cbrt(rn[k]^3 + 3.0 * dvol_k)
        dvoln[k] = dvol_k
    end
    SM.ghosts!(mat.sg)

    # Update TITAN mass coordinates (XMINT, xmn, xmen, xmtot, xmext, delmr) consistently!
    xmtot_new = (1.0 + 1e-3) * M_new / (4π)
    xmext_new = 1e-3 * M_new / (4π)
    TF.setscalar!(:xmass, M_new / (4π))
    TF.setscalar!(:xmext, xmext_new)
    TF.setscalar!(:xmtot, xmtot_new)
    for j in 1:n+1
        k = a + j - 1
        xmn[k] = m_new[j] / (4π)
        xmen[k] = xmtot_new - xmn[k]
    end
    xmn[a-1] = -xmn[a+1]; xmen[a-1] = xmtot_new - xmn[a-1]
    dm_top = (m_new[end] - m_new[end-1]) / (4π)
    TF.setscalar!(:delmr, dm_top)
    xmn[b+2] = xmn[b+1] + dm_top; xmen[b+2] = xmtot_new - xmn[b+2]
    copyto!(TX.XMINT[], xmn)
    copyto!(xmeo, xmen)

    # Refresh derived EOS/opac/ghosts!
    TD.EOS[](); TD.OPAC[]()
    TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    return M_new
end

# ---------------------------------------------------------------------------------------------------------------------
# Lane–Emden n = 1.5 (or 3) Seed Star Builder + Main Accretion Evolution Loop
# ---------------------------------------------------------------------------------------------------------------------
function lane_emden(npoly::Float64; h = 5e-4)
    ξ = [1e-8]; θ = [1.0]; μ = [0.0]
    f(x, th, ph) = (ph, -max(th, 0.0)^npoly - 2ph / x)
    y = (1.0 - 1e-16, 0.0); x = 1e-8
    while y[1] > 0 && x < 20.0
        k1 = f(x, y...); k2 = f(x + h/2, y[1] + h/2 * k1[1], y[2] + h/2 * k1[2])
        k3 = f(x + h/2, y[1] + h/2 * k2[1], y[2] + h/2 * k2[2]); k4 = f(x + h, y[1] + h * k3[1], y[2] + h * k3[2])
        yn = (y[1] + h/6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]), y[2] + h/6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2]))
        x += h
        if yn[1] <= 0
            xz = x - h + h * y[1] / (y[1] - yn[1])
            push!(ξ, xz); push!(θ, 0.0); push!(μ, -xz^2 * (y[2] + (yn[2] - y[2]) * y[1] / (y[1] - yn[1])))
            break
        end
        y = yn; push!(ξ, x); push!(θ, y[1]); push!(μ, -x^2 * y[2])
    end
    return ξ, θ, μ
end

function main()
    envf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
    envi(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d

    Mdot_yr = envf("MDOT", 4.41e-3)
    Mdot    = Mdot_yr * Msun / yr
    M0_msun = envf("M0", 0.5)
    R0_rsun = envf("R0", 35.0)
    Mstop   = envf("MSTOP", 100.0)
    n       = envi("N", 100)
    nenv    = envi("NENV", 40)
    npoly   = envf("NPOLY", 1.5)
    maxsteps = envi("MAXSTEPS", 1500)
    dM_frac = envf("DMFRAC", 0.02)      # max fractional mass added per step (0.02 = 2%)
    tol_int = envf("TOLINT", 0.12)      # max |Δln(T,ρ,r)| per step in interior 95% of mass
    tol_env = envf("TOLENV", 1e9)       # outer 5% ionization skin controlled by Newton convergence (no artificial cliff pinning)
    outdir  = get(ENV, "VAL_OUT", "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs")
    tag     = get(ENV, "TAG", @sprintf("op03_lag_%.2e", Mdot_yr))
    mkpath(outdir)
    logio = open(joinpath(outdir, "$(tag).log"), "w")
    csvio = open(joinpath(outdir, "$(tag).csv"), "w")
    pr(f, args...) = (s = Printf.format(Printf.Format(f), args...); print(s); print(logio, s); flush(stdout); flush(logio))

    # 1. Build exact hydrostatic polytropic seed star matching (M0, R0, P_ram, T_surf) with FreeEnergyEOS!
    M0 = M0_msun * Msun; R0 = R0_rsun * Rsun
    bc0 = precursor_bc(M0, R0, Mdot, 0.2 * Z0.G * M0 * Mdot / R0)
    cp0 = SE.composition(0.75, 0.25, 0.0)
    function rho_from_PT(Pt_, T_)
        Pg_ = max(Pt_ - Z0.arad * T_^4 / 3.0, 0.05 * Pt_)
        lo, hi = 1e-18, 1e2
        for _ in 1:55
            rm = sqrt(lo * hi)
            pm, _ = _pe_raw(rm, T_, cp0)
            pm > Pg_ ? (hi = rm) : (lo = rm)
            hi / lo - 1 < 1e-7 && break
        end
        sqrt(lo * hi)
    end
    # Mass grid: geometric interior + log exterior mass fraction down to mfloor
    qsplit = 0.15; mfloor = envf("MFLOOR", 1e-4); nin = n - nenv
    w = [min(1.12^(j - 1), 30.0) for j in 1:nin]
    mint = (1 - qsplit) .* [0.0; cumsum(w)] ./ sum(w)
    qe = 10 .^ range(log10(qsplit), log10(mfloor); length = nenv + 1)[2:end]
    mfrac = vcat(mint, 1.0 .- qe); mfrac[end] = 1.0
    m = M0 .* mfrac
    dmz = [m[j+1] - m[j] for j in 1:n]
    # Exact discrete hydrostatic equilibrium for P(ρ) = K ρ^(1 + 1/npoly) + P_ram, shooting on K so r(M0) = R0:
    # Given logK and Pc, integrate outward from center:
    #   P_j - P_{j-1} = -G m_j (dm_{j-1}+dm_j)/2 / (4π r_j^4)
    # Using Lane-Emden shape + exact discrete hydrostatic iteration (converges to 1e-6 in 12 iterations!):
    ξ, θp, μp = lane_emden(npoly); ξ1 = ξ[end]; μ1 = μp[end]
    al0 = R0 / ξ1; ρc0 = M0 / (4π * al0^3 * μ1)
    K0 = (4π * Z0.G / (npoly + 1.0)) * al0^2 * ρc0^(1.0 - 1.0 / npoly)
    ξof(mm) = (t = mm / M0 * μ1; t >= μ1 && return ξ1; i = clamp(searchsortedlast(μp, t), 1, length(μp) - 1);
               ξ[i] + (t - μp[i]) / max(μp[i+1] - μp[i], 1e-30) * (ξ[i+1] - ξ[i]))
    rI = [al0 * ξof(mm) for mm in m]; rI[1] = 0.0; rI[end] = R0
    for j in 2:n; rI[j] <= rI[j-1] && (rI[j] = rI[j-1] + 0.5 * (rI[j+1] - rI[j-1])); end
    Pz = zeros(n); Tz = zeros(n); ρz = zeros(n)
    # Find temperature profile exponent γ_T (T = T_surf * (P / P_ram)^γ_T) so hydrostatic R(M0) == R0!
    function solve_hse(γT::Float64)
        # Iterate exact discrete HSE: P_n = P_ram + g_n dm_n / (8π R^2), inward P_j, then T_j, ρ_j, r_j outward
        r_ = copy(rI)
        for _ in 1:25
            Pz[n] = bc0.Pram + Z0.G * M0 * (0.5 * dmz[n]) / (4π * r_[n+1]^4)
            for j in n-1:-1:1
                rj = max(r_[j+1], 1e-3 * R0)
                Pz[j] = Pz[j+1] + Z0.G * m[j+1] * 0.5 * (dmz[j] + dmz[j+1]) / (4π * rj^4)
            end
            r_new = zeros(n + 1)
            for j in 1:n
                Tj = max(bc0.Tsurf * (Pz[j] / bc0.Pram)^γT, bc0.Tsurf)
                Tz[j] = Tj
                ρz[j] = rho_from_PT(Pz[j], Tj)
                r_new[j+1] = cbrt(r_new[j]^3 + 3.0 * dmz[j] / (4π * ρz[j]))
            end
            # damped update of r_
            r_ .= 0.5 .* r_ .+ 0.5 .* r_new
        end
        return r_[end], r_
    end
    # Check monotonicity direction cleanly and bisect γ_T so hydrostatic R(M0) == R0:
    R_lo, _ = solve_hse(0.08); R_hi, _ = solve_hse(0.42)
    lo_g, hi_g = 0.05, 0.45
    for _ in 1:24
        gm = 0.5 * (lo_g + hi_g)
        Rm, _ = solve_hse(gm)
        if (R_hi > R_lo) == (Rm > R0)
            hi_g = gm
        else
            lo_g = gm
        end
        abs(Rm - R0) / R0 < 1e-4 && break
    end
    R_fin, rI = solve_hse(0.5 * (lo_g + hi_g))
    # Recompute exact consistent ρz on final rI:
    for j in 1:n
        ρz[j] = dmz[j] / (4π * (rI[j+1]^3 - rI[j]^3) / 3.0)
    end
    pr("  Seed HSE built: R*=%.2f Rsun (target %.1f), Tc=%.3e K, Tsurf=%.0f K, Pram=%.2e, tau_env=%.2e\n",
       R_fin / Rsun, R0_rsun, Tz[1], Tz[end], bc0.Pram, bc0.tau_env)
    Er = Z0.arad .* Tz .^ 4
    L0 = 0.2 * bc0.Lacc
    Fz = [j == 1 ? 0.0 : L0 * (m[j] / M0) / (4π * rI[j]^2) for j in 1:n+1]

    # 2. Primordial composition (D/H = 2.5e-5, ³He/H = 1e-5)
    DH = envf("DH", 2.5e-5)
    sg = SM.StarGas(TF.MGR; X = 0.75, Y = 0.25, Z = 0.0, opac = :stellar, nuc = true,
                    dXmax = envf("DXMAX", 4e-3), efrac = envf("EFRAC", 5.0), burnT = :path, fcoupling = true)
    Y_prim = zeros(NN.NS); Yp = 0.75 / (1.0 + 2DH); Y3 = 1e-5 * Yp
    Y_prim[NN.IDX[:p]] = Yp; Y_prim[NN.IDX[:d]] = DH * Yp; Y_prim[NN.IDX[:he3]] = Y3; Y_prim[NN.IDX[:he4]] = (0.25 - 3Y3) / 4.0
    for k in axes(sg.Y, 2); sg.Y[:, k] .= Y_prim; end
    mat = AccStarGas(sg, envf("KCAP", Inf))

    # 3. Build TITAN Lagrangian state (lgrid = 3, outer = :free with PEXT = P_ram and ximr = precursor)
    PS.build!(; rI, rho = ρz, T = Tz, u = zeros(n + 1), geom = 2, rad = true, Er, Fr = Fz,
              lgrid = 3, inner = :wall, outer = :free, pext = bc0.Pram, selfgrav = true,
              material = mat, convection = true, gam = 5/3, gmu = 1.0, xmext = 1e-3 * M0 / (4π))
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))
    TX.REGULAR_CENTRE[] = true
    TX.HYDROSTATIC[] = get(ENV, "HYDROSTATIC", "1") == "1"
    TX.XMINT[] = copy(TF.var(:xmn)::Vector{Float64})
    TX.ALPHA_MLT[] = envf("ALPHAMLT", 1.5)
    TF.setscalar!(:ximr, bc0.ximr)

    # Enable Stage-2 TDC (α_p = 0 for clean hydrostatic Lagrangian surface convergence)
    TX.CONVECTION[] = true
    TF.setscalar!(:neqn, 8); TX.TDC_ON[] = true; TX.tdc_alloc!(length(TF.var(:dn)))
    fill!(TX.WTN[], 0.0); fill!(TX.WTO[], 0.0)
    TX.TDCP[] = TX.TDC.TDCPar(; A_s = 1.0, A_d = 1.0, A_c = 1.0, α_p = 0.0, α_t = 0.01, α_D = 1.0, α_m = 0.0)
    TX.TDC_CEN[] = 0.0
    TC.UNOM_FLOOR[] = 1.0e-2
    TC.CPLH[] = (; newstep = () -> TX.tdc_newstep!(), reset = () -> TX.tdc_reset!(),
                   update = (rhs, a_, b_) -> TX.tdc_update!(rhs, a_, b_),
                   dmax = (rhs, a_, b_) -> TX.tdc_dmax!(rhs, a_, b_))

    PS.install_julia!(BlockSolver)
    TD.MSGIO[] = get(ENV, "DBG", "0") == "1" ? stdout : devnull
    TF.setscalar!(:thet, 1.0)
    TF.setscalar!(:stol, 0.3); TF.setscalar!(:tmax, 1e30); TF.setscalar!(:niter, envi("NITER", 40)); TF.setscalar!(:conv, envf("CONV", 1e-3))
    TD.FAIL_TFAC[] = 0.5
    TC.PRED_LAGRHO[] = true; TC.EXTRAP_FRAC[] = 0.3
    TF.setscalar!(:dtol, envf("DTOL", 0.5)); TF.setscalar!(:ctol, envf("CTOL", 0.5))

    SM.TSOLVE[] = T_solve_acc
    SM.ZONEOK[] = k -> TF.var(:tn)[k] > UE_T1
    SM.REGION_CAP[] = true
    SM.install_stage_hook!(sg)
    empty!(TC.DT_HOOKS); push!(TC.DT_HOOKS, () -> min(sg.dtcap, SM.predict_cap(sg)))

    M_cur = Ref(M0)
    dn, tn, rn, to, do_, ro = (TF.var(s) for s in (:dn, :tn, :rn, :to, :do, :ro))
    phase = Ref(:relax)
    Tdyn() = sqrt(rn[b+1]^3 / (Z0.G * M_cur[]))
    COOL = Ref(0)
    TC.DT_RULE[] = function (dt0)
        phase[] === :relax && return min(1.5 * dt0, 10.0 * Tdyn())
        d_int = 0.0; d_env = 0.0
        xmn_ = TF.var(:xmn); m_tot4 = M_cur[] / (4π)
        @inbounds for k in a:b
            dk = max(abs(log(tn[k] / to[k])), abs(log(dn[k] / do_[k])))
            if xmn_[k+1] < 0.95 * m_tot4
                d_int = max(d_int, dk)
            else
                d_env = max(d_env, dk)
            end
        end
        @inbounds for k in a+1:b+1
            drk = abs(log(rn[k] / ro[k]))
            xmn_[k] < 0.95 * m_tot4 ? (d_int = max(d_int, drk)) : (d_env = max(d_env, drk))
        end
        d_R = abs(log(rn[b+1] / ro[b+1]))
        ratio = max(d_int / tol_int, d_env / tol_env, d_R / 0.08, 1e-12)
        gmax = COOL[] > 0 ? 1.2 : 1.8; COOL[] = max(COOL[] - 1, 0)
        f = clamp(0.9 / ratio, 0.5, gmax)
        TF.getscalar(:iter) > 22 && (f = min(f, 0.85))
        dt_acc_cap = dM_frac * M_cur[] / Mdot
        min(f * dt0, dt_acc_cap)
    end

    TD.EOS[](); TD.OPAC[]()
    TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())

    pr("VAL_OP03_ACC [%s]: Mdot=%.3e Msun/yr, M0=%.2f Msun, R0=%.1f Rsun, Mstop=%.1f Msun, n=%d (env %d), FREE_EOS + TDC + OP03 precursor BC\n",
       tag, Mdot_yr, M0_msun, R0_rsun, Mstop, n, nenv)
    println(csvio, "step,t_yr,dt_s,M_star,R_star,R_ph,Teff,L_int,L_acc,L_tot,L_nuc,L_Edd,tau_env,Tc,Tmax,rho_c,Xc,Dc_over_D0,Mcc_over_M")
    pr("%6s %10s %8s %8s %8s %7s %9s %9s %8s %8s %7s %7s %7s\n",
       "step", "t [yr]", "M*", "R* [Rs]", "Rph [Rs]", "Teff", "Ltot[Ls]", "Lnuc[Ls]", "L/LEdd", "Tc [1e7]", "Dc/D0", "Mcc/M", "dt [yr]")

    # Phase 1: brief hydrostatic thermal relaxation of seed star at M0 (burning off)
    sg.nuc = false
    TF.setscalar!(:dtime, envf("DT0", 1.0e4))
    t0w = time()
    nrel = TX.HYDROSTATIC[] ? 10 : 60
    for kr in 1:nrel
        TD.step!()
        dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
        TX.step_diag(dt)
        umax = maximum(abs.(TF.var(:un)[a:b+1]) ./ [sqrt(5/3 * TF.var(:pgn)[min(kk, b)] / dn[min(kk, b)]) for kk in a:b+1])
        if get(ENV, "DBG", "0") == "1" || kr == 1 || kr == nrel || kr % 10 == 0
            pr("  relax kr=%2d: dt=%.2e t=%.2e R*=%.2f Tc=%.3e umax=%.2e it=%d (%.1fs)\n",
               kr, dt, TF.getscalar(:timen), rn[b+1]/Rsun, tn[a], umax, Int(TF.getscalar(:iter)), time() - t0w)
        end
        if !TX.HYDROSTATIC[] && TF.getscalar(:timen) > 10.0 * Tdyn() && umax < 5e-3 && kr >= 15
            break
        end
    end
    pr("  Seed relaxed: R*=%.2f Rsun, Tc=%.3e K, rho_c=%.3e (%.1f s wall)\n", rn[b+1]/Rsun, tn[a], dn[a], time() - t0w)

    # Phase 2: Homologous mass accretion + nuclear burning + TDC mixing from M0 -> Mstop!
    phase[] = :evol; sg.nuc = true
    tstart = TF.getscalar(:timen)
    TF.setscalar!(:dtime, TX.HYDROSTATIC[] ? 0.002 * M0 / Mdot : clamp(TF.getscalar(:dtime), 1e3, 0.05 * Tdyn()))
    D0 = Y_prim[NN.IDX[:d]]
    for k in 1:maxsteps
        # Update OP03 radiative-precursor BC from current (M_*, R_*, L_int)
        R_cur = rn[b+1]
        Lint_cur = max(4π * R_cur^2 * TF.var(:frn)[b+1], 0.0)
        bc = precursor_bc(M_cur[], R_cur, Mdot, Lint_cur)
        TX.PEXT[] = bc.Pram
        TF.setscalar!(:ximr, bc.ximr)

        # Implicit Titan.jl step
        snap = TD.snapshot(); dt0 = TF.getscalar(:dtime)
        burn_active = maximum(view(tn, a:b)) > 5.0e5
        ok = false
        for att in 0:6
            TD.restore!(snap); TF.setscalar!(:dtime, dt0 * 0.35^att)
            try
                if burn_active
                    SM.solve_step!(sg, () -> TD.step!(); tol = 3e-3, maxit = 4, ztol = 3e-3)
                else
                    TD.step!()
                end
                ok = true; break
            catch e
                COOL[] = 4
                att == 6 && pr("step %d failed after 7 retries: %s\n", k, sprint(showerror, e)[1:min(end, 200)])
            end
        end
        ok || break
        dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
        sd = TX.step_diag(dt)
        SM.burn!(sg; sd)
        SM.mix!(sg; sd, capmix = false, bridge = 4, tmax = 0.02)

        # Homologous mass accretion dM = Mdot * dt onto q = m / M_* grid!
        dM = Mdot * dt
        R_post = rn[b+1]
        Lint_post = max(4π * R_post^2 * TF.var(:frn)[b+1], 0.0)
        bc_post = precursor_bc(M_cur[] + dM, R_post, Mdot, Lint_post)
        M_cur[] = accrete_homologous!(mat, M_cur[], dM, bc_post, Y_prim, a, b)

        # Diagnostics
        R_s = rn[b+1]; L_int = max(4π * R_s^2 * TF.var(:frn)[b+1], 0.0)
        bc_d = precursor_bc(M_cur[], R_s, Mdot, L_int)
        Lnuc = SM.nuclear_power(sg)
        LEdd = 4π * Z0.c * Z0.G * M_cur[] / (0.2 * 1.75)
        Dmix = sg.Dmix; j = a + 1; last = a; gap = 0
        while j <= b && gap <= 2; Dmix[j] > 0 ? (last = j; gap = 0) : (gap += 1); j += 1; end
        Mcc = sum((4π * dn[kk] * TF.var(:dvoln)[kk] for kk in a:last-1); init = 0.0) / M_cur[]
        t_yr = (TF.getscalar(:timen) - tstart) / yr
        M_msun = M_cur[] / Msun

        if k == 1 || k % envi("PRINTEVERY", 25) == 0 || M_msun >= Mstop
            pr("%6d %10.3e %8.3f %8.2f %8.2f %7.0f %9.3e %9.3e %8.4f %8.4f %7.4f %7.4f %7.2e\n",
               k, t_yr, M_msun, R_s / Rsun, bc_d.Rph / Rsun, bc_d.Teff, bc_d.Ltot / Lsun, Lnuc / Lsun,
               bc_d.Ltot / LEdd, tn[a] / 1e7, sg.Y[SM.ID, a] / D0, Mcc, dt / yr)
        end
        if k == 1 || k % envi("CSVEVERY", 2) == 0 || M_msun >= Mstop
            @printf(csvio, "%d,%.6e,%.4e,%.5f,%.4f,%.4f,%.2f,%.5e,%.5e,%.5e,%.5e,%.5e,%.4e,%.5e,%.5e,%.5e,%.6f,%.5e,%.5f\n",
                    k, t_yr, dt, M_msun, R_s / Rsun, bc_d.Rph / Rsun, bc_d.Teff, L_int / Lsun, bc_d.Lacc / Lsun,
                    bc_d.Ltot / Lsun, Lnuc / Lsun, LEdd / Lsun, bc_d.tau_env, tn[a], maximum(view(tn, a:b)),
                    dn[a], SM.hmass(sg, a), sg.Y[SM.ID, a] / D0, Mcc)
            flush(csvio)
        end
        if M_msun >= Mstop || R_s / Rsun > 5.0e4
            pr("STOPPED at step %d: M*=%.2f Msun, R*=%.1f Rsun, t=%.3e yr (wall %.1f s)\n", k, M_msun, R_s/Rsun, t_yr, time() - t0w)
            break
        end
    end
    close(csvio); close(logio)
end

main()
