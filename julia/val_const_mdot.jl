#!/usr/bin/env julia
"""
    val_const_mdot.jl

Controlled constant-Mdot spherical RHD protostellar accretion validation driver on `popiii-val`,
built to reproduce Omukai & Palla (2003, OP03) and Hosokawa, Omukai & Yorke (2012, HOY12)
quantitatively using Titan.jl's unified FreeEnergyEOS microphysics + Stage-2 TDC + NucNet + StellarOpacity.

Physics stack:
  - Unified C² FreeEnergyEOS: self-generated C² B-spline LTE table (X=0.75, Y=0.25) below 3e5 K
    blended C²-quintically into exact FreeEnergyEOS.eos_full_ion(ρ, T, comp) above 9e5 K
    (variable NucNet X/Y/Z, Fermi–Dirac e± pairs, PC00+KSDT14 Coulomb, exact atomic energy zero).
  - StellarOpacity (κ_R, κ_P; H⁻, CIA, Kramers, OPAL, Thomson, electron conduction) with live per-zone X.
  - NucNet 19-species nuclear burning (D/H = 2.5e-5, ³He/H = 1e-5, pp, 3α, CNO) + TDC turbulent mixing
    + exact conservative donor-cell nuclear advection with pristine primordial inflow at r_out.
  - Full 1D spherical RHD (r, m, ρ, u, T, E_r, F_r, w): fixed logarithmic Eulerian mesh from regular
    centre (r_1 = 0, r_2 = RIN) to ROUT, leobc = 2 constant-Mdot supersonic inflow outer BC
    (phir0 = -Mdot/(4π r_out²) strictly constant; ur, dr, phir2 updated smoothly with enclosed mass),
    XMZERO centre-referenced mass coordinate, ECOMPAT first-law gas row, TR-BDF2 with error control.

Usage:
  FE_EOS_CACHE=/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls \
  MDOT=4.41e-3 M0=0.5 R0=35.0 MSTOP=150 NZ=220 ROUT=1e5 TAG=op03_fid \
  julia -t 4 --project=/zpool/nvme/data/tabel_scratch/titan_popiii/subconv/jlenv val_const_mdot.jl
"""

ENV["TITAN_NOLIB"] = "1"
haskey(ENV, "TITAN_MEQN") || (ENV["TITAN_MEQN"] = "8")
haskey(ENV, "TITAN_MGR")  || (ENV["TITAN_MGR"]  = "500")
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")
haskey(ENV, "DTERR") || (ENV["DTERR"] = "0.03")
haskey(ENV, "DTERRTRACE") || (ENV["DTERRTRACE"] = "0")
haskey(ENV, "DTERRBE") || (ENV["DTERRBE"] = "1")

using Printf, ForwardDiff, Serialization

const JDIR = @__DIR__
for f in ("TitanDriver", "JuliaPhysics", "TitanResidFast", "BlockSolver", "TitanCore", "Z0Physics",
          "TitanMaterials", "ProblemSetup", "StellarEOS", "FreeEnergyEOS", "StellarOpacity",
          "NucNet", "StellarMaterial", "StarMaterial", "GridEq", "Trbdf2", "NewtonKrylov")
    include(joinpath(JDIR, f * ".jl"))
end

const TD = TitanDriver; const TF = TD.TF; const TX = TitanResidFast; const PS = ProblemSetup; const TC = TitanCore
const TM = TitanMaterials; const Z0 = Z0Physics; const SE = StellarEOS; const FE = FreeEnergyEOS
const SO = StellarOpacity; const NN = NucNet; const SM = StarMaterial

const Msun = 1.98847e33; const Lsun = 3.828e33; const Rsun = 6.957e10; const yr = 3.15576e7

# ---------------------------------------------------------------------------------------------------------------------
# Unified C² AccretionGas Material: FreeEnergyEOS (C² B-spline LTE + eos_full_ion) + StellarOpacity + NucNet + TDC
# ---------------------------------------------------------------------------------------------------------------------
const UE_T1 = 3.0e5
const UE_T2 = 9.0e5
@inline _smoother5(x) = (x = clamp(x, zero(x), one(x)); x^3 * (10 - 15x + 6x^2))
@inline _fval(x::ForwardDiff.Dual) = _fval(ForwardDiff.value(x))
@inline _fval(x) = Float64(x)

# Preload our self-generated C² B-spline LTE table once!
const _FE_TAB = FE.lte_table()

mutable struct AccretionGas{T} <: TM.Material
    sg::SM.StarGas{T}
    kcap::Float64
end

"Exact unified C² FreeEnergyEOS (p, e) at (ρ, T, comp) in StarGas atomic energy zero (neutral H = 0, neutral He = 0)."
@inline function _acc_pe_raw(ρ, T, cp)
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
        u1, p1 = FE.lte_eos(ρ, T)
        e1 = u1 + (0.75 / FE.mH) * FE.χH
        u2, p2 = FE.eos_full_ion(ρ, T, cp; pairs = true, coulomb = true)
        e2 = u2 + ΔH
        return p1 + w * (p2 - p1), e1 + w * (e2 - e1)
    end
end

"Chain-rule lifted _acc_pe (zero nested 48-wide Dual compilation inside EOS!)."
function _acc_pe(ρ, T, cp)
    if !(ρ isa ForwardDiff.Dual || T isa ForwardDiff.Dual)
        return _acc_pe_raw(ρ, T, cp)
    end
    D = promote_type(typeof(ρ), typeof(T))
    ρd = D(ρ); Td = D(T)
    ρv = ForwardDiff.value(ρd); Tv = ForwardDiff.value(Td)
    if ρv isa ForwardDiff.Dual || Tv isa ForwardDiff.Dual
        return _acc_pe_raw(ρd, Td, cp)
    end
    ρ0 = Float64(ρv); T0 = Float64(Tv)
    p0, e0 = _acc_pe_raw(ρ0, T0, cp)
    J = ForwardDiff.jacobian(v -> begin; p_, e_ = _acc_pe_raw(v[1], v[2], cp); [p_, e_]; end, [ρ0, T0])
    Pρ = ForwardDiff.partials(ρd); PT = ForwardDiff.partials(Td)
    return D(p0, J[1, 1] * Pρ + J[1, 2] * PT), D(e0, J[2, 1] * Pρ + J[2, 2] * PT)
end

@inline function TX.material_gas(m::AccretionGas, ρ, T)
    cp = SM.zonecomp(m.sg)
    p, e = _acc_pe(ρ, T, cp)
    return p, e, sqrt(5.0 / 3.0 * max(p, 1e-30 * one(p)) / ρ)
end

@inline function _acc_ne(ρ, T, cp)
    Tv = _fval(T)
    if Tv >= UE_T2
        return ρ * cp.νe
    elseif Tv <= UE_T1
        yH = FE.lte_composition(ρ, T)
        z1, z2 = FE.lte_composition_he(ρ, T)
        return (ρ * 0.75 / FE.mH) * max(yH[1] + z1 + 2z2, zero(yH[1]))
    else
        w = _smoother5(log(T / UE_T1) / log(UE_T2 / UE_T1))
        yH = FE.lte_composition(ρ, T)
        z1, z2 = FE.lte_composition_he(ρ, T)
        ne1 = (ρ * 0.75 / FE.mH) * max(yH[1] + z1 + 2z2, zero(yH[1]))
        ne2 = ρ * cp.νe
        return ne1 + w * (ne2 - ne1)
    end
end

@inline function TX.material_opac(m::AccretionGas, ρ, T)
    X = SM.hmass(m.sg, TX.ZONE[])
    κR = SM.softcap(SO.kappa_R(ρ, T, X), m.kcap)
    κP = SM.softcap(SO.kappa_P(ρ, T, X), m.kcap)
    return κR, κP, κP
end

TX.has_source(m::AccretionGas) = m.sg.nuc
TX.material_source(m::AccretionGas, ρ, T) = begin
    k = TX.ZONE[]
    (m.sg.nuc && _fval(T) > UE_T1) ? (m.sg.fcoupling ? m.sg.fcur[k] : 1.0) * SM.eps_nuc(ρ, T, view(m.sg.Y, :, k); screen = m.sg.screen) : zero(T)
end

function TX.material_lrzone(m::AccretionGas, k, ρ, T, x, rs, dual::Bool)
    cp = SM.zonecomp(m.sg, k)
    ρv = _fval(ρ); Tv = _fval(T)
    f1(v) = begin; p_, e_ = _acc_pe_raw(v[1], v[2], cp); [p_, e_]; end
    p0, u0 = _acc_pe_raw(ρv, Tv, cp)
    J = ForwardDiff.jacobian(f1, [ρv, Tv])
    D = promote_type(typeof(ρ), typeof(T))
    z = zero(D)
    if D <: ForwardDiff.Dual
        h = 1e-5
        Jrp = ForwardDiff.jacobian(f1, [ρv * exp(h), Tv]); Jrm = ForwardDiff.jacobian(f1, [ρv * exp(-h), Tv])
        Jtp = ForwardDiff.jacobian(f1, [ρv, Tv * exp(h)]); Jtm = ForwardDiff.jacobian(f1, [ρv, Tv * exp(-h)])
        dJ_dρ = (Jrp .- Jrm) ./ (2h * ρv)
        dJ_dT = (Jtp .- Jtm) ./ (2h * Tv)
        ρd = D(ρ); Td = D(T)
        Pρ = ForwardDiff.partials(ρd); PT = ForwardDiff.partials(Td)
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

function TM.fill_eos!(m::AccretionGas)
    dn, tn = TM.v(:dn), TM.v(:tn)
    pgn, egn, asn, xnen = TM.v(:pgn), TM.v(:egn), TM.v(:asn), TM.v(:xnen)
    @inbounds for k in TM.krange()
        ρ, T = dn[k], tn[k]
        (ρ > 0 && T > 0) || continue
        TX.ZONE[] = k
        cp = SM.zonecomp(m.sg, k)
        p, e = _acc_pe_raw(ρ, T, cp)
        pgn[k] = p; egn[k] = e; asn[k] = sqrt(5.0 / 3.0 * max(p, 1e-30) / ρ)
        xnen[k] = _acc_ne(ρ, T, cp)
    end
    nothing
end

function TM.fill_opac!(m::AccretionGas)
    dn, tn = TM.v(:dn), TM.v(:tn)
    chifn, xken, xkpn, xnen = TM.v(:chifn), TM.v(:xken), TM.v(:xkpn), TM.v(:xnen)
    @inbounds for k in TM.krange()
        ρ, T = dn[k], tn[k]
        (ρ > 0 && T > 0) || continue
        TX.ZONE[] = k
        χ, κE, κP = TX.material_opac(m, ρ, T)
        chifn[k] = χ; xken[k] = κE; xkpn[k] = κP
        xnen[k] = _acc_ne(ρ, T, SM.zonecomp(m.sg, k))
    end
    nothing
end

# T_of_e solver for StarMaterial.burn!/mix! using our exact unified FreeEnergyEOS
function T_solve_acc(ρ, e, T0, comp)
    lT = log(clamp(T0, 100.0, 1e9))
    err = Inf
    for _ in 1:40
        T = exp(lT)
        _, ev = _acc_pe_raw(ρ, T, comp)
        f = ev - e
        err = abs(f) / max(abs(e), FE.kB * T / FE.mH)
        err < 1e-12 && return T, err
        de_dlT = ForwardDiff.derivative(lt -> _acc_pe_raw(ρ, exp(lt), comp)[2], lT)
        (!isfinite(de_dlT) || de_dlT <= 0) && break
        δ = clamp(-f / de_dlT, -0.4, 0.4)
        lT = clamp(lT + δ, log(50.0), log(3e9))
    end
    T = exp(lT)
    _, ev = _acc_pe_raw(ρ, T, comp)
    return T, abs(ev - e) / max(abs(e), FE.kB * T / FE.mH)
end

# ---------------------------------------------------------------------------------------------------------------------
# Lane–Emden Polytrope (n = 3 or n = 1.5) Seed Core + Free-Fall Envelope Builder on Fixed Log Mesh
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

# Conservative upwind nuclear advection with pristine primordial inflow at outer boundary (j = n)
function advect_nuc!(sg::SM.StarGas, xme_prev::Vector{Float64}, Y_prim::Vector{Float64}, a::Int, b::Int)
    xme = TF.var(:xmen)
    n = b - a + 1
    F = [xme[a+j-1] - xme_prev[a+j-1] for j in 1:n+1]      # > 0 outward, < 0 inward
    M = [xme[a+j-1] - xme[a+j] for j in 1:n]
    M0 = [xme_prev[a+j-1] - xme_prev[a+j] for j in 1:n]
    lo = zeros(n); up = zeros(n); di = zeros(n)
    for j in 1:n
        di[j] = M[j] + max(F[j+1], 0.0) + max(-F[j], 0.0)
        j > 1 ? (lo[j] = -max(F[j], 0.0)) : (di[j] -= max(F[j], 0.0))
        j < n && (up[j] = -max(-F[j+1], 0.0))
    end
    # Note at j = n: inflow max(-F[n+1], 0) enters RHS as max(-F[n+1], 0) * Y_prim (no up[n] term)!
    cp_ = zeros(n); den = zeros(n)
    den[1] = di[1]; cp_[1] = up[1] / den[1]
    for j in 2:n; den[j] = di[j] - lo[j] * cp_[j-1]; cp_[j] = up[j] / den[j]; end
    Y = sg.Y
    fin_out = max(-F[n+1], 0.0)
    ys = Vector{Vector{Float64}}(undef, n)
    ys[1] = (M0[1] .* Y[:, a]) ./ den[1]
    for j in 2:n
        rhs_j = M0[j] .* Y[:, a+j-1]
        j == n && fin_out > 0 && (rhs_j .+= fin_out .* Y_prim)
        ys[j] = (rhs_j .- lo[j] .* ys[j-1]) ./ den[j]
    end
    for j in n-1:-1:1; ys[j] = ys[j] .- cp_[j] .* ys[j+1]; end
    for j in 1:n; Y[:, a+j-1] .= max.(ys[j], 0.0); end
    # Pin outer ghosts to pristine primordial composition; inner ghosts reflect centre
    for k in max(a-2, 1):a-1; Y[:, k] .= Y[:, a]; end
    for k in b+1:min(b+2, size(Y, 2)); Y[:, k] .= Y_prim; end
    copyto!(xme_prev, xme)
    nothing
end

"Locate accretion shock radius R_star, enclosed mass M_star, photosphere R_ph, and interior L_max."
function acc_diagnostics(mat::AccretionGas, sd, a::Int, b::Int, Mdot::Float64)
    rn, dn, tn, un, frn, ern, pgn, chifn, asn = (TF.var(s) for s in (:rn, :dn, :tn, :un, :frn, :ern, :pgn, :chifn, :asn))
    xmen = TF.var(:xmen); xmtot = TF.getscalar(:xmtot)
    m_enc(k) = 4π * (xmtot - xmen[k])      # enclosed mass at interface k [g]
    # Photosphere τ_R = 2/3 integrated from outside in:
    tau = zeros(b + 2)
    tt = 0.0
    for k in b:-1:a
        tt += chifn[k] * dn[k] * (rn[k+1] - rn[k])
        tau[k] = tt
    end
    kph = a
    for k in b:-1:a
        if tau[k] >= 2.0 / 3.0
            kph = k; break
        end
    end
    Rph = rn[kph+1]
    Lph = 4π * Rph^2 * max(frn[kph+1], 1e-30)
    Lsurf = 4π * rn[b+1]^2 * max(frn[b+1], 1e-30)
    Teff = (Lph / (4π * Rph^2 * Z0.σSB))^0.25
    # Accretion shock radius R_star: outermost interface where gas is hydrostatic / subsonic inside the infall shock
    ksh = a + 1
    for k in b:-1:a+1
        cs = max(asn[k-1], 1e4)
        if abs(un[k]) < 0.5 * cs && dn[k-1] > 1e-13
            ksh = k; break
        end
    end
    Rstar = rn[ksh]
    Mstar = m_enc(ksh)
    Mgrid = m_enc(b + 1)
    # Interior maximum luminosity (rad + conv) inside R_star:
    Lmax = 0.0; Tmax = 0.0; umax_cs = 0.0
    for k in a:max(ksh-1, a)
        Lk = 4π * rn[k+1]^2 * (frn[k+1] + sd.fcv[k+1])
        Lk > Lmax && (Lmax = Lk)
        tn[k] > Tmax && (Tmax = tn[k])
        mc = abs(un[k+1]) / max(asn[k], 1e4)
        mc > umax_cs && (umax_cs = mc)
    end
    Lacc = Z0.G * Mstar * Mdot / max(Rstar, 1e8)
    Lnuc = SM.nuclear_power(mat.sg)
    κ_es = 0.2 * (1.0 + 0.75)
    LEdd = 4π * Z0.c * Z0.G * Mstar / κ_es
    # Convective core mass fraction (inside R_star)
    D = mat.sg.Dmix; j = a + 1; last = a; gap = 0
    while j < ksh && gap <= 2
        D[j] > 0 ? (last = j; gap = 0) : (gap += 1); j += 1
    end
    Mcc = (m_enc(last)) / max(Mstar, 1e20)

    # Supergiant envelope pulsation & density-inversion diagnostics (inside R_ph):
    wtn = TX.WTN[]; par = TX.TDCP[]
    max_dlnrho = -Inf; r_inv = 0.0; max_Gam = 0.0; max_PtP = 0.0; Ekin_star = 0.0
    k_lim = clamp(max(ksh, kph), a + 2, b - 1)
    for k in a:k_lim
        rc_k = 0.5 * (rn[k] + rn[k+1])
        Pr_k = ern[k] / 3.0
        Pt_k = par.α_p * dn[k] * 0.5 * (wtn[k]^2 + wtn[k+1]^2)
        Ptot_k = max(pgn[k] + Pr_k, 1e-30)
        PtP = Pt_k / Ptot_k
        PtP > max_PtP && (max_PtP = PtP)
        mk_g = 4π * dn[k] * (rn[k+1]^3 - rn[k]^3) / 3.0
        k < ksh && (Ekin_star += 0.5 * mk_g * 0.5 * (un[k]^2 + un[k+1]^2))
        if k > a && k < k_lim
            rc_p = 0.5 * (rn[k+1] + rn[k+2])
            dlr = log(rc_p / rc_k)
            if dlr > 0 && dn[k] > 1e-13 && dn[k+1] > 1e-13
                g_rho = log(dn[k+1] / dn[k]) / dlr
                if g_rho > max_dlnrho
                    max_dlnrho = g_rho; r_inv = rn[k+1]
                end
            end
        end
        m_k = m_enc(k + 1)
        if m_k > 0.01 * Mstar
            Lk = 4π * rn[k+1]^2 * max(frn[k+1], 0.0)
            Gam_k = chifn[k] * Lk / (4π * Z0.c * Z0.G * m_k)
            Gam_k > max_Gam && (max_Gam = Gam_k)
        end
    end
    isfinite(max_dlnrho) || (max_dlnrho = 0.0)
    u_ph_kms = un[clamp(kph + 1, a, b + 1)] / 1e5

    return (; Mstar = Mstar / Msun, Mgrid = Mgrid / Msun, Rstar = Rstar / Rsun, Rph = Rph / Rsun,
              Teff, Lsurf = Lsurf / Lsun, Lph = Lph / Lsun, Lmax = Lmax / Lsun, Lacc = Lacc / Lsun,
              Lnuc = Lnuc / Lsun, LEdd = LEdd / Lsun, Tc = tn[a], Tmax, ρc = dn[a],
              Xc = SM.hmass(mat.sg, a), Dc = mat.sg.Y[SM.ID, a], Mcc, umax_cs, ksh = ksh - a + 1, kph = kph - a + 1,
              max_dlnrho, r_inv = r_inv / Rsun, max_Gam, max_PtP, u_ph_kms, Ekin_star)
end

"Dump full radial structure profile (for density inversion, super-Eddington, TDC, and pulsation eigenmode analysis)."
function dump_profile(path::String, mat::AccretionGas, sd, a::Int, b::Int, D0::Float64)
    rn, dn, tn, un, frn, ern, pgn, chifn, asn = (TF.var(s) for s in (:rn, :dn, :tn, :un, :frn, :ern, :pgn, :chifn, :asn))
    xmen = TF.var(:xmen); xmtot = TF.getscalar(:xmtot)
    wtn = TX.WTN[]; par = TX.TDCP[]; Dmix = mat.sg.Dmix
    tau = zeros(b + 2); tt = 0.0
    for k in b:-1:a; tt += chifn[k] * dn[k] * (rn[k+1] - rn[k]); tau[k] = tt; end
    open(path, "w") do io
        println(io, "zone,r_Rsun,dr_Rsun,m_Msun,rho,T,u_kms,cs_kms,P_gas,P_rad,P_turb,beta,kappa_R,tau_R,Gamma_Edd,L_rad_Lsun,L_conv_Lsun,w_conv_kms,Dmix,nabla,nabla_ad,dlnrho_dlnr,X_H,D_over_D0")
        for k in a:b
            j = k - a + 1
            rc = 0.5 * (rn[k] + rn[k+1]); dr = rn[k+1] - rn[k]
            m_k = 4π * (xmtot - xmen[k+1])
            Pr = ern[k] / 3.0; Pg = pgn[k]
            Pt = par.α_p * dn[k] * 0.5 * (wtn[k]^2 + wtn[k+1]^2)
            Ptot = max(Pg + Pr, 1e-30); beta = Pg / Ptot
            Lrad = 4π * rn[k+1]^2 * frn[k+1]
            Lconv = 4π * rn[k+1]^2 * sd.fcv[k+1]
            Gam = m_k > 0 ? chifn[k] * max(Lrad, 0.0) / (4π * Z0.c * Z0.G * m_k) : 0.0
            # Local nabla and exact nabla_ad (gas + rad) via 2x2 AD on _acc_pe_raw:
            nab = NaN; dlnrho = NaN
            if k < b
                Ptot_p = max(pgn[k+1] + ern[k+1] / 3.0, 1e-30)
                dlP = log(Ptot_p / Ptot)
                abs(dlP) > 1e-12 && (nab = log(tn[k+1] / tn[k]) / dlP)
                rc_p = 0.5 * (rn[k+1] + rn[k+2])
                dlnrho = log(dn[k+1] / dn[k]) / log(rc_p / rc)
            end
            cp = SM.zonecomp(mat.sg, k)
            J = ForwardDiff.jacobian(v -> begin
                p_, e_ = _acc_pe_raw(v[1], v[2], cp)
                [p_ + Z0.arad * v[2]^4 / 3.0, e_ + Z0.arad * v[2]^4 / v[1]]
            end, [dn[k], tn[k]])
            pρ, pT, uρ, uT = J[1, 1], J[1, 2], J[2, 1], J[2, 2]
            cP_tot = uT + (tn[k] / dn[k]^2) * pT^2 / max(pρ, 1e-30)
            nad = (Ptot / max(tn[k], 1.0)) * pT / max(dn[k]^2 * pρ * cP_tot, 1e-30)
            @printf(io, "%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.5e,%.5e,%.6e,%.6e,%.6e,%.5f,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.6f,%.5e\n",
                    j, rc / Rsun, dr / Rsun, m_k / Msun, dn[k], tn[k], un[k+1] / 1e5, asn[k] / 1e5,
                    Pg, Pr, Pt, beta, chifn[k], tau[k], Gam, Lrad / Lsun, Lconv / Lsun,
                    wtn[k+1] / 1e5, Dmix[k+1], isfinite(nab) ? nab : 0.0, isfinite(nad) ? nad : 0.0,
                    isfinite(dlnrho) ? dlnrho : 0.0, SM.hmass(mat.sg, k), mat.sg.Y[SM.ID, k] / D0)
        end
    end
    nothing
end

"Configure TitanResidFast.GRIDEQ[] for the Dorfi-Drury implicit moving mesh (tracks shock + photosphere smoothly)."
function grideq_update!(a::Int, b::Int, k_rel::Int; κ_target = 30.0, κ_ramp0 = 2.0, n_ramp = 60,
                        fv = 0.5, fvr = 10.0, w_tau = 1.0, w_chi = 3.0, w_u = 1.0, fdyn = 0.2,
                        pr_dec = 6.0, pm_dec = 6.0, amix = 0.4)
    rn, un, asn = TF.var(:rn), TF.var(:un), TF.var(:asn)
    n = b - a + 1; mgr = length(rn)
    κ = κ_ramp0 * (κ_target / κ_ramp0)^clamp(k_rel / max(n_ramp, 1), 0.0, 1.0)
    tauv = zeros(mgr); tfl = zeros(mgr); sx = zeros(mgr)
    flc = 1.0 / min(fv * κ_target / κ, fvr)
    @inbounds for j in 1:n
        k = a + j - 1
        rc = 0.5 * (rn[k] + rn[k+1])
        um = 0.5 * (abs(un[k]) + abs(un[k+1]))
        tc = rc / max(um + asn[k], 1.0)
        dxj = j >= 2 ? log(rn[k+1] / rn[k]) : 1.0
        tfl[k] = dxj * tc * flc
        tauv[k] = κ * tc
        sx[k] = w_u > 0 ? w_u * ((un[k+1] - un[k]) / max(um + asn[k], 1.0))^2 : 0.0
    end
    fr = (pr_dec / log(10))^2; fm = (pm_dec / log(10))^2
    TX.GRIDEQ[] = (; k1 = a + 1, k2 = b + 1, w = (1.0, 1.0, 1.0, 1.0, 1.0), alph = 1.5, gmax = 100.0,
                   tauv, tfl, flc, wtau = w_tau, tlo = 0.1, thi = 1e4, kap = κ, fdyn, sx,
                   wχ = w_chi, gp = (0.0, 0.0), fr, fm, amix)
    nothing
end

# ---------------------------------------------------------------------------------------------------------------------
# Main Driver
# ---------------------------------------------------------------------------------------------------------------------
function main()
    envf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
    envi(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d

    Mdot_yr = envf("MDOT", 4.41e-3)
    Mdot    = Mdot_yr * Msun / yr
    M0_msun = envf("M0", 0.5)
    R0_rsun = envf("R0", 35.0)
    Mstop   = envf("MSTOP", 120.0)
    rin     = envf("RIN", 0.2) * Rsun
    rout    = envf("ROUT", 1.0e5) * Rsun
    nz      = envi("NZ", 220)
    Tinflow = envf("TINFLOW", 300.0)
    npoly   = envf("NPOLY", 1.5)            # n=1.5 adiabatic/convective protostellar seed (or 3.0)
    maxsteps = envi("MAXSTEPS", 10000)
    outdir  = get(ENV, "VAL_OUT", "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs")
    tag     = get(ENV, "TAG", @sprintf("mdot_%.2e", Mdot_yr))
    mkpath(outdir)
    logpath = joinpath(outdir, "$(tag).log")
    csvpath = joinpath(outdir, "$(tag).csv")
    logio = open(logpath, "w")
    csvio = open(csvpath, "w")
    pr(f, args...) = (s = Printf.format(Printf.Format(f), args...); print(s); print(logio, s); flush(stdout); flush(logio))

    # 1. Build fixed logarithmic radial mesh from regular centre (r_1 = 0, r_2 = rin .. rout)
    lr = collect(range(log(rin), log(rout); length = nz))
    rI = vcat(0.0, exp.(lr))
    n = length(rI) - 1   # nz zones

    # 2. Construct seed protostar (polytrope M0, R0) + steady free-fall inflow envelope at constant Mdot
    M0 = M0_msun * Msun; R0 = R0_rsun * Rsun
    ξ, θp, μp = lane_emden(npoly); ξ1 = ξ[end]; μ1 = μp[end]
    al = R0 / ξ1
    ρc0 = M0 / (4π * al^3 * μ1)
    K0 = (4π * Z0.G / (npoly + 1.0)) * al^2 * ρc0^(1.0 - 1.0 / npoly)
    θof(x) = (x >= ξ1 && return 0.0; i = clamp(searchsortedlast(ξ, x), 1, length(ξ) - 1);
              θp[i] + (x - ξ[i]) / (ξ[i+1] - ξ[i]) * (θp[i+1] - θp[i]))

    ρz = zeros(n); Tz = zeros(n); uI = zeros(n + 1); Fz = zeros(n + 1)
    μ_mol = 0.62 * SE.amu
    Lacc0 = Z0.G * M0 * Mdot / R0
    Teff0 = (Lacc0 / (4π * R0^2 * Z0.σSB))^0.25
    # Invert FreeEnergyEOS P(ρ, T) = P_target cleanly for core zones
    function T_from_P(ρ_, Pt_)
        lo, hi = 200.0, 2e8
        for _ in 1:60
            m = sqrt(lo * hi)
            p_, _ = _acc_pe_raw(ρ_, m, SE.composition(0.75, 0.25, 0.0))
            p_ + Z0.arad * m^4 / 3 > Pt_ ? (hi = m) : (lo = m)
            hi / lo - 1 < 1e-6 && break
        end
        sqrt(lo * hi)
    end
    for j in 1:n
        rc = cbrt(0.5 * (rI[j]^3 + rI[j+1]^3))
        if rc < R0
            th = max(θof(rc / al), 1e-4)
            ρ_ = ρc0 * th^npoly
            P_ = K0 * ρ_^(1.0 + 1.0 / npoly)
            ρ_ff = Mdot / (4π * rc^2 * sqrt(2 * Z0.G * M0 / rc))
            if ρ_ > 2.0 * ρ_ff
                ρz[j] = ρ_
                Tz[j] = max(T_from_P(ρ_, P_), Teff0)
            else
                ρz[j] = ρ_ff
                Tz[j] = max(Teff0 * sqrt(R0 / (2rc)), Tinflow)
            end
        else
            vff = sqrt(2 * Z0.G * M0 / rc)
            ρz[j] = Mdot / (4π * rc^2 * vff)
            # Smooth radiative equilibrium envelope temperature above Tinflow:
            Tz[j] = max(Teff0 * sqrt(R0 / (2rc)), Tinflow)
        end
    end
    for j in 2:n+1
        r = rI[j]
        if r > R0
            w = _smoother5((r - R0) / (0.3 * R0))
            uI[j] = -w * sqrt(2 * Z0.G * M0 / r)
            Fz[j] = Lacc0 / (4π * r^2)
        else
            uI[j] = 0.0
            Fz[j] = Lacc0 * (r / R0)^3 / (4π * max(r, rin)^2)
        end
    end
    # Outer inflow BC state at r_out:
    vff_out = sqrt(2 * Z0.G * M0 / rout)
    ρ_out = Mdot / (4π * rout^2 * vff_out)
    u_out = -vff_out
    uI[end] = u_out

    # 3. Build AccretionGas material with primordial D/H = 2.5e-5, ³He/H = 1e-5
    DH = envf("DH", 2.5e-5)
    sg = SM.StarGas(TF.MGR; X = 0.75, Y = 0.25, Z = 0.0, opac = :stellar, nuc = true,
                    dXmax = envf("DXMAX", 4e-3), efrac = envf("EFRAC", 5.0), burnT = :path, fcoupling = true)
    Y_prim = zeros(NN.NS)
    Yp = 0.75 / (1.0 + 2DH)
    Y_prim[NN.IDX[:p]] = Yp
    Y_prim[NN.IDX[:d]] = DH * Yp
    Y3 = 1e-5 * Yp
    Y_prim[NN.IDX[:he3]] = Y3
    Y_prim[NN.IDX[:he4]] = (0.25 - 3Y3) / 4.0
    for k in axes(sg.Y, 2); sg.Y[:, k] .= Y_prim; end
    mat = AccretionGas(sg, envf("KCAP", Inf))

    # 4. Initialize TITAN state (lgrid = 2 fixed Eulerian, outer = :inflow leobc = 2, regular centre)
    Er = [begin
              rc = cbrt(0.5 * (rI[j]^3 + rI[j+1]^3))
              rc <= R0 ? Z0.arad * Tz[j]^4 : max(Lacc0 / (4π * rc^2 * Z0.c), Z0.arad * 20.0^4)
          end for j in 1:n]
    av = (ql0 = 0.0, ql1 = envf("AVQL1", 0.15), cq1 = 0.0, cq2 = envf("AVCQ2", 1.0))
    PS.build!(; rI, rho = ρz, T = Tz, u = uI, geom = 2, rad = true, Er, Fr = Fz,
              lgrid = 2, inner = :wall, outer = :inflow, inflow = (ρ_out, Tinflow, u_out),
              selfgrav = true, material = mat, convection = true, av, xmext = 10.0 * Mstop * Msun / (4π))
    a = Int(TF.getscalar(:ngrs)); b = Int(TF.getscalar(:ngre))

    # Set XMZERO (centre-referenced mass coordinate for 16-digit relative precision everywhere!)
    let xmen = TF.var(:xmen), xmeo = TF.var(:xmeo), dn = TF.var(:dn), dv = TF.var(:dvoln)
        m = zeros(length(xmen))
        for k in a+1:b+2; m[k] = m[k-1] + dn[k-1] * dv[k-1]; end
        ε = 0.5 * m[a+1]
        for k in a:b+2; xmen[k] = ε - m[k]; end
        xmen[a-1] = 2xmen[a] - xmen[a+1]
        xmeo .= xmen
        TF.setscalar!(:xmtot, ε); TF.setscalar!(:xmext, ε - m[b+2])
    end

    # Fixed Eulerian RHD + Stage-2 TDC + ECOMPAT + VEFIMP + FLUXTIE + SLABFLUX
    TX.REGULAR_CENTRE[] = true
    empty!(TX.XMINT[])
    TX.GASROW[] = true
    TX.AV_ZONE[] = 2
    TX.AV_TENSION[] = false
    TX.ECOMPAT[] = true
    TX.TIEFRAME[] = true
    TX.SLABFLUX[] = true
    TX.SLABTIE[] = true
    TD.VEF_ITER[] = true
    TD.VEF_SPH[] = true
    TD.VEF_BG[] = true
    TD.FLUXTIE[] = (0.1, 1.0)
    TX.ALPHA_MLT[] = envf("ALPHAMLT", 1.5)
    TX.TDC_ON[] = true
    TF.setscalar!(:neqn, 8)
    TX.tdc_alloc!(length(TF.var(:dn)))
    TX.TDCP[] = TX.TDC.TDCPar(; A_s = 1.0, A_d = 1.0, A_c = 1.0, α_p = envf("ALPHAP", 1.0), α_t = 0.01, α_D = 1.0, α_m = 0.0)
    TX.TDC_CEN[] = 0.0
    TC.UNOM_FLOOR[] = 1.0e-2
    TC.CPLH[] = (; newstep = () -> TX.tdc_newstep!(), reset = () -> TX.tdc_reset!(),
                   update = (rhs, a_, b_) -> TX.tdc_update!(rhs, a_, b_),
                   dmax = (rhs, a_, b_) -> TX.tdc_dmax!(rhs, a_, b_))

    PS.install_julia!(BlockSolver)
    if get(ENV, "NKRYLOV", "1") != "0"
        NewtonKrylov.configure!(get(ENV, "NKRYLOV", "1"); activate = true)
    end
    TD.MSGIO[] = devnull
    TF.setscalar!(:thet, 1.0)
    TF.setscalar!(:stol, 0.3); TF.setscalar!(:tmax, 1e30); TF.setscalar!(:niter, 25); TF.setscalar!(:conv, envf("CONV", 1e-3))
    TD.FAIL_TFAC[] = 0.5; Trbdf2.TFAC[] = 0.5
    TC.PRED_LAGRHO[] = true; TC.EXTRAP_FRAC[] = 0.3; TC.ER_TAUMIN[] = 3000.0

    # Optional GRIDEQ=1 Dorfi-Drury implicit moving mesh (tracks accretion shock & photosphere smoothly!)
    use_grideq = get(ENV, "GRIDEQ", "0") == "1"
    geq_kap    = envf("GEQKAP", 30.0)
    geq_kap0   = envf("GEQKAP0", 2.0)
    geq_nramp  = envi("GEQRAMPN", 60)
    geq_fmove  = envf("GEQFMOVE", 0.4)
    geq_wu     = envf("GEQWU", 1.0)
    geq_wtau   = envf("GEQWTAU", 1.0)
    geq_wchi   = envf("GEQWCHI", 3.0)
    geq_ctol   = envf("GEQCTOL", 0.3)

    # Wire StarMaterial nuclear burning / mixing hooks to our unified AccretionGas
    SM.TSOLVE[] = T_solve_acc
    SM.ZONEOK[] = k -> TF.var(:tn)[k] > UE_T1
    SM.REGION_CAP[] = true
    SM.install_stage_hook!(sg)
    nuc_cap() = min(sg.dtcap, SM.predict_cap(sg))
    empty!(TC.DT_HOOKS); push!(TC.DT_HOOKS, nuc_cap)
    empty!(Trbdf2.DT_BOUNDS); push!(Trbdf2.DT_BOUNDS, nuc_cap)
    if use_grideq
        fmcap = function ()
            TX.GRIDEQ[] === nothing && return Inf
            rn_ = TF.var(:rn); ro_ = TF.var(:ro); f = 0.0
            for kk in a+2:b
                wd = min(ro_[kk] - ro_[kk-1], ro_[kk+1] - ro_[kk])
                wd > 0 && (f = max(f, abs(rn_[kk] - ro_[kk]) / wd))
            end
            f > 0 ? geq_fmove / f * TF.getscalar(:dtime) : Inf
        end
        push!(TC.DT_HOOKS, fmcap)
        push!(Trbdf2.DT_BOUNDS, fmcap)
    end

    # Relaxation rule (Backward Euler for initial ~3 tdyn to damp polytrope acoustic modes cleanly);
    # once relaxed, TR-BDF2's Hosea-Shampine embedded error controller (DTERR) drives dt automatically!
    tdyn0 = sqrt(R0^3 / (Z0.G * M0))
    relaxed = Ref(false)
    TC.DT_RULE[] = function (dt0)
        relaxed[] && return dt0 * 1.3   # diagnostic only under DTERR
        min(1.35 * dt0, 0.5 * tdyn0)
    end

    # Initial dt: 0.02 * seed dynamical time
    TF.setscalar!(:dtime, envf("DT0", 0.02 * tdyn0))

    TD.EOS[](); TD.OPAC[]()
    TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
    xme_prev = copy(TF.var(:xmen)::Vector{Float64})
    phir0_const = -Mdot / (4π * rout^2)

    # Checkpoint save / restart support
    ckdir = joinpath(outdir, "ck")
    save_name = get(ENV, "SAVE", "")
    rest_name = get(ENV, "RESTART", "")
    ckevery   = envi("CKEVERY", 250)
    profevery = envi("PROFEVERY", 0)
    profmasses = haskey(ENV, "PROFMASS") && !isempty(ENV["PROFMASS"]) ? sort(parse.(Float64, split(ENV["PROFMASS"], ","))) : Float64[]
    k_start = 0; k_geq0 = Ref(-1)

    function save_ck(name::String, k_cur::Int)
        mkpath(ckdir)
        TF.save_state(joinpath(ckdir, "$(name).tstate.gz"))
        serialize(joinpath(ckdir, "$(name).aux"), (; Y = copy(sg.Y), Yprev = copy(sg.Yprev),
                  fscale = copy(sg.fscale), fslope = copy(sg.fslope), Dmix = copy(sg.Dmix),
                  dtcap = sg.dtcap, xme_prev = copy(xme_prev), wtn = copy(TX.WTN[]), wto = copy(TX.WTO[]),
                  k = k_cur, relaxed = relaxed[], k_geq0 = k_geq0[]))
    end
    if !isempty(rest_name)
        TF.load_state!(joinpath(ckdir, "$(rest_name).tstate.gz"))
        aux = deserialize(joinpath(ckdir, "$(rest_name).aux"))
        sg.Y .= aux.Y; sg.Yprev .= aux.Yprev; sg.fscale .= aux.fscale; sg.fslope .= aux.fslope
        sg.Dmix = copy(aux.Dmix); sg.dtcap = aux.dtcap
        copyto!(xme_prev, aux.xme_prev); copyto!(TX.WTN[], aux.wtn); copyto!(TX.WTO[], aux.wto)
        k_start = aux.k; relaxed[] = aux.relaxed; k_geq0[] = aux.k_geq0
        # Re-apply runtime tolerances / dt cap on restart
        TF.setscalar!(:stol, 0.3); TF.setscalar!(:tmax, 1e30); TF.setscalar!(:niter, 25); TF.setscalar!(:conv, envf("CONV", 1e-3))
        haskey(ENV, "DT0") && TF.setscalar!(:dtime, envf("DT0", TF.getscalar(:dtime)))
        TD.EOS[](); TD.OPAC[]()
        TC._ghost_and_derived!(:update, () -> TD.EOS[](), () -> TD.OPAC[]())
        pr("RESTARTED from %s at step %d (t=%.3e yr, dt=%.2e s)\n", rest_name, k_start, TF.getscalar(:timen)/yr, TF.getscalar(:dtime))
    end

    pr("VAL_CONST_MDOT [%s]: Mdot=%.3e Msun/yr, M0=%.2f Msun, R0=%.1f Rsun, Mstop=%.1f Msun, nz=%d, r=%.2f..%.1e Rsun, FREE_EOS + TDC(αp=%.1f) + DTERR(%s) + NKRYLOV + GRIDEQ=%d\n",
       tag, Mdot_yr, M0_msun, R0_rsun, Mstop, n, rin/Rsun, rout/Rsun, TX.TDCP[].α_p, ENV["DTERR"], use_grideq ? 1 : 0)
    println(csvio, "step,t_yr,dt_s,M_star,M_grid,R_star,R_ph,Teff,L_surf,L_ph,L_max,L_acc,L_nuc,L_Edd,Tc,Tmax,rho_c,Xc,Dc_over_D0,Mcc_over_M,umax_cs,ksh,kph,max_dlnrho_dlnr,r_inv_Rsun,max_Gamma_Edd,max_Pt_over_P,u_ph_kms,Ekin_star_erg")
    pr("%6s %10s %8s %8s %8s %8s %7s %9s %9s %8s %8s %7s %7s %7s\n",
       "step", "t [yr]", "M*", "Mgrid", "R* [Rs]", "Rph [Rs]", "Teff", "Lph [Ls]", "Lnuc", "L/LEdd", "Tc [1e7]", "dlnρ_m", "Γ_max", "dt [s]")

    STEP_CORE = () -> (relaxed[] ? Trbdf2.step!() : TD.step!())
    D0 = Y_prim[NN.IDX[:d]]
    t0w = time()
    for step_i in 1:maxsteps
        k = k_start + step_i
        if !relaxed[] && (TF.getscalar(:timen) >= 3.0 * tdyn0 || step_i >= 30)
            relaxed[] = true
        end
        if use_grideq && relaxed[]
            if k_geq0[] < 0
                k_geq0[] = k
                TC.GRIDMOVE[] = true; TF.setscalar!(:ctol, geq_ctol)
            end
            grideq_update!(a, b, k - k_geq0[]; κ_target = geq_kap, κ_ramp0 = geq_kap0, n_ramp = geq_nramp,
                           w_u = geq_wu, w_tau = geq_wtau, w_chi = geq_wchi)
        end
        # Solve coupled hydro + radiation + TDC (+ nuclear level factor once Tmax > 5e5 K)
        snap = TD.snapshot(); dt0 = TF.getscalar(:dtime)
        burn_active = maximum(view(TF.var(:tn), a:b)) > 5.0e5
        ok = false
        for att in 0:6
            TD.restore!(snap); TF.setscalar!(:dtime, dt0 * 0.4^att)
            try
                if burn_active
                    SM.solve_step!(sg, STEP_CORE; tol = 3e-3, maxit = 4, ztol = 3e-3)
                else
                    STEP_CORE()
                end
                ok = true; break
            catch e
                att == 6 && pr("step %d failed after 7 retries: %s\n", k, sprint(showerror, e)[1:min(end, 200)])
            end
        end
        ok || break
        dt = TF.getscalar(:timen) - TF.getscalar(:timeo)
        sd = TX.step_diag(dt)
        # Advect nuclear species conservatively (with pristine primordial inflow at r_out), then burn & mix!
        advect_nuc!(sg, xme_prev, Y_prim, a, b)
        SM.burn!(sg; sd)
        SM.mix!(sg; sd, capmix = false, bridge = 4, tmax = 0.02)

        # Gently update outer inflow (ur, dr, phir2) to track growing M_grid while keeping phir0 STRICTLY constant!
        M_tot = 4π * (TF.getscalar(:xmtot) - TF.var(:xmen)[b+1])
        vff_new = sqrt(2 * Z0.G * max(M_tot, M0) / rout)
        ur_new = -vff_new
        dr_new = -phir0_const / vff_new
        _, er_new = _acc_pe_raw(dr_new, Tinflow, SE.composition(0.75, 0.25, 0.0))
        TF.setscalar!(:ur, ur_new); TF.setscalar!(:dr, dr_new); TF.setscalar!(:tr, Tinflow)
        TF.setscalar!(:phir0, phir0_const); TF.setscalar!(:phir2, phir0_const * er_new)
        TF.var(:dn)[b+1] = dr_new; TF.var(:dn)[b+2] = dr_new
        TF.var(:tn)[b+1] = Tinflow; TF.var(:tn)[b+2] = Tinflow
        TF.var(:un)[b+1] = ur_new; TF.var(:un)[b+2] = ur_new
        for kk in b+1:min(b+2, size(sg.Y, 2)); sg.Y[:, kk] .= Y_prim; end

        # Re-fill EOS/opac with updated composition & outer ghost
        TD.EOS[](); TD.OPAC[]()

        diag = acc_diagnostics(mat, sd, a, b, Mdot)
        t_yr = TF.getscalar(:timen) / yr
        if step_i == 1 || step_i % envi("PRINTEVERY", 25) == 0 || diag.Mstar >= Mstop
            pr("%6d %10.3e %8.3f %8.3f %8.2f %8.2f %7.0f %9.3e %9.3e %8.4f %8.4f %7.2f %7.2f %7.2e\n",
               k, t_yr, diag.Mstar, diag.Mgrid, diag.Rstar, diag.Rph, diag.Teff, diag.Lph, diag.Lnuc,
               (diag.Lph) / max(diag.LEdd, 1.0), diag.Tc / 1e7, diag.max_dlnrho, diag.max_Gam, dt)
        end
        if step_i == 1 || step_i % envi("CSVEVERY", 5) == 0 || diag.Mstar >= Mstop
            @printf(csvio, "%d,%.6e,%.4e,%.5f,%.5f,%.4f,%.4f,%.2f,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.6f,%.5e,%.5f,%.4f,%d,%d,%.4f,%.4f,%.4f,%.4e,%.3f,%.5e\n",
                    k, t_yr, dt, diag.Mstar, diag.Mgrid, diag.Rstar, diag.Rph, diag.Teff, diag.Lsurf, diag.Lph,
                    diag.Lmax, diag.Lacc, diag.Lnuc, diag.LEdd, diag.Tc, diag.Tmax, diag.ρc, diag.Xc,
                    diag.Dc / D0, diag.Mcc, diag.umax_cs, diag.ksh, diag.kph,
                    diag.max_dlnrho, diag.r_inv, diag.max_Gam, diag.max_PtP, diag.u_ph_kms, diag.Ekin_star)
            flush(csvio)
        end
        if (profevery > 0 && (step_i == 1 || step_i % profevery == 0)) || (!isempty(profmasses) && diag.Mstar >= profmasses[1])
            pdir = joinpath(outdir, "profiles"); mkpath(pdir)
            ppath = joinpath(pdir, @sprintf("%s_s%05d_M%.1f.csv", tag, k, diag.Mstar))
            dump_profile(ppath, mat, sd, a, b, D0)
            while !isempty(profmasses) && diag.Mstar >= profmasses[1]; popfirst!(profmasses); end
        end
        if !isempty(save_name) && (step_i % ckevery == 0 || diag.Mstar >= Mstop)
            save_ck(save_name, k)
        end
        if diag.Mstar >= Mstop
            pr("Target M_star = %.2f Msun reached at step %d (t = %.3e yr, wall = %.1f s)\n", diag.Mstar, k, t_yr, time() - t0w)
            break
        end
    end
    !isempty(save_name) && save_ck(save_name, k_start + maxsteps)
    close(csvio); close(logio)
end

main()
