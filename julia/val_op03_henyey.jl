#!/usr/bin/env julia
"""
    val_op03_henyey.jl

Full-mass (0.5 -> 100+ Msun) hydrostatic Henyey + Sugimoto (1970) homologous-accretion
solver built directly on Titan.jl's unified microphysics:
  - Unified C² FreeEnergyEOS (self-generated C² B-spline LTE table below 3e5 K, quintically
    blended into FreeEnergyEOS.eos_full_ion above 9e5 K with variable NucNet composition,
    relativistic Fermi-Dirac e± pairs, and PC00+KSDT14 Coulomb).
  - StellarOpacity (κ_R) + 19-species NucNet (D/H = 2.5e-5, pp, 3α, CNO) + Kippenhahn-Weigert MLT.
  - Exact OP03 / Stahler (1986) shock outer BC (BCMODE=shock) or cold-disk photospheric BC (BCMODE=phot).
  - Exact Sugimoto (1970) / OP03 eq. (4) homologous accretion entropy advection in q = m / M_*(t).

Purpose:
  1. Reproduce the full-mass (0.5 -> 100 Msun) Omukai & Palla (2003, OP03) Mdot_crit ≈ 4e-3 Msun/yr
     bifurcation sweep (subcritical KH contraction to ZAMS vs supercritical Eddington expansion).
  2. Cross-validate directly against Titan.jl's full 1D Eulerian RHD runs (val_const_mdot.jl)
     at Mdot = 4.41e-3, 6e-3, 3e-2, 0.1, 0.3 Msun/yr and Titan.jl's ZAMS benchmarks (s5_run.jl).

Usage:
  MDOT=4.41e-3 M0=0.5 R0=25.0 MSTOP=100.0 TAG=op03_hen_4e3 \
  julia --project=/zpool/nvme/data/tabel_scratch/titan_popiii/subconv/jlenv val_op03_henyey.jl
"""

ENV["TITAN_NOLIB"] = "1"
haskey(ENV, "FE_EOS_CACHE") || (ENV["FE_EOS_CACHE"] = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/fe_lte_table_X75_Y25_v1.jls")

using Printf, LinearAlgebra, SparseArrays, ForwardDiff

const JDIR = @__DIR__
for f in ("Z0Physics", "StellarEOS", "FreeEnergyEOS", "StellarOpacity", "NucNet")
    include(joinpath(JDIR, f * ".jl"))
end

const Z0 = Z0Physics; const SE = StellarEOS; const FE = FreeEnergyEOS
const SO = StellarOpacity; const NN = NucNet

const G     = Z0.G
const c_cgs = Z0.c
const arad  = Z0.arad
const σSB   = Z0.σSB
const Msun  = 1.98847e33
const Lsun  = 3.828e33
const Rsun  = 6.957e10
const yr    = 3.15576e7

# ---------------------------------------------------------------------------------------------------------------------
# 1. Unified C² FreeEnergyEOS + Thermodynamics (ρ, e, cP, δ, ∇ad) from (P_tot, T, comp)
# ---------------------------------------------------------------------------------------------------------------------
const UE_T1 = 3.0e5
const UE_T2 = 9.0e5
@inline _smoother5(x) = (x = clamp(x, zero(x), one(x)); x^3 * (10 - 15x + 6x^2))

const _FE_TAB = FE.lte_table()

@inline function _pe_gas(ρ::Float64, T::Float64, cp)
    ΔH = (cp.x[1] / cp.m[1]) * FE.χH
    if T >= UE_T2
        u, p = FE.eos_full_ion(ρ, T, cp; pairs = true, coulomb = true)
        return p, u + ΔH
    elseif T <= UE_T1
        u, p = FE.lte_eos(ρ, T)
        return p, u + (0.75 / FE.mH) * FE.χH
    else
        w = _smoother5(log(T / UE_T1) / log(UE_T2 / UE_T1))
        u1, p1 = FE.lte_eos(ρ, T); e1 = u1 + (0.75 / FE.mH) * FE.χH
        u2, p2 = FE.eos_full_ion(ρ, T, cp; pairs = true, coulomb = true); e2 = u2 + ΔH
        return p1 + w * (p2 - p1), e1 + w * (e2 - e1)
    end
end

@inline function _pe_gas_dual(ρ, T, cp)
    Tv = ForwardDiff.value(T)
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

function eos_from_PT(Ptot::Float64, T::Float64, cp; ρ_guess::Float64 = 1e-6)
    Pr = arad * T^4 / 3.0
    Pg_tgt = max(Ptot - Pr, 1e-4 * Ptot)
    # Fast Newton in ln ρ with bisection safeguard
    lρ = log(clamp(ρ_guess, 1e-18, 1e5))
    ok = false
    for _ in 1:25
        ρ = exp(lρ)
        p, _ = _pe_gas(ρ, T, cp)
        f = log(max(p, 1e-40)) - log(Pg_tgt)
        if abs(f) < 1e-10; ok = true; break; end
        # finite-diff d(ln p)/d(ln ρ)
        h = 1e-4
        pp, _ = _pe_gas(ρ * exp(h), T, cp)
        dpdlr = max((log(max(pp, 1e-40)) - log(max(p, 1e-40))) / h, 0.2)
        lρ = clamp(lρ - clamp(f / dpdlr, -1.5, 1.5), log(1e-20), log(1e6))
    end
    if !ok
        lo, hi = 1e-20, 1e6
        for _ in 1:60
            rm = sqrt(lo * hi)
            pm, _ = _pe_gas(rm, T, cp)
            pm > Pg_tgt ? (hi = rm) : (lo = rm)
            hi / lo - 1 < 1e-9 && break
        end
        lρ = log(sqrt(lo * hi))
    end
    ρ = exp(lρ)
    # Exact local thermodynamic derivatives via central differences in ln ρ, ln T (smooth C²!)
    h = 1e-4
    p0, u0 = _pe_gas(ρ, T, cp)
    pp, up = _pe_gas(ρ * exp(h), T, cp); pm, um = _pe_gas(ρ * exp(-h), T, cp)
    ptp, utp = _pe_gas(ρ, T * exp(h), cp); ptm, utm = _pe_gas(ρ, T * exp(-h), cp)
    pρ = (pp - pm) / (2h * ρ)
    pT = (ptp - ptm) / (2h * T)
    uT = (utp - utm) / (2h * T)
    Pt = p0 + Pr
    Ptρ = max(pρ, 1e-30)
    PtT = max(pT + 4.0 * arad * T^3 / 3.0, 1e-30)
    cv = max(uT + 4.0 * arad * T^3 / ρ, 1e6)
    χρ = (ρ / Pt) * Ptρ
    χT = (T / Pt) * PtT
    δ  = χT / max(χρ, 1e-6)
    cP = cv + (Pt / (ρ * T)) * χT * δ
    ∇ad = clamp((Pt * δ) / max(ρ * T * cP, 1e-30), 0.05, 0.45)
    β = clamp(p0 / Pt, 1e-6, 1.0)
    return (; ρ, Pg = p0, Pr, Pt, u = u0, cP, δ, ∇ad, β)
end

# ---------------------------------------------------------------------------------------------------------------------
# 2. Composition helpers, Nuclear Burning (19-species NucNet), and Kippenhahn-Weigert MLT
# ---------------------------------------------------------------------------------------------------------------------
const IH = NN.IDX[:p]; const ID = NN.IDX[:d]
const HEL = [NN.IDX[s] for s in (:he3, :he4, :li7, :be7, :b8)]
const METALS = [i for i in 1:NN.NS if !(i == IH || i == ID || i in HEL)]

function xyz_of(Y)
    X  = NN.AA[IH] * Y[IH] + NN.AA[ID] * Y[ID]
    Yh = 0.0; for i in HEL; Yh += NN.AA[i] * Y[i]; end
    Z  = 0.0; for i in METALS; Z += NN.AA[i] * Y[i]; end
    Z < 1e-18 && (Z = 0.0)
    return X, Yh, Z
end

comp_of(Y) = (t = xyz_of(Y); SE.composition(t[1], t[2], t[3], SE.ZMIX))

const IHE3 = NN.IDX[:he3]; const ILI7 = NN.IDX[:li7]

function eps_nuc(ρ::Float64, T::Float64, Y, dt::Float64 = 0.0)
    T < 3.0e5 && return 0.0
    pl = NN.Plasma(ρ, Y)
    k = NN.rates_k(T, ρ, pl, NN.ye_of(Y); screen = :chugunov)
    e = 0.0
    @inbounds for (ir, r) in enumerate(NN.NET)
        ev = NN.event(r, k[ir], Y)
        # Step-averaged fuel depletion cap for fast light-element burning (D, ³He, ⁷Li) when dt >> τ_burn:
        if dt > 0.0 && ev > 0.0
            for isp in (ID, IHE3, ILI7)
                if Y[isp] > 0.0 && (isp in r.rin)
                    rate = ev / max(Y[isp], 1e-30)
                    x = rate * dt
                    if x > 1e-4
                        ev *= (-expm1(-min(x, 700.0))) / x
                    end
                    break
                end
            end
        end
        e += (NN.NA * NN.MeV * (r.Q - r.Enu)) * ev
    end
    return e
end

"Exact Kippenhahn-Weigert MLT temperature gradient ∇ and convective velocity v_conv [cm/s]."
function mlt_grad(P::Float64, T::Float64, ρ::Float64, κ::Float64, cP::Float64, δ::Float64,
                  ∇ad::Float64, L::Float64, r::Float64, m::Float64; α_mlt::Float64 = 1.5)
    ∇rad = (3.0 * κ * L * P) / (16π * arad * c_cgs * G * max(m, 1e20) * T^4)
    if ∇rad <= ∇ad || L <= 0.0
        return ∇rad, ∇rad, 0.0, false
    end
    g = G * m / max(r^2, 1e10)
    HP = P / max(ρ * g, 1e-30)
    ℓ = α_mlt * HP
    U = (3.0 * arad * c_cgs * T^3 / max(cP * ρ^2 * κ * ℓ^2, 1e-30)) * sqrt(8.0 * HP / max(g * δ, 1e-30))
    W = (∇rad - ∇ad) / max(U^2, 1e-60)
    z = W < 1.0 ? W / 2.0 : cbrt((8.0 / 9.0) * W)
    for _ in 1:15
        f  = 1.125 * z^3 + z^2 + 2.0 * z - W
        df = 3.375 * z^2 + 2.0 * z + 2.0
        dz = -f / df
        z = max(z + dz, 0.0)
        abs(dz) <= 1e-12 * max(z, 1.0) && break
    end
    y = z * U
    ∇ = clamp(∇ad + y * (y + 2.0 * U), ∇ad, ∇rad)
    v_conv = sqrt(g * δ / (8.0 * HP)) * ℓ * y
    return ∇, ∇rad, v_conv, true
end

# ---------------------------------------------------------------------------------------------------------------------
# 3. OP03 / Stahler (1986) Shock Outer Boundary Condition (and Cold-Disk Photospheric BC)
# ---------------------------------------------------------------------------------------------------------------------
"Integrate exact 16-substep hydrostatic + MLT atmosphere from τ=2/3 (P_ph, T_ph) down to P_N = P_ph + g*(0.5 dm_N)/(4π R^2)."
function integrate_atm_to_PN(M::Float64, R::Float64, L_int::Float64, P_ph::Float64, T_ph::Float64,
                             dm_half::Float64, cp, X::Float64, α_mlt::Float64)
    g = G * M / R^2
    P_N = P_ph + g * dm_half / (4π * R^2)
    P_N <= P_ph * (1.0 + 1e-9) && return P_N, T_ph
    n_sub = 16
    lP0 = log(P_ph); lP1 = log(P_N)
    dlP = (lP1 - lP0) / n_sub
    lT = log(T_ph)
    ρg = 1e-9
    for j in 1:n_sub
        lPa = lP0 + (j - 1) * dlP
        Pa = exp(lPa); Ta = exp(lT)
        eoa = eos_from_PT(Pa, Ta, cp; ρ_guess = ρg); ρg = eoa.ρ
        κa = SO.kappa_R(eoa.ρ, Ta, X)
        ∇a, _, _, _ = mlt_grad(Pa, Ta, eoa.ρ, κa, eoa.cP, eoa.δ, eoa.∇ad, L_int, R, M; α_mlt)
        # Midpoint RK2 step in ln P (cap ∇ <= 0.45 so Prad never exceeds P in super-Eddington layers)
        ∇a = clamp(∇a, 0.0, 0.45)
        Pb = exp(lPa + 0.5 * dlP); Tb = exp(lT + 0.5 * dlP * ∇a)
        eob = eos_from_PT(Pb, Tb, cp; ρ_guess = ρg); ρg = eob.ρ
        κb = SO.kappa_R(eob.ρ, Tb, X)
        ∇b, _, _, _ = mlt_grad(Pb, Tb, eob.ρ, κb, eob.cP, eob.δ, eob.∇ad, L_int, R, M; α_mlt)
        ∇b = clamp(∇b, 0.0, 0.45)
        lT += dlP * ∇b
    end
    return P_N, exp(lT)
end

function outer_bc(M::Float64, R::Float64, L_int::Float64, Mdot::Float64, X::Float64; bcmode::Symbol = :shock)
    vff = sqrt(2.0 * G * M / R)
    Pram = Mdot * vff / (4π * R^2)
    Lacc = G * M * Mdot / R
    g = G * M / R^2
    cp0 = SE.composition(X, 1.0 - X, 0.0)
    if bcmode === :shock
        Ltot = max(L_int, 1e-2 * Lsun) + Lacc
        Teff_sh = (Ltot / (4π * R^2 * σSB))^0.25
        if bcmode === :shock && Mdot > 0.0
            n_env = 40
            r_e = exp.(range(log(R), log(30.0 * R); length = n_env + 1))
            tau_sph = 0.0; tau_lin = 0.0
            Rph = R
            tau_arr = zeros(n_env + 1)
            for i in n_env:-1:1
                rm = sqrt(r_e[i] * r_e[i+1]); dr = r_e[i+1] - r_e[i]
                ρm = Mdot / (4π * rm^2 * sqrt(2.0 * G * M / rm))
                Tm = max(Teff_sh * sqrt(R / rm) * 0.84, 300.0)
                κm = SO.kappa_R(ρm, Tm, X)
                tau_sph += κm * ρm * (R / rm)^2 * dr
                tau_lin += κm * ρm * dr
                tau_arr[i] = tau_lin
            end
            if tau_lin >= 2.0 / 3.0
                for i in n_env:-1:1
                    if tau_arr[i] >= 2.0 / 3.0
                        f = (2.0 / 3.0 - tau_arr[i+1]) / max(tau_arr[i] - tau_arr[i+1], 1e-30)
                        Rph = r_e[i+1] + clamp(f, 0.0, 1.0) * (r_e[i] - r_e[i+1])
                        break
                    end
                end
            end
            Tsurf = Teff_sh * (0.5 + 0.75 * min(tau_sph, 1.5))^0.25
        else
            tau_sph = 0.0; tau_lin = 0.0; Rph = R; Tsurf = Teff_sh
        end
        Pr_s = arad * Tsurf^4 / 3.0
        rhs = (2.0 / 3.0) * g
        lo_p = max(1.01 * Pr_s, Pram + 1.01 * Pr_s, 1e-3); hi_p = max(1e12, 1e6 * lo_p)
        f_sh(Pt_) = (Pt_ - Pram - 0.5 * Pr_s) * max(SO.kappa_R(eos_from_PT(Pt_, Tsurf, cp0).ρ, Tsurf, X), 1e-5) - rhs
        if f_sh(lo_p) >= 0.0
            Pext = lo_p
        else
            for _ in 1:40
                pm = sqrt(lo_p * hi_p)
                f_sh(pm) > 0.0 ? (hi_p = pm) : (lo_p = pm)
                hi_p / lo_p - 1.0 < 1e-6 && break
            end
            Pext = sqrt(lo_p * hi_p)
        end
        Teff_obs = (Ltot / (4π * Rph^2 * σSB))^0.25
        return (; Pext, Tsurf, Rph, Teff = Teff_obs, Ltot, Lacc, tau_env = tau_lin, Pr_s)
    else
        # Cold-disk / photospheric accretion (HOY12 §2.1): L = 4π R_*^2 σ T_eff^4
        # Exact grey Eddington photosphere (Kippenhahn & Weigert eq. 11.18):
        #   P_gas(τ=2/3) = (2/3)(g/κ_R)(1 - Γ) => (P_tot - 0.5 P_rad) * κ_R(ρ, T_eff) = (2/3) g
        Ltot = max(L_int, 1e-2 * Lsun)
        Teff = (Ltot / (4π * R^2 * σSB))^0.25
        Tsurf = Teff
        Pr_s = arad * Tsurf^4 / 3.0
        rhs = (2.0 / 3.0) * g
        lo_p = max(1.01 * Pr_s, 1e-3); hi_p = max(1e12, 1e6 * Pr_s)
        f_ph(Pt_) = (Pt_ - 0.5 * Pr_s) * max(SO.kappa_R(eos_from_PT(Pt_, Tsurf, cp0).ρ, Tsurf, X), 1e-5) - rhs
        if f_ph(lo_p) >= 0.0
            Pext = lo_p
        else
            for _ in 1:40
                pm = sqrt(lo_p * hi_p)
                f_ph(pm) > 0.0 ? (hi_p = pm) : (lo_p = pm)
                hi_p / lo_p - 1.0 < 1e-6 && break
            end
            Pext = sqrt(lo_p * hi_p)
        end
        return (; Pext, Tsurf, Rph = R, Teff, Ltot, Lacc, tau_env = 0.0, Pr_s)
    end
end

# ---------------------------------------------------------------------------------------------------------------------
# 4. Exact 4-Equation Henyey Block-Tridiagonal Solver (Sugimoto Homologous Accretion!)
# ---------------------------------------------------------------------------------------------------------------------
"""
Evaluate the 4N Henyey residual vector F(x) for state x = [ln r_2, ln P_1, ln T_1, L_2/L_sc, ...]:
  - Eq 1 (k=1..N): Mass conservation  ln[(r_{k+1}^3 - r_k^3)/3] - ln[dm_k / (4π ρ_k)] = 0
  - Eq 2 (k=1..N-1): HSE             (P_{k+1} - P_k)/P̄_k + G m_{k+1} dm̄_k / (4π r_{k+1}^4 P̄_k) = 0
         (k=N):      Outer P BC       ln P_N - ln[P_ext + G M_* (dm_N/2) / (4π R_*^4)] = 0
  - Eq 3 (k=1..N-1): Energy transp.  (ln T_{k+1} - ln T_k) - ∇_{k+1/2} (ln P_{k+1} - ln P_k) = 0
         (k=N):      Outer T BC       ln T_N - ln T_surf(M, R_*, L_{N+1}) = 0
  - Eq 4 (k=1..N):   Energy eq.      [(L_{k+1} - L_k)/dm_k - ε_nuc,k - ε_g,k] / (L_sc / M_*) = 0
"""
function henyey_residual!(F::Vector{Float64}, x::Vector{Float64}, st)
    N = st.N; M = st.M; Mdot = st.Mdot; dt = st.dt; Lsc = st.Lsc
    q = st.q; dm = st.dm; qm = st.qm
    lnP_old = st.lnP_old; lnT_old = st.lnT_old
    Y = st.Y; cps = st.cps; ρ_hint = st.ρ_hint; α_mlt = st.α_mlt; bcmode = st.bcmode

    # Unpack state variables
    r   = st.r_buf;   P   = st.P_buf;   T   = st.T_buf;   L   = st.L_buf
    ρ   = st.ρ_buf;   cP  = st.cP_buf;  δ   = st.δ_buf;   ∇ad = st.nad_buf
    κ   = st.κ_buf;   εn  = st.εn_buf;  vc  = st.vc_buf;  isconv = st.conv_buf

    r[1] = 0.0; L[1] = 0.0
    @inbounds for k in 1:N
        i0 = 4 * (k - 1)
        r[k+1] = exp(clamp(x[i0 + 1], 10.0, 42.0))
        P[k]   = exp(clamp(x[i0 + 2], -15.0, 50.0))
        T[k]   = exp(clamp(x[i0 + 3], 5.5, 22.0))
        # For photospheric BC, store outer luminosity as ln(L_surf / Lsc) so L_surf > 0 strictly!
        L[k+1] = (k == N && bcmode === :phot) ? exp(clamp(x[i0 + 4], -30.0, 15.0)) * Lsc : x[i0 + 4] * Lsc
        eo = eos_from_PT(P[k], T[k], cps[k]; ρ_guess = ρ_hint[k])
        ρ[k] = eo.ρ; cP[k] = eo.cP; δ[k] = eo.δ; ∇ad[k] = eo.∇ad
        Xk = NN.AA[IH] * Y[IH, k] + NN.AA[ID] * Y[ID, k]
        κ[k] = SO.kappa_R(ρ[k], T[k], Xk)
        εn[k] = st.burn_on ? eps_nuc(ρ[k], T[k], view(Y, :, k), dt) : 0.0
    end

    # Outer BC at k = N (R_* = r[N+1], L_int = L[N+1])
    Xsurf = NN.AA[IH] * Y[IH, N] + NN.AA[ID] * Y[ID, N]
    bc = outer_bc(M, r[N+1], L[N+1], Mdot, Xsurf; bcmode)

    @inbounds for k in 1:N
        i0 = 4 * (k - 1)
        # Eq 1: Mass conservation
        dV = (r[k+1]^3 - r[k]^3) / 3.0
        F[i0 + 1] = log(max(dV, 1e10)) - log(dm[k] / (4π * ρ[k]))

        # Eq 2 & 3: HSE and Energy Transport
        if k < N
            rk = r[k+1]; mk = q[k+1] * M; dmb = 0.5 * (dm[k] + dm[k+1])
            Pb = sqrt(P[k] * P[k+1])
            F[i0 + 2] = (P[k+1] - P[k]) / Pb + G * mk * dmb / (4π * rk^4 * Pb)

            # Interface averages for exact Kippenhahn-Weigert MLT gradient (continuous weight 1/(1+max(-∇rad,0)) for L < 0!)
            Tb   = sqrt(T[k] * T[k+1])
            ρb   = sqrt(ρ[k] * ρ[k+1])
            κb   = sqrt(κ[k] * κ[k+1])
            cPb  = 0.5 * (cP[k] + cP[k+1])
            δb   = 0.5 * (δ[k] + δ[k+1])
            nadb = 0.5 * (∇ad[k] + ∇ad[k+1])
            ∇, ∇rad, vck, cvk = mlt_grad(Pb, Tb, ρb, κb, cPb, δb, nadb, L[k+1], rk, mk; α_mlt)
            vc[k] = vck; isconv[k] = cvk
            w_neg = 1.0 / (1.0 + max(-∇rad, 0.0))
            F[i0 + 3] = w_neg * ((log(T[k+1]) - log(T[k])) - ∇ * (log(P[k+1]) - log(P[k])))
        else
            # Outer mechanical & thermal BC at center of outermost shell k = N via exact 16-substep HSE+MLT atmosphere!
            P_N_tgt, T_N_tgt = integrate_atm_to_PN(M, r[N+1], max(L[N+1], 1e-2 * Lsun),
                                                   bc.Pext, bc.Tsurf, 0.5 * dm[N], cps[N], Xsurf, α_mlt)
            F[i0 + 2] = log(P[N]) - log(max(P_N_tgt, 1e-10))
            F[i0 + 3] = log(T[N]) - log(max(T_N_tgt, 300.0))
            vc[N] = 0.0; isconv[N] = false
        end

        # Eq 4: Energy conservation with Sugimoto (1970) / OP03 eq. (4) homologous accretion!
        # ε_g = -T (∂s/∂t)_q + (Mdot/M) q T (∂s/∂q)_t
        # where T Δs = cP T Δ(ln T) - (δ P / ρ) Δ(ln P) = cP T [Δ(ln T) - ∇ad Δ(ln P)]
        if dt > 0.0
            dlnT_dt = (log(T[k]) - lnT_old[k]) / dt
            dlnP_dt = (log(P[k]) - lnP_old[k]) / dt
            εg_loc  = -cP[k] * T[k] * (dlnT_dt - ∇ad[k] * dlnP_dt)
            # Upwind homologous advection from shell k+1 -> k (for k = N, fresh gas enters at surface entropy -> 0 gradient)
            if k < N && Mdot > 0.0
                dqk = qm[k+1] - qm[k]
                dlnT_dq = (log(T[k+1]) - log(T[k])) / dqk
                dlnP_dq = (log(P[k+1]) - log(P[k])) / dqk
                εg_nh = (Mdot / M) * qm[k] * cP[k] * T[k] * (dlnT_dq - ∇ad[k] * dlnP_dq)
            else
                εg_nh = 0.0
            end
            εg = εg_loc + εg_nh
        else
            εg = 0.0
        end
        F[i0 + 4] = ((L[k+1] - L[k]) - dm[k] * (εn[k] + εg)) / Lsc
    end
    return bc
end

"Exact 8-diagonal colored finite-difference Jacobian + sparse LU Newton solver for the 4N Henyey system."
function solve_henyey!(x::Vector{Float64}, st; maxit::Int = 30, tol::Float64 = 1e-6)
    nvar = length(x)
    F0 = zeros(nvar); Fp = zeros(nvar); xp = copy(x)
    # Each equation in block k depends ONLY on blocks k-1, k, k+1 (and outer BC k=N on R_*, L_*):
    # Bandwidth <= 12 -> 12 colored perturbation sweeps evaluate the ENTIRE sparse Jacobian in 12 residual calls!
    I_idx = Int[]; J_idx = Int[]; V_idx = Float64[]
    bc = nothing
    for it in 1:maxit
        bc = henyey_residual!(F0, x, st)
        err = maximum(abs, F0)
        if get(ENV, "DBG", "0") == "1" && (it <= 5 || it % 5 == 0)
            im = argmax(abs.(F0))
            @printf("    it=%2d err=%.3e (blk=%d eq=%d) R*=%.2f Tc=%.3e Lint=%.3e\n",
                    it, err, (im-1)÷4+1, (im-1)%4+1, st.r_buf[end]/Rsun, st.T_buf[1], st.L_buf[end]/Lsun)
        end
        if err < tol
            # Update ρ_hint
            st.ρ_hint .= st.ρ_buf
            return true, it, err, bc
        end
        empty!(I_idx); empty!(J_idx); empty!(V_idx)
        # 12-color compression: perturb variable j = c, c+12, c+24, ... simultaneously!
        # Note: outer BC depends on r[N+1] (var 4N-3) and L[N+1] (var 4N), which are in block N, so k=N couples ONLY to block N (and k=N-1)!
        # Block k eqs couple to blocks k-1, k, k+1 => max column separation in any row is <= 11 => 12 colors exact!
        for col_c in 1:12
            copyto!(xp, x)
            jlist = col_c:12:nvar
            hvec = [j % 4 == 0 ? 1e-7 * max(abs(x[j]), 1e-4) : 1e-6 * max(abs(x[j]), 1.0) for j in jlist]
            for (idx, j) in enumerate(jlist); xp[j] += hvec[idx]; end
            henyey_residual!(Fp, xp, st)
            for (idx, j) in enumerate(jlist)
                hj = hvec[idx]
                blk = (j - 1) ÷ 4 + 1
                ilo = max(1, 4 * (blk - 2) + 1)
                ihi = min(nvar, 4 * (blk + 1))
                @inbounds for i in ilo:ihi
                    v = (Fp[i] - F0[i]) / hj
                    if v != 0.0
                        push!(I_idx, i); push!(J_idx, j); push!(V_idx, v)
                    end
                end
            end
        end
        Jmat = sparse(I_idx, J_idx, V_idx, nvar, nvar)
        dx = try
            -(Jmat \ F0)
        catch
            return false, it, err, bc
        end
        any(!isfinite, dx) && return false, it, err, bc
        # Damped step cap on log variables (|Δln r, ln P, ln T| <= 0.5, |Δln L_surf| <= 0.8) + backtracking line search
        max_log = 0.0
        @inbounds for k in 1:st.N
            i0 = 4 * (k - 1)
            max_log = max(max_log, abs(dx[i0+1]), abs(dx[i0+2]), abs(dx[i0+3]))
        end
        if st.bcmode === :phot
            max_log = max(max_log, abs(dx[4*st.N]) / 1.6)
        end
        α = min(1.0, 0.5 / max(max_log, 1e-30))
        norm0 = norm(F0)
        accepted = false
        for _ in 1:12
            @. xp = x + α * dx
            henyey_residual!(Fp, xp, st)
            if isfinite(norm(Fp)) && (norm(Fp) < (1.0 - 1e-5 * α) * norm0 || norm(Fp) < tol)
                copyto!(x, xp); st.ρ_hint .= st.ρ_buf
                accepted = true; break
            end
            α *= 0.5
        end
        if !accepted
            # Take a gentle damped step if close, else bail
            @. x = x + 0.25 * min(1.0, 0.2 / max(max_log, 1e-30)) * dx
            st.ρ_hint .= st.ρ_buf
        end
    end
    bc = henyey_residual!(F0, x, st)
    return maximum(abs, F0) < 10 * tol, maxit, maximum(abs, F0), bc
end

# ---------------------------------------------------------------------------------------------------------------------
# 5. Burn + Homologous Composition Advection + Convective Mixing
# ---------------------------------------------------------------------------------------------------------------------
function evolve_composition!(st, dt::Float64, Y_prim::Vector{Float64})
    N = st.N; Y = st.Y; M = st.M; Mdot = st.Mdot; dm = st.dm; q = st.q
    # 1. Nuclear burning per zone (NucNet.burn where T > 5e5 K)
    if st.burn_on
        @inbounds for k in 1:N
            st.T_buf[k] > 5.0e5 || continue
            Ynew, _, _ = NN.burn(view(Y, :, k), st.ρ_buf[k], st.T_buf[k], dt; screen = :chugunov)
            Y[:, k] .= Ynew
        end
    end
    # 2. Conservative upwind homologous accretion advection of abundances Y_i(q):
    #    d(dm_k Y_k)/dt = Mdot * (q_{k+1} Y_{k+1} - q_k Y_k)  (with Y_{N+1} = Y_prim)
    #    Implicit backward Euler (unconditionally positive and conservative!):
    if Mdot > 0.0 && dt > 0.0
        # Sweep inward from k = N down to 1 (exact implicit upwind!)
        Y_up = copy(Y_prim)
        for k in N:-1:1
            fin  = (Mdot * dt / dm[k]) * q[k+1]
            fout = (Mdot * dt / dm[k]) * q[k]
            # Y_new * (1 + fout) = Y_old + fin * Y_up? Wait: dm_k also grew by (Mdot*dt/M)*dm_k = (fin - fout)*dm_k!
            # Exact finite-volume conservation on expanding mass grid:
            #   dm_k(t+dt) Y_k^new - dm_k(t) Y_k^old = Mdot*dt * (q_{k+1} Y_{k+1}^new - q_k Y_k^new)
            #   Since dm_k(t+dt) = dm_k(t) * (1 + dM/M) and dm_k(t)*(dM/M) = Mdot*dt*(q_{k+1} - q_k):
            #   => (1 + fin) Y_k^new = Y_k^old + fin * Y_{k+1}^new  (EXACTLY conserves Σ dm_k Y_k + dM Y_prim!)
            denom = 1.0 + fin
            @inbounds for s in 1:NN.NS
                Y[s, k] = (Y[s, k] + fin * Y_up[s]) / denom
                Y_up[s] = Y[s, k]
            end
        end
    end
    # 3. Convective mixing across contiguous convective zones (instantaneous mixing where v_conv > 0 & τ_mix < dt)
    k = 1
    while k <= N
        if k < N && st.conv_buf[k]
            k0 = k
            while k < N && st.conv_buf[k]; k += 1; end
            k1 = k  # zones k0..k1 are convectively coupled!
            mtot = sum(dm[k0:k1])
            for s in 1:NN.NS
                ybar = sum(dm[j] * Y[s, j] for j in k0:k1) / mtot
                for j in k0:k1; Y[s, j] = ybar; end
            end
        end
        k += 1
    end
    # Refresh zone compositions
    for k in 1:N; st.cps[k] = comp_of(view(Y, :, k)); end
end

# ---------------------------------------------------------------------------------------------------------------------
# 6. Main Driver
# ---------------------------------------------------------------------------------------------------------------------
function main()
    envf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
    envi(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d

    Mdot_yr = envf("MDOT", 4.41e-3)
    Mdot    = Mdot_yr * Msun / yr
    M0_msun = envf("M0", 0.5)
    R0_rsun = envf("R0", 25.0)
    Mstop   = envf("MSTOP", 100.0)
    N       = envi("N", 120)
    nenv    = envi("NENV", 45)
    α_mlt   = envf("ALPHAMLT", 1.5)
    bcmode  = Symbol(get(ENV, "BCMODE", "shock"))
    dM_frac = envf("DMFRAC", 0.015)
    maxsteps = envi("MAXSTEPS", 1500)
    outdir  = get(ENV, "VAL_OUT", "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs")
    tag     = get(ENV, "TAG", @sprintf("op03_hen_%.2e", Mdot_yr))
    mkpath(outdir)
    csvio = open(joinpath(outdir, "$(tag).csv"), "w")
    logio = open(joinpath(outdir, "$(tag).log"), "w")
    pr(f, args...) = (s = Printf.format(Printf.Format(f), args...); print(s); print(logio, s); flush(stdout); flush(logio))

    # Fixed relative mass coordinate q = m / M_*(t) down to atmosphere base qfloor = 1e-5
    qsplit = 0.15; qfloor = envf("MFLOOR", 1e-5); nin = N - nenv
    w = [min(1.10^(j - 1), 25.0) for j in 1:nin]
    qint = (1.0 - qsplit) .* [0.0; cumsum(w)] ./ sum(w)
    qe = 10 .^ range(log10(qsplit), log10(qfloor); length = nenv + 1)[2:end]
    q = vcat(qint, 1.0 .- qe); q[end] = 1.0
    qm = [0.5 * (q[k] + q[k+1]) for k in 1:N]

    # Primordial abundances (X=0.75, Y=0.25, D/H=2.5e-5, 3He/H=1e-5)
    DH = envf("DH", 2.5e-5)
    Y_prim = zeros(NN.NS); Yp = 0.75 / (1.0 + 2DH); Y3 = 1e-5 * Yp
    Y_prim[IH] = Yp; Y_prim[ID] = DH * Yp; Y_prim[NN.IDX[:he3]] = Y3; Y_prim[NN.IDX[:he4]] = (0.25 - 3Y3) / 4.0
    Y = repeat(Y_prim, 1, N)
    cps = [comp_of(view(Y, :, k)) for k in 1:N]

    # Build exact Lane-Emden polytrope seed star of exact (M0, R0) on mass grid q:
    M0 = M0_msun * Msun; R0 = R0_rsun * Rsun
    dm0 = [M0 * (q[k+1] - q[k]) for k in 1:N]
    npoly = envf("NPOLY", 1.5)
    ξ = [1e-8]; θp = [1.0]; μp = [0.0]; h_le = 1e-3; y_le = (1.0 - 1e-16, 0.0); x_le = 1e-8
    while y_le[1] > 0 && x_le < 25.0
        f_le(xx, th, ph) = (ph, -max(th, 0.0)^npoly - 2ph / xx)
        k1 = f_le(x_le, y_le...); k2 = f_le(x_le + h_le/2, y_le[1] + h_le/2*k1[1], y_le[2] + h_le/2*k1[2])
        k3 = f_le(x_le + h_le/2, y_le[1] + h_le/2*k2[1], y_le[2] + h_le/2*k2[2]); k4 = f_le(x_le + h_le, y_le[1] + h_le*k3[1], y_le[2] + h_le*k3[2])
        yn = (y_le[1] + h_le/6*(k1[1]+2k2[1]+2k3[1]+k4[1]), y_le[2] + h_le/6*(k1[2]+2k2[2]+2k3[2]+k4[2]))
        x_le += h_le
        if yn[1] <= 0
            xz = x_le - h_le + h_le * y_le[1] / (y_le[1] - yn[1])
            push!(ξ, xz); push!(θp, 0.0); push!(μp, -xz^2 * (y_le[2] + (yn[2] - y_le[2]) * y_le[1] / (y_le[1] - yn[1])))
            break
        end
        y_le = yn; push!(ξ, x_le); push!(θp, y_le[1]); push!(μp, -x_le^2 * y_le[2])
    end
    ξ1 = ξ[end]; μ1 = μp[end]; al = R0 / ξ1
    ξof(qq) = (t = qq * μ1; t >= μ1 && return ξ1; i = clamp(searchsortedlast(μp, t), 1, length(μp)-1);
               ξ[i] + (t - μp[i]) / max(μp[i+1] - μp[i], 1e-30) * (ξ[i+1] - ξ[i]))
    r0 = [al * ξof(qq) for qq in q]; r0[1] = 0.0; r0[end] = R0
    let qs = 1e-4, δs = R0 - al * ξof(1.0 - qs), pexp = 1.0 / (npoly + 1.0)
        for k in 2:N
            qext = 1.0 - q[k]
            qext < qs && (r0[k] = R0 - δs * (qext / qs)^pexp)
        end
    end
    for k in 2:N; r0[k] <= r0[k-1] && (r0[k] = 0.5 * (r0[k-1] + r0[k+1])); end
    ρ0 = [dm0[k] / (4π * (r0[k+1]^3 - r0[k]^3) / 3.0) for k in 1:N]
    # Exact discrete HSE pressures P0 inward from photospheric BC + atmosphere integration at R0 (for clean pre-relaxation):
    L_init = 4π * R0^2 * σSB * 5000.0^4
    bc0 = outer_bc(M0, R0, L_init, 0.0, 0.75; bcmode = :phot)
    P_N0, T_N0 = integrate_atm_to_PN(M0, R0, max(L_init, 1e-2 * Lsun), bc0.Pext, bc0.Tsurf, 0.5 * dm0[N], cps[N], 0.75, α_mlt)
    P0 = zeros(N)
    P0[N] = P_N0
    for k in N-1:-1:1
        P0[k] = P0[k+1] + G * (q[k+1] * M0) * 0.5 * (dm0[k] + dm0[k+1]) / (4π * r0[k+1]^4)
    end
    # Exact T0(ρ0, P0) via bisection on FreeEnergyEOS, smoothly anchored to T_N0 at surface:
    T0 = zeros(N)
    for k in 1:N
        lo_t, hi_t = 1000.0, 5e8
        for _ in 1:50
            tm = sqrt(lo_t * hi_t)
            pm = _pe_gas(ρ0[k], tm, cps[k])[1] + arad * tm^4 / 3.0
            pm > P0[k] ? (hi_t = tm) : (lo_t = tm)
            hi_t / lo_t - 1.0 < 1e-6 && break
        end
        T0[k] = max(sqrt(lo_t * hi_t), bc0.Tsurf)
    end
    T0[N] = T_N0
    Lsc0 = 4π * c_cgs * G * M0 / (0.2 * 1.75)   # L_Edd(M0)
    L0 = [L_init * q[k] for k in 1:N+1]

    st = (;
        N, q, qm, dm = copy(dm0),
        M = M0, Mdot, dt = 0.0, Lsc = Lsc0,
        lnP_old = log.(P0), lnT_old = log.(T0),
        Y, cps, ρ_hint = copy(ρ0), α_mlt, bcmode, burn_on = false,
        r_buf = copy(r0), P_buf = copy(P0), T_buf = copy(T0), L_buf = copy(L0),
        ρ_buf = copy(ρ0), cP_buf = zeros(N), δ_buf = zeros(N), nad_buf = zeros(N),
        κ_buf = zeros(N), εn_buf = zeros(N), vc_buf = zeros(N), conv_buf = falses(N)
    )
    # Use a mutable wrapper for scalar fields
    M_ref = Ref(M0); Mdot_ref = Ref(Mdot); dt_ref = Ref(0.0); Lsc_ref = Ref(Lsc0); burn_ref = Ref(false); bc_ref = Ref(:phot)
    make_st() = (; N, q, qm, dm = st.dm, M = M_ref[], Mdot = Mdot_ref[], dt = dt_ref[], Lsc = Lsc_ref[],
                   lnP_old = st.lnP_old, lnT_old = st.lnT_old, Y = st.Y, cps = st.cps,
                   ρ_hint = st.ρ_hint, α_mlt, bcmode = bc_ref[], burn_on = burn_ref[],
                   r_buf = st.r_buf, P_buf = st.P_buf, T_buf = st.T_buf, L_buf = st.L_buf,
                   ρ_buf = st.ρ_buf, cP_buf = st.cP_buf, δ_buf = st.δ_buf, nad_buf = st.nad_buf,
                   κ_buf = st.κ_buf, εn_buf = st.εn_buf, vc_buf = st.vc_buf, conv_buf = st.conv_buf)

    x = zeros(4 * N)
    for k in 1:N
        i0 = 4 * (k - 1)
        x[i0 + 1] = log(r0[k+1]); x[i0 + 2] = log(P0[k]); x[i0 + 3] = log(T0[k])
        x[i0 + 4] = (k == N) ? log(max(L0[k+1] / Lsc0, 1e-12)) : L0[k+1] / Lsc0
    end
    # Quick diagnostic helper if DBG=1
    if get(ENV, "DBG", "0") == "1"
        Ftmp = zeros(4 * N); henyey_residual!(Ftmp, x, make_st())
        im = argmax(abs.(Ftmp))
        @printf("  DBG init max|F|=%.3e at var=%d (blk=%d, eq=%d) Tc=%.3e P1=%.3e Tsurf=%.1f\n",
                abs(Ftmp[im]), im, (im-1)÷4+1, (im-1)%4+1, T0[1], P0[1], T0[N])
    end

    pr("VAL_OP03_HENYEY [%s]: Mdot=%.3e Msun/yr, M0=%.2f Msun, R0=%.1f Rsun, Mstop=%.1f Msun, N=%d, BC=%s, FREE_EOS + StellarOpacity + NucNet + MLT\n",
       tag, Mdot_yr, M0_msun, R0_rsun, Mstop, N, string(bcmode))
    println(csvio, "step,t_yr,dt_yr,M_star,R_star,R_ph,Teff,L_int,L_acc,L_tot,L_nuc,L_Edd,tau_env,Tc,rho_c,Xc,Dc_over_D0,Mcc_over_M")
    pr("%5s %10s %8s %8s %8s %7s %9s %9s %9s %7s %8s %7s %7s\n",
       "step", "t [yr]", "M*", "R* [Rs]", "Rph [Rs]", "Teff", "Lint[Ls]", "Ltot[Ls]", "Lnuc[Ls]", "L/LEdd", "Tc [1e7]", "Dc/D0", "Mcc/M")

    # Pre-relax the Lane-Emden seed star (3 gentle Mdot=0 photospheric steps, burn=off) onto the exact FreeEnergyEOS+MLT adiabat
    # so EVERY accretion rate (phot OR shock) starts from the EXACT same converged hydrostatic seed star!
    bc_ref[] = :phot; Mdot_ref[] = 0.0; burn_ref[] = false
    for rstep in 1:3
        st.lnP_old .= [x[4*(k-1)+2] for k in 1:N]
        st.lnT_old .= [x[4*(k-1)+3] for k in 1:N]
        dt_ref[] = 5.0 * yr
        ok_r, _, err_r, _ = solve_henyey!(x, make_st(); maxit = 35, tol = 1e-6)
        !ok_r && @warn "Pre-relaxation step $rstep did not fully converge (err=$err_r)"
    end
    bc_ref[] = bcmode
    if bcmode === :shock
        x[4*N] = exp(x[4*N])   # convert ln(L_surf/Lsc) -> linear L_surf/Lsc for shock BC
    end
    Mdot_ref[] = Mdot

    t_yr = 0.0
    dt_cur = 0.002 * M0 / Mdot
    D0 = Y_prim[ID]
    burn_ref[] = true
    t0w = time()

    for step in 1:maxsteps
        M_old = M_ref[]
        x_save = copy(x); ρ_save = copy(st.ρ_hint); Y_save = copy(st.Y)
        lnP_prev = [x[4*(k-1)+2] for k in 1:N]
        lnT_prev = [x[4*(k-1)+3] for k in 1:N]
        st.lnP_old .= lnP_prev; st.lnT_old .= lnT_prev

        ok = false; bc = nothing; it_used = 0; err_used = NaN
        dt_try = min(dt_cur, dM_frac * M_old / Mdot)
        for att in 0:8
            dt_now = dt_try * 0.45^att
            M_new = M_old + Mdot * dt_now
            M_ref[] = M_new; dt_ref[] = dt_now
            Lsc_ref[] = 4π * c_cgs * G * M_new / (0.2 * 1.75)
            for k in 1:N; st.dm[k] = M_new * (q[k+1] - q[k]); end
            # Rescale L unknowns to new Lsc
            copyto!(x, x_save); copyto!(st.ρ_hint, ρ_save)
            scale_L = (4π * c_cgs * G * M_old / (0.2 * 1.75)) / Lsc_ref[]
            for k in 1:N
                if k == N && bcmode === :phot
                    x[4*(k-1)+4] += log(scale_L)
                else
                    x[4*(k-1)+4] *= scale_L
                end
            end
            ok, it_used, err_used, bc = solve_henyey!(x, make_st(); maxit = 35, tol = 1e-5)
            if ok
                dt_cur = dt_now
                break
            end
        end
        if !ok
            pr("STOPPED (convergence failure) at step %d: M*=%.3f Msun, R*=%.1f Rsun, err=%.2e (%.1fs)\n",
               step, M_old / Msun, exp(x_save[4*(N-1)+1]) / Rsun, err_used, time() - t0w)
            break
        end

        # Evolve nuclear abundances + homologous accretion + convective mixing
        evolve_composition!(make_st(), dt_cur, Y_prim)
        t_yr += dt_cur / yr

        # Adaptive dt growth based on max change in (ln R_*, ln T_c, ln P_c, iter)
        dR = abs(x[4*(N-1)+1] - x_save[4*(N-1)+1])
        dTc = abs(x[3] - lnT_prev[1])
        dPc = abs(x[2] - lnP_prev[1])
        dmax = max(dR / 0.06, dTc / 0.05, dPc / 0.08, 1e-6)
        grow = it_used <= 8 ? 1.4 : (it_used <= 16 ? 1.15 : 0.85)
        dt_cur = min(dt_cur * clamp(0.9 / dmax, 0.5, grow), dM_frac * M_ref[] / Mdot)

        # Diagnostics
        R_s = st.r_buf[N+1]
        L_int = st.L_buf[N+1]
        L_nuc = sum(st.dm[k] * st.εn_buf[k] for k in 1:N)
        LEdd = Lsc_ref[]
        # Convective core mass fraction
        kc = 0
        while kc < N && st.conv_buf[kc + 1]; kc += 1; end
        Mcc = kc > 0 ? q[kc + 1] : 0.0
        M_msun = M_ref[] / Msun
        X_c = NN.AA[IH] * st.Y[IH, 1] + NN.AA[ID] * st.Y[ID, 1]

        if step == 1 || step % envi("PRINTEVERY", 15) == 0 || M_msun >= Mstop || R_s / Rsun > 4000.0
            pr("%5d %10.3e %8.3f %8.2f %8.2f %7.0f %9.3e %9.3e %9.3e %7.4f %8.4f %7.4f %7.4f\n",
               step, t_yr, M_msun, R_s / Rsun, bc.Rph / Rsun, bc.Teff, L_int / Lsun, bc.Ltot / Lsun,
               L_nuc / Lsun, bc.Ltot / LEdd, st.T_buf[1] / 1e7, st.Y[ID, 1] / D0, Mcc)
        end
        @printf(csvio, "%d,%.6e,%.4e,%.5f,%.4f,%.4f,%.2f,%.5e,%.5e,%.5e,%.5e,%.5e,%.4e,%.5e,%.5e,%.6f,%.5e,%.5f\n",
                step, t_yr, dt_cur / yr, M_msun, R_s / Rsun, bc.Rph / Rsun, bc.Teff, L_int / Lsun, bc.Lacc / Lsun,
                bc.Ltot / Lsun, L_nuc / Lsun, LEdd / Lsun, bc.tau_env, st.T_buf[1], st.ρ_buf[1], X_c, st.Y[ID, 1] / D0, Mcc)
        flush(csvio)

        if M_msun >= Mstop || R_s / Rsun > 4000.0
            pr("DONE at step %d: M*=%.2f Msun, R*=%.1f Rsun, Teff=%.0f K, L/LEdd=%.4f, t=%.3e yr (wall %.1f s)\n",
               step, M_msun, R_s / Rsun, bc.Teff, bc.Ltot / LEdd, t_yr, time() - t0w)
            break
        end
    end
    close(csvio); close(logio)
end

main()
