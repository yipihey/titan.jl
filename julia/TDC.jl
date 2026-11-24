"""
    TDC

Pure, number-type-generic functions of the sub-shock convection work (docs/POPIII_CONVECTION.md, D1/M2; the time-dependent
convection (TDC) of M3 goes in here too). Nothing here reads global state: the arguments may be Float64 or ForwardDiff Duals
(of any nesting), so the Newton Jacobian sees every dependence.
"""
module TDC

export gamma1_inv, gamma1, blob_sup, blob_f
export TDCPar, TDCGeom, TDCBack, AC_CGS, fb, tdc_ay, tdc_local, tdc_flux, tdc_pturb, tdc_dmix, tdc_mlt_w, tdc_mlt_zeta, kw_mlt_flux
export tdc_F!, tdc_rows!, tdc_exchange, tdc_gas_exchange, tdc_budget_check, tdc_facew

"""
    gamma1_inv(α, δ, ∇ad) = 1/Γ₁ = α − δ ∇_ad

with α = (∂ln ρ/∂ln P)_{T,x} and δ = −(∂ln ρ/∂ln T)_{P,x}. [derivation] A parcel moving adiabatically has dln T = ∇_ad dln P, so
dln ρ = α dln P − δ dln T = (α − δ ∇_ad) dln P ≡ dln P/Γ₁.
"""
@inline gamma1_inv(α, δ, nad) = α - δ * nad
@inline gamma1(α, δ, nad) = inv(gamma1_inv(α, δ, nad))

"""
    blob_f(ℓ, a, t_chem) = ℓ/(ℓ + t_chem a)

Fraction of equilibrium chemistry in a convective blob: t_turn/(t_turn + t_chem) with the turnover time ℓ/a. t_chem = 0 (equilibrium
zone) gives 1, t_chem = ∞ (frozen composition) gives 0. For t_chem = Inf the result is the Float64 zero (no Inf·Dual).
"""
@inline function blob_f(ℓ, a, t_chem)
    isinf(t_chem) && return zero(ℓ * a)
    return ℓ / (ℓ + t_chem * a)
end

"""
    (sup, δ_b) = blob_sup(ΔlnP, ΔlnT, Δlnμ, α_fr, δ_fr, φ_fr, Γ1_eq, Γ1_fr, δ_eq, f, ΔlnP_hs)

Superadiabaticity of the buoyancy of a blob across one interface, replacing ∇ − ∇_ad (design D1). Differences are taken outwards
(zone k−1 → k). Background, from the frozen-composition EOS of the layer and the composition difference,

    Δln ρ = α_fr ΔlnP − δ_fr ΔlnT + φ_fr Δlnμ,

and the blob, which keeps the pressure of the surroundings and moves adiabatically with blended thermodynamics

    dln ρ_b = ΔlnP/Γ₁,b,   1/Γ₁,b = f/Γ₁,eq + (1−f)/Γ₁,fr,   δ_b = f δ_eq + (1−f) δ_fr.

The blob is lighter than its surroundings, hence unstable, if A = Δln ρ − Δln ρ_b > 0, so

    A   = α_fr ΔlnP − δ_fr ΔlnT + φ_fr Δlnμ − ΔlnP/Γ₁,b,        sup = A/(δ_b |ΔlnP_hs|).

[derivation] In a hydrostatic layer ΔlnP = ΔlnP_hs < 0, ΔlnT = ∇ ΔlnP, and with the frozen Γ₁ (f = 0, 1/Γ₁ = α − δ∇_ad):
A = −δ ΔlnP (∇ − ∇_ad) + φ Δlnμ, so sup = (∇ − ∇_ad) − (φ/δ) ∇_μ (|ΔlnP_hs| = −ΔlnP) with ∇_μ = Δlnμ/ΔlnP: the Ledoux form of the production
∇ − ∇_ad. Where ΔlnP is far from ΔlnP_hs (behind the accretion shock, ΔlnP → 0) the buoyancy A is still the true density contrast, and the
hydrostatic |ΔlnP_hs| is only the (non-zero) normalisation that turns it into the dimension-less superadiabaticity the KW flux needs.

Frozen blobs (f = 0) never touch the equilibrium arguments: pass anything (NaN is fine) for Γ1_eq, δ_eq then.
"""
@inline function blob_sup(ΔlnP, ΔlnT, Δlnμ, α_fr, δ_fr, φ_fr, Γ1_eq, Γ1_fr, δ_eq, f, ΔlnP_hs)
    if iszero(f)                                   # frozen: structural (not a threshold), and keeps NaN equilibrium data out
        invG = inv(Γ1_fr); δb = δ_fr
    else
        invG = f / Γ1_eq + (1 - f) / Γ1_fr; δb = f * δ_eq + (1 - f) * δ_fr
    end
    A = α_fr * ΔlnP - δ_fr * ΔlnT + φ_fr * Δlnμ - ΔlnP * invG
    return A / (δb * abs(ΔlnP_hs)), δb
end


# =====================================================================================================================
# M3a: the time-dependent convection (TDC) model, docs/POPIII_CONVECTION.md D2 and "M3a".
#
# Unknown w >= 0 [cm/s] on the interfaces (staggered like u). Index convention = TITAN's: zone j (1..nz) lies between
# interfaces j and j+1 (1..nz+1); interface k is between zones k-1 and k; its control volume (CV) spans the two zone
# centres k-1 .. k, mass m_k = (mz[k-1] + mz[k])/2 with mz = rho*dvol (the momentum row's mk). The CV faces are the zone
# centres. dmdt[k] is the mass flux at interface k with TITAN's sign: POSITIVE = INWARD (towards lower index). The flux
# through the zone-centre face of zone j is Phi_j = (dmdt[j] + dmdt[j+1])/2 (as in the momentum row).
# =====================================================================================================================

"radiation constant a_rad * c [erg cm^-3 K^-4 cm s^-1] = 4 sigma_SB (cgs)"
const AC_CGS = 7.565733250033928e-15 * 2.99792458e10

"""
    TDCPar(; A_s=1, A_d=1, A_c=1, α_p=2/3, α_t=0.01, α_D=1, α_m=0)

Coefficients of the TDC model. A_s source, A_d dissipation, A_c enthalpy-flux prefactor (A_s = A_d = A_c = 1: steady local
root = KW/MLT exactly). α_p turbulent pressure P_t = α_p rho w^2, α_t turbulent-energy diffusion (overshoot), α_D species
mixing D_mix = α_D ℓ w, α_m eddy viscosity (not used by the rows; reserved, 0). Defaults and the mapping to RSP/MESA/Kuhfuss
constants: docs/POPIII_CONVECTION.md "M3a".
"""
Base.@kwdef struct TDCPar{T<:Real}
    A_s::T = 1.0
    A_d::T = 1.0
    A_c::T = 1.0
    α_p::T = 2 / 3
    α_t::T = 0.01
    α_D::T = 1.0
    α_m::T = 0.0
end

"""
    fb(w, F) = w + F − √(w² + F²)

Fischer–Burmeister function: fb = 0 ⇔ w ≥ 0, F ≥ 0, w F = 0. For F ≥ 0 the cancellation-free form w − w²/(F + √(w²+F²)) is
used. At exactly (0,0) the kink is resolved by the element (1,1) of the generalised Jacobian (returns w + F).
"""
@inline function fb(w, F)
    s2 = w * w + F * F
    iszero(s2) && return w + F
    s = sqrt(s2)
    return F >= 0 ? w - w * w / (F + s) : w + F - s
end

"""
    (a, U) = tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κ, T, ac = AC_CGS)

a = ℓ√(gδ/(8H_P)) [cm/s], U = 3 a_rad c T³/(c_P ρ² κ ℓ²) √(8H_P/(gδ)) (dimensionless radiative-loss parameter of KW).
"""
@inline function tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κ, T, ac = AC_CGS)
    a = ℓ * sqrt(g * δ / (8H_P))
    U = 3ac * T^3 / (c_P * ρ^2 * κ * ℓ^2) * sqrt(8H_P / (g * δ))
    return a, U
end

"""
    tdc_local(w, a, U, sup, ℓ, par) = (A_s a² (sup − 2U w/a) − A_d w²)/ℓ

The local (source − dissipation − radiative loss) part of Dw/Dt, [cm/s²]. Steady root for A_s = A_d: w = a ζ.
"""
@inline tdc_local(w, a, U, sup, ℓ, par::TDCPar) = (par.A_s * a^2 * (sup - 2U * w / a) - par.A_d * w^2) / ℓ

"""
    tdc_flux(w, a, U, sup, ℓ, ρ, c_P, T, H_P, par) = A_c ρ c_P T w (ℓ/2H_P) Ỹ,  Ỹ = sup − 2U w/a

Convective (enthalpy) flux [erg cm⁻² s⁻¹], negative in overshoot (sup < 0, w > 0).
"""
@inline tdc_flux(w, a, U, sup, ℓ, ρ, c_P, T, H_P, par::TDCPar) = par.A_c * ρ * c_P * T * w * (ℓ / (2H_P)) * (sup - 2U * w / a)

"turbulent pressure of a zone from the two interface values: α_p ρ ē, ē = (w_k² + w_{k+1}²)/2"
@inline tdc_pturb(w1, w2, ρ, par::TDCPar) = par.α_p * ρ * 0.5 * (w1^2 + w2^2)
"species mixing diffusivity D_mix = α_D ℓ w [cm² s⁻¹]"
@inline tdc_dmix(w, ℓ, par::TDCPar) = par.α_D * ℓ * w

"""
    ζ = tdc_mlt_zeta(U, sup) = √(U² + sup) − U  (written as sup/(√(U²+sup)+U), cancellation-free), 0 for sup ≤ 0.

Reference (MLT/KW) quantities for sup > 0 [a branch on the sign of the *input*, the reference only]:
`tdc_mlt_w(a,U,sup) = a ζ`, `kw_mlt_flux(ρ,c_P,T,g,δ,ℓ,H_P,U,sup) = ρ c_P T √(gδ) ℓ² ζ³/(4√2 H_P^{3/2})`.
"""
@inline tdc_mlt_zeta(U, sup) = sup > 0 ? sup / (sqrt(U^2 + sup) + U) : zero(sup + U)
@inline tdc_mlt_w(a, U, sup) = a * tdc_mlt_zeta(U, sup)
@inline kw_mlt_flux(ρ, c_P, T, g, δ, ℓ, H_P, U, sup) = ρ * c_P * T * sqrt(g * δ) * ℓ^2 * tdc_mlt_zeta(U, sup)^3 / (4 * sqrt(2.0) * H_P * sqrt(H_P))

"""
    TDCGeom(nz, r, Az, mzn, mzo, dmdt, dVdt)

Geometry/mass of the staggered grid (nz zones, nz+1 interfaces). r: interface radii (nz+1); Az: zone-centre area (nz; 4π r_c²
for μ=2, 1 for plane); mzn, mzo: zone masses ρ dvol at the new and old time level (nz); dmdt: mass flux at the interfaces (nz+1),
TITAN's sign (positive inward), the same ṁ as the momentum row; dVdt: Lagrangian volume rate of change of each zone (nz), for the
P_t dV work in `tdc_exchange` (zeros for a static grid).
"""
struct TDCGeom
    nz::Int
    r::AbstractVector
    Az::AbstractVector
    mzn::AbstractVector
    mzo::AbstractVector
    dmdt::AbstractVector
    dVdt::AbstractVector
end

"""
    TDCBack(a, U, sup, ell, rhoz, ellz)

Background: per interface (nz+1) a, U, sup (= ∇−∇_ad Ledoux/blob form, `blob_sup`) and ℓ; per zone (nz) ρ and ℓ (zone-centred,
used by the diffusion flux ρℓ w ∇w and by P_t).
"""
struct TDCBack
    a::AbstractVector
    U::AbstractVector
    sup::AbstractVector
    ell::AbstractVector
    rhoz::AbstractVector
    ellz::AbstractVector
end

"""
    tdc_facew(Φ, wlo, whi, cen) 

Face value of w at a zone centre between the interfaces (wlo = lower index, whi = upper). Φ ≥ 0 is inward flow (upper donor):
first-order upwind (donor) value, blended with the centred value by `cen` ∈ [0,1] (cen = 0: pure upwind, cen = 1: centred,
2nd order; the θ-weighting of the momentum row). No limiter (the FB row keeps w ≥ 0 but a centred face can overshoot).
"""
"""Default upwind decision of the zone-centre face j (local zone index): the sign of the mass flux itself. The TITAN residual passes
`TitanResidFast.decide(:up, 30, …)` instead, so that the donor side is recorded at the base state and held during the Jacobian / FD
J·v evaluations, like every other upwind and limiter branch of the residual (`advect!`)."""
@inline _upraw(j, x::Bool) = x
@inline function tdc_facew(Φ, wlo, whi, cen, isup::Bool = Φ >= 0)
    don = isup ? whi : wlo
    return cen == 0 ? don : (1 - cen) * don + cen * 0.5 * (wlo + whi)
end

"""
    tdc_F!(F, w, wo, G, B, dt, par; cen = 0.0)

F_k [cm/s], interior interfaces k = 2..nz (F[1], F[nz+1] set to w, the Dirichlet rows). With m_n, m_o the CV masses, Φ_hi/Φ_lo the
mass fluxes through the upper/lower zone-centre faces (inward positive), wb the face values, A the zone-centre areas and
G_j = α_t ρ_j ℓ_j w̄_j (w_{j+1} − w_j)/Δr_j, w̄_j = (w_j + w_{j+1})/2:

    F_k = [ (m_n w_k − m_o wo_k) − Δt (Φ_hi wb_hi − Φ_lo wb_lo) − w_k ((m_n − m_o) − Δt (Φ_hi − Φ_lo)) ] / m_n
          − Δt [ (A_hi G_hi − A_lo G_lo)/m_n + α_t ℓ_k (g_lo² + g_hi²)/2 ]  − Δt R_loc(w_k),    g_j = (w_{j+1} − w_j)/Δr_j.

The first bracket is the interface-CV advection exactly as in the GAS MOMENTUM row (conservative form, same ṁ and face values)
minus w_k times the discrete continuity residual, i.e. the advective form m_o(w_n − w_o) − Δt[Φ_hi (wb_hi − w_k) − Φ_lo (wb_lo − w_k)]
when continuity closes. The diffusion is (1/ρ)∇·(α_t ρ ℓ w ∇w) + α_t ℓ |∇w|² of the doc.
"""
function tdc_F!(F, w, wo, G::TDCGeom, B::TDCBack, dt, par::TDCPar; cen = 0.0, up = _upraw)
    nz = G.nz; r = G.r; dm = G.dmdt
    F[1] = w[1]; F[nz+1] = w[nz+1]
    @inbounds for k in 2:nz
        mn = 0.5 * (G.mzn[k-1] + G.mzn[k]); mo = 0.5 * (G.mzo[k-1] + G.mzo[k])
        Φh = 0.5 * (dm[k] + dm[k+1]); Φl = 0.5 * (dm[k-1] + dm[k])
        wbh = tdc_facew(Φh, w[k], w[k+1], cen, up(k, Φh >= 0)); wbl = tdc_facew(Φl, w[k-1], w[k], cen, up(k - 1, Φl >= 0))
        adv = ((mn * w[k] - mo * wo[k]) - dt * (Φh * wbh - Φl * wbl) - w[k] * ((mn - mo) - dt * (Φh - Φl))) / mn
        gh = (w[k+1] - w[k]) / (r[k+1] - r[k]); gl = (w[k] - w[k-1]) / (r[k] - r[k-1])
        Gh = par.α_t * B.rhoz[k] * B.ellz[k] * 0.5 * (w[k] + w[k+1]) * gh
        Gl = par.α_t * B.rhoz[k-1] * B.ellz[k-1] * 0.5 * (w[k-1] + w[k]) * gl
        dif = (G.Az[k] * Gh - G.Az[k-1] * Gl) / mn + par.α_t * B.ell[k] * 0.5 * (gl^2 + gh^2)
        F[k] = adv - dt * dif - dt * tdc_local(w[k], B.a[k], B.U[k], B.sup[k], B.ell[k], par)
    end
    return F
end

"""
    tdc_rows!(φ, w, wo, G, B, dt, par; cen = 0.0)

The complementarity rows φ_k = fb(w_k, F_k) (interior k), φ_1 = w_1, φ_{nz+1} = w_{nz+1} (w = 0 boundaries). `φ` must have the
element type of `w` (ForwardDiff). The Jacobian is tridiagonal.
"""
function tdc_rows!(φ, w, wo, G::TDCGeom, B::TDCBack, dt, par::TDCPar; cen = 0.0)
    F = similar(φ)
    tdc_F!(F, w, wo, G, B, dt, par; cen = cen)
    nz = G.nz
    φ[1] = w[1]; φ[nz+1] = w[nz+1]
    @inbounds for k in 2:nz; φ[k] = fb(w[k], F[k]); end
    return φ
end

"""
    (Q, parts) = tdc_exchange(w, wo, G, B, dt, par; cen = 0.0)

Realised exchange between the gas and the turbulent energy e = w², per interface CV (all nz+1 CVs; the two boundary CVs hold half
a zone and e = 0):

    Q_k = (m_n e_k − m_o eo_k)/Δt + (Fe_hi − Fe_lo) + (Fd_hi − Fd_lo) + P_dV,k

with the zone-centre (outward positive) fluxes Fe_j = −Φ_j wb_j² (advective, the same face value as the row) and
Fd_j = −α_t A_j ρ_j ℓ_j w̄_j (e_{j+1} − e_j)/Δr_j (diffusive), both zero outside the domain, and P_dV,k = (Pt_{k−1} dVdt_{k−1} + Pt_k dVdt_k)/2
with Pt_j = α_p ρ_j ē_j (keyword `Ptc` (M3b): the zone P_t the momentum row used, e.g. θ-centred between the old and new level, in the P_dV term only). Q_k > 0: e gains net from the gas (production − dissipation, as realised by the row, including its
discretisation error). Telescoping: Σ_k Q_k Δt = Δ Σ_k m e + Δt Σ_j Pt_j dVdt_j exactly (`tdc_budget_check`). Returns Q and a
NamedTuple of the pieces.
"""
function tdc_exchange(w, wo, G::TDCGeom, B::TDCBack, dt, par::TDCPar; cen = 0.0, Ptc = nothing, up = _upraw)
    nz = G.nz; r = G.r; dm = G.dmdt; T = promote_type(eltype(w), eltype(wo), eltype(G.mzn))
    Fe = zeros(T, nz); Fd = zeros(T, nz); Pt = zeros(T, nz)
    @inbounds for j in 1:nz
        Φ = 0.5 * (dm[j] + dm[j+1])
        wb = tdc_facew(Φ, w[j], w[j+1], cen, up(j, Φ >= 0))
        Fe[j] = -Φ * wb^2
        Fd[j] = -par.α_t * G.Az[j] * B.rhoz[j] * B.ellz[j] * 0.5 * (w[j] + w[j+1]) * (w[j+1]^2 - w[j]^2) / (r[j+1] - r[j])
        Pt[j] = tdc_pturb(w[j], w[j+1], B.rhoz[j], par)
    end
    Q = zeros(T, nz + 1); dedt = zeros(T, nz + 1); pdv = zeros(T, nz + 1)
    Pc = Ptc === nothing ? Pt : Ptc          # M3b: the θ-centred zone P_t of the momentum row (default: the new-level P_t, as M3a)
    @inbounds for k in 1:nz+1
        mn = (k > 1 ? G.mzn[k-1] : 0.0) + (k <= nz ? G.mzn[k] : 0.0); mo = (k > 1 ? G.mzo[k-1] : 0.0) + (k <= nz ? G.mzo[k] : 0.0)
        mn *= 0.5; mo *= 0.5
        dedt[k] = (mn * w[k]^2 - mo * wo[k]^2) / dt
        fhe = k <= nz ? Fe[k] + Fd[k] : zero(T); fle = k > 1 ? Fe[k-1] + Fd[k-1] : zero(T)
        pdv[k] = 0.5 * ((k > 1 ? Pc[k-1] * G.dVdt[k-1] : zero(T)) + (k <= nz ? Pc[k] * G.dVdt[k] : zero(T)))
        Q[k] = dedt[k] + (fhe - fle) + pdv[k]
    end
    return Q, (Fe = Fe, Fd = Fd, Pt = Pt, dedt = dedt, pdv = pdv)
end

"""
    tdc_gas_exchange(Q) → gas-zone source (nz): zone j gets −(Q_j + Q_{j+1})/2, the two boundary CVs (k = 1, nz+1) entirely to
their only zone, so that Σ_j gas_j = −Σ_k Q_k (energy conservation).
"""
function tdc_gas_exchange(Q)
    nz = length(Q) - 1
    g = [-0.5 * (Q[j] + Q[j+1]) for j in 1:nz]
    g[1] -= 0.5 * Q[1]; g[nz] -= 0.5 * Q[nz+1]
    return g
end

"""
    (resid, scale) = tdc_budget_check(Q, w, wo, G, B, dt, par; cen = 0.0)

resid = Σ_k Q_k Δt − [Σ_k (m_n e)_k − Σ_k (m_o e)_k] − Δt Σ_j Pt_j dVdt_j (all the e fluxes telescope away; the domain is closed,
nothing leaves through the boundary CVs, whose e = 0 and Dirichlet); scale = Σ_k (m_n e + m_o e) + Δt Σ|Q_k|. Round-off if
the advective and diffusive fluxes telescope.
"""
function tdc_budget_check(Q, w, wo, G::TDCGeom, B::TDCBack, dt, par::TDCPar; cen = 0.0)
    nz = G.nz
    mn = [0.5 * ((k > 1 ? G.mzn[k-1] : 0.0) + (k <= nz ? G.mzn[k] : 0.0)) for k in 1:nz+1]
    mo = [0.5 * ((k > 1 ? G.mzo[k-1] : 0.0) + (k <= nz ? G.mzo[k] : 0.0)) for k in 1:nz+1]
    EN = sum(mn .* w .^ 2); Eo = sum(mo .* wo .^ 2)
    pdv = sum(tdc_pturb(w[j], w[j+1], B.rhoz[j], par) * G.dVdt[j] for j in 1:nz)
    resid = dt * sum(Q) - (EN - Eo) - dt * pdv
    return resid, EN + Eo + dt * sum(abs, Q)
end

end # module
