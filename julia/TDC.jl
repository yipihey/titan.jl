"""
    TDC

Pure, number-type-generic functions of the sub-shock convection work (docs/POPIII_CONVECTION.md, D1/M2; the time-dependent
convection (TDC) of M3 goes in here too). Nothing here reads global state: the arguments may be Float64 or ForwardDiff Duals
(of any nesting), so the Newton Jacobian sees every dependence.
"""
module TDC

export gamma1_inv, gamma1, blob_sup, blob_f

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
A = −δ ΔlnP (∇ − ∇_ad) + φ Δlnμ, so sup = (∇ − ∇_ad) + (φ/δ) ∇_μ with ∇_μ = Δlnμ/ΔlnP: the Ledoux form of the production
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

end # module
