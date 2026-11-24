"""
    LResp

M6 (docs/POPIII_CONVECTION.md): the linear response of the network's own chemistry to a convective blob, as pure number-type-generic functions (Float64 or ForwardDiff Duals of any
nesting; nothing here reads global state). A blob displaced at angular frequency ω has δρ, δT e^{iωt}; the coupled species x (the solver's x_H⁺, x_H₂) follow

    iω δx = S_ρ δρ + S_T δT + J δx   ⇒   δx = (iωI − J)⁻¹ (S_ρ δρ + S_T δT) ≡ X_ρ δρ + X_T δT,        J = ∂S/∂x, S = dx/dt (the CK rates).

With the EOS derivatives at fixed x and the composition derivatives p_x, e_x, the dynamic (complex) derivatives are  P_ρ^ω = P_ρ + p_x·X_ρ,  P_T^ω = P_T + p_x·X_T, e_ρ^ω, e_T^ω likewise, and
c_P, δ, ∇_ad, α, Γ₁ follow from the usual formulas of the complex derivatives (the adiabat de = (P/ρ²)dρ closes the system for δT/δρ). The blob uses the real parts (the part of the
response in phase with the displacement). ω → 0 is the equilibrium response of the network (latent heat), ω → ∞ the frozen composition.
"""
module LResp

export CC, lresp, lresp_frozen, lresp_derivs, lresp_tau, ARAD_LR

"Radiation constant [erg cm⁻³ K⁻⁴] as `conv_thermo` (PopIIIMaterial)."
const ARAD_LR = 7.565723e-15

"Complex number over any real type (also Duals): only what the response needs."
struct CC{R<:Real}
    re::R
    im::R
end
CC(a::Real, b::Real) = (p = promote(a, b); CC{typeof(p[1])}(p[1], p[2]))
CC(a::Real) = CC(a, zero(a))
Base.:+(a::CC, b::CC) = CC(a.re + b.re, a.im + b.im)
Base.:-(a::CC, b::CC) = CC(a.re - b.re, a.im - b.im)
Base.:-(a::CC) = CC(-a.re, -a.im)
Base.:*(a::CC, b::CC) = CC(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
Base.:*(a::Real, b::CC) = CC(a * b.re, a * b.im)
Base.:*(a::CC, b::Real) = b * a
Base.:+(a::Real, b::CC) = CC(a + b.re, b.im)
Base.:+(a::CC, b::Real) = b + a
Base.:-(a::CC, b::Real) = CC(a.re - b, a.im)
Base.:-(a::Real, b::CC) = CC(a - b.re, -b.im)
function Base.:/(a::CC, b::CC)
    n = b.re * b.re + b.im * b.im
    CC((a.re * b.re + a.im * b.im) / n, (a.im * b.re - a.re * b.im) / n)
end
Base.:/(a::CC, b::Real) = CC(a.re / b, a.im / b)
Base.:/(a::Real, b::CC) = CC(a) / b

"""
    lresp_derivs(d, ω) -> (Pρ, PT, eρ, eT)   (complex, radiation included)

`d` = (; ρ, T, p, u, dp::NTuple{4}, du::NTuple{4}, dS1::NTuple{4}, dS2::NTuple{4}) with the gas pressure p, the gas internal energy per mass u and the rates S₁, S₂ of the coupled species
(x₁ = x_H⁺, x₂ = x_H₂) and their derivatives with respect to (ρ, T, x₁, x₂) (the CK EOS `cpl_g` rows, the neutral hydrogen following x₁, x₂ as 1 − x₁ − 2x₂). `ω = nothing`: frozen composition.
"""
function lresp_derivs(d, ω)
    ρ, T = d.ρ, d.T
    a4 = ARAD_LR * T^4
    Pρ = d.dp[1]; PT = d.dp[2] + 4ARAD_LR * T^3 / 3
    eρ = d.du[1] - a4 / ρ^2; eT = d.du[2] + 4ARAD_LR * T^3 / ρ
    if ω === nothing
        return CC(Pρ), CC(PT), CC(eρ), CC(eT)
    end
    m11 = -d.dS1[3]; m12 = -d.dS1[4]; m21 = -d.dS2[3]; m22 = -d.dS2[4]          # A = iω I + M, M = −J
    Dr = m11 * m22 - m12 * m21 - ω^2; Di = ω * (m11 + m22); det = CC(Dr, Di)
    function resp(s1, s2)                                                          # A⁻¹ s = adj(A) s / det A
        n1 = CC(m22 * s1 - m12 * s2, ω * s1); n2 = CC(-m21 * s1 + m11 * s2, ω * s2)
        n1 / det, n2 / det
    end
    X1ρ, X2ρ = resp(d.dS1[1], d.dS2[1]); X1T, X2T = resp(d.dS1[2], d.dS2[2])
    Pρ_ = Pρ + d.dp[3] * X1ρ + d.dp[4] * X2ρ; PT_ = PT + d.dp[3] * X1T + d.dp[4] * X2T
    eρ_ = eρ + d.du[3] * X1ρ + d.du[4] * X2ρ; eT_ = eT + d.du[3] * X1T + d.du[4] * X2T
    Pρ_, PT_, eρ_, eT_
end

function _thermo(d, dv)
    Pρ, PT, eρ, eT = dv
    ρ, T = d.ρ, d.T; P = d.p + ARAD_LR * T^4 / 3
    δ = (T / ρ) * PT / Pρ                                    # (∂ ln ρ/∂ ln T)_P = −P_T T/(ρ P_ρ): δ = T P_T/(ρ P_ρ)
    cP = eT - (eρ - P / ρ^2) * PT / Pρ
    α = P / (ρ * Pρ)                                         # (∂ ln ρ/∂ ln P)_T
    nad = P * δ / (T * ρ * cP)
    Γ = 1 / (α - δ * nad)
    (; gam1 = Γ.re, del = δ.re, cp = cP.re, nad = nad.re, alpha = α.re)
end
"""
    lresp(d, ω) -> (; gam1, del, cp, nad, alpha)   real parts of the dynamic Γ₁(ω), δ(ω), c_P(ω), ∇_ad(ω), α(ω)  (see the module docstring); (the complex values: `_thermo(d, lresp_derivs(d, ω))` is not exported; the unit test forms them from `lresp_derivs`).
"""
lresp(d, ω) = _thermo(d, lresp_derivs(d, ω))
"The frozen-composition limit ω → ∞ (the CK EOS at fixed x₁, x₂, the fast species following ρ, T)."
lresp_frozen(d) = _thermo(d, lresp_derivs(d, nothing))
"""
    lresp_tau(d) -> (τ_max, τ_min)  relaxation times 1/|λ| of the eigenvalues of J = ∂S/∂x (2×2, real or complex pair); Float64 diagnostics (CONVDIAG ωτ_max).
"""
function lresp_tau(d)
    a, b, c, e = d.dS1[3], d.dS1[4], d.dS2[3], d.dS2[4]
    tr = a + e; dt = a * e - b * c; disc = tr^2 / 4 - dt
    if disc >= 0
        s = sqrt(disc); l1 = abs(tr / 2 + s); l2 = abs(tr / 2 - s)
    else
        l1 = l2 = sqrt(dt)
    end
    (inv(min(l1, l2)), inv(max(l1, l2)))
end

end # module
