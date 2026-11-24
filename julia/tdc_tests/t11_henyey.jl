include(joinpath(@__DIR__, "harness.jl"))
open_log("t11_henyey"; dir = normpath(joinpath(@__DIR__, "..", "..", "docs", "logs", "subconv", "henyey")))
lp("T11: Henyey-Vardya-Bodenheimer (1965) thin-eddy radiative losses in tdc_ay (y = 1/3, nu = 8: MESA defaults)")
par = TDCPar(A_s = 1.0, A_d = 1.0, A_c = 1.0, α_t = 0.0)
# fixed physical state; ω = κρℓ is varied through κ
ρ = 3e-9; c_P = 2.5e8; T = 4e4; g = 1e5; δ = 1.1; H_P = 7e9; ℓ = 0.7 * H_P
σ = 5.670374419e-5   # Stefan-Boltzmann [cgs]
check("AC_CGS = 4 σ_SB (a_rad c)", abs(AC_CGS / (4σ) - 1) < 1e-8, @sprintf("ratio-1 = %.1e", AC_CGS / (4σ) - 1))
κof(ω) = ω / (ρ * ℓ)
# (i) thick limit
lp("(i) omega >> 1: U_Henyey/U_KW (y=1/3) = omega^2/(omega^2+3) -> 1")
let worst = 0.0, worst2 = 0.0, mono = true, prev = 0.0
    for ω in 10.0 .^ (1:6)
        a, U = tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κof(ω), T); _, Uk = tdc_ay_kw(ℓ, g, δ, H_P, c_P, ρ, κof(ω), T)
        r = U / Uk; lp(@sprintf("   omega=1e%d  U_H/U_KW = %.12f   1-ratio = %.3e (3/omega^2 = %.3e)", Int(log10(ω)), r, 1 - r, 3 / ω^2))
        worst = max(worst, abs(r - ω^2 / (ω^2 + 3)) / r); worst2 = max(worst2, abs((1 - r) - 3 / (ω^2 + 3)))
        mono &= r > prev; prev = r
    end
    check("U_H/U_KW = ω²/(ω²+3) to round-off, ω=1e1..1e6", worst <= 1e-13 && worst2 <= 1e-15 && mono, @sprintf("max rel dev %.1e, monotone %s", worst, mono))
    _, U1 = tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κof(1e6), T); _, Uk1 = tdc_ay_kw(ℓ, g, δ, H_P, c_P, ρ, κof(1e6), T)
    check("ω=1e6: U_H == U_KW to 3e-12 (KW is the ω→∞ limit)", abs(U1 / Uk1 - 1) <= 3e-12, @sprintf("%.2e", U1 / Uk1 - 1))
    yh = 3 / (4π^2)   # Ludwig+99 / MESA-doc value: thick limit = U_KW/(3y)
    _, U2 = tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κof(1e6), T; y = yh)
    check("y = 3/(4π²): thick limit U_KW/(3y)", abs(U2 / (Uk1 / (3yh)) - 1) <= 1e-10, @sprintf("%.2e", U2 / (Uk1 / (3yh)) - 1))
end
# (ii) thin limit
lp("(ii) omega << 1: U -> a_rad c T^3 kappa l/(c_P a) = 4 sigma T^3 kappa l/(c_P a)  <=>  Gamma_thin = 2 v t_loss/l, t_loss = c_P/(16 sigma T^3 kappa)")
let worst = 0.0, wg = 0.0, wd = 0.0
    for ω in 10.0 .^ (-1:-1:-8), y in (1 / 3, 3 / (4π^2))
        κ = κof(ω); a, U = tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κ, T; y = y)
        Uthin = 4σ * T^3 * κ * ℓ / (c_P * a)
        worst = max(worst, abs(U * (1 + y * ω^2) / Uthin - 1))          # exact: U = U_thin/(1+yω²)
        w = 0.37a; Γ = w / (2U * a); tloss = c_P / (16σ * T^3 * κ)   # eddy efficiency from the Riccati ↔ Γ map (M7.2)
        wg = max(wg, abs(Γ / (1 + y * ω^2) / (2w * tloss / ℓ) - 1))
        dU = ForwardDiff.derivative(k -> tdc_ay(ℓ, g, δ, H_P, c_P, ρ, k, T; y = y)[2], κ)
        h = 1e-6 * κ; fd = (tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κ + h, T; y = y)[2] - tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κ - h, T; y = y)[2]) / 2h
        wd = max(wd, abs(dU / fd - 1))
    end
    check("U = U_thin/(1+yω²), U_thin = 4σT³κℓ/(c_P a); ω=1e-1..1e-8, y∈{1/3,3/4π²}", worst <= 1e-14, @sprintf("max rel dev %.1e", worst))
    check("Γ = (1+yω²)·2 v t_loss/ℓ, t_loss = c_P/(16σT³κ) (thin: Γ_thin = 2 v t_loss/ℓ)", wg <= 1e-14, @sprintf("max rel dev %.1e", wg))
    check("ForwardDiff dU/dκ == centred FD (thin→thick)", wd <= 1e-6, @sprintf("%.1e", wd))
    _, U0 = tdc_ay(ℓ, g, δ, H_P, c_P, ρ, 0.0, T)
    check("κ = 0: U = 0 finite (KW U diverges)", U0 == 0.0, "")
end
# (iii) steady TDC == Henyey-MLT. Two independent references:
#  (a) given ∇ (sup = ∇−∇_ad, what TITAN supplies): Ludwig+99 App. A.1–A.2 solved by Newton on v WITHOUT using U:
#      v² = f1 Λ² gδ (∇−∇e)/H_P,  Γ = (∇−∇e)/(∇e−∇ad) = ρ c_P v τ (1+f4/τ²)/(f3 σ T³),  τ = ω,  ∇−∇e = sup Γ/(1+Γ).
#  (b) given ∇_rad: the Cox&Giuli 14.81–14.82 cubic as coded in MESA calc_MLT (a0 depends on ω); its (v, ∇) must satisfy our steady relation.
function henyey_v_given_grad(sup, ρ, c_P, T, g, δ, H_P, ℓ, κ, y)
    ω = κ * ρ * ℓ; f1 = 1 / 8; f3 = 8 / y; f4 = 1 / y; σ = 5.670374419e-5
    cΓ = ρ * c_P * ω * (1 + f4 / ω^2) / (f3 * σ * T^3)
    F(v) = v^2 - f1 * ℓ^2 * g * δ / H_P * sup * (cΓ * v) / (1 + cΓ * v)
    v = sqrt(f1 * ℓ^2 * g * δ / H_P * sup)                              # upper bound (Γ → ∞)
    for _ in 1:200; dv = -F(v) / ForwardDiff.derivative(F, v); v += dv; abs(dv) < 1e-15v && break; end
    return v, cΓ * v
end
function mesa_cubic(supr, ρ, c_P, T, g, δ, H_P, ℓ, κ, y)   # supr = ∇_rad − ∇_ad
    ac = AC_CGS; P = ρ * g * H_P; Q = δ; α = ℓ / H_P; ω = κ * ρ * ℓ
    ff1 = 1 / 8; ff2 = 1 / 2; ff3 = 8 / y; ff4 = 1 / y
    f4w = ff4 / ω^2 + 1; a0 = (3 / 16) * ff2 * ff3 / f4w
    A = 4c_P * sqrt(ff1 * P * Q * ρ) * α * ω * f4w / (ff3 * ac * T^3)
    B3 = A^2 / a0 * supr
    lo, hi = 0.0, max(1.0, cbrt(B3) + 1, B3 * a0 + 1); for _ in 1:300; mid = (lo + hi) / 2; (a0 * mid^3 + mid^2 + mid - a0 * B3 > 0) ? (hi = mid) : (lo = mid); end
    Γ = (lo + hi) / 2
    for _ in 1:4; f = a0 * Γ^3 + Γ^2 + Γ - a0 * B3; Γ -= f / (3a0 * Γ^2 + 2Γ + 1); end
    a = α * sqrt(Q * P / (8ρ)); v = a * Γ / A                       # C&G 14.86b
    Zeta = Γ^3 / B3
    return v, Γ, (1 - Zeta) * supr                                  # Y = ∇ − ∇ad of the solution (C&G 14.79)
end
let ma = 0.0, mfa = 0.0, mb = 0.0, mfb = 0.0, mΓ = 0.0, na = 0, nb = 0
    for y in (1 / 3, 3 / (4π^2)), ω in 10.0 .^ (-4:1.0:6), s in 10.0 .^ (-6:1.0:1), ℓfac in (0.3, 1.0)
        ℓ = ℓfac * H_P; κ = ω / (ρ * ℓ)
        a, U = tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κ, T; y = y)
        # (a)
        w = tdc_mlt_w(a, U, s); v, Γ = henyey_v_given_grad(s, ρ, c_P, T, g, δ, H_P, ℓ, κ, y)
        ma = max(ma, abs(w / v - 1)); mΓ = max(mΓ, abs(w / (2U * a) / Γ - 1))
        Fm = ρ * c_P * T * v * (ℓ / (2H_P)) * (s * Γ / (1 + Γ)); Fk = kw_mlt_flux(ρ, c_P, T, g, δ, ℓ, H_P, U, s)
        Fs = par.A_c * ρ * c_P * T * w * (ℓ / (2H_P)) * (par.A_d / par.A_s) * (w / a)^2
        mfa = max(mfa, abs(Fk / Fm - 1), abs(Fs / Fm - 1)); na += 1
        # (b)
        vb, Γb, Y = mesa_cubic(s, ρ, c_P, T, g, δ, H_P, ℓ, κ, y)
        Y > 0 || continue
        wb = tdc_mlt_w(a, U, Y); mb = max(mb, abs(wb / vb - 1))
        Fb = ρ * c_P * T * vb * (ℓ / (2H_P)) * (Y * Γb / (1 + Γb)); mfb = max(mfb, abs(kw_mlt_flux(ρ, c_P, T, g, δ, ℓ, H_P, U, Y) / Fb - 1)); nb += 1
    end
    check("(a) steady TDC a ζ(U_H, sup) == Henyey-MLT v (Ludwig A.1-2, independent of U), $na cases (ω 1e-4..1e6, sup 1e-6..10, y∈{1/3,3/4π²}, ℓ∈{0.3,1}H_P)", ma <= 1e-9, @sprintf("max rel err %.2e", ma))
    check("(a) Γ = w/(2Ua) == Ludwig Γ(ω, v)", mΓ <= 1e-9, @sprintf("max rel err %.2e", mΓ))
    check("(a) F_c: kw_mlt_flux(U_H) and tdc_flux steady form == f2 ρ c_P v T Λ(∇−∇e)/H_P", mfa <= 1e-9, @sprintf("max rel err %.2e", mfa))
    check("(b) MESA/C&G cubic (a0(ω)): its (v, ∇) satisfies a ζ(U_H, ∇−∇ad) = v, $nb cases", mb <= 1e-8, @sprintf("max rel err %.2e", mb))
    check("(b) F_c flux of the cubic solution == kw_mlt_flux(U_H, ∇−∇ad)", mfb <= 1e-8, @sprintf("max rel err %.2e", mfb))
end
# Newton root of the local equation with the Henyey U (as T9), and harness long-time with a thin-eddy U
let maxerr = 0.0
    for ω in 10.0 .^ (-3:1.0:4), sup in 10.0 .^ (-6:1.0:1)
        a, U = tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κof(ω), T); wref = tdc_mlt_w(a, U, sup); w = 1.7wref
        for it in 1:80
            f(x) = tdc_local(x, a, U, sup, ℓ, par); dw = -f(w) / ForwardDiff.derivative(f, w); w += dw; abs(dw) < 1e-15w && break
        end
        maxerr = max(maxerr, abs(w - wref) / wref)
    end
    check("local root == a ζ with Henyey U (ω 1e-3..1e4)", maxerr <= 1e-12, @sprintf("max rel err %.2e", maxerr))
end
# the thin limit keeps convection alive: ζ(U_H) finite as ω→0, while KW U → ∞ kills it
let sup = 1.0
    a, U = tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κof(1e-6), T); _, Uk = tdc_ay_kw(ℓ, g, δ, H_P, c_P, ρ, κof(1e-6), T)
    ζH = tdc_mlt_zeta(U, sup); ζK = tdc_mlt_zeta(Uk, sup)
    check("ω=1e-6, sup=1: ζ_Henyey ≈ √sup, ζ_KW suppressed", abs(ζH / sqrt(sup) - 1) < 1e-2 && ζK < 1e-6 * sqrt(sup), @sprintf("ζ_H/√sup = %.6f, ζ_KW/√sup = %.2e (U_KW = %.1e, U_H = %.1e)", ζH / sqrt(sup), ζK / sqrt(sup), Uk, U))
end
lp("T11 failures: ", NFAIL[]); close_log(); exit(NFAIL[] == 0 ? 0 : 1)
