# Upper bound on H2 made via the H− channel when free electrons are NOT in Saha equilibrium.
#   H + e → H− + γ   (k7, Galli & Palla 1998)
#   H− + H → H2 + e  (k8, associative detachment)
#   H− + γ → H + e   (photodetachment by the quasi-star's radiation, diluted blackbody W·B(T_rad))
# Steady-state H−; H2 destruction ignored (so x_H2 is an upper bound).
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; const Q = QuasiStarInfall
using Printf
const kB = 1.380649e-16; const h = 6.62607015e-27; const c = 2.99792458e10; const mH = 1.6735575e-24
k7(T) = 3.0e-16 * (T / 300)^0.95 * exp(-T / 9320)
k8(T) = 1.3e-9
"Photodetachment rate (s⁻¹) in a diluted blackbody, using the John (1988) σ(H−) already in the code."
function kpd(Trad, W)
    s = 0.0
    λs = exp.(range(log(0.126), log(1.6418), length = 4000))      # μm, threshold 1.6419 μm
    for i in 1:length(λs)-1
        λ = 0.5(λs[i] + λs[i+1]) * 1e-4; dλ = (λs[i+1] - λs[i]) * 1e-4
        ν = c / λ; dν = c / λ^2 * dλ
        B = 2h * ν^3 / c^2 / expm1(h * ν / (kB * Trad))
        s += 4π * B / (h * ν) * Q.σHm_bf(λ * 1e4) * dν
    end
    return W * s
end
@printf("check: undiluted photodetachment at 6000 K = %.2e s⁻¹ (fit 0.11 T^2.13 e^(-8823/T) = %.2e)\n",
        kpd(6000.0, 1.0), 0.11 * 6000^2.13 * exp(-8823 / 6000))

M = 1e6Msun; R = 8.42e15; L = LEdd_es(M); TL = (L / (4π * σSB * R^2))^0.25
println("\nPhotosphere (extended branch, M_BH/M = 0.1): n_H ≈ 2e13, T ≈ T_rad ≈ 7300 K, W = 1/2")
nH = 2e13; T = 7300.0; kp = kpd(T, 0.5)
for xe in (:LTE, 1e-3, 1e-2, 0.1, 0.5)
    ne = xe == :LTE ? Q.saha(0.76 \ 1 * nH * mH / 0.76 * 1, T).ne : xe * nH   # LTE uses Saha at ρ = n_H m_H / X
    nHm = k7(T) * nH * ne / (k8(T) * nH + kp)
    rate = k8(T) * nH * nHm                                              # H2 formation, cm⁻³ s⁻¹
    @printf("  x_e = %-5s: k_pd/(k8 n_H) = %.1e, H2 formation time n_H/rate = %.1e s\n",
            string(xe), kp / (k8(T) * nH), nH / rate)
end
println("  (equilibrium H2 fraction at 7300 K is ~1e-8; formation catalysts cannot push it above that)")

println("\nInfall, Ṁ = 1 Msun/yr, Lucy temperature, upper bound on x_H2 accumulated during fall from r to R")
Mdot = 1Msun / yr
for rr in (2.0, 3.0, 10.0, 30.0, 100.0)
    r = rr * R
    W = 0.5 * (1 - sqrt(1 - (R / r)^2))
    T = TL * W^0.25                                  # optically thin Lucy limit
    n = 0.76 * rho_infall(Mdot, M, r) / mH
    tff = (2 / 3) * r / vff(M, r)                    # time to fall from r to ≈ R (r ≫ R)
    kp = kpd(TL, W)
    out = String[]
    for xe in (1e-3, 1e-2, 0.1, 0.5)
        ne = xe * n; n0 = (1 - xe) * n
        nHm = k7(T) * n0 * ne / (k8(T) * n0 + kp)
        xH2 = k8(T) * n0 * nHm * tff / n
        push!(out, @sprintf("x_e=%.0e→x_H2≤%.1e", xe, xH2))
    end
    @printf("  r=%5.0fR: T_gas=%5.0f K n_H=%.1e W=%.1e k_pd=%.1e s⁻¹ vs k8 n_H=%.1e s⁻¹ | %s\n",
            rr, T, n, W, kp, k8(T) * n, join(out, "  "))
end
