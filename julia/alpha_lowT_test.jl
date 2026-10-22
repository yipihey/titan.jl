# alpha_lowT_test.jl — every term of ContinuumOpacity.alpha is finite, non-negative and continuous in T from 10 K
# to 1e5 K (so RadTransfer needs no temperature floor), and the means are finite in the cold envelope.
#   julia --project=env alpha_lowT_test.jl
using Printf
include(joinpath(@__DIR__, "ContinuumOpacity.jl")); const CO = ContinuumOpacity
pass = Ref(0); fail = Ref(0)
chk(name, ok) = (ok ? (pass[] += 1) : (fail[] += 1); @printf("%-82s %s\n", name, ok ? "PASS" : "FAIL"))
# compositions as smooth functions of T (alpha is linear in each density, so only the coefficients' T dependence is probed)
comp(nH, fH2, xHp, xe, xHm) = (nH = nH * (1 - fH2 - xHp), nH2 = fH2 * nH / 2, nHe = nH / 12, nHp = xHp * nH, ne = xe * nH, nHm = xHm * nH,
                               nHep = 0.1 * xHp * nH / 12, nHepp = 0.01 * xHp * nH / 12)
cases = (("cold envelope n=1e6, f_H2 1e-3", comp(1e6, 1e-3, 1e-4, 1e-4, 1e-12)),
         ("molecular core n=1e10", comp(1e10, 0.5, 1e-10, 4e-10, 1e-14)),
         ("dense n=1e16", comp(1e16, 0.5, 1e-10, 4e-10, 1e-14)),
         ("ionising n=1e20", comp(1e20, 0.05, 1e-3, 1e-3, 1e-14)),
         ("hot ionised n=1e12", comp(1e12, 1e-6, 0.9, 0.9, 1e-16)))
Ts = 10 .^ range(1.0, 5.0; length = 801)
νs = (1e12, 3e12, 1e13, 3e13, 1e14, 3e14, 1e15, 3e15)
ε = 1e-6
for (label, s) in cases
    allfinite = true; cont_ok = true; worst = 0.0; worstinfo = ""
    for ν in νs, T in Ts
        a, sc, p = CO.alpha(ν, T; s...)
        allfinite &= isfinite(a) && isfinite(sc) && a >= 0 && sc >= 0 && all(v -> isfinite(v) && v >= 0, values(p))
        a2, _, p2 = CO.alpha(ν, T * (1 + ε); s...)
        for k in keys(p)
            r = abs(p2[k] - p[k]) / max(p[k], 1e-300)
            p[k] > 1e-60 && r > worst && (worst = r; worstinfo = @sprintf("%s at ν=%.0e T=%.0f", k, ν, T))
        end
        r = abs(a2 - a) / max(a, 1e-300); a > 1e-60 && r > worst && (worst = r; worstinfo = @sprintf("total at ν=%.0e T=%.0f", ν, T))
    end
    chk("$label: all parts finite and ≥ 0, 10 K–1e5 K, 1e12–3e15 Hz", allfinite)
    chk(@sprintf("%s: continuous in T (max relative change for ΔT/T = 1e-6: %.1e, %s)", label, worst, worstinfo), worst < 1e-3)
end
# means at cold temperatures
for (nH, T) in ((1e6, 10.0), (1e8, 30.0), (1e10, 100.0), (1e12, 300.0))
    s = comp(nH, 0.5, 1e-10, 4e-10, 1e-14); ρ = nH * 1.6735575e-24 * (1 + 4 / 12)
    m = CO.means(ρ, T; s...)
    chk(@sprintf("means finite and positive at n_H=%.0e, T=%.0f K: κ_P=%.3e κ_R=%.3e", nH, T, m.κP, m.κR), all(isfinite, (m.κP, m.κR, m.κR_abs)) && m.κP > 0 && m.κR > 0)
end
# the previously defective pieces specifically
chk("John 1988 H⁻ ff held at its 1400 K value below 1400 K (continuous, positive)",
    CO.k_Hm_ff(1.0, 10.0) == CO.k_Hm_ff(1.0, 1400.0) > 0 && isapprox(CO.k_Hm_ff(1.0, 1399.0), CO.k_Hm_ff(1.0, 1401.0); rtol = 1e-2))
chk("LTE H₂⁺ estimate finite and ≤ min(n_H, n_H⁺) down to 10 K", all(T -> (n = CO.FF.nH2p_lte(T, 1e10, 1e2); isfinite(n) && 0 <= n <= 1e2), (10.0, 100.0, 1000.0, 3150.0, 1e4)))
@printf("\n%d passed, %d failed\n", pass[], fail[])
