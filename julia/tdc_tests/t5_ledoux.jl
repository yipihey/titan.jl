include(joinpath(@__DIR__, "harness.jl"))
open_log("t5_ledoux")
lp("T5: Ledoux — frozen blob (f = 0) stable, equilibrium blob (f = 1) unstable; units ℓ = a = 1")
lp("    NO semiconvection model: a layer with sup(f=0) < 0 < sup(f=1) is convective or not purely according to f.")
# a Schwarzschild-unstable (∇ = 0.45 > ∇_ad = 0.4 of the frozen gas) layer with a stabilising composition gradient
ΔlnP = -0.1; ∇ = 0.45; ΔlnT = ∇ * ΔlnP; Δlnμ = -0.01
αfr = 1.0; δfr = 1.0; φfr = 1.0; Γfr = 1 / 0.6; Γeq = 1 / 0.9; δeq = 1.5
Δlnρ = αfr * ΔlnP - δfr * ΔlnT + φfr * Δlnμ      # M6: blob_sup takes the actual Δlnρ and the blob's Γ₁,b, δ_b
sup0, _ = TDC.blob_sup(ΔlnP, Δlnρ, Γfr, δfr, ΔlnP)
sup1, _ = TDC.blob_sup(ΔlnP, Δlnρ, Γeq, δeq, ΔlnP)
supS = ∇ - 0.4
lp(@sprintf("    ΔlnP = %g ∇ = %g ∇_ad(frozen) = 0.4 Δlnμ = %g: Schwarzschild ∇−∇_ad = %+g;  sup(f=0, frozen Ledoux) = %+.5f;  sup(f=1, equilibrium, Γ1_eq = 1/0.9, δ_eq = 1.5) = %+.5f", ΔlnP, ∇, Δlnμ, supS, sup0, sup1))
check("sup(f=0) < 0 < sup(f=1)", sup0 < 0 < sup1)
par = TDCPar(α_t = 0.0)
nz = 10; G = make_grid(collect(range(0.0, 10.0, length = nz + 1)))
Bf = make_back(G; a = 1.0, U = 0.0, sup = sup0, ell = 1.0)
Bg = make_back(G; a = 1.0, U = 0.0, sup = sup1, ell = 1.0)
# f = 0 from w = 0: exactly zero
let w = zeros(nz + 1); bad = 0
for s in 1:200
    w, info = newton_w(w, G, Bf, 0.37, par); bad += !info.ok
end
check("f=0, w=0 start: w stays exactly 0 (200 steps)", all(w .== 0.0) && bad == 0, @sprintf("max w = %g", maximum(w)))
end
# f = 0 from a seeded w: decays to exactly 0 in finite time, Riccati: dw/dt = -(a²|sup| + w²)/ℓ
wc = sqrt(abs(sup0)); w0 = 0.6
tstar = atan(w0 / wc) / wc
let
lp(@sprintf("\n    decay of a seeded w0 = %g: dw/dt = −(a²|sup| + w²)/ℓ, w = w_c tan(atan(w0/w_c) − w_c t/ℓ), w_c = a√|sup| = %.5f, analytic t* = %.5f ℓ/a (= (ℓ/w_c) atan(w0/w_c))", w0, wc, tstar))
lp(@sprintf("    %-9s %-12s %-12s %-10s %s", "dt", "t_hit (BE)", "t*", "rel diff", "w after hit"))
ok = true
for dt in (0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001)
    w = fill(w0, nz + 1); w[1] = 0.0; w[end] = 0.0; t = 0.0; thit = NaN; zero_after = true
    for s in 1:round(Int, 3 * tstar / dt + 5)
        w, info = newton_w(w, G, Bf, dt, par); t += dt
        if isnan(thit) && w[5] == 0.0; thit = t; elseif !isnan(thit); zero_after &= all(w .== 0.0); end
    end
    # BE time to zero by linear interpolation of the last positive step is not needed; the first step at which w is exactly 0
    lp(@sprintf("    %-9g %-12.5f %-12.5f %-10.2e %s", dt, thit, tstar, abs(thit - tstar) / tstar, zero_after ? "exactly 0 thereafter" : "NOT zero"))
    ok &= zero_after && abs(thit - tstar) <= 4dt
end
check("f=0 seeded: reaches exactly 0 in finite time (first order: |t_hit - t*| <= 4 dt of the analytic t*, ~2.5 dt observed), stays 0", ok)
end
# BE reproduces the tan solution before the hit
let dt = 0.001; w = fill(w0, nz + 1); w[1] = 0.0; w[end] = 0.0; t = 0.0; emax = 0.0
    for s in 1:round(Int, 0.9 * tstar / dt)
        w, _ = newton_w(w, G, Bf, dt, par); t += dt
        emax = max(emax, abs(w[5] - wc * tan(atan(w0 / wc) - wc * t)) / w0)
    end
    lp(@sprintf("    max |w_BE − w_tan|/w0 for t < 0.9 t*, dt = 1e-3: %.2e", emax)); check("BE follows the tan solution (dt = 1e-3)", emax < 2e-3)
end
# f = 1: convection, w -> tanh solution
let ws = sqrt(sup1); dt = 0.02; w = zeros(nz + 1); t = 0.0; emax = 0.0
    for s in 1:round(Int, 8 / (ws * dt))
        w, _ = newton_w(w, G, Bg, dt, par); t += dt
        emax = max(emax, abs(w[5] - ws * tanh(ws * t)) / ws)
    end
    lp(@sprintf("\n    f=1 from w=0: w_s = a√sup = %.5f, final w[5]/w_s = %.8f, max |w − tanh|/w_s = %.2e (BE, dt = 0.02)", ws, w[5] / ws, emax))
    check("f=1: convection grows to the tanh/MLT amplitude", abs(w[5] / ws - 1) < 1e-6 && emax < 2e-2)
end
lp("T5 failures: ", NFAIL[]); close_log(); exit(NFAIL[] == 0 ? 0 : 1)
