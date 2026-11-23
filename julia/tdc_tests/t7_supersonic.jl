include(joinpath(@__DIR__, "harness.jl"))
open_log("t7_supersonic")
lp("T7: supersonic inflow through an unstable profile (sup = 1 everywhere, w_s = a√sup = 1, ℓ = a = 1, g = 8, δ = H_P = ρ = c_P = T = 1), w = 0 inflow, no gate")
lp("    estimate: w(x) ≈ w_s tanh(w_s d/(u ℓ)) → w_s (w_s d/ℓ)/u ≪ w_s for u ≫ w_s d/ℓ;  F_c/F_MLT = (w/w_s)(Ỹ/sup)·... ≈ w/w_s (U = 0)")
g = 8.0; HP = 1.0; ℓ = 1.0
Fmlt = kw_mlt_flux(1.0, 1.0, 1.0, g, 1.0, ℓ, HP, 0.0, 1.0)
lp(@sprintf("    F_MLT (sup=1, U=0) = %.6f; F_c(w_s) = %.6f", Fmlt, tdc_flux(1.0, 1.0, 0.0, 1.0, ℓ, 1.0, 1.0, 1.0, HP, TDCPar())))
lp(@sprintf("\n    %-6s %-6s %-5s %-8s | %-12s %-12s %-12s %-12s | %s", "u/w_s", "Ld/ℓ", "α_t", "CFL-dt", "max w/w_s", "w_s·tanh est", "max F_c/F_MLT", "est", "Newton: max its, fail"))
okall = true; rows = []
for ur in (10.0, 30.0, 100.0), Ld in (1.0, 3.0), α in (0.0, 0.1)
    nz = 100; xs = collect(range(0.0, Ld, length = nz + 1)); h = xs[2]
    par = TDCPar(α_t = α); G = make_grid(xs; mdot = ur)
    B = make_back(G; a = 1.0, U = 0.0, sup = 1.0, ell = 1.0)
    for cfl in (5.0, Ld / h * 1.0)
        dt = cfl * h / ur
        w, info = evolve(zeros(nz + 1), G, B, dt, 3000, par; stop = 1e-13)
        wmax = maximum(w); est = tanh(Ld / ur)
        Fc = maximum(tdc_flux(w[k], 1.0, 0.0, 1.0, ℓ, 1.0, 1.0, 1.0, HP, par) for k in 1:nz+1) / Fmlt
        lp(@sprintf("    %-6g %-6g %-5g %-8s | %-12.5e %-12.5e %-12.5e %-12.5e | %d, %d", ur, Ld, α, cfl == 5.0 ? "CFL 5" : "CFL nz", wmax, est, Fc, est, info.maxit, info.nbad))
        push!(rows, (ur, Ld, α, cfl, wmax, est, Fc, info.maxit, info.nbad))
        cfl == 5.0 && (global okall &= wmax < 1.05 * est + 1e-12 && wmax <= est * 1.0001 && info.nbad == 0 && Fc < 1.001 * est && wmax / 1.0 < 0.3)
    end
end
write_csv("t7_supersonic.csv", "u_over_ws,Ld,alpha_t,cfl,max_w_over_ws,tanh_estimate,max_Fc_over_Fmlt,max_newton_its,newton_fail", rows)
check("w <= w_s tanh(w_s Ld/(u ℓ)) (≪ w_s for u ≫ w_s Ld/ℓ), F_c/F_MLT <= same, no gate needed, no Newton failure", okall)
lp("T7 failures: ", NFAIL[]); close_log(); exit(NFAIL[] == 0 ? 0 : 1)
