include(joinpath(@__DIR__, "harness.jl"))
open_log("t9_mlt")
lp("T9: steady local MLT limit (A_s=A_d=A_c=1)")
par = TDCPar(A_s = 1.0, A_d = 1.0, A_c = 1.0, α_t = 0.0)
# (1) pure local equation: Newton root of R_loc(w) = 0 vs a ζ, over a grid (sup, U)
maxerr = 0.0; n = 0
for U in 10.0 .^ (-3:0.25:3), sup in 10.0 .^ (-8:0.5:2)
    a = 1.7e6; ℓ = 3e9
    wref = tdc_mlt_w(a, U, sup)
    w = wref * (1.7)                       # start off the root
    for it in 1:60
        f(x) = tdc_local(x, a, U, sup, ℓ, par); d = ForwardDiff.derivative(f, w)
        dw = -f(w) / d; w += dw; abs(dw) < 1e-15 * w && break
    end
    global maxerr = max(maxerr, abs(w - wref) / wref); global n += 1
end
check("local root == a*zeta, $n cases (U 1e-3..1e3, sup 1e-8..1e2)", maxerr <= 1e-12, @sprintf("max rel err %.2e", maxerr))
# (2) flux: F_c(w = a zeta) == KW MLT flux, with physical inputs. Ỹ = sup − 2U w/a cancels catastrophically when U² >> sup
#     (Ỹ = ζ² ≈ (sup/2U)²), so the direct evaluation has relative error ~ eps·4U²/sup; checked against that bound,
#     and to 1e-12 where the bound permits. The steady-root identity Ỹ = (A_d/A_s)(w/a)² gives the cancellation-free value.
let maxf = 0.0, maxw = 0.0, m = 0, worst = 0.0, maxfs = 0.0, mtight = 0
    for ℓfac in (0.5, 1.0, 2.0), κ in 10.0 .^ (-3:1.0:6), sup in 10.0 .^ (-6:1.0:1)
        ρ = 3e-9; c_P = 2.5e8; T = 4e4; g = 1e5; δ = 1.1; H_P = 7e9; ℓ = ℓfac * H_P
        a, U = tdc_ay(ℓ, g, δ, H_P, c_P, ρ, κ, T)
        w = tdc_mlt_w(a, U, sup)
        Fc = tdc_flux(w, a, U, sup, ℓ, ρ, c_P, T, H_P, par)
        Fk = kw_mlt_flux(ρ, c_P, T, g, δ, ℓ, H_P, U, sup)
        err = abs(Fc - Fk) / Fk; bound = 1e-12 + 4 * eps() * 4U^2 / sup * 4
        worst = max(worst, err / bound); m += 1
        U^2 / sup <= 1e3 && (maxf = max(maxf, err); mtight += 1)
        Fs = par.A_c * ρ * c_P * T * w * (ℓ / (2H_P)) * (par.A_d / par.A_s) * (w / a)^2     # steady-root form of Ỹ
        maxfs = max(maxfs, abs(Fs - Fk) / Fk)
        maxw = max(maxw, abs(tdc_local(w, a, U, sup, ℓ, par)) * ℓ / (a^2 * sup))
    end
    check("F_c(a zeta) == KW MLT flux (direct form), $mtight cases with U²/sup<=1e3", maxf <= 1e-12, @sprintf("max rel err %.2e", maxf))
    check("F_c direct form within the cancellation bound eps*16U²/sup, all $m cases", worst <= 1, @sprintf("max err/bound %.2e", worst))
    check("F_c steady-root form == KW MLT flux, all $m cases", maxfs <= 1e-12, @sprintf("max rel err %.2e; max |R_loc|ℓ/(a²sup) %.1e", maxfs, maxw))
end
# (3) harness: uniform background, no advection/diffusion, long-time solution, from small seed and from a large seed
for (U, sup) in ((1e-3, 0.3), (1.0, 0.05), (30.0, 1e-3), (1e3, 1e-2), (0.1, 5.0)), w0f in (1e-3, 4.0)
    nz = 12; r = collect(range(0.0, 12.0, length = nz + 1))
    G = make_grid(r); a = 1.0; ℓ = 1.0
    B = make_back(G; a = a, U = U, sup = sup, ell = ℓ)
    wref = tdc_mlt_w(a, U, sup)
    w = fill(w0f * wref, nz + 1); w[1] = 0.0; w[end] = 0.0
    dt = 3.0 * ℓ / (2U * a + 2wref)                     # a few relaxation times (a huge dt only adds cancellation noise in F)
    w, info = evolve(w, G, B, dt, 4000, par; stop = 1e-15)
    err = maximum(abs.(w[2:end-1] .- wref)) / wref
    check(@sprintf("harness long-time U=%g sup=%g seed=%g wref", U, sup, w0f), err <= 1e-12 && info.nbad == 0,
          @sprintf("steps %d  max rel err %.2e  maxNewton %d", info.steps, err, info.maxit))
end
lp("T9 failures: ", NFAIL[]); close_log(); exit(NFAIL[] == 0 ? 0 : 1)
