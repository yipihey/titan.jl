include(joinpath(@__DIR__, "harness.jl"))
open_log("mesa_compare")
lp("MESA TDC pointwise ODE tests (Jermyn et al. 2023, Sect. 3, eqs 8-13) vs the implicit FB solver of julia/TDC.jl")
lp("Variable map [derived]: MESA w_J = sqrt(3/2) w (w_J^2 = e_t, w = our velocity = MLT v_c).  Our local equation")
lp("   Dw/Dt = (A_s a² (sup − 2U w/a) − A_d w²)/ℓ   <=>   2 Dw_J/Dt = ξ0 + ξ1 w_J + ξ2 w_J²  with")
lp("   ξ0 = 2√(3/2) A_s a² sup/ℓ,   ξ1 = −4 A_s U a/ℓ,   ξ2 = −2√(2/3) A_d/ℓ   (J² = ξ1² − 4 ξ0 ξ2)")
lp("   (MESA: ξ2 = −(8√2/(3√3)) α_D/(α h) = −2.177 α_D/ℓ  =>  A_d = (4/3) α_D;  ξ0 ∝ α_S => A_s = (4/3) α_S)")
xi(a, U, sup, ℓ, As, Ad) = (2sqrt(1.5) * As * a^2 * sup / ℓ, -4As * U * a / ℓ, -2sqrt(2 / 3) * Ad / ℓ)
"MESA eqs (12)/(13) literally (w_J, frozen coefficients, exact over the step). tanh branch for w0 below the plateau, coth above; tan branch (13) for J² < 0, clipped at 0."
function mesa_step(wJ0, t, ξ0, ξ1, ξ2)
    J2 = ξ1^2 - 4ξ0 * ξ2
    if J2 > 0
        J = sqrt(J2); x0 = (-2ξ2 * wJ0 - ξ1) / J             # tanh(λ/4) = x0
        if abs(x0) < 1
            λ = 4atanh(x0); return -(J * tanh((λ + J * t) / 4) + ξ1) / (2ξ2)
        else
            λ = 4atanh(1 / x0); return -(J / tanh((λ + J * t) / 4) + ξ1) / (2ξ2)    # coth branch (w0 above the plateau)
        end
    else
        # (13) re-derived [derived]: w − c = q' tan θ, c = −ξ1/(2ξ2), q' = |J|/(2|ξ2|), θ' = −|J|/4 (ξ2 < 0);
        # w = −[ξ1 + |J| tan(θ0 − |J| t/4)]/(2ξ2); reaches 0 at t* and stays 0 (MESA clips there)
        Jm = sqrt(-J2); c = -ξ1 / (2ξ2); qp = Jm / (2abs(ξ2)); θ0 = atan((wJ0 - c) / qp)
        tstar = 4 * (θ0 - atan(-c / qp)) / Jm
        return t >= tstar ? zero(wJ0) : -(ξ1 + Jm * tan(θ0 - Jm * t / 4)) / (2ξ2)
    end
end
rk4(w, T, f; n = 20000) = (h = T / n; for i in 1:n; k1 = f(w); k2 = f(w + h / 2 * k1); k3 = f(w + h / 2 * k2); k4 = f(w + h * k3); w += h / 6 * (k1 + 2k2 + 2k3 + k4); end; w)
cases = [
    ("A growth,  U=0.5 sup=+0.8 A_s=A_d=1  w0 small", 1.0, 0.5, 0.8, 1.0, 1.0, 1.0, 0.01),
    ("B overshoot of plateau, same, w0 = 3", 1.0, 0.5, 0.8, 1.0, 1.0, 1.0, 3.0),
    ("C RSP-standard A_s=A_d=4/3, U=0, sup=+0.3", 1.0, 0.0, 0.3, 1.0, 4 / 3, 4 / 3, 0.05),
    ("D alpha_D=2 (A_d = 8/3, A_s = 4/3), U=0.2, sup=0.5", 1.0, 0.2, 0.5, 1.0, 4 / 3, 8 / 3, 0.05),
]
const OK = Ref(true)
for (name, a, U, sup, ℓ, As, Ad, w0) in cases
    par = TDCPar(A_s = As, A_d = Ad, α_t = 0.0)
    ξ0, ξ1, ξ2 = xi(a, U, sup, ℓ, As, Ad)
    wplat = tdc_mlt_w(a, U, sup)            # for A_s = A_d only
    # plateau of MESA (12) as dt -> inf, mapped back
    J = sqrt(ξ1^2 - 4ξ0 * ξ2); wJp = -(J + ξ1) / (2ξ2); wp = wJp / sqrt(1.5)
    wp_ours = (-As * U * a + sqrt((As * U * a)^2 + As * Ad * a^2 * sup)) / Ad         # root of A_d w² + 2 A_s U a w − A_s a² sup = 0
    lp("\n", name)
    lp(@sprintf("    ξ0 = %.5f ξ1 = %.5f ξ2 = %.5f J = %.5f;  plateau MESA (12) mapped back %.12f   our root %.12f   (A_s=A_d: a ζ = %.12f)", ξ0, ξ1, ξ2, J, wp, wp_ours, wplat))
    OK[] &= abs(wp - wp_ours) <= 1e-13 * wp
    lp(@sprintf("    %-10s %-12s %-14s %-14s %-14s", "T (=τ units)", "n BE steps", "w_BE(T)", "w_exact(T)", "rel err  (MESA (12) one step == exact)"))
    τ = ℓ / wp
    for Tτ in (0.1, 1.0, 5.0), n in (1, 4, 16, 64)
        T = Tτ * τ; f(w) = tdc_local(w, a, U, sup, ℓ, par)
        wex = rk4(w0, T, f)
        wm = mesa_step(sqrt(1.5) * w0, T, ξ0, ξ1, ξ2) / sqrt(1.5)
        wb = w0; for i in 1:n; wb = local_step(wb, a, U, sup, ℓ, T / n, par); end
        lp(@sprintf("    %-10g %-12d %-14.9f %-14.9f %-10.3e   [MESA(12) vs RK4: %.2e]", Tτ, n, wb, wex, abs(wb - wex) / wex, abs(wm - wex) / wex))
        OK[] &= abs(wm - wex) / wex < 1e-9
    end
    # long implicit step: lands on the plateau regardless of w0 (as MESA), any dt
    wl = local_step(w0, a, U, sup, ℓ, 1e3 * τ, par)
    lp(@sprintf("    one BE step of 1e3 τ: %.12f (plateau %.12f)  rel %.2e  [BE is not exactly at the plateau after one finite step: error = (w−w0)/(dt R') ~ τ/dt = 1e-3 order]", wl, wp_ours, abs(wl - wp_ours) / wp_ours))
    OK[] &= abs(wl - wp_ours) <= 10e-3 * wp_ours
end
check("MESA (12) == exact ODE (RK4) and plateau == our steady root; BE converges (tables)", OK[])
# ---- decay to zero in finite time (MESA (13)) — independent derivation of the tan solution
lp("\nDecay (sup < 0): our Riccati  dw/dt = −(A_s a² |sup| + 2A_s U a w + A_d w²)/ℓ.  Exact: w = w_c tan(arctan(w0 ..)) generalised with U:")
for (U, sup, w0) in ((0.0, -0.3, 1.0), (0.5, -0.3, 1.0), (2.0, -0.05, 0.5))
    a = 1.0; ℓ = 1.0; par = TDCPar(α_t = 0.0)
    # complete the square: A_d (w + A_s U a/A_d)² + A_s a²|sup| − A_s² U² a²/A_d ; q ≡ A_s a²|sup|/A_d − (A_s U a/A_d)² (>0 for tan)
    As = par.A_s; Ad = par.A_d; c = As * U * a / Ad; q = As * a^2 * abs(sup) / Ad - c^2
    if q > 0
        sq = sqrt(q); tstar = (atan((w0 + c) / sq) - atan(c / sq)) * ℓ / (Ad * sq)
        wex(t) = sq * tan(atan((w0 + c) / sq) - Ad * sq * t / ℓ) - c
        # RK4 check of the formula, then BE hit time
        chk = rk4(w0, 0.5 * tstar, w -> tdc_local(w, a, U, sup, ℓ, par)); fm = wex(0.5 * tstar)
        lp(@sprintf("    U=%g sup=%g w0=%g: analytic t* = %.6f ℓ/a (formula vs RK4 at t*/2: %.2e)", U, sup, w0, tstar, abs(chk - fm) / fm))
        OK[] &= abs(chk - fm) / fm < 1e-9
        ξ0, ξ1, ξ2 = xi(a, U, sup, ℓ, As, Ad)
        for fr in (0.3, 0.9)
            wmx = mesa_step(sqrt(1.5) * w0, fr * tstar, ξ0, ξ1, ξ2) / sqrt(1.5); wrk = rk4(w0, fr * tstar, w -> tdc_local(w, a, U, sup, ℓ, par))
            lp(@sprintf("      MESA (13) at %.1f t*: %.9f  RK4 %.9f  rel %.2e", fr, wmx, wrk, abs(wmx - wrk) / wrk)); OK[] &= abs(wmx - wrk) / wrk < 1e-9
        end
        OK[] &= mesa_step(sqrt(1.5) * w0, 1.1 * tstar, ξ0, ξ1, ξ2) == 0.0
        for dt in (0.1, 0.01, 0.001)
            wb = w0; t = 0.0; thit = NaN
            while wb > 1e-20 && t < 3tstar
                wb = local_step(wb, a, U, sup, ℓ, dt, par); t += dt
            end
            lp(@sprintf("      BE dt=%-6g  w hits 0 (|w|<1e-20, Newton leaves a −1e-30 residue) at t = %.5f  (t*: %.5f, diff/dt = %.2f)  w = %g", dt, t, tstar, (t - tstar) / dt, wb))
            OK[] &= abs(wb) < 1e-14 && abs(t - tstar) < 4dt
        end
    else
        sc = sqrt(-q); y0 = w0 + c
        tstar = (atanh(sc / c) - atanh(sc / y0)) * ℓ / (Ad * sc)            # coth solution y = sc coth(u0 + A_d sc t/ℓ), w = y − c  [derived]
        wex2 = t -> sc / tanh(atanh(sc / y0) + Ad * sc * t / ℓ) - c
        chk = rk4(w0, 0.5 * tstar, w -> tdc_local(w, a, U, sup, ℓ, par)); fm = wex2(0.5 * tstar)
        lp(@sprintf("    U=%g sup=%g w0=%g: q = %.3f <= 0 (J²>0 in MESA's variables; both roots of the Riccati are negative): coth solution, t* = %.6f (formula vs RK4 at t*/2: %.2e)", U, sup, w0, q, tstar, abs(chk - fm) / fm))
        OK[] &= abs(chk - fm) / fm < 1e-9
        ξ0, ξ1, ξ2 = xi(a, U, sup, ℓ, As, Ad)
        wmx = mesa_step(sqrt(1.5) * w0, 0.5 * tstar, ξ0, ξ1, ξ2) / sqrt(1.5)
        lp(@sprintf("      MESA (12) literal at t*/2: %.9f vs exact %.9f (rel %.2e);  at 1.1 t* (12) gives %.5f < 0: MESA must clip to w=0 (not stated in the paper)", wmx, fm, abs(wmx - fm) / fm, mesa_step(sqrt(1.5) * w0, 1.1 * tstar, ξ0, ξ1, ξ2) / sqrt(1.5)))
        OK[] &= abs(wmx - fm) / fm < 1e-9
        for dt in (0.1, 0.01, 0.001)
            wb = w0; t = 0.0
            while wb > 1e-20 && t < 3tstar; wb = local_step(wb, a, U, sup, ℓ, dt, par); t += dt; end
            lp(@sprintf("      BE dt=%-6g  w reaches 0 at t = %.5f  (t*: %.5f, diff/dt = %.2f)  w = %g", dt, t, tstar, (t - tstar) / dt, wb))
            OK[] &= abs(wb) < 1e-14 && abs(t - tstar) < 4dt
        end
    end
end
check("decay: MESA (13) == exact ODE; finite-time zero reached by the FB step, first order in dt", OK[])
lp("MESA compare failures: ", NFAIL[]); close_log(); exit(NFAIL[] == 0 ? 0 : 1)
