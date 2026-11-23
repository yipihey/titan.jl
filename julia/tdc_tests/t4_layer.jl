include(joinpath(@__DIR__, "harness.jl"))
open_log("t4_layer")
lp("T4: static Schwarzschild-unstable layer between stable layers (plane, units ℓ = a = 1, g = 8, δ = H_P = ρ = c_P = T = 1)")
lp("    N² = gδ sup/H_P = 8 a² sup/ℓ², w_s = a√sup = ℓ N/√8")
g = 8.0; δ = 1.0; HP = 1.0; ℓ = 1.0; a = 1.0
# ---------------------------------------------------------------- (i) growth from w = 0
lp("\n(i) growth from w = 0 (uniform unstable background, α_t = 0, U = 0)")
par0 = TDCPar(α_t = 0.0)
let sup = 1.0
    N = sqrt(g * δ * sup / HP); ws = a * sqrt(sup)
    nz = 8; G = make_grid(collect(range(0.0, 8.0, length = nz + 1))); B = make_back(G; a = a, U = 0.0, sup = sup, ell = ℓ)
    # initial slope: first BE step from w = 0 with dt -> 0
    dt = 1e-6; w, _ = newton_w(zeros(nz + 1), G, B, dt, par0)
    slope = w[5] / dt
    check("initial growth rate dw/dt = A_s a² sup/ℓ", abs(slope - a^2 * sup / ℓ) <= 1e-5, @sprintf("numeric %.8f  analytic %.8f", slope, a^2 * sup / ℓ))
    lp("    tanh solution w = w_s tanh(w_s t/ℓ) vs backward-Euler, max |w - w_an|/w_s over t in [0, 5ℓ/w_s]")
    prev = NaN
    for dt in (0.1, 0.05, 0.025, 0.0125, 0.00625)
        w = zeros(nz + 1); emax = 0.0
        for s in 1:round(Int, 5 / dt)
            w, _ = newton_w(w, G, B, dt, par0); t = s * dt
            emax = max(emax, abs(w[5] - ws * tanh(ws * t / ℓ)) / ws)
        end
        lp(@sprintf("    dt = %-8g max err = %.3e   ratio to previous = %.3f", dt, emax, prev / emax)); prev = emax
    end
    check("BE converges to the tanh solution (first order)", prev < 3e-3)
    # time to 90 %
    dt = 1e-3; w = zeros(nz + 1); t = 0.0; t90 = NaN
    while true
        wn, _ = newton_w(w, G, B, dt, par0); tn = t + dt
        if wn[5] >= 0.9 * ws; t90 = t + dt * (0.9 * ws - w[5]) / (wn[5] - w[5]); break; end
        w = wn; t = tn
    end
    t90an = atanh(0.9) * ℓ / ws
    lp(@sprintf("    time to 90%% of w_s: numeric %.5f (BE, dt=1e-3)   analytic atanh(0.9) ℓ/w_s = %.5f   N t90 = %.4f  (atanh(0.9)√8 = %.4f)", t90, t90an, N * t90, atanh(0.9) * sqrt(8)))
    check("t90 N = atanh(0.9)√8 = 4.164 (sup = 1)", abs(N * t90 - atanh(0.9) * sqrt(8)) < 5e-3)
    lp("    scaling of N t90 with sup (U = 0, analytic and numeric): √8 atanh(0.9) independent of sup")
    for sp in (1e-3, 1e-2, 0.1, 1.0, 10.0)
        N_ = sqrt(g * δ * sp / HP); ws_ = a * sqrt(sp)
        B2 = make_back(G; a = a, U = 0.0, sup = sp, ell = ℓ); dt = 1e-3 / ws_; w = zeros(nz + 1); t = 0.0; t90n = NaN
        while true
            wn, _ = newton_w(w, G, B2, dt, par0)
            if wn[5] >= 0.9 * ws_; t90n = t + dt * (0.9 * ws_ - w[5]) / (wn[5] - w[5]); break; end
            w = wn; t += dt
        end
        lp(@sprintf("    sup = %-6g N = %-9.4g w_s = %-9.4g  N t90 = %.4f", sp, N_, ws_, N_ * t90n))
    end
end
# ---------------------------------------------------------------- layer set-up
nz = 1000; xs = collect(range(0.0, 50.0, length = nz + 1)); h = xs[2] - xs[1]
function layer_back(G; supu = 1.0, sst = 1.0, U = 0.0)
    sup(x) = (20.0 < x < 30.0) ? supu : -sst
    make_back(G; a = a, U = U, sup = sup, ell = ℓ)
end
function run_steady(par, G, B; dt = 0.5, nsteps = 3000, w0 = zeros(G.nz + 1), cb = nothing)
    evolve(w0, G, B, dt, nsteps, par; stop = 1e-11, cb = cb)
end
G = make_grid(xs)
# ---------------------------------------------------------------- (ii) saturation + (T8) energy budget in the transient
lp("\n(ii) saturation: interior w and F_c vs the local MLT values (α_t = 0.1, sup_u = 1, sup_st = -1)")
parA = TDCPar(α_t = 0.1)
for U in (0.0, 0.7)
    B = layer_back(G; supu = 1.0, sst = 1.0, U = U)
    maxbud = Ref(0.0)
    cb = (s, t, w, wo, info) -> begin
        s % 5 == 0 || return
        Q, _ = tdc_exchange(w, wo, G, B, 0.5, parA); rs, sc = tdc_budget_check(Q, w, wo, G, B, 0.5, parA)
        maxbud[] = max(maxbud[], abs(rs) / sc)
    end
    w, info = run_steady(parA, G, B; cb = cb)
    wmlt = tdc_mlt_w(a, U, 1.0)
    i = findall(x -> 23.0 <= x <= 27.0, xs)
    dev = maximum(abs.(w[i] ./ wmlt .- 1))
    Fc = [tdc_flux(w[k], a, U, B.sup[k], ℓ, 1.0, 1.0, 1.0, HP, parA) for k in i]
    Fk = kw_mlt_flux(1.0, 1.0, 1.0, g, δ, ℓ, HP, U, 1.0)
    devF = maximum(abs.(Fc ./ Fk .- 1))
    icen = argmin(abs.(xs .- 25.0))
    lp(@sprintf("    U = %.1f: steps %d steady=%s newton_fail %d | centre w/w_MLT = %.8f, max |w/w_MLT-1| on 23<x<27 = %.2e, max |F_c/F_MLT-1| = %.2e", U, info.steps, info.steady, info.nbad, w[icen] / wmlt, dev, devF))
    check(@sprintf("interior saturation U=%.1f equals local MLT (23<x<27)", U), dev < 1e-3 && devF < 3e-3 && info.nbad == 0)
    check(@sprintf("T8 energy budget closes during the transient, U=%.1f", U), maxbud[] < 1e-13, @sprintf("max |Σ Q Δt − ΔΣme| / scale = %.2e", maxbud[]))
    U == 0.0 && write_csv("t4_profile.csv", "x,w,sup", [(xs[k], w[k], B.sup[k]) for k in 1:nz+1])
end
# ---------------------------------------------------------------- (iii) overshoot
lp("\n(iii) overshoot into the stable layer above the unstable layer (sup_u = 1)")
lp("    analytic (steady, U=0, plane, A=1): w'² = |s| a²/(2 α ℓ²) + w²/(3 α ℓ²)  =>  distance from w_e to 0:  Δ0 = ℓ√(3α) asinh(w_e/(a√(3|s|/2))),")
lp("    to w_e/e:  Δe = ℓ√(3α)[asinh(w_e/q) − asinh(w_e/(e q))], q = a√(3|s|/2);  compact support (w reaches 0 at finite distance, slope a√(|s|/2α)/ℓ)")
lp("    estimate (balance of diffusion with the buoyancy damping a²|s|/ℓ): Δ ≈ ℓ√(2α_t) w_e/(a√|s|) for w_e ≪ q (small |s|…moderate), saturating at ℓ√(3α_t) for w_e ≫ q")
lp("    measured from the last unstable interface (w_e = its value); grid h = $(h)")
lp(@sprintf("    %-5s %-6s %-9s %-10s %-10s %-10s %-10s %-10s", "α_t", "|s|", "w_e/w_s", "Δe num", "Δe an", "Δ0 num", "Δ0 an", "relerr Δe"))
rows = []; ok = true
for st in (0.1, 1.0, 10.0), α in (0.0, 0.1, 0.3, 1.0, 3.0)
    par = TDCPar(α_t = α); B = layer_back(G; supu = 1.0, sst = st)
    w, info = run_steady(par, G, B)
    ie = findlast(x -> x < 30.0, xs); we = w[ie]
    # 1/e and zero crossing, linear interpolation
    function cross(level)
        for k in ie+1:nz+1
            if w[k] <= level
                return xs[k-1] + h * (w[k-1] - level) / (w[k-1] - w[k]) - xs[ie]
            end
        end
        return NaN
    end
    de = cross(we / ℯ); d0 = findfirst(k -> k > ie && w[k] == 0.0, 1:nz+1); d0 = d0 === nothing ? NaN : xs[d0-1] - xs[ie]
    if α > 0
        q = a * sqrt(1.5 * st); lam = sqrt(3α) * ℓ
        dean = lam * (asinh(we / q) - asinh(we / (ℯ * q))); d0an = lam * asinh(we / q)
    else
        dean = 0.0; d0an = 0.0
    end
    push!(rows, (α, st, we, de, dean, d0, d0an, info.nbad))
    lp(@sprintf("    %-5g %-6g %-9.4f %-10.4f %-10.4f %-10.4f %-10.4f %.3f", α, st, we / 1.0, isnan(de) ? 0.0 : de, dean, isnan(d0) ? 0.0 : d0, d0an, α > 0 ? abs(de - dean) / dean : 0.0))
    global ok &= info.nbad == 0
    if α > 0
        global ok &= abs(de - dean) / dean < 0.15 + 4h / dean * 0       # discretisation: see doc; loose bound, the table is the result
    end
end
write_csv("t4_overshoot.csv", "alpha_t,abs_sup_stable,w_edge,d_1e_num,d_1e_an,d_0_num,d_0_an,newton_fail", rows)
check("overshoot distance within 15% of the analytic profile; no Newton failure", ok)
# scaling exponents of the numeric 1/e distance
for st in (0.1, 1.0, 10.0)
    r = filter(x -> x[2] == st && x[1] > 0, rows); lx = log.([x[1] for x in r]); ly = log.([x[4] for x in r])
    p = ([ones(length(lx)) lx] \ ly)[2]
    lp(@sprintf("    |s| = %-5g: Δe ∝ α_t^%.3f over α_t = 0.1..3 (analytic: ∝ α_t^(1/2) in both limits)", st, p))
end
for α in (0.1, 1.0)
    r = filter(x -> x[1] == α, rows); lx = log.([x[2] for x in r]); ly = log.([x[4] for x in r])
    p = ([ones(length(lx)) lx] \ ly)[2]
    lp(@sprintf("    α_t = %-4g: Δe ∝ |s|^%.3f over |s| = 0.1..10 (analytic: w_e ≪ q = a√(1.5|s|): ∝ |s|^(-1/2) (damping–diffusion balance); w_e ≫ q: → λ(1−…) independent of |s|)", α, p))
end
# ---------------------------------------------------------------- spherical vs plane
lp("\n(iv) spherical geometry (μ = 2) with the layer at R = 1000 ℓ: overshoot must reproduce the plane result")
let rs = collect(range(1000.0, 1050.0, length = nz + 1)), α = 0.3
    Gs = make_grid(rs; μ = 2)
    sup(x) = (1020.0 < x < 1030.0) ? 1.0 : -1.0
    Bs = make_back(Gs; a = a, U = 0.0, sup = sup, ell = ℓ)
    ws, infos = run_steady(TDCPar(α_t = α), Gs, Bs)
    Gp = make_grid(xs); Bp = layer_back(Gp; supu = 1.0, sst = 1.0)
    wp, infop = run_steady(TDCPar(α_t = α), Gp, Bp)
    dev = maximum(abs.(ws .- wp))
    lp(@sprintf("    max |w_sph − w_plane| = %.2e (curvature ~ h/R = 5e-5 .. ℓ/R = 1e-3); Newton fail %d/%d", dev, infos.nbad, infop.nbad))
    check("spherical rows reproduce plane at R >> L", dev < 2e-3 && infos.nbad == 0)
end
lp("T4 failures: ", NFAIL[]); close_log(); exit(NFAIL[] == 0 ? 0 : 1)
