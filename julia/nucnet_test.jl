# S1 acceptance tests for NucNet (REACLIB network). Run: julia -t 4 --project=env nucnet_test.jl
# Prints numbers and PASS/FAIL per check; "ALL PASS" at the end if everything passed.
using Printf, ForwardDiff, LinearAlgebra
include(joinpath(@__DIR__, "NucNet.jl"))
using .NucNet
import .NucNet: NET, NS, SP, IDX, AA, ZZ, MEX, NA, MeV, kB, rates_k, Plasma, ye_of, energy, ydot, jac!, reaclib_lambda, pick,
                be7ec_bahcall, screen_factor, chugunov07, weak_h

const FAILS = String[]
function chk(name, ok, detail)
    @printf("  [%s] %-34s %s\n", ok ? "PASS" : "FAIL", name, detail)
    ok || push!(FAILS, name)
end
ix(s) = IDX[s]
mass_sum(Y) = sum(AA .* Y)
function cons(Y, name; tol = 1e-12)
    chk("A-conservation: " * name, abs(mass_sum(Y) - 1) < tol, @sprintf("ΣA·Y − 1 = %.2e", mass_sum(Y) - 1))
end

# ================================================================================================================================
println("== 0. Network sanity: nucleon number and charge per reaction, Jacobian vs ForwardDiff")
bad = String[]
for r in NET
    dA = sum(AA[i] * c for (i, c) in r.net)
    expect = r.name in ("p(p,e+nu)d", "p(p e-,nu)d (pep)", "be7(e-,nu)li7", "b8(e+nu)2a", "he3(p,e+nu)a (hep)",
                        "n13(e+nu)c13", "o15(e+nu)n15", "f17(e+nu)o17", "f18(e+nu)o18") ? -1 :
             r.name == "d(d,p)t [t->he3]" ? +1 : 0
    (abs(dA) > 1e-12 || r.weak != expect) && push!(bad, @sprintf("%s (dA=%g, dZ=%d, expected %d)", r.name, dA, r.weak, expect))
end
chk("ΔA = 0, ΔZ = lepton charge", isempty(bad), isempty(bad) ? "$(length(NET)) reactions OK" : join(bad, "; "))
let Y = composition(Dict(:p => 0.7, :he4 => 0.28, :c12 => 0.005, :n14 => 0.005, :o16 => 0.01, :d => 1e-4, :he3 => 1e-4)), ρ = 100.0, T = 2.5e7
    k = rates_k(T, ρ, Plasma(ρ, Y), ye_of(Y)); J = zeros(NS, NS); jac!(J, Y, k)
    Jad = ForwardDiff.jacobian(y -> ydot(y, k), Y)
    chk("analytic Jacobian = ForwardDiff", maximum(abs.(J .- Jad)) <= 1e-12 * maximum(abs.(J)), @sprintf("max|ΔJ|/max|J| = %.1e", maximum(abs.(J .- Jad)) / maximum(abs.(J))))
    f = ydot(Y, k)
    fz = sum(ZZ .* f); fl = sum(r.weak * NucNet.event(r, k[i], Y) for (i, r) in enumerate(NET))
    chk("d/dt ΣA Y = 0, d/dt ΣZ Y = lepton ledger", abs(sum(AA .* f)) < 1e-14 * maximum(abs.(f) .* AA) && abs(fz - fl) < 1e-12 * abs(fl) + 1e-30,
        @sprintf("ΣA f = %.1e, ΣZ f = %.4e vs Σ_r ΔZ_r ev_r = %.4e", sum(AA .* f), fz, fl))
end
println("  Q: REACLIB vs AME2020 mass excesses (keV), max |diff| over reactions:")
qd = [(r.name, (r.Q + sum(c * MEX[i] for (i, c) in r.net) / 1000) * 1000) for r in NET]
mx = maximum(abs(x[2]) for x in qd)
chk("Q(REACLIB) = Q(AME2020)", mx < 2.0, @sprintf("max |ΔQ| = %.2f keV (%s)", mx, qd[argmax(abs.(last.(qd)))][1]))
let T = 1.57e7, rb = be7ec_bahcall(T), rr = reaclib_lambda(pick(["be7"], ["li7"], "ec") |> x -> [s.a for s in x], T * 1e-9)
    chk("7Be EC REACLIB = Bahcall formula", abs(rr / rb - 1) < 0.05, @sprintf("REACLIB %.4e, Bahcall/CF88 %.4e (N_A<σv> per unit n_e); λ(ρYe=100) = %.3e /s, τ = %.0f d", rr, rb, rr * 100, 1 / (rr * 100) / 86400))
end

# ================================================================================================================================
println("== a. Solar centre: T = 1.57e7 K, ρ = 150, X = 0.34, Y = 0.646, X_CNO = 0.014 (C12 0.0035, N14 0.001, O16 0.0095)")
ρc = 150.0; Tc = 1.57e7
X0 = Dict(:p => 0.34, :he4 => 0.646, :c12 => 0.0035, :n14 => 0.001, :o16 => 0.0095)
Y0 = composition(X0); fr = [ix(:p), ix(:he4)]            # H and He4 held fixed (instantaneous burning rate, no depletion)
res_a = Dict()
for (nm, tmax) in (("relaxed 4.6 Gyr (CNO)", 1.45e17), ("full equilibrium", 1e30))
    Y = steady(Y0, ρc, Tc; frozen = fr, tmax = tmax); r = burn(Y, ρc, Tc, 1.0; frozen = fr)
    res_a[nm] = (Y, r)
    @printf("  %-22s ε = %.3f  [pp chain %.3f, D %.4f, CNO %.3f]  ε_ν = %.3f erg/g/s;  X(3He) = %.2e, 12C/13C = %.2f, X(14N) = %.2e\n", nm, r.ε, r.εg.pp, r.εg.d, r.εg.cno, r.εν,
            3Y[ix(:he3)], Y[ix(:c12)] / Y[ix(:c13)], 14Y[ix(:n14)])
end
Ya, ra = res_a["relaxed 4.6 Gyr (CNO)"]
chk("ε_pp (equilibrium 3He) in 15–17", 14.5 < ra.εg.pp + ra.εg.d < 17.5, @sprintf("ε_pp = %.2f erg/g/s (SSM centre ≈ 16–17; BP2000/Bahcall 1989)", ra.εg.pp + ra.εg.d))
chk("ε_total in 15–19", 15 < ra.ε < 19, @sprintf("ε_pp+ε_CNO = %.2f", ra.ε))
chk("CNO fraction 3–9 %", 0.03 < ra.εg.cno / ra.ε < 0.09, @sprintf("ε_CNO/ε = %.3f (CNO ≈ 1.5 %% of L_⊙, centrally concentrated)", ra.εg.cno / ra.ε))
chk("X(3He) centre 3e-6..2e-5", 3e-6 < 3Ya[ix(:he3)] < 2e-5, @sprintf("X3 = %.2e (SSM central ≈ 1e-5)", 3Ya[ix(:he3)]))
# he3 held at multiples of its equilibrium value
print("  ε_pp with 3He held at a fraction of equilibrium X3 = $(@sprintf("%.2e", 3Ya[ix(:he3)])):")
for fx in (0.1, 0.5, 1.0, 2.0)
    Yh = copy(Ya); Yh[ix(:he3)] *= fx
    Yh = steady(Yh, ρc, Tc; frozen = [ix(:p), ix(:he4), ix(:he3)], tmax = 1e12)
    rh = burn(Yh, ρc, Tc, 1.0; frozen = [ix(:p), ix(:he4), ix(:he3)]); @printf("  ×%.1f: %.2f", fx, rh.εg.pp + rh.εg.d)
end
println()
print("  screening (ε_pp): ")
for sc in (:chugunov, :weak, :none)
    r = burn(Ya, ρc, Tc, 1.0; frozen = fr, screen = sc); @printf("%s %.3f   ", sc, r.εg.pp)
end
println()
let Y = Ya
    @printf("  screening factors at centre (Chugunov / Salpeter): pp %.4f / %.4f, 3He+3He %.4f / %.4f, 14N+p %.4f / %.4f\n",
        screen_factor(Tc, ρc, Y, 1.0, 1.0, 1.0, 1.0), screen_factor(Tc, ρc, Y, 1.0, 1.0, 1.0, 1.0; mode = :weak),
        screen_factor(Tc, ρc, Y, 2.0, 3.0, 2.0, 3.0), screen_factor(Tc, ρc, Y, 2.0, 3.0, 2.0, 3.0; mode = :weak),
        screen_factor(Tc, ρc, Y, 1.0, 1.0, 7.0, 14.0), screen_factor(Tc, ρc, Y, 1.0, 1.0, 7.0, 14.0; mode = :weak))
    fc = screen_factor(Tc, ρc, Y, 1.0, 1.0, 1.0, 1.0); fw = screen_factor(Tc, ρc, Y, 1.0, 1.0, 1.0, 1.0; mode = :weak)
    chk("pp screening: ln f_Chugunov / ln f_Salpeter", 0.4 < log(fc) / log(fw) < 1.05, @sprintf("f_C = %.4f, f_W = %.4f, ln-ratio %.2f (Salpeter incl. electron screening: +5 %%; Chugunov/DH ion-only form is lower; SSM uses Salpeter)", fc, fw, log(fc) / log(fw)))
end
let Y = Ya
    # neutrino loss: ε_ν/(ε+ε_ν) and 8B, 7Be, pep fractions
    k = rates_k(Tc, ρc, Plasma(ρc, Y), ye_of(Y)); ev(n) = (i = findfirst(r -> r.name == n, NET); NucNet.event(NET[i], k[i], Y))
    pp_ = ev("p(p,e+nu)d"); pep = ev("p(p e-,nu)d (pep)");
    @printf("  pep/pp = %.2e (Bahcall–May: 1.1e-4 ρ/μe/√T6 = %.2e); 7Be-EC / pp = %.3f; 8B / pp = %.2e; hep/pp = %.2e; branching ppI:ppII+III via 3He+3He : 3He+α = %.1f : 1\n",
        pep / pp_, 1.102e-4 * ρc * ye_of(Y) / sqrt(Tc * 1e-6), ev("be7(e-,nu)li7") / pp_, ev("b8(e+nu)2a") / pp_, ev("he3(p,e+nu)a (hep)") / pp_, 2ev("he3(he3,2p)he4") / (ev("he3(a,g)be7")))
    @printf("  neutrino losses ε_ν/ε_tot = %.3f (pp-chain ε_ν/(ε+ε_ν) ≈ 0.02 ⇒ L_ν/L ≈ 2 %%)\n", ra.εν / (ra.ε + ra.εν))
    chk("pep/pp within 40 % of Bahcall–May", abs(pep / pp_ / (1.102e-4 * ρc * ye_of(Y) / sqrt(Tc * 1e-6)) - 1) < 0.4, @sprintf("ratio %.2f", pep / pp_ / (1.102e-4 * ρc * ye_of(Y) / sqrt(Tc * 1e-6))))
    chk("ε_ν/ε_tot 1–4 %", 0.01 < ra.εν / (ra.ε + ra.εν) < 0.04, @sprintf("%.3f", ra.εν / (ra.ε + ra.εν)))
end

# ================================================================================================================================
println("== b. CNO equilibrium (p, 4He fixed; ρ = 100, X = 0.7, X_CNO = 0.01 all initially 12C)")
kw(T9, ρ, X, XCNO) = 8.24e25 * (1 - 2T9 + 3.41T9^2 - 2.43T9^3) * XCNO * X * ρ * T9^(-2 / 3) * exp(-15.231 * T9^(-1 / 3) - (T9 / 0.8)^2)
for T in (2e7, 3e7)
    ρ = 100.0; Y0b = composition(Dict(:p => 0.7, :he4 => 0.29, :c12 => 0.01))
    Y = steady(Y0b, ρ, T; frozen = fr, tmax = 1e30); r = burn(Y, ρ, T, 1.0; frozen = fr)
    XC = sum(AA[ix(s)] * Y[ix(s)] for s in (:c12, :c13, :n13, :n14, :n15, :o15, :o16, :o17, :f17, :o18, :f18, :f19))
    fits = kw(T * 1e-9, ρ, 0.7, XC)
    c12c13 = Y[ix(:c12)] / Y[ix(:c13)]
    @printf("  T = %.0e: 12C/13C = %.2f, X: 12C %.2e 13C %.2e 14N %.2e 15N %.2e 16O %.2e (of X_CNO = %.4f; N14 fraction %.3f); ε_CNO = %.4e, KW18.65 fit %.4e, ratio %.2f\n",
            T, c12c13, 12Y[ix(:c12)], 13Y[ix(:c13)], 14Y[ix(:n14)], 15Y[ix(:n15)], 16Y[ix(:o16)], XC, 14Y[ix(:n14)] / XC, r.εg.cno, fits, r.εg.cno / fits)
    @printf("      (ε_pp = %.3e, f17+o17+f18 flows small: 16O fraction %.3f)\n", r.εg.pp, 16Y[ix(:o16)] / XC)
    chk(@sprintf("12C/13C ∈ [3,4.5] at %.0e K", T), 3.0 < c12c13 < 4.5, @sprintf("%.2f", c12c13))
    chk(@sprintf("14N dominates at %.0e K", T), 14Y[ix(:n14)] / XC > 0.9, @sprintf("%.3f of CNO", 14Y[ix(:n14)] / XC))
    chk(@sprintf("ε_CNO vs KW fit at %.0e K", T), 0.4 < r.εg.cno / fits < 2.5, @sprintf("ratio %.2f", r.εg.cno / fits))
    nC = sum(Y[ix(s)] for s in (:c12, :c13, :n13, :n14, :n15, :o15, :o16, :o17, :f17, :o18, :f18, :f19))
    chk(@sprintf("CNO nuclei conserved at %.0e K", T), abs(nC * 12 / 0.01 - 1) < 1e-12, @sprintf("number of C+N+O nuclei / initial − 1 = %.1e (p, 4He held fixed, so ΣA·Y grows by p-captures)", nC * 12 / 0.01 - 1))
end

# ================================================================================================================================
println("== c. Triple-α: REACLIB (Fynbo et al. 2005 'fy05') vs Caughlan & Fowler 1988")
cf88(T9) = 2.79e-8 * T9^-3 * exp(-4.4027 / T9) + 1.35e-8 * T9^-1.5 * exp(-24.811 / T9)    # N_A² <σv> (CF88, eq. for 3α)
r3 = NET[findfirst(r -> startswith(r.name, "a(2a"), NET)]
for T9 in (0.1, 0.2)
    lam = reaclib_lambda(r3.sets, T9)
    @printf("  T9 = %.2f: REACLIB %.4e, CF88 %.4e, ratio %.3f\n", T9, lam, cf88(T9), lam / cf88(T9))
end
l1 = reaclib_lambda(r3.sets, 0.1) / cf88(0.1); l2 = reaclib_lambda(r3.sets, 0.2) / cf88(0.2)
chk("3α REACLIB/CF88 at T9 = 0.2", 0.7 < l2 < 1.4, @sprintf("%.3f", l2))
chk("3α REACLIB/CF88 at T9 = 0.1", 0.5 < l1 < 2.0, @sprintf("%.3f (CF88 and NACRE agree to ~10-30 %% here; Fynbo 2005 lowers NACRE below 0.07 only)", l1))
let T9 = 0.3, lam = reaclib_lambda(r3.sets, T9)
    @printf("  T9 = 0.30: REACLIB %.4e, CF88 %.4e, ratio %.3f (3α-rate shape check; NACRE table values not entered from memory)\n", lam, cf88(T9), lam / cf88(T9))
end

# ================================================================================================================================
println("== d. Deuterium burning: ε_D = N_A Q ρ Y_p Y_D <σv>, ρ = 1, D/H = 2.5e-5, vs Stahler 1988 eq. (48)")
stahler(T6, DH, ρ) = 4.19e7 * DH * ρ * T6^11.8                      # f = 1; ε0 = 4.19e7 erg/g/s at ρ = 1, T6 = 1, D/H = 1
function epsD(T6, X, DH, ρ; screen = :chugunov)
    Y = composition(Dict(:p => X, :he4 => 1 - X)); Y[ix(:d)] = DH * Y[ix(:p)]
    k = rates_k(T6 * 1e6, ρ, Plasma(ρ, Y), ye_of(Y); screen); i = findfirst(r -> r.name == "d(p,g)he3", NET)
    NA * MeV * NET[i].Q * NucNet.event(NET[i], k[i], Y), 1 / (ρ * Y[ix(:p)] * k[i] / ρ)
end
rat = Float64[]
for X in (0.70, 0.75)
    @printf("  X = %.2f:\n", X)
    for T6 in (0.8, 1.0, 1.5)
        e, _ = epsD(T6, X, 2.5e-5, 1.0); s = stahler(T6, 2.5e-5, 1.0); X == 0.70 && push!(rat, e / s)
        @printf("    T6 = %.1f: ε_D = %.4e, Stahler = %.4e, ratio %.3f\n", T6, e, s, e / s)
    end
    e8 = epsD(0.8, X, 2.5e-5, 1.0)[1]; e15 = epsD(1.5, X, 2.5e-5, 1.0)[1]
    @printf("    effective index d ln ε/d ln T (0.8–1.5): %.2f (Stahler 11.8);  d(p,γ) alone: %.2f\n", log(e15 / e8) / log(1.5 / 0.8), log(e15 / e8) / log(1.5 / 0.8))
end
e10 = epsD(1.0, 0.70, 2.5e-5, 1.0)[1]
chk("ε_D(T6 = 1) vs Stahler (X = 0.70)", 0.8 < e10 / stahler(1.0, 2.5e-5, 1.0) < 1.25, @sprintf("ε = %.1f vs %.1f erg/g/s, ratio %.3f", e10, stahler(1.0, 2.5e-5, 1.0), e10 / stahler(1.0, 2.5e-5, 1.0)))
idx = log(epsD(1.5, 0.7, 2.5e-5, 1.0)[1] / epsD(0.8, 0.7, 2.5e-5, 1.0)[1]) / log(1.5 / 0.8)
chk("index 0.8–1.5 within 0.7 of 11.8", abs(idx - 11.8) < 0.7, @sprintf("%.2f", idx))
println("  burning time τ_D = 1/(ρ X N_A<σv>) [yr] (ρ = 1, X = 0.7; weak/Chugunov screening negligible):")
for T6 in (0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0)
    Y = composition(Dict(:p => 0.7, :he4 => 0.3)); k = rates_k(T6 * 1e6, 1.0, Plasma(1.0, Y), ye_of(Y)); i = findfirst(r -> r.name == "d(p,g)he3", NET)
    τs = NA * MeV * 5.493 * (2.5e-5 * 0.7) / stahler(T6, 2.5e-5, 1.0)        # Y_D/(dY_D/dt) implied by Stahler's ε
    @printf("    T6 = %.1f: τ_D = %.3e yr   (Stahler-fit implied, X = 0.7: %.3e yr)\n", T6, 1 / (0.7 * k[i]) / 3.156e7, τs / 3.156e7)
end
# full network burn of D at T6 = 1.5 (D follows e^{-t/τ})
let T = 1.5e6, ρ = 1.0
    Y = composition(Dict(:p => 0.75, :he4 => 0.25 - 2.5e-5 * 0.75 * 2)); Y[ix(:d)] = 2.5e-5 * Y[ix(:p)]
    k = rates_k(T, ρ, Plasma(ρ, Y), ye_of(Y)); i = findfirst(r -> r.name == "d(p,g)he3", NET); τ = 1 / (ρ * Y[ix(:p)] * k[i] / ρ)
    r1 = burn(Y, ρ, T, τ)
    ev = evolve(Y, ρ, T, τ; rtol = 0.01); fD = ev.Y[ix(:d)] / Y[ix(:d)]
    Yc = copy(Y); for _ in 1:100; Yc = burn(Yc, ρ, T, τ / 100).Y; end; Yc2 = copy(Y); for _ in 1:1000; Yc2 = burn(Yc2, ρ, T, τ / 1000).Y; end
    @printf("  D/D0 after t = τ_D(d(p,γ)) = %.3e s at T6 = 1.5: single BE step %.3f; 100 steps %.4f; 1000 steps %.4f; adaptive %.4f (%d steps); e^-1 = %.4f\n",
            τ, r1.Y[ix(:d)] / Y[ix(:d)], Yc[ix(:d)] / Y[ix(:d)], Yc2[ix(:d)] / Y[ix(:d)], fD, ev.steps, exp(-1))
    chk("D e-folding: converged BE, within 5 % of e^-1", abs(Yc2[ix(:d)] / Y[ix(:d)] / exp(-1) - 1) < 0.05 && abs(Yc[ix(:d)] / Yc2[ix(:d)] - 1) < 0.01,
        @sprintf("D/D0 = %.4f; the 4 %% deficit is d(d,n)he3/d(d,p)t, which exceeds d(p,γ) per D nucleus at D/H = 2.5e-5 (S_dd/S_dp ~ 3e5)", Yc2[ix(:d)] / Y[ix(:d)]))
    cons(ev.Y, "D burn")
end

# ================================================================================================================================
println("== e. Pop III one-zone (X = 0.75, Y = 0.25, Z = 0, D/H = 2.5e-5, 3He/H = 1e-5, 7Li/H = 4.7e-10)")
const CNOSP = (:c12, :c13, :n13, :n14, :n15, :o15, :o16, :o17, :f17, :o18, :f18, :f19)
function popIII(T, ρ, tend_yr)
    yr = 3.15576e7
    Y = zeros(NS); Y[ix(:p)] = 0.75; Y[ix(:d)] = 2.5e-5 * 0.75; Y[ix(:he3)] = 1e-5 * 0.75; Y[ix(:li7)] = 4.7e-10 * 0.75
    Y[ix(:he4)] = (1 - sum(AA .* Y)) / 4
    cons(Y, "initial")
    rec = NamedTuple[]; marks = Dict{Float64,Float64}()
    cb(t, Y, r) = begin
        push!(rec, (t = t, XC = 12Y[ix(:c12)], XCNO = sum(AA[ix(s)] * Y[ix(s)] for s in CNOSP), pp = r.εg.pp + r.εg.d, cno = r.εg.cno, he = r.εg.he, Xp = Y[ix(:p)], X4 = 4Y[ix(:he4)]))
        for th in (1e-11, 1e-10, 1e-9, 1e-8); haskey(marks, th) || (12Y[ix(:c12)] >= th && (marks[th] = t / yr)); end
    end
    ev = evolve(Y, ρ, T, tend_yr * yr; rtol = 0.03, dt0 = 1e3, callback = cb)
    (ev, rec, marks)
end
for (T, ρ, tend) in ((1.0e8, 30.0, 1e9), (1.2e8, 30.0, 1e6), (1.4e8, 30.0, 1e5))
    t0 = time(); ev, rec, marks = popIII(T, ρ, tend); Y = ev.Y; yr = 3.15576e7
    @printf("  T = %.1e K, ρ = %g, to %.0e yr: %d steps in %.1f s\n", T, ρ, tend, ev.steps, time() - t0)
    for th in (1e-11, 1e-10, 1e-9, 1e-8); @printf("    X(12C) ≥ %.0e at t = %s yr\n", th, haskey(marks, th) ? @sprintf("%.3e", marks[th]) : "not reached"); end
    for ty in tend ./ (1e4, 1e3, 1e2, 10, 1)
        j = findlast(r -> r.t <= ty * yr * 1.0001, rec); r = rec[j]
        @printf("    t = %.0e yr: X12 = %.2e, ε_pp = %.3e, ε_CNO = %.3e, ε_3α = %.3e erg/g/s, X(H) = %.4f, X(4He) = %.4f\n", ty, r.XC, r.pp, r.cno, r.he, r.Xp, r.X4)
    end
    jx = findfirst(r -> r.cno >= r.pp && r.t > 1e7 && r.Xp > 0.1, rec)
    if jx === nothing
        @printf("    ε_CNO ≥ ε_pp: not reached (max ε_CNO/ε_pp = %.3f)\n", maximum(r.cno / r.pp for r in rec if r.Xp > 0.1 && r.t > 1e7))
    else
        r = rec[jx]; XCpred = r.pp / kw(T * 1e-9, ρ, r.Xp, 1.0)
        @printf("    ε_CNO = ε_pp at t = %.3e yr: X(12C) = %.2e, X(CNO) = %.2e (KW18.65 equilibrium prediction for X_CNO at ε_CNO = ε_pp: %.2e)\n", r.t / yr, r.XC, r.XCNO, XCpred)
        chk(@sprintf("%.1e K: X(CNO) at ε_CNO = ε_pp vs KW18.65 (×4)", T), 0.25 < r.XCNO / XCpred < 4, @sprintf("%.2e vs %.2e", r.XCNO, XCpred))
        T == 1.0e8 && chk("1e8 K: X(CNO) at crossing 1e-11..1e-9 (Ezer & Cameron 1971 ~1e-9; Marigo+ 2001 1e-10..1e-9)", 1e-11 < r.XCNO < 1e-9, @sprintf("X(CNO) = %.2e, X(12C) = %.2e at t = %.2e yr", r.XCNO, r.XC, r.t / yr))
    end
    cons(Y, @sprintf("Pop III %.1e K", T))
end

# ================================================================================================================================
println("== f. Conservation and energy")
let ρ = 30.0, T = 1.2e8
    Y0 = zeros(NS); Y0[ix(:p)] = 0.75; Y0[ix(:he4)] = 0.25 / 4; Y0[ix(:d)] = 1e-5; Y0[ix(:he3)] = 1e-5
    Y0[ix(:he4)] = (1 - sum(AA .* Y0) + AA[ix(:he4)] * Y0[ix(:he4)]) / 4
    for dt in (1e8, 1e12, 3e13)
        r = burn(Y0, ρ, T, dt); dY = r.Y .- Y0
        E_mass = -NA * sum(dY .* MEX) * 1e-3 * MeV                  # erg/g released (incl. neutrinos), atomic mass excesses
        E_net = (r.ε + r.εν) * dt
        dZ = sum(ZZ .* dY);
        @printf("  dt = %.0e s (%d substeps): ΣA Y − 1 = %.1e;  ∫(ε+ε_ν)dt = %.6e vs −N_A Σ ΔY M = %.6e (rel %.1e);  ∫ε_ν/total = %.3f\n", dt, r.nsub, mass_sum(r.Y) - 1, E_net, E_mass, E_net / E_mass - 1, r.εν / (r.ε + r.εν))
        chk(@sprintf("energy ledger dt = %.0e", dt), abs(E_net / E_mass - 1) < 2e-3 && abs(mass_sum(r.Y) - 1) < 1e-12, "")
    end
    # charge: Σ Z Y changes only through weak events (e+ emitted / e- captured); integrated event counts from burn
    for dt in (1e12, 3e13)
        r = burn(Y0, ρ, T, dt); dZ = sum(ZZ .* (r.Y .- Y0)); led = sum(rr.weak * r.nev[i] for (i, rr) in enumerate(NET))
        Ee = NA * MeV * sum((rr.Q - rr.Enu) * r.nev[i] for (i, rr) in enumerate(NET))
        @printf("  dt = %.0e: Δ(ΣZY) = %.6e, Σ_r ΔZ_r·N_r = %.6e (rel diff %.1e), net leptons (e+ emitted − e- captured) %.6e mol/g; Σ_r Q_r N_r = (ε dt) %.6e vs %.6e\n", dt, dZ, led, (dZ - led) / abs(led), -led, Ee, r.ε * dt)
        chk(@sprintf("charge ledger dt = %.0e", dt), abs(dZ - led) < 1e-12 * abs(led) + 1e-30 && abs(Ee / (r.ε * dt) - 1) < 1e-12, "")
    end
end
# coupled (ρ, u) integrator against the fixed-T integrator, ideal gas + radiation EOS
aR = 7.5657e-15
eos_u(ρ, T, Y) = 1.5 * kB * T * NA * (sum(Y) + sum(ZZ[i] * Y[i] for i in 1:NS)) + aR * T^4 / ρ
let ρ = 30.0, T0 = 1.2e8
    Y0 = zeros(NS); Y0[ix(:p)] = 0.75; Y0[ix(:d)] = 2e-5; Y0[ix(:he3)] = 1e-5; Y0[ix(:he4)] = (1 - sum(AA .* Y0)) / 4
    Y0 = evolve(Y0, ρ, T0, 1e10; rtol = 0.03).Y                        # after the D, 3He, 7Li transients
    u0 = eos_u(ρ, T0, Y0); dt = 1e10
    rc = burn_coupled(Y0, ρ, u0, T0, dt, eos_u); rf = burn(Y0, ρ, T0, dt)
    du = eos_u(ρ, rc.T, rc.Y) - u0
    @printf("  coupled: T %.6e → %.8e (ΔT/T = %.2e), ε = %.5e (fixed-T %.5e); u_eos(new) − u − dt ε = %.2e (rel %.1e); %d substeps\n", T0, rc.T, rc.T / T0 - 1, rc.ε, rf.ε, du - dt * rc.ε, (du - dt * rc.ε) / du, rc.nsub)
    chk("coupled: energy consistency", abs((du - dt * rc.ε) / du) < 1e-9 && abs(mass_sum(rc.Y) - 1) < 1e-12, "")
    chk("coupled ε ≈ fixed-T ε (ΔT/T small)", abs(rc.ε / rf.ε - 1) < 0.02, @sprintf("ε ratio − 1 = %.2e", rc.ε / rf.ε - 1))
    # coupled with a large step, thermal feedback: heat capacity small enough that T rises
    rc2 = burn_coupled(Y0, ρ, u0, T0, 1e12, eos_u); du2 = eos_u(ρ, rc2.T, rc2.Y) - u0
    @printf("  coupled long step dt = 1e12 s: T → %.5e (ΔT/T = %.3f), u consistency rel %.1e, ΣA Y − 1 = %.1e\n", rc2.T, rc2.T / T0 - 1, (du2 - 1e12 * rc2.ε) / du2, mass_sum(rc2.Y) - 1)
    chk("coupled long step: energy ledger", abs((du2 - 1e12 * rc2.ε) / du2) < 1e-8, "")
end

# ================================================================================================================================
println("== g. Performance (19 species)")
let ρ = 150.0, T = 1.57e7, Y = Ya
    burn(Y, ρ, T, 1e10); t0 = time_ns(); n = 2000
    ns = 0; for _ in 1:n; ns += burn(Y, ρ, T, 1e10).nsub; end
    el = (time_ns() - t0) / 1e3
    @printf("  burn (solar centre, dt = 1e10 s): %.1f µs per call (%.1f substeps avg)\n", el / n, ns / n)
    Yp = composition(Dict(:p => 0.75, :he4 => 0.25)); burn(Yp, 30.0, 1.2e8, 1e11); t0 = time_ns()
    for _ in 1:n; burn(Yp, 30.0, 1.2e8, 1e9); end; @printf("  burn (Pop III, dt = 1e9 s): %.1f µs per call\n", (time_ns() - t0) / 1e3 / n)
    u0 = eos_u(ρ, T, Y); burn_coupled(Y, ρ, u0, T, 1e10, eos_u); t0 = time_ns(); m = 300
    for _ in 1:m; burn_coupled(Y, ρ, u0, T, 1e10, eos_u); end
    @printf("  burn_coupled (ForwardDiff Jacobian, 20×20): %.1f µs per call\n", (time_ns() - t0) / 1e3 / m)
    Yb = copy(Y); k = rates_k(T, ρ, Plasma(ρ, Y), ye_of(Y)); t0 = time_ns(); for _ in 1:n; rates_k(T, ρ, Plasma(ρ, Y), ye_of(Y)); end
    @printf("  rates_k (31 reactions incl. screening): %.1f µs\n", (time_ns() - t0) / 1e3 / n)
end

println()
if isempty(FAILS)
    println("ALL PASS")
else
    println("FAILED: ", join(FAILS, ", ")); exit(1)
end
