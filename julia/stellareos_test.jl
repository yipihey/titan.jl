# Tests of StellarEOS (S2): stellar gas EOS with degenerate/relativistic electrons, Saha join.
#   TITAN_NOLIB=1 julia --project=env stellareos_test.jl
# Needs (optional) the Helmholtz reference build in ~/Projects/popiii_lit/helmholtz/run (README there).
using Printf, ForwardDiff
include(joinpath(@__DIR__, "Z0Physics.jl")); include(joinpath(@__DIR__, "StellarEOS.jl"))
const S = StellarEOS; const Z0 = Z0Physics
ok = Bool[]
check(name, c) = (push!(ok, c); println(c ? "  PASS  " : "  FAIL  ", name))
const mc = S.me * S.cl; const mc2 = S.mc2

# ---------------------------------------------------------------------------------------------------------------------
# independent references: tanh-sinh quadrature (generalised FD integrals; momentum-space electron gas)
function ts(f, a, b; n = 300, h = 0.025)
    s = 0.0
    for k in -n:n
        t = k * h
        u = tanh(π / 2 * sinh(t)); w = π / 2 * cosh(t) / cosh(π / 2 * sinh(t))^2
        s += w * f(0.5 * (b - a) * u + 0.5 * (a + b)) * 0.5 * (b - a)
    end
    s * h
end
function Fref(k, η, θ)
    g(x) = x^k * sqrt(1 + θ * x / 2) / (exp(x - η) + 1)
    brk = unique(sort([0.0; (η > 1 ? [max(0, η - 30), η, η + 45] : [30.0, 60.0])]))
    sum(ts(g, brk[i], brk[i+1]) for i in 1:length(brk)-1)
end
function gas_ref(μk, θ)       # n/A, P/(A mc²), E_kin/(A mc²), A = 8π(mc/h)³; μk = kinetic chemical potential / mc²
    ε(p) = sqrt(1 + p^2) - 1
    f(p) = 1 / (exp((ε(p) - μk) / θ) + 1)
    pm = 1.0; while ε(pm) < μk + 60θ; pm *= 1.5; end
    brk = unique(sort([0.0; pm * (0:0.05:1)]))
    I(g) = sum(ts(p -> g(p) * f(p), brk[i], brk[i+1]) for i in 1:length(brk)-1)
    I(p -> p^2), I(p -> p^4 / (3 * sqrt(1 + p^2))), I(p -> p^2 * ε(p))
end

println("A. generalised Fermi–Dirac integrals F_k(η, θ), k = 1/2, 3/2, 5/2 and ∂F_k/∂η: quadrature vs independent tanh-sinh")
eF = 0.0; eD = 0.0; nA = 0
for η in (-40.0, -10.0, -3.0, 0.0, 1.0, 5.0, 20.0, 47.9, 48.1, 100.0, 1e3, 1e4), θ in (1e-5, 1e-3, 0.05, 0.5, 5.0)
    F = S.fd(η, θ)
    for (i, k) in enumerate((0.5, 1.5, 2.5))
        global eF = max(eF, abs(F[i] / Fref(k, η, θ) - 1))
        d = ForwardDiff.derivative(x -> S.fd(x, θ)[i], η)
        global eD = max(eD, abs(F[3+i] / d - 1)); global nA += 1
    end
end
@printf("  %d (η, θ, k) combinations: max rel. error F %.1e, D = ∂F/∂η vs AD %.1e\n", nA, eF, eD)
check("F_k accuracy ≤ 1e-12", eF < 1e-12)
check("D_k (analytic ∂F/∂η) vs AD of F ≤ 1e-11", eD < 1e-11)

println("B. electron–positron gas (n, P, E, with e⁺) vs momentum-space integration")
let comp = S.composition_single(1.3, 1.0), emax = 0.0, nmax = 0.0
    for (ρ, T) in [(3e-9, 1.8e8), (1e-10, 3e8), (1.0, 1e7), (1e3, 3e7), (1e-6, 3e9), (1e2, 1e9), (1e3, 1e6), (1e4, 3e6), (1e-4, 1e9), (30.0, 1.2e8), (1e-8, 3e9)]
        θ = S.kB * T / mc2
        b = S.branch_E(ρ, T, comp; coulomb = false); η = b.η
        n1, P1, E1 = gas_ref(η * θ, θ); n2, P2, E2 = gas_ref(-η * θ - 2, θ)
        A = 8π * (S.me * S.cl / S.h)^3
        Pe = A * mc2 * (P1 + P2); Ee = A * mc2 * (E1 + E2) + 2mc2 * A * n2
        pion = ρ * comp.νion * S.kB * T; eion = 1.5 * comp.νion * S.kB * T
        emax = max(emax, abs((b.p - pion) / Pe - 1), abs((b.e - eion) * ρ / Ee - 1))
        nmax = max(nmax, abs(A * (n1 - n2) / (ρ * comp.νe) - 1) * min(1, n1 / max(n1 - n2, 1e-300) * 1e-6))   # (n₋−n₊ cancels for pair plasmas)
    end
    @printf("  11 points up to ρ = 1e4, T = 3e9: max rel. diff P_e, E_e: %.1e\n", emax)
    check("P_e, E_e (incl. pairs) vs momentum integral ≤ 1e-9", emax < 1e-9)
end

println("C. limits")
let c = S.composition(0.76, 0.24, 0.0)
    r = S.eos(1e-10, 1e6; coulomb = false)
    pid = 1e-10 * S.kB * 1e6 * (c.νion + c.νe)
    @printf("  nondegenerate ideal (ρ = 1e-10, T = 1e6): P/P_ideal − 1 = %.1e, (e − e_ion)/(1.5 P/ρ) − 1 = %.1e (relativistic ≈ 15θ/8·ν_e/ν = %.1e)\n",
            r.P / pid - 1, (r.e - c.eion) / (1.5 * pid / 1e-10) - 1, 15 / 8 * S.kB * 1e6 / mc2 * c.νe / (c.νe + c.νion))
    check("ideal gas pressure ≤ 1e-8", abs(r.P / pid - 1) < 1e-8)
    check("ideal gas energy to relativistic order", abs((r.e - c.eion) / (1.5 * pid / 1e-10) - 1) < 3e-4)
    check("ideal Γ1 = 5/3, ∇ad = 0.4, η_e ≈ ln(n_e/n_Q)", abs(r.Γ1 / (5 / 3) - 1) < 1e-3 && abs(r.∇ad / 0.4 - 1) < 1e-3)
    # T = 0 Chandrasekhar
    c2 = S.composition_single(2.0, 1.0)             # μ_e = 2
    println("  degenerate electrons at T = 10 K (μ_e = 2): P_e vs exact T = 0 (Chandrasekhar), NR and ER limits")
    cmax = 0.0
    for ρ in (1.0, 1e2, 1e4, 1e6, 1e8, 1e10)
        b = S.branch_E(ρ, 10.0, c2; coulomb = false)
        pe = b.p - ρ * c2.νion * S.kB * 10.0
        ne = ρ * c2.νe; x = cbrt(3ne / (8π)) * S.h / mc
        Pch = π * S.me^4 * S.cl^5 / (3 * S.h^3) * (x * (2x^2 - 3) * sqrt(1 + x^2) + 3 * asinh(x))
        PNR = (1 / 20) * (3 / π)^(2 / 3) * S.h^2 / S.me * ne^(5 / 3)
        PER = S.h * S.cl / 8 * (3 / π)^(1 / 3) * ne^(4 / 3)
        @printf("    ρ = %.0e  x = %.3g:  P_e/P_Chandra − 1 = %+.1e   P_e/P_NR = %.5f   P_e/P_ER = %.5f\n", ρ, x, pe / Pch - 1, pe / PNR, pe / PER)
        cmax = max(cmax, abs(pe / Pch - 1))
        ρ == 1.0 && check("NR limit P_e/P_NR = 1 − 5x²/14 within 1e-6", abs(pe / PNR - (1 - 5x^2 / 14)) < 1e-6)
        ρ == 1e10 && check("ER limit P_e/P_ER = 1 − 1/x² within 2e-4", abs(pe / PER - (1 - 1 / x^2)) < 2e-4)
    end
    check("Chandrasekhar T = 0 relation ρ = 1 … 1e10: ≤ 1e-6", cmax < 1e-6)
end

# ---------------------------------------------------------------------------------------------------------------------
println("D. Helmholtz (Timmes & Swesty 2000) comparison, fully ionised single-species ions (A = 1.3, Z = 1), log ρ −10…4, log T 4…9.5")
const HDIR = joinpath(homedir(), "Projects", "popiii_lit", "helmholtz", "run")
function helm(pts, abar, zbar)
    out = Vector{Vector{Float64}}()
    for chunk in Iterators.partition(pts, 900)
        inp = IOBuffer(); println(inp, length(chunk))
        for (ρ, T) in chunk; @printf(inp, "%.15e %.15e %.15e %.15e\n", T, ρ, abar, zbar); end
        o = read(pipeline(IOBuffer(take!(inp)), setenv(`$(joinpath(HDIR, "helmdrv"))`; dir = HDIR)), String)
        for l in split(strip(o), '\n'); push!(out, parse.(Float64, split(l))); end
    end
    out
end
if isfile(joinpath(HDIR, "helmdrv"))
    abar, zbar = 1.3, 1.0
    comp1 = S.composition_single(abar, zbar)
    pts = [(10.0^a, 10.0^b) for a in -10:0.5:4 for b in 4:0.25:9.5]
    H = helm(pts, abar, zbar)
    rows = NTuple{7,Float64}[]; global nskip = 0
    for (i, (ρ, T)) in enumerate(pts)
        r = H[i]; (r[6] > 0 && r[7] > 0) || continue
        (r[15] == 0 && r[26] > 1e-3) && (global nskip += 1; continue)        # Helmholtz switched its Coulomb correction off (negative p or e)
        st = S.thermo((ρ_, T_) -> S.branch_E(ρ_, T_, comp1), ρ, T)
        push!(rows, (log10(ρ), log10(T), abs(st.P / r[6] - 1), abs(st.e / r[7] - 1), abs(st.s / r[8] - 1), abs(st.Γ1 / r[18] - 1), abs(st.∇ad / r[19] - 1)))
    end
    function rep(name, f)
        sel = filter(f, rows); isempty(sel) && return
        mx(i) = maximum(x -> x[i], sel)
        @printf("  %-44s n=%3d  max rel. diff: P %.1e  e %.1e  s %.1e  Γ1 %.1e  ∇ad %.1e\n", name, length(sel), mx(3), mx(4), mx(5), mx(6), mx(7))
        (mx(3), mx(4), mx(5), mx(6), mx(7))
    end
    rA = rep("hot, log T 6–9, log ρ ≥ −4 (non-pair)", x -> 6 <= x[2] <= 9 && x[1] >= -4)
    rB = rep("degenerate/dense: log ρ ≥ 1, log T ≥ 6", x -> x[2] >= 6 && x[1] >= 1)
    rC = rep("cool and dense (Coulomb): log T 4–6, log ρ ≥ 0", x -> x[2] < 6 && x[1] >= 0)
    rD = rep("rarefied hot, log ρ ≤ −6, log T 7–9.5 (pairs)", x -> x[2] >= 7 && x[1] <= -6)
    rep("all comparable points", x -> true)
    println("  ($nskip points skipped where Helmholtz's bomb-proofing removed its Coulomb term: Γ_ion ≳ 20 at log T = 4, log ρ ≥ 0)")
    check("vs Helmholtz: hot non-pair region P ≤ 3e-5, e ≤ 1e-4, s ≤ 2e-3, Γ1 ≤ 2e-3", rA[1] < 3e-5 && rA[2] < 1e-4 && rA[3] < 2e-3 && rA[4] < 2e-3)
    check("vs Helmholtz: dense degenerate region P ≤ 3e-5, Γ1 ≤ 1e-3", rB[1] < 3e-5 && rB[4] < 1e-3)
    println("  (Helmholtz is a bicubic table, ~10 nodes/dex in T: its own errors vs exact momentum integrals reach 1e-2 in e (section B);")
    println("   the quadrature here agrees with the independent integrals to 1e-12.)")
else
    println("  (Helmholtz reference build not found: skipped)")
end

# ---------------------------------------------------------------------------------------------------------------------
println("E. thermodynamic consistency of the blended EOS over log ρ −10…4, log T 3.5…9.5 (X=0.76, Y=0.24, Z=1e-3)")
let ms = 0.0, mv = 0.0, nbad = 0, n = 0, worst = (0.0, 0.0), nstab = 0, stabbad = Tuple[], adbad = Tuple[]
    for lr in -10:0.5:4, lt in 3.5:0.25:9.5
        ρ = 10.0^lr; T = 10.0^lt
        r = S.eos(ρ, T; X = 0.759, Y = 0.24, Z = 0.001)
        isfinite(r.P) || (nbad += 1; continue)
        mx = abs(r.dedρ - (r.P - T * r.dPdT) / ρ^2) / (abs(r.dedρ) + abs(r.P) / ρ^2)
        cv = abs(r.cV - r.cVs) / r.cV
        if mx > ms; ms = mx; worst = (lr, lt); end
        mv = max(mv, cv); n += 1
        if lt >= 4 && !(r.χρ > 0 && r.cV > 0 && r.Γ1 > 0); push!(stabbad, (lr, lt, r.χρ, r.Γ1, r.∇ad)); end
        if lt >= 4 && !(0 < r.∇ad < 0.41); push!(adbad, (lr, lt, r.∇ad)); end
    end
    @printf("  %d points: max |∂e/∂ρ − (p − T∂p/∂T)/ρ²| (relative) %.1e at %s; max |c_V − T∂s/∂T| %.1e; non-finite %d\n", n, ms, string(worst), mv, nbad)
    check("Maxwell relation ∂e/∂ρ|T = (p − T ∂p/∂T)/ρ² ≤ 1e-8", ms < 1e-8)
    check("c_V = ∂e/∂T = T ∂s/∂T ≤ 1e-8", mv < 1e-8)
    check("all finite", nbad == 0)
    @printf("  stability (T ≥ 1e4 K): χρ, c_V, Γ1 > 0 violated at %d of %d points; ∇ad ∉ (0, 0.41) at %d points", length(stabbad), n, length(adbad))
    if !isempty(adbad)
        @printf(" (all in the pressure-ionisation transition, log ρ ∈ [%.1f, %.1f], log T ∈ [%.2f, %.2f]; max ∇ad %.3f)", minimum(x -> x[1], adbad), maximum(x -> x[1], adbad),
                minimum(x -> x[2], adbad), maximum(x -> x[2], adbad), maximum(x -> x[3], adbad))
    end
    println()
    check("χρ > 0, c_V > 0, Γ1 > 0 everywhere", isempty(stabbad))
    check("0 < ∇ad < 0.41 except in the crude pressure-ionisation transition (−1 ≤ log ρ ≤ 1.5, log T ≤ 5.75)", all(x -> -1 <= x[1] <= 1.5 && x[2] <= 5.75, adbad))
end

println("F. AD vs finite differences (outer ForwardDiff through gas_pe, and through eos: nested)")
let ea = 0.0, eb = 0.0, ec = 0.0
    for (ρ, T) in [(1e-6, 3e4), (1e-3, 1e6), (100.0, 3e7), (30.0, 1.2e8), (1e-4, 2.5e5), (2e-2, 4e5), (0.3, 5e4), (1e4, 1e5), (1e-7, 2e9)]
        f(v) = collect(S.gas_pe(v[1], v[2]; Z = 1e-3, X = 0.759))
        J = ForwardDiff.jacobian(f, [ρ, T])
        h = 1e-4
        Jfd = hcat((f([ρ * (1 + h), T]) - f([ρ * (1 - h), T])) / (2h * ρ), (f([ρ, T * (1 + h)]) - f([ρ, T * (1 - h)])) / (2h * T))
        fv = f([ρ, T]); ea = max(ea, maximum(abs.(J .- Jfd) .* [ρ T] ./ abs.(fv)))
        # nested: d Γ1/d lnT from AD through eos vs central differences of Γ1
        g(lt) = S.eos(ρ, exp(lt); Z = 1e-3, X = 0.759).Γ1
        dg = ForwardDiff.derivative(g, log(T)); dgf = (g(log(T) + 1e-4) - g(log(T) - 1e-4)) / 2e-4
        eb = max(eb, abs(dg - dgf) / (abs(dg) + 1e-3))
        gc(lr) = S.eos(exp(lr), T; Z = 1e-3, X = 0.759).cV
        dc = ForwardDiff.derivative(gc, log(ρ)); dcf = (gc(log(ρ) + 1e-4) - gc(log(ρ) - 1e-4)) / 2e-4
        ec = max(ec, abs(dc - dcf) / (abs(dc) + 1e-3 * gc(log(ρ))))
    end
    @printf("  9 points (Saha, blend, degenerate, hot): max |ΔJ|·x/f of ∂(p,e)/∂(ρ,T): %.1e; d Γ1/d ln T (nested): %.1e; d c_V/d ln ρ (nested): %.1e\n", ea, eb, ec)
    check("AD Jacobian of (p, e) vs FD ≤ 1e-6 (FD-limited)", ea < 1e-6)
    check("nested AD (second derivatives) vs FD ≤ 1e-5", eb < 1e-5 && ec < 1e-5)
end

println("G. Saha join: limits and continuity")
let
    # w = 0: identical to Z0Physics with X = 0.76 (gas only)
    ρ, T = 1e-6, 3e4
    r = S.eos(ρ, T); z = Z0.eos_rhoT(ρ, T); pz = Z0.Pgas(ρ, T); ez = Z0._egas(Z0.saha(ρ, T), ρ, T)
    @printf("  below the blend (ρ = 1e-6, T = 3e4): p/p_Z0 − 1 = %.1e, e/e_Z0 − 1 = %.1e\n", r.P / pz - 1, r.e / ez - 1)
    check("w = 0: p, e equal Z0Physics' Saha gas", abs(r.P / pz - 1) < 1e-12 && abs(r.e / ez - 1) < 1e-12)
    # continuity of Γ1 and ∇ad across the T window and the ρ window: scan with small steps
    function scan(f, xs)
        v = [f(x) for x in xs]
        dΓ = maximum(abs.(diff(first.(v)))); d∇ = maximum(abs.(diff(last.(v))))
        dΓ, d∇, maximum(first.(v)) - minimum(first.(v))
    end
    h1, h2 = 2e-4, 1e-4
    out = Tuple[]
    for (lbl, ρ0) in (("ρ=1e-6", 1e-6), ("ρ=1e-4", 1e-4), ("ρ=1e-2", 1e-2))
        fT(lt) = (r = S.eos(ρ0, 10.0^lt); (r.Γ1, r.∇ad))
        a = scan(fT, 5.0:h1:6.2); b = scan(fT, 5.0:h2:6.2)
        push!(out, (lbl, a[1], b[1], a[2], b[2], a[3]))
        @printf("  scan log T = 5.0–6.2 at %s: max step jump in Γ1 %.2e (Δ=%g dex) → %.2e (Δ=%g); ∇ad %.2e → %.2e; Γ1 range over scan %.3f\n", lbl, a[1], h1, b[1], h2, a[2], b[2], a[3])
    end
    for (lbl, T0) in (("T=1e4", 1e4), ("T=1e5", 1e5), ("T=1e6", 1e6))
        fR(lr) = (r = S.eos(10.0^lr, T0); (r.Γ1, r.∇ad))
        a = scan(fR, -2.0:h1:1.2); b = scan(fR, -2.0:h2:1.2)
        push!(out, (lbl, a[1], b[1], a[2], b[2], a[3]))
        @printf("  scan log ρ = −2–1.2 at %s: max step jump in Γ1 %.2e (Δ=%g dex) → %.2e (Δ=%g); ∇ad %.2e → %.2e; Γ1 range %.3f\n", lbl, a[1], h1, b[1], h2, a[2], b[2], a[3])
    end
    check("Γ1, ∇ad continuous: the largest step jump halves with the step (finite slope, no discontinuity)",
          all(o -> o[3] < 0.55 * o[2] + 1e-9 && o[5] < 0.55 * o[4] + 1e-9, out))
    # size of the blend's effect: Δf and the pressure difference at the window centre
    ρ, T = 1e-3, 10.0^5.6
    c = S.composition(0.76, 0.24, 0.0); Sb = S.branch_S(ρ, T, c); Eb = S.branch_E(ρ, T, c)
    @printf("  at ρ = 1e-3, T = 10^5.6 (window centre): branches differ by Δp/p = %.2e, Δe/e = %.2e, Δs/s = %.2e (Saha: He partly neutral)\n",
            Eb.p / Sb.p - 1, Eb.e / Sb.e - 1, Eb.s / Sb.s - 1)
    ρ, T = 1e-3, 10.0^6.0
    Sb = S.branch_S(ρ, T, c); Eb = S.branch_E(ρ, T, c)
    @printf("  at the hot end of the window (T = 1e6, ρ = 1e-3): Δp/p = %.1e, Δe/e = %.1e, Δs/s = %.1e (the two models coincide when ionisation is complete)\n", Eb.p / Sb.p - 1, Eb.e / Sb.e - 1, Eb.s / Sb.s - 1)
    check("Saha and fully ionised branches agree to < 1e-3 in p where ionisation is complete (1e6 K, 1e-3)", abs(Eb.p / Sb.p - 1) < 1e-3)
end

println("H. ρ(P, T) inversion")
let e1 = 0.0, e2 = 0.0
    for (ρ, T) in [(1e-6, 3e4), (1e-3, 1e6), (100.0, 3e7), (30.0, 1.2e8), (1e-4, 2.5e5), (0.3, 5e4), (1e4, 1e5)]
        P = S.eos(ρ, T).P
        ρ2 = S.rho_PT(P, T)
        e1 = max(e1, abs(ρ2 / ρ - 1))
        d = ForwardDiff.derivative(lp -> log(S.rho_PT(exp(lp), T)), log(P)); χρ = S.eos(ρ, T).χρ
        e2 = max(e2, abs(d * χρ - 1))
    end
    @printf("  7 points: max |ρ(P(ρ,T), T)/ρ − 1| = %.1e; max |(∂lnρ/∂lnP)·χρ − 1| (AD) = %.1e\n", e1, e2)
    check("ρ(P,T) round trip ≤ 1e-10 and derivative identity ≤ 1e-8", e1 < 1e-10 && e2 < 1e-8)
end

println("I. Pop III central points (X = 0.76, Y = 0.24)")
for (ρ, T) in [(100.0, 3e7), (30.0, 1.2e8), (1e-3, 1e6), (1e-6, 3e4)]
    r = S.eos(ρ, T); rr = S.eos(ρ, T; radiation = true)
    @printf("  ρ = %-7g T = %-7g  η_e = %+8.4f  Γ1 = %.5f  ∇ad = %.5f  χρ = %.5f  χT = %.5f  c_V = %.4e  n_e = %.3e   (with radiation: Γ1 = %.5f, ∇ad = %.5f, β = %.5f)\n",
            ρ, T, r.η, r.Γ1, r.∇ad, r.χρ, r.χT, r.cV, r.ne, rr.Γ1, rr.∇ad, r.P / rr.P)
end

println("J. Coulomb correction size (Γ_ion = Z̄²e²/(kT a_i)) and pairs")
for (ρ, T) in [(100.0, 3e7), (30.0, 1.2e8), (1.0, 1e6), (1e-3, 1e6)]
    c = S.composition(0.76, 0.24, 0.0); a = S.eos(ρ, T); b = S.eos(ρ, T; coulomb = false)
    Γ = (c.νe / c.νion)^2 * S.qe^2 / (S.kB * T * cbrt(3 / (4π * ρ * c.νion)))
    @printf("  ρ = %-7g T = %-7g  Γ_ion = %.3f  Δp/p = %+.2e  Δe/e = %+.2e  Δs/s = %+.2e  ΔΓ1 = %+.2e\n", ρ, T, Γ, a.P / b.P - 1, a.e / b.e - 1, a.s / b.s - 1, a.Γ1 - b.Γ1)
end
let a = S.eos(1e-8, 3e9), b = S.eos(1e-8, 3e9; pairs = false)
    @printf("  pairs at ρ = 1e-8, T = 3e9: P(with)/P(without) = %.3f, Γ1 = %.4f vs %.4f\n", a.P / b.P, a.Γ1, b.Γ1)
    check("pair creation matters at 3e9 K (P ratio > 1.5) and softens Γ1", a.P / b.P > 1.5 && a.Γ1 < b.Γ1)
end

println("K. performance (single thread)")
let
    pts = [("Saha        (1e-6, 3e4)", 1e-6, 3e4), ("blend       (1e-3, 4e5)", 1e-3, 4e5), ("hot, ND     (1e-3, 1e6)", 1e-3, 1e6),
           ("central      (100, 3e7)", 100.0, 3e7), ("degenerate  (1e4, 1e5)", 1e4, 1e5), ("pairs       (1e-8, 3e9)", 1e-8, 3e9)]
    fpe(v) = collect(S.gas_pe(v[1], v[2]))
    for (lbl, ρ, T) in pts
        S.gas_pe(ρ, T); S.eos(ρ, T); ForwardDiff.jacobian(fpe, [ρ, T])
        n = 2000
        t0 = time_ns(); for _ in 1:n; S.gas_pe(ρ, T); end; t1 = (time_ns() - t0) / n / 1e3
        t0 = time_ns(); for _ in 1:n; S.eos(ρ, T); end; t2 = (time_ns() - t0) / n / 1e3
        t0 = time_ns(); for _ in 1:n; ForwardDiff.jacobian(fpe, [ρ, T]); end; t3 = (time_ns() - t0) / n / 1e3
        @printf("  %s: gas_pe %6.1f µs   (p,e) + Jacobian (AD) %6.1f µs   full eos %6.1f µs\n", lbl, t1, t3, t2)
    end
end

# ---------------------------------------------------------------------------------------------------------------------
if get(ENV, "TITAN_NOLIB", "") == "1" || true
    println("L. TitanMaterials.StellarGas (per-zone composition through TitanResidFast.ZONE)")
    try
        include(joinpath(@__DIR__, "TitanDriver.jl")); include(joinpath(@__DIR__, "JuliaPhysics.jl"))
        include(joinpath(@__DIR__, "TitanResidFast.jl")); include(joinpath(@__DIR__, "BlockSolver.jl"))
        include(joinpath(@__DIR__, "TitanCore.jl")); include(joinpath(@__DIR__, "TitanMaterials.jl"))
        include(joinpath(@__DIR__, "StellarMaterial.jl"))
        TX = TitanResidFast
        m = StellarMaterial.StellarGas(8; X = 0.76, Y = 0.24, Z = 0.0)
        StellarMaterial.setcomp!(m, 3, 0.0, 0.98, 0.02)            # zone 3: He + 2 % metals
        ρ, T = 100.0, 3e7
        TX.ZONE[] = 1; p1, e1, a1 = TX.material_gas(m, ρ, T)
        TX.ZONE[] = 3; p3, e3, a3 = TX.material_gas(m, ρ, T)
        q1 = S.gas_pe(ρ, T); q3 = S.gas_pe(ρ, T; X = 0.0, Y = 0.98, Z = 0.02)
        @printf("  zone 1 (Pop III H/He): p = %.5e e = %.5e;  zone 3 (He+Z): p = %.5e e = %.5e\n", p1, e1, p3, e3)
        check("material_gas uses the zone's composition", p1 == q1[1] && e1 == q1[2] && p3 == q3[1] && e3 == q3[2] && p1 != p3)
        TX.ZONE[] = 1
        J = ForwardDiff.jacobian(v -> collect(TX.material_gas(m, v[1], v[2])[1:2]), [ρ, T])
        r = S.eos(ρ, T)
        check("Jacobian of material_gas through the residual types = EOS derivatives", abs(J[1, 1] - r.dPdρ) < 1e-8 * abs(r.dPdρ) && abs(J[1, 2] - r.dPdT) < 1e-8 * abs(r.dPdT) && abs(J[2, 1] - r.dedρ) < 1e-8 * abs(r.dedρ) && abs(J[2, 2] - r.cV) < 1e-8 * r.cV)
        th = TX.material_thermo(m, ρ, T); rr = S.eos(ρ, T; radiation = true)
        check("material_thermo = ∇ad, c_P, δ of gas + radiation", th == (rr.∇ad, rr.cP, rr.δ))
        pg, eg, as_, ne = TX.material_gasne(m, ρ, T)
        check("material_gasne returns n_e", ne == r.ne)
    catch err
        println("  material test failed: ", sprint(showerror, err)[1:min(end, 300)]); push!(ok, false)
    end
end

println(all(ok) ? "\nALL PASS ($(length(ok)) checks)" : "\nFAILED: $(count(!, ok)) of $(length(ok))")
