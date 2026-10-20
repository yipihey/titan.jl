# Population III project, S3: acceptance tests of StellarOpacity (Rosseland/Planck means, Z = 0, OPAL + own low-T means).
#   julia -t 4 --project=env stellaropacity_test.jl
# Needs the run-time table julia/stellar/kaptable_stellar.bin.gz (stellaropacity_build.jl) and, for the comparison with the
# raw data, ~/Projects/popiii_lit/opal/ (OPAL GN93hz, Cassisi conduction tables).
using Printf, ForwardDiff, Random, Statistics
const J = @__DIR__
include(joinpath(J, "PrimordialChem.jl")); include(joinpath(J, "ContinuumOpacity.jl")); include(joinpath(J, "Z0Physics.jl"))
include(joinpath(J, "StellarOpacity.jl")); include(joinpath(J, "StellarOpacityLow.jl"))
include(joinpath(J, "stellaropacity_data.jl"))
const SO = StellarOpacity; const SL = StellarOpacityLow; const PC = PrimordialChem; const Z0 = Z0Physics
results = Tuple{String,Bool,String}[]
check(name, ok, detail) = (push!(results, (name, ok, detail)); @printf("%-34s %s  %s\n", name, ok ? "PASS" : "FAIL", detail))
Random.seed!(20261020)

tabs, lTo, lRo = read_gn93(OPALFILE)
t = SO.table()
lTf(i) = SO.LT0 + SO.DLT * (i - 1); lRf(j) = SO.LR0 + SO.DLR * (j - 1)
iT(lT) = round(Int, (lT - SO.LT0) / SO.DLT) + 1; iR(lR) = round(Int, (lR - SO.LR0) / SO.DLR) + 1

# 1. reproduction of the OPAL grid points -------------------------------------------------------------------------------------
let
    worst = 0.0; n = 0; where_ = ""
    for (ix, X) in enumerate(SO.XNODES), (i, lT) in enumerate(lTo), (j, lR) in enumerate(lRo)
        a = tabs[X][i, j]; isnan(a) && continue
        lT < 4.1 - 1e-9 && continue
        (abs(lT - 8.7) < 1e-9 && lR <= -7) && continue                  # placeholders in OPAL, masked
        v = SO.log_kappa_rad(lT, lR, X); e = abs(v - a); n += 1
        e > worst && (worst = e; where_ = @sprintf("X=%.2f lT=%.2f lR=%.1f", X, lT, lR))
    end
    check("opal_grid_points", worst <= 1e-6, @sprintf("%d nodes (log T 4.1–8.7), max |Δ log κ| = %.2e at %s", n, worst, where_))
    # reference points, including the one named in the work package
    s = ""; ok = true
    for (lT, lR, X) in ((7.0, -3.0, 0.7), (6.0, -2.0, 0.5), (5.0, 0.0, 0.35), (4.5, -4.0, 0.0), (8.3, -4.0, 1.0), (4.2, -6.0, 0.9))
        a = tabs[X][findfirst(x -> abs(x - lT) < 1e-9, lTo), findfirst(x -> abs(x - lR) < 1e-9, lRo)]
        v = SO.log_kappa_rad(lT, lR, X); ok &= abs(v - a) <= 1e-6
        s *= @sprintf(" (%.1f,%.0f,%.2f): %.3f=%.3f", lT, lR, X, v, a)
    end
    check("opal_reference_points", ok, "table vs raw:" * s)
end

# 2. off-node accuracy of the fine-grid bicubic against OPAL's own Hermite interpolant ---------------------------------------
let
    X = 0.7; A = tabs[X]; e = Float64[]
    for _ in 1:2000
        lT = 4.2 + 2.7 * rand(); lR = -7.5 + 8.0 * rand()               # rows 4.0–7.0 are complete up to log R = 1
        rows = [herm(lRo, A[i, :], lR) for i in eachindex(lTo)]
        push!(e, SO.log_kappa_rad(lT, lR, X) - herm(lTo, rows, lT))
    end
    check("opal_offnode_interp", maximum(abs, e) < 5e-3, @sprintf("%d random points, X=0.7: rms %.1e, max %.1e dex vs OPAL's cubic Hermite", length(e), sqrt(mean(abs2, e)), maximum(abs, e)))
end

# 3. C¹ continuity: jumps of the gradient across cell edges, blend centres and X nodes -------------------------------------------
let
    f(lT, lR, X) = SO.log_kappa_R(lT, lR, X)
    grad(lT, lR, X) = ForwardDiff.gradient(v -> f(v[1], v[2], v[3]), [lT, lR, X])
    jumpT = 0.0; jumpR = 0.0; jumpX = 0.0; wT = ""
    δ = 1e-9
    for X in (0.05, 0.3, 0.62, 0.75, 0.97), lR in (-12.3, -6.1, -3.3, -1.1, 0.4, 2.2, 3.7), lT in (2.1, 3.0, 3.8, 3.95, 4.1, 4.5, 5.5, 6.0, 6.5, 7.0, 8.0, 8.4, 8.7, 9.2)
        # edges in log T at grid nodes, at the blend knots and in log R
        g1 = grad(lT - δ, lR, X); g2 = grad(lT + δ, lR, X); d = abs(g1[1] - g2[1]); d > jumpT && (jumpT = d; wT = @sprintf("lT=%.2f lR=%.1f X=%.2f", lT, lR, X))
        jumpR = max(jumpR, abs(grad(lT, lR - δ, X)[2] - grad(lT, lR + δ, X)[2]))
    end
    for lT in (3.9, 5.1, 6.3, 7.7), lR in (-9.4, -4.4, -0.7, 0.6), X in SO.XNODES[2:end-1]
        jumpX = max(jumpX, abs(grad(lT, lR, X - δ)[3] - grad(lT, lR, X + δ)[3]))
    end
    for lT in (3.9, 5.1, 6.3, 7.7), lR in (-8.0, -7.75, 0.0, 1.0, 1.25, -15.0 + 0.25)         # across log R nodes, OPAL's edges (−8, +1)
        jumpR = max(jumpR, abs(grad(lT, lR - δ, 0.6)[2] - grad(lT, lR + δ, 0.6)[2]))
    end
    check("c1_continuity", max(jumpT, jumpR, jumpX) < 1e-5,
          @sprintf("max jump d lnκ: in log T %.1e (%s), in log R %.1e, in X %.1e (at the X nodes)", jumpT, wT, jumpR, jumpX))
    # no kinks in the blend zones: second differences of the fine grid relative to the pure-OPAL region
    A = t.lkR[:, :, 6]; ratio = Dict{String,Float64}()
    for (nm, a, b) in (("own↔OPAL blend 3.8–4.1", 3.8, 4.1), ("OPAL 4.1–5.5", 4.1, 5.5), ("anchor ramp 5.5–6.5", 5.5, 6.5), ("OPAL 6.5–8.4", 6.5, 8.4))
        m = 0.0
        for i in 2:SO.NT-1, j in iR(-8):iR(1)
            a <= lTf(i) < b && (m = max(m, abs(A[i+1, j] - 2A[i, j] + A[i-1, j]) / SO.DLT^2))
        end
        ratio[nm] = m
    end
    @printf("   max |d² log κ/d log T²| at X = 0.7 (log R −8…1): %s\n", join([@sprintf("%s: %.0f", k, v) for (k, v) in sort(collect(ratio))], "; "))
end

# 4. AD vs finite differences -----------------------------------------------------------------------------------------------------
let
    worst = 0.0; w = ""
    for _ in 1:60
        lT = 2.0 + 7.0 * rand(); lR = -14 + 18 * rand(); X = 0.75 * rand()
        T = 10^lT; ρ = 10^(lR + 3 * (lT - 6))
        for fn in (SO.kappa_R, SO.kappa_P, SO.kappa_R_rad)
            F(v) = log(fn(exp(v[1]), exp(v[2]), v[3]))
            v0 = [log(ρ), log(T), X]
            g = ForwardDiff.gradient(F, v0)
            ek(k) = [m == k ? 1e-6 : 0.0 for m in 1:3]
            fd = [(F(v0 .+ ek(k)) - F(v0 .- ek(k))) / 2e-6 for k in 1:3]
            e = maximum(abs.(g .- fd) ./ (1 .+ abs.(fd)))
            e > worst && (worst = e; w = @sprintf("%s lT=%.2f lR=%.1f X=%.2f", nameof(fn), lT, lR, X))
        end
    end
    check("ad_vs_fd", worst < 1e-4, @sprintf("60 random points × (κ_R, κ_P, κ_R,rad): max |AD − FD|/(1+|FD|) of dlnκ/d(lnρ, lnT, X) = %.1e (%s)", worst, w))
end

# 5. electron-scattering limit ------------------------------------------------------------------------------------------------------
let
    ok = true; s = ""
    for X in (0.0, 0.35, 0.7, 0.77), lT in (6.0, 7.0, 8.0, 9.0)
        T = 10^lT; ρ = 1e-12 * (T / 1e6)^3 * 1e3            # log R = −9
        κ = SO.kappa_R(ρ, T, X); ref = 0.2 * (1 + X) * SO.thermal_factor(T)
        ok &= abs(κ / ref - 1) < 0.02
        X == 0.7 && (s *= @sprintf(" T=1e%.0f: %.4f", lT, κ / (0.2 * (1 + X))))
    end
    κ0 = SO.kappa_R(1e-20, 1e6, 0.7)
    check("electron_scattering_limit", ok && abs(κ0 / 0.34 - 1) < 0.01, @sprintf("κ/[0.2(1+X)] at X=0.7, log R = −9:%s (thermal factor 1/(1+(T/4.5e8)^0.95)); T=1e6, ρ=1e-20: %.4f", s, κ0 / 0.34))
end

# 6. low-T means vs M3c: Rosseland/Planck at T = 1500–6000 K, n_H = 1e14–1e21 --------------------------------------------------------------
let
    # (n_H, T, κ_P, κ_R,ext without Li) of the M3c table (NOTES M3c "First means")
    m3c = ((1e14, 1500.0, 4.5e-8, 1.2e-9), (1e16, 2000.0, 5.0e-6, 7.4e-6), (1e18, 3000.0, 3.4e-4, 5.6e-5), (1e20, 4000.0, 2.3e-2, 3.1e-3), (1e21, 6000.0, 6.4, 4.9))
    okm = true; okt = true; line1 = ""; line2 = ""
    for (nH, T, kP, kR) in m3c
        ρ = nH * PC.μH
        # (a) the M3c physics without Li electrons, our solver
        m = SL.means_at(ρ, T, 0.75; li = false)
        okm &= abs(m.κP / kP - 1) < 0.15 && abs(m.κR / kR - 1) < 0.1
        line1 *= @sprintf("\n     n_H=%.0e T=%4.0f: κ_P %.2e (M3c %.1e), κ_R %.2e (M3c %.1e)  [own solver, no Li]", nH, T, m.κP, kP, m.κR, kR)
        # (b) the table (X = 0.75, with Li) against the direct calculation
        d = SL.means_at(ρ, T, 0.75); tab = SO.kappas(ρ, T, 0.75)
        eR = tab[1] / d.κR - 1; eP = tab[2] / d.κP - 1
        okt &= abs(eR) < 0.06 && abs(eP) < 0.06
        line2 *= @sprintf("\n     n_H=%.0e T=%4.0f: κ_R ρ table %.3e direct %.3e (%+.1f%%), κ_P %+.1f%%", nH, T, tab[1] * ρ, d.κR * ρ, 100eR, 100eP)
    end
    check("m3c_lowT_means_nolithium", okm, "recomputed M3c table (no Li) with the new solver within 15 %/10 %:" * line1)
    check("m3c_lowT_table_vs_direct", okt, "table at X = 0.75 (Li on) vs direct means, ≤ 6 %:" * line2)
    # composition agrees with PrimordialChem.y_saha at X = 0.75 (same Saha chemistry)
    ρ = 1e18 * PC.μH; T = 3000.0; y = PC.y_saha(ρ, T); s = SL.equilibrium(ρ, T, 0.75; li = false)
    check("equilibrium_vs_ysaha", abs(2 * s.nH2 / s.nHn / (2 * y[3]) - 1) < 1e-6, @sprintf("f_H2 %.6f vs PrimordialChem.y_saha %.6f (n_H = 1e18, 3000 K, X = 0.75)", 2 * s.nH2 / s.nHn, 2 * y[3]))
    # fast implementation equals ContinuumOpacity.means on the same composition (He Rayleigh aside)
    ρ = 1e18 * PC.μH; T = 2000.0; s = SL.equilibrium(ρ, T, 0.75); mc = ContinuumOpacity.means(ρ, T; nH = s.nH, nH2 = s.nH2, nHe = s.nHe, nHp = s.nHp, ne = s.ne, nHm = s.nHm)
    m = SL.means_at(ρ, T, 0.75)
    check("lowT_vs_ContinuumOpacity", abs(m.κP / mc.κP - 1) < 1e-3 && abs(m.κR / mc.κR - 1) < 0.02, @sprintf("κ_P ratio %.5f, κ_R ratio %.4f (He Rayleigh)", m.κP / mc.κP, m.κR / mc.κR))
end

# 7. comparison with the Z0Physics table (X = 0.76, continuum only) -------------------------------------------------------------------
let
    z0 = Z0.table(); X = 0.76; rows = String[]; okP = true; okR = true
    ratios = Dict{Float64,Vector{Float64}}()
    for lT in 4.0:0.5:8.0
        rR = Float64[]; rP = Float64[]
        for lR in -8:1.0:1
            T = 10^lT; ρ = 10^(lR + 3 * (lT - 6))
            push!(rR, SO.kappa_R_rad(ρ, T, X) / Z0.kap(z0, ρ, T)); push!(rP, SO.kappa_P(ρ, T, X) / Z0.kapP(z0, ρ, T))
        end
        push!(rows, @sprintf("log T=%.1f: κ_R(OPAL)/κ_R(Z0) median %.2f (min %.2f, max %.2f); κ_P ours/Z0 median %.3f (min %.3f, max %.3f)", lT, median(rR), minimum(rR), maximum(rR), median(rP), minimum(rP), maximum(rP)))
        lT >= 4.5 && lT <= 7.0 && (okP &= all(abs.(rP[1:5] .- 1) .< 0.01))          # log R ≤ −4: non-degenerate, same continuum physics
        lT >= 5.0 && lT <= 5.0 && (okR &= maximum(rR) > 1.5)
    end
    check("vs_Z0Physics_table", okP && okR, "X = 0.76, log R −8…1:\n     " * join(rows, "\n     "))
end

# 8. conduction -----------------------------------------------------------------------------------------------------------------------
let
    cond = read_cond(CONDFILE)
    # independent evaluation at table nodes: κ = 16σT³/(3ρK)
    lρ = 5.0; lT = 7.0; K = 10^cond[1.0][round(Int, (lρ + 6) / 0.25) + 1, round(Int, (lT - 3) / 0.25) + 1]
    ref = 16 * 5.670374e-5 * 1e7^3 / (3 * 1e5 * K); v = SO.kappa_cond(1e5, 1e7, 1.0)
    # solar-centre-like and degenerate He-core-like points
    sun = (SO.kappa_cond(150.0, 1.5e7, 0.7) / SO.kappa_R_rad(150.0, 1.5e7, 0.7), SO.kappa_R(150.0, 1.5e7, 0.7) / SO.kappa_R_rad(150.0, 1.5e7, 0.7))
    core = (SO.kappa_R_rad(1e6, 1e7, 0.0), SO.kappa_cond(1e6, 1e7, 0.0), SO.kappa_R(1e6, 1e7, 0.0))
    check("conduction", abs(v / ref - 1) < 1e-8 && sun[1] > 100 && sun[2] > 0.99 && core[3] < core[1] && core[2] < core[1],
          @sprintf("H node ρ=1e5 T=1e7: κ_cond %.4e (raw %.4e); solar-centre-like κ_cond/κ_rad = %.0f, κ_R/κ_rad = %.4f; He core ρ=1e6, 1e7 K: κ_rad %.2e, κ_cond %.2e, total %.2e", v, ref, sun[1], sun[2], core...))
end

# 9. coverage, clamping, Planck ≤ ... , finite everywhere ----------------------------------------------------------------------------------
let
    ok = all(isfinite, t.lkR) && all(isfinite, t.lkP) && minimum(t.lkR) > -45 && maximum(t.lkR) < 15
    counts = [count(==(UInt8(k)), t.flag[:, :, 6]) for k in 1:4]
    check("table_finite_coverage", ok, @sprintf("log κ_R ∈ [%.1f, %.1f], log κ_P ∈ [%.1f, %.1f]; X=0.7 node counts: OPAL data %d, extended %d, blend %d, own means %d", minimum(t.lkR), maximum(t.lkR), minimum(t.lkP), maximum(t.lkP), counts...))
end

# 10. performance ---------------------------------------------------------------------------------------------------------------------
let
    N = 200_000; ρs = 10 .^ (-12 .+ 8 * rand(N)); Ts = 10 .^ (3.5 .+ 4 * rand(N)); Xs = 0.7 * rand(N)
    s = 0.0; SO.kappa_R(1e-5, 1e5, 0.7)
    t1 = @elapsed for i in 1:N; s += SO.kappa_R(ρs[i], Ts[i], Xs[i]); end
    t2 = @elapsed for i in 1:N; s += SO.kappa_P(ρs[i], Ts[i], Xs[i]); end
    g(v) = SO.kappa_R(v[1], v[2], v[3]); v0 = [1e-5, 1e5, 0.7]; ForwardDiff.gradient(g, v0)
    t3 = @elapsed for i in 1:20_000; ForwardDiff.gradient(g, [ρs[i], Ts[i], Xs[i]]); end
    check("performance", t1 / N < 5e-6, @sprintf("κ_R %.2f µs, κ_P %.2f µs, κ_R with gradient (ρ, T, X) %.2f µs per call", 1e6t1 / N, 1e6t2 / N, 1e6t3 / 20_000))
end

nf = count(r -> !r[2], results)
println(nf == 0 ? "ALL PASS ($(length(results)) checks)" : "FAILED: $nf of $(length(results))")
exit(nf == 0 ? 0 : 1)
