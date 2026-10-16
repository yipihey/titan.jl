push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; const Q = QuasiStarInfall
using ForwardDiff, Printf, Random
tab = load_table()
@assert !any(isnan, tab.lkP) "v2 table with Planck mean not loaded"
dlnk(f, ρ, T) = ForwardDiff.derivative(t -> log(f(tab, ρ, exp(t))), log(T))
println("== 1. Derivative continuity across table cell edges (d ln κ_R / d ln T at T-edge ± 1e-9 dex) ==")
jb = Float64[]; jc = Float64[]
rng = MersenneTwister(1)
for _ in 1:400
    i = rand(rng, 10:length(tab.lT)-10); j = rand(rng, 10:length(tab.lR)-10)
    lT = tab.lT[i]; lR = tab.lR[j] + 0.37 * (tab.lR[2] - tab.lR[1])
    f(lt, which) = (T = 10^lt; ρ = 10^(lR + 3(lt - 6)) ; which(tab, ρ, T))
    for (arr, which) in ((jb, Q.kap_bilinear), (jc, Q.kap))
        d(lt) = ForwardDiff.derivative(x -> log10(which(tab, 10^(lR + 3(lT - 6)), 10^x)), lt)
        a = d(lT - 1e-9); b = d(lT + 1e-9)
        push!(arr, abs(a - b) / max(abs(a), abs(b), 1e-3))
    end
end
@printf("  bilinear: median jump %.2e, 99th pct %.2e, max %.2e (relative)\n", sort(jb)[200], sort(jb)[396], maximum(jb))
@printf("  bicubic : median jump %.2e, 99th pct %.2e, max %.2e (relative)\n", sort(jc)[200], sort(jc)[396], maximum(jc))

println("== 2. Interpolation accuracy vs direct calculation at random off-grid points ==")
for (lTlo, lThi, lab) in ((3.5, 3.9, "3.5<logT<3.9 (H− regime)"), (3.9, 4.6, "3.9<logT<4.6 (H/He ionization)"), (4.6, 7.5, "4.6<logT<7.5"))
    eR = Float64[]; eRb = Float64[]; eP = Float64[]
    for _ in 1:60
        lT = lTlo + (lThi - lTlo) * rand(rng); lR = -11.0 + 9.5 * rand(rng)
        T = 10^lT; ρ = 10^(lR + 3(lT - 6))
        kR, kP = Q.kappa_means(ρ, T)
        push!(eR, abs(Q.kap(tab, ρ, T) / kR - 1)); push!(eRb, abs(Q.kap_bilinear(tab, ρ, T) / kR - 1))
        push!(eP, abs(Q.kapP(tab, ρ, T) / kP - 1))
    end
    @printf("  %-32s κ_R bicubic: med %.1e max %.1e | bilinear: med %.1e max %.1e | κ_P bicubic: med %.1e max %.1e\n",
            lab, sort(eR)[30], maximum(eR), sort(eRb)[30], maximum(eRb), sort(eP)[30], maximum(eP))
end

println("== 3. Planck vs Rosseland at representative points ==")
for (ρ, T) in ((4e-11, 6241.0), (1e-10, 1e4), (1e-9, 3e4), (1e-8, 3e5), (1e-16, 6000.0))
    @printf("  ρ=%.0e T=%.0e: κ_R=%.3e κ_P=%.3e (κ_P/κ_R=%.2f)\n", ρ, T, Q.kap(tab, ρ, T), Q.kapP(tab, ρ, T), Q.kapP(tab, ρ, T) / Q.kap(tab, ρ, T))
end

println("== 4. AD derivative of κ_R vs FD of the *interpolant* (consistency) ==")
for (ρ, T) in ((4e-11, 6241.0), (1e-10, 1.1e4), (1e-9, 3.3e4))
    g = ForwardDiff.gradient(v -> log(Q.kap(tab, exp(v[1]), exp(v[2]))), [log(ρ), log(T)])
    hh = 1e-6
    fd = [(log(Q.kap(tab, ρ * exp(hh), T)) - log(Q.kap(tab, ρ * exp(-hh), T))) / 2hh,
          (log(Q.kap(tab, ρ, T * exp(hh))) - log(Q.kap(tab, ρ, T * exp(-hh)))) / 2hh]
    @printf("  ρ=%.0e T=%.1e: ∂lnκ/∂lnρ AD=%.6f FD=%.6f | ∂lnκ/∂lnT AD=%.5f FD=%.5f\n", ρ, T, g[1], fd[1], g[2], fd[2])
end
t = @elapsed for _ in 1:20000; Q.kap(tab, 4e-11, 6241.0); end
@printf("== 5. kap bicubic: %.2f μs/call\n", 1e6 * t / 20000)
