# G2Opacity.means_fast vs means (1500 points): accuracy on a (ρ, T) grid in LTE, timing, ForwardDiff vs finite differences.
include(joinpath(@__DIR__, "load.jl")); include(joinpath(@__DIR__, "mats.jl"))
include(joinpath(JD, "G2Opacity.jl")); const G2O = G2Opacity
using ForwardDiff
kw(s, T) = (; n1 = s.nH, n2 = G2O.n2_lte(T, s.nH, s.ne, s.nH + s.nHe + s.nH2), nH2 = s.nH2, nHe = s.nHe, nHep = s.nHep, nHepp = s.nHepp, nHp = s.nHp, ne = s.ne, nHm = s.nHm, nHmS = s.nHm)
t0 = time(); G2O.fastgrid(); @printf("fast tables built in %.1f s, %d nodes\n", time() - t0, length(G2O.fastgrid().ν))
errs = [Float64[] for _ in 1:3]; worst = [(0.0, 0.0, 0.0) for _ in 1:3]
for lρ in -14.0:0.5:-1.0, lT in 3.0:0.05:5.45
    ρ = 10^lρ; T = 10^lT; s = net_lte(ρ, T)
    for hb in (true, false)
        a = G2O.means(ρ, T; kw(s, T)..., hbfP = hb); b = G2O.means_fast(ρ, T; kw(s, T)..., hbfP = hb)
        for q in 1:3
            e = abs(log(b[q] / a[q])); push!(errs[q], e); e > worst[q][1] && (worst[q] = (e, lρ, lT))
        end
    end
end
for (q, nm) in enumerate(("κ_E", "κ_P", "κ_R"))
    v = sort(errs[q]); @printf("%s: |Δ ln κ| median %.2e, 95%% %.2e, 99%% %.2e, max %.2e at log ρ %.1f log T %.2f\n", nm, v[end÷2], v[round(Int, 0.95end)], v[round(Int, 0.99end)], worst[q]...)
end
for (lρ, lT) in ((-6.66, 3.525), (-6.42, 3.98), (-4.0, 4.4))
    ρ = 10^lρ; T = 10^lT; s = net_lte(ρ, T); k = kw(s, T)
    G2O.means_fast(ρ, T; k..., hbfP = false)
    n = 200; t = time_ns(); for _ in 1:n; G2O.means_fast(ρ, T; k..., hbfP = false); end; tf = (time_ns() - t) / n / 1e3
    f(x) = collect(G2O.means_fast(ρ, x[1]; merge(k, (; ne = x[2]))..., hbfP = false))
    D = ForwardDiff.jacobian(f, [T, s.ne])
    t = time_ns(); for _ in 1:n; ForwardDiff.jacobian(f, [T, s.ne]); end; tj = (time_ns() - t) / n / 1e3
    h = 1e-6; Dfd = hcat((f([T * (1 + h), s.ne]) - f([T * (1 - h), s.ne])) / (2h * T), (f([T, s.ne * (1 + h)]) - f([T, s.ne * (1 - h)])) / (2h * s.ne))
    @printf("lρ %.2f lT %.2f: means_fast %.1f µs, with 2-partial Jacobian %.1f µs; max |AD/FD − 1| %.2e\n", lρ, lT, tf, tj, maximum(abs.(D ./ Dfd .- 1)))
end
