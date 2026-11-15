# Density dependence of the H2 line cooling per molecule: statistical equilibrium (H2Levels) against LTE, and the run's own zones
#   julia --project=<env> h2crit_check.jl <run directory (RUNREC)> <log n_c of the frame to check>
using Printf
include(joinpath(@__DIR__, "H2Levels.jl")); include(joinpath(@__DIR__, "runrec.jl"))
L = H2Levels.load(); yHe = 0.0823
permol(p) = H2Levels.cooling(L, p, 1.0)
println("A. optically thin, x(H2) = 1e-3, x(H+) = 1e-4: cooling per H2 molecule / LTE value")
@printf("%7s %11s |", "T [K]", "LTE erg/s"); for ln in 0:8; @printf(" n=1e%d ", ln); end; @printf("| n(half)\n")
for T in (200., 300., 500., 1000., 2000.)
    lte = permol(H2Levels.lte(L, T)); @printf("%7.0f %11.3e |", T, lte); nh = NaN; prev = 0.0
    f(n) = permol(H2Levels.populations(L, T, n * (1 - 2e-3), 1e-3n, yHe * n, 1e-4n)) / lte
    for ln in 0:8; @printf(" %6.4f", f(10.0^ln)); end
    lo, hi = 0.0, 9.0; for _ in 1:40; mid = (lo + hi) / 2; f(10^mid) < 0.5 ? (lo = mid) : (hi = mid); end
    @printf(" | %.2e\n", 10^lo)
end
K = lineage(ARGS[1]); x = K[argmin([abs(log10(y.ρ[1] / mH) - parse(Float64, ARGS[2])) for y in K])]; n = length(x.ρ)
@printf("\nB. run %s, step %d, log n_c %.2f: line cooling per molecule in the run's zones\n", ARGS[1], x.step, log10(x.ρ[1] / mH))
xH2 = sp(x, :H2); xHp = sp(x, :Hp)
@printf("%5s %9s %8s %9s | %10s %10s %10s | %8s %8s\n", "zone", "n_H", "T", "x(H2)", "run", "SE thin", "LTE", "run/SE", "run/LTE")
for j in unique(round.(Int, range(1, n, length = 26)))
    nH = x.ρ[j] * 0.75 / mH; nH2 = xH2[j] * nH; nH2 > 0 || continue
    run = x.Λext[j] / nH2
    se = permol(H2Levels.populations(L, x.T[j], nH * (1 - 2xH2[j] - xHp[j]), nH2, yHe * nH, xHp[j] * nH)); lte = permol(H2Levels.lte(L, x.T[j]))
    @printf("%5d %9.2e %8.1f %9.2e | %10.3e %10.3e %10.3e | %8.3f %8.3f\n", j, nH, x.T[j], xH2[j], run, se, lte, run / se, run / lte)
end
