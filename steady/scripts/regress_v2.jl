push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; using Printf
tab = load_table(); M = 1e6Msun; L = LEdd_es(M)
ref = Dict((0.1, 0.0) => 8.470e15, (0.01, 0.0) => 1.349e16, (0.1, 1.0) => (8.273e14, 8.420e15), (0.01, 1.0) => (6.519e14, 1.347e16))
for fBH in (0.1, 0.01)
    p = QSParams(M = M, MBH = fBH * M, L = L, tab = tab)
    s = solve_R(p; lRmin = 15.6, lRmax = 16.4, nscan = 33)
    @printf("fBH=%.2f Ṁ=0: R=%.4e (v1 %.4e, rel %.2e) Teff=%.0f\n", fBH, s.R, ref[(fBH, 0.0)], s.R / ref[(fBH, 0.0)] - 1, s.bc.Tb)
    # both branches at Ṁ = 1
    for (lo, hi, k) in ((14.6, 15.2, 1), (15.8, 16.3, 2))
        s = solve_R(p; Mdot = 1.0Msun / yr, lRmin = lo, lRmax = hi, nscan = 25)
        r = ref[(fBH, 1.0)][k]
        @printf("fBH=%.2f Ṁ=1 branch %d: R=%.4e (v1 %.4e, rel %.2e) Tb=%.0f\n", fBH, k, s.R, r, s.R / r - 1, s.bc.Tb)
    end
end
