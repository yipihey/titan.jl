include(joinpath(@__DIR__, "load.jl")); include(joinpath(@__DIR__, "mats.jl"))
include(joinpath(JD, "G2Opacity.jl")); const G2O = G2Opacity
function tm(f, n)
    f(); t = time_ns(); for _ in 1:n; f(); end; (time_ns() - t) / n / 1e6
end
for (lρ, lT) in ((-6.657, 3.525), (-6.418, 3.983), (-4.0, 4.4), (-9.0, 3.6))
    ρ = 10^lρ; T = 10^lT; s = net_lte(ρ, T); n2 = G2O.n2_lte(T, s.nH, s.ne, s.nH + s.nHe + s.nH2)
    a = tm(() -> net_means(ρ, T, s), 20)
    b = tm(() -> G2O.means(ρ, T; n1 = s.nH, n2, nH2 = s.nH2, nHe = s.nHe, nHep = s.nHep, nHepp = s.nHepp, nHp = s.nHp, ne = s.ne, nHm = s.nHm, nHmS = s.nHm, hbfP = true), 20)
    c = tm(() -> tab_kap(ρ, T), 2000); d = tm(() -> EQ.eos(ρ, T), 2000); e = tm(() -> net_eos(ρ, T, s), 2000)
    f = tm(() -> net_lte(ρ, T), 50)
    @printf("lρ %.2f lT %.2f: nlte_means-equivalent %.2f ms, G2Opacity.means %.2f ms, table κ %.4f ms, table EOS %.4f ms, network EOS %.4f ms, network Saha solve %.3f ms\n", lρ, lT, a, b, c, d, e, f)
end
