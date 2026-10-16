push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall
tab = load_table()
M = 1e6Msun; L = LEdd_es(M)
for fBH in (0.1, 0.01), ε in (1e-4, 1e-3, 1e-2)
    p = QSParams(M = M, MBH = fBH*M, L = L, ε_in = ε, tab = tab)
    s = solve_R(p; lRmin = 15.6, lRmax = 16.4, nscan = 33)
    s === nothing && (println("fBH=$fBH ε=$ε: no solution"); continue)
    pr = s.prof; g = G*M/s.R^2
    imax = argmax(pr[:v] ./ pr[:cs])
    println("fBH=$fBH ε=$ε: R=$(round(s.R,sigdigits=4)) cm = $(round(s.R/Rsun,sigdigits=3)) Rsun, Teff=$(round(s.bc.Tb)) K, log g=$(round(log10(g),digits=3)), " *
            "max Mach=$(round(maximum(pr[:v] ./ pr[:cs]),sigdigits=3)) at T=$(round(pr[:T][imax],sigdigits=3)), max Γ=$(round(maximum(pr[:Γ]),sigdigits=3)), Γκ_phot=$(round(s.bc.Γκ,sigdigits=3)), Pgas_b=$(round(s.bc.Pgas_b,sigdigits=3))")
end
println("MESA-QUEST (paper): Teff = 6241 K, log g = -0.003, R = 1.15e16 cm (from L,Teff)")
