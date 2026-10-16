# Refine the saddle-node (fold) where the extended and compact branches merge.
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall
using Printf
tab = load_table(); M = 1e6Msun; L = LEdd_es(M)
lRs = collect(15.2:0.01:16.25)
open(joinpath(@__DIR__, "..", "output", "fold.csv"), "w") do out
    println(out, "fBH,Mdot,nroots,lR1,lR2")
    for (fBH, mds) in ((0.1, [6.0, 7.0, 8.0, 9.0]), (0.01, [12.0, 15.0, 20.0, 25.0]))
        p = QSParams(M = M, MBH = fBH * M, L = L, tab = tab)
        for md in mds
            sg = [shoot_sign(10^lR, p; Mdot = md * Msun / yr)[1] for lR in lRs]
            idx = [i for i in 1:length(lRs)-1 if sg[i] != 0 && sg[i+1] != 0 && sg[i] != sg[i+1]]
            r = [0.5(lRs[i] + lRs[i+1]) for i in idx]
            @printf(out, "%g,%g,%d,%s,%s\n", fBH, md, length(r), length(r) > 0 ? @sprintf("%.3f", r[1]) : "NaN", length(r) > 1 ? @sprintf("%.3f", r[end]) : "NaN")
            @printf("fBH=%.2f Ṁ=%5.1f: %d roots %s\n", fBH, md, length(r), join(r, ", ")); flush(stdout); flush(out)
        end
    end
end
