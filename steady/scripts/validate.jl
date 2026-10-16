push!(LOAD_PATH, joinpath(@__DIR__, "..", "src")); using QuasiStarInfall; const Q = QuasiStarInfall
tab = load_table()
M = 1e6Msun; MBH = 0.1M; L = LEdd_es(M)
p = QSParams(M = M, MBH = MBH, L = L, tab = tab)
println("L = ", L/Lsun, " Lsun")
# inspect sign scan for no infall
for lR in 15.0:0.25:17.0
    s, mrin, bc = shoot_sign(10^lR, p)
    println("logR=$lR sign=$s m(rin)/M=$(round(mrin/M,sigdigits=4)) Teff=$(bc===nothing ? NaN : round(bc.Tb)) Pb=$(bc===nothing ? NaN : round(bc.Pb,sigdigits=3)) Γκ=$(bc===nothing ? NaN : round(bc.Γκ,sigdigits=3))")
end
