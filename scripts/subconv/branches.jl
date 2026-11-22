using Serialization
f = ARGS[1]; want = parse.(Int, split(ARGS[2], ","))
io = open(f)
while !eof(io)
    s = try deserialize(io) catch; break end
    s.step in want || continue
    n = length(s.ρ); T = s.T; u = s.u; cs = sqrt.(s.pg ./ s.ρ)
    uc = 0.5 .* (u[1:end-1] .+ u[2:end])
    br = [T[k] > 3e5 ? "H" : (s.eq[k] > 0.5 ? (T[k] > 1e4 ? "Q+" : "Q") : "N") for k in 1:n]
    # shock: most negative du between zones
    ks = argmin(diff(uc))
    println("step $(s.step) t=$(s.t)  nz=$n  counts H=$(count(==("H"),br)) Q+(table,T>1e4)=$(count(==("Q+"),br)) Q=$(count(==("Q"),br)) N=$(count(==("N"),br)); shock near zone $ks")
    for k in max(1,ks-8):min(n,ks+3)
        println("  k=$k r=$(round(s.r[k]/6.96e10,sigdigits=4)) Rsun T=$(round(T[k],sigdigits=3)) rho=$(round(s.ρ[k],sigdigits=3)) nH=$(round(s.ρ[k]/1.67e-24*0.76,sigdigits=3)) |u|/cs=$(round(abs(uc[k])/cs[k],sigdigits=2)) branch=$(br[k])")
    end
end
