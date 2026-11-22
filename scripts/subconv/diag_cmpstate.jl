# relative state differences of two snapshot streams at a few common steps: julia diag_cmpstate.jl <A.jls> <B.jls>
using Serialization
rd(f) = (d = Dict{Int,Any}(); open(f) do io; while !eof(io); s = deserialize(io); d[s.step] = s; end; end; d)
A = rd(ARGS[1]); B = rd(ARGS[2])
for st in sort(collect(intersect(keys(A), keys(B))))[[1,2,5,10,20,40]]
    a = A[st]; b = B[st]
    f(x) = maximum(abs.((x .- b.T) ./ b.T))
    println(st, " t ", a.t, " ", b.t, " dt ", a.dt, " ", b.dt, " maxrelT ", maximum(abs.(a.T ./ b.T .- 1)), " maxrel rho ", maximum(abs.(a.ρ ./ b.ρ .- 1)), " maxrel u ", maximum(abs.(a.u .- b.u)) / maximum(abs.(a.u)), " nz ", length(a.T), length(b.T))
end
