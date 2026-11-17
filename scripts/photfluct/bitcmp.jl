# bitcmp.jl <dirA> <dirB> : compare snapshots at common steps (max |rel diff| of r, T, ρ, E_r, F_r, t)
using Serialization, Printf
isdefined(Main, :NucNet) || (@eval Main module NucNet end)
include(joinpath(@__DIR__, "..", "..", "julia", "runrec.jl"))
A = Dict(x.step => x for x in lineage(ARGS[1])); B = Dict(x.step => x for x in lineage(ARGS[2]))
rd(a, b) = length(a) == length(b) ? maximum(abs.(a .- b) ./ max.(abs.(a), abs.(b), 1e-300)) : Inf
for s in sort(collect(intersect(keys(A), keys(B))))
    a = A[s]; b = B[s]
    @printf("step %d  t %s  dt %s  r %.2e T %.2e rho %.2e Er %.2e Fr %.2e\n", s, a.t == b.t ? "==" : @sprintf("%.3e", (a.t - b.t) / a.t), a.dt == b.dt ? "==" : "!=",
            rd(a.r, b.r), rd(a.T, b.T), rd(a.ρ, b.ρ), rd(a.er, b.er), rd(a.fr, b.fr))
end
