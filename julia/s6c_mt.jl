# S6b Part 1: M_pt(t) of a sink run (lib/popiii/s6c_<TAG>.csv): point mass and sink radius at given times after the restart.
#   julia s6c_mt.jl <TAG> [t list in yr]
using Printf
tag = ARGS[1]; ts = length(ARGS) > 1 ? parse.(Float64, ARGS[2:end]) : [10, 30, 100, 300, 1000, 3000, 10000, 20000]
L = readlines(joinpath(@__DIR__, "..", "lib", "popiii", "s6c_$(tag).csv"))[2:end]
rows = [parse.(Float64, split(l, ",")[1:7]) for l in L if !isempty(strip(l))]
t0 = rows[1][2]
println(tag, ": point mass of the sink vs time after the first row (t0 = ", t0, " yr)")
for T in ts
    i = findfirst(r -> r[2] - t0 >= T, rows)
    i === nothing && continue
    @printf("  t = %7g yr: M_pt = %9.3f Msun  r_s = %8.1f AU   (step %d)\n", T, rows[i][6], rows[i][5], rows[i][1])
end
@printf("  last row: t = %.1f yr, M_pt = %.3f Msun, step %d\n", rows[end][2] - t0, rows[end][6], rows[end][1])
