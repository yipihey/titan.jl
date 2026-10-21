# S6a: summary of the structure dumps lib/popiii/s6b_<tag>_prof_<M>.csv (T_max and its mass coordinate, convective intervals, D depletion).
#   julia s6_profsum.jl TAG
using DelimitedFiles, Printf
tag = ARGS[1]; dir = joinpath(@__DIR__, "..", "lib", "popiii")
for f in sort(filter(x -> startswith(x, "s6b_$(tag)_prof_"), readdir(dir)))
    path = joinpath(dir, f); hdr = readline(path); D = readdlm(path, ','; skipstart = 2)
    m = D[:, 1]; T = D[:, 4]; conv = D[:, 8]; fD = D[:, 10]; r = D[:, 2]
    iT = argmax(T); iv = String[]; inr = false; m0 = 0.0
    for i in eachindex(m)
        if conv[i] > 0 && !inr; inr = true; m0 = i == 1 ? 0.0 : m[i-1]; end
        if conv[i] == 0 && inr; inr = false; push!(iv, @sprintf("%.3f-%.3f", m0, m[i-1])); end
    end
    inr && push!(iv, @sprintf("%.3f-1", m0))
    println(f, ": ", hdr[1:min(end, 120)])
    @printf("   T_max = %.3e K at m/M = %.3f (r = %.2f Rsun); convective m/M: %s; X_D/X_D0 min %.3f at m/M = %.3f\n", T[iT], m[iT], r[iT], join(iv, " "), minimum(fD), m[argmin(fD)])
end
