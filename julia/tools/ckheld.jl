# Compare two .held dumps (HELDDUMP=1): which module-level globals differ.   julia ckheld.jl a.held b.held
using Serialization, Printf
A, B = deserialize(ARGS[1]), deserialize(ARGS[2]); res = []
for k in sort(collect(union(keys(A), keys(B))))
    (haskey(A, k) && haskey(B, k)) || (push!(res, (k, Inf, 0, "missing in " * (haskey(A, k) ? "B" : "A"))); continue)
    a, b = A[k], B[k]
    length(a) == length(b) || (push!(res, (k, Inf, 0, "length $(length(a)) vs $(length(b))")); continue)
    m = 0.0; km = 0
    for i in eachindex(a)
        (isnan(a[i]) && isnan(b[i])) && continue; d = abs(a[i] - b[i]); s = max(abs(a[i]), abs(b[i])); r = s == 0 ? 0.0 : d / s
        r > m && (m = r; km = i)
    end
    m > 0 && push!(res, (k, m, km, ""))
end
sort!(res; by = x -> -x[2])
for (k, m, i, msg) in res; @printf("%-60s %.3e idx %d %s\n", k, m, i, msg); end
println(length(res), " globals differ of ", length(A))
