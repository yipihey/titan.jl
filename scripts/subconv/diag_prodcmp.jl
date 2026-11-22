# Step 0: bit-identity of replay snapshots vs production snapshots at common steps. julia diag_prodcmp.jl <prodsnaps.jls> <replaysnaps.jls>
using Serialization
function rd(f); d = Dict{Int,Any}(); open(f) do io; while !eof(io); s = deserialize(io); d[s.step] = s; end; end; d; end
function main()
A = rd(ARGS[1]); B = rd(ARGS[2]); common = sort(collect(intersect(keys(A), keys(B))))
nf = 0; nb = 0; bad = Int[]
for st in common
    a = A[st]; b = B[st]
    for k in keys(a)
        (k === :cd || !haskey(b, k)) && continue
        nf += 1
        if !isequal(a[k], b[k]); nb += 1; push!(bad, st); println("DIFF step $st field $k"); end
    end
end
println("production steps ", length(A), " (", minimum(keys(A)), "..", maximum(keys(A)), "); replay steps ", length(B), " (", minimum(keys(B)), "..", maximum(keys(B)), "); common ", length(common),
        "; fields compared ", nf, "; fields differing ", nb, "; distinct steps differing ", length(unique(bad)))
println("fields of replay not in production: ", setdiff(keys(first(values(B))), keys(first(values(A)))), "; production-only: ", setdiff(keys(first(values(A))), keys(first(values(B)))))

end
main()