# cmp_snaps.jl <dirA> <dirB> [after_step] <prefix> : bitwise (isequal, all fields but :cd) comparison of the snapshots of two RUNREC dirs that have the SAME step
# in common (steps > after_step); prints "<prefix>.snaps_compared N", "<prefix>.fields_compared N", "<prefix>.fields_diff N" and DIFFERENT lines on stderr.
# (scripts/subconv/compare_snaps.jl needs equal snapshot counts; a restarted run has snapshots of other steps too.)
using Serialization
function rd(dir)
    v = Dict{Int,Any}()
    for f in sort(filter(x -> occursin(r"^snaps_\d+_\d+\.jls$", x), readdir(dir)))
        open(joinpath(dir, f)) do io; while !eof(io); s = deserialize(io); v[s.step] = s; end; end
    end
    v
end
A = rd(ARGS[1]); B = rd(ARGS[2]); after = parse(Int, ARGS[3]); pre = ARGS[4]
common = sort([k for k in keys(A) if haskey(B, k) && k > after])
nf = 0; nb = 0
for k in common
    for f in keys(A[k])
        f === :cd && continue
        global nf += 1
        if !(haskey(B[k], f) && isequal(A[k][f], B[k][f])); global nb += 1; println(stderr, "DIFFERENT step ", k, " field ", f); end
    end
end
println(pre, ".snaps_compared ", length(common)); println(pre, ".fields_compared ", nf); println(pre, ".fields_diff ", nb)
