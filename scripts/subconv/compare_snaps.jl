# Bit-identity of two RUNREC snapshot streams: julia compare_snaps.jl <dirA> <dirB>
# Every field of A's snapshots (all solver fields) must be isequal (bitwise for Float64; NaN == NaN, -0.0 != 0.0) to B's; fields only B has (cd) are listed, not compared.
using Serialization
function rd(dir)
    f = filter(x -> occursin(r"^snaps_\d+_\d+\.jls$", x), readdir(dir)); length(f) == 1 || error("expected one snaps file in $dir: $f")
    v = Any[]; open(joinpath(dir, f[1])) do io; while !eof(io); push!(v, deserialize(io)); end; end; v
end
function main()
    A = rd(ARGS[1]); B = rd(ARGS[2]); ok = true
    length(A) == length(B) || (println("different number of snapshots: ", length(A), " vs ", length(B)); ok = false)
    nf = 0; nb = 0
    for (a, b) in zip(A, B)
        for k in keys(a)
            k === :cd && continue
            nf += 1
            if !isequal(a[k], b[k]); ok = false; nb += 1; println("DIFFERENT step ", a.step, " field ", k); end
        end
    end
    println("snapshots compared: ", min(length(A), length(B)), " (steps ", join([a.step for a in A], ","), "); fields compared: ", nf, "; fields with differences: ", nb,
            "; extra fields in B: ", join(setdiff(keys(B[1]), keys(A[1])), ","))
    println(ok ? "BIT-IDENTICAL" : "NOT IDENTICAL")
end
main()
