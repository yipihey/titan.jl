# Compare two RUNREC snapshot streams step by step (bit-identity checks of default-off switches).
#   julia --project=<env> snapcmp.jl <snaps_a.jls> <snaps_b.jls>
using Serialization, Printf
include(joinpath(@__DIR__, "runrec.jl"))
A = readstream(ARGS[1]); B = readstream(ARGS[2])
ib = Dict(x.step => x for x in B)
const FIELDS = (:r, :u, :ρ, :T, :e, :pg, :er, :fr, :chi, :kp, :Y)
nsame = 0; ncmp = 0
for x in A
    y = get(ib, x.step, nothing); y === nothing && continue
    global ncmp += 1
    d = Dict{Symbol,Float64}()
    for f in FIELDS
        a = getproperty(x, f); b = getproperty(y, f)
        size(a) == size(b) || (d[f] = Inf; continue)
        d[f] = maximum(abs.(a .- b) ./ (abs.(a) .+ abs.(b) .+ 1e-300); init = 0.0)
    end
    idt = x.t == y.t && x.dt == y.dt && all(v -> v == 0, values(d))
    idt && (global nsame += 1)
    @printf("step %d: t %s dt %s  max rel diff %s\n", x.step, x.t == y.t ? "=" : "≠", x.dt == y.dt ? "=" : "≠",
            join(["$f $(round(d[f], sigdigits = 2))" for f in FIELDS], " "))
end
println("$nsame of $ncmp common snapshots bit-identical")
