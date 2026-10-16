"""
    PortHarness

In-flight tier-(a) checks for Phase 5 translations. `wrap!(:update, TitanCore.update!)`
replaces `TitanDriver.ROUTINES[:update]` by a function that, on every call:
  1. snapshots all COMMON blocks (S0),
  2. runs the Fortran routine and records every COMMON variable (A) and the blocks (SA),
  3. restores S0, runs the Julia routine and records every variable (B),
  4. accumulates max |A−B| / max|A| per variable (integers: 1 if any element differs),
  5. restores SA, so the run continues on the Fortran trajectory (each call is isolated).
`FOLLOW_JULIA[] = true` instead continues on the Julia result (cumulative test).
"""
module PortHarness

using Printf
import ..TitanDriver
const TD = TitanDriver
const TF = TD.TF

const VARS = [(b, var) for b in sort(collect(keys(TF.BLOCKS))) for var in TF.BLOCKS[b] if var.T in (Float64, Int32)]
capture() = Dict{Symbol,Array}(var.name => copy(TF.var(var.name)) for (_, var) in VARS)

const STATS = Dict{Symbol,Dict{Symbol,Float64}}()     # routine => var => max relative difference
const WRITES = Dict{Symbol,Set{Symbol}}()             # routine => variables the Fortran actually changed
const NCALL = Dict{Symbol,Int}()
const FOLLOW_JULIA = Ref(false)

function vdiff(a::Array{Float64}, b::Array{Float64})
    s = 0.0; m = 0.0
    @inbounds for i in eachindex(a)
        s = max(s, abs(a[i] - b[i])); m = max(m, abs(a[i]))
    end
    s == 0 ? 0.0 : s / max(m, floatmin())
end
vdiff(a::Array{Int32}, b::Array{Int32}) = a == b ? 0.0 : 1.0

wrap!(s::Symbol, jfun::Function) = (TD.ROUTINES[s] = make_compare(s, () -> TD.fcall(Symbol(s, "_")), jfun))

"Comparison closure for any Fortran/Julia pair (e.g. routines called from inside other routines)."
function make_compare(s::Symbol, fort::Function, jfun::Function)
    st = get!(STATS, s, Dict{Symbol,Float64}()); wr = get!(WRITES, s, Set{Symbol}())
    return function ()
        S0 = TD.snapshot(); P0 = capture()
        fort(); A = capture(); SA = TD.snapshot()
        TD.restore!(S0)
        jfun(); B = capture()
        for (n, a) in A
            d = vdiff(a, B[n]); st[n] = max(get(st, n, 0.0), d)
            a != P0[n] && push!(wr, n)
        end
        NCALL[s] = get(NCALL, s, 0) + 1
        FOLLOW_JULIA[] || TD.restore!(SA)
        nothing
    end
end

"Print, per wrapped routine, the variables with the largest differences and any written-but-differing ones."
function report(io = stdout; top = 8, tol = 1e-12)
    for s in sort(collect(keys(STATS)))
        st = STATS[s]; wr = WRITES[s]
        bad = sort([(d, n) for (n, d) in st if d > 0]; rev = true)
        nbad = count(x -> x[1] > tol, bad)
        @printf(io, "%-8s calls %6d | Fortran writes %3d variables | %3d differ > %.0e | worst: %s\n",
                s, get(NCALL, s, 0), length(wr), nbad, tol,
                isempty(bad) ? "none (identical)" : join([@sprintf("%s %.1e", n, d) for (d, n) in bad[1:min(top, end)]], ", "))
    end
end

end # module
