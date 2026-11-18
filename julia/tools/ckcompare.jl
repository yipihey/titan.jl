# Compare two checkpoints (<base>.tstate.gz + <base>.aux) of popiii_m2c.jl: max relative difference per quantity.
#   julia ckcompare.jl baseA baseB [ngrs-range]     (stdlib only)
using Serialization, Printf
fmt(x) = @sprintf("%.3e", x)
function readstate(path)
    io = IOBuffer(read(`gzip -dc $path`)); String(read(io, 8)) == "TITANST1" || error("not TITANST1")
    n = read(io, Int64); st = Dict{Symbol,Array}()
    for _ in 1:n
        name = Symbol(String(read(io, read(io, Int64)))); code = read(io, UInt8)
        dims = Tuple(read(io, Int64) for _ in 1:read(io, Int64))
        T = code == 0x01 ? Float64 : code == 0x02 ? Int32 : NTuple{80,UInt8}
        st[name] = reshape(collect(reinterpret(T, read(io, prod(dims; init = 1) * sizeof(T)))), dims)
    end
    st
end
function rel(x::AbstractArray{<:Real}, y::AbstractArray{<:Real})
    size(x) == size(y) || return (Inf, 0)
    m = 0.0; km = 0
    for i in eachindex(x)
        (isnan(x[i]) && isnan(y[i])) && continue
        d = abs(x[i] - y[i]); s = max(abs(x[i]), abs(y[i])); r = s == 0 ? 0.0 : d / s
        r > m && (m = r; km = i)
    end
    (m, km)
end
flat(x::Number) = [Float64(x)]
flat(x::AbstractArray{<:Number}) = vec(Float64.(x))
flat(x::AbstractArray) = reduce(vcat, [flat(e) for e in x]; init = Float64[])
flat(x::Tuple) = reduce(vcat, [flat(e) for e in x]; init = Float64[])
flat(x::NamedTuple) = reduce(vcat, [flat(e) for e in values(x)]; init = Float64[])
flat(x::Dict) = reduce(vcat, [flat(x[k]) for k in sort(collect(keys(x)); by = string)]; init = Float64[])
flat(x::Bool) = [Float64(x)]
flat(x) = Float64[]
function walk!(out, pa, pb, name)
    if pa isa NamedTuple
        for k in keys(pa); hasproperty(pb, k) && walk!(out, getfield(pa, k), getfield(pb, k), "$name.$k"); end
    else
        fa, fb = flat(pa), flat(pb); length(fa) == length(fb) || (push!(out, (name, Inf, 0)); return)
        isempty(fa) || push!(out, (name, rel(fa, fb)...))
    end
end
function main(A, B)
    sa, sb = readstate(A * ".tstate.gz"), readstate(B * ".tstate.gz"); res = []
    for k in sort(collect(keys(sa)))
        eltype(sa[k]) === Float64 || continue
        r, i = rel(sa[k], sb[k]); push!(res, ("tstate." * string(k), r, i))
    end
    try   # (aux of older checkpoints holds instances of the driver's own types: not readable here -> tstate only)
        aa, ab = deserialize(A * ".aux"), deserialize(B * ".aux"); walk!(res, aa, ab, "aux")
    catch e
        println("(aux not compared: ", first(sprint(showerror, e), 80), ")")
    end
    sort!(res; by = x -> -x[2])
    prim = ("rn", "dn", "tn", "egn", "ern", "frn", "un", "xmen", "pgn", "timen", "timeo", "dtime", "jstep")
    println("PRIMARY state (relative difference, linear index = zone for 1-D arrays):")
    for (n, r, i) in sort(filter(x -> (startswith(x[1], "tstate.") && x[1][8:end] in prim) || x[1] in ("aux.mat.x", "aux.k", "aux.mat.Tref", "aux.mat.Y"), res); by = x -> x[1])
        println("  ", rpad(n, 24), fmt(r), "  idx ", i)
    end
    println("ALL (sorted):")
    println("max relative difference per quantity (nonzero only), sorted; index = linear index in the array (zone = first index for COMMON arrays)")
    for (n, r, i) in res; r > 0 && println(rpad(n, 40), fmt(r), "  idx ", i); end
    println("overall max: ", isempty(res) ? 0.0 : res[1][2], " (", isempty(res) ? "" : res[1][1], "); quantities compared: ", length(res), ", identical: ", count(x -> x[2] == 0, res))
end

main(ARGS[1], ARGS[2])
