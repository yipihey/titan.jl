# Turn a `julia --trace-compile=<file> --trace-compile-timing popiii_m2c.jl ...` trace into src/pc.jl, the precompile statements that
# TitanStack applies when it is built (statements that do not resolve inside the package, e.g. those naming the driver's own
# functions in Main, are skipped when the package is built).
#   julia julia/TitanStack/trace_to_pc.jl <trace.jl> [<trace2.jl> ...]  > julia/TitanStack/src/pc.jl
seen = Set{String}()
for f in ARGS, l in eachline(f)
    s = replace(l, r"^#=.*?=#\s*" => "")
    startswith(s, "precompile(") || continue
    occursin("Main.", s) && continue
    s in seen && continue
    push!(seen, s); println(s)
end
