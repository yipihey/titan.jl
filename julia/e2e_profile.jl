# E2E-A: statistical profile of a driver run, attributed to the source file of the innermost known frame (EOS, opacity, network, residual, solver ...).
#   E2E_SCRIPT=s5_run.jl [driver environment ...] TITAN_NOLIB=1 julia -t 4 --project=env e2e_profile.jl
# The driver runs inside Profile.@profile (include); output: table on stdout and lib/popiii/perf/<E2E_PROFOUT>_profile.txt. Wall-clock sampling of the main thread only (interval 5 ms).
using Profile, Printf
script = joinpath(@__DIR__, ENV["E2E_SCRIPT"])
Profile.init(n = 10^7, delay = 0.005)
@profile try include(script) catch e; println("driver stopped: ", sprint(showerror, e)[1:min(end, 200)]) end
data, lidict = Profile.retrieve()
cats = ["StellarEOS", "StellarOpacity", "ContinuumOpacity", "NucNet", "StarMaterial", "PopIIIMaterial", "ChemNet", "AtomicHHe", "H2Lines", "H2Levels", "BlockSolver", "TitanResidFast", "TitanCore", "TitanDriver", "Trbdf2", "Z0Physics", "TitanMaterials", "ProblemSetup"]
counts = Dict{String,Int}(); total = 0; fd = 0
# the samples are separated by zeros; frames listed leaf first
i = 1; n = length(data)
stack = UInt64[]
function attribute(stack)
    for ip in stack
        for fr in get(lidict, ip, [])
            f = string(fr.file)
            for c in cats
                occursin(c * ".jl", f) && return c
            end
        end
    end
    for ip in stack, fr in get(lidict, ip, [])
        occursin(basename(script), string(fr.file)) && return "driver script (" * basename(script) * ")"
    end
    nothing                          # idle threads (scheduler, GC wait): not counted
end
for k in 1:n
    if data[k] == 0
        if !isempty(stack)
            c = attribute(stack); c === nothing || (global total += 1; counts[c] = get(counts, c, 0) + 1)
        end
        empty!(stack)
    else
        push!(stack, data[k])
    end
end
# leaf function names within the module named by E2E_PROFDETAIL (default TitanDriver)
det = Dict{String,Int}(); dmod = get(ENV, "E2E_PROFDETAIL", "TitanDriver")
let st = UInt64[]
    for k in 1:n
        if data[k] == 0
            for ip in st
                hit = false
                for fr in get(lidict, ip, [])
                    if occursin(dmod * ".jl", string(fr.file)); key = string(fr.func) * ":" * string(fr.line); det[key] = get(det, key, 0) + 1; hit = true; break; end
                end
                hit && break
            end
            empty!(st)
        else
            push!(st, data[k])
        end
    end
end
out = String[]
push!(out, @sprintf("profile of %s: %d samples (main thread, ~5 ms)", ENV["E2E_SCRIPT"], total))
for (c, v) in sort(collect(counts); by = x -> -x[2])
    push!(out, @sprintf("  %-50s %6d  %5.1f %%", c, v, 100 * v / max(total, 1)))
end
push!(out, "leaf functions in " * dmod * ".jl:")
for (c, v) in first(sort(collect(det); by = x -> -x[2]), 12); push!(out, @sprintf("  %-50s %6d", c, v)); end
println(join(out, "\n"))
mkpath(joinpath(@__DIR__, "..", "lib", "popiii", "perf"))
write(joinpath(@__DIR__, "..", "lib", "popiii", "perf", get(ENV, "E2E_PROFOUT", "profile") * "_profile.txt"), join(out, "\n") * "\n")
