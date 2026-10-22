"""
    Perf

Lightweight performance recording for the Pop III runs (no dependencies). Every driver records the same way, so runs
can be compared later when looking for room for improvement.

    Perf.reset!(run = "z19_e2e_collapse"; stage = "collapse", meta...)    # start of a run
    Perf.@section :residual begin ... end    # inclusive wall time, calls, allocated bytes per named section
    Perf.count!(:newton_iter, 3)             # integer counters (iterations, failures, retries, ...)
    Perf.sample!(step = k, t = t, dt = dt, zones = n, newton = it)   # optional per-step row (keep it sparse: every N steps)
    Perf.write()                             # lib/popiii/perf/<run>_{sections,counters,samples}.csv + <run>_meta.txt
    Perf.summary()                           # table to stdout

Sections may nest; times are inclusive. With threads, call `@section` only around serial code or around the whole
threaded loop (the accumulators are not thread-safe).
"""
module Perf

using Printf

mutable struct Sec
    calls::Int; ns::UInt64; bytes::Int; first_ns::UInt64     # first_ns: the very first call (includes JIT compilation of everything it runs); also part of ns
end
const SECS = Dict{Symbol,Sec}()
const ORDER = Symbol[]
const COUNTS = Dict{Symbol,Int}()
const SAMPLES = Vector{Vector{Pair{Symbol,Any}}}()
const META = Dict{String,String}()
const RUN = Ref("run")
const T0 = Ref(time_ns())
const DIR = Ref(joinpath(@__DIR__, "..", "lib", "popiii", "perf"))

gitrev() = try
    strip(read(`git -C $(@__DIR__) rev-parse --short HEAD`, String))
catch
    "?"
end

"Start a run: clears all accumulators and records machine/run metadata (extra keyword pairs are stored as metadata)."
function reset!(run::AbstractString = "run"; kw...)
    empty!(SECS); empty!(ORDER); empty!(COUNTS); empty!(SAMPLES); empty!(META)
    RUN[] = run; T0[] = time_ns()
    META["run"] = run
    META["date"] = string(Base.Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
    META["host"] = gethostname()
    META["cpu"] = Sys.CPU_NAME * " ($(Sys.CPU_THREADS) hw threads)"
    META["julia"] = string(VERSION)
    META["julia_threads"] = string(Threads.nthreads())
    META["git"] = gitrev()
    for (k, v) in kw; META[string(k)] = string(v); end
    nothing
end
"Add or overwrite a metadata entry."
meta!(k, v) = (META[string(k)] = string(v); nothing)

"Accumulate one timed call of section `name`."
function add!(name::Symbol, ns, bytes)
    s = get!(SECS, name) do
        push!(ORDER, name); Sec(0, 0, 0, 0)
    end
    s.calls == 0 && (s.first_ns = ns)
    s.calls += 1; s.ns += ns; s.bytes += bytes
    nothing
end

"""
    @section :name expr

Time `expr` (wall clock and allocated bytes) and accumulate under `:name`; returns the value of `expr`.
"""
macro section(name, ex)
    quote
        local _b0 = Base.gc_bytes()
        local _t0 = time_ns()
        local _v = $(esc(ex))
        add!($(esc(name)), time_ns() - _t0, Base.gc_bytes() - _b0)
        _v
    end
end

"Increment counter `name` by `n`."
count!(name::Symbol, n::Integer = 1) = (COUNTS[name] = get(COUNTS, name, 0) + n; nothing)

"Append a per-step sample (keyword pairs; columns are the union over all samples)."
sample!(; kw...) = (push!(SAMPLES, [k => v for (k, v) in kw]); nothing)

"Total wall time since `reset!` [s]."
elapsed() = (time_ns() - T0[]) * 1e-9

"Write the CSV files and the metadata to `dir` (default lib/popiii/perf). Returns the file prefix."
function write(; dir = DIR[])
    mkpath(dir)
    pre = joinpath(dir, RUN[])
    META["wall_s"] = @sprintf("%.3f", elapsed())
    open(pre * "_meta.txt", "w") do io
        for k in sort(collect(keys(META))); println(io, k, " = ", META[k]); end
    end
    open(pre * "_sections.csv", "w") do io
        println(io, "section,calls,total_s,mean_ms,frac_of_wall,alloc_MB,first_call_s,mean_ms_excl_first")
        tw = elapsed()
        for n in ORDER
            s = SECS[n]
            @printf(io, "%s,%d,%.6f,%.6f,%.5f,%.3f,%.6f,%.6f\n", n, s.calls, s.ns * 1e-9, s.ns * 1e-6 / max(s.calls, 1), s.ns * 1e-9 / tw, s.bytes / 2^20, s.first_ns * 1e-9,
                    (s.ns - s.first_ns) * 1e-6 / max(s.calls - 1, 1))
        end
    end
    open(pre * "_counters.csv", "w") do io
        println(io, "counter,value")
        for k in sort(collect(keys(COUNTS))); println(io, k, ",", COUNTS[k]); end
    end
    if !isempty(SAMPLES)
        cols = Symbol[]
        for r in SAMPLES, (k, _) in r; k in cols || push!(cols, k); end
        open(pre * "_samples.csv", "w") do io
            println(io, join(cols, ","))
            for r in SAMPLES
                d = Dict(r)
                println(io, join([haskey(d, c) ? string(d[c]) : "" for c in cols], ","))
            end
        end
    end
    pre
end

"Print a summary table (sections sorted by total time, counters) to `io`."
function summary(io::IO = stdout)
    tw = elapsed()
    @printf(io, "perf %s: wall %.1f s, %s threads, git %s\n", RUN[], tw, get(META, "julia_threads", "?"), get(META, "git", "?"))
    @printf(io, "  %-24s %10s %12s %10s %8s %10s %10s\n", "section", "calls", "total [s]", "mean [ms]", "% wall", "alloc [MB]", "1st [s]")
    for n in sort(ORDER; by = n -> -SECS[n].ns)
        s = SECS[n]
        @printf(io, "  %-24s %10d %12.3f %10.3f %8.1f %10.1f %10.3f\n", n, s.calls, s.ns * 1e-9, s.ns * 1e-6 / max(s.calls, 1), 100 * s.ns * 1e-9 / tw, s.bytes / 2^20, s.first_ns * 1e-9)
    end
    for k in sort(collect(keys(COUNTS))); @printf(io, "  counter %-20s %d\n", k, COUNTS[k]); end
    nothing
end

end # module
