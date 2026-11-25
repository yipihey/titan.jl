# inproc.jl <jobfile> : several PopIIIRun runs in ONE julia process (scripts/dev/gates.sh --inproc; docs/POPIII_DEVTOOLS.md).
# jobfile line: "<run.out> KEY=VAL KEY=VAL ... -- <ARGS...>" (written by scripts/proto/gate_run.sh with GATE_PRINT=<jobfile>).
# Per job the process environment is reset to the one the process started with, the words are applied, then the stack defaults of
# popiii_stack.jl (as a fresh process would), stdout/stderr go to <run.out>, ARGS is set and the loader popiii_m2c.jl runs
# (first job: loads the modules or the TitanStack image, then PopIIIRun.main(ARGS)); later jobs call Main.PopIIIRun.main(ARGS) only.
# A run that ends in an error is recorded in its run.out and the next job starts (state leaks are then the question, see the docs).
# Appends "===== EXIT <rc> wall <s> =====" to every run.out like gate_run.sh. NOTE: with TITAN_PKG the image's baked load-time keys
# are those of the FIRST job; a later job that needs other baked keys (MGROUPS, HOT, ...) must go in another process.
const BASE_ENV = Dict{String,String}(ENV)
const TREE = haskey(ENV, "INPROC_TREE") ? ENV["INPROC_TREE"] : error("INPROC_TREE (the tree under test) is not set")
const LOADER = joinpath(TREE, "julia", "popiii_m2c.jl")
function setenv_job!(words)
    for k in collect(keys(ENV)); haskey(BASE_ENV, k) || delete!(ENV, k); end
    for (k, v) in BASE_ENV; ENV[k] = v; end
    for w in words; i = findfirst('=', w); ENV[w[1:i-1]] = w[i+1:end]; end
end
jobs = [split(l) for l in readlines(ARGS[1]) if !isempty(strip(l))]
first_job = true
for j in jobs
    out = String(j[1]); sep = findfirst(==("--"), j)
    words = String.(j[2:sep-1]); jargs = String.(j[sep+1:end])
    setenv_job!(words)
    if first_job
        include(joinpath(TREE, "julia", "popiii_stack.jl"))     # defaults for unset keys, as popiii_m2c.jl does
    else
        for (k, v) in Main.POPIII_STACK; haskey(ENV, k) || (ENV[k] = v); end
    end
    t0 = time(); rc = 0
    open(out, "w") do io
        redirect_stdio(stdout = io, stderr = io) do
            try
                empty!(ARGS); append!(ARGS, jargs)
                if first_job
                    Base.include(Main, LOADER)               # loads the stack, then PopIIIRun.main(ARGS)
                else
                    Core.eval(Main, :(PopIIIRun.main(ARGS)))
                end
            catch e
                rc = 1; showerror(stdout, e, catch_backtrace()); println()
            end
        end
        println(io, "===== EXIT ", rc, " wall ", round(time() - t0, digits = 3), " =====")
    end
    global first_job = false
    println("inproc: ", out, " rc ", rc, " wall ", round(time() - t0, digits = 1), " s")
end
