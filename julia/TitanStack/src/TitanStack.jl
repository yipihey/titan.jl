# TitanStack: the Pop III stack's module files (julia/*.jl) as one precompiled package, so that a run does not compile them again.
# The module files stay where they are and are included here in the order popiii_m2c.jl includes them; `popiii_m2c.jl` loads this
# package instead of including them when TITAN_PKG is set (default off). Nothing about the physics or numerics differs.
#
# What the image freezes: all load-time state of the modules (`const X = get(ENV, ...)`, conditional includes, ...). envguard.jl finds
# the environment keys that matter mechanically from the sources and records their values at build time (BAKED); `__init__` compares
# them with the loading process's environment and refuses to load on any difference (the driver then runs on the include path).
module TitanStack

const SRC = normpath(@__DIR__, "..", "..")      # julia/: where the module files are

# Names that the precompile statements (pc.jl: a --trace-compile trace) use unqualified; without them the statements that mention
# ForwardDiff / Printf / StaticArrays types fail to resolve and are skipped (about half of them).
import ForwardDiff, Printf, StaticArrays, Serialization, LinearAlgebra, Statistics, Libdl
const StaticArraysCore = parentmodule(StaticArrays.SArray)
const DiffResults = ForwardDiff.DiffResults
include("envguard.jl")

# --- conditional parts, exactly the conditions of popiii_m2c.jl (the guard records MGROUPS, HOT and LINEGPU through these lines) -----
const WITH_LINEGPU = haskey(ENV, "LINEGPU")                                          # LINEGPU=cuda|metal|cpu includes h2lines_gpu.jl (KernelAbstractions): not part of the image
WITH_LINEGPU && error("TitanStack: LINEGPU is set; the formal-solution GPU path (h2lines_gpu.jl) is not part of the package image, use the include path")
const WITH_MGROUPS = haskey(ENV, "MGROUPS") && ENV["MGROUPS"] != "0"                 # JuliaSolver.jl: the pivoted banded solver for the wider block system
const WITH_HOT = get(ENV, "HOT", "0") == "1"                                         # HOT=1: stellar EOS, nuclear network, convective mixing in the hot core

"The module files, in the order of popiii_m2c.jl."
const MODULE_FILES = String[
    "TitanDriver", "JuliaPhysics", "TitanResidFast", "BlockSolver", "TitanCore", "Z0Physics", "TitanMaterials",
    "ProblemSetup", "PrimordialChem", "SCVH", "PopIIIEqEOS", "ContinuumOpacity", "G2Opacity", "ChemNet", "ChemRates",
    "PopIIIMaterial", "H2Lines", "StellarOpacity", "IonRay",
    (WITH_MGROUPS ? ["JuliaSolver"] : String[])...,
    (WITH_HOT ? ["StellarEOS", "NucNet", "StellarMaterial", "StarMaterial"] : String[])...,
    "Trbdf2", "NewtonKrylov", "StageCons",
    "H2Levels"]      # (popiii_m2c.jl includes it into Main after the driver's line functions; here it is part of the image, Main.H2Levels is bound to it)

# --- the guard: keys and the values they have now --------------------------------------------------------------------------------
"Reviewed load-time uses of the thread count that an `__init__` of the same file recomputes (file, symbol, ENV key that feeds it): TitanResidFast.MAT_THREADS (MATTHREADS, default the thread count), H2Levels.POPWS_POOL and ChemNet.LNKMEMO (one cache per thread id). `scan` checks the `__init__`s exist."
const REREAD = Tuple{String,String,Union{String,Nothing}}[("TitanResidFast.jl", "nthreads", "MATTHREADS"), ("H2Levels.jl", "maxthreadid", nothing), ("ChemNet.jl", "maxthreadid", nothing)]   # (symbols as strings: this file is scanned too)
"""Keys the scan reaches only through its by-name approximation, which are read when a value is first needed, not at load: EQEOS_CACHE (the
path of the equilibrium-table cache, `PopIIIEqEOS.cache_path()`, called from the lazy `table()`). Each has a check that is run at build
time, after the includes, and fails the build if the key could have been used at load time."""
const RUNTIME_ONLY = Dict{String,Function}("EQEOS_CACHE" => () -> PopIIIEqEOS.TAB[] === nothing)     # the table was not built while loading
const SCAN = scan(vcat([joinpath(SRC, f * ".jl") for f in MODULE_FILES], [@__FILE__]); dynamic_include_ok = Set(["TitanStack.jl"]), skip = Set(["envguard.jl"]), reread = REREAD)
isempty(SCAN.problems) || error("TitanStack: the environment scan failed, the image would not be safe:\n" * join(SCAN.problems, "\n"))
"Every key (or pseudo-key) that the image's load-time state depends on, with the value it was built with."
const BAKED = Pair{String,Union{String,Nothing}}[k => current(k) for k in SCAN.keys if !haskey(RUNTIME_ONLY, k)]

# --- the modules ---------------------------------------------------------------------------------------------------------------
for f in MODULE_FILES
    include(joinpath(SRC, f * ".jl"))
end
for (k, ok) in RUNTIME_ONLY
    ok() || error("TitanStack: $k is read at load time after all (the build-time check of RUNTIME_ONLY failed): remove it from RUNTIME_ONLY")
end
"Names of the modules this image defines (popiii_m2c.jl binds them in Main)."
const MODULES = Symbol[Symbol(f) for f in MODULE_FILES]

# --- precompile statements from a --trace-compile run of the production configuration (generate with julia/TitanStack/trace_to_pc.jl) ---
const PC_APPLIED = Ref((0, 0))
let ok = 0, bad = 0, f = joinpath(@__DIR__, "pc.jl")
    if isfile(f)
        include_dependency(f)
        for l in eachline(f)
            try
                (Core.eval(@__MODULE__, Meta.parse(l)) === true) ? (ok += 1) : (bad += 1)
            catch
                bad += 1
            end
        end
    end
    PC_APPLIED[] = (ok, bad)
end

"The keys that differ now (empty: the image can be used)."
check() = mismatches(BAKED)

function __init__()
    # the directories TitanFortran creates when it is loaded (done at build time otherwise)
    for d in ("jl_runs", "p2_runs", "p3_runs", "p4_runs", "p4_legs"); mkpath(joinpath(TitanDriver.TitanFortran.ROOT, "lib", d)); end
    bad = check()
    isempty(bad) || error(mismatch_message(bad, length(BAKED)))
end

end # module
