# Precompiled start (TitanStack)

`julia/TitanStack/` is a package whose module includes the Pop III module files from where they are (`julia/*.jl`, in the order of
`popiii_m2c.jl`) so that Julia caches their compiled code. `popiii_m2c.jl` uses it when `TITAN_PKG` is set; unset (default), the
driver includes the files into Main exactly as before. Nothing about physics, numerics or defaults differs between the two paths.

## Use

```bash
# once per environment (never into the shared production environment: copy it first)
julia --project=$JLENV -e 'using Pkg; Pkg.develop(path = "julia/TitanStack")'
# (re)build the image with the environment of the runs it will serve (also happens automatically on the first `using` after a code change)
env <run settings> julia --project=$JLENV julia/TitanStack/build.jl
# a run
env TITAN_PKG=1 ... julia -t 8 --project=$JLENV popiii_m2c.jl lines 128
```

`TITAN_PKG=1`: use the image; if it was built with a different load-time environment, say why and run on the include path.
`TITAN_PKG=strict`: stop instead of falling back.

## The environment guard

An image freezes everything the modules compute when they are loaded: `const X = get(ENV, ...)`, conditional includes (HOT, MGROUPS,
LINEGPU), defaults read once. `src/envguard.jl` finds these mechanically from the sources at build time:

1. every file the package includes is parsed (the module files and the files they include in turn, resolved from the AST);
2. load-time code = everything outside function / macro / lambda bodies; each `ENV` use in it (`get`, `haskey`, `ENV[..]`, ...) gives a key;
   an `ENV` use with a non-literal key, or the whole `ENV`, fails the build;
3. every identifier in load-time code that names a function of the scanned files adds the keys read in that function's body, transitively
   (matching by name only: over-strict);
4. `homedir()` adds `HOME`. Other value-fixing calls in load-time code (`Threads.nthreads()`, `maxthreadid()`, `threadid()`, `time`, `getpid`,
   `ccall`, ...) fail the build, except at the reviewed sites in `TitanStack.REREAD`: the three per-thread tables sized at load
   (`TitanResidFast.MAT_THREADS`, `H2Levels.POPWS_POOL`, `ChemNet.LNKMEMO`; the build process has one thread) are recomputed in a
   one-line `__init__` of the same file, which the scan checks exists. These `__init__`s also run on the include path, where they
   recompute the value the module body just computed.

The value (or absence) of every key is stored in `TitanStack.BAKED`; `__init__` compares them with the environment of the loading
process and throws, naming the key and both values. The conditional includes are the conditions of `popiii_m2c.jl`; `LINEGPU` set
at build time fails the build (the GPU formal solution is not part of the image).

`EQEOS_CACHE` (the equilibrium-table cache path differs per run directory) is read when the table is first needed
(`PopIIIEqEOS.cache_path()`, was `const CACHE`), so that runs with different cache files share one image. The by-name closure of
step 3 would still list it; `TitanStack.RUNTIME_ONLY` removes it and the build checks that the table was not built while loading.

`H2Levels` is part of the image (`import ..H2Lines`, was `import Main.H2Lines`); the driver binds `Main.H2Levels` and `Main.H2Lines`
to the image's modules, which is where `H2Lines` / `PopIIIMaterial` look them up.

One image per configuration: the newest build wins when Julia loads the package, so a second configuration needs its own depot
(`JULIA_DEPOT_PATH=<dir>:~/.julia`) or a rebuild with `build.jl`.

## Regenerating the precompile statements

`julia --trace-compile=trace.jl --trace-compile-timing popiii_m2c.jl lines 128` (include path, production settings), then
`julia julia/TitanStack/trace_to_pc.jl trace.jl > julia/TitanStack/src/pc.jl`.
