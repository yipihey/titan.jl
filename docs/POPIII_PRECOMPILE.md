# Precompiled start (TitanStack)

`julia/TitanStack/` is a package whose module includes the Pop III module files from where they are (`julia/*.jl`, in the order of
`popiii_m2c.jl`) so that Julia caches their compiled code. `popiii_m2c.jl` uses it when `TITAN_PKG` is set; unset (default), the
driver includes the files into Main exactly as before. The driver itself (`julia/PopIIIRun.jl`, below) is one of the module files. Nothing about physics, numerics or defaults differs between the two paths.

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

## The driver is in the image (PopIIIRun)

`popiii_m2c.jl` is a thin loader: it includes the module files (or loads the image), then `PopIIIRun.main(ARGS)`.

- `julia/PopIIIRun.jl`: module `PopIIIRun`, the last module file of TitanStack (`GridEq.jl`, which the GRIDEQ driver code needs, is now a module
  file of the stack too). It holds every function and struct definition of the former script and of `grideq_driver.jl`, the driver's constant
  Refs and counters, the driver's few load-time constants, and a `global` declaration with an initial value for every driver variable that its
  functions read. `run_julia!`'s step callback, formerly an anonymous function of the script, is the named function `output_step` (it captured
  nothing but globals); `pchip`, `interp`, `finy`, `vfin`, `Tfin`, `xefin`, `xinit`, `ytup`, `Bν_`, `J_cmb`, `J_zero`, `coolnet` (helpers and
  closures the script defined inside its IC blocks and that the material and `iceulbuild!` code is specialised on) are named module functions
  with the same bodies; with `HOT=1` the stellar-core helpers (`e_blend`, `T_blend`, `hotburn!`, `SM`, ...) are module level too.
- `julia/popiii_run_body.jl`: the remaining statements, in the original order, evaluated by `main` in module `PopIIIRun` with `Base.include`
  when the run starts. Everything that reads ENV or ARGS, builds the IC and sets solver switches is here, so none of it is ever baked into an
  image. `julia/grideq_driver.jl` is now only the GRIDEQ settings (it is included by the body at the old place).
- The script's top-level `const X = get(ENV, ...)` (`TIMING`, `PRINTEVERY`, `CKEVERY`, `ENVEUL`, `XADV2`, `GRIDEQP`, `GEQW`, ... 38 + 8 names) are
  now assigned by the body to globals of the module. A global that functions read is declared with a type where simple (`TIMING::Bool`) and an
  initial value: reads stay type-stable, and the first assignment at run time is a plain store (a bare `global x` declaration makes it a
  binding change that invalidates every cached method reading `x`; a binding created at run time does the same, which is why the helpers
  above were moved). Constants without ENV (`Ref`s, counters, `Msun`) are `const` in the module. No load-time ENV constant was kept except
  `HOTON` (`HOT`), a key the guard has anyway (the HOT includes depend on it).
- The guard scans `PopIIIRun.jl`'s load-time code but not its function and lambda bodies (`TitanStack.RUNTIME_FILES`, `scan(...; runtime)`): the
  bodies read ENV with non-literal keys (the checkpoint environment, the RUNREC meta file) and run only from the body, never while a module loads.
- Main globals of the old script are now `PopIIIRun.<name>` (`HELDDUMP=1` lists them under the old `Main.` prefix). Nothing in `julia/` or
  `scripts/` includes `popiii_m2c.jl` or reads its Main globals (checked by grep); the shell launchers `scripts/popiii_run/iceulrun.sh` and
  `scripts/g2audit/{bitrun,exprun,exprun2}.sh` call it as a command, which is unchanged.
- The statements the body still evaluates and compiles at run time: the IC builder's helpers defined inside `begin`/`if` blocks (`readfig2`,
  `fincurve`, `rhs_hs`, ...), the closures stored in Refs (`TX.CPL_ZONE[]`, `TD.STEPFN[]` wrappers), the HOT block's lambdas, and anything a
  `RESTART` or a diagnostic switch defines. What is left of the compile time with `TITAN_PKG=1` is listed in the report of this step.

## Regenerating the precompile statements

A run on the image with `--trace-compile` lists what still compiles; `trace_to_pc.jl` turns it into statements, which are added to `src/pc.jl`
(union, order kept) and applied at the next build. Three rounds were enough: the trace of the first run on a new image shows the driver's own
specialisations (they are only named like that once the driver is in the module), the next one what the new statements did not cover.
`julia --trace-compile=trace.jl --trace-compile-timing popiii_m2c.jl lines 128` (TITAN_PKG=strict, production settings), then
`julia julia/TitanStack/trace_to_pc.jl trace.jl >> new.jl; awk '!seen[$0]++' julia/TitanStack/src/pc.jl new.jl > pc.tmp; mv pc.tmp julia/TitanStack/src/pc.jl`.
A trace of the include path names the driver's functions `Main.PopIIIRun...`, which `trace_to_pc.jl` drops: use image traces.

### Measured (23 Nov 2026, Julia 1.12.6, 8 threads, NUMA node 3, grey stack defaults from the IC, `lines 128`)

| | include path (= cc351b8) | image, driver in Main (cc351b8) | image, driver in the image |
|---|---|---|---|
| 3 steps, whole run | 382.5 s | 92.2 s | 67.7 s |
| time to the first step (whole run minus steps 2-3) | 349 s | 73 s | 50 s |
| 100 steps | 664.1 s | 363.5 s | 338.7 s |
| MGROUPS=4, first step | 635.1 s | 260.1 s | 231.5 s |
| image build (grey / MGROUPS=4) | - | 249 s | 283 s / 286 s |

The include path of this commit is byte-identical in output to cc351b8 (step lines, all non-timing stdout, decompressed tstate, 3 and 100
steps; MG 3 steps; RESTART of the 3-step checkpoint for 2 steps also identical to the image run, including the `.aux`). The image run is
identical to both. The `.aux` files differ between any two runs, including two runs of cc351b8 itself: `TD.FXEDD` has an uninitialised tail
that is written to the checkpoint (not fixed here).

Compile still traced with `TITAN_PKG=1` (sum of `--trace-compile-timing`, inclusive and overlapping: 77 s before, 41 s now): about 22 s are
statements flagged `# recompile` (cached code that Julia invalidated when ForwardDiff / StaticArrays methods were added at load; the same
statements as before this step), 12 s are specialisations on closures the body still creates at run time (`PopIIIRun.var"#NNN#MMM"`: the
HOT block's lambdas, the `TX.CPL_ZONE[]`-type closures, generators in the IC), 7 s are `Printf.format` specialisations for the format
strings of `pr`. Nothing in Main compiles any more.
