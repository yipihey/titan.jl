# titan.jl

A verified, incremental port of **TITAN**, the 1D implicit adaptive-grid radiation-hydrodynamics code of
Gehmeyr & Mihalas (1994, *Physica D* 77, 320), from Fortran 77 to Julia. The scientific target is dynamical
models of quasi-stars under infall: accretion shocks with radiative precursors on radiation-dominated
envelopes, with a central black-hole luminosity. See `docs/PORTING_NOTES.md` for the full record.

## Status (Phase 4)

| Piece | State | Verification |
|---|---|---|
| TITAN 1.2 Fortran, gfortran 13 | builds; 3 one-line build fixes (see git history) | all 9 benchmark legs reproduce the 1995 step counts where the run ends on time |
| `libtitan.so` + Julia COMMON-block wrapper | done | Julia-driven runs byte-identical to the standalone executable |
| Time-step control flow (`step.f`) | Julia (`TitanDriver.step!`) | byte-identical on all 9 legs |
| EOS / opacity (ideal gas, benchmark opacities) | Julia | tier (a) ≤ 3e-16 |
| **All residuals + Jacobian** (7 equations, ForwardDiff, stencil colouring) | Julia (`TitanResidFast`) | Jacobian vs analytic ≤ 2e-11 of row max |
| Linear solve | Julia static block-pentadiagonal LU with row equilibration (`BlockSolver.Penta`) | vs Fortran matslv 7e-11 |
| Full Julia path (`jfast`) vs Fortran (both with bookkeeping fixes) | SOD / SUBC / SUPC: identical step counts, observables ≤ 1e-12 | tier (c) ACCEPT |
| Speed (steady state, excluding JIT) | Julia 4.1–5.6 ms/step vs Fortran 6.0–6.7 ms/step | 2-core VM, ±20 % |
| Still Fortran | `peq` (Eddington factors), `maxdel`, `update`, `timstp`, `reset`, `totnrg`, I/O | — |

**Bugs found in TITAN 1.2**, fixed in the Julia driver. Both fixes are on by default; set `TITAN_BITWISE=1` to reproduce 1995:
1. **Energy diagnostic:** `totnrg` runs after `timstp`, so each step's boundary work is booked with the *next* dt,
   which is 20× too large before a divergence retry. This caused the reported 2.6–4 % "energy errors".
2. **Clock:** after a Newton-divergence retry `reset` rewinds the time, so the retried dt is never added and all later times are early (0.2–0.5 %).
3. **`lopac = 5` opacity derivatives** are wrong below ~3×10⁵ K. No benchmark uses this option.

**Genuine non-conservation** of the internal-energy scheme: 1.5–3×10⁻³, converging roughly as first order in the zone count and concentrated at
shock formation and wall reflection. A **total-energy formulation** (`TitanResidFast.ETOTAL[] = true`) conserves to 1e-13 and
removes SOD wall heating. It is still a prototype: SUPC takes 2.4× more steps.

## Layout
```
fortran/            TITAN 1.2 sources (pristine import + build-fix commit), input decks titan.in01..10
lib/                Makefile → libtitan.so, titan.x, titan_trace.x;  src/: entry points, eos/opac hooks, NCELL override
julia/              TitanFortran (ccall + COMMON map), TitanDriver (step control), TitanResidFast (residuals + AD Jacobian),
                    BlockSolver, DampedNewton (experimental), JuliaPhysics, Metrics, tierc.jl (acceptance), phase3_run.jl (runner)
scripts/            benchmark runners and sweeps (paths relative to the repo)
reference_runs/     final-time profiles of the 9 legs and the observable reference tables (metrics_reference*.csv)
docs/PORTING_NOTES.md   architecture map, build notes, verification policy, Phase 1–4 results
```

## Quick start
```bash
cd lib && make -j1                       # needs gfortran, objcopy
cd ../julia && julia --project=env -e 'using Pkg; Pkg.instantiate()'
julia --project=env phase3_run.jl SUPC 9 std     # Julia driver, Fortran routines
julia --project=env phase3_run.jl SUPC 9 jfast   # all residuals + solver in Julia
```
Problem numbers: SOD 3, BLST 4, SDOV 5, HEAT 6, COOL 7, SUBC 8, SUPC 9, RLST 10 (see `scripts/legs.sh`). Run outputs go to the git-ignored `lib/p*_runs/`.

## Verification policy (summary)
The port is **not** required to be bitwise identical. Parity mode, which calls glibc libm through `ccall`, is kept only to debug translations.
1. **Per routine:** relative difference ≤ 1e-12 on identical inputs.
2. **One step from an identical state:** within the Newton tolerance.
3. **Physical observables at fixed times:** shock and contact positions 0.5 % of the length scale, temperatures 1 %, step count 3 %.

## Licence
The TITAN sources in `fortran/` are © 1995 University of Illinois. They are redistributed under the permission notice in `fortran/README`,
which must be retained. The Julia code and build scaffolding are new work.
