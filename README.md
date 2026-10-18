# titan.jl

A verified, incremental port of **TITAN**, the 1D implicit adaptive-grid radiation-hydrodynamics code of
Gehmeyr & Mihalas (1994, *Physica D* 77, 320), from Fortran 77 to Julia. The scientific target is dynamical
models of quasi-stars under infall: accretion shocks with radiative precursors on radiation-dominated
envelopes, with a central black-hole luminosity. See `docs/PORTING_NOTES.md` for the full record.

**The book:** https://yipihey.github.io/titan.jl/ documents every experiment: the steady quasi-star model (`steady/`,
with an interactive WebAssembly port), the port of TITAN, and the quasi-star runs. Source in `report/`.

**Phase 7 (quasi-star physics):** Z = 0 Saha material, MLT convective flux, quasi-star set-up (`julia/QuasiStarSetup.jl`,
`julia/qstar_run.jl`). The hydrostatic test holds the interior to 0.5 % and validates the MLT flux. With four numerical
defects fixed (artificial-viscosity length, gas energy by cancellation, Eddington-factor source for tabulated opacities,
TITAN's gas + radiation energy row; P7.5) the infall runs are stable from 3e-3 to 3 M⊙/yr, to 2.6e8 s (run R2).
P7.6 adds an inner mass sink (and fixes a sign error in TITAN 1.2's inner mass flux), second-order TR-BDF2 time stepping
(`julia/Trbdf2.jl`), and an adaptive grid that works on this problem; the adaptive runs resolve a hot cooling layer
behind the shock that the Eulerian grid misses, which changes the late evolution. P7.7: MLT in non-hydrostatic layers,
a positivity-preserving Eddington-factor solution, `TITAN_MGR` (array size), zone-count convergence at 220/314/440 zones,
and runs to 6e8 s that settle within 1.4e-3 R and 50 K of the steady model. P7.8: emergent spectra of TITAN snapshots
(`julia/Spectrum.jl`, `spectra_run.jl`, `shock_budget.jl`, `line_profile.jl`): static LTE continuum with Unsöld–Lucy
correction, the shock's X-rays through the infall column, the Balmer jump, and bracketed Hα/Hβ profiles of the infall.
P7.9: non-LTE hydrogen in the infall (`julia/NLTEH.jl`, `nlte_run.jl`): the star keeps the inner infall ionized and
puts Hα into strong emission.
P7.10: a wide-domain run with the outer boundary at 30 R (W30; `QS_XSCALE`, `MLT_NOINFALL`), whose photosphere matches
A220's and whose outer infall is close to dilute radiative equilibrium; the full non-LTE hydrogen line spectrum
(`nlte_spectrum.jl`, 23 lines, Balmer decrement 5.1) with a per-line energy check (`nlte_linepower.jl`). The non-LTE
calculation on W30 is still to do.

## Status (Phase 6)

**SOD, SUBC and SUPC now run in pure Julia, with no Fortran at run time.** `libtitan.so` is never loaded (`TITAN_NOLIB=1`).
The initial states were produced once by TITAN's own set-up routines and are stored in `reference_runs/init/*.tstate.gz`.

| Piece | State | Verification |
|---|---|---|
| TITAN 1.2 Fortran, gfortran 13 | builds; 3 one-line build fixes (see git history) | all 9 benchmark legs reproduce the 1995 step counts where the run ends on time |
| `libtitan.so` + Julia COMMON wrapper | kept as reference and for problem set-up | Julia-driven runs byte-identical to the standalone executable |
| Every routine of the time step | Julia: driver, residuals + ForwardDiff Jacobian, block solver, EOS/opacity, peq, eddfac, update, maxdel, timstp, totnrg, reset | each checked call by call against Fortran over whole runs (`PortHarness`): identical or ≤ 1e-12, apart from 1-ulp `pow` cases |
| State | Julia-owned arrays (`TitanFortran` `:julia` backend), same names and shapes as the COMMON blocks | — |
| Pure-Julia runs vs. fixed Fortran reference | SUBC 587/587 and SUPC 713/713 steps; SOD 455 vs 452 | tier (c) ACCEPT on all three |
| Speed (steady state, excluding JIT) | 2.4 (SOD), 4.3–5.0 (SUPC/SUBC) ms/step vs Fortran 6.4 | 2-core VM |
| Still Fortran | problem set-up (`start*.f`, `gridinit.f`) and the original output files | to be replaced by Julia set-ups (quasi-star first) |

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
TITAN_NOLIB=1 julia --project=env phase6_run.jl pure SUPC   # pure Julia from the stored initial state
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
