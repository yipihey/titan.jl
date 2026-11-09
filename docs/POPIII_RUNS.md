# Running the Pop III protostar from scratch

How to run `popiii_m2c.jl` from the cloud to the accreting protostar with the current physics, record it, and see where the
time goes. First run with this procedure: `full1` (started 2026-11-09).

## Procedure

`scripts/popiii_run/fullrun.sh <run> [threads]` (needs `JLENV`, a Julia project with the repository's packages; data go to
`$DATA/runs/<run>`, default `DATA=/zpool/nvme/data/tabel_scratch/titan_popiii`). It writes `procedure.log` with every phase
start, its full settings and how it ended.

| Phase | Mesh | Settings | Ends |
|---|---|---|---|
| A: collapse | Lagrangian, rezoned every 20 steps (`REZONE=1`) | `phaseA.env`: the production physics without the accretion-only Δt rules (`CDYN=1000`, no `DTENV`) | at n_c = 2.55e23 cm⁻³ (log 23.41, where the earlier history switched), final checkpoint `z19<run>A` |
| B: accretion | remapped onto the fixed accretion mesh at the restart (`ENVEUL=0,0.03,0.1,1000`; `REZONE=0`) | `phaseB.env`: the production flags | segments of 100 000 steps, each from the previous segment's final checkpoint `z19<run>B<n>`; `ENVMS0` is read from the remap line in the log |

A segment that ends any other way than `:maxsteps` stops the chain; `procedure.log` says how.

Notes:
- `NLOPAC=1` from the cloud needed a guard (5337e9f): before RAD's opacity tables exist the non-LTE factors are skipped.
- Phase A uses the collapse solver controls `NITER=10 STOL=0.1 NTRY=10` (those of the old collapse). With the accretion
  controls (`NITER=25 STOL=0.3 NTRY=25 FAILTFAC=0.5`) 61 % of TR stages diverge from step 1 (15 % with the collapse
  controls), and `full1` died at step 3375 (log n_c 5.43): BDF2 failed on every step from 3363, Δt fell from 3e10 to 350 s,
  then NaN everywhere and a NaN Jacobian in the chemistry. The checkpoint carries the controls, so a restart must set them.
- Open: the divergences that remain have Newton corrections of 1e17–1e21 in the central zones 1–3 (dm 3e-7 M☉).
- `MSG=1` now prints Newton messages from the start (they were silent until RAD switched on).
- `fullrun.sh <run> [checkpoint | B:<checkpoint>] [threads]` resumes phase A or B; every attempt saves under a fresh name.
- `m*`/`R*` in the `hot:` line and perf.csv are meaningful only after the remap (they take the first interface with
  u < −2 km/s, which in the cloud is the infall front).
- Checkpoints of a run live in `runs/<run>/ck` (`CKDIR`), apart from the shared `ck/`.

## Recording (`RUNREC=<dir>`)

Per process segment (start step, pid): `snaps_<k0>_<pid>.jls`, one snapshot every `SNAPEVERY` steps, appended and never
rewritten, and `meta_<k0>_<pid>.txt` with the settings the code reads, the git commit and whether `julia/` was modified, host
and threads. A snapshot holds r, u, m, ρ, T, e, p, E_r, F, χ, κ_P, Eddington factors, equilibrium flags, line cooling, all
network species, the lagged chemistry heat, the stellar composition of the hot zones, L_nuc, L_lines, chemistry fallbacks and
the star-light T_eff. About 85 KB per snapshot at 200 zones; `SNAPEVERY=20`. `SNAPS=0` (the old in-memory history is
rewritten in full at every snapshot and stored in each checkpoint).

`perf.csv`, one row per `PRINTEVERY` steps (needs `TIMING=1`): wall seconds in the hydro+radiation solve, chemistry, H₂
lines and the rest; GC seconds and GB allocated; cumulative Newton divergences and repeated steps, TR-BDF2 steps and stage
failures, BE fallbacks, ALI calls, chemistry fallbacks; threads, zones, log n_c, T_c, m*, R*, fastest infall and the variable
and zone that limit Δt.

## Looking at a run

- `julia --project=<env> julia/runrec.jl runs/<run> [every]` assembles the surviving branch over all segments (a restart
  from an earlier checkpoint re-runs steps; the abandoned ones are dropped by a backward scan in time) and writes
  `export/frames.f32`, `frames.json` and `track.csv`. `runrec.jl --old <snaps file> <outdir> [every]` converts an old
  `SNAPS=1` history the same way (four species only, no κ_P).
- `python scripts/popiii_movie.py runs/<run> [--x r|m] [--fps 24]`: nine profile panels and the track, fixed axes, mp4.
- The interactive explorer (artifact "Pop III Run Explorer") reads the same export, base64-encoded, with the run's perf.csv;
  it overlays the reference history matched by protostar mass (by central density before the protostar).
