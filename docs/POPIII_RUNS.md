# Running the Pop III protostar from scratch

How to run `popiii_m2c.jl` from the cloud to the accreting protostar with the current physics, record it, and see where the
time goes. Runs with this procedure: `full1` (2026-11-09, earlier phase A settings, stopped at log n_c 16.9) and `full2`
(2026-11-09, the phase A settings below).

## Procedure

`scripts/popiii_run/fullrun.sh <run> [threads]` (needs `JLENV`, a Julia project with the repository's packages; data go to
`$DATA/runs/<run>`, default `DATA=/zpool/nvme/data/tabel_scratch/titan_popiii`). It writes `procedure.log` with every phase
start, its full settings and how it ended.

| Phase | Mesh | Settings | Ends |
|---|---|---|---|
| A: collapse | Lagrangian, rezoned every 20 steps (`REZONE=1`) | `phaseA.env`: the production physics without the accretion-only Δt rule `DTENV`, with `CDYN=0.01 ADVREL=1 EPGRAD=thermal TCOOLTH=1` (see "Stability of the collapse") | at n_c = 2.55e23 cm⁻³ (log 23.41, where the earlier history switched), final checkpoint `z19<run>A` |
| B: accretion | remapped onto the fixed accretion mesh at the restart (`ENVEUL=0,0.03,0.1,1000`; `REZONE=0`) | `phaseB.env`: the production flags | segments of 100 000 steps, each from the previous segment's final checkpoint `z19<run>B<n>`; `ENVMS0` is read from the remap line in the log |

A segment that ends any other way than `:maxsteps` stops the chain; `procedure.log` says how.

Notes:
- `NLOPAC=1` from the cloud needed a guard (5337e9f): before RAD's opacity tables exist the non-LTE factors are skipped.
- Phase A uses the collapse solver controls `NITER=10 STOL=0.1 NTRY=10` (those of the old collapse). With the accretion
  controls (`NITER=25 STOL=0.3 NTRY=25 FAILTFAC=0.5`) 61 % of TR stages diverge from step 1 (15 % with the collapse
  controls), and `full1` died at step 3375 (log n_c 5.43): BDF2 failed on every step from 3363, Δt fell from 3e10 to 350 s,
  then NaN everywhere and a NaN Jacobian in the chemistry. The checkpoint carries the controls, so a restart must set them.
- Phase A runs without `ETOTAL=1` (and `ETASW`, `EXCAP`). With it, `full1` died twice at the same physical time
  (t = 3.309e14 s, log n_c 5.43, steps 3375 and 5196 with different solver controls): in the cold supersonic outer infall
  (n ~ 1e4, u ≈ −4 km/s, kinetic ≫ thermal energy) two zones cooled in the hydro step from ~35 K to 17 K, below T_CMB,
  until Newton failed. Without it the same checkpoint passes. Phase B keeps it.
- `MSG=1` now prints Newton messages from the start (they were silent until RAD switched on).
- `fullrun.sh <run> [checkpoint | B:<checkpoint>] [threads]` resumes phase A or B; every attempt saves under a fresh name.
- `m*`/`R*` in the `hot:` line and perf.csv are meaningful only after the remap (they take the first interface with
  u < −2 km/s, which in the cloud is the infall front).
- Checkpoints of a run live in `runs/<run>/ck` (`CKDIR`), apart from the shared `ck/`.

## Stability of the collapse (phase A)

Three defects showed in `full1` between log n_c 5 and 17. Each has its own option; all are off by default, so phase B and
the accretion production runs are unchanged.

1. **Newton divergences (one per ~6 steps), `ADVREL=1`.** With frozen interior masses the interfaces move with the gas
   only at convergence. In an unconverged iterate the continuity row carries mass r² u_rel ρ across the interfaces, but
   that mass carried no energy and no momentum. The energy row then holds e × (mass defect), and with e = 0 for ionised
   gas neutral gas has e ≈ −1e13 erg/g ≈ 400 c_v T at 250 K: a zone gaining 1e-3 of its mass in the linearised step is
   heated by 40 %. Above a critical Δt the Jacobian is singular at the predicted state (one singular value 1e-15, a
   zone-to-zone velocity mode of the central zones) and Newton diverges from its first iteration. With `ADVREL=1` the
   residual mass flux advects energy and momentum; the converged solution is the same. Divergences: 0 in all tests.
2. **Δt, `CDYN=0.01`.** The divergences had been the de-facto Δt control. Without them `STOL=0.1` alone lets the central
   density change ~10 % per step and the collapse comes late. Offset of the collapse time against `full1` at equal
   log n_c: `CDYN=1000` +46 000 yr, `0.03` +7 000 yr (T_c within 1 %), `0.01` +1 950 yr (T_c within 0.1 %), `0.003`
   +165 yr (to log n_c 10). `CDYN=0.01` (Δt ≤ 0.01 × the shortest √(r³/Gm)) takes ~125 steps per dex of central density.
   The cloud phase needs it too. `full1` ran its first 3000 steps with `CDYN=1000` (and `ETOTAL=1` and the accretion
   solver controls, with 61 % of the TR stages diverging) and collapsed about twice too fast: log n_c 16.9 at 10.77 Myr,
   where the October collapse that the production history starts from (`C_dyn = 0.01`) reached log n_c 23.4 at 21.86 Myr.
   At 8.9 Myr `full1` had log n_c 4.21 and T_c 186 K; the October run 3.73 and 151 K; `full2` 3.71 and 154 K (step 1000
   in both). Which of the three settings sped `full1` up was not separated. The offsets above are measured from a common
   checkpoint of `full1` at log n_c 5.06, so they stand, but `full1` itself is not a valid collapse history.
3. **Ringing of the core while it turns molecular (log n_c 10.5–13), `EPGRAD=thermal`.** Once the H₂ lines are thick the
   cooling depends on the escape probability, whose Sobolev velocity gradient was taken from zone-local quantities
   (|Δu|/Δr across the zone and |ū|/r_c). In the subsonic core (|u| ~ 0.01 km/s, v_th ≈ 3 km/s) both are noise: a 10 m/s
   wave (Mach 0.003) at zone 1 changed the escape-probability cooling by 40 % from one step to the next while the
   cooling from the full line transfer moved by 6 %. A zone that contracts faster then cools more, so the waves grow
   (std of ln ρ over the 20 central zones 0.1–0.4, velocities to 1.7 km/s); the periodic full line-transfer call adds
   10–18 % jumps on top. With `EPGRAD=thermal` both gradients are taken over the distance in which the line-of-sight
   velocity changes by one thermal width, which is what a line photon samples. Std of ln ρ then stays 1e-4 at log n_c
   10.4 and 4e-3 at 16.5, with the same result at `CDYN` 0.003, 0.01 and 0.03. Explicit against implicit cooling, the
   interval of the full line-transfer calls, rezoning and Δt alone were each tested and are not the cause.
   `TCOOLTH=1` (cooling-time rule from the thermal energy instead of the energy above cold H₂) made no difference
   once the gradient was fixed; it is kept as the physically meant rule.
4. **Chemistry failure at log n_c 16.9** (`ChemNet.step: Newton failed`, where `full1` stopped): H⁺ ≈ 4e-12 in fast
   equilibrium flips in its 10th digit between Newton iterations while T and every other abundance are converged, so the
   test (1e-10 per species) was never met at any sub-step. A stagnated iterate is now accepted after 60 full Newton steps
   if it is converged to 1e-6 in the abundances and 1e-8 in ln T. Steps that converged before are unchanged bit for bit.
   This one is not an option: it applies to every run.

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
  it overlays another run matched by protostar mass (by central density before the protostar).
