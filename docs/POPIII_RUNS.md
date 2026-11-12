# Running the Pop III protostar from scratch

How to run `popiii_m2c.jl` from the cloud to the accreting protostar with the current physics, record it, and see where the
time goes. Runs of 2026-11-09: `full1` (total-energy row in the collapse, stopped at log n_c 16.9; not a valid collapse),
`full2` (single-stage phase A; valid below log n_c ≈ 14 only, see "Radiation must be switched on by a rebuild") and
`full3` (the staged procedure below).

## Procedure

`scripts/popiii_run/fullrun.sh <run> [threads]` (needs `JLENV`, a Julia project with the repository's packages; data go to
`$DATA/runs/<run>`, default `DATA=/zpool/nvme/data/tabel_scratch/titan_popiii`). It writes `procedure.log` with every phase
start, its full settings and how it ended.

| Stage | What runs | Settings | Ends |
|---|---|---|---|
| A1: cloud collapse | Lagrangian mesh, 200 zones. Hydrodynamics, chemistry, H₂ line transfer, RA04 CIE continuum cooling. No radiation field (`RAD=0`) | `phaseA.env` + `stageA1.env` (`CDYN=0.01 ADVREL=1 EPGRAD=thermal TCOOLTH=1`, see "Stability of the collapse") | n_c = 3.5e18 cm⁻³ (log 18.54) |
| A2: radiation, core formation | `RAD=1 RADFROM=<A1 checkpoint>`: E_r and F_r are built from the hydro state and the StellarOpacity continuum replaces the CIE escape-factor term. The future core is pre-refined while it is still smooth (`REFINE=1e-3,0.1,0.03`: 1e-3 to 0.1 M☉ to Δm/m ≤ 0.03) | `phaseA.env` + `stageA2.env` (`DTRADMAX=1500`, `CDYN=0.01`) | T_c = 5e4 K (the core has formed) |
| A3: core growth | Lagrangian mesh with a moving refinement window ahead of the star (`REZWIN=1.3,4,0.03 REZNMAX=340`), dynamical Δt limit relaxed (`CDYN=1000 DTRADMAX=0`): with `CDYN=0.01` the core's own density holds Δt at ~90 s | `phaseA.env` + `stageA3.env` | M* = 0.05 M☉ since 11 Nov, 0.5 M☉ before (`MSTOP`: mass inside the first interface falling faster than 2 km/s) |
| B: accretion | remapped onto the fixed accretion mesh at the restart (`ENVEUL=0,0.03,0.1,1000`; `REZONE=0`) | `phaseB.env`: the production flags | segments of 100 000 steps, each from the previous segment's final checkpoint `z19<run>B<n>`; `ENVMS0` is read from the remap line in the log |

Stages A1–A3 are the sequence that produced the production history in October (S6-A/B/C in `POPIII_NOTES.md`: `z19base` →
`RADFROM=z19s9500` with `REFINE` → accretion with `REZWIN` → `z19n1`), with the present code and options. The October
hand-over was at log n_c 23.407 (M* ≈ 0.5 M☉); the stage now ends on the mass, which is what the remap needs.

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
  until Newton failed. Without it the same checkpoint passes. Phase B keeps it. The row has no gravity work
  on Lagrangian zones (see "The total-energy row on the Lagrangian mesh").
- `MSG=1` now prints Newton messages from the start (they were silent until RAD switched on).
- `fullrun.sh <run> <stage>:<checkpoint> [threads]` resumes: `A1:`, `A2:`, `A3:` or `B:` continue (or start) that stage from
  a checkpoint, `R:` starts A2 with the rebuild from a hydro-only checkpoint. Every attempt saves under a fresh name
  (`z19<run>_<stage>_<n>`).
- `m*`/`R*` in the `hot:` line and perf.csv are meaningful only after the remap (they take the first interface with
  u < −2 km/s, which in the cloud is the infall front).
- Checkpoints of a run live in `runs/<run>/ck` (`CKDIR`), apart from the shared `ck/`.
- `CHEMFBDENSE=1` (default, new on 9 Nov): a zone at p ≥ `PSWITCH` that is still on the network (T < `PSWTMIN`) and whose
  network step fails takes the equilibrium-table step instead, in every phase; counted in `nchemfb` (perf.csv). Before, that
  fallback existed only in phase B (`CHEMFALLBACK`), and both `full3` and `full4` died in stage A3: `full3` at step 11935,
  `ChemNet.step: Newton failed` in the H₂ dissociation layer just below `PSWTMIN` (n_H 4.5e18 cm⁻³, T 4346 K, Δt 515 s,
  twelve halvings without convergence); `full4` between steps 11700 and 11800 in the same way. Reproduced from the A2
  checkpoint (`t_a3b`) and run through the failure step with the fallback (`t_a3c`). Both restarted stages
  (`z19full3_A3_2`, `z19full4_A3_1`) are past it; the fallback is taken a few times per 100 steps (38 times by step 13000 in
  `full3`, 18 by step 12400 in `full4`). `CHEMFBDENSE=0` restores the old behaviour.
- A run that ends on an exception prints the first 400 characters of the message. That used a byte slice, which threw a
  `StringIndexError` on a message with a multi-byte character at byte 400 and hid `full4`'s failure; now `first(msg, 400)`.

## Stability of the collapse (phase A)

Three defects showed in `full1` between log n_c 5 and 17. Each has its own option. `ADVREL`, `CDYN=0.01` and `TCOOLTH`
are set in `phaseA.env` only. `EPGRAD=thermal` is the code default since 9 Nov 2026 (`EPGRAD=zone` gives the old method).

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
   where the October history that production starts from reached log n_c 23.4 at 21.86 Myr.
   At 8.9 Myr `full1` had log n_c 4.21 and T_c 186 K; the October run 3.73 and 151 K; `full2` 3.71 and 154 K (step 1000
   in both). The offsets above are measured from a common checkpoint of `full1` at log n_c 5.06, so they stand, but
   `full1` itself is not a valid collapse history. The cause is `ETOTAL=1` on the Lagrangian mesh (next section).
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
   `full2`, on its own history from the cloud, passes log n_c 10–12.5 with std of ln ρ 1e-4 to 6e-4.
   In the accretion phase the choice has no effect: the H₂ lines of the infalling envelope are thin (x(H₂) ≤ 1e-3), so
   the escape probability is 1 whatever the gradient. From an `fx3` checkpoint at 99.5 M☉ (step 187 550), 400 steps with
   `zone` and with `thermal` gave the same line cooling in every zone to three decimals in the log and temperatures equal
   to 5e-7 dex. `fx3` takes the new default at its next segment (`z19fx3d`).
4. **Chemistry failure at log n_c 16.9** (`ChemNet.step: Newton failed`, where `full1` stopped): H⁺ ≈ 4e-12 in fast
   equilibrium flips in its 10th digit between Newton iterations while T and every other abundance are converged, so the
   test (1e-10 per species) was never met at any sub-step. A stagnated iterate is now accepted after 60 full Newton steps
   if it is converged to 1e-6 in the abundances and 1e-8 in ln T. Steps that converged before are unchanged bit for bit.
   This one is not an option: it applies to every run.

## Radiation must be switched on by a rebuild (why `full2` is wrong above log n_c ≈ 14)

`RAD=1` does two things in the driver: it removes the RA04 CIE continuum cooling ("the continuum is in the transfer"), and
it applies the radiation settings to a state that already has a radiation field. The field itself is only created by
`RADFROM=<hydro-only checkpoint>` (`radbuild!`). `phaseA.env` as first written (5934c9e) set `RAD=1` from the initial
condition with no rebuild, so `full1`, `full2` and the `t_e_*` tests had **no continuum cooling and no radiation at all**
(no `RAD: rebuilt` line in the log; T_rad and τ are NaN and L = 0 in every recorded zone). Nothing warned.

Effect on `full2`, T_c against the October collapse at equal central density:

| log n_c | ≤ 14 | 15 | 16 | 17 | 19 | 19.9 | 22.2 |
|---|---|---|---|---|---|---|---|
| T_c `full2` / October | 1.002–1.005 | 1.016 | 1.07 | 1.13 | 1.20 | 1.8 | 1.6 (51 165 against 31 726 K) |

Below log n_c ≈ 14 the continuum is negligible and the run is valid; above, the centre is too hot and H₂ too dissociated.
`full2` was stopped at log n_c 22.45. Results obtained below log n_c 14 or independent of the continuum stand: the
total-energy row finding, `ADVREL`, the Δt convergence, `EPGRAD=thermal`, the chemistry round-off fix.

The driver now refuses both ways of getting there: `RAD=1` from the initial condition, and `RAD=1 RESTART=<checkpoint>`
of a checkpoint without a radiation field, stop with an error that names `RADFROM`. Checked on 9 Nov: a `RAD=1` start
from the initial condition stops with that error, and a restart of a production checkpoint (`z19fx3c_auto1`) runs as before.

A second defect of the single-stage phase A: after the core forms, `CDYN=0.01` limits Δt to 0.01 × the dynamical time
of the core itself (87 s at log n_c 22.45). October relaxed the limit to `CDYN=1000` for the core-growth stage.

Timing of the cloud collapse: `t_e_old` (the three new options off) follows `full2` to 2e-4 Myr up to log n_c 6.3
(23.1323 against 23.1325 Myr), so the options do not change the collapse time. Both are 1.5 Myr later than the October
collapse, the difference arising between log n_c 4 and 6. Its cause is not known (October code and settings).

## When to switch the continuum transport on (`full4` against `full3`, 9 Nov)

Between log n_c ≈ 14 and 18.5 the October sequence cools the gas with the RA04 fit of the H₂ collision-induced emission
times its escape factor, and only then builds the radiation field. The alternative is to build the field while the cloud
is still thin and let the transfer with our own continuum opacity do the cooling throughout. `full4` is that run: a
branch of `full3` at step 9500 (log n_c 14.007), `fullrun.sh full4 R:z19full4_A1_s9500` with `A2X="DTRADMAX=0"`.

- `A2X="<settings>"` (environment of `fullrun.sh`) appends settings to the rebuild stage; they override `stageA2.env`.
- `DTRADMAX` is a cap on Δt in seconds (1500 by default with `RAD=1`), not a ratio. At log n_c 14 the step is ~1e6 s, so
  an early rebuild needs `DTRADMAX=0`. With it the run went from log n_c 14 to 18.6 with 3 Newton divergences and 4
  failed trapezoidal stages in total, all recovered.
- The rebuild in a thin cloud: τ_R(centre) = 3.7e-3, E_r/aT⁴ = 0.004 at the centre, no photosphere, L_ph = 0.
  For comparison `full3`'s rebuild at log n_c 18.55: τ_R(centre) = 528, photosphere at zone 80, L_ph = 9.74e34 erg/s
  (October: 790, zone 81, 1.113e35 erg/s, with T_c 3639 K against 3523 K now).

T_c at equal central density:

| log n_c | 14.5 | 15.0 | 15.5 | 16.0 | 16.5 | 17.0 | 17.5 | 18.0 | 18.546 |
|---|---|---|---|---|---|---|---|---|---|
| `full4` (transport from 14) [K] | 1891.5 | 2020.3 | 2133.7 | 2212.2 | 2314.6 | 2528.7 | 2795.9 | 3107.5 | 3523.8 |
| `full3` (RA04 × escape factor) [K] | 1895.3 | 2031.2 | 2157.3 | 2247.3 | 2324.2 | 2526.9 | 2793.7 | 3105.2 | 3518.9 |
| ratio | 0.998 | 0.995 | 0.989 | 0.984 | 0.996 | 1.001 | 1.001 | 1.001 | 1.001 |

The transport run is at most 1.6 % cooler (log n_c 16) and the two agree to 0.1 % from log n_c 17, where the centre is
optically thick and neither loses energy. This is the sign and size a static comparison of the two coolants gives: thin
emission of our opacity, 4σT⁴κ_P ρ at equilibrium composition, over the RA04 fit is 0.95 (log n 12, 1100 K), 1.11 (13),
1.15 (14), 1.29 (15), 1.28 (16), 1.23 (17), 1.04 (18), 1.10 (18.5); the RA04 escape factor is 1.0 to log n 15, 0.85 at 16,
4.7e-3 at 17. (Only the centre is compared here.)

Cost: the pre-refinement to 326 zones comes with the rebuild, so the early rebuild carries the finer mesh through four
more decades of density, and a full line-transfer solve takes 23 s at 326 zones against 7.6 s at 200.

## The H₂ line cooling follows the level populations (checked 9 Nov)

`H2Levels.jl` solves statistical equilibrium for 63 levels (v = 0–2, J = 0–20) with collisions with H (Wrathmall, Gusdorf
& Flower 2007), ortho- and para-H₂ (Le Bourlot et al. 1999), He and H⁺ (Gerlich 1990), no electron collisions, and the
radiative rates from the line transfer (`NLTE=hyb`); the gas cooling is the net collisional exchange at those populations
(`LINECOLL=1`). Cooling per molecule over its LTE value, optically thin, x(H₂) = 1e-3, x(H⁺) = 1e-4:

| T [K] | n = 1e2 | 1e3 | 1e4 | 1e5 | 1e6 | 1e7 | n at half LTE [cm⁻³] |
|---|---|---|---|---|---|---|---|
| 200 | 0.092 | 0.358 | 0.794 | 0.965 | 0.996 | 1.000 | 2.0e3 |
| 300 | 0.067 | 0.306 | 0.662 | 0.910 | 0.988 | 0.999 | 3.3e3 |
| 500 | 0.044 | 0.180 | 0.469 | 0.808 | 0.944 | 0.990 | 1.2e4 |
| 1000 | 0.012 | 0.060 | 0.214 | 0.582 | 0.921 | 0.991 | 6.6e4 |
| 2000 | 0.010 | 0.067 | 0.367 | 0.844 | 0.982 | 0.998 | 1.8e4 |

In the run (`full3`, step 8180, log n_c 8.01; `julia/h2crit_check.jl <run dir> 8.0`) the recorded line cooling per molecule
over the local statistical-equilibrium thin value is 1.002–1.004 below n_H = 1e4 cm⁻³, 0.999 at 1.2e5, 0.996 at 3e6 and
0.991 in the core (7.7e7 cm⁻³, 630 K); the run's value includes the line transfer, the comparison value does not.
Relative to LTE it is 0.001 at n = 2.3, 0.03 at 45, 0.33 at 1.5e3 (244 K), 0.84 at 1.2e5 (349 K), 0.96 at 3e6 (507 K),
0.989 in the core. Newer H–H₂ rates (Lique 2015) are in the Cloudy data set next to the ones used and have not been
evaluated.

The recorded `Λext` is per volume (erg cm⁻³ s⁻¹); the export label said erg/s/g until 9 Nov.

## The total-energy row on the Lagrangian mesh (why `full1`'s envelope differs)

`full1` and `full2` differ from the first Myr and out to the edge of the cloud (in `full1` the shells at 2–4e3 M☉ are
2–3× colder at 6 Myr and fall at Mach 1.1 instead of below 0.5). Single-variable tests from the initial condition, 700 steps
each (6.7 Myr), with the energy budget of Lagrangian shells (`de + ∫p dV`, what the gas lost other than by expansion):

| test | setting | result |
|---|---|---|
| `t_e_old` | `phaseA.env` with `ADVREL=0 EPGRAD=zone TCOOLTH=0` | shell budget equal to `full2` in every printed digit: the three new options do not act on the cloud |
| `t_e_c3` | `CDYN=0.003` | equal to `full2` (shell budgets within 3 %, log n_c 3.577 and T_c 155.6 K at 6.94 Myr): converged in Δt |
| `t_e_dt` | `CDYN=1000 STOL=0.3`, accretion solver controls | 0.25 dex per step, unresolved, and slower than `full2`, not faster |
| `t_e_etot` | `ETOTAL=1 ETASW=0.6 EXCAP=1e8` | reproduces `full1`: the shells lose 1.3–6.4e10 erg/g (thermal energy 5–10e10) |
| `t_e_eg` | the same with `ETOTGRAVGRID=1` | equal to `full2` (budgets within 5 % or 0.01e10 erg/g, except the edge shell) |

The gravity work of the total-energy row (`ETOT_GRAV`) is written for the fixed accretion mesh: it acts on the mass that
crosses a face, r² u_rel ρ. On Lagrangian zones u_rel = 0, the row has no gravity work at all, and the kinetic energy
that gravity gives the gas is taken from its thermal energy. In the envelope the dark matter dominates gravity
(M_DM/M_gas = 15–23 outside 2500 M☉ with the read-off halo profile), so the loss is large, and it equals the work of
gravity ∫ G(m + M_DM)/r² |dr|: at 6300 and 7900 M☉ the loss beyond `full2`'s is 3.7 and 2.4e10 erg/g, the work of gravity
3.8 and 2.5e10. The cooled envelope then falls in under the external
pressure. The cooling below T_CMB in the supersonic infall that killed `full1` at log n_c 5.43 is probably the same
loss (not tested separately).

`ETOTGRAVGRID=1` (default off) adds the work on mass carried by moving faces, u_grid g (m_{j-1} + m_j)/2 per face, the
work of the momentum row's force; it is zero on a fixed mesh. With it the total-energy and the internal-energy rows give
the same cloud. The edge shell still differs (+0.33e10 against −0.06e10 erg/g): probably because the boundary flux of the row
uses the pressure of the last zone where the free-surface momentum row uses `PEXT`. Not fixed; phase A does not use `ETOTAL`.

Independent check of `full2`: its shell losses equal the optically thin H₂ cooling (Galli & Palla 1998, from the
recorded T, ρ, x(H₂)) in the outer shells: −0.06, −0.09, −0.10, −0.42, −2.00e10 erg/g against −0.06, −0.09, −0.15,
−0.53, −2.12 at log m = 3.97, 3.9, 3.8, 3.6, 3.4. (Further in the fit is used below its range and overestimates.)

Phase B is not affected by this defect: its mesh is fixed in every zone (checked on the `fx3` history, all 286
interfaces), so u_rel = u and the gravity work is complete.

Open: in the `fx3` history the gas at 1–6 pc cooled by 9–22 % during the first 500 yr of the accretion phase (steps
23 250–25 350; e.g. 271 → 212 K at log r = 8.0 R☉) at constant density and velocity, and then stayed within 0.1 %
for the next 5400 yr. H₂ cooling is three orders of magnitude too slow for that. Those steps were run with October
code and settings; `full2` will show whether the present code does it after the remap.

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
- `julia -t 8 --project=<env> julia/snapspec.jl runs/<run> [stepmod]` computes emergent spectra of the snapshots (below) and
  writes `export/spec.json` and `spec.f32`; the explorer shows them under "Emergent spectrum".

## What the optical depth is

The τ in the exports, the explorer and the driver (`logtau`, `r_ph`, `T_eff_lum`, `LTAUCUT`) is the grey continuum
optical depth

    τ_R(r_k) = Σ_{j ≥ k} χ_j ρ_j Δr_j ,

summed from the outer mesh boundary inwards, with τ_R = 0 at the boundary: gas beyond the mesh does not count. χ is the
Rosseland mean of the extinction (absorption plus Rayleigh and Thomson scattering) from the equilibrium opacity table,
multiplied in the network zones by the lagged non-LTE ratio for the zone's own composition (`nlte_kappas`). The photosphere
is the outermost interface with τ_R > 2/3 (none while the whole mesh is thinner than that: until step ~9760 in `full4`); `T_eff_lum` is (L/4πr²σ)^¼ there with the run's L.

It is not a line optical depth and not a depth at any one wavelength. The H₂ lines have their own transfer, line by line;
the continuum enters it only through the merging of the zones deeper than τ_R = `LTAUCUT` (10). The monochromatic depth
τ_ν = ∫(κ_ν + σ_ν) ρ dr and the radius where it reaches 2/3 come with the spectra.

## Emergent spectra (`snapspec.jl`)

Post-processing of the recorded snapshots; nothing feeds back into a run. For each selected snapshot (`step % stepmod == 0`,
default 160, plus the last) the continuum L_ν at 216 wavelengths from 0.05 to 1000 μm and the luminosity of each of the 231
H₂ lines, both at the outer mesh boundary ("edge": what leaves the computed volume) and through the run's photosphere
(τ_R = 2/3), and the radius where τ_ν = 2/3 at each wavelength. About 3 s per snapshot on 8 threads; results are cached in
`runs/<run>/spec/`.

- Continuum: static spherical transfer (tangent rays, Feautrier, accelerated Λ-iteration; a direct solve where Rayleigh
  scattering dominates) on the snapshot's ρ, T_gas and network composition. Opacity: collision-induced absorption of H₂,
  H⁻ bound–free and free–free, H bound–free and H⁺ free–free, Rayleigh and Thomson scattering (coherent, isotropic).
  No atomic lines, no helium opacity.
- Emission is thermal at T_gas except in the bound–free continua, which emit by radiative attachment (H⁻) and recombination
  (H) of the zone's own electrons. LTE emission instead made the collisionally ionised gas outside the photosphere of the
  accreting star 16 times brighter than the run in the Balmer continuum. These continua are averaged over each 4.7 %
  frequency band: in cold gas the recombination emission falls off above a threshold within much less than a band, and
  sampled at the grid frequency alone, the solver and the optically thin volume sum disagreed by a factor 6 and the
  spectrum peaked at 0.36 μm instead of 1.5 μm.
- Absorption by the excited levels of H is approximate: the run's three-level populations are not recorded, so the smaller
  of the Saha (from the zone's n_e n_p) and Boltzmann (from the zone's n = 1) populations is used.
- H₂ lines: the run's line transfer (ALI) on the snapshot. No continuum extinction of the lines. The run cools by the
  collisional loss at the ALI populations (`LINECOLL=1`), which where lines are thick is 0.53 to 0.56 of the formal line
  flux; the exported line luminosities are the formal ones scaled to the collisional total (`lscale`), so the sticks add up
  to what the gas loses.
- Zones deeper than τ_R = 300 are replaced by a thermal inner boundary.

The temperature structure is the grey run's, so ∫L_ν dν need not equal the run's luminosity; both are exported and the
ratio is the check.

| Run, phase | ∫L_ν dν / run's L, edge | photosphere |
|---|---|---|
| `full4`, core formation (steps 10240–12000; T_eff ≈ 2050 K, L_edge ≈ 250 L☉) | 0.79–0.82 | 0.80–0.82 |
| `full4`, first frame after the rebuild (step 9760) | 0.87 (109 against 126 L☉) | |
| reference, 146 frames (0.5 to 97 M☉) | median 1.19; within 0.5–2 in 131 | median 0.93; within 0.5–2 in 137 |
| reference, last frame (step 139400, T_eff 6051 K) | 1.01 (3.56e6 against 3.53e6 L☉) | 0.79 (2.66e6 against 3.35e6 L☉) |

Optically thin frames agree with the direct volume sum of the emissivity to 2–3 % (2.89e30 against 2.96e30 erg/s at step
1180 of `full4`). The H₂ line total reproduces the run's `L_lines` to 0.990–1.026 in `full4`.

Where it does not hold: the first ~10 000 steps of the reference's accretion phase (steps < 34 000). The photosphere sits
inside the accretion shock, where τ_R goes from below 2/3 to above 300 within one zone, so the "photosphere" spectrum is
that of the first thick zone (69 times the run's L at step 28800); and the run's flux at the edge is a light-travel
transient there (1.3e9 L☉ at step 24000, negative at step 33600), which a static solution cannot have.

What the spectra show for the collapse: until the continuum transport is on, the continuum is 0.1–0.2 % of the H₂ line
cooling (6e-4 against 0.3–0.6 L☉), and still 0.227 against 516 L☉ at step 9280, just below log n_c 14.

Open: one zone where the run's κ_P is 2.5e-14 against the spectral 3.35e-12 cm² g⁻¹ (energetically negligible).

## The envelope oscillation at 101 M☉ and the stop of `fx3` (9 Nov)

*Mesh tests of 10 Nov show that the amplitude and period of this oscillation are not resolved on the production mesh:
see "Is the envelope cycle at 104 M☉ resolved?" below.*

The production line `fx3` stopped at step 226129, M* = 101.3 M☉ (R* 19.96 R☉, T_c 6.06e7 K, L_nuc 2.8e38 erg/s):
`TITAN step01: failure to converge after 25 attempts`. Continued as `z19fx3d` from the final checkpoint with `NTRY=60`
(now also in `phaseB.env`).

**The stop was the retry limit, not a state without a solution.** The 25 halvings took Δt from 3.5e5 s to 0.02 s with
the last corrections of one zone stuck at −46 % in T and +34 % in E_r for every Δt between 2700 s and 1 s, shrinking only
below 0.2 s. The zone (mesh zone 144, 1789 R☉) had just changed from the network to the table at T = 7600 K with
T_rad = 5416 K; on the table κ_P = 9 cm² g⁻¹ at ρ = 7.4e-12 g cm⁻³, so gas and radiation equilibrate in 1/(c κ_P ρ) ≈ 0.5 s
and no larger step has a nearby solution. Restarts with the unchanged code pass: from the final checkpoint the first step
(Δt = 7.7e-3 s after a restart) changes the zone by 24 % and Δt is back at 5.9e5 s within 40 steps (`t_stop1`); from the
checkpoint 129 steps earlier 400 steps run with 54 divergences, none beyond 10 attempts (`t_stop0`). The path is chaotic,
so the failing step itself does not recur.

**What the envelope does** (old-format history `snaps/snaps_abn02z19_lines_n200_fx3.jls`; tables by
`fx3stop.jl`/`fx3prof.jl`, kept in the session scratchpad, not in the repository):

- From 93 to 100.9 M☉ the flow is steady. L_ph rises from 3.2e6 to 3.6e6 L☉, which is 0.92 to 0.97 of the
  electron-scattering Eddington luminosity (3.77e6 L☉ at 101 M☉ for κ = 0.35 cm² g⁻¹). With the run's opacity Γ = χL/4πGmc
  is 1.03–1.07 locally (235–250 R☉ and at the H ionization zone, 1050–1800 R☉). Ṁ is 7.5e-3 M☉/yr at 10⁴ R☉ and
  6.0–6.7e-3 M☉/yr onto the star.
- The accretion shock stands at about 115 R☉, far above the hydrostatic surface at 21 R☉. Between them lies a hot
  settling layer (u ≈ −20 km/s, T 1e5–3e5 K, ρ 1e-9–7e-9 g cm⁻³, Γ = 0.98, of order 1e-3 M☉).
- Between steps 216900 and 217050 (100.96 M☉) that layer falls onto the star: infall reaches −230 to −250 km/s at
  25–30 R☉, the density there drops 15-fold, and a dense shell (ρ 5e-11 g cm⁻³, u = +1.8 km/s, Γ = 1.18) appears at
  545 R☉. The Γ > 1 shell moves outwards, 513 → 1556 R☉ in 11 yr.
- Once it reaches the photosphere (step ~218500, 101.04 M☉) the region 300–3000 R☉ oscillates with a period of a year or
  less: Γ up to 2.8 at 1300–1900 R☉, Ṁ at 10³ R☉ between −0.26 and +0.16 M☉/yr, gas piling up at 1000–2500 R☉ (ρ up to
  1e-11 against 2e-13 g cm⁻³ in the steady state), the photosphere hopping between zones 141 and 151 (1450–2400 R☉) with
  T_eff 4200–6300 K, and one mesh zone carrying an optical depth of 30 next to one of 0.2. Ṁ onto the star is
  2.6–5.7e-3 M☉/yr. Ṁ at 10⁴ and 10⁵ R☉ has not changed (7.50e-3, 7.34e-3 M☉/yr).

**Not established:** that the star is ending its accretion. The inflow at 10⁴ R☉ is untouched after 45 yr, no energy
budget of the envelope has been done, and the photosphere is not resolved during the oscillation (Δln r = 0.07 against a
jump of 30 in τ per zone). Earlier lines without the star-light and non-LTE terms passed this mass (`_q9`, `_a8`), others
broke down at 107–113 M☉ (`nlo8r`, `_gf1`, `_g9`). Omukai & Palla (2003, ApJ 589, 677) find for steady rates above
4e-3 M☉/yr that the luminosity reaches the Eddington limit before hydrogen ignition and the star expands, possibly
reversing the flow near 100 M☉; this run is above that rate at this mass, but the comparison has not been made
quantitatively.

## From the initial conditions to 0.5 M☉: `full3` and `full4` (9–10 Nov)

Both from-scratch runs finished the collapse (stages A1–A3) and started the accretion phase (`z19full3B1`, `z19full4B1`)
on the evening of 9 Nov. Each died once at the start of stage A3 (17:02 and 17:23 in `procedure.log`) on a failed chemistry
network step in the H₂ dissociation layer; that is the failure fixed by `CHEMFBDENSE=1` (see Procedure). After the 17:27
restart with the fallback, stage A3 ran to the remap without stopping (the fallback was taken 63 times in `full3` and 37
times in `full4`).

| | `full3` | `full4` | October run |
|---|---|---|---|
| step at the remap (M* ≥ 0.5 M☉) | 29773 | 29632 | 23250 |
| M* there | 0.5133 M☉ | 0.5118 M☉ | |
| log n_c | 23.377 | 23.368 | 23.374 |
| T_c | 118 100 K | 116 800 K | 123 800 K |
| Newton divergences in A3 | 1431 | 1479 | |

The central temperature is 5–6 % below October's at the same central density. The cause has not been traced. The remap
onto the accretion mesh itself lowers T_c by 4 % (122 717 → 117 753 K in `full3`) and log n_c from 23.403 to 23.375.

**Where the wall time went** (`perf.csv`, 4 threads; fractions of the summed step time):

| stage | steps | zones | wall | hydro + radiation | chemistry | lines |
|---|---|---|---|---|---|---|
| A1, to log n_c 17.9 (`full3`) | 9900 | 200 | 0.9 h | 10 % | 17 % | 73 % |
| A2, to T_c 5e4 K (`full3`) | 1600 | 260 | 0.6 h | 11 % | 47 % | 41 % |
| A2 from log n_c 14 (`full4`) | 1700 | 260 | 1.1 h | 7 % | 27 % | 66 % |
| A3, to 0.5 M☉ (`full3`) | 17900 | 340 | 5.9 h | 15 % | 27 % | 57 % |
| A3 (`full4`) | 18200 | 340 | 6.0 h | 15 % | 29 % | 56 % |
| B, first 1800 steps (`full3`) | 1800 | 370 | 1.7 h | 18 % | 18 % | 64 % |

The line transfer is the largest cost in every stage. Stage A3 (the adiabatic core growing from formation to 0.5 M☉) is
two thirds of the wall time to the remap. Switching the continuum transport on at log n_c 14 instead of 17.9 (`full4`)
costs half an hour more in A2 and nothing afterwards.

**Accretion phase so far** (10 Nov 01:00; `full3` at 0.66 M☉, `full4` at 0.69 M☉): about 600–1000 steps per hour, with
700 divergences each in the first 1300–1800 steps and up to 28 halvings on one step (`NTRY=60`). The failing corrections
are E_r and F_r in the first zones outside the stellar surface (zones 113–125), where one zone has τ ≈ 900 and T = 8500 K
and the next has τ ≈ 0.008 and T = 3100 K with T_rad = 930 K, and where the zone flag changes from table to network. When Δt
falls to seconds there L_ph falls from about 3000 to about 10 L☉. This is the same unresolved surface that the production
line has at 101 M☉ (previous section).

**Same mass growth as October; a swollen, dim transient from 0.6 to 1.2 M☉** (10 Nov 04:10; an earlier version of this paragraph said the
mass grew at the October pace per step, which was not checked and is wrong). Time is counted from M* = 0.52 M☉; R*, r_ph
and L_ph are medians over seven snapshots (`cmpacc.jl` on the recorded histories):

| M* | | October lineage | `full3` | `full4` |
|---|---|---|---|---|
| 0.60 M☉ | t [yr] / steps | 1.46 / 400 | 1.32 / 1280 | 1.22 / 860 |
| | R* [R☉] | 71 | 82 | 82 |
| 0.80 M☉ | t / steps | 5.2 / 500 | 5.0 / 1900 | 5.4 / 1520 |
| | R* | 75 | 121 | 125 |
| 1.00 M☉ | t / steps | 9.6 / 600 | 9.7 / 2240 | 10.1 / 1880 |
| | R* / r_ph | 84.5 / 95 | 149 / 159 | 154 / 163 |
| 1.20 M☉ | t / steps | 14.6 / 700 | 13.4 / 3720 | 13.9 / 3100 |
| | R* / L_ph [L☉] | 90 / 4500 | 125 / 4600 | 132 / 4450 |
| 1.5 M☉ | t / steps | 21.4 / 750 | 21.4 / 3800 | |
| | R* / L_ph | 108 / 4650 | 117 / 5000 | |
| 2.1 M☉ | t / steps | 34.6 / 950 | 34.8 / 3860 | |
| | R* / L_ph | 129 / 5500 | 125 / 5700 | |
| 2.5 M☉ | t / steps | 45.3 / 1050 | 45.5 / 3920 | |
| | R* / r_ph | 129 / 145 | 132 / 145 | |

The mass grows at the same rate in time (0.05 M☉/yr averaged from 0.52 M☉; 0.023–0.033 M☉/yr through 10⁴ R☉) in all
three. The difference is a transient between 0.6 and 1.2 M☉: the star swells to 1.8 times the October radius while
little luminosity leaves it (the luminosity at the τ = 2/3 interface is between 10 and 2000 L☉ and flickers as that
interface hops across the accretion shock; about 1500 L☉ leaves the shock), and the solver stalls at the shock
(3000–3700 steps where October took 700). At about 1.2 M☉ the photospheric luminosity reaches the October value, the
radius falls back, and from there `full3` follows October: at 2.1 M☉ the time agrees to 0.3 yr in 35, R* to 3 % and L_ph
to 3 %, and 1.2 → 2.5 M☉ took 200 steps. The radius has not been compared with published tracks.

The October run crossed this mass range on 3–4 Nov, before the star-light set (`STARLIGHT`, `SL3`, `NLOPAC`, `SLBLEND`,
`SLEB`, `LYABB`), the flux tie at the outer photosphere (`FLUXTIE`) and the envelope Δt rule (`DTENV`) existed; these
were developed between 36 and 100 M☉ and `phaseB.env` applies them from the remap on. Tests from the remap checkpoint
`z19full3_A3_2` (`t_<tag>` run directories), one switch at a time:

| test | result |
|---|---|
| `FLUXTIE` removed | worse: Δt falls to 1e-8 s at 0.534 M☉, failure to converge after 60 attempts at step 30448 |
| `DTENV` removed | identical to `full3` for 600 steps; after 1200 steps 0.563 M☉ at 0.7 yr with 299 divergences (`full3`: 0.564 M☉, 334) |
| star-light set off (`STARLIGHT=0 SL3=0 NLOPAC=0 SLBLEND=0 SLEB=0 LYABB=0`), with or without `FLUXTIE` | fails on the first step (20 divergences, then a non-finite Jacobian) |
| `LYABB=0` | after 1100 steps 0.562 M☉ at 0.7 yr, 281 divergences: as `full3` |
| `NLOPAC=0 SLBLEND=0` | after 1200 steps 0.568 M☉ at 0.9 yr, 190 divergences: as `full3` |

No single switch removes the stalls of the first year (the tests cover 0.50–0.57 M☉, where the transient begins; none
was run to 1.2 M☉). The flux tie is needed. What produces the swollen, dim state between 0.6 and 1.2 M☉ is not known.

**`full3` crash at 1.12 M☉** (10 Nov 02:38, step 33216): `PopIIIEqEOS.eos: rho = 9.2e-21 below the table`. The run had
stalled at the accretion shock (0.3 yr in 600 steps, Δt limited by T in zone 141, where τ falls from 575 to 0.1 across one
zone). One step needed 18 halvings, which brought Δt to 4.6e-3 s. At t = 7.4e14 s the spacing of Float64 is 0.125 s, so
the step had `timen == timeo`. `mix!` divides the zone mass by the step duration, the compositions of the convective core
came back NaN, and the next equation-of-state call failed. Fix: `mix!` returns without mixing when the step duration is
not positive (the same guard `burn_zone` has since 81105fc). `full3` continues from the autosave `z19full3B1_auto0` as
segment `z19full3B2`; the final checkpoint `z19full3B1` was written after the failure and must not be used. `full4`
failed the same way at 03:13 (step 33178, 1.30 M☉; its process predated the guard) and continues from `z19full4B1_auto0`
as `z19full4B2`. The guard removes the crash, not the stall: by 04:00 the restarted segments had met steps of zero
duration 6 and 8 times and passed them.

## The spurious hot root of the chemistry step: `CHEMSUB` (10 Nov)

`full3` (step 34704, 2.6 M☉) and `full4` (step ~34400, 4.3 M☉) both stopped with `failure to converge after 60 attempts`.
Cause, traced in `full3`: the photosphere zone in the accretion shock (zone 157, r = 256 R☉, n = 9.3e13 cm⁻³, p = 100–250,
below `PSWITCH`) left the table at step 34702; the hydro step 34704 heated it to 6965 K and the chemistry step then took it
to 2.6e6 K, fully ionized, +5.5e14 erg/g. The network integrator (`ChemNet.jl`) is one backward-Euler step in (x, ln T) at
fixed ρ and u with no accuracy control, and the three-level atom's Peebles factor C uses the electron density of the start
of the step (C ≈ 1.7e-3 there). The ionization avalanche has an e-folding time of about 6 s against Δt = 1.1e4 s, and the
single step converged on a root that the resolved evolution never reaches. The energy went into E_r and Newton failed on
the following steps.

`CHEMSUB=fT,fx` (in `phaseA.env` and `phaseB.env` as `CHEMSUB=1.3,3`): a three-level network zone whose single step changes
T by more than the factor fT, or x_e by more than fx (and by more than 1e-4), is integrated in halved sub-steps
(recursively, to depth 16), each with its own start-of-step electron density. With it the same zones (n 5e13–2.4e14,
T 5300–5500 K at the start of the step) end at x_e 2e-4–1.6e-3 with T 4–5 % lower and no step heats (`t_cs1`, steps
34700–34900; `t_cs4` for `full4`, smooth to 4.26 M☉). Where the single step still converges on the right branch it
overshoots x_e by 2–5×. Both chains continue from autosaves before the failure (`z19full3B3` from step 33000,
`z19full4B3` from step 34000).

Not removed by it: the rough spell of `full3` at 1.12 M☉ (steps 33100–33300, Δt 1e1–1e4 s, limited by T in the shock
zone 141) recurs with a similar number of divergences (170 against 161 by step 33200).

Cost on the production state (99.5 M☉): zones 148/149 at the ionization front (1.9e3 R☉, n 7.7e10, x_e 0.46) are
sub-stepped to depth 16 on almost every step, with the same result as the single step (the hydro step cools the zone to
~4800 K, the chemistry step heats it back to ~7000 K); the step costs 0.53 s against 0.49 s. The trigger could skip zones
where the sub-stepped and single results agree; not done.

`CHEMSUB` in phase A is exercised first by `full6`.

## What heats the infalling envelope, and the tension viscosity (10 Nov)

Production state at 99.46 M☉; rates in erg g⁻¹ s⁻¹:

| r [R☉] | T (T_rad) [K] | x_e | κ_E cE | −κ_P 4σT⁴ | artificial viscosity | recombination-cycle heat (`GASQ`) | p dV |
|---|---|---|---|---|---|---|---|
| 2.9e3 | 5419 (3444) | 1.7e-2 | +2.27e6 | −3.22e6 | +6.56e5 | +2.84e5 | +2.8e4 |
| 1e4 | 5418 (1824) | 7.1e-3 | +1.47e5 | −2.04e5 | +5.05e4 | +7.5e3 | +4.5e3 |
| 2.9e4 | 5477 (1043) | 2.1e-3 | +8.9e3 | −1.19e4 | +2.74e3 | +1.2e2 | +7.9e2 |
| 9.7e4 (front) | 4671 (569) | 7e-5 | +52 | −63 | +80 | 0.05 | +110 |
| 2e5 | 1812 (399) | 2e-5 | net −1.5 | | +16.7 | | +13; H₂ lines −18.9 |
| 3e5 | 1401 | | net −0.44 | | +5.4 | | +5.1; H₂ lines −6.4 |
| 1e6 | 716 | | | | +0.11 | | +0.26; H₂ lines −0.18 |

Inside ~5e4 R☉ the gas is heated by absorption of the star's continuum by H⁻ (photodetachment; κ_E/κ_P ≈ 56 because
the diluted light is far hotter in colour than T_rad) and cooled by H⁻ radiative attachment at the gas temperature. The
electrons (x_e 1e-3–1e-2) come from the three-level atom: Lyα is trapped, n = 2 is photoionized by the Balmer continuum.
Direct Lyα bound–bound heating is ≲ 1 %, the recombination-cycle heat 3–9 %, Lyman–Werner heating negligible (H₂ is gone,
which also removes the molecular cooling). Outside the front compression and viscosity heat and H₂ lines cool.

**The viscosity term is numerical.** With `AVTENSION=1` (TITAN's default) a zone whose volume shrinks (div u < 0) keeps its
viscous stress while it is stretched radially (∂u/∂r − u/r > 0), which is the state of smooth free fall; with `AVZONE=2`
that supplied 20–25 % of the heating at 3e3–5e4 R☉ and as much as p dV out to 3e5 R☉. `AVTENSION=0` keeps the viscosity
only where the zone is also compressed radially, that is in shocks. Test `t_avt0` from `z19fx3b` (99.46 → 100.22 M☉, 3000
steps), against production at the same time:

| r [R☉] | T, `AVTENSION=0` | T, production |
|---|---|---|
| 300 | 26647 | 27074 |
| 1e3 | 10521 | 10392 |
| 3e3 | 4174 | 5388 |
| 1e4 | 4367 | 5398 |
| 3e4 | 5059 | 5488 |
| 1e5 | 2003 | 3923 |
| 2e5 | 1469 | 1706 |
| 3e5 | 1234 | 1324 |
| 1e6 | 683 | 685 |

T_rad is unchanged (3231 against 3261 K at 3e3 R☉), x(H⁺) is 3.6e-3 against 1.6e-2 at 3e3 R☉ and 1.5e-3 against 6.8e-3 at
1e4 R☉, the infall slightly faster (−110 against −102 km/s at 3e3 R☉). Δt was 1.09e6 s with one divergence; production
took 14 650 steps (mean 1.93e5 s) for the same interval.

**`AVTENSION=0` is the standard from 10 Nov** (`phaseA.env`, `phaseB.env`; decided with T. Abel: the viscosity is there for
the shocks). `full3` and `full4` keep `AVTENSION=1` through their `run.env`. The test continues past the mass where `fx3`
stopped as `t_avt1`. In the collapse phase the setting is first used by `full6`.

This corrects the 6 Nov note "Precursor temperature" in `POPIII_GOAL.md`: the mechanism (photodetachment heating with
H⁻ above Saha) stands, but about 1000 K of the 5400 K plateau at 3e3–1e4 R☉, and the extent of the warm region beyond
5e4 R☉, were viscous heating.

## Table entry with T_gas = T_rad (`PSWRAD`) and a finer envelope mesh: tests from 99.46 M☉ (10 Nov)

The zone on which `fx3` stopped at 101.3 M☉ had entered the table by the τ rule at 7600 K with T_rad = 5416 K (see "The
envelope oscillation at 101 M☉"). `PSWRAD=ε[,phard]` (default off): a network zone enters the table only if
|T_rad/T − 1| ≤ ε, unless p ≥ phard × `PSWITCH` (default 1: only entries below the pressure switch are tested). Staying
and leaving are unchanged.

Tests from `z19fx3b`, all with `ENVMS0=0.4737` and `AVTENSION=1` (started before the standard changed):

| test | change | state after ~1000 steps |
|---|---|---|
| `t_pr1` | `PSWRAD=0.05` | Δt at the production mean (1.7e5 s), 0 divergences; the test has not yet rejected an entry |
| `t_ms1` | `ENVSPLIT=300,3000,3` (286 → 360 zones; old zones 115–150 split in three) | mean Δt 5.4e4 s, limited by T in the split zones 201–205; 0 divergences |
| `t_mp1` | both | as `t_ms1` |

None reached 100.96 M☉ (onset of the oscillation) before the viscosity standard changed; `t_pr1` was stopped at
100.16 M☉ with 0 divergences. The same switches on the `AVTENSION=0` line are in the next section.
A restart from a checkpoint written with `ENVSPLIT` must not set it again.

## The `AVTENSION=0` line from 100.2 to 101 M☉: the same oscillation, and what a step costs (10 Nov)

*Mesh tests of 10 Nov show that the amplitude and period of this oscillation are not resolved on the production mesh:
see "Is the envelope cycle at 104 M☉ resolved?" below.*

`z19avt0` (the `AVTENSION=0` branch of `z19fx3b`, 100.23 M☉) was continued in four variants, all with `ENVMS0=0.4737` on
4 threads. Numbers from each test's `perf.csv`, state at 08:00 on 10 Nov:

| test | added to `AVTENSION=0` | 100.24 → 100.63 M☉ (smooth) | from 100.63 M☉ (oscillating) |
|---|---|---|---|
| `t_avt1` | nothing | 1500 steps, mean Δt 9.9e5 s, 0 divergences | to 101.08 M☉: 9100 steps, mean Δt 2.5e5 s, 946 divergences, 60 chemistry fallbacks, 42 model yr per wall hour (31 in the last 0.15 M☉) |
| `t_ac1` | `CHEMSUB=1.3,3` | 3900 steps, mean Δt 3.7e5 s, 1 divergence | to 100.77 M☉: 2700 steps, mean Δt 2.6e5 s, 35 divergences, 0 fallbacks, 54 yr/h (`t_avt1` over the same mass: 168 divergences, 46 yr/h) |
| `t_ap1` | `CHEMSUB`, `PSWRAD=0.05` | identical to `t_ac1` line by line to 100.75 M☉ (step 195400): the test never rejects an entry. Stopped. | |
| `t_am1` | `CHEMSUB`, `ENVSPLIT=300,3000,3` | 5900 steps for 0.07 M☉: mean Δt 4.8e4 s, limited by T in the split zones (200–204), 67 fallbacks. 10 yr/h. | not reached |

- **The oscillation is not caused by the tension viscosity.** Without it the envelope leaves the steady state at
  100.63 M☉ (`fx3`: 100.96 M☉). `t_avt1` at 100.87 M☉: the ionized, optically thick layer between 25 and 2000 R☉ has
  Γ = 0.93–1.24, the inflow is brought to rest near 750 R☉, and the gas between 1200 and 2070 R☉ moves outwards at 10 to
  26 km/s. `fx3` at 102.55 M☉: Γ = 0.98–1.00 from 25 to 740 R☉ and a dense cold shell at 2800–3900 R☉. L is
  3.4–3.8e6 L☉ against an electron-scattering Eddington luminosity of 3.82e6 L☉. In both lines the mass of the star keeps
  growing at the rate it had before, 6.2e-3 M☉/yr (`t_avt1` 100.63 → 101.08 M☉ in 71 yr; `fx3` 101.3 → 102.65 M☉ in
  220 yr). No energy budget of the envelope has been made; that the accretion is ending is still not established.
- Before the oscillation the hydrogen ionization front is one mesh zone wide (zones 148 → 149: x(H⁺) 0.84 → 0.012,
  τ 0.74 → 0.05).
- The steady state does not depend on Δt: `t_avt1` (Δt 1.09e6 s) and `t_ac1` (4.3e5 s) at 100.30 M☉ agree to about 1 %
  in the profiles (4 % at zone 148).
- `CHEMSUB` removes most divergences and all chemistry fallbacks at this mass, and costs a factor 2.7 in Δt in the smooth
  phase (cause not traced). `PSWRAD` and the split mesh bring nothing here.

**What sets the step.** In the smooth phase Δt sits on the `DTENV=1000,1e5,0.2` cap: 0.2 of the shortest free-fall time
(r³/Gm)^½ outside 1000 R☉, 1.09e6 s at 100 M☉. In the oscillation Δt is set by the 30 % limit (`STOL`) on the temperature
change of the zones the ionization front crosses (`dtlim` T@143–151), not by the divergences. The free-fall time at
1500 R☉ is 1.4e7 s and the oscillation period is of that order, so the oscillating envelope is followed with about 50
steps per dynamical time. One step costs 0.55–0.67 s (calm: 0.23 s hydrodynamics and transport solve, 0.27 s chemistry,
serial because the star-light ray passes the zones in order, 0.07 s lines; rough: 0.42/0.31/0.30 s).

**Projection.** 300 M☉ at 6e-3 M☉/yr is 33 000 yr more. On the `DTENV` cap that is 9.5e5 steps, about 140 h. At the mean
Δt of the oscillating phase (2.5e5 s) it is 4e6 steps, 600–800 h. `fx3` itself is making 16.5 yr/h at 102.6 M☉ (220 yr in
13.35 h), which would be 2000 h. `t_px1` (the `fx3` state at 102.65 M☉ continued with `AVTENSION=0 CHEMSUB=1.3,3`) showed
that the difference to the 30–50 yr/h of the `AVTENSION=0` line was the cost per step, not Δt; see the next section, which
replaces this projection.

**`CHEMROBUST=1` with `CHEMSUB`** (`t_fb8`, from `z19full4B3_auto1` at 4.8 M☉): 0 fallbacks against about 5 per step in
`full4`, Δt 9e4 s limited by the T change of zone 141, about 25 times less model time per step. This is the result of
6 Nov (`POPIII_GOAL.md`) again; `CHEMSUB` does not change it. The fallbacks of `full4` are at zones 144–149, the shock
layer just outside the photosphere (168–200 R☉, τ 2–0.1, T 5700–4300 K, n ≈ 1e14 cm⁻³).

## Deep in the oscillation: the line transfer, the production switch, and the step at the front (10 Nov)

**`t_px1`** (the `fx3` state at step 272020, 102.65 M☉, continued with `AVTENSION=0 CHEMSUB=1.3,3`, 4 threads): 1780 steps
to 102.73 M☉ with 289 divergences and 11 chemistry fallbacks, mean Δt 2.5e5 s. `fx3` with its own flags over the same
stretch: 15–47 divergences per 100 steps (`t_px1`: 16), mean Δt about 3.5e5 s. So the switch can be made in place at this
mass, and by itself it changes neither Δt nor the pace. Stopped by hand at step 273800.

**Where the time went: the ALI line-transfer call.** In the last 13.7 h segment of `fx3` 38 % of the wall time was
`t_lines`, and 50–240 s of about 240 s per 100 steps towards its end (`t_px1`: 177–228 s of 263 s). The oscillation leaves
a cold molecular shell at 2700–6500 R☉ (T 1300–1800 K, x(H₂) 0.14–0.32, falling at 46–90 km/s). There 100–113 H₂ lines are
optically thick, the velocity span puts each on the maximum of 200 frequency bins, and the adaptive refresh interval of
the full solution (`ALIINT`) drops to 1–8 steps. One call cost 10–14 s at 4 threads. The lines cannot be dropped there:
they are 2–20 % of the net continuum heating in the shell.

**The line solver was made cheaper without changing a bit of its result** (commit 0222371, `H2Lines.jl`): the ray geometry
(path lengths, direction cosines) is built once per call instead of once per line and ray, `transfer_line` keeps its work
arrays for the whole line, the exponential of each ray segment is taken once, and the lines are distributed over threads
one task per line (`@threads :greedy`, since the thick lines are neighbours in the list). All floating-point expressions
are in the original order. One call at step 272600 of `t_px1` (107 thick lines, 200 bins, 2 iterations):

| threads | before | after |
|---|---|---|
| 4 | 14.5 s | 8.3 s |
| 8 | | 3.3 s |
| 16 | 7.1–7.9 s | 1.8 s |

Allocation per call 16.5 GB → 0.32 GB (the allocation was what kept the old version from scaling with threads). Checked:
the edge and photosphere line luminosities, the level populations and the total and collisional cooling are identical
(sha1 of the arrays, and the sums to 17 digits) at steps 272400 and 272600, in the default mode, with `LINESRC=zone` and
with `LINEQUAD=p2`; and a 300-step run from the production state (`t_pt16b` against `t_pt16`, 16 threads) is bit-identical
in all 15 snapshots (T, ρ, u, E_r, species, line cooling) and in every non-timing column of `perf.csv`. In that pair the
wall time per 100 steps went from 89 to 59 s and from 126 to 83 s.

**Production switched to the standard flags at step 273550 (102.69 M☉, 08:41 on 10 Nov).** The segment `z19fx3d`
(`AVTENSION=1`, no `CHEMSUB`, 4 threads) was stopped by hand; `z19fx3e` continues from `z19fx3d_auto1` with
`AVTENSION=0 CHEMSUB=1.3,3 NTRY=60`, 16 threads, the new line solver, and `RUNREC` into `runs/fx3` (snapshots every 50
steps, `perf.csv`). The first 600 steps (to 102.73 M☉): 0.59 s per step, mean Δt 3.3e5 s, 63 model yr per wall hour
against 16.5 before; 82 divergences, 4 fallbacks. Per 100 steps: hydrodynamics and transport solve 22–58 s, chemistry
17–27 s, lines 4–10 s, garbage collection 7–13 s. The solve is now the largest item (0.2–0.6 s per step for about 160
zones, including the repeated steps); it has not been profiled.

**Projection (replaces the one above).** 33 000 yr at 63 yr/h is about 520 h, three weeks, if the pace of the first 600
steps holds; on the `DTENV` cap it would be about 140 h.

**The step at the ionization front** (`t_px1`, steps 272020–273000, `DTDIAG=1`):

- Δt is always limited by the 30 % rule on T or ρ in zones 137–155, the zones the front, the shock and the photosphere
  pass through (1200–5000 R☉).
- Distribution of Δt: below 1e4 s 4 % of the steps, 1e4–3e4 s 11 %, 3e4–1e5 s 18 %, 1e5–3e5 s 35 %, 3e5–1e6 s 31 %, above
  1 %. The steps below 1e5 s are a third of all steps and advance 5 % of the model time.
- The mesh there is fixed in radius with Δln r = 0.07–0.10, zones 150–300 R☉ wide. Gas falling at 50–100 km/s crosses a
  zone in about 1.6e6 s; the `DTENV` cap is 1.08e6 s. The run is at an advective Courant number of 0.15–0.5 in these zones
  even on the cap, so the implicit solver is not being used for steps far above the flow time there.
- The front is one interface wide: T 10 000–12 900 K and x(H⁺) ≈ 1 on the inner side, 1830–2070 K and x(H⁺) = 2–3e-5 on
  the outer (steps 272660, 272700, 272980). The photosphere is in the same zone and moves between 1450 and 2400 R☉ in
  whole-zone jumps within about 7 yr. When the front crosses a zone, that zone goes from neutral to ionized in one or two
  steps, which is what the 30 % rule reacts to.

So time adaptivity is already there globally (`STOL`), and zone-local sub-stepping of chemistry and temperature is what
`CHEMSUB` does. Sub-cycling the hydrodynamics and radiation of the front zones alone could save at most the third of the
steps that are short, and the radiation couples those zones implicitly to all others. What is missing is spatial
resolution that follows the front, so that no zone has to change its state in one step. The code has TITAN's adaptive
grid equation in a window on the fixed mesh (`GRIDWIN=r1,r2,…`). It was tried once, on 8 Nov, together with a threefold
`ENVSPLIT` (`gw5`, from 107.5 M☉ of an earlier line: `ENVSPLIT=1000,8000,3 GRIDWIN=500,8000,1e5,0.3,1,1.5,0.3,0.3,30`);
Δt collapsed to 2e-7 s and the test was stopped without a diagnosis. The static split alone (`t_am1` above) was 5–7 times
slower. `GRIDWIN` alone was tested on 10 Nov (next section).

## The adaptive grid in a window (`t_gw1`–`t_gw3`, 10 Nov)

Strategy. Most of the mass crosses the accretion shock, so the mesh there has to be close to Eulerian; a Lagrangian mesh
loses its zones into the star. The fixed log mesh is that frame. What it lacks is resolution that follows the two thin
layers (shock, ionization front with the photosphere). `GRIDWIN=r1,r2,tau,wd,wt` lets the interfaces between r1 and r2
follow TITAN's grid equation (equal arc length in ln r, √wd ln ρ, √wt ln T, time constant tau) with both window ends
fixed, so mass still flows through.

Tests from `z19avt0` (100.23 M☉, step 189000), 1400 steps, window 800–6000 R☉, against the same steps on the fixed mesh:

| run | tau [s], wd, wt | model yr | mean Δt [s] | steps with Δt < 1e5 s | r_ph mean ± sd [R☉] | L_ph [L☉] | T_eff [K] | Δln r at the photosphere | div. |
|---|---|---|---|---|---|---|---|---|---|
| `t_ac1` (fixed mesh) | — | 16.6 | 3.73e5 | 0.29 | 1788 ± 72 | 3.581e6 | 5943 | 0.071 | 0 |
| `t_gw1` | 1e7, 0.3, 1 | 11.6 | 2.64e5 | 0.39 | 1758 ± 111 | 3.583e6 | 5999 | 0.065 | 2 |
| `t_gw2` | 1e5, 0.3, 1 | 6.6 | 1.47e5 | 0.49 | 1764 ± 77 | 3.591e6 | 5987 | 0.066 | 2 |
| `t_gw3` | 1e5, 0, 10 | 4.9 | 1.14e5 | 0.06 | 1734 ± 38 | 3.571e6 | 6027 | 0.041 | 4 |

The grid moves without the Δt collapse of 8 Nov. Luminosity and T_eff agree to 0.5 % and 1.4 %. With the weight on
temperature (`t_gw3`) the photosphere zone is 1.7 times finer and the photosphere's scatter halves: the front no longer
jumps by whole zones. It does not buy time step: Δt falls in proportion to the width of the front zone, and in `t_gw3` it
is the 30 % temperature rule that sets nearly every step. That rule is applied to the state after the hydrodynamic step
and before the chemistry step (`CHEMSUB`): the hydrodynamic step moves the front zone's temperature at the old
composition, and the chemistry step brings it back. So the limit acts on an intermediate state of the operator split,
not on the change over the step. Ways out, none adopted: apply the 30 % test in network zones after the chemistry step;
`BALSPLIT=1` (the ionization balance inside the hydrodynamic step; it drove an edge flicker on 8 Nov); a chemistry sweep
that follows the flow through the front. `GRIDWIN` stays off in production.

## Is the envelope cycle at 104 M☉ resolved? (`t_pd1`–`t_pd8`, `t_pm1`–`t_pm5`, 10 Nov)

Question: `fx3` shows the photosphere swinging between 1360 and 7300 R☉ (T_eff 1900–8200 K) with a period growing from 3
to 7 yr over 272 yr of record. Is that physical? All tests start from one frozen production state
(`t_ref6/ck/m2c_z19fx3e_auto1`, step 297500, 103.80 M☉, t0 = 2.18640451e7 yr) with the production flags
(`AVTENSION=0 CHEMSUB=1.3,3`), 8 threads. Diagnostics per snapshot by `pulse2.jl`, events by `events.py`, the energy
budget by `bounce.py` (session scratchpad, not in the repository).

**Answer so far: the loading of the envelope is robust, what happens to the loaded shell is not resolved on any mesh
tried.** This replaces the description of the oscillation in the two sections above ("The envelope oscillation at
101 M☉", "The `AVTENSION=0` line from 100.2 to 101 M☉") where they differ.

What holds on every mesh and step limit:

- With L at 40 R☉ pinned at 3.82e6 L☉ = 0.98 of the electron-scattering Eddington luminosity, the infalling gas is
  stopped at 300–3000 R☉ and about 0.9–1.0 M☉ collects there within 0.8 yr (0.1 → 0.9 M☉ between t0 + 1.0 and 1.8 yr in
  every run).
- The energy budget of the envelope (40–10⁴ R☉) closes to 0.4 % of L on the production mesh (`t_pd1`).
- The cycle does not need the H/H⁻ opacity peak: with the flux-mean opacity capped at electron scattering below 5e4 K
  (`CHICAP=0.35`, `t_pd5`) it persists (6.4–8.4 yr, ±20 km/s).
- The star grows at 5.88e-3 M☉/yr (`t_pd1`, 16.9 yr), 5.89e-3 (`t_pd2`, 18.7 yr), 6.32e-3 (`t_pd6`, 18.6 yr).

What differs. A loaded shell does one of two things: a soft rebound (velocity at 1000 R☉ stays above about −60 km/s,
photosphere stays outside 1400 R☉), or it falls nearly freely (−165 to −177 km/s at 1000 R☉; free fall is 199 km/s)
through its own photosphere onto the hot layer and bounces. Collapse events (inflow at 300 R☉ above 0.5 M☉/yr), times in
yr after t0, rates in M☉/yr:

| run | mesh | t | shell mass [M☉] | r_ph min [R☉] | Ṁ max at 300 / 100 / 30 R☉ | L_ph min–max [L☉] |
|---|---|---|---|---|---|---|
| `t_pd1` | 286 zones, fixed (production) | none in 16.9 yr (max Ṁ at 300 R☉ 0.23) | 0.99 | 1478 | | |
| `t_pd2` | same, `STOL=0.1` | 15.82 | 0.97 | 1229 | 1.3 / 0.24 / 0.04 | 5.9e4–1.4e7 |
| `t_pd6` | `ENVSPLIT=300,8000,2` (334) | 1.83 | 0.92 | 869 | 4.1 / 0.54 / 0.08 | 5.5e4–2.5e7 |
| | | 4.92 | 0.96 | 262 | 10 / 9.5 / 1.0 | 1.8e4–2.4e7 |
| `t_pd4` | `ENVSPLIT=300,8000,3` (382) | 1.88 | 0.95 | 302 | 11 / 6.5 / 0.79 | 1.2e4–3.2e7 |
| | | 4.51 (run at 4.55) | 0.68 | 362 | 30 / 52 / 0.36 | |
| `t_pd7` | `ENVSPLIT=800,8000,3` (350) | 1.87 (run ended 2.00) | 0.95 | 369 | 7.8 / 3.9 / 0.54 | |
| `t_pm5` | 286, moving, `gmax` 5 | none in 4.65 yr (max 0.15) | 0.96 | 1379 | | |
| `t_pm4` | 334, moving, `gmax` 2.5 | 1.80 | 0.88 | 1195 | 1.7 / 0.11 / 0.02 | 2.1e5–2.5e7 |
| `t_pm3` | 334, moving, `gmax` 5 | 1.76 (stopped at 2.4) | 0.83 | 1219 | 1.1 / 0.07 / 0.01 | |

Median Δln r of the photosphere zone: 0.075 (`t_pd1`), 0.035 (`t_pd6`), 0.025 (`t_pd4`), 0.024 (`t_pm3`), 0.031
(`t_pm4`), 0.036 (`t_pm5`). The outcome does not converge monotonically with the mesh and is not set by the resolution at
the photosphere alone (`t_pm3` has the photosphere resolution of `t_pd4` and the partial collapse of `t_pd6`). On the
production mesh the photosphere sits on one zone for most of the cycle (d ln r_ph / d ln L_ph = −0.07 to 0.12 on the
fixed 286-zone runs against 0.29–0.47 on the finer or moving ones), so the photosphere and T_eff excursions of `fx3` are
not to be trusted.

**The bounce conserves energy and is a compression rebound** (`t_pd4`, first event; control volume 18–10⁴ R☉, where only
0.0015 M☉ crosses the inner boundary; with the boundary at 25 or 40 R☉ the donor-cell estimate of the advected energy
fails, not the code). Between t0 + 1.65 and 1.97 yr E_grav falls by 20.3e47 erg, E_rad rises by 19.4e47, E_gas by 0.96e47,
E_kin by 0.6e47; the star supplies 1.07e47, 0.33e47 leaves as light. Total energy against the integrated boundary terms:
1.1e46 erg, 0.5 % of the energy converted. So 96 % of the gravitational energy released is stored as trapped radiation and
the rebound is its re-expansion. Profiles: free-falling neutral gas (−226 to −343 km/s, 2400–4700 K, about 2 M☉/yr)
outside 300 R☉, a shock at 228–294 R☉, a hot pile-up at 50–170 R☉ (ρ 1.5–2.6e-7 g cm⁻³, 2.4e5–6.8e5 K); then an outward
shock (294 → 540 R☉ in 0.03 yr, +12 to +120 km/s behind it, local Γ 2–13, L up to 4.7e7 L☉). L inside 100 R☉ stays at
3.7–3.8e6 L☉. The flash at the photosphere (to 3.2e7 L☉) carries 1.6e47 erg, 8 % of the stored energy.

**After a deep event** (`t_pd6`, 18.6 yr): the bounce of t0 + 4.9 yr drives gas outwards at up to 84 km/s at 5000 R☉
(escape speed there 89 km/s). From 6.8 to 12.1 yr mass leaves through 10⁴ R☉ (up to 0.29 M☉/yr) and the mass between 40
and 10⁴ R☉ falls from 0.98 to 0.44 M☉. From 8.9 to 15 yr the shell at 300–3000 R☉ is nearly empty (0.04–0.08 M☉) and the
star accretes steadily at 5.5–6e-3 M☉/yr. The gas returns from 12.6 yr and the shell is loaded again at 18 yr. The cycle
of the finer mesh is therefore longer and irregular, not a 3 yr limit cycle.

**The moving mesh** (`GRIDWIN=40,8000,1e5,0,10 UPFREEZE=4`, without and with `ENVSPLIT=300,8000,2`: `t_pm1`, `t_pm2`)
stalled at t0 + 0.47 and 0.34 yr: about 60 zones with Δln r < 1e-3 collected on the thin radiative shock near 5060 R☉ and
left Δln r = 0.13–0.26 inside 1000 R☉. A ninth `GRIDWIN` parameter `gmax` caps the resolved gradient (zones at most
√(1 + gmax²) finer than in smooth gas); with 5 or 2.5 the runs proceed at Δt of order 1e4–2e5 s outside the events (`t_pm3`–`t_pm5`).

**Cost.** `t_pd4` took 8700 steps for 4.55 yr (Δt 24–300 s inside a collapse, 2.8 s per step at 8 threads); the
production mesh takes about 600 steps per yr at these masses. A resolved cycle is roughly ten times more expensive per
model year. `t_pd8` (`ENVSPLIT=300,8000,6`, about 520 zones) was started to see whether the first event at six-fold
resolution matches `t_pd4`.

**Switches added for these tests (all default off; with all off the code is bit-identical over 600 steps from the frozen
state, `t_bi2` against `t_ref6`):**

- `UPFREEZE=n`: from Newton iteration n on, each interface keeps its upwind direction. Removes a two-cycle of the donor
  choice at interfaces of a moving window where the flow relative to the grid is near zero.
- `STOLCHEM=1`: the step-size rule is evaluated again after the chemistry step, so that it acts on the change over the
  whole step and not on the state between the two halves of the operator split. Under test, not adopted; it ends the
  stall of the resolved-shock runs (see "The mesh of the baseline").
- `CHICAP=κ[,Tmax]`: caps the flux-mean opacity below Tmax. Diagnostic, not physics.
- `NEWTFLIP=1`: prints which branch decisions of the residual changed between Newton iterations.
- `GRIDWIN=…,gmax`: see above.

**Not established:** which outcome (soft rebound, deep collapse and ejection to beyond 10⁴ R☉) the converged solution
has, and how often; whether the deep collapse depends on the shell recombining and losing its opacity (it looks like a
threshold; not tested); whether 1D is adequate for a radiation-supported shell above a lighter layer. None of the runs
shows the accretion ending.

## The cycle has two modes; coarser meshes, tension and the work integral (10 Nov, evening)

Decision of 10 Nov (user): keep the coarse mesh for one full solution to the goal, document the cycle as unresolved, do
resolution studies from restarts at selected masses afterwards, then the final runs. Open: which coarse mesh.

**Two modes of the cycle on the production mesh** (`fx3`, 447 yr of snapshots from 102.7 to 105.4 M☉; cycle = upward
crossing of the mass between 300 and 3000 R☉ through 0.5 M☉):

| | period | photosphere maximum | T_eff minimum | steps per model yr |
|---|---|---|---|---|
| short | 2.3–3.9 yr | 2400–3100 R☉ | 4000–5400 K | 90–140 |
| long | 5.5–10.2 yr (median 7.8) | 4600–7300 R☉ | 1900–3200 K | 125–295 |

Production was in the short mode until 103.85 M☉ (three isolated long cycles before) and has stayed in the long mode
since (27 cycles in a row). The mass between 40 and 10⁴ R☉ grows: 0.57 M☉ (20-yr mean, 190 yr before the frozen
state at 103.80 M☉), 0.94, 1.01, 1.49 M☉ (257 yr after), about 2.1e-3 M☉/yr, while the star grows 5.97e-3 M☉/yr.
Not shown: that the mode change is caused by this growth.

**A comparison of cost or of cycle properties is valid only between runs in the same mode.** The tests from the frozen
short-mode state (`t_co2`, `t_av05`, `t_av20`, `t_st5`) stayed in the short mode for 31–36 yr while production went to
the long mode after 8.9 yr; their 100–125 steps per model yr against production's ~200 is that difference, not a gain.
The `ndiv` column of `perf.csv` is cumulative per process. Single realizations differ (null runs `t_nul1`, `t_nul2`
with `STOL=0.2999`), so a difference needs a null run before it is attributed to a change.

**Tests from a frozen long-mode state** (`frozen_long/m2c_z19fx3f_auto0`, 105.37 M☉, t = 2.18643043e7 yr; production
flags; `t_L0` without changes reproduces production row by row for 400 steps). `ENVMERGE=r1,r2,k` merges groups of k
zones with centres in [r1, r2] R☉ on restart (conserves mass, gas and radiation energy; commits eedb250, 67bc61f).
Tests with 6 threads, production with 16, so seconds are not comparable one to one.

| run | change | zones | model yr | steps per yr | s per step | s per model yr | cycle |
|---|---|---|---|---|---|---|---|
| `fx3` | (production) | 286 | 16.0 | 244 | 0.59 | 145 | long mode, P 7.9 yr, photosphere to 6600 R☉ |
| `t_Lav` | `AVTENSION=0.05` | 286 | 11.1 | 171 | 0.97 | 166 | long mode, same amplitude (mass maximum 1.48–1.50 M☉) |
| `t_Lco` | `ENVMERGE=1e4,1e12,2` | 226 | 20.7 | 140 | 0.68 | 95 | long mode, P 8.1 yr |
| `t_Lm2` | `ENVMERGE=30,1e12,2` | 182 | 34.7 | 107 | 0.54 | 58 | long mode, P 7.9 and 9.5 yr, photosphere to 7000 R☉ |
| `t_t_Lm3c` | `ENVMERGE=30,1e12,3` | 146 | 105.0 | 37 | 0.27 | 10 | **none after the first cycle** |

On the 146-zone mesh (Δln r 0.18 from 30 to 1000 R☉, 0.3 outside 8800 R☉) the envelope goes through one long cycle and
then stays put for the remaining 97 yr: photosphere 6610–6910 R☉, T_eff 3000–3090 K, L 3.4–3.6e6 L☉, inflow at 30 R☉
5.5e-3 M☉/yr, the mass between 300 and 3000 R☉ wandering between 0.15 and 0.32 M☉ with a period near 12 yr. Production
over the 257 yr before: growth 5.98e-3 M☉/yr, mean L 3.70e6 L☉, mean T_eff 4894 K. So the cycle's amplitude grows with
resolution (146 zones: none; 182 and 286: the bounce; 2× and 3× finer: deep collapses onto the star), and none of
these is converged. The star's own mesh and the accretion shock (inside 30 R☉) are the same in all of them.
`t_c146` continues the 146-zone run (8 threads) to see whether the steady state lasts as the mass grows.

**Production moved to the 146-zone mesh (10 Nov, 22:55; user decision: the one coarse full solution runs on 146
zones).** The run directory is `runs/t_c146` (link `runs/c146`). Lineage: `fx3` checkpoint `z19fx3f_auto0` (step
343100, 105.37 M☉, 286 zones) → `t_t_Lm3c` with `ENVMERGE=30,1e12,3` (to step 347050) → `z19c146` restarted without
the switch (to step 387050) → `z19c146b`, … in segments of 100 000 steps (8 threads, standard flags `AVTENSION=0
CHEMSUB=1.3,3`, from 11 Nov also `ETASWTH=1e-3`, see the section on the stop at 117.6 M☉; chain script `chain_c146.sh` in the session scratchpad; the restart at step 387050 continues the record
without a jump). The 286-zone run `fx3` was stopped by hand at step 359500 (105.86 M☉); its last checkpoint is kept
as `frozen_fx3_286/m2c_z19fx3f_auto1` for the resolution studies and for a restart of the 286-zone line.

State of the 146-zone line over its first 1488 yr (steps 347050–388200, 105.93 → 113.35 M☉):

| quantity | value |
|---|---|
| cost (8 threads) | 27.6 steps and 4.8 s of wall time per model yr (0.18 s per step), Δt 1.13–1.16e6 s, 8 divergences in 40 000 steps |
| growth of the star | 5.02e-3 M☉/yr in the mean, falling from 5.5e-3 to 4.6e-3 |
| photosphere | on one mesh interface (5638 R☉, Δln r 0.27) in every snapshot; L 3.52–3.78e6 L☉; T_eff from L and that radius 3330–3391 K |
| centre | T_c 7.09e7 → 8.97e7 K, linear in time (1.3e4 K per yr) |
| burning marks (0.1 % of the median L_ph) | pp chain from 112.33 M☉ (step 381750); deuterium already on at the start of the record; CNO 8e-12 of L |

**The envelope is not stationary on either mesh: it gains mass.** Mean inflow through each fixed interface, from the
change of the enclosed mass over the record (`mstore.jl`, `menv.jl` in the session scratchpad):

| | 286 zones (`fx3`, 533 yr, 102.66 → 105.83 M☉) | 146 zones (`t_c146`, 1488 yr) |
|---|---|---|
| inflow at 1e4–1e5 R☉ | 8.0–8.1e-3 M☉/yr | 6.9–7.0e-3 M☉/yr |
| growth of the star (inside 30 R☉) | 5.94e-3 M☉/yr | 5.02e-3 M☉/yr |
| growth of the envelope (30 to 1.0–1.2e4 R☉) | 2.1e-3 M☉/yr (0.56 → 1.68 M☉) | 1.9e-3 M☉/yr (1.72 → 4.60 M☉) |

So about a quarter of the inflow stays in the envelope, on both meshes at the same rate, and the star receives the
rest. On the 286-zone mesh this mass is spread over 1450–6200 R☉ and takes part in the cycle; on the 146-zone mesh
2.7 of the 2.9 M☉ gained went into the two zones on either side of the photosphere (4290–7540 R☉: 0.54 → 1.13 and
0.72 → 2.86 M☉). The 146-zone mesh therefore gives the same mean mass budget as the 286-zone mesh without the cycle,
and its "steady" photosphere is a growing reservoir held on one interface. Whether that reservoir keeps growing,
drains, or collapses is open, and neither mesh resolves it; that is what the resolution studies from restarts are
for.

**The supply falls with enclosed mass.** Mean inflow through fixed interfaces of the 146-zone record against the mass
enclosed now: 6.6e-3 M☉/yr at 120 M☉ (2.8e5 R☉), 4.9e-3 at 135, 3.1e-3 at 158, 2.4e-3 at 177, 1.3e-3 at 229, 7.1e-4
at 304 M☉ (7.6e6 R☉ = 0.17 pc, free-fall time from rest 67 kyr). The falling growth rate of the star follows the
falling supply (inflow at 1e5 R☉ fell from 7.6e-3 to 6.7e-3 M☉/yr over the record). Reaching 300 M☉ takes of order
60–120 kyr of model time, 3.4–6.7 days of wall time at the present cost.

Tension 0.05 and 0.2 do not damp the cycle in either mode and do not lower the cost; `STOL=0.5` (`t_st5`) gives no gain.

**Work integral** (`t_pd1`, a snapshot every 2 steps; ∮P d(1/ρ) on Lagrangian shells with P = p_gas + f_Edd E_rad;
scripts `workint.jl`, `worksum.py` in the session scratchpad). The moving layer is one persistent shell of 0.93–0.95 M☉:
the same mass shells stay in it over a cycle, the innermost 0.013–0.02 M☉ falls onto the star per 3-yr cycle and
0.03–0.1 M☉ joins from 5000–30,000 R☉. Every shell goes from fully ionized (0.7–1.7e5 K at compression) to fully
recombined (2000–3500 K at maximum radius) each cycle.

| cycle | positive work | by | negative work | by | net | Σ dm max(u²/2) |
|---|---|---|---|---|---|---|
| short, 1.4–4.5 yr | 8.49e46 erg | inner 0.39 M☉ | −3.82e46 | outer 0.54 M☉ | 4.67e46 | 7.8e46 |
| short, 4.5–8.1 yr | 5.41e46 | inner 0.55 | −2.86e46 | outer 0.40 | 2.55e46 | 7.9e46 |
| long, 8.1–16.7 yr | 6.76e46 | inner 0.44 | −3.64e46 | outer 0.49 | 3.11e46 | 9.1e46 |

Radiation pressure does 97–100 % of the positive work. Reading: the hydrogen-ionization opacity valve in the limit
where the whole shell ionizes and recombines (a relaxation oscillation, not a linear mode). Limits: the loops do not
close exactly (up to 40 % in specific volume for the driving shells), the mesh is unresolved, the phase of the luminosity
through the shell was not computed, and He is neutral by construction in the network although the driving shells reach
0.7–1.7e5 K, so a He⁺ contribution cannot appear.

Literature for the comparison, each found on SciX as refereed (existence and status checked, abstracts read for the
first five; quoted model numbers not checked against the papers): Moriya & Langer 2015 (2015A&A...573A..18M),
Nakauchi, Inayoshi & Omukai 2020 (2020ApJ...902...81N), Yoon & Cantiello 2010 (2010ApJ...717L..62Y), Heger et al.
1997 (1997A&A...327..224H), Lobel et al. 2003 (2003ApJ...583..923L), Inayoshi, Hosokawa & Omukai 2013
(2013MNRAS.431.3036I), Inayoshi et al. 2013 (2013ApJ...769L..20I), Cantiello et al. 2026 (2026ApJ..1000L...4C),
Grassitelli et al. 2021 (2021A&A...647A..99G), Yadav et al. 2018 (2018MNRAS.475.4881Y), Grott, Chernigovski & Glatzel
2005 (2005MNRAS.360.1532G), Glatzel & Kiriakidis 1993 (1993MNRAS.263..375G).

## The stop at 117.6 M☉: the total-energy row cools transparent infalling gas (11 Nov)

**What happened.** The 146-zone production run stopped at 00:36:53 on 11 Nov at step 415553 (117.62 M☉) with a NaN
Jacobian. The star was fine (T_c 1.045e8 K, ρ_c 30.3 g/cm³, photosphere at 7536 R☉). The NaN came from one zone at
2.6 pc that had cooled to 17.3 K, far below the CMB (54.5 K).

**Two causes, one behind the other.**

1. *Overflow.* `u_H2(T)` (vibrational energy of H₂, `ChemNet.jl`) was written with `exp(6100/T)`; its derivative
   overflows below 17.3 K and its value below 8.6 K. It is now written with `exp(−6100/T)` for 6100/T ≥ 100
   (bit-identical value and derivative from 17.5 K to 1e9 K).
2. *Why the zone was at 17 K.* With only the overflow fixed (test `t_nan7`) the same zone kept cooling at 0.41 K per 100
   steps, reached 11.7 K and Newton failed again at step 416953. The cooling is not physical. The total-energy row
   (`ETOTAL=1`) minus the kinetic-energy equation minus the radiation-energy row leaves on the gas

       −u [∂(fE)/∂r + (3f−1)E/r + χρF/c],

   which is zero only if the static radiation momentum equation holds. In the transparent outer zones that equation is
   replaced by the closure F = hcE (`FLUXTIE`), f ≈ 1/3, χ ≈ 0, and the bracket is ∂(E/3)/∂r. Wherever E_r has a
   gradient, gas moving at 0.2–0.3 km/s exchanges energy with a radiation field it does not interact with. On the
   146-zone mesh E_r has a slow "front" near 1.5–3 pc (a plateau of 1.87 E_CMB inside, E_CMB outside), and the inflow
   through it loses energy. Check on `t_nan7` (`spur.jl` in the session scratchpad; last 20 snapshots):

   | zone (snapshot index) | measured m·de/dt | predicted from the term | ratio |
   |---|---|---|---|
   | four cooling zones, 1.4–3.5 pc | −2.6e34, −1.85e35, −2.9e35, −1.3e35 erg/s | | 0.98, 1.00, 1.00, 1.02 |
   | nearly all zones 128–145 | | | 0.9–1.3 |

   E_th/E_r of the gas is 1.4e-2 at 0.23 pc, 1e-3 near 0.75 pc and 1e-5–1e-6 beyond 2.5 pc, so a relative error of 1e-5
   in the radiation terms is the whole thermal energy there.

**Fix used to continue (a new switch, default off).** `ETASWTH=x`: a zone also uses the internal-energy (gas-only) row
where 3p/2ρ < x (3p/2ρ + E_r/ρ + K/m), evaluated on the old state like `ETASW` (`TX.ETASW_TH`, `dualswitch` in
`TitanResidFast.jl`; needs `ETASW > 0`). Production runs with `ETASWTH=1e-3` since the restart: zones ≥ 132 of 146
(r ≳ 0.6 pc) are on the gas row. Test `t_nan8` (from the frozen checkpoint `frozen_c146_nan/m2c_z19c146b_auto1`, 3500
steps to step 419000, 118.12 M☉): no divergence, the zones outside 0.6 pc change by 1e29–3e32 erg/s instead of
1e34–3e35, the coldest zone stays at 17.64–17.65 K, and log n_c, T_c, M*, R* and min u agree with `t_nan7` in all printed
digits at steps 416000 and 416900. Commit 508105c. Production restarted at 01:30:23 from `z19c146b_auto1` as
`z19c146c` (54 min down; it starts at step 415550, three steps before the failure). It passed both failure points: at step
417650 no divergence, the coldest zone at 17.65 K.

**What this does not repair.**

- Zones 126–131 (0.1–0.6 pc, E_th/E_r 4e-2 to 2e-3) are still on the total row and get the same term, there a net
  heating of about 1e33 erg/s per zone (4e-3 K/yr at 250–400 K). `ETASWTH=1e-2` would move the limit to 0.27 pc. The
  proper cure is to make the total-energy row consistent with the momentum row that is solved: where `FLUXTIE` replaces
  the radiation momentum equation, subtract u times its residual. Not done (a change of method).
- The zones that were cooled stay cold: 132, 42.7, 17.6 and 215 K at 1.4–3.5 pc (they were 206, 226, 272 and 342 K at the mesh switch).
  The radiation field reheats them at ~3e-5 K/yr only. This gas is 2–5 pc from the star and does not reach it during the
  run; the state was not altered.
- The radiation field outside ≈ 0.3 pc is not the physical one on either mesh. With f ≈ 1/3 and no opacity the moment
  equations leave E_r flat. 146 zones: a plateau of 1.87 E_CMB (T_r 64 K) carrying 7e39 erg/s of the luminosity (L
  falls from 1.525e40 to 8.3e39 erg/s across the front); the physical field is E_CMB + L/4πr²c, 6 % above E_CMB at
  1 pc. 286 zones (`fx3`, steps 273550–359550): a depression instead, E_r/E_CMB falling from 0.92–0.96 to 0.59–0.67
  between 1 and 3 pc; there the gas temperature at 1–5 pc changed by less than 20 K over those 86 000 steps (162 → 160 K
  at 1.0 pc, 215 → 223 K at 2.3 pc, 325 → 344 K at 4.2 pc), so the 286-zone run did not have the runaway cooling. The
  note of 5 Nov in `POPIII_GOAL.md` (F/cE > 1, E 12 % below the CMB) is probably the same defect.
- The from-IC runs `full3`, `full4`, `full6` run without `ETASWTH` (their `phaseB.env` was read before the switch
  existed). Their outer zones were not checked for this.
- Untested: whether the same term is behind the cooling below T_CMB that killed `full1` in phase A with `ETOTAL=1`
  (lines 33–37 and "The total-energy row on the Lagrangian mesh" above attribute it to the missing gravity work).

**Envelope reservoir at the stop** (146 zones, steps 347050–412500, 105.93 → 117.18 M☉; continues the table above):
inflow through 1.2e4 R☉ fell from 7.4e-3 to 6.1e-3 M☉/yr and the growth of the star from 5.5e-3 to 4.2e-3 M☉/yr; the
envelope gains 1.9e-3 M☉/yr throughout (1.75 → 6.26 M☉), 26 % of the inflow at the start and 31 % at the end. 4.2 M☉
of it is in the zone just outside the photosphere and 1.4 M☉ in the zone inside it; that shell warmed from 9900 to
15 800 K.

## The cold zone outside the photosphere and `ETASWTIE` (11 Nov)

The same defect acts at the photosphere of the 146-zone run. The first zone outside it (7536–10 180 R☉) is at 164–168 K
between 17 600–18 400 K inside and 5070 K outside: 237 K at step 355 000, 210 K at 395 000, 197 K at 428 000, 168 K at
449 000. It has been cold since the mesh switch, so it is not a consequence of `ETASWTH`. On this mesh the photosphere
is one interface: the zone inside has an optical depth of 2620, the flux-mean optical depth from the outer boundary to
its outer face is 0.010, so the `FLUXTIE` weight is 1 at every interface from there outwards and every zone outside the
photosphere that is on the total-energy row carries the −u × (moment-equation residual) term.

`ETASWTIE=w` (default 0 = off, commit b8c7d8d): a zone also keeps the internal-energy row when the `FLUXTIE` weight at its
inner interface is ≥ w. This selects the first law by cause. Tests from a frozen copy of production at step 449 550
(122.2 M☉, 8 threads):

| test | setting | steps | first zone outside the photosphere | divergences | Δt [s] |
|---|---|---|---|---|---|
| `t_row1` | `ETASWTH=1e-3` (production) | 750 | 168 K | 0 | 1.085e6 |
| `t_row2` | `ETASWTH=0.1` | 750 | 1100–2200 K, fluctuating | 328 | 6.7e5 |
| `t_row3` | `ETASWTH=1e-3 ETASWTIE=0.5` | 1500 | 5240 K | 0 | 1.085e6 |
| `t_row4` | `ETASWTIE=0.5` | 1500 | 5240 K | 0 | 1.085e6 |
| `t_row5` | `ETASWTH=1e-3 ETASWTIE=0.5` | 20 000 | 5229 K | 0 | 1.074e6 |

`t_row2` is rough because the zone changes row from step to step (its thermal fraction crosses the threshold); with the
selection fixed by the closure weight the roughness is absent. `t_row3` and `t_row4` agree in every zone from the
photosphere to the edge, so on this mesh `ETASWTIE` covers what `ETASWTH` covers. `t_row5` ended at step 469 500
(681 yr) with 124.613 M☉; production at the same step has 124.614 M☉, the same Δt (1.0743e6 s) and a minimum velocity
of −197.1 against −196.9 km/s.

Temperatures at step 469 500, `t_row5` against production (zone outer radius):

| r [R☉] | 1.0e4 | 1.4e4 | 1.9e4 | 4.6e4 | 6.2e4 | 8.3e4 | 1.1e5 | 1.5e5 | 2.1e5 | 3.8e5 | 9.3e5 | 2.3e6 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `ETASWTIE=0.5` | 5229 | 5620 | 5535 | 5537 | 4466 | 3165 | 2341 | 1884 | 1618 | 1237 | 884 | 606.7 |
| production | 164 | 5077 | 5396 | 5557 | 5691 | 5699 | 5017 | 3207 | 2075 | 1530 | 980 | 633.1 |

So the artefact has both signs: it cools the zone at the photosphere and it heats the infalling gas between 5e4 and
a few 1e6 R☉, where the physical rates are small (the budget table in "What heats the infalling envelope": ±50–100 erg g⁻¹ s⁻¹
at the front). The warm plateau of the 146-zone production extends to 1.1e5 R☉; with the internal-energy row the front is
at about 6e4 R☉. The zone at 6.2e4 R☉ was still cooling at the end of `t_row5` (4542 → 4466 K in the last 2000 steps), the
zones outside it changed by ≤ 1.2 % in that interval. The radiation energy density from the photosphere to 0.1 pc is 2–8 % higher. Beyond 0.6 pc
the two runs are identical (both use the internal-energy row there), including the zones at 132, 42.7 and 17.76 K.

Against the last snapshot of the 286-zone run (`fx3`, step 359 550, 105.9 M☉, one phase of a pulsating photosphere), at
fixed radii:

| r [R☉] | 3e4 | 6e4 | 1e5 | 2e5 | 3e5 | 1e6 |
|---|---|---|---|---|---|---|
| `fx3`, 286 zones | 3402 | 2104 | 1703 | 1336 | 1159 | 707 |
| production, 146 zones (step 479 500) | 5481 | 5710 | 4916 | 1984 | 1609 | 882 |
| `t_row5`, 146 zones | 5574 | 3892 | 2291 | 1528 | 1277 | 807 |

The internal-energy row moves the coarse mesh towards the fine one; a difference remains (mass, phase and resolution
are not separated here). The residual of a discretised moment equation grows with the zone size, which is why the
146-zone run shows the artefact far more strongly than the 286-zone run did.

Production runs with `ETASWTH=1e-3 ETASWTIE=0.5` since 11 Nov 06:06 (segment `z19c146e`, from step 516 000, 129.76 M☉;
decision of T. Abel, see the end of the next section). `full3`, `full4` and `full6` use neither switch.

## Are the cold outer zones physical? Thermal balance at 0.1–5 pc, radiative pumping and HD (11 Nov)

Question: the cooling time at large radii could be a little shorter than the infall time, so the cold zones might be real; if so, a
heating term (collisional de-excitation of H₂ whose levels are populated by the star's radiation) must not be missed, and HD would
take the gas lower still. Checked with `outerbudget2.jl` / `outerbudget3.jl` (scratch) on the first 146-zone snapshot (step 347050,
before the artefact acted) and the latest one (step 504050).

**Balance.** Between 0.1 and 3 pc the gas is close to a balance of compression heating (−p ∇·u) and H₂ line cooling: at step 347050
the cooling time is 1.4–3 infall times (r/|u|) and the line cooling is 0.3–0.8 of the compression heating, at 190–380 K. So the
regime is marginal, as suspected, but the balance sits at 240–450 K (temperature at which the code's own H₂ cooling equals the
compression heating), not at tens of K.

**The cold zones are not on that balance.** In the zones at 1.9, 2.6 and 3.5 pc (now 132, 42.7 and 17.8 K; 206, 226 and 272 K at
the mesh switch) the H₂ cooling is 0.15, 2.5e-4 and < 1e-6 of the compression heating. Nothing in the model cools gas at 43 K or
18 K; both are below `TCMB` = 54.5 K. Together with the earlier evidence (the cooling equalled the −u × moment-residual term to 2 %,
the 286-zone run kept these zones within 20 K, the cooling stopped when the row was switched) they are numerical. They stay cold
because the compression heating time there is ~1e7 yr.

**Radiative pumping.** The gas term of the line transfer is the net collisional exchange with the levels (`LINECOLL=1`), so heating
by collisional de-excitation of radiatively populated levels is included for the H₂ line radiation itself. Not included: the line
transfer has I = 0 at the outer boundary (no CMB in the lines) and no continuum (no diluted star light). Evaluated here by adding
J = B_ν(T_CMB) + W B_ν(T_eff) to the statistical equilibrium (optically thin, W from the photosphere): the exchange changes by
< 0.5 % in every zone above 130 K. In the 42.7 K zone it becomes a net heating of 1.7e-29 erg cm⁻³ s⁻¹ (compression: 5.4e-27); the
star adds 2 % to the CMB term (dilution W ~ 1e-6 at 0.1 pc to 1e-9 at 3 pc). The term matters only for gas near or below T_CMB, and for HD.

**HD (estimate, not a run).** D⁺/D and HD/D from k50–k53 of ChemistryKernels.jl with the snapshot's H⁺ and H₂, D/H = 2.5e-5, HD
cooling from `HDlte`/`HDlow` minus its value at T_CMB. On the warm branch y_HD = 1–2.5e-7 and HD adds 20–45 % to the H₂ cooling inside 1 pc
(70 % at 187 K, less than 15 % beyond 2.5 pc); the balance temperature drops by 2–34 K (e.g. 333 → 299 K at 0.3 pc, 266 → 236 K at 0.57 pc). With all deuterium
in HD a second, cold balance exists at 58–64 K in the zones between 0.2 and 2.6 pc, reached only from below an unstable threshold
at 77–114 K; with HD formation limited to one infall time it does not exist. So in this envelope HD is a temperature correction of at most about 10 %, not a runaway, unless something first takes the gas below ~100 K. An artefact that does exactly that would, with HD in
the network, be locked in by real cooling: one more reason to remove the artefact before adding HD.

**Decisions (T. Abel, 11 Nov 05:45): all four.** (1) production with `ETASWTIE=0.5`; (2) reset the artefact-cooled zones at that
restart; (3) D, D⁺, HD and HD cooling behind a default-off switch, side run first; (4) the CMB and the diluted continuum as incident
radiation of the line transfer. (1) and (2) are done; (4) is in the code as the default-off switch `LINEINC` (next section); (3) is in the code as the
default-off switch `HDCHEM` (the section after it). Neither switch is on in production.

**The reset (`TRESET`, commit 397801f).** History of the outer zones on the 146-zone mesh, first snapshot (step 347 050, the mesh
switch) → step 511 900:

| zone | r_out [pc] | T at the switch [K] | T at step 511 900 [K] |
|---|---|---|---|
| 133 | 1.04 | 188.7 | 186.9 |
| 134 | 1.41 | 214.5 | 194.1 |
| 135 | 1.90 | 206.1 | 132.0 |
| 136 | 2.57 | 225.7 | 42.7 |
| 137 | 3.47 | 271.9 | 17.9 |
| 138 | 4.68 | 341.7 | 215.2 |
| 139 | 6.32 | 389.6 | 403.0 |

Densities and H₂ fractions did not change, and all temperatures have been constant since `ETASWTH=1e-3` came in (step ≈ 427 000).
The artefact moved outward zone by zone, so besides the three cold zones it had also taken 127 K from zone 138 (still falling when
the row was switched) and 20 K from zone 134 (within ~700 yr, far shorter than its thermal time). Five zones were therefore reset,
to their temperatures at the mesh switch, the last state before the artefact acted on this mesh (not to a computed balance
temperature): `TRESET=134:214.51,135:206.11,136:225.73,137:271.87,138:341.70`. `TRESET=k:T,…` (with `RESTART`) sets the gas
temperature of the listed zones once after the checkpoint is read; density, velocity, radiation field and composition are untouched
and the held line-cooling rate of the zone is cleared. It is a repair, not conservative: it added 3.11e46 erg of gas energy
(1065 M☉ of gas in the five zones). Zones 130, 131, 139 and 140 (3–4 % warmer than at the switch) and 133 (1 % cooler) were left alone.

Test `t_rst1` (frozen production state at step 512 500, 129.4 M☉, `ETASWTH=1e-3 ETASWTIE=0.5` and the `TRESET` above, 3000 steps
alongside production): 0 divergences, the same time step (1.053e6 s), M* = 129.71 M☉ at step 515 500 in both, the reset zones
hold their temperatures to 0.1 K, zone 105 at 5233 K (20 100 K inside it, 5632 K outside). Production was then stopped at step
516 000 and restarted from `z19c146d_auto0` as `z19c146e` with both switches and the one-time `TRESET` (later segments restart
without it).

## Incident radiation in the H₂ line transfer (`LINEINC`, 11 Nov)

Until now the line transfer had no incident radiation: lines were formed against a dark sky and a dark star. So gas colder than the
CMB still lost energy through its lines, and the level populations never saw the star's continuum. `LINEINC` (default `0`: the code
path is bit-identical to before) adds both:

| `LINEINC` | incident field at each line frequency |
|---|---|
| `0` (default) | none |
| `cmb` | B_ν(T_CMB) from all directions (`TCMB`, 54.5 K in this run) |
| `1` | B_ν(T_CMB), and B_ν(T_eff) inside the solid angle of the photosphere (R_ph and T_eff of the continuum transport, the same numbers `STARLIGHT` uses) |

**How.** `H2Lines.incident` builds, per line, the background mean intensity J_bg = B_c + W(r)(B_s − B_c) with the dilution factor
W = ½[1 − √(1 − R_ph²/r²)] (zone value: mean of the two faces; W = ½ on the photosphere). Thin lines (escape probability and the
thin branch of the ALI loop) add β·J_bg to the mean intensity of the rate equations. Thick lines solve the formal solution for the
difference I − B_bg along every ray (`transfer_line(…; bg)`), with B_bg = B_s on the part of a ray that left the photosphere and B_c
elsewhere; rays that hit the photosphere restart there (the star is opaque), and J_bg is added back. Nothing inside the photosphere
is changed. The gas term stays the collisional exchange at the converged populations (`LINECOLL=1`), which is now negative
(heating) wherever the radiation holds the levels above the kinetic temperature. Between two ALI calls the rate of a zone is
normally the EP rate scaled by the last ALI/EP ratio; with a background that ratio is meaningless where absorption is a large part
of the EP balance, so zones where the EP net and the EP emission differ by more than half are held at their ALI value until the
next ALI call (`popiii_m2c.jl`, `linecool_full`). `LINEINCDIAG=1` prints T_bg, T_star, R_ph, the photosphere node and the number of
held zones at every ALI call.

**Checks (off line, snapshots of `c146` at 130.5 M☉ and of `full3` at 9.85 M☉).**

- Switch off: Λ, the luminosity, the line table and the populations are `==` to the previous code in both modes.
- Thermodynamic equilibrium: gas at T0 inside a cavity at T0 (with or without a star at T0) exchanges nothing. Residual
  |ΣVΛ_coll| relative to the unirradiated value ≤ 4.8e-8 on `c146` (all lines thin) and ≤ 4.1e-4 on `full3` (27–92 thick lines at
  300–2500 K), the same in the EP and ALI modes.
- The floor: with the zones outside 0.5 pc set to 40 K, they still cool without the field (+6.8e27 erg/s) and are all heated with the
  CMB (−1.96e29 erg/s; 15 of 15 zones on `c146`, 46 of 46 on `full3`).
- Size at the present state. `c146`: every H₂ line is optically thin; the collisional exchange changes by −0.18 % in total
  (7.423e34 → 7.410e34 erg/s), by 0.16–0.21 % in the molecular zones at 0.003–0.07 pc, by at most 0.1 % outside 0.07 pc and by
  12 % in the first zone outside the photosphere, which holds no molecules (2y_H₂ = 2e-7). `full3`: −0.28 % in total.

**Side runs from the frozen production state at step 526 000 (130.83 M☉), 3000 steps ≈ 100 yr each.** `t_inc0` (code before the
change), `t_inc2p` (new code, switch off), `t_inc2` (`LINEINC=1`).

- `t_inc2p` = `t_inc0` in T, u, ρ, composition, Λ and E_r in all 60 snapshots, bit for bit.
- `t_inc2` against `t_inc0` at step 529 000: the same M* (131.147 M☉), time step (1.0472e6 s) and R*, T_c equal to 1e-6; 0 divergences in both; the gas is
  warmer by at most 1.05 K (0.05 %), at the H₂ formation front (zones 112–114, 2000–4100 K), by 0.01–0.6 K between 0.003 and 0.04 pc
  and by less than 0.001 K outside 0.1 pc; total line cooling −0.16 %. Star: T_eff 3092–3094 K at 7536 R☉.

So at 130 M☉ the incident field is a small correction, as estimated; its job is the floor (no line cooling below T_CMB) and the
pumping once the star is hotter or the gas colder. Production still runs with `LINEINC` unset.

**What it does not do.** The energy the lines take from the continuum is not debited from the grey radiation field (at most the
0.2 % above). Continuum absorption between the photosphere and the line-forming zones is neglected (τ = 0.03 on `c146`). The star is
B_ν(T_eff) at every line frequency, although on the 146-zone mesh the photospheric zone is at 20 100 K. In the formal solution the edge of the stellar disc
is resolved to one zone width. Held zones have no temperature slope between ALI calls. `snapspec.jl` does not pass the incident
field yet (spectra of a `LINEINC` run would be computed without it).

## Deuterium and HD (`HDCHEM`, 11 Nov)

The network had no deuterium chemistry (deuterium existed only in the nuclear network). `HDCHEM=1` (default off: the code path is
bit-identical to before) adds D, D⁺ and HD and optically thin HD cooling, ported from ChemistryKernels.jl.

**Method.** Each zone carries two extra numbers, x_D⁺ and x_HD per H nucleus, appended to the zone's abundance vector behind the
species of the H network, so the advection and every remap move them with the gas. x_D = y_D − x_D⁺ − x_HD, with y_D = D/H from `DH`
(2.5e-5, the variable the nuclear network uses). Deuterium is a trace: it does not enter the H network, the electron balance or the
EOS. After the chemistry–temperature step of a zone, the three-state linear system D ⇌ D⁺ ⇌ HD ⇌ D takes one backward-Euler step at
the zone's new temperature and composition, with the rates k50–k56 of ChemistryKernels.jl (H⁺ + D ⇌ H + D⁺, H₂ + D⁺ ⇌ HD + H⁺,
H₂ + D ⇌ HD + H, D + H⁻ → HD + e), D ionisation and recombination with the H rates, and HD photodissociation with the zone's
Lyman–Werner rate. The step is evaluated as the stationary state of the chain with a restart to the old state at rate 1/Δt, which
uses sums of positive products only: it is exact for any rate × Δt and conserves y_D. Table zones and hot zones take the stationary
state. HD cooling is n_HD [W(T) − W(T_CMB)] with W = Λ_LTE/(1 + Λ_LTE/(n_HI Λ_low)) (the `HDlte` and `HDlow` fits of
ChemistryKernels.jl). It sits inside the implicit chemistry–temperature step and uses the x_HD of the start of the step.

| switch | meaning |
|---|---|
| `HDCHEM=1` | D, D⁺, HD and HD cooling on (default `0`) |
| `HDINIT=eq` (default) or `HDINIT=<x_HD>` | state of zones that have no deuterium entries yet (restart from a run without HD): the stationary state, or the given x_HD with D⁺/D = H⁺/H |
| `DH` | D/H by number, 2.5e-5 |

Snapshots of an `HDCHEM` run carry a field `hd = (; Dp, HD, Λ)` (x_D⁺, x_HD, Λ_HD per volume); `runrec.jl` reads runs with and
without it. A run with the switch off that restarts from an `HDCHEM` checkpoint drops the two entries.

**Offline checks** (`julia/hdchem_check.jl`). The ported rates agree with ChemistryKernels.jl to 6e-14 over 10–1e6 K, the cooling fits to 2e-13.
The chain step against a BigFloat solve of the backward-Euler system, 20 000 random cases (rates 1e-30–1e10 s⁻¹ including zeros,
Δt 1e-10–1e20 s): largest relative error 6.7e-16, |Σp − 1| ≤ 2.2e-16. The long-step limit is the stationary state, and repeated
steps converge to the matrix exponential. The cooling vanishes at T = T_CMB and heats below it.

**Side runs** (3000 steps = 96 yr from step 526 000, 130.83 M☉, production flags; all ran clean to 131.147 M☉ with no divergence
and no measurable cost: 582–594 s with the switch on, 583 s with it off).

- `t_hd0`, patched code with the switch off: every snapshot equal bit for bit to the unpatched reference `t_inc0`.
- `t_hd1b`, `HDCHEM=1 HDINIT=eq`: HD radiates 6.50e33 erg/s, 9.3 % of the H₂ lines (6.98e34 erg/s). Largest temperature change
  0.18 %: −4.6 K at 4085 K, −3.7 K at 2015 K, −1.0 K at 1314 K, below 0.1 K outside 0.07 pc (the run is far shorter than the cooling time out there, which is 1.4–3 infall times).
  Λ_HD/Λ_H₂ by radius: 1–2 % at 2e-4–0.03 pc, 3–10 % at 0.04–0.17 pc, 14 % at 0.23 pc, 23–36 % at 0.3–0.8 pc, 72 % at 1.04 pc
  (187 K), 5–28 % at 1.4–3.5 pc, below 2 % outside 4.7 pc. x_HD is 4–8e-8 at 0.003–0.17 pc and peaks at 2.45e-7 at 1.04 pc: at most 1 % of
  the deuterium is in HD, because the gas is only 1e-3 molecular out there. In the hot zones all deuterium is D⁺.
- `t_hd2b`, `HDINIT=1e-9`: inside 0.1 pc the same state as the equilibrium start within the 96 yr; outside 0.13 pc x_HD stays near
  the seed (the relaxation time there is longer than the run, see the table).
- `t_hd3`, `t_hd4`: restart from the `HDCHEM` checkpoint with the switch on (entries carried) and off (entries dropped), both clean.
- `t_hd5`, `HDCHEM=1 LINEINC=1` together: clean to the same step and mass (620 s). The two changes add: outside the star the
  temperature differs from the sum of the two single changes by at most 0.003 K (−3.5 K at 4085 K, −2.7 K at 2015 K, −0.8 K at
  1314 K); H₂ line cooling 6.969e34 erg/s (−0.44 %), HD 6.505e33 erg/s.

**Which start is right.** The relaxation time of the deuterium chain (the slowest non-zero eigenvalue of the rate matrix, `relax.jl`)
against the infall time r/|u|, last snapshot of `t_hd2b`:

| r [pc] | n_H [cm⁻³] | T [K] | stationary x_HD | t_relax [yr] | t_infall [yr] | ratio |
|---|---|---|---|---|---|---|
| 0.0034 | 7.1e7 | 2012 | 8.1e-8 | 7e-5 | 1.6e2 | 5e-7 |
| 0.021 | 5.6e6 | 945 | 5.5e-8 | 1e-2 | 3.3e3 | 3e-6 |
| 0.094 | 3.9e5 | 479 | 4.4e-8 | 7.4 | 5.3e4 | 1e-4 |
| 0.23 | 4.5e4 | 328 | 5.1e-8 | 1.5e3 | 2.2e5 | 7e-3 |
| 0.42 | 1.0e4 | 260 | 1.2e-7 | 2.0e4 | 5.3e5 | 0.04 |
| 1.04 | 1.2e3 | 187 | 2.5e-7 | 1.3e5 | 1.9e6 | 0.07 |
| 1.9 | 2.9e2 | 206 | 1.9e-7 | 1.7e5 | 5.0e6 | 0.03 |
| 4.7 | 37 | 342 | 6.6e-8 | 1.1e5 | 2.0e7 | 0.005 |
| 11.5 | 5.1 | 601 | 2.0e-8 | 7e4 | 1.1e8 | 6e-4 |

The chain relaxes at the sum of the forward and backward rates, 14 to more than 1e6 times faster than the gas falls in (the formation time
alone, 1e7 yr at 1 pc, is the wrong measure). So the infalling gas is in deuterium equilibrium everywhere, `HDINIT=eq` is the
state to restart with, and the numbers of `t_hd1b` are the relevant ones. What 96 yr cannot show is the thermal response at
0.2–3 pc, where HD adds 14–72 % to the cooling and the cooling time is 1.4–3 infall times (3e5 yr and more): that needs a long
stretch of the run.

**What it does not do.** Chemistry deuterium is not coupled to the nuclear network's deuterium (burning in the star does not lower
y_D here; the star's zones are fully ionised, so their deuterium is D⁺ either way). HD–H collisions only, no HD–H₂ or HD–He. The
rate fits are not in exact detailed balance. HD lines are optically thin and are not pumped by the star or the CMB beyond the
W(T_CMB) term. HD is not in the EOS. HD photodissociation uses the zone's H₂ Lyman–Werner rate, with the H₂ shielding. In the
5000 K zones the fits are used at the edge of their range (the temperature is clamped to 10–1e6 K). The explorer does not show
x_HD yet.

## Cloud mass, the initial force balance and dark matter (`full5`, `full6`, 10 Nov)

`MCLOUD` (default 1e4 M☉) sets the gas mass on the mesh. The ABN02 read-off n(r) ends at 93 pc, which encloses 1.5e4 M☉;
beyond, the builder continues the power law. T(M) and v(M) reach 2.2e4 M☉. `MCLOUD=2e4` puts the edge at 166.9 pc
(T_edge 786 K, ρ_c/ρ_edge 1.6e5, 200 shells with mass ratio 1.12), with the outer quarter by mass on the extrapolated
density.

Dark matter is included (`DM=1`, default): the enclosed dark mass enters the hydrostatic interior and the momentum
equation (`EXTMASS`). M_DM/M_gas is ≈ 77 inside 0.01 pc (a slope −1 continuation inwards, not a measured profile), 6.5 at
0.5 pc, 3 at 2–4 pc and 15–22 outside 16 pc. The initial transient is not missing gravity. The model is hydrostatic only
inside the splice at 4.77 pc; outside, the read-off profile is used as it is. Pressure gradient over gravity (gas + dark
matter) in the first snapshot with `MCLOUD=2e4`: 0.985 at 3.5 pc, 0.79 at 7.1 pc, 0.075 at 16 pc, 0.40 at 27.6 pc, 1.00 at
49 pc, 1.58 at 130 pc. So the gas between 7 and 30 pc falls (with the imposed −0.3 to −1 km/s), and the extrapolated outer
part is over-pressured and expands against the constant-pressure boundary. Measured radial profiles of the ABN02 halo
(gas, temperature, velocity, enclosed dark mass) to the virial radius would remove both the inner dark-matter guess and
the outer extrapolation.

`full5` (`MCLOUD=2e4`, radiation from log n_c 14, `AVTENSION=1`) was stopped at step 3200 (log n_c 4.45) when the viscosity standard
changed; up to there log n_c followed `full3` to 0.005. `full6` is the same run with `AVTENSION=0`
(`A1X="NSTOP=1e14" A2X="DTRADMAX=0"`, `runs/full6/run.env` = `MCLOUD=2e4`).

`full6` to 08:00 on 10 Nov, the first from-IC run with `AVTENSION=0` and `CHEMSUB=1.3,3`: stage A1 to log n_c 14 in
9171 steps (40 min, 0 divergences), stage A2 to T_c = 5e4 K at step 10739 (55 min, 3 divergences, n_c 4.9e22 cm⁻³), 0
chemistry fallbacks in both; stage A3 (core growth) running, 0.004 M☉ at step 11400 with 0 divergences. Core formation at
t = 2.27067718e7 yr (`full4`, `MCLOUD=1e4` and `AVTENSION=1`: 2.33773179e7 yr).

`full3` and `full4` in phase B (both `AVTENSION=1` by their `run.env`, with `CHEMSUB`): `full3` pass 3 took about 1650
divergences between 1.1 and 1.25 M☉ (steps 33000–35000, Δt 10–1000 s; the one-zone photosphere and shock layer at
160–185 R☉, zones 141–146) and has run freely since, 3.2 M☉ at step 35400 with Δt 6e6 s. `full4` is at 7.8 M☉ (step
35500) with 59 divergences and about 5 chemistry fallbacks per step (see the `CHEMROBUST` note above).

New in `fullrun.sh`: `A1X="<settings>"` (appended to stage A1, as `A2X` for A2) and `<run directory>/run.env` (appended
to the settings of every stage and restart; read when the chain starts).

## The entropy staircase behind the shock in `full6`, and the zone split (10 Nov)

With `AVTENSION=0` the settled gas behind the accretion shock of `full6` (stage A3, 0.04–0.08 M☉) shows steps in entropy
and density. The question was whether a reduced tension viscosity would smooth them. Tests from `full6`'s checkpoint at
step 16020 (0.0359 M☉, 336 zones), 1500 steps each (`AVTENSION=f`: the stress of the tension
zones, div u < 0 with the traceless ∂u/∂r − u/r > 0, is kept times f; `radial`: it is kept only where ∂u/∂r < 0 in the
plain sense; both new, and the standard 0 and 1 are unchanged):

| test | `AVTENSION` | div. | M★ at step 17500 | T at 2 R★ [K] | q/p min | q/p max | largest entropy step, new zones | roughness |
|---|---|---|---|---|---|---|---|---|
| `t_sa0` | 0 (= `full6`) | 0 | 0.0636 | 2470–2720 | 0 | 10–45 | 0.13 / 0.14 | 0.038 / 0.039 |
| `t_sa1` | 1 | 149 | 0.0584 | 3100–11400 | −9.5 … −63 | 0.9–2.5 | 0.40 / 0.37 | 0.153 / 0.134 |
| `t_saf` | 0.25 | 8 | 0.0949 | 3000–8400 | −9 … −16 | 2.6–3.9 | 0.23 / 0.25 | 0.084 / 0.060 |
| `t_sar` | radial | 17 | 0.0537 | 2500–3080 | 0 … −0.2 (spikes −54, −84) | 6–32 | 0.16 / 0.11 | 0.050 / 0.054 |

(Entropy step: |Δln s| between neighbouring zones accreted during the test, at two epochs; s is a p/ρ^γ proxy, not an
adiabatic invariant in partly ionized gas, so the numbers compare runs and are not absolute.) Any tension makes it worse:
the stretched infall ahead of the shock is heated by the viscosity (T at 20 R☉ at M★ = 0.0567: 2899 K with 0, 8791 K with
0.25, 12999 K with 1), the shock spreads, the solver diverges, and the steps behind the shock grow. `AVTENSION=0` stays.

The steps have a strict period of four zones, which is one zone of the initial mesh after its two splits (`REZWIN`), and
the infall carries the pattern before it reaches the shock: in `t_sa0` at step 17500 the difference of ln ρ between
neighbours is −0.13 to −0.14 across the edge of each block of four and −0.056 to −0.081 inside (smooth would be about
−0.07 everywhere). The seed is the split itself. `Rezone.apply` gives the two halves of a split zone the limited slope
`0.5·minmod(Δln ρ to the outer zone, Δln ρ to the inner zone)` per zone index. The split window moves outward in mass, so
a zone is split when its inner neighbour already has been; that neighbour's centre is 0.6–0.75 of a zone away, the
difference to it is that much smaller, minmod takes it, and the halves come out too flat with the missing jump left at
the block edge. A sequence of such splits of a power law (ρ ∝ r^−2.2, mass ratio 0.1078) gives Δln ρ = −0.0515 inside a
block and −0.127 between blocks, where the exact values on that mesh are −0.072, −0.070, −0.069 and −0.071.

`REZSLOPE=lnr` (default off) takes the slopes of ln ρ, e and ln E_r per unit of ln r and multiplies by the distance of
the half centres. In the same sequence it gives −0.0705, −0.0684, −0.0667 inside and −0.0759 between blocks (errors of
ln ρ ten times smaller). In a run (`t_sml` against the control `t_sm0`, both with `REZNMAX=420` from step 16020; the two
blocks split by step 16680, still in free fall): inside the blocks −0.071 to −0.090 against −0.054 to −0.083, between the
blocks −0.106 against −0.131, and to the half-split zone outside −0.135 against −0.177. The sawtooth is reduced, not
removed, on a profile with curvature. Whether the entropy steps behind the shock go with it is open until those zones
have passed the shock.

**Outcome (13:20 on 10 Nov): the split rule does not change the entropy steps behind the shock.** `t_sml` and `t_sm0`
ran 6000 steps, to 0.22 M☉. The gas accreted after 0.1 M☉ was split after step 16020, so with the two different rules.
Compared at M★ = 0.2228 (`stcmp.jl`), `t_sml` / `t_sm0` / `full6`: roughness of ln s over the 25 zones accreted since
0.1 M☉ 0.037 / 0.038 / 0.035, over the 8 zones since 0.16 M☉ 0.048 / 0.050 / 0.045; largest step 0.141 / 0.138 / 0.136;
mean step across block edges 0.066 / 0.068 / 0.068 and inside blocks 0.032 / 0.031 / 0.032. The steps at the block edges
are twice those inside with either rule, so the density sawtooth of the split is not what sets them. What does is open.

A statement made here earlier was wrong: that with `REZNMAX=340` reached at step 16500 nothing beyond 0.24 M☉ is split.
The window keeps splitting as the star grows (the zone at 0.24 M☉ went from 26.9e-3 to 6.7e-3 M☉ between steps 16500
and 18000; by step 24000 everything inside 1.0 M☉ has been split twice), so all the gas that crosses the shock in stage
A3 arrives in 2.8 % zones in blocks of four. Stage A3 ends at 0.5 M☉ with the remap onto the fixed mesh. Two other
attempts did not work: splitting 40 or more zones in one pass by mass
(`t_saz`, `t_saz2`: the radiation rows diverge far outside the split range; `t_saz3`, `REZWIN=1.3,4,0.02` with
`REZNMAX=420`: 60 divergences in 1500 steps) and splitting by Δr/r as well (`REZDR=0.2`, `t_sad`: 44 divergences in 1100
steps, Δt 2e2–1e3 s, and the surface found at 3 R☉ and 0.027 M☉ at step 16500 where the control has 7.4 R☉ and 0.044 M☉).

## Where the time of a step goes: a profile of production (10 Nov)

Per step, from `perf.csv`:

| run | threads | zones | solve | chemistry | H₂ lines | garbage collection | allocated |
|---|---|---|---|---|---|---|---|
| `fx3` (103.6 M☉, 3000 steps) | 16 | 286 | 0.43 s | 0.14–0.24 s | 0.085 s | 0.087 s | 0.4 GB |
| `full4` (8.1 M☉) | 4 | 370 | 0.22 s | 0.86 s | 1.95 s | 0.29 s | 2.7 GB |
| `full6` (0.22 M☉, stage A3) | 4 | 340 | 0.18 s | 0.44 s | 0.69 s | 0.16 s | 1.1 GB |

In production the implicit solve is the largest item; in the runs from the initial conditions the line transfer is.

A sampling profile of 600 steps from the production checkpoint `z19fx3e_auto1` (`t_prf`, 16 threads, 5 ms sampling,
`Profile` around the whole run; flat and tree listings in `runs/t_prf/prof.txt` and `prof.txt.tree`):

- **80 % of all thread samples are idle.** With 16 threads the run uses about three cores' worth.
- Main task, 101.7k samples: the implicit solve 49.6k (`TitanDriver.jl:188`), the chemistry and line closure 46.1k (`:189`), of
  which 27.4k waiting for the threaded chemistry and 14.7k waiting for the line threads.
- **Solve** (serial): the Jacobian assembly is 27k of the 49.6k, the residual 16.7k. The assembly is one residual sweep in
  `ForwardDiff.Dual` numbers with 42 partials (6 colours × 7 equations, `_assemble!(…, Val{42})`). In the first stage of the
  step the material calls inside it (the equation of state, the nuclear rates, the thermodynamic chain of the network gas,
  `_material!` lines 474–479) are 7.0k of 11.7k assembly samples: functions of (ρ, T) of one zone, evaluated with 42 partials.
- **Chemistry threads** (84.7k samples): 18.3k are the start-up build of the star-light photo-rate table in
  `ChemNet.photorates_t`, not a per-step cost (measured afterwards: the table is built once per process, 4.8 s for H⁻, H and
  Li, and a lookup takes 1–15 µs; in the first step all 16 chemistry threads built the same table at once). Of the other 66k:
  the network Jacobian 7.2k; the non-LTE opacity means 8.0k, with `kcia` (`searchsortedfirst`) 8.1k. That is 0.55 CPU-s per
  step, for which the main task waited 0.23 s; evenly spread over 16 threads it would be 0.034 s.
- **Line threads** (160k samples): `H2Lines.transfer_line` (the exponentials; 12.6k in a two-index `setindex!`),
  `H2Levels.populations` 29.8k with `gth` allocating a fresh array per call (5.2k).
- `exp` is 5.2 % of all samples.

A change at rounding level in the Jacobian only, not made: evaluate the material functions with two partials and apply the
chain rule. A GPU does not fit the problem (a few hundred zones, small dense blocks, a time-sequential implicit solve), and
single precision cannot hold it (Δt/t is 1e-10 to 1e-12 at t = 2e7 yr).

### Bit-identical changes, measured (10 Nov)

All of these leave every number of a run unchanged. The check is a 600-step run from a frozen copy of the production
checkpoint `z19fx3e_auto1` (step 297550, 16 threads) against the same run on commit 3564bda (`t_ref6`): all 12 snapshots
bit-identical (T, ρ, u, E_r, F_r, r, χ, species, line cooling), all non-timing columns of `perf.csv` identical, the logs
identical apart from timing strings and the checkpoint name. The same holds for a 300-step pair in stage A3 from the
`full6` checkpoint (`t_a3r`, `t_a3d`, 4 threads, 15 snapshots).

| change | switch | effect |
|---|---|---|
| Chemistry zones handed out one at a time from a shared counter instead of `Threads.@threads` blocks. The expensive zones (sub-cycled, fallback) are neighbours, so a few tasks carried most of the work. | `CHEMSCHED=zone` (default), `chunk` | chemistry 0.29 → 0.13 s per step |
| Material calls of the residual and Jacobian sweep (`_material!`) on tasks, one zone at a time, each task with its own zone cell | `MATTHREADS=<n>` (default 0: serial) | solve 0.31 → 0.22 s |
| MLT coefficients (∇_ad, c_P, δ, κ_R at the interfaces) evaluated in a threaded pre-pass | with `MATTHREADS` | solve 0.22 → 0.18 s |
| Work arrays in `H2Levels.populations` (`PopWS`): no rate-matrix allocations, collision coefficients and Boltzmann factors kept over the ortho:para iteration. Checked directly as well: 4000 random cases bit-identical, 216 → 29 KB and 620 → 500 µs per call. | automatic with `LINESCHED=zone` (default) | lines 0.076 → 0.051 s, allocation 0.32 → 0.26 GB |
| Level populations distributed one zone at a time | `LINESCHED=zone` (default), `chunk` | none measurable |
| Star-light rate of the chemistry step computed without the unused heating integral (`ChemNet.photorate`) | none | none measurable (2 % of the chemistry samples) |
| Photo-rate table built once under a lock | none | 4.8 s × 15 threads of start-up work |

Per step over steps 297700–298100 (the runs shared the machine with each other, load 40–45: read the numbers to 10 %):

| | wall | solve | chemistry | H₂ lines | garbage collection | allocated |
|---|---|---|---|---|---|---|
| commit 3564bda (`t_ref6`) | 0.67 s | 0.30 s | 0.29 s | 0.076 s | 0.096 s | 0.32 GB |
| new defaults (`t_finC`, without `PopWS`) | 0.51 s | 0.30 s | 0.13 s | 0.072 s | 0.088 s | 0.32 GB |
| all, `MATTHREADS=16` (`t_finD`) | 0.37 s | 0.18 s | 0.14 s | 0.051 s | 0.085 s | 0.26 GB |
| stage A3, 4 threads, 3564bda (`t_a3r`) | 1.22 s | 0.17 s | 0.60 s | 0.41 s | 0.042 s | 0.31 GB |
| stage A3, 4 threads, all, `MATTHREADS=4` (`t_a3d`) | 0.78 s | 0.11 s | 0.24 s | 0.39 s | 0.032 s | 0.25 GB |

Production: 1.8 times faster per step. `MATTHREADS` is off by default because only the Pop III material has been checked
for calls from several tasks; `phaseA.env` and `phaseB.env` set it (the number of tasks is the smaller of `MATTHREADS` and
the thread count). What is left in production: the serial Newton iteration (assembly and block solve), the slowest single
chemistry zone of a step, and garbage collection, which is now a quarter of the wall time.

## Threads, pinning and a 100-zone mesh (11 Nov)

**No run is pinned.** The launchers start `julia -t N` and nothing in the code sets an affinity; the 8-thread production
process may run on all 128 cores (`taskset -cp`: 0–127). An unpinned process also creates 64 idle OpenBLAS threads
(80 threads in all; with an affinity mask of 8 cores it creates 20). The machine has 16 NUMA nodes of 8 cores, each node one
shared L3 cache.

**The 146-zone production line ran without `MATTHREADS`.** The threaded material sweep of 10 Nov (section above) is off
by default; `phaseA.env` and `phaseB.env` set it for the runs from the initial conditions, but the environment of the
`fx3` line, which the `c146` chain and the test launcher inherit, never had it. In production the main thread has
therefore used 8 times the CPU time of each of the 8 workers. The switch changes no number (checked again below).

**A 100-zone mesh.** From the frozen 146-zone checkpoint (`frozen_c146_e`, step 526000, 130.83 M☉) with the one-time
`ENVMERGE='0.1,10,2;11,30,2;30,1500,2'` (`t_z100m`): the stellar interior 57 → 29 zones, 11–30 R☉ 16 → 8, 30–1500 R☉
21 → 11; the two zones at the surface, the photosphere and everything outside 1500 R☉ are unchanged. Merging inside the
star breaks the force balance: for about 400 steps the star rings (Δt 40–300 s, infall to 700 km/s, 6 divergences), then
Δt is back at 1.05e6 s with R* = 10.4662 R☉. At step 528000 against the 146-zone run at the same step: T_c 1.3852e8
against 1.3892e8 K, L at the photosphere 4.670e6 against 4.684e6 L☉, T_eff 3091 against 3093 K, nuclear luminosity
3.27e6 against 3.56e6 L☉. That is 1500 steps (50 yr) after the merge, so the interior has not relaxed and the 8 % in
the nuclear luminosity is not yet a statement about the mesh. Checkpoint: `frozen_z100/m2c_z19z100m`.

**Timing.** 3000 steps from a frozen checkpoint with `ETASWTH=1e-3 ETASWTIE=0.5 HDCHEM=1 LINEINC=1`; wall time over the
last 2500 steps (start-up excluded); the four component columns are seconds per 100 steps. Pinned means `numactl
--cpunodebind=n --membind=n taskset -c <cores>` with as many cores as threads. Within each mesh all runs give identical
numbers (every non-timing column of `perf.csv`), whatever the thread count, the placement and `MATTHREADS`; no
divergences. These are single runs on a shared machine, and the load changed between the rounds, so compare within a
block.

100 zones, `MATTHREADS` = thread count (load average 50–56):

| run | threads | placement | s per step | solve | chemistry | H₂ lines | GC |
|---|---|---|---|---|---|---|---|
| `t_tm8p` | 8 | pinned, one node | 0.0765 | 3.7 | 1.7 | 1.6 | 0.7 |
| `t_tm8u` | 8 | not pinned | 0.0851 | 4.1 | 2.0 | 1.7 | 1.2 |
| `t_tm4p` | 4 | pinned, alone on the node | 0.0966 | 4.2 | 2.3 | 2.5 | 0.7 |
| `t_tm4u` | 4 | not pinned | 0.0970 | 4.3 | 2.4 | 2.4 | 0.8 |
| `t_tmA4a`, `b` | 4 | two runs on one node | 0.103, 0.102 | 4.4 | 2.5–2.6 | 2.7 | 0.8 |
| `t_tm2p` | 2 | pinned, alone on the node | 0.136 | 5.5 | 3.5 | 4.1 | 0.8 |
| `t_th8p3` (no `MATTHREADS`) | 8 | pinned, one node | 0.107 | 6.3 | 2.0 | 1.9 | 0.8 |
| `t_th4p3` (no `MATTHREADS`) | 4 | pinned, alone on the node | 0.123 | 6.4 | 2.5 | 2.8 | 0.7 |

100 zones, no `MATTHREADS` (load average 25–38):

| run | threads | placement | s per step | solve | chemistry | H₂ lines | GC |
|---|---|---|---|---|---|---|---|
| `t_th1p` | 1 | pinned | 0.198 | 6.6 | 5.6 | 7.0 | 1.1 |
| `t_th2p` | 2 | pinned, alone on the node | 0.144 | 6.3 | 3.5 | 4.1 | 0.8 |
| `t_thB21`–`4` | 2 | four runs on one node | 0.149–0.152 | 6.6–6.7 | 3.6–3.7 | 4.1 | 0.9 |
| `t_th4p` | 4 | pinned, alone on the node | 0.117 | 6.2 | 2.4 | 2.6 | 0.7 |
| `t_thA4a`, `b` | 4 | two runs on one node | 0.125, 0.126 | 6.5 | 2.6 | 2.8–2.9 | 0.8 |
| `t_th4u`, `t_th4u2` | 4 | not pinned | 0.118, 0.122 | 6.5–6.6 | 2.4–2.5 | 2.3–2.4 | 0.9–1.0 |
| `t_th8p` | 8 | pinned, one node | 0.103 | 6.2 | 1.9 | 1.7 | 0.7 |
| `t_th8u`, `t_th8u2` | 8 | not pinned | 0.130, 0.112 | 7.8, 6.8 | 2.6, 2.1 | 1.8, 1.7 | 2.1, 1.2 |

146 zones, from `frozen_c146_e`, steps 526500–529000 (`t_tc*` at load average 25–35):

| run | threads | placement | `MATTHREADS` | s per step | solve | chemistry | H₂ lines | GC |
|---|---|---|---|---|---|---|---|---|
| `t_tc8p` | 8 | pinned, one node | 8 | 0.0995 | 5.5 | 2.3 | 1.6 | 0.9 |
| `t_tc4p` | 4 | pinned, half a node | 4 | 0.122 | 6.4 | 2.9 | 2.4 | 0.8 |
| `t_tc8p0` | 8 | pinned, one node | off | 0.151 | 10.3 | 2.5 | 1.7 | 0.9 |
| `t_hd5` (earlier) | 8 | not pinned | off | 0.160 | 10.9 | 2.8 | 1.6 | 1.5 |

Production itself (146 zones, 8 threads, not pinned, no `MATTHREADS`, both switches off) runs at 0.148 s per step.

- `MATTHREADS` takes 29 % off a step at 8 threads and 22 % at 4 threads on 100 zones, and 34 % at 8 threads on
  146 zones. Without it the solve is 6.2–6.4 s per 100 steps on 100 zones whatever the thread count (10.3 s on 146).
- With `MATTHREADS`, 4 threads cost 26 % more wall time than 8 on 100 zones (23 % on 146), 2 threads 77 %.
- Pinning makes no measurable difference at 4 threads. At 8 threads the unpinned runs were 9–26 % slower than the pinned
  one (11 % with `MATTHREADS`), with more time in garbage collection and in the solve; three runs are not enough to give
  that a number.
- Two 4-thread runs on one node each lose 6 % against one alone. Steps per second per 8-core node on 100 zones with
  `MATTHREADS`: 13.1 with one 8-thread run, 19.5 with two 4-thread runs.
- A 100-zone run on 4 pinned threads (0.097–0.103 s per step) is 1.45–1.5 times faster than production as it runs now
  on 8; the 146-zone mesh on 4 pinned threads with `MATTHREADS=4` (0.122) is 1.2 times faster.
- `LINEINC=1` adds 0.5–0.65 s per 100 steps to the line transfer on 146 zones (`t_inc2`, `t_hd5` against `t_inc0`,
  `t_hd1b`), 3–4 % of a step; `HDCHEM=1` adds nothing measurable.

These numbers are for the present state (131 M☉, Δt 1.05e6 s, line transfer cheap). In the runs from the initial
conditions the line transfer is the largest item (the table in the section above), so they depend more on the thread
count; that was not measured here.

## A 64-zone baseline from the initial conditions, 128-zone continuations, and HD in the collapse (11 Nov)

From 11 Nov 11:00 the programme runs as four chains, each on 4 threads pinned to half a NUMA node with `MATTHREADS=4`.
All earlier chains are stopped; their records and checkpoints stay.

| run | what | zones | start | switches | cpus |
|---|---|---|---|---|---|
| `base64` | the new baseline, from the initial conditions | 64 | IC | `HDCHEM=1 LINEINC=1`; phase B with `ETASWTH=1e-3 ETASWTIE=0.5` | 8–11 |
| `nohd64` | twin of `base64` without HD | 64 | IC | `LINEINC=1` | 32–35 |
| `c128` | production toward 300 M☉ (continues `t_c146`) | 128 | 141.5 M☉ | `HDCHEM=1 LINEINC=1 ETASWTH=1e-3 ETASWTIE=0.5` | 16–19 |
| `f128` | `full6` continued | 128 | 14.3 M☉ | the same | 24–27 |

What was stopped, and why:

- `full3`, `full4` (10.53 and 10.88 M☉): `AVTENSION=1`, `MCLOUD=1e4`, phase B without `ETASWTH`/`ETASWTIE`. Their
  questions (when to switch the transport on, the collapse procedure) are answered in the sections above.
- `full6` (14.33 M☉, 380 zones, 16 threads, 0.095 M☉ per wall hour): the only run with the 2e4 M☉ cloud and
  `AVTENSION=0` from the start. Kept as `f128`.
- `t_c146` (141.52 M☉, 146 zones, 8 threads): the production line. Kept as `c128`.

### The launcher

`scripts/popiii_run/fullrun.sh` now reads, from the run directory:

- `launch.env`: `NZ` (shells of the initial model, default 200), `TH` (threads), `PIN` (a prefix command such as
  `numactl --physcpubind=8-11 --membind=1`), `REFINE` (`off` skips the shell refinement at the rebuild), `A1X`, `A2X`;
- `phaseA.env`, `stageA1.env`, `stageA2.env`, `stageA3.env`, `phaseB.env`: used instead of the defaults next to the
  script when present, so a run with its own mesh keeps its complete settings in its directory;
- `run.env` as before (it overrides the stage files).

`B1X="<settings>"` in the environment is passed to the first phase-B segment only: one-time switches of a restart.
The first line a chain writes to `procedure.log` gives shells, threads, pinning and the commit.

### `base64`

- Stage A1 with 64 shells to log n_c 14 (`RAD=0 NSTOP=1e14`), then the rebuild with continuum transport
  (`DTRADMAX=0`, as `full4`), no `REFINE`. The first attempt also ran stage A3 without `REZWIN` (64 Lagrangian
  zones to 0.5 M☉); that state could not be handed to phase B, see "The hand-over from 64 shells" below.
- Phase B: `ENVEUL=0,0.33,0.7,3000 ENVCORE=0.1,0.33,10`, i.e. Δln r = 0.33 from 0.1 to 3000 R☉ growing to 0.7 outside
  (`meshdesign.jl` logic: 64 zones to the cloud edge). About 14 zones lie inside 10 R☉. This is the coarse end of the
  planned resolution study and not a mesh for the stellar interior.
- A1 took 16 910 steps in 26 minutes (0.09 s per step). The rebuild carried the HD state (x_HD = 2.4e-5 in the
  molecular core, all deuterium in HD).

### The hand-over from 64 shells (11 Nov, not finished)

- **Unrefined stage A3 cannot be remapped.** With 64 Lagrangian shells the star at 0.5 M☉ is followed by one shell
  from 34.9 to 2884 R☉ holding 0.16 M☉. `envbuild!` interpolates ρ, T and E_r in log r between zone centres, so the
  Eulerian mesh receives a smooth ramp across that shell in place of a star edge and free fall. `base64` failed at the
  remap; `nohd64` passed it and then took up about 2.5 M☉, which is an artefact (enclosed mass at fixed radii before and
  after, `mprof.jl`). Both phase-B segments were discarded.
- **Stage A3 is therefore refined ahead of the star also for 64 shells:** `REZWIN=1.3,4,0.03 REZNMAX=340` as in the
  200-shell runs, with `REZQMAX=3.5` (the default 2.5 forbids the first splits of shells this coarse) and
  `DTSTOP=1e-7` in the run's `stageA3.env`. The lower `DTSTOP` is needed because each Newton divergence right after a
  bulk split cuts Δt by a factor 20–30; three or four in a row took Δt below the phase-A value 1e-2 s at step 19 022
  and ended the second attempt, although the step then converges.
- Tests of the window fraction from the A2 checkpoint (`t_ra1`–`t_ra3`): with `frac = 0.08` the stage reaches
  0.5174 M☉ at step 19 125 on 172 zones with 620 divergences; `frac = 0.03` gives about 300 zones and far fewer
  divergences and is what the two chains use.
- **A zone-flag error in `envbuild!`** stopped the first remap of a refined state onto 64 zones (`t_rb1`:
  "rho below the table" in the EOS fill). `PS.build!` filled the EOS while the per-zone table/network flag and
  composition were still on the old zone indices, so outer zones of the new mesh carried the flag of zones inside the
  old star whenever the source has more zones in the star than the target has in all. The flags are now set on the
  new indices before the build, as the `ENVMERGE` rebuild already did. Hand-overs that worked before are unchanged.
- With the fix the remap of the 172-zone test state onto 64 zones completes (`t_rb2`). Enclosed mass [M☉]:

  | radius [R☉] | 100 | 300 | 1e3 | 3e3 | 1e4 | 3e4 |
  |---|---|---|---|---|---|---|
  | source, step 19 120 | 0.477 | 0.485 | 0.574 | 0.865 | 1.962 | 5.080 |
  | after the remap, step 19 140 | 0.508 | 0.544 | 0.621 | 0.909 | 1.978 | 5.074 |
  | step 19 300, 6.9 yr later | 0.911 | 1.052 | 1.095 | 1.292 | 2.264 | 5.288 |

  The interpolating remap is not conservative: it adds 0.06 M☉ (12 % of the star) inside 300 R☉ and 3 % to the grid
  mass (20 000 → 20 596 M☉, mostly far out). The inflow afterwards, about 0.07 M☉/yr, is comparable to what the
  source did before the remap (0.1 M☉/yr), so it is not an obvious artefact. About half of the steps diverge
  (limiter at zones 24–26, the surface is zone 22; 151 divergences by step 19 400). This is not yet a clean hand-over.
- **A mass-conserving remap, `ENVCONS=1` (default off).** The zone densities of the new mesh are taken from the
  source's enclosed mass (a monotone cubic of ln m against ln r through the source's interfaces) in place of the
  log-interpolated density at the zone centres. T, E_r, the composition and the velocity are interpolated as before.
  Test `t_rb3`, the same hand-over as `t_rb2`: grid mass 20 000.0 → 20 000.0 M☉, M_shell0 0.4806 M☉; the largest
  change against the interpolated density is a factor 0.40 in the zone just outside the star (102–142 R☉).
  Enclosed mass [M☉] at step 19 140, the first snapshot after the remap:

  | radius [R☉] | 100 | 300 | 1e3 | 3e3 | 1e4 | 3e4 | 1e7 | 1e8 |
  |---|---|---|---|---|---|---|---|---|
  | interpolating (`t_rb2`) | 0.508 | 0.544 | 0.621 | 0.909 | 1.978 | 5.074 | 382.0 | 1233 |
  | `ENVCONS=1` (`t_rb3`) | 0.487 | 0.505 | 0.575 | 0.860 | 1.938 | 5.040 | 384.8 | 1240 |

  What it changes and what it does not, over 575 steps (41.5 yr; `t_rb2` ran 275 steps, 13.7 yr):
  - The accretion afterwards is the same. At 13.7 yr after the remap the star has 1.43 M☉ in `t_rb2` and 1.40 M☉ in
    `t_rb3` (interpolated between rows). `t_rb3` reaches 2.97 M☉ after 41.5 yr, 0.060 M☉/yr; `full6` (380 zones, no
    HD) went from 0.46 to 3.20 M☉ in 53.6 yr, 0.051 M☉/yr. So the inflow after the hand-over is the collapse, not the
    remap.
  - Fewer divergences early (63 by 11.4 yr and 104 by 22.7 yr, against 151 by 13.7 yr), and larger steps. They come
    back in bursts each time the photosphere moves out by one zone: 105 between steps 19 400 and 19 500, 96 between
    19 600 and 19 700.
  - The photosphere sits on mesh interfaces spaced by a factor 1.39 (Δln r = 0.33): R* goes 73, 102, 142, 197, 274,
    382, 531 R☉ and is 531 R☉ at 2.97 M☉, where `full6` has 147 R☉ at 3.2 M☉. This is the resolution of the
    64-zone mesh at the star's surface, not the remap, and it is the open problem of the 64-zone baseline.
- Cost of the refined stage A3 in `base64`: by step 26 000 (300 zones, 0.045 M☉) the ALI refresh interval had fallen
  to 1, driven by one zone at 85 R☉ (T 2150 K, 0.16 % of the line cooling) whose ALI/EP factor changes by 12–17 % per
  step. Each ALI call took 2.5 s, so the lines were 80 % of 1.6 s per step. `LINEFAST=1` was switched on in `base64`
  at step 26 000 for this reason (see the `LINEFAST` section; the twin from that checkpoint is bitwise equal). `nohd64` has it from
  step 18 500 (its row at step 18 700 equals the earlier trajectory).

### The mesh of the baseline: 64 fixed zones cannot carry it, and what a resolved shock costs (11 Nov, open)

Tom's reading of `base64` in stage A3: far too many zones at small enclosed mass, where nothing varies, and far too few
at the accretion shock. The record agrees. At step 50 000 (0.44 M☉, 340 zones) 31 zones lie inside 1e-3 M☉, 103 inside
1e-2 M☉ and 196 inside 0.1 M☉, and the shock is one zone from 55.8 to 347 R☉ (Δln r = 1.83). After the hand-over the
fixed 64-zone mesh has Δln r = 0.33 where the stellar surface is, and R★ sits on mesh interfaces: `base64` stayed at
200.952 R☉ from 0.75 to 1.70 M☉ and then hopped inward to 144.371 R☉; `nohd64c` is at 737.996 R☉ from 16.8 to
18.4 M☉. Both chains run on with tens of thousands of divergences (52 000 and, in the two segments, 85 000, by the evening of
11 Nov). They are not a baseline.

**An early hand-over removes the over-resolved core.** Stage A3 refines in mass and is where the zones at small enclosed
mass come from. Handing over to the radial mesh at 0.045 M☉ (HD; `t_eh1`, from `base64`'s checkpoint at step 26 000)
gives M★(t) within 1–2 % of `base64`'s own A3 + B and takes about 17 times less wall time than the rest of stage A3.
All tests below hand over early: at 0.045 M☉ with HD, at 0.164 M☉ without (from `nohd64`'s checkpoint at step 18 500),
with the conservative remap (`ENVCONS=1`) in the tests without HD.

**The reference is `full6`** (no HD; 380 fixed zones, Δln r = 0.03 from about 10 to 1000 R☉; 0.056 M☉/yr). Its shock
is one 3 % zone and R★ moves in steps of 3 %. It flashes and diverges from 0.5 to 2.05 M☉ (3404 divergences, R★
between 130 and 199 R☉), and then runs from 2.05 to 10.4 M☉ with no divergence at all, Δt 3–6.5e6 s, 200 steps from
4.59 to 5.76 M☉. R★ there: 147.2 R☉ at 3.20 M☉, 156.3 at 3.93, 166.0 at 4.59, 171.1 at 5.20, 176.3 at 5.76, 181.7 at
6.29, 187.3 at 6.78–7.25, 193.0 at 7.7–8.1, 205.0 at 9.7–10.5. It is a reference only for runs without HD: the HD
runs accrete at 0.008–0.010 M☉/yr and the radius follows the rate. There is no fine-mesh HD reference.

**A moving window that resolves the shock** (`GRIDWIN=1,10000,1e6,0.05,10,1.5,0.3,0.3,40` with `UPFREEZE=4`; the mesh
follows the temperature gradient, the shock is spread over about 8 zones of Δln r = 0.006):

| test | zones | physics | result |
|---|---|---|---|
| `t_gw64e` | 64 + window | no HD | R★ 481 R☉ at 3.12 M☉ (`full6`: 147). Too coarse outside the window. |
| `t_eh5` | 124 + window | no HD | 5–13 % above `full6` in the quiet phases, a swelling episode at 5.2–5.8 M☉ (to 414 R☉) that `full6` does not have; 230.7 R☉ at 12.0 M☉. |
| `t_eh7` | 124 + window, gmax 10 | no HD | clearly worse than `t_eh5`; stopped. |
| `t_eh6` | 206 + window | no HD | within 1–2 % of `full6` from 3.8 to 6.1 M☉ (157.8 R☉ at 3.73, 164.6 at 4.45, 171.6 at 5.07, 176.2 at 5.57, 180.3 at 6.11). Then it swells: 219 R☉ at 6.28 M☉, 287 at 6.64, 335 at 6.90, where `full6` has 182–187. |
| `t_st1`, `t_st4` | as `t_eh6`, `STOL=0.15` from 3.73 M☉ | no HD | the same swelling at the same mass: 184 R☉ at 6.06 M☉, 254 at 6.34, 292 at 6.52. |
| `t_stc1` | `t_eh6` from 6.89 M☉ with `STOLCHEM=1` | no HD | leaves the swollen state; 204.9 R☉ at 8.41 M☉ (`full6`: 193.0–198.9). See 3. below. |
| `t_eh3` | 124 + window | HD | 57.9 R☉ at 1.19 M☉, 61.3 at 2.00. |
| `t_eh4`, `t_eh8` | 206 + window | HD | 52.5 R☉ at 0.96 M☉, 54.2 at 1.34 M☉; then stuck (3. below). |
| `t_stc2` | `t_eh8` from 1.34 M☉ with `STOLCHEM=1` | HD | accretes again; 54.1 R☉ at 1.356 M☉. |

So 206 zones with the window agree with the fixed 380-zone mesh over 3.8–6.1 M☉, to within the 3 % by which the
reference is quantised. Three things stand against using it as the baseline as it is.

1. *Steps.* From 4 to 5.7 M☉ `t_eh6` took 2100 steps with 742 divergences (mean Δt 4.6e5 s) where `full6` took about
   200 with none (mean Δt 4.6e6 s): ten times the steps, four to five times the wall time (4 against 16 threads).
   Each divergence cuts Δt by 8 and the step then grows back to 1.5–2e6 s at most.
2. *The cause of the divergences is not found.* In a typical failed step the largest Newton correction is in the
   radiation energy at the upstream edge of the shock; it falls for five iterations, jumps at the sixth, falls again and
   then grows by 1.8 per iteration. Tested and without effect: `STOL=0.15` (fewer divergences for 1200 steps, then the
   same rate), the time constant of the mesh (1e5 and 1e7 s instead of 1e6), and holding the van Leer limiter decisions
   from iteration 5 on (a trial switch, 1101 decisions held in 60 steps, the divergence count unchanged to the digit;
   removed again). A smaller gmax (10) is worse. `NEWTFLIP` lines for the upwind direction show the raw decision,
   not the held one, so they report flips after iteration `UPFREEZE` that the solver does not make.
3. *Beyond 6.1 M☉ (no HD) and 1.34 M☉ (HD) the step falls to the time the gas needs to cross one zone at the shock.*
   This is an artefact of where the step rule is evaluated. The first zone behind the dissociation front (5500 K,
   optical depth 0.05 across the zone, x(H⁺) ≈ 1e-3) is heated by the transport step at fixed composition and cooled
   back by the chemistry step, every step: `chemsub k 118 … T 7204 → 5541` in `t_eh6`, with Δt = 2.2e4 s set by the
   30 % rule on that intermediate temperature and equal to the crossing time of the zone (0.6 % of 348 R☉ at 72 km/s).
   `t_st4` sits at 1.1–1.7e4 s in the same way, `t_eh8` at 1.2e3 s for 4400 steps (0.65 % of 54 R☉ at 95 km/s is
   2.6e3 s; its rule is 15 %). In the quiet phase of `t_eh6` the same zones existed (5700–6550 K behind a 2470 K zone
   at step 26 000) and the step was 37 crossing times, so the layer by itself is not the trigger; the gas ahead of the
   shock is three times less dense in the swollen state (5e-11 against 1.6e-10 g/cm³).
   With the rule applied after the chemistry step (`STOLCHEM=1`) both stalls end:
   - `t_stc2`, from `t_eh8` at step 43 500 (HD, 1.340 M☉): Δt goes from 1.2e3 s to 2–6e5 s at once and the star
     accretes again (1.356 M☉, 54.1 R☉ after 400 steps), with about 70 divergences per 100 steps.
   - `t_stc1`, from `t_eh6` at step 28 500 (6.89 M☉, 335 R☉): for 1500 steps the solver diverges on every second step
     (Δt 4e3–2e4 s) while R★ falls to 215 R☉ and then moves between 215 and 280 R☉; at step 30 000 the step returns
     to 0.3–2.5e6 s and R★ settles: 253 R☉ at 7.02 M☉, 221 at 7.35, 213 at 7.63, 207.2 at 8.00, 204.9 at 8.41, where
     `full6` has 193.0 at 7.7–8.1 and 198.9 at 8.55, so 3–7 % above the reference and closing. It takes 100 steps and
     80 divergences per 0.065 M☉ there, against 0.4–0.5 M☉ per 100 steps in `full6`. `t_stc3` continues it.

   So the radius between 6.1 and 7 M☉ on this mesh depends on the step rule (335 R☉ and rising with the rule on the
   intermediate state, 215–280 R☉ with it on the full step, 182–187 R☉ in `full6`), and the split between transport
   and chemistry is not converged in the relaxation layer behind a resolved shock. Whether a swelling episode at
   6.1 M☉ is real is not established. `t_eh5` has a similar one at 5.2–5.8 M☉ and left it without help; `full6` has
   them only below 2.05 M☉ and at 10.47 M☉.

**A fixed mesh that is fine only where the surface lives.** `ENVEUL` takes two more entries (3fc1c41): the factor per
zone by which the spacing grows above rb and the factor by which it shrinks above the core boundary.
`ENVCORE=0.1,0.3,3 ENVEUL=0,0.03,0.7,1000,1.15,0.7` gives 224 zones: 12 inside 3 R☉, 176 at Δln r = 0.03 from 3 to
1000 R☉ (the range R★ has covered in all runs, 8 to 740 R☉), 36 outside. It is `full6`'s mesh where the shock and
the surface are, with the coarse core and envelope of the 64-zone runs. `t_eh9` (no HD, to be compared with `full6`)
and `t_eh10` (HD) started on it at 19:24 and 19:47 on 11 Nov. First rows of `t_eh9`: 103.0 R☉ at 0.700 M☉ and 138.4 at
0.959 (`full6`: 102.5 at 0.70, 142.8 at 0.98). The comparison that counts is the quiet phase above 3 M☉.

*Without HD the fixed 224-zone mesh reproduces `full6`* (`t_eh9`, 21:00 on 11 Nov). It flashes up to 2.0 M☉ as `full6`
does (5456 failed solves from the hand-over at 0.164 M☉; `full6` 3404 from 0.5 to 2.05 M☉), is quiet from 2.26 M☉ on, and has had no failed solve since.
Both runs take 600 steps from 3.0 to 6.5 M☉ with Δt = 3–6 × 10⁶ s. R★ sits on 3 % interfaces in both runs (the two
ladders are offset), so the `full6` column is interpolated in mass between its rows:

| M★ [M☉] | R★ `t_eh9` [R☉] | R★ `full6` [R☉] | ratio |
|---|---|---|---|
| 2.26 | 138.4 | 134.4 | 1.03 |
| 3.02 | 151.2 | 144.0 | 1.05 |
| 3.74 | 160.4 | 153.9 | 1.04 |
| 4.43 | 170.2 | 163.6 | 1.04 |
| 5.06 | 175.3 | 169.9 | 1.03 |
| 5.66 | 180.6 | 175.4 | 1.03 |
| 6.22 | 186.0 | 181.1 | 1.03 |
| 6.75 | 191.6 | 186.9 | 1.03 |
| 7.25 | 191.6 | 187.3 | 1.02 |

The fixed band mesh is 2–5 % larger, about one interface, and the difference shrinks with mass. The runs differ in more
than the mesh: hand-over at 0.164 instead of 0.5 M☉, the conservative remap, 12 instead of 43 zones inside 3 R☉ and
a coarse envelope. Which of these carries the 3 % is not separated.

*With HD the early radius is not settled between meshes* (`t_eh10`, fixed 224, against `t_eh4`, 206 with the window,
and `t_eh3`, 124 with the window; all from the same hand-over checkpoint, M★(t) the same within a few per cent):

| M★ [M☉] | `t_eh10` | `t_eh4` | `t_eh3` |
|---|---|---|---|
| 0.30 | 56.2 | 52.6 | 63.0 |
| 0.40 | 75.6 | 52.8–54.4 | 64.4 |
| 0.45 | 75.6–82.6 | 49.5–51.4 | 56.1 |
| 0.50 | 61.4 | 49.0–49.3 | 65.8 |

`t_eh10` flashes through this range (Δt 1–10⁴ s, an excursion to 82.6 R☉ at 0.40–0.47 M☉ and back to 60 R☉ at
0.52 M☉); `t_eh4` crosses it smoothly at Δt 10⁵–10⁶ s. This is the phase in which `full6` flashed as well and in which
the failed-step fraction does not discriminate between meshes; there is no fine-mesh HD reference. The numbers to
compare later are the window runs past the flashing: 49.6–52.7 R☉ from 0.6 to 1.0 M☉ and 55.7 R☉ at 1.63 M☉
(`t_eh4`, `t_eh8`, `t_stc2`), 57.6–59.6 R☉ from 1.2 to 1.5 M☉ and 61.3 R☉ at 2.0 M☉ (`t_eh3`).

*Past the flashing the HD runs agree* (22:00 on 11 Nov). `t_eh10` flashed until 0.63 M☉ (step 33 600, 4553 failed
solves) and has been quiet since: 70 failed solves in the next 2400 steps, Δt 10⁶–10⁷ s.

| M★ [M☉] | fixed 224 `t_eh10` | window 206 `t_eh4`/`t_eh8`/`t_stc2` | window 124 `t_eh3` |
|---|---|---|---|
| 0.6 | 57.9 | 49.6–49.9 | |
| 0.8 | 54.5 | 51.4–51.5 | 84.7 |
| 1.0 | 56.2 | 52.6–52.7 | 60.9–62.5 |
| 1.3 | 57.9 | 54.0–54.3 | 57.6–58.8 |
| 1.6 | 59.6 | 54.8–79.4 | 59.9–60.3 |
| 2.0 | 61.4 | | 61.2–61.3 |
| 2.4 | 63.3 | | |

From 1.3 M☉ the fixed mesh and the 124-zone window run agree within 3 %, one interface; the 206-zone window run is
7 % smaller at 1.0–1.3 M☉. The 20–50 % of the first table belongs to the flashing excursion, not to the quiet branch.
The baseline `b224` is a second realization of the same flashing (different step sequence from 26 500 on) and reads the
same R★(M★) up to 0.3 M☉ (17.2–18.2, 37.1, 56.2 R☉ at 0.1, 0.2, 0.3 M☉).

`t_eh9` meanwhile reached `full6`'s slow phase: 197.3 R☉ at 8.0 M☉ (`full6` 193.0), 203.2 at 9.0 (198.9), 209.3 from
10.0 to 10.55 M☉ (205.0, then 211.2 from 10.57), with Δt falling to 1.8 × 10⁵ s at 10.5 M☉ as in `full6` and six
failed solves since 2.26 M☉.

*The window mesh with `STOLCHEM=1` rejoins `full6`.* `t_stc3` continues `t_stc1` from 8.41 M☉: 203.4 R☉ at 9.0 M☉
(`full6` 198.9), 204.0–205.9 at 10.0 (205.0), 207.5–209.4 at 11.0 (205.0–211.2), 211.3 at 11.52 (211.2). From 10 M☉ on
the two agree within 1 %. Its steps are still short (Δt 4 × 10⁴–4 × 10⁵ s) and it fails about 35 solves per 100
steps; `full6` itself is slow there (6 × 10⁴–2 × 10⁵ s past 10.47 M☉). `t_stc2`, the HD window run with the same
rule, ran its 6000 steps from 1.34 to 1.63 M☉ at Δt 2–6 × 10⁵ s instead of 1.2 × 10³ s, with 3657 failed solves.

**Decided on 11 Nov (Tom).**

- `base64` and `nohd64c` are stopped (20:01; 1.78 M☉ at step 128 600 and 18.7 M☉ at step 138 200). Their stage-A3
  checkpoints `z19base64_A3_3f` (0.045 M☉) and `z19nohd64_A3_3f` (0.164 M☉) are the starting points of the new baseline.
- The early hand-over is the procedure: `stageA3.env` ends the stage at `MSTOP=0.05`. Tested from 0.045 M☉ with HD on
  three meshes (`t_eh1`, `t_eh4`, `t_eh10`); without HD only from 0.164 M☉.
- New runs take `LINEFAST=1 LINEXC2=16 EPTAB=1` (in `phaseA.env` and `phaseB.env`) and the conservative remap
  `ENVCONS=1` (`phaseB.env`). `c128` and `f128` keep the settings they were started with: the previous `phaseB.env` is
  copied into their run directories.
- The mesh of the baseline waited for `t_eh9` above 3 M☉ (agreed plan: if the run without HD reproduces `full6`,
  the HD baseline starts on the fixed 224-zone mesh and the window mesh is kept for the resolution studies). With the
  table above the chain `b224` was started at 21:00 on 11 Nov from `z19base64_A3_3f`: HD, CMB, 4 threads on cpus 8–11,
  mesh and settings of `t_eh10` plus `LINEXC2=16 EPTAB=1`, all in its own `phaseB.env`. The default `phaseB.env` still
  has the 380-zone mesh `ENVEUL=0,0.03,0.1,1000`. `t_eh10` keeps running beside it as a second realization of the
  flashing phase (it differs from `b224` by `LINEXC2`, a 10⁻⁹ change in the line luminosity).

**Where this stands.** The baseline `b224` runs on the fixed 224-zone mesh, which inherits `full6`'s one-zone shock
and 3 % steps in R★ and reproduces `full6` within 2–5 % without HD. With HD the flashing below 0.63 M☉ makes excursions
of 20–50 % in R★ that the 206-zone window run does not have; past it the three HD meshes agree within 3–7 %. The window mesh resolves the shock, needs
`STOLCHEM=1` to run at all past the first swollen state, costs about ten times the steps, and its divergences are not
understood; with that rule it agrees with `full6` within 1 % from 10 to 11.5 M☉. `STOLCHEM` and `GRIDWIN` are off in
all production chains.

### A mesh that follows radius and mass: `GRIDMAP` (12 Nov, in test)

**The fixed 224-zone mesh after the first hours of `b224`.**

- *The stall of the step rule is there too.* `t_eh10` stopped advancing at 2.40 M☉ (63.3 R☉) and `b224` was held at
  1.337–1.348 M☉ from step 34 800 to 38 200 (22:15 to 23:53 on 11 Nov) with Δt 2.4 × 10⁴ s and no failed solves: the
  same artefact as on the window mesh (zones behind the dissociation front, `STOL` on the intermediate temperature).
  `t_stc4` continues `t_eh10`'s checkpoint with `STOLCHEM=1`: 2.40 → 5.06 M☉ in 4800 steps at Δt 2–5 × 10⁶ s with 117
  failed solves, R★ 67.1 R☉ at 3.2 M☉, 65.2 at 3.5–3.8, 67.1 at 4.5–5.1. `t_eh9` (no HD) stalled the same way at
  10.9–11.05 M☉ (215.6 R☉) and was stopped; `t_stc3` (window mesh, `STOLCHEM=1`) ended at 12.53 M☉, 220.7 R☉.
- *`b224` left its stall by itself, through a second swelling*: 95.8 R☉ at 1.357 M☉ (step 38 800), back to 59.6 R☉ at
  1.53 M☉, then 61.4 at 1.8–2.0, 63.3 at 2.1–2.8 and 65.2 at 3.05 M☉ (00:53, 12 Nov), on `t_eh10`'s quiet branch and past
  the mass where `t_eh10` stalled. `STOLCHEM` is still off in `b224`.
- *Where the zones are* (`meshcensus.jl`). `b224` has 1–4 zones per half decade of radius beyond 3000 R☉, with
  Δ ln m of 0.23–0.70 per zone: the infalling envelope that feeds the star for the next 10³ yr is described by a
  handful of zones. `c128` carries a 23.5 M☉ shell near 5600 R☉ with Δ ln ρ = 7.7 across one zone; `f128` has
  Δ ln r = 0.24 inside 10 R☉. Each fixed mesh is good in the range it was built for and poor elsewhere, and that
  range moves: R★ changes by a factor 20 and M★ by 10⁴ over a run.

**The rule** (`julia/GridMap.jl`, no dependence on the rest of the code). A fixed number of zones N; every zone covers
the same arc length in three coordinates,

    X  = ln(r + r_c) / h_r                 radius (r_c: the radius that encloses f_c M★, so the uniform centre is not resolved in ln r)
    Q1 = ln(m + f_c M★) / h_1              enclosed mass: stellar interior, dense shells, the envelope once it outweighs the star
    Q2 = asinh((m − M★)/(ε M★)) / h_2      mass counted from the surface in both directions: outer layers of the star, the infall

so radius sets the spacing where the enclosed mass does not change (the free-fall region between the shock and the
massive envelope) and mass sets it where it does, and the hand-over between the two follows M★(t). M★ is the mass
inside the point where the infall has dropped to 0.2 of its fastest value, interpolated between interfaces so that
it is continuous in time. Two limits: no zone narrower than `dmin` in ln r (0.02), neighbouring zones differ by at
most the factor `g` (1.2). The three widths (defaults h_r 0.3, h_1 0.3, h_2 0.4, ε 10⁻⁵, f_c 10⁻³) are scaled by one
factor λ, found by bisection, so that the zones add up to N. One call costs 23 ms for N = 200.

**How the mesh moves.** `GRIDMAP=N[,hr,h1,h2,eps,fc,dmin,g,tau,fmove,nrelax,ctol,ftar]` (with `ENVEUL`, default off).
At the hand-over the phase-B mesh is built from the rule instead of the `ENVCORE`/`ENVEUL` spacing (same remap,
`ENVCONS=1`). After every step, once the chemistry and `STOLCHEM` have fixed the next Δt, the driver computes the
target mesh of the present state and sets one velocity per interface (`TitanResidFast.GRIDVEL`): the mesh closes
1/`nrelax` (1/3) of its distance from the target in the next step, at most `fmove` (0.4) zone widths. The grid row of
the residual is r_new − r_old − w Δt, so the velocities are data for the Newton solve and add nothing to the
Jacobian; everything else is TITAN's moving-grid form as for `GRIDWIN` (u_rel = u − ṙ, fluxes through the moving
faces, lab-frame gravity work, abundances advected with the mass that crosses each interface). Δt is limited so that
the target itself moves at most `ftar` (0.2) zone widths per step, judged from its move in the step just taken.
The first version without this limit (`t_mm1`: at most 0.2 zone widths per step and not faster than in 10⁶ s) fell 35
zone widths behind its target once the steps reached 3 × 10⁶ s.

**First test** (`t_mm2`: `GRIDMAP=199 STOLCHEM=1`, otherwise `b224`'s settings, from the same checkpoint
`z19base64_A3_3f`; 200 zones).

- With the switch off the patched sources are bit-identical to the previous ones over 400 steps from a 224-zone
  checkpoint (`t_bi3`/`t_bi4`: all fields of all 20 snapshots equal).
- The total mass on the grid is constant to 1 part in 10¹² over 7240 steps (19 994.670 688 69 → 19 994.670 688 67 M☉; the fixed mesh of
  `t_eh10` does the same).
- The mesh stays within 0.01–0.6 zone widths of its target. At 0.906 M☉ (step 33 260): 34 zones between 31.6 and
  100 R☉ with eight zones of Δ ln r = 0.020 across the accretion shock (the infall goes from −71 to −1 km/s between
  68 and 57 R☉; on the fixed mesh this is one 3 % zone), 8–10 zones per half decade everywhere from 0.3 to 3 × 10⁹ R☉
  (6–7 between 3 × 10⁵ and 10⁷ R☉), and Δ ln m ≤ 0.18 per zone outside 1 R☉.
- Through the flashing phase it follows `b224`: the same swelling to 100–107 R☉ at 0.50–0.59 M☉ (`b224` 99–111) and the
  collapse at 0.60–0.62 M☉; times at fixed mass agree within 0.5 yr up to 0.5 M☉. Steps to 0.6 M☉: 7000 with 3468 failed
  solves (`b224` 7100 with 3628, `t_eh10` 7400 with 4517): no saving there.
- Quiet branch: 58.1 R☉ at 0.815 M☉ and 55.9 at 0.906 (`t_eh10` 54.5 at 0.8, 56.2 at 0.9–1.0).
- *Open: a swelling at 0.91 M☉ that the fixed meshes do not show there.* After steps of 10⁷ s a step failed down to
  Δt = 4 s (step 33 200); a layer of 2.6–4.8 × 10⁴ K appeared at the shock, which moved from 56 to 126 R☉ within a year
  and came back with three smaller swings (up to 5 km/s outflow at the surface); about 430 failed solves until it settled at
  65 R☉, 0.926 M☉ (step 34 360). `b224` did the same at 1.35 M☉ and not here. Whether the trigger is the mesh motion
  at large Δt is not known.
- *Open: the infall is faster than on the fixed mesh.* Mass inside 10⁴ R☉, both runs starting from 0.839 M☉ at
  t₀ + 30.1 yr (t₀ = 10 178 900 yr): 1.002 M☉ at 57.4 yr and 1.171 at 85.4 yr in `t_mm2` against 0.991 and 1.153 in `b224`
  (interpolated; `t_eh10` as `b224`), i.e. 6 % more mass through 10⁴ R☉ in 55 yr; inside 1000 R☉ the difference is
  1.3 % at 85 yr. Between 3 × 10³ and 3 × 10⁴ R☉ `b224` has zones of
  Δ ln r 0.21–0.49 (0.69 beyond 10⁵ R☉) and `t_mm2` 0.13–0.14. `t_fe1` (fixed mesh as `b224` but
  Δ ln r ≤ 0.13 outside 1000 R☉, 310 zones, `STOLCHEM=1`, started 00:48 on 12 Nov) will tell whether the difference is
  the resolution of the envelope or the moving mesh.
- Not tested: a restart of a fixed-mesh checkpoint with `GRIDMAP`; `mat.Λext`, `Tref` and the per-zone equilibrium
  flags are not advected with the mesh.

`GRIDMAP` is off in all production chains; adopting it for the baseline is Tom's decision.

### `c128` and `f128`: a restart with a one-time merge

Both start from frozen copies of the last checkpoints (`frozen_c146_f/m2c_z19c146f_auto0`, 146 zones;
`frozen_full6_380/m2c_z19full6B1_auto1`, 380 zones) with `ENVMERGE` in `B1X`:

- `c128`: `ENVMERGE=11,30,2;30,1500,2` (146 → 138 → 128 zones). No `MCLOUD` (the restart check needs the settings of
  the original run), `ENVMS0=0.4737`.
- `f128`: `ENVMERGE=0.06,1e12,3` (380 → 128 zones, every three zones merged), `MCLOUD=2e4`, `ENVMS0=5.014e-01`.

Mass, radiation energy and gas energy are conserved by the merge to the printed digits. HD starts from its equilibrium
(`hd_pad!`). To undo HD and the CMB term: copy the frozen checkpoint into a new run's `ck/` and start
`fullrun.sh <run> B:<checkpoint>` with a `run.env` without `HDCHEM`/`LINEINC`.

First hour:

- `c128`: 0.113 s per step (146 zones on 8 threads: 0.148), no divergences. R* went from 10.47 to 9.66 R☉ at the
  merge: one interface of the mesh.
- `f128`: 1.2–2 M☉ per wall hour against 0.095 for `full6`. R* 231 → 245.5 R☉ at the merge, 268.8 R☉ from step
  83 000. From step 81 600 the run diverges about 50 times per 100 steps, in bursts (none between steps 83 100 and 83 300), with
  Δt between 7e2 and 9e5 s. The limiter
  is the temperature of the zone just inside the photosphere. The mesh has Δr/r = 0.086 there (0.03 in `full6`): the
  shock is one zone (31 900 K, 8.8e-8 g cm⁻³ → 5600 K, 7.8e-10), four zones at 5100–5600 K follow, then 2750 K, and
  the photosphere moves between two of them. `c128` has the same Δr/r at its surface but its photosphere is far out
  in the infall (1e4 R☉). Not fixed; the run advances. A restart with fewer merged zones between 150 and 700 R☉ is
  the remedy to try if it stalls.

### HD puts the collapse on a colder branch

At the same central density (log n_c = 7.80):

| run | shells | HD | `LINEINC` | t [yr] | T_c [K] | x_H2 | x_e |
|---|---|---|---|---|---|---|---|
| `base64` | 64 | on | on | 1.0005e7 | 220.5 | 4.04e-3 | 5.75e-10 |
| `nohd64` | 64 | off | on | 1.7458e7 | 675.5 | 9.06e-4 | 6.19e-9 |
| `b64n` (control, stopped) | 64 | off | off | 1.7454e7 | 675.5 | 9.06e-4 | 6.19e-9 |
| `full6` | 200 | off | off | 2.2691e7 | 627.3 | 9.20e-4 | 4.23e-9 |

- `base64` and `b64n` differ only in `HDCHEM=1 LINEINC=1`; `nohd64` shows that `LINEINC` alone changes nothing here.
  The cold core is the HD.
- In `base64` at that time Λ_HD/Λ_H₂ is 0.2 in the core, 1.2 at n = 8e6, 3.4 at 7e5, 4.7 at 1.3e5 (T = 159 K,
  x_HD = 3.7e-7), 2.1 at 1.2e4, 0.8 at 2.5e3 cm⁻³: HD is the main coolant between about 3e3 and 1e7 cm⁻³ and holds the
  gas at 159–185 K there. Without HD: 333 K at 3.7e4 and 495 K at 1e6 cm⁻³ (`b64n`). The mass above 3e7 cm⁻³ is 4 M☉
  with HD and 23 M☉ without. (`hdcore.jl <run A> <run B> <log n_c>` in the session scratch directory prints these.)
- Later the tracks approach each other: T_c = 482 K at log n_c 10.0 and 626 K at 10.64 with HD, 862 K at 9.98 and
  964 K at 10.57 without.
- Resolution alone (`b64n` against `full6`, both without HD): T_c +8 %, time of collapse −23 %.
- The literature has this branch: Ripamonti (2007, MNRAS 376, 709) finds HD cooling to 50–100 K in objects below
  ~3e5 M☉ at z_vir = 20 where H₂ alone stays above 200 K; McGreer & Bryan (2008, ApJ 685, 8) find HD equal to or above
  H₂ cooling in the lowest-mass halos and between 1e2 and 1e6 cm⁻³ in pre-ionized ones, with lower stellar masses.
  Whether a cloud enters it depends on how long the gas stays near 1e3–1e5 cm⁻³. Our 1D cloud takes 1e7 yr to
  collapse, which gives HD that time; the result is therefore tied to the initial model and is not a general
  statement about halos of this mass.

Consequences:

- `base64` is a different star from the `c128`/`f128` lineage, which collapsed without HD and received it at 141.5 and
  14.3 M☉. Its core is colder, so a lower accretion rate is expected at first.
- A resolution comparison needs the same physics on both sides: `nohd64` against `full6`/`f128` (no HD in the
  collapse), `base64` against future finer runs with HD.

## A faster formal solution for the H₂ lines: `LINEFAST=1` (11 Nov)

Where the line time goes (4 threads, the `NLTE=hyb` production setting):

- Every step makes one escape-probability call. Every `ALIINT` steps (20 at the cap, where the runs sit almost always)
  a full ALI call follows: two iterations, each a formal solution of every thick line (τ_c > 0.02) on
  `clamp(2 v_span / (0.5 min v_th) + 1, 30, 200)` frequency bins.
- Before this change the lines took 0.10–0.13 s of the 0.26–0.49 s per step in `base64`, `nohd64` and `f128`
  (24–31 %; the chemistry takes 0.17–0.24 s), and 0.023 of 0.113 s in `c128`, which has no thick H₂ lines.
- In the accretion phase almost all of it is the formal solution: 7.1 CPU-s of a 1.85 s ALI call on `f128` step 87 000
  (83 reduced zones, 55 thick lines, 198 bins, |u| up to 114 km/s). Every bin of every segment of every ray was
  integrated, although a segment only has opacity within a few thermal widths of its own velocity.

`LINEFAST=1` (default off, commit a614b3e) keeps rays, bins and quadratures and reorganises the solution
(`transfer_line_fast!` in `julia/H2Lines.jl`, used with `LINESRC=zone`):

1. Per segment only the bins inside the resonance window |v − μu| ≤ v_th √(38 + ln max(τ₀, 1)) are integrated; outside
   it the optical depth of the segment is below e⁻³⁸ and the intensity is unchanged to round-off. Segments without
   opacity are skipped.
2. The loop order is ray → segment → bin, with a branch-free exponential (within 1 ulp of `Base.exp`, 1.3 against
   3.8 ns per element) so that the bin loop vectorises.
3. Bin sums are formed first and the angle quadratures once per line; the work arrays come from a per-thread pool.

ALI call, wall seconds on 4 threads, reference → `LINEFAST=1`:

| snapshot | regime | thick lines, bins | reference | fast |
|---|---|---|---|---|
| `f128` 87 000 | accretion | 55, 198 | 1.85 | 0.14 |
| `base64` 12 000 | early collapse | 3, 40 | 0.071 | 0.076 |
| `base64` 15 000 | collapse | 23, 32 | 0.123 | 0.081 |
| `base64` 17 500 | late collapse | 217, 31 | 0.384 | 0.135 |
| `base64` 18 400 | late collapse | 216, 33 | 0.347 | 0.122 |
| `base64` 22 000 | stage A3, 250 zones | 217, 44 | 2.08 | 0.55 |

Agreement with the reference on these cases and on two without the incident background (`LINEINC` off):

- populations: every one above 1.4e-7 is bitwise equal; the rest agree to 1.2e-14;
- −∇·F per zone to 5e-14 of its maximum, L(r) and the luminosity per line to 2.5e-15;
- the ALI iteration counts and convergence histories are the same.

With `LINEFAST` unset the original code runs and its output is bitwise what it was before the change.

Twin runs, `LINEFAST=0` against `1`, 4 pinned threads, same steps:

| | steps | wall s/step | lines s/step | chemistry s/step | divergences |
|---|---|---|---|---|---|
| stage A3 from `base64`'s A2 checkpoint (`t_lfa0`, `t_lfa1`) | 18 600–19 000 | 0.294 → 0.259 | 0.072 → 0.045 | 0.170, 0.162 | 5, 5 |
| `f128` from `frozen_f128_lf` (`t_lfb0`, `t_lfb1`) | 91 600–92 100 | 0.405 → 0.345 | 0.102 → 0.042 | 0.202, 0.204 | 99, 99 |

- All snapshots of each pair (25 and 30) are bitwise equal in t, Δt, r, T, ρ, u, e, the composition and L_lines. The
  reason: with `LINECOLL=1` the gas term is the collisional exchange at the ALI populations, and the populations that
  carry it are bitwise the same.
- The mean ALI call went from 0.66 to 0.26 s (A3) and from 1.32 to 0.16 s (`f128`). The step gains 12 and 15 %,
  because what remains of the line cost (0.04 s per step) is the escape-probability call of every step, and the
  chemistry is now the largest part.

A third twin where the lines dominate: stage A3 of `base64` from its step-26 000 checkpoint (298 zones, ALI
refresh every step; `t_lfc0`, `t_lfc1`), steps 26 100–26 200: wall 1.64 → 0.67 s per step, lines 1.29 → 0.32,
chemistry 0.18 in both, ALI call 2.5 → 0.5 s, no divergences. The ten common snapshots are bitwise equal in the
state; one entry of the line heating term differs by 6e-30 of its maximum in one snapshot. `base64` continues from
this checkpoint with `LINEFAST=1` in its `run.env`.

One error was found by the snapshot comparison and fixed before the twin runs: the reused work array kept the
outer-node intensity of the previous ray, which showed only as a 7e-8 difference in L(r) while the populations agreed
to 9e-15.

Not done: changes of the method (one ALI iteration per call, a longer refresh interval) and work on the
escape-probability call. The benchmark tools are in the session scratch directory (`lb/linebench.jl`, `lb/cmp2.jl`,
`lb/twin.jl`, `lftest.sh`).

## More on the line cost: refresh cap, resonance window, screening table and a GPU formal solution (11 Nov)

After `LINEFAST=1` the lines take 20 % of the step in `c128` (0.113 s per step), 10 % in `f128` (0.340 s), 13 % in
`nohd64c` (0.183 s) and 24 % in stage A3 of `base64` (0.288 s). The two regimes differ:

- **Phase B.** The formal solution is at most about 4 % of the line time, under 1 % of the step. What remains is the
  escape-probability call that every step makes. Nothing that speeds up the formal solution pays here.
- **Stage A3** (refined Lagrangian mesh, about 300 zones, more than 100 thick lines, refresh interval often 1). The
  ALI calls were 2036 s of the 8668 s wall time of `base64`'s last A3 segment. This is where the options below pay.

Four switches, all default off or default unchanged, none switched on in a running chain:

| switch | commit | what it does | measured |
|---|---|---|---|
| `ALICAP=<n>` | b1b766e | cap of the ALI refresh interval, separate from `NALI` (which the restart check pins) | `c128` and `f128` sit at the cap of 20 with drift 0.000, `nohd64c` 98 % of the time; a higher cap is safe there and saves under 1 % of the step |
| `LINEXC2=<x>` | 50c85cc | the resonance window is \|v − μu\| ≤ v_th √(x + ln max(τ₀, 1)); default 38 as before | `LINEXC2=16`: formal solution −20 % (`base64` step 50 000: 0.250 → 0.198 s); line luminosity changes by 2.9e-9 (`base64`) and 7.3e-10 (`f128`) |
| `EPTAB=1` | 1ba004f | screening factors of the escape-probability call from a table built once per call | bitwise equal; the call goes 0.037 → 0.022 s (`base64`), 0.032 → 0.022 (`f128`), 0.017 → 0.011 (`nohd64c`), no gain in `c128` |
| `LINEGPU=cuda\|metal\|cpu` | a5eeef4 | formal solution of all thick lines of an ALI iteration in one batch on a GPU | see below |

### `LINEGPU`: the formal solution through KernelAbstractions

`julia/H2LinesKA.jl` is a port of `transfer_line_fast!` to KernelAbstractions.jl. All thick lines of one ALI iteration
are solved in three kernel launches: the sweeps along the rays (one work item per line, ray and chunk of frequency
bins), the sum over chunks, and the angle quadrature per line and zone. The flux quadrature and everything that
carries ν₀ dv/c or r is done on the host in Float64.

- The kernels are generic in the float type and hold no backend-specific code. What reaches the device is normalised
  on the host (sources by a power of two near the line's largest, lengths by a power of two near the outer radius),
  so Float32 does not underflow on Planck-like source values.
- In Float32 the device accumulates (S − I)ψ per ray segment and the host forms the intensity differences by a
  running sum in Float64, so no difference of nearly equal numbers is taken in single precision.
- `LINEGPU=cuda` or `metal` selects the backend; `LINEGPUFLOAT=32|64`, default 32. Metal has no Float64. `LINEGPU=cpu`
  runs the same kernels on the KernelAbstractions CPU backend; it is a check, five times slower than the native path.
- `popiii_m2c.jl` includes `h2lines_gpu.jl` only when `LINEGPU` is set, and the run then needs a Julia environment
  with KernelAbstractions and CUDA or Metal next to ForwardDiff and StaticArrays. The production environment is
  unchanged. It acts only with `LINEFAST=1 LINESRC=zone`.
- With `LINEGPU` unset, ALI and EP results on `f128` step 122 500 and `base64` step 50 000 are bitwise equal to the
  source before the switch.

Measured on the RTX A6000 through the switch, against the 4-thread CPU path:

| | CPU | CUDA Float32 | CUDA Float64 |
|---|---|---|---|
| ALI call, `base64` step 50 000 (127 thick lines) | 0.30 s | 0.060–0.075 s | 0.11–0.13 s |
| ALI call, `f128` step 122 500 (14 thick lines) | 0.074 s | 0.04–0.07 s | 0.078 s |
| net line cooling per zone, largest difference / peak | | 1.1e-8, 4.5e-9 | 8e-17 |
| populations above 1e-12, largest relative difference | | 7.3e-7, 2.4e-7 | 7e-16 |

Twin run in the driver, stage A3 from `base64`'s step-26 000 checkpoint, steps 26 200–27 500, 4 pinned threads
(`t_lg0` CPU, `t_lg1` `LINEGPU=cuda` Float32): 0.743 → 0.405 s per step, lines 0.46 → 0.13 s per step. Δt agrees to
three or four digits throughout, M* to all printed digits, and R* to 5e-4 at the last step (17.774 against
17.783 R☉). Device memory of the batch: 373 MB for 92 reduced zones and 127 lines.

Limits of what was tested:

- **Metal has not been run.** There is no Apple machine here. `julia/h2lines_ka_selftest.jl` checks a backend without
  any run data: a synthetic accreting envelope (opaque star, shock, free fall, sub-thermal sources, line-centre
  optical depths 0.6 to 3e4), Float32 on the device against Float64 on the CPU backend, with and without the incident
  field. On CUDA it passes with errors of 5e-7 to 1.6e-6 in J̄, 5e-7 in Λ* and 7e-8 of the peak in the net cooling.
  On a Mac: `julia -t 4 --project=<env with KernelAbstractions and Metal> julia/h2lines_ka_selftest.jl metal`.
- The workspace cache is process-global and not thread-safe; the driver calls it from one task.
- The twin covers 1300 steps of one stage. Phase B was not run with it, since there is nothing to gain there.

### A cooling-share floor in stage A3 (`ALIWMIN=1e-2`): not adopted

In stage A3 the refresh interval is driven to 1 by a single zone ahead of the shock that carries 0.16 % of the line
cooling. Raising the floor below which a zone is ignored in the refresh test from 1e-3 to 1e-2 of the total makes
the stage 1.75 times faster over a 0.7 yr twin, but the total line luminosity then differs by 10–15 %, in the two
zones just ahead of the shock. The floor stays at 1e-3.

### Where the chemistry time goes (12 Nov)

Profiled (`Profile`, 300 steps from `fz_b224`, step 84 500, 4.22 M☉, the stalled state before `STOLCHEM`): chemistry 70 % of the
wall time (hydro 12, chemistry 32–37, lines 4 s per 100 steps). Counters (`CHEMSTATS=1`, default off): 38 834 network steps,
55 736 attempts, 925 935 Newton iterations — 24 per step, against 3.2 in `c128` (182 M☉). About 17 000 attempts fail and each runs
to `maxit` = 60 before the step is halved; that is nearly all the iterations.

- Bit-identical (0c442fa, in production code): the reverse-rate power law below T_u is computed once (it was recomputed with two
  `lnKeq` on every call); the frequency-only terms of `nlte_means` are tabulated once. A memoised `lnKeq` was slower and dropped.
- `CHEMSTAG=8 CHEMGIVEUP=15` (57d6f45, default off): accept a round-off-stagnated attempt from iteration 8; abandon an attempt that is
  still far from converged at iteration 15. Newton iterations −61 %, chemistry time −24 %, total −18 % (two runs of each variant in
  parallel on idle nodes: 90.8/90.5 against 67.6/70.3 s per 200 steps). Differences to the default (T ≤ 2.5 % in the two zones at the
  dissociation front, abundances above 1e-8 up to 50 %) are those of a round-off-level perturbation (null test: `CHEMSTAG` alone):
  that layer is chaotic. Not yet in any chain (decision).
- With those on, `nlopc` (2–3 × `nlte_means` × 1000 frequencies per network zone per step) is half the chemistry time. Reusing the
  factors while all inputs stay within 1 % saved nothing (2 % of zone-steps qualify); not committed. Next candidates: fewer frequencies
  (accuracy test against the 1000-point means), or a fixed refresh interval.
- Machine-load scatter of whole-run timings is ~30 % between nodes; compare variants only in parallel runs, repeated.
