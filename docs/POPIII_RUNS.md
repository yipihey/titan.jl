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
| A3: core growth | Lagrangian mesh with a moving refinement window ahead of the star (`REZWIN=1.3,4,0.03 REZNMAX=340`), dynamical Δt limit relaxed (`CDYN=1000 DTRADMAX=0`): with `CDYN=0.01` the core's own density holds Δt at ~90 s | `phaseA.env` + `stageA3.env` | M* = 0.5 M☉ (`MSTOP`: mass inside the first interface falling faster than 2 km/s) |
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
  whole step and not on the state between the two halves of the operator split. Under test, not adopted.
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
CHEMSUB=1.3,3`, chain script `chain_c146.sh` in the session scratchpad; the restart at step 387050 continues the record
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
