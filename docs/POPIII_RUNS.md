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
