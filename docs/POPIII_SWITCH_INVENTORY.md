# titan.jl run-time switch inventory (branch `popiii-s6`, HEAD 326df84, 14 Nov 2026)

Read-only audit. Nothing in the repository or the run directories was changed or run.

## 0. Method and scope

- **What counts as a switch.** Every `ENV["X"]`, `get(ENV, "X", …)` and `haskey(ENV, "X")` in `julia/*.jl` and `julia/z0/*.jl`.
  That gives **393 distinct names**. The production driver does not use shared helpers. `envf`/`envi` exist only in
  `s5_run.jl:20-21`.
- **Production code path.** These are the files that `popiii_m2c.jl` includes. Some are included only conditionally: StellarEOS,
  NucNet, StellarMaterial and StarMaterial with `HOT=1`, `h2lines_gpu.jl` with `LINEGPU`, and `LyaTransfer.jl` through
  PopIIIMaterial. **273 of the names** are read on this path. Sections 1–13 cover all 273. The other 120 are read only by
  test, benchmark or old-stage scripts and are listed in section 14.
- **"Prod B"** is the effective value in the current production segment `z19t_c128cB5`, taken from its recorded `env`
  lines (`runs/t_c128c/meta_002165050_1512373.txt`). That file is the union of `runs/t_c128c/phaseB.env`, `run.env` and the
  `fullrun.sh` COMMON settings. Where a key appears twice, the later value wins (`CDYN=1000 … CDYN=3e5` gives 3e5). "–" means
  not set, so the code default applies. Note that the run keeps **its own** `phaseB.env`. It is older than the repo copy and
  lacks `EPTAB=1 LINEXC2=16 ENVCONS=1`, which the repo `scripts/popiii_run/phaseB.env` has.
- **"Prod A"** is phase A as the repo stage files define it: `phaseA.env` plus `stageA1/A2/A3.env`. `t_c128c` itself never ran
  phase A. It continued from the c128 checkpoint `z19c128B16` (procedure.log, 13 Nov 22:04), so phase A is the procedure that
  produced that history, not something re-run in t_c128c.
- **Origin** is the oldest commit from `git log -S'"NAME"'`. I searched on the quoted name to avoid false hits on short names
  such as `RAD`, `DH` and `IC`. **`1889108` (2026-10-20) is the root import commit.** For those switches the real origin
  predates the repository history, so it is shown as "root".
- **Verdict legend.** **R** = likely removable once the named general method is in. **N** = needed physics or needed
  infrastructure. **U** = unclear (quoted reasons). **K** = keep as tooling (I/O, diagnostic, performance; not a patch).
- **General methods referred to.** I call them **G1–G7** throughout:
  - **G1** Fully coupled implicit chemistry: `COUPLED` in every zone and for every species, including He, D and the
    three-level heats, so there is no split step and no lagged source.
  - **G2** One non-LTE material description everywhere: EOS and opacity from the populations the network carries. There is no
    LTE table and no table/network switch.
  - **G3** Error-controlled time stepping: an embedded local-error estimate of TR-BDF2 on all unknowns, including species,
    over the whole step. It replaces caps, stol, retry limits and halving factors.
  - **G4** Proper radiation transport and boundary: a consistent moment closure or multigroup, with the CMB entering only as
    the incoming boundary/background. It replaces CMB floors, flux ties and Eddington-factor clamps.
  - **G5** One consistent conservative energy formulation. A single discretisation keeps total energy conservative and
    recovers the gas first law, with no dual-energy row selection and no zero-point dependence.
  - **G6** A globalised nonlinear solver: line search or trust region on the full Newton system, with exact Jacobians. It
    replaces clamps, damping heuristics, decoupling and give-up rules.
  - **G7** One mesh from the IC to the end: a fixed log-Eulerian mesh or GRIDSEG. It replaces the stage chain, the remaps and
    rezoning. GRIDSEG is named in `docs/POPIII_GOAL.md` but is **not in the code**: there is no `GRIDSEG` string in `julia/`.

---

## 1. Run control, I/O and setup (category I/O unless noted)

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| SAVE | unset | popiii_m2c.jl:2296 | checkpoint name (end + every CKEVERY) | z19t_c128cB5 | per stage | root | K |
| CKEVERY | 500 | popiii_m2c.jl:1623 | checkpoint interval (steps) | – | – | root | K |
| CKDIR | lib/popiii/ck | popiii_m2c.jl:582 | checkpoint directory | run/ck | run/ck | root | K |
| RESTART | unset | popiii_m2c.jl:1289 | restart from checkpoint, with a CKENV consistency check (popiii_m2c.jl:597-605) | z19t_c128cB4_auto0 | per stage | root | K |
| RADFROM | unset | popiii_m2c.jl:1285 | build the radiation field from a hydro-only checkpoint (stage A2) | – | A2 | f75ab23 2026-10-31 | R (G7) |
| STARFROM | unset | popiii_m2c.jl:1283 | start from an S5 star state (`S5DIR`) | – | – | f4baf46 2026-11-01 | K |
| S5DIR | /tmp/s5 | popiii_m2c.jl:935 | directory of the S5 states for STARFROM | – | – | f4baf46 2026-11-01 | K |
| POPIII_OUT | lib/popiii | popiii_m2c.jl:46 | log/snapshot directory | run dir | run dir | c765bfd 2026-11-08 | K |
| POPIII_LIT | ~/Projects/popiii_lit | SCVH.jl:62 (also PopIIIEqEOS:121, ContinuumOpacity:71, H2Lines:32, NucNet:31) | literature-data directory | – | – | 07b074c 2026-11-10 | K |
| NUCNET_REACLIB | $POPIII_LIT/reaclib/… | NucNet.jl:31 | REACLIB file | – | – | root | K |
| RUNREC | "" | popiii_m2c.jl:311 | run recording (snaps_*.jls, perf.csv, meta) | run dir | run dir | 5337e9f 2026-11-09 | K |
| TAG | "" | popiii_m2c.jl:47 | log/snapshot file-name suffix | _t_c128c | yes | root | K |
| SNAPS | 0 | popiii_m2c.jl:436 | in-memory profile snapshots every 0.25 dex n_c | 0 | 0 | root | K |
| SNAPEVERY | 0 | popiii_m2c.jl:437 | snapshot every n steps (RUNREC) | 50 | 20 | root | K |
| PRINTEVERY | 500 | popiii_m2c.jl:1623 | status-line interval | 100 | 100 | f75ab23 2026-10-31 | K |
| MAXSTEPS | 200000 | popiii_m2c.jl:2056 | steps per process (segment length) | 100000 | 1e7 | root | R (G7: one run) |
| TSTOP | 6000 | popiii_m2c.jl:577 | stop at central T | 2e8 | 2e8, A2: 5e4 | root | R (G7) |
| NSTOP | 1e18 | popiii_m2c.jl:577 | stop at central n | 1e27 | A1: 3.5e18 | root | R (G7) |
| MSTOP | 0 | popiii_m2c.jl:1621 | stop at M★ (A3→B hand-over) | – | A3: 0.05 | 35ef587 2026-11-09 | R (G7) |
| DTSTOP | 0 | popiii_m2c.jl:1622 | stop when Δt collapses below this (protects checkpoints) | 1e-9 | 1e-2 | 4b90c91 2026-11-01 | K |
| TITAN_NOLIB | 0 | TitanFortran.jl:161 | pure-Julia backend, no libtitan | 1 | 1 | root | K |
| TITAN_MGR | unset | TitanFortran.jl:99 | grid array size (needs NOLIB) | 600 | 600 | root | K |
| TITAN_BITWISE | 0 | TitanDriver.jl:100 | bitwise comparison mode vs Fortran (FIX_CLOCK/FIX_ENERGY off) | – | – | root | K |
| MSG | 0 | popiii_m2c.jl:279 | Newton messages to stdout | 1 | 1 | f75ab23 2026-10-31 | K |
| DEBUG | 0 | popiii_m2c.jl:2300 | print backtrace on exception | – | – | root | K |
| Z0_TABLE | "" | Z0Physics.jl:17 | `ext` = extended Z=0 opacity table | – | – | root | U |
| REFINE | unset | popiii_m2c.jl:1288 | static pre-refinement mlo,mhi,frac at the A2 rebuild | – | A2: 1e-3,0.1,0.03 (fullrun.sh) | f42916f 2026-11-01 | R (G7) |

## 2. Problem set-up and initial conditions (physics model)

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| IC | p100 | popiii_m2c.jl:61 | initial condition (`abn02z19` = ABN02 z=19 cloud) | abn02z19 | abn02z19 | root | N |
| MCLOUD | 1e4 (z19) | popiii_m2c.jl:133 | cloud mass [M☉] | – | – | root | N |
| DM | 1 | popiii_m2c.jl:116 | ABN02 dark-matter potential on | – | – | root | N |
| LRJOIN | auto | popiii_m2c.jl:165 | join radius of the hydrostatic interior with the ABN02 n(r) | – | – | root | N |
| DLNJOIN | 0.15 | popiii_m2c.jl:169 | tolerance [dex] of the automatic join | – | – | root | N |
| XE0 | 1e-4 (z19) | popiii_m2c.jl:205 | initial x_e ("[assumed; ABN02 does not show it]") | – | – | root | N |
| K5 | pss | popiii_m2c.jl:227 | choice of rate k5 | – | – | root | N |
| CHEM | rhfc02 | popiii_m2c.jl:244 | material: RHFC02 PrimordialGas or ChemNet `net*` | netmodernli | netmodernli | root | N |
| RHOEQ | 5e-11 | popiii_m2c.jl:246 | ρ of the equilibrium switch (PrimordialGas only; unused with netmodernli) | – | – | root | R (G2) |
| EQEOS | scvh | popiii_m2c.jl:246 | equilibrium EOS of PrimordialGas (unused with net*) | – | – | root | U |
| CHEM_YHE | 1/12 | ChemNet.jl:32 | He/H by number (0 = pure-H benchmarks) | – | – | e96c96c 2026-11-04 | N |
| DH | 2.5e-5 | popiii_m2c.jl:651, PopIIIMaterial.jl:1238 | D/H for NucNet and HD chemistry | – | – | root | N |
| TCMB | 0 | popiii_m2c.jl:39 | T_CMB: cooling floor `cool_floor` (PrimordialChem.jl:313, Λ(T)−Λ(T_CMB)), outer radiation BC `ximr`, starlight E_cmb, LINEINC background, HD | 54.5 | 54.5 | root | N (value); the subtractive floor is a patch, see §12 |
| HDINIT | eq | PopIIIMaterial.jl:1239 | initial HD/D⁺ state at HDCHEM start | – | – | d97be5f 2026-11-11 | K |

## 3. Physics-model switches (adding or choosing physics)

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| RAD | 0 | popiii_m2c.jl:43 | TITAN radiation E_r, F_r on; RA04 CIE cooling off ("must never be set from the IC without the rebuild", fullrun.sh) | 1 | A1: 0, A2/A3: 1 | f75ab23 2026-10-31 | N (the staged on/off is R, G7) |
| LTRAN | 1 | popiii_m2c.jl:733 | 1 = formal solution for f, g; 2 = Eddington approximation | – | – | f75ab23 2026-10-31 | N |
| HOT | 0 | popiii_m2c.jl:22 | stellar EOS, nuclear network and MLT in the hot core (T > HOTT1) | 1 | 1 | ee1bedc 2026-11-01 | N |
| HOTT1 / HOTT2 | 3e5 / 9e5 | popiii_m2c.jl:665 | blend window from the PopIII EOS to the stellar EOS (smoothstep, PopIIIMaterial.jl:139-142) | – | – | ee1bedc 2026-11-01 | U (a blend between two EOS; G2 would make it one EOS) |
| ALPHAMLT | 1.5 | popiii_m2c.jl:739 | mixing-length α | 1.5 | 1.5 | ee1bedc 2026-11-01 | N |
| COLDMLT | 0 | popiii_m2c.jl:740 | MLT also in hot table zones above this T; turns on MLT_NOINFALL (no MLT flux where infall > 30 km/s, TitanResidFast.jl:635) | 1e4 | 1e4 | 452ec51 2026-11-01 | U |
| STARLIGHT | 0 | popiii_m2c.jl:403 | star light below 13.6 eV in network zones (H⁻ detachment, LW, optional H(n≥2)) | 1 | 1 | 66737c6 2026-11-05 | N (physics), but its J_ν = W B(T_col) closure is a patch, see §5 |
| SLHM | 1 | popiii_m2c.jl:405 | H⁻ photodetachment by star light | – | – | 66737c6 2026-11-05 | N |
| SLH | 0 | popiii_m2c.jl:406 | H(n≥2) photoionisation via network ("energy-unsafe") | – | – | 66737c6 2026-11-05 | R (superseded by SL3) |
| SLSHIELD | static | popiii_m2c.jl:407 | LW shielding: static / sobolev / none / off | sobolev | sobolev | 66737c6 2026-11-05 | N |
| SL3 | 0 | popiii_m2c.jl:251 | three-level H atom (Peebles C, Sobolev Lyα) in lit network zones | 1 | 1 | d7373a7 2026-11-06 | N |
| NLOPAC | 0 | popiii_m2c.jl:422 | network zones: non-LTE opacity = LTE table × lagged non-LTE ratio | 1 | 1 | 72a27c2 2026-11-06 | R (G2) |
| NLOPACJ | 0 | popiii_m2c.jl:425 | κ_E weighted with B_ν(T_c) of star light instead of B_ν(T_gas) | – | – | 69bef61 2026-11-06 | R (G4: multigroup) |
| LYABB | 0 | PopIIIMaterial.jl:745 | bound–bound energy exchange of the three-level atom ("missing from tl3_heat's ionization bookkeeping") | 1 | 1 | 7e763d6 2026-11-08 | N (a missing term; would be native in G1+G2) |
| TL3RECKT | 0 | PopIIIMaterial.jl:696 | recombination kinetic energy ⟨hν⟩(T)−χ₂ instead of 0.75 kT (otherwise "cycle heat does not vanish in LTE") | – | – | 8807830 2026-11-14 | N (correctness fix; should be default) |
| LOCALJ | 0 | PopIIIMaterial.jl:711 | network photo-processes see W B(T_c) instead of the CMB-only `m.env` | – | – | 3b8cc2b 2026-11-14 | N (physics gap: G4 per-zone J_ν) |
| H2PPD | 0 | ChemRates.jl:108 | H₂⁺ photodissociation by detailed balance | – | – | 8807830 2026-11-14 | N |
| HECHEM | 0 | PopIIIEqEOS.jl:41, ChemNet.jl:74 | He, He⁺, He²⁺ as network species; He table EOS | 1 | – | 8a7fbd8 2026-11-12 | N |
| HDCHEM | 0 | PopIIIMaterial.jl:1237 | D⁺/HD trace network and HD cooling | 1 | – | d97be5f 2026-11-11 | N |
| CHEMCK | 0 | ChemNet.jl:548 | reduced network (H₂, H⁺ slow; fast species algebraic) as in ChemistryKernels | 1 | – | 6bf97a8 2026-11-13 | N (model reduction underlying COUPLED) |
| CHEMCKSWEEPS | 3 | ChemNet.jl:548 | Gauss–Seidel sweeps for fast-species steady state | – | – | 6bf97a8 2026-11-13 | U |
| NLTE | 1 | popiii_m2c.jl:441 | H₂ line source function: 1/0/ep/ali/hyb | hyb | hyb | root | N |
| NALI | 20 | popiii_m2c.jl:441 | ALI refresh interval in `hyb` | – | – | root | U (lag scheme) |
| NBIN | 30 | popiii_m2c.jl:576 | frequency bins of line transfer | – | – | root | N |
| LINECOLL | 0 | popiii_m2c.jl:448 | gas term = net collisional exchange at ALI populations instead of −∇·F of the formal solution | 1 | 1 | f6f2bad 2026-11-08 | N |
| LINEINC | 0 | popiii_m2c.jl:456 | H₂ levels see the incident CMB plus diluted photosphere ("otherwise … line cooling can take the gas below T_CMB") | 1 | – | b6fa7fd 2026-11-11 | N |
| LINESRC | linear | H2Lines.jl:238 | zone-constant line source | zone | zone | f6f2bad 2026-11-08 | N |
| LINEQUAD | mu | H2Lines.jl:242 | flux quadrature trapezoid in μ² (p²) | p2 | p2 | f6f2bad 2026-11-08 | N |
| LINECORE | 1 | popiii_m2c.jl:493 | number of transparent merged-core zones for impact parameters | 8 | 8 | f6f2bad 2026-11-08 | N |
| LTAUCUT | 10 | popiii_m2c.jl:493 | H₂ line transfer only outside continuum τ (commit "outside the continuum photosphere (10x)") | – | – | f0bee17 2026-10-31 | U |
| EPGRAD | thermal | H2Lines.jl:602 | Sobolev gradient over one thermal width (fix for core ringing, POPIII_RUNS.md:76-95) | – (=thermal) | thermal | 73491cc 2026-11-09 | N (now the default) |
| LYART + LYATAU, LYAGAMMA, LYAVMAX, LYAH, LYANCORE, LYASUBTOP, LYANSUB, LYAMAXIT, LYATOL, LYANG | 0; 20, 1, 2e5, 0.2, 6, 30, 24, 40, 1e-3, 4 | LyaTransfer.jl:26-28, 36, 69, 98, 168-169 | Lyα transfer for the three-level atom replacing Sobolev escape (Voigt, ALI+Ng, sublayers) | – | – | d565e21 2026-11-08 | N (off; physics option) |
| RAY, RAYTMIN, RAYXMIN | 0, 1.5e4, 1e-6 | popiii_m2c.jl:376, 378, 399 | photon-conserving ionising ray from the photosphere (M3); a zone the ray can ionise by > RAYXMIN stays on the network | – | – | e96c96c 2026-11-04 | N (the H II region is the expected end of accretion, POPIII_STRATEGY.md) |
| SHOCKRAD, SHOCKDU | 0, unset | popiii_m2c.jl:270, 276 | sub-grid radiative shock | – | – | a8eb0ae 2026-11-04 | U |
| INTEGRATOR | trbdf2 | popiii_m2c.jl:32 | TR-BDF2 or TITAN θ scheme | – | – | root | N |
| THETA | 0.55 | popiii_m2c.jl:291 | θ of TITAN scheme ("0.5 lets the centre ring") | – | – | root | U |
| DXMAX, EFRAC | 4e-3, 5 | popiii_m2c.jl:689 | HOT burning: composition cap and energy fraction of the S5 split | – (active, HOT=1) | – | root | U (G1/G3 would absorb it) |
| SPLITTOL, SPLITIT, ZSPLITTOL | 3e-3, 4, 3e-3 | popiii_m2c.jl:692 | HOT: tolerance and iterations of the burn/hydro source iteration | – (active) | – | root | R (G1: nuclear sources in Newton) |
| S5NUM, S5TOL | 0, 0.03 | popiii_m2c.jl:695, 701 | S5 quasi-static numerics (predictor, EXTRAP_FRAC, own Δt rule) | – | – | f4baf46 2026-11-01 | R (G3) |

## 4. Table ↔ network switching rules (physics patch)

These rules all decide, per zone and per step, whether the gas is the LTE equilibrium table (SCVH/Saha; latent heat inside
Newton) or the ideal-gas ChemNet network. The decision is in PopIIIMaterial.jl:1033-1060.

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| PSWITCH | 1e4 | popiii_m2c.jl:258 | table where p_gas ≥ PSWITCH ("EOS validity, RHFC02's 1e4 dyn cm⁻²") | 300 | 300 | root | R (G2) |
| PSWHYS | 1 | popiii_m2c.jl:846 | stay on table until p < PSWHYS × pswitch | 0.3 | 0.3 | c458b57 2026-11-01 | R (G2) |
| PSWTMIN | 0 | popiii_m2c.jl:847 | table also needs T ≥ PSWTMIN | 4500 | 4500 | a8eb0ae 2026-11-04 | R (G2); Tom set 0 in z100 (POPIII_RUNS.md:2049) |
| PSWTEX | 0 | popiii_m2c.jl:848 | table zone stays on table until T < PSWTEX at any p | 5000 | 5000 | 21cac65 2026-11-05 | R (G2) |
| PSWSTAY | 1 | popiii_m2c.jl:430 | 0 drops the PSWTEX stay rule (only τ rule) | – (=1) | – | 69bef61 2026-11-06 | R (G2) |
| PSWTAU | 0 | popiii_m2c.jl:848 | table where T ≥ PSWTEX and τ_R from surface ≥ PSWTAU (half to stay) | 3 | 3 | ef988af 2026-11-05 | R (G2) |
| PSWAGREE | 0 | popiii_m2c.jl:423 | table below pswitch only where table and three-level steady state agree to ε | – | – | 72a27c2 2026-11-06 | R (G2) |
| PSWNH | 0 | PopIIIMaterial.jl:505 | table also where n_H ≥ n (non-ideal regime) | – | – | 8807830 2026-11-14 | U (non-ideal EOS physics; G2 needs non-ideal network EOS) |
| PSWRAMP | 0 | PopIIIMaterial.jl:498 | gradual table→network hand-over, one zone per n steps | – | – | 8807830 2026-11-14 | R (G2) |
| PSWRAD | 0 | PopIIIMaterial.jl:550 | enter table only if \|T_rad/T−1\| ≤ ε | – | – | 6478132 2026-11-10 | R (G2) |
| (hard-wired) PSWTAUDEEP | 3.0 | PopIIIMaterial.jl:521 | table always at τ ≥ 3 in PSWAGREE path | n/a | n/a | 72a27c2 / 69bef61 | R (G2) |
| CHEMFBDENSE | 1 | PopIIIMaterial.jl:585 | failed network step at p ≥ pswitch → table-equilibrium fallback | – (=1) | – (=1) | ce95c1b 2026-11-09 | R (G1/G6) |
| CHEMTMIN | 0 | PopIIIMaterial.jl:948 | chemistry step returning T < floor counts as failed → table fallback | – | – | deddb60 2026-11-12 | R (G1/G6) |

## 5. Opacity / star-light closures in network zones (physics patch)

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| NLOIMPT | 0 | PopIIIMaterial.jl:372 | non-LTE opacity = OPAC(ρ,T_iterate) × extrapolated lagged ratio | 1 | – | 8807830 2026-11-14 | R (G2) |
| NLOFIX | 0 | PopIIIMaterial.jl:426 | floor χ of LTE-table network zones at σ_T n_e/ρ of the zone's own electrons | 1 | – | 34450da 2026-11-13 | R (G2) |
| NLOTMAX | Inf | PopIIIMaterial.jl:433 | above T, plain table opacity instead of the non-LTE ratio | – | – | c3de545 2026-11-14 | R (G2) |
| NLODAMP | 1 | PopIIIMaterial.jl:260 | under-relaxation in time of lagged NLO factors ("Test switch") | – | – | f4c0389 2026-11-12 | R (G1/G2) |
| CPLNLOX | 1 | PopIIIMaterial.jl:1352 | coupled zones: opacity follows the iterate's x_H⁺ (d ln κ/d ln x) | – (=1) | – | f6292e4 2026-11-13 | R (G2: opacity from populations) |
| SLBLEND | 0 | popiii_m2c.jl:424 | star-light colour temperature and dilution blended smoothly in τ (0.3…3) | 1 | 1 | 72a27c2 2026-11-06 | R (G4) |
| SLEB | 0 | popiii_m2c.jl:427 | colour temperature = local T_r at τ_k + 2/3 (Eddington–Barbier) | 1 | 1 | 69bef61 2026-11-06 | R (G4) |
| SLWFLUX | 0 | PopIIIMaterial.jl:435 | dilution W from max(c(E_r−E_cmb), F_r) | 1 | – | 3dbf916 2026-11-13 | R (G4) |
| CHICAP | 0 | PopIIIMaterial.jl:390 | cap χ ("Diagnostic switch … Not physics") | – | – | 271e101 2026-11-10 | K (diag) |
| OPACSENS | 1,1 | popiii_m2c.jl:725 | κ_cond, κ_rad sensitivity factors (test) | – | – | 9f44a74 2026-11-02 | K (diag) |

## 6. Operator-splitting bookkeeping (numerics workaround)

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| COUPLED | 0 | popiii_m2c.jl:1086, TitanFortran.jl:104 | x_H⁺, x_H₂ as Newton unknowns 8, 9 (implies CHEMCK) | 1 | – | ab256c2 2026-11-13 | N (it is G1, partial) |
| CPLBAND | "" (all) | PopIIIMaterial.jl:1359 | which zones are coupled (hss=f / r=r1,r2 / z=k1,k2); the rest stay split | r=1500,1e6 | – | e982bb0 2026-11-13 | R (G1 everywhere) |
| CPLBADDELTA | 100 | TitanResidFast.jl:273 | decouple a zone for the rest of a failed step if \|δ ln x\| > 100 | – (=100) | – | f72ba29 2026-11-13 | R (G6) |
| CPLCLAMP | 1 | TitanResidFast.jl:222 | clamp of species correction per Newton iteration (ln x) | 3 | – | fe4b3a1 2026-11-13 | R (G6) |
| CPLTHETA | 1 | TitanResidFast.jl:228 | backward-Euler chemistry also in the TR stage | – (=1) | – | ec9c58d 2026-11-13 | U (L-stability need; G3 would choose) |
| CPLIMPQ | 0 | TitanResidFast.jl:220 | three-level heats cq, yq of the iterate in the rows (experimental) | – | – | f6292e4 2026-11-13 | R→should become default under G1 (docstring: "~300 K below the split scheme's … unexplained") |
| CPLLOCAL | 1 | popiii_m2c.jl:2039 | initial species guess from the zone's species rows solved alone | – (=1) | – | ec9c58d 2026-11-13 | U (predictor; G6) |
| ECHEMSPLIT | 0 | popiii_m2c.jl:1079 | hydro advects e−Σx E0; chemical energy moves with species | 1 | – | fab1a39 2026-11-13 | R (G1/G5) |
| XADV2 | = SL3 | popiii_m2c.jl:1032 | explicit 2nd-order (van Leer) correction to species advection | – (=1 via SL3) | – (=1) | 89c218f 2026-11-06 | R (G1: species advected in the Newton with e's scheme) |
| XADVU | 0 | popiii_m2c.jl:1070 | 1: diagnose the e/x mismatch after advection; 2: reset e (not conservative) | – | – | 883b4de 2026-11-13 | R |
| BALSPLIT | 0 | popiii_m2c.jl:412 | balanced splitting (Speth et al. 2013) of chemistry heat | 0 | 0 | c559871 2026-11-08 | R (G1) |
| CHEMSUB | 0 | PopIIIMaterial.jl:860 | sub-step a three-level chemistry step changing T > fT or n_e > fx | 1.3,3 (**set but bypassed**: skipped when CHEMCK && CHEMCKHSS>0, PopIIIMaterial.jl:1123) | 1.3,3 (active) | 6478132 2026-11-10 | R (G1) |
| CHEMSUBMAX | 16 | PopIIIMaterial.jl:861 | max bisection depth of CHEMSUB | – | – | 6478132 2026-11-10 | R |
| CHEMCKHSS | 0 | ChemNet.jl:564 | H⁺ quasi-steady state where it relaxes fast; self-consistent n_e in C | 100 | – | 79ee3ca 2026-11-13 | R (G1) |
| CHEMCKADAPT | "" | ChemNet.jl:555 | adaptive BE sub-steps inside CHEMCK step | – | – | dc6548c 2026-11-13 | R (G1/G3) |
| CHEMSUBCYC, CHEMSUBMAX_SC | 0, 200000 | ChemNet.jl:414-415 | sub-cycled chemistry (Anninos 1997) | – | – | deb6f13 / 59ac3c2 2026-11-13 | R |
| IMPCOOL / IMPALWAYS | 1 / 0 | popiii_m2c.jl:235, 237 | H₂ line cooling Λ_j (T/T_j)^s implicit in chemistry; slope by a second transfer call | – (=1) | – (=1) | root / 73491cc | R (G1: Λ in Newton) |
| LINEFRAC | 0 | popiii_m2c.jl:1659 | reuse Λ until LINEFRAC × t_cool elapsed | 0 | 0 | 225ef23 2026-11-04 | R |
| ALIWMIN, ALIGROW, ALISHRINK, ALICAP | 0, 0.01, 0.05, NALI | popiii_m2c.jl:441 | adaptive ALI refresh interval in NLTE=hyb | 1e-3, 0.1, 0.5, – | same | 38b8e81 2026-10-31; ALICAP b1b766e 2026-11-11 | U (lagged-line scheme; R if lines enter Newton) |

## 7. Chemistry-Newton workarounds (numerics workaround)

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| CHEMROBUST | 0 | ChemNet.jl:384 | non-finite Jacobian = failed attempt; clip trace Li | – | – | 7d91552 2026-11-06 | R (G6) |
| CHEMTRCLIP | 0 | ChemNet.jl:393 | per-species positivity instead of common step length | – | – | 7930234 2026-11-13 | R (G6 / log variables) |
| CHEMLOGX | 0 | ChemNet.jl:401 | Newton in ln x | – | – | 06db35d 2026-11-13 | U (a reformulation, possibly the general fix) |
| CHEMSTAG | 0 | ChemNet.jl:521 | accept stagnated iterate from iteration n | – | – | 57d6f45 2026-11-12 | R (G6) |
| CHEMGIVEUP | 0 | ChemNet.jl:527 | abandon attempt early at iteration n | – | – | 57d6f45 2026-11-12 | R (G6/G3) |
| (hard-wired) stagnation acceptance after 60 iterations | – | POPIII_RUNS.md:96-100 ("not an option: it applies to every run") | | on | on | 73491cc | R (G6) |

## 8. Radiation numerics (numerics workaround / physics fix)

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| FLUXTIE | 0,0 (off) | popiii_m2c.jl:1735; TitanDriver.jl:87 | thin zones: momentum row → F = h c E (upwinded), weight by τ between τ1, τ2 | 0.1,1 | 0.1,1 | a85f13b 2026-11-08 | R (G4) |
| EDDGMU0 | 0 | popiii_m2c.jl:1741; TitanCore.jl:101 | outer g integrated to μ = 0 | 1 | 1 | 6fcf304 2026-11-05 | N (a quadrature correctness fix; should be default) |
| EDDSCAT | 0 | popiii_m2c.jl:1734; TitanCore.jl:96 | formal solution: all χ above κ_P is scattering | 1 | 1 | 6fcf304 2026-11-05 | R (G4/G2: true absorption/scattering split from the material) |
| EDDCHI | 0 | popiii_m2c.jl:1744 | extinction ≤ κ_P + σ_e in formal solution | – | – | 6fcf304 2026-11-05 | R |
| EDDREAL, EDDTAU, EDDRELAX | 0, 0.3, 1 | popiii_m2c.jl:1635-1636 | enforce f ≥ (F/cE)² above τ = EDDTAU, under-relaxed | – | – | ef988af / 2ae79bc 2026-11-05 | R (G4) |
| ERTAUCONV | 1.0 (**on by default**) | popiii_m2c.jl:734; TitanCore.jl:91 | causality floor E_r ≥ max(F,0)/c, 0.1E in thin zones (τ < value) | 0 (off) | 0 (off) | 4b90c91 2026-11-01 | R (G4) |
| ERTAUMIN | 1.0 | popiii_m2c.jl:735; TitanCore.jl:191 | zones with τ_surface < value excluded from E_r Δt criterion | 3000 | 3000 | f75ab23 2026-10-31 | R (G3) |
| ERTHICK | unset | popiii_m2c.jl:1652; TitanCore.jl:189 | zones thicker than value excluded from E_r Δt criterion | – | – | 225ef23 2026-11-04 | R (G3) |
| EXCAP | unset | popiii_m2c.jl:1658; TitanResidFast.jl:488 | scale gas–radiation exchange so c κ_E ρ Δt ≤ EXCAP | 1e8 | – | 225ef23 2026-11-04 | U (round-off conditioning; a scaled/1-T formulation would remove it) |
| DTRADMAX | 1500 (**on by default**) | popiii_m2c.jl:742 | Δt cap for lagged Eddington factors | 0 (off) | A2: 1500; A1/A3: 0 | f75ab23 2026-10-31 | R (G3/G4) |

## 9. Energy formulation (numerics workaround / physics fix)

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| ETOTAL | 0 | popiii_m2c.jl:1657 | TITAN conservative total-energy row; gas row kept where e/(e+KE) < ETASW | 1 | – (off: broke collapse, POPIII_RUNS.md:33-38) | 225ef23 2026-11-04 | U (keep conservation; replace the selection, G5) |
| ETASW | 0.1 | popiii_m2c.jl:1657 | dual-energy threshold | 0.6 | – | 225ef23 2026-11-04 | R (G5) |
| ETASWTH | unset | popiii_m2c.jl:1654; TitanResidFast.jl:68 | also gas row where 3p/2ρ < f × (3p/2ρ + E_r/ρ + K/m) | 1e-3 | – | 508105c 2026-11-11 | R (G5/G4) |
| ETASWTIE | unset | popiii_m2c.jl:1655; TitanResidFast.jl:76 | also gas row where FLUXTIE weight ≥ value | 0.5 | – | b8c7d8d 2026-11-11 | R (G4: with no FLUXTIE it has no purpose) |
| ETOTGRAVGRID | 0 | popiii_m2c.jl:1656; TitanResidFast.jl:1325 | gravity work on moving-face mass in total-energy row | – | – | 2395554 2026-11-09 | N if ETOTAL is kept on Lagrangian zones (missing term) |
| ADVREL | 0 | TitanResidFast.jl:134 | unconverged relative mass flux carries momentum and energy | – (fixed mesh: inert?) | 1 | 73491cc 2026-11-09 | N (a consistency fix; should be default) |
| TRESET | unset | popiii_m2c.jl:1527 | one-time set of zone T at restart (repair; "not conservative") | – | – | 397801f 2026-11-11 | R (G4/G5) |

## 10. Hydro / artificial viscosity / advection (numerics)

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| AVZONE | 0 | popiii_m2c.jl:277 | viscous length = AVZONE × zone width | 2 | 2 | root | N (AV choice) |
| AVQL1, AVCQ2 | 0.15, 1.0 | popiii_m2c.jl:266 | viscous length ql1 r̄; quadratic coefficient | – | – | root | N |
| AVTENSION | 1 (TITAN) | popiii_m2c.jl:271; TitanResidFast.jl:318 | 0 = no viscous tension in zones that shrink but stretch radially | 0 | 0 | a8eb0ae 2026-11-04 | U (POPIII_RUNS.md:507 "The viscosity term is numerical") |
| AVCAP | 0 | popiii_m2c.jl:275 | cap viscous length at AVCAP × r̄ | – | – | a8eb0ae 2026-11-04 | R |
| UPFREEZE | 0 | TitanResidFast.jl:502 | freeze upwind direction from Newton iteration n | – | – | 271e101 2026-11-10 | R (G6: smooth upwinding) |

## 11. Mesh and remap (numerics; G7 targets)

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| ENVEUL | unset | popiii_m2c.jl:956 | remap to a fixed log mesh (lgrid=2) at restart; also turns on GASROW, INNER_SINK | 0,0.03,0.1,1000 | – | a8eb0ae 2026-11-04 | R (G7: start on that mesh) |
| ENVCORE | 0.1,0.08,10 | popiii_m2c.jl:959 | first interface, Δln r, core boundary for a regular centre | – (default) | – | a8eb0ae 2026-11-04 | N (mesh parameters) |
| ENVCONS | 0 | popiii_m2c.jl:960 | mass-conserving remap | – (not in run's phaseB.env) | – | 20ad4c2 2026-11-11 | R (G7) |
| ENVSINKTAU | 1e6 | popiii_m2c.jl:961 | relaxation time of the inner sink (r_in > 0) | – (centre: inert) | – | a8eb0ae 2026-11-04 | U |
| ENVMS0 | unset | popiii_m2c.jl:1294 | reference shell mass for the sink at restart | 0.4737 | – | a8eb0ae 2026-11-04 | R (G7) |
| ENVDT0 | 100 / 1e4 | popiii_m2c.jl:1274, 1366, 1451 | first Δt after remap/split/merge | – | – | a8eb0ae 2026-11-04 | R (G3/G7) |
| ENVSPLIT, ENVMERGE | unset | popiii_m2c.jl:1497-1498 | one-time split/merge of fixed-mesh zones at restart | – | – | 4095a0a 2026-11-08 / eedb250 2026-11-10 | R (G7) |
| GRIDMAP (+GMTAU, GMWIN, GMSMOOTH, GMJUMP) | unset/1/0 | popiii_m2c.jl:976, 986, 999, 1011, 984 | one moving mesh following radius and mass (arc length) | – | – | f96d334, 15dcf8f, adf0bdf, 54af9e5, 0b98e65 (2026-11-12) | U (candidate G7) |
| GRIDWIN | unset | popiii_m2c.jl:1537 | Dorfi–Drury window on the fixed mesh | – | – | f383e69 2026-11-08 | U |
| REZONE, REZEVERY, REZJSPLIT, REZJMERGE, REZDMMIN, REZQMAX, REZNMAX, REZFRAC, REZKEEPC | 0, 20, 0.3, 0.05, 1e-7, 2.5, 260, 0.01, 5 | popiii_m2c.jl:842-845 | conservative mass rezoning (Rezone.jl) | REZONE=0 (off) | 1; REZJSPLIT=Inf; A3: REZNMAX=340 | c81412a 2026-10-31 | R (G7) |
| REZWIN | unset | popiii_m2c.jl:849 | moving refinement window ahead of the star | – | A3: 1.3,4,0.03 | c458b57 2026-11-01 | R (G7) |
| REZDR, REZPSPLIT, REZSLOPE | 0, 1, "" | popiii_m2c.jl:855, 888, 856 | split by Δr/r; pressure-consistent split; slopes per ln r | – | – | a8eb0ae / f0bee17 / 74c679c | R (G7) |

## 12. Time-step control

| Switch | Default | Where | What it does | Prod B | Prod A | Origin | V |
|---|---|---|---|---|---|---|---|
| CDYN | 0.01 | popiii_m2c.jl:292; TitanCore.jl:41 | Δt ≤ CDYN × min (r³/Gm)^½ over the mesh | 3e5 | A1/A2: 0.01; A3: 1000 | root | R (G3) |
| DTENV | unset | popiii_m2c.jl:1485 | Δt ≤ C × dynamical time over r1…r2 only | – (removed 14 Nov) | – | 4095a0a 2026-11-08 | R (G3) |
| STOL | (TITAN) | popiii_m2c.jl:1626 | TITAN max relative change per step | 0.3 | 0.1 | root | R (G3) |
| STOLCHEM | 0 | popiii_m2c.jl:1630 | re-evaluate stol after chemistry (whole split step) | 1 | – | 271e101 2026-11-10 | R (G1+G3) |
| CPLXSTOL | ln 3 | popiii_m2c.jl:1087 | Δt cap: ln x_H⁺ change ≤ value in coupled zones | 3 | – | ab256c2 2026-11-13 | R (G3) |
| NITER | (TITAN) | popiii_m2c.jl:1625 | Newton iterations per attempt | 25 | 10 | root | R (G3/G6) |
| NTRY | (TITAN) | popiii_m2c.jl:1631 | attempts (halvings) before giving up | 60 | 10 | a2d829a 2026-11-04 | R (G3/G6) |
| FAILTFAC | (TITAN) | popiii_m2c.jl:1632 | Δt factor after failed attempt | 0.5 | – | root | R (G3) |
| DAMPNEWTON, DAMPALPHA, DAMPTRIG | 0, 0.3, 0.1 | Trbdf2.jl:70-72 | damped Newton for stage solves | – | – | a85f13b 2026-11-08 | R→G6 (a crude version of it) |
| DTCOOL | 1.0 | popiii_m2c.jl:236 | Δt ≤ DTCOOL × t_cool when implicit cooling (else 0.1 t_cool hard-wired, popiii_m2c.jl:2150) | – (=1) | – | root | R (G3) |
| TCOOLTH | 0 | popiii_m2c.jl:242 | t_cool from 3p/2ρ instead of e + UOFF_TH | – | 1 | 73491cc 2026-11-09 | R (G3) |
| DTMAXTEST | unset | popiii_m2c.jl:2029 | test: Δt ≤ value (Δt → 0 limit) | – | – | fd64528 2026-11-13 | K (test) |
| (see §8) DTRADMAX, ERTAUMIN, ERTHICK | | | | | | | |

## 13. Performance (claimed bit-identical or round-off only) and diagnostics

**Performance (K).**

| Switch | Default | Where | Claim | Prod B | Prod A | Origin |
|---|---|---|---|---|---|---|
| MATTHREADS | 0 | TitanResidFast.jl:561 | "bit-identical to the serial sweep" | 4 | 16 (the run's phaseB.env has 16, run.env overrides to 4) | f57b27e 2026-11-10 |
| CHEMTHREADS | 0 | PopIIIMaterial.jl:840 | zones independent | 1 | 1 | ef4ba2c 2026-11-09 |
| CHEMSCHED | zone | PopIIIMaterial.jl:847 | "the results are the same" | – | – | f57b27e 2026-11-10 |
| LINESCHED | zone | H2Lines.jl:604 | scheduling | – | – | f57b27e 2026-11-10 |
| LINEFAST | 0 | H2Lines.jl:431 | "agree to round-off, not to the last bit" | 1 | 1 | a614b3e 2026-11-11 |
| LINEXC2 | 38 | H2Lines.jl:432 | 16 "changes the cooling by ≲ 3e-9" (not bit-identical) | – (38) | 16 | 50c85cc 2026-11-11 |
| EPTAB | 0 | H2Lines.jl:629 | "the same to the last bit" | – | 1 | 1ba004f 2026-11-11 |
| LINEGPU, LINEGPUFLOAT | unset, 32 | h2lines_gpu.jl:45, 35 | GPU formal solution (Float32 by default: not identical) | – | – | a5eeef4 2026-11-11 |
| THCHAIN | 1 | popiii_m2c.jl:635 | "same values and same first derivatives" | – (=1) | – | ef4ba2c 2026-11-09 |
| CPLJLAG | 3e-3,0.1 | PopIIIMaterial.jl:1393 | inexact Newton, "converged state is unchanged" | – (on) | – | e982bb0 2026-11-13 |

**Diagnostics (K, all off in production except DTDIAG and TIMING).** Name: default, location, origin.

- Chemistry: `CHEMDIAG` 0 popiii_m2c.jl:1617 (c458b57); `CHEMNANDBG` 0 ChemNet.jl:509 (7d91552); `CHEMTRACE` 0 ChemNet.jl:531 (7d91552);
  `CHEMSTATS` 0 ChemNet.jl:514 (d0f34ea); `CHEMREPLAY` 0 PopIIIMaterial.jl:867 (13a6a96); `CHEMREPLAYTRACE` 0 PopIIIMaterial.jl:928
  (acc42e0); `CHEMZPROF` 0 PopIIIMaterial.jl:952 (6fbdd4b); `TL3DBG` 0 PopIIIMaterial.jl:640 (d7373a7); `ZONEDBG` 0 PopIIIMaterial.jl:642 (03c8987).
- Coupled solve: `CPLDBG` 0 TitanResidFast.jl:222 (ec9c58d); `CPLDBGZ`/`CPLDBGN` 0/20 PopIIIMaterial.jl:1353 (fd64528); `CPLJCHECK`
  popiii_m2c.jl:2042 (ab256c2); `CPLPROF` 0 popiii_m2c.jl:1985 and `CPLPROFILE` popiii_m2c.jl:1995 (e982bb0); `CPLSCAN` popiii_m2c.jl:2053 and `CPLSCANNLO`
  0 popiii_m2c.jl:1967 (8807830).
- Newton/Jacobian: `NEWTKX` 0 popiii_m2c.jl:1661 (225ef23); `NEWTSEC` popiii_m2c.jl:1775 (ef988af); `NEWTDBG` popiii_m2c.jl:1829 (a8eb0ae);
  `NEWTFLIP` 0 popiii_m2c.jl:1836 (271e101); `JACCHK` popiii_m2c.jl:1673 (21cac65); `JACDUMP` popiii_m2c.jl:1719 (73491cc); `DAMPTRACE` 0 Trbdf2.jl:74 (a85f13b);
  `TSTPDBG` 0 TitanCore.jl:195 (5267f83); `STEPKX` 0 popiii_m2c.jl:2183 (225ef23).
- Radiation: `EDDDBG` popiii_m2c.jl:1755 (6fcf304); `ERDIAG` popiii_m2c.jl:2187 (f75ab23); `ERDIAGZ` popiii_m2c.jl:2192 (4b90c91); `FLUXTIEDIAG`
  TitanCore.jl:797 (a85f13b).
- Lines: `LLDIAG` 0 popiii_m2c.jl:544 (a41cbba); `LINECONS` 0 popiii_m2c.jl:552 (f6f2bad); `LINEINCDIAG` 0 popiii_m2c.jl:456 (b6fa7fd).
- Zones/profiles: `ZDIAG` popiii_m2c.jl:2059 (c559871); `SLZ` popiii_m2c.jl:2081 (4095a0a); `IFDBG`/`IFDBGN` 4 popiii_m2c.jl:2085/1986 (5267f83);
  `PROFAT` popiii_m2c.jl:2197 (root); `UDIAG` 0 popiii_m2c.jl:2218 (root); `FORCECHECK` 0 popiii_m2c.jl:282 (root); `REZDBG` 0 popiii_m2c.jl:875 (c81412a);
  `REZNULL` 0 popiii_m2c.jl:849 (f42916f); `RAYTEST` popiii_m2c.jl:379 (e96c96c).
- Audits: `TABAUDIT` popiii_m2c.jl:2026 (c3de545; values 1–4); `TABAUDNH` 17 popiii_m2c.jl:1907 (8807830).
- Timing: `DTDIAG` 0 popiii_m2c.jl:2175 (root), **on (1) in production**; `TIMING` 0 popiii_m2c.jl:1617 (f0bee17), **on (1) in production**.

---

## 14. Switches outside the production path (120; test, benchmark and old-stage scripts)

These are not read by `popiii_m2c.jl` or its includes. They have no bearing on the production run. Listed for completeness, by file.

- `s5_run.jl` (S5 isolated stars): ATM BURNT DBGTO DETAIL DTWHY DUMPNOCORE DUMPREJ FCOUPLING HSDISC ITERTEST JACTEST LAGRHO NANREPORT NOCONV NORULE PREDICT PROFILE_RELAX PROFILE_X RCAP REGIONCAP RELAX RESDUMP VSPLIT ZONEDUMP, plus shared with s4_d.jl / popiii_tests.jl / s5_b.jl / s4_b.jl: ALPHA GROW KCAP TEND STOPZAMS TRACE INT N NENV DT0 CONV EVERY MFLOOR M OUT PROFILE.
- `s4_d.jl`, `s4_b.jl`, `s4_a.jl`, `s4_star.jl`, `s5_b.jl`: CAPDIAG DTMAX ITMAX MIX MIXTEST NUCSTEP PROFILE_AT RSUN SMAXF GEOM NUC TMAX RHO NSTEPS S4DBG T0.
- `popiii_m1.jl`, `popiii_m2_onezone.jl`, `popiii_tests.jl`, `popiii_regimes.jl`: DECADES DUMPAT DUMPRHO FT Q RLIM XIMAX DLNRHO NMAX POPIII_PROJECT POPIII_SERIAL LI.
- Quasi-star line (`qstar_run.jl`, `QuasiStarSetup.jl`, `stepdiag.jl`, `newtrace.jl`, `infall_tests.jl`, `qstar_analyze.jl`, `test_run_convention.jl`): QS_SET QS_GASROW QS_CONV QS_COMPTON QS_TRBDF2 QS_MSG QS_META QS_SAVE QS_FREEZE QS_DTRIG QS_MLT_OLD QS_GRIDTAU QS_GRIDW QS_XSCALE SCAN RUN_CHECK.
- Phase/porting harness (`etotal_run.jl`, `phase3/5/6_run.jl`, `phase5_check.jl`, `newton_run.jl`, `energy_ledger.jl`, `order_test.jl`, `z0_jac_check.jl`, `av_test.jl`, `compton_box.jl`): TITAN_SET DN_SET ETOT_VERBOSE RUNTAG FOLLOW_JULIA PARITY TITAN_NCELL AMP CADV NLIST NREF SCHEMES AVC QL1P QL1S THET NOCOMPTON.
- Analysis and spectra (`closed_form.jl`, `line_profile.jl`, `nlte_spectrum.jl`, `lrd_metrics.jl`, `snapspec.jl`, `snapfull.jl`, `SpecFull.jl`, `ionfront_*.jl`, `hdchem_check.jl`, `stellareos_test.jl`, `tierc.jl`): CF_FUV CF_L CF_NLTE CF_STAR CF_TEFF LP_NLTE LP_NZ LP_TAG LP_XI LRD_R LRD_STAR SPECMAX HDLINES OUTDIR SPEC CKPATH EKPATH HELMDIR TIERC_REF.

---

## 15. Patch and workaround dossier: symptom, then the general replacement

Only switches categorised as physics patch, numerics workaround or time-step control are covered. The symptoms quote the
docstring or doc. "→" gives the candidate general replacement.

### 15.1 Table/network switching

- **PSWITCH=300 / PSWHYS=0.3** (c458b57). PSWHYS: "Warm gas of the accretion layer sitting on pswitch otherwise flips every
  step between table and network, which differ by ~10 % in T at fixed e" (PopIIIMaterial.jl:470-471). → **G2.** The ~10 % jump
  exists only because two material models meet. **R**
- **PSWTMIN=4500** (a8eb0ae, "Accretion on the fixed mesh: get past the 1.5 Msun stop"). "At 7000-9500 K the network's H
  ionization lags the table's Saha by ×5, and switching there drops T by up to 16 %" (PopIIIMaterial.jl:1033-1034). On
  12 Nov Tom set `PSWTMIN=0` for z100. In those runs the floor kept an H₂-dissociation zone on the network and made the
  failed solves (POPIII_RUNS.md:2049-2058). → **G2.** **R**
- **PSWTEX=5000 / PSWSTAY** (21cac65). "The network has no photoionization by the local radiation field, so at 8500 K it holds
  x(H+) well below Saha … T 8509 -> 13938 K in one step at the surface of a 28 M_sun protostar" (PopIIIMaterial.jl:476-478).
  → **G2 plus a per-zone J_ν (G4).** Then the network itself reaches Saha where J = B. **R**
- **PSWTAU=3** (ef988af, "get past the 29 M_sun slowdown and the 39.5 M_sun stop"). "The network has no radiation field: at
  the outer edge of an opaque layer … coronal x(H+) ~ 0.05 against 0.99, Newton heats the zone by a factor of 2"
  (PopIIIMaterial.jl:490-492). → **G2+G4.** **R**
- **PSWAGREE, PSWRAD, PSWRAMP, PSWNH** (all off). These are further entry/exit rules. PSWRAD covers zone 144 at 101.3 M☉,
  which "entered the table at T 7600 K with T_rad 5416 K … the step had no solution near the old state"
  (PopIIIMaterial.jl:545-547). PSWNH covers the ideal-network EOS departing from SCVH at n_H ≳ 1e17. → **G2.** PSWNH also
  needs a non-ideal network EOS (occupation probabilities, pressure ionisation). **R** for PSWAGREE, PSWRAD and PSWRAMP;
  **U** for PSWNH, because it marks real missing physics.
- **CHEMFBDENSE (default on)** (ce95c1b). This is the table fallback after a failed network Newton at p ≥ pswitch: "full3 at
  step 11935 … n_H 4.5e18, T 4346 K, Δt 515 s, twelve halvings without convergence" (POPIII_RUNS.md:45-52). → **G1+G6**:
  dissociation inside the hydro Newton. **R**
- **CHEMTMIN** (deddb60). A step returning T < floor counts as failed: "zone 69 … went 1516 K -> 0.3 K in the sub-stepped
  chemistry" (PopIIIMaterial.jl:942-946). → **G1/G6.** **R**

### 15.2 Opacity and star-light closures

- **NLOPAC** (72a27c2). The opacity is the LTE table × a lagged non-LTE ratio. Its follow-ups patch the lag:
  - **NLODAMP**: "flips between a hot opaque and a cool transparent state every 20–40 steps through that one-step lag"
    (PopIIIMaterial.jl:257-258).
  - **CPLNLOX**: "The lagged factors alone drove a two-step limit cycle (±30 % in T …, Δt held at ~2e3 s)"
    (PopIIIMaterial.jl:377-379).
  - **NLOIMPT** (on in production since 14 Nov 06:12): "OPAC's κ_P rises ×7 … while the continuum exponent q sees none of
    it … zones 81/82 at 312.9 M☉, Δt pinned at ~2e3 s" (PopIIIMaterial.jl:367-371; procedure.log 05:03 and 06:12).
  - **NLOFIX**: "χ 3e-9 cm²/g at 21 pc in c128 … against electron scattering on the zone's own x_e" (PopIIIMaterial.jl:409-411).
  - **NLOTMAX**: the ratio "came out 0.01–0.03, i.e. κ_R far below electron scattering" in ionised gas (PopIIIMaterial.jl:428-432).

  → **G2**: opacity computed from the zone's own populations at the Newton iterate, with no LTE base table and no lagged
  ratio. All six would go together. **R**
- **SLBLEND, SLEB, SLWFLUX** (on). The star light is closed as J_ν = W B_ν(T_col), with T_col and W derived from E_r − E_cmb
  and the τ = 2/3 zone.
  - SLEB: "the τ = 2/3 zone hopped between 122 and 130 … T_eff … swung 5580 ↔ 7100 K" (PopIIIMaterial.jl:525-529).
  - SLWFLUX: "E_r − E_cmb went negative in one zone and W = 0 there" (PopIIIMaterial.jl:784).
  - Several `max(er, Ecmb)` clamps (PopIIIMaterial.jl:773-799) are the hard-wired E_r ≈ E_cmb special cases.

  → **G4**: multigroup J_ν per zone, or at least a two-component field (star + CMB background as a boundary). The network
  would then read J_ν directly. This also covers **LOCALJ** and **NLOPACJ**. **R**

### 15.3 Operator splitting

- **CPLBAND=r=1500,1e6** (e982bb0). Only zones between 1500 R☉ and 1e6 R☉ are coupled; the rest use the split CHEMCK +
  CHEMCKHSS + ECHEMSPLIT. → **G1 in every zone.** Its satellites go with it:
  - **CPLBADDELTA**: "δ ln x = −1e76 in every retry" (TitanResidFast.jl:268-271).
  - **CPLCLAMP=3**: no docstring beyond "the clamp of a species correction per Newton iteration" (TitanResidFast.jl:217).
  - **CPLXSTOL=3**: "the guard against the hot, fully ionized root of an implicit step in the avalanche zones"
    (popiii_m2c.jl:1083-1084).
  - **CPLTHETA**: TR stage "x_γ ≈ −x_n: unreachable in the logarithmic species variables".

  CPLBADDELTA and CPLCLAMP are **R** under **G6**. CPLXSTOL is **R** under **G3**. CPLTHETA is **U**: an L-stable choice may
  be needed, so G3 should own it.
- **CPLIMPQ** (off). The three-level heats are booked lagged even in coupled zones: "with CPLNLOX its Δt → 0 limit was
  ~300 K below the split scheme's in the pre-shock band, unexplained" (TitanResidFast.jl:218-219). **U**: this is the open
  defect a full G1 must resolve.
- **ECHEMSPLIT** (on; fab1a39) and **XADV2** (on via SL3) and **XADVU**. "With x_e ~ 0.1–0.5 the first-order scheme carried a
  different amount of ionization energy across each interface than e did: up to ~1 eV per H a step" (popiii_m2c.jl:1028-1029).
  Also "zone 127 … T(e, x_new) = 13 800 K for the hydro's 5595 K" (popiii_m2c.jl:1063-1065). → **G1** with species advected
  inside the same Newton by the same interface scheme as e. **R**
- **BALSPLIT** (off, superseded): "first network zone outside the photosphere was heated 6227 → 8097 K in every hydro step
  … held Δt at the STOL limit (~2e3 s)" (TitanResidFast.jl:186-188). → **G1.** **R**
- **CHEMSUB / CHEMSUBMAX** (set to 1.3,3 but **inactive in phase B**, PopIIIMaterial.jl:1123). They counter the spurious
  hot root of a single backward-Euler step: "Newton ended on one at 2.6e6 K, fully ionized … 5.5e14 erg/g from nothing"
  (PopIIIMaterial.jl:854-857; POPIII_RUNS.md:458).
  - **CHEMCKHSS=100** (on): "electron avalanche e-folding in seconds within steps of 1e5–1e6 s … CHEMSUB tree … bisected to
    ~10 s (48 % of the chemistry time in 4 zones)" (ChemNet.jl:559-562).
  - **CHEMCKADAPT** and **CHEMSUBCYC** are alternatives for the same problem.

  → **G1** with a globalised Newton (G6) and error control (G3). The avalanche is then resolved by Δt control on species,
  not by sub-stepping a lagged split step. **R**
- **IMPCOOL / DTCOOL / TCOOLTH / LINEFRAC / ALI*** (line-cooling lag). Λ is held fixed or linearised (slope from a second
  transfer call, clamped to [−10, 30], popiii_m2c.jl:2170) over a step, and Δt is capped at 0.1 or 1 × t_cool. TCOOLTH: "the
  explicit Λ was used up to Δt ≈ 1.2 thermal cooling times … the core rings" (popiii_m2c.jl:238-241). → **G1** with Λ (or at
  least its EP part) in the Newton, plus **G3**. **R**. ALI refresh adaptivity is **U** because the ALI cost is real.

### 15.4 Chemistry Newton

- **CHEMROBUST, CHEMTRCLIP, CHEMSTAG, CHEMGIVEUP**, and the hard-wired stagnation acceptance. Symptoms:
  - "H⁺ ≈ 4e-12 … flips in its 10th digit" (POPIII_RUNS.md:96-100).
  - "b224 … 24 Newton iterations on average … many running to maxit = 60" (ChemNet.jl:517-519).
  - "~17 000 of ~56 000 attempts failed, each after 60 iterations" (ChemNet.jl:524-525).

  → **G6**, with scaled log variables (**CHEMLOGX**, which is **U**: it may itself be the general fix) and a residual-based
  convergence test. **R**

### 15.5 Radiation

- **FLUXTIE=0.1,1** (a85f13b). "With f → 1 the two moment equations decouple and r²E is pinned only by the outer boundary:
  F/cE reached 1.2 … and 6–25 (E_r → 1e-7 at step 87450)" (TitanDriver.jl:90-93). → **G4.** **R**

  FLUXTIE has knock-on patches:
  - **ETASWTIE=0.5**: total row minus KE minus radiation row "leaves −u × (moment-equation residual) on the gas"
    (TitanResidFast.jl:77-80).
  - **ETASWTH=1e-3**: "inflow against the gradient cools the gas: 0.12 K/yr at 2.6 pc … to 17 K and then 11.7 K, where the
    run stopped twice" (TitanResidFast.jl:69-73; POPIII_RUNS.md:901-925).
  - **TRESET** repairs zones cooled by this artefact.

  All four go together. **R**
- **ERTAUCONV** (default on, production off): "F/cE up to 440 seen … Δt → 1e-245 s" (TitanCore.jl:86-90). → **G4.** **R**
- **EDDGMU0** (on): "E 12 % below the CMB" because g stops at the last ray μ ≈ 0.33 (TitanCore.jl:98-100). This is a
  correctness fix of the outer boundary; it should become the default. **N**
- **EDDSCAT** (on): "chi emits at S ~ B: 50-1000 times the energy rows' emission" (TitanCore.jl:93-95). It treats all
  non-Planck extinction as scattering, which is a closure assumption. → **G2/G4**: a true absorption/scattering split from
  the material. **R**
- **EDDREAL, EDDCHI** (off): further Eddington-factor guards. **R**
- **ERTAUMIN=3000, ERTHICK, DTRADMAX**: zones are excluded from the E_r Δt rule, or Δt is capped for lagged Eddington factors
  ("steps much longer than the light-crossing time … fail in Newton (Eddington factors lagged over the step)",
  popiii_m2c.jl:736-737). → **G3** (error estimate excluding nothing) **+ G4.** **R**
- **EXCAP=1e8** (on): "c κ_E ρ Δt … can be 1e20 and the rows carry it as round-off" (TitanResidFast.jl:485-486). **U**: this
  is a conditioning device, and it conserves energy (the same factor in both rows). A scaled unknown (T_rad − T) or a
  one-temperature limit would remove it.

### 15.6 Energy formulation

- **ETOTAL=1 / ETASW=0.6** (phase B only). Phase A cannot use it: "two zones cooled in the hydro step from ~35 K to 17 K,
  below T_CMB" (POPIII_RUNS.md:33-38). The row also lacks gravity work on Lagrangian zones (**ETOTGRAVGRID**, off). → **G5**:
  one energy formulation valid on every mesh and in every phase, which removes ETASW, ETASWTH, ETASWTIE and ETOTGRAVGRID.
  The switch itself is **U**, because conservation is wanted.
- **ADVREL** (phase A only): "a zone that gains 1e-3 of its mass in the linearised step is heated by 40 %"
  (TitanResidFast.jl:129-130). This is a consistency fix. Keep it as the default. **N**

### 15.7 Hydro

- **AVTENSION=0** (on both phases). "With l = AV_ZONE·Δr and Δr ~ r … it can reach ~10× gravity and accelerate an
  infalling shell past the escape speed" (TitanResidFast.jl:319-321). **U**: an AV formulation choice that should be fixed
  once and not be a switch. A tensor/compression-only viscosity is the standard choice.
- **UPFREEZE, AVCAP** (off). **R**

### 15.8 Time step

- **CDYN** (0.01 in A1/A2, 1000 in A3, 3e5 in B), **DTENV** (removed), **STOL** 0.1/0.3, **NITER** 10/25, **NTRY** 10/60,
  **FAILTFAC 0.5**, **STOLCHEM**, **CPLXSTOL**, **DTCOOL**, **DTRADMAX**, **ERTAUMIN**. Per-phase hand-tuned values:
  - The A-phase solver controls differ because "with the accretion controls … 61 % of TR stages diverge from step 1"
    (POPIII_RUNS.md:29-32).
  - NTRY went 25 → 60 after the fx3 retry-limit stop (POPIII_RUNS.md:323-325).
  - CDYN 1e5 → 3e5 on 14 Nov gave 3× Δt, then a stall at 312.9 M☉ (procedure.log 04:30, 05:03).

  → **G3**: one local-error tolerance for the whole run, with a PI step controller. **R**

### 15.9 Mesh and staging

- The stage chain is RAD=0 → RADFROM rebuild → REZWIN core growth → ENVEUL remap, with MSTOP/NSTOP/TSTOP hand-overs,
  REFINE, REZONE*, ENVCONS, ENVMS0, ENVDT0, ENVSPLIT and ENVMERGE. → **G7**: one mesh from the IC, which is the stated
  next goal (POPIII_GOAL.md:3-10). **R**. RAD from the IC needs a radiation field built at t = 0; today
  "RAD=1 must never be set from the IC without the rebuild" (fullrun.sh).

---

## 16. Ranked: patches most worth retiring first

1. **The whole table/network switch** (PSWITCH, PSWHYS, PSWTMIN, PSWTEX, PSWSTAY, PSWTAU, the hard-wired PSWTAUDEEP, and the
   off-by-default PSWAGREE, PSWRAD and PSWRAMP). Six of them are active in production, and these rules caused most of the
   documented stops: 1.5, 28, 39.5 and 101.3 M☉, and the 312.9 M☉ stall "after zone 81 left the table" (procedure.log
   05:03). Every interface heuristic in §15.1–15.2 exists because two material models meet. Retiring it needs G2.
2. **The lagged non-LTE opacity stack** (NLOPAC, NLOIMPT, NLOFIX, NLOTMAX, CPLNLOX, NLODAMP). Three of these were added in
   the last 48 h, each for a new limit cycle or a stall. They are the purest symptom patches: opacity from the iterate's own
   populations removes all six.
3. **CPLBAND**: make COUPLED global. This removes the split CHEMCK path in the remaining zones, with CHEMCKHSS, ECHEMSPLIT,
   XADV2/XADVU, STOLCHEM and the inactive CHEMSUB. It is the single largest simplification, and COUPLED already exists.
4. **FLUXTIE + ETASWTIE + ETASWTH + TRESET.** One closure, plus the three patches for the gas-energy artefact it causes in
   transparent zones (the 117.6 M☉ stop).
5. **The time-step caps and per-phase solver controls** (CDYN, STOL, NITER, NTRY, FAILTFAC, CPLXSTOL, DTCOOL, DTRADMAX,
   ERTAUMIN). On 100 zones from the IC there cannot be per-stage tuning; G3 is a prerequisite for "one native run".
6. **The star-light J_ν closure** (SLBLEND, SLEB, SLWFLUX, the max(E_r, E_cmb) clamps), plus LOCALJ and NLOPACJ. These are
   the E_r ≈ E_cmb special cases; they need a per-zone spectrum (G4).
7. **Chemistry fallbacks** (CHEMFBDENSE, CHEMTMIN, the hard-wired `u < u(10 K, x)` error at PopIIIMaterial.jl:1119). These
   hide non-convergence by swapping the material model.
8. **Coupled-Newton guards** (CPLBADDELTA, CPLCLAMP; CHEMROBUST, CHEMTRCLIP, CHEMSTAG, CHEMGIVEUP; the 60-iteration
   stagnation acceptance). They become unnecessary with a globalised Newton (G6).
9. **The line-cooling lag** (IMPCOOL slope from a second transfer call, the [−10, 30] clamp, the 0.1 t_cool explicit cap,
   TCOOLTH, LINEFRAC). Put Λ_EP into the Newton.
10. **The dual-energy machinery** (ETOTAL/ETASW, ETOTGRAVGRID). Phase A and phase B use different energy rows today, which
    blocks a single run.
11. **EDDSCAT, EDDREAL, EDDCHI and ERTAUCONV.** These are Eddington-factor guards; only EDDGMU0 should survive, as the default.
12. **The stage and remap switches** (RADFROM, REFINE, REZWIN, REZONE*, ENVEUL remap, ENVCONS, ENVMS0, ENVDT0, MSTOP). They
    vanish by construction in the one-run goal (G7).
13. **EXCAP**. It is low risk but hides conditioning, so reformulate rather than cap.

Make these the defaults rather than retiring them: **TL3RECKT** (an LTE-consistency fix, currently off), **EDDGMU0**,
**ADVREL**, **LYABB** and **EPGRAD=thermal** (already the default). They are correctness fixes, not patches.

## 17. Groups that one general mechanism replaces together

| General mechanism | Switches it retires |
|---|---|
| **G2: one non-LTE material** (EOS and opacity from the network populations everywhere; non-ideal terms added to the network EOS) | PSWITCH PSWHYS PSWTMIN PSWTEX PSWSTAY PSWTAU PSWAGREE PSWRAD PSWRAMP (PSWNH needs non-ideal EOS) PSWTAUDEEP; NLOPAC NLOIMPT NLOFIX NLOTMAX NLODAMP CPLNLOX; CHEMFBDENSE CHEMTMIN; RHOEQ; HOTT1/HOTT2 blend (U); EDDSCAT |
| **G1: fully coupled chemistry** (COUPLED in every zone, every species, three-level heats implicit, Λ_lines in Newton, species advected in Newton) | CPLBAND CPLIMPQ (→ on) CPLTHETA (→ G3) CPLLOCAL; ECHEMSPLIT XADV2 XADVU BALSPLIT; CHEMSUB CHEMSUBMAX CHEMCKHSS CHEMCKADAPT CHEMSUBCYC CHEMSUBMAX_SC; STOLCHEM; IMPCOOL IMPALWAYS LINEFRAC TCOOLTH; SPLITTOL SPLITIT ZSPLITTOL DXMAX EFRAC (nuclear sources likewise) |
| **G3: error-controlled time stepping** | CDYN DTENV STOL NITER NTRY FAILTFAC CPLXSTOL DTCOOL DTRADMAX ERTAUMIN ERTHICK ENVDT0 S5NUM S5TOL; the 0.1 t_cool cap; per-stage A/B controls |
| **G4: proper radiation transport and boundary** (multigroup or a consistent closure; CMB only as the incoming boundary/background; per-zone J_ν for the network) | FLUXTIE ETASWTIE ETASWTH TRESET; ERTAUCONV EDDREAL EDDTAU EDDRELAX EDDCHI (EDDSCAT with G2); SLBLEND SLEB SLWFLUX LOCALJ NLOPACJ; the max(E_r, E_cmb) clamps and the TCMB `cool_floor` subtraction |
| **G5: one conservative energy formulation** | ETOTAL ETASW ETASWTH ETOTGRAVGRID (with G4: ETASWTIE TRESET); EXCAP (U) |
| **G6: globalised Newton** | CPLBADDELTA CPLCLAMP; CHEMROBUST CHEMTRCLIP CHEMSTAG CHEMGIVEUP (CHEMLOGX possibly part of it); DAMPNEWTON DAMPALPHA DAMPTRIG; UPFREEZE; the 60-iteration stagnation acceptance |
| **G7: one mesh from the IC** (fixed log mesh / GRIDSEG, not yet in code) | RADFROM REFINE REZONE REZEVERY REZJSPLIT REZJMERGE REZDMMIN REZQMAX REZNMAX REZFRAC REZKEEPC REZWIN REZDR REZPSPLIT REZSLOPE; ENVEUL-as-remap ENVCONS ENVMS0 ENVSPLIT ENVMERGE ENVSINKTAU; MSTOP NSTOP TSTOP as hand-overs; MAXSTEPS segmenting |

Open questions I could not settle from the code and docs:

- With COUPLED=1, the docstring says it "excludes ECHEMSPLIT" (popiii_m2c.jl:1082). Production runs both, which the
  CPLBAND path allows (popiii_m2c.jl:2035).
- ADVREL is unset in phase B. Whether it matters on the fixed mesh (lgrid=2, no frozen XMINT) is not stated.
- EQEOS, Z0_TABLE, SHOCKRAD, GRIDWIN and GRIDMAP are off in production; I have no evidence on whether they are still wanted.
