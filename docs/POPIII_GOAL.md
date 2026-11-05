# Goal: a first star to 300 M☉ or to the end of its own accretion

Set 2026-11-04. Measurements below are from run `_c2` unless stated otherwise.

## The goal

One 1D spherical calculation, continued from the collapse of the z = 19 metal-free cloud (`IC=abn02z19`), on one fixed
logarithmic mesh from the centre through the star, the shock and the whole infall. It ends at the first of:

| End | Measured how |
|---|---|
| **A. 300 M☉** | mass inside the accretion shock reaches 300 M☉ |
| **B. Accretion halted** | inflow through the shock stays low enough that M\*/Ṁ exceeds the star's remaining life (about 2 Myr), or the gas outside the star moves outward at all radii up to the ionization front or the Eddington-driven surface, for longer than 10 % of the star's age |
| **C. A documented failure** | a limit we can name and could not get past, with the mass reached and why |

A and B are the science results. C counts as an honest outcome only when the cause is identified.

The result to report is M\*(t), R\*(t), L(t), Ṁ(t), the end mass and the process that set it, plus the velocity, density and
ionization of the infall at each stage, because that is what later sets the observable line emission.

## What it will not be

It is the spherical answer, not the mass of real first stars. There is no rotation, disk, fragmentation, magnetic field or 3D shadowing. Ionizing light
reaches all directions, which makes feedback stronger than behind a disk's shadow, while disk photoevaporation (the end McKee & Tan 2008 find) cannot happen. It is
one cloud, so it says nothing about the spread between clouds.

## Is it realistic?

**Reservoir.** Enclosed mass and free-fall time at step 30580 (a rough forecast, since real arrival depends on how the infall evolves):

| M [M☉] | r [R☉] | t_ff [yr] | Ṁ between rows [M☉/yr] |
|---|---|---|---|
| 1.5 | 1.1e3 | 1.7 | |
| 10 | 7.6e4 | 360 | 0.025 |
| 50 | 4.6e5 | 2.5e3 | 0.019 |
| 106 | 1.25e6 | 7.6e3 | 0.011 |
| 155 | 2.3e6 | 1.6e4 | 0.006 |
| 206 | 3.8e6 | 2.9e4 | 0.004 |
| 312 | 8.4e6 | 7.7e4 | 0.002 |
| 513 | 2.3e7 | 2.7e5 | 0.001 |

So 300 M☉ is available in about 10⁵ yr, far less than the stellar lifetime.

**What published work leads us to expect.** Omukai & Palla 2003 find Ṁ_crit ≈ 4×10⁻³ M☉/yr: below it, growth continues to about 300 M☉; above it, the star
reaches the Eddington limit before CN ignition, at about 90 M☉ in the case they discuss. Our forecast Ṁ stays above 4×10⁻³ until about 200 M☉, which covers the
mass range where they find trouble, so **ending near 100–200 M☉ through end B is a real possibility**, which is why the goal allows it. Hosokawa & Omukai 2009 describe adiabatic
accretion, swelling, contraction and main sequence phases at lower masses for their accretion rates; ZAMS arrival at about 40 M☉ at their fiducial rate.
Omukai & Inutsuka 2002 and McKee & Tan 2008 (`popiii_lit/feedback`) give the ionization conditions. These are the comparisons, not targets.

**Time steps.** From 0.66 to 1.47 M☉, 0.8 e-folds of mass took 4200 steps, about 5300 per e-fold. From 1.5 to 300 M☉ is 5.3 e-folds, so about 3×10⁴ steps at 2–3 s each
is **about a day of wall time if the steps keep scaling with M/Ṁ**. Steps must reach about 10⁸ s from today's 10⁶ s. This is the optimistic case, and is not yet demonstrated: `_c3`, on the finer mesh,
averaged 640 s per step and would take about 100× longer. A run needing more than about two weeks fails the budget gate below.

## What stands between here and there

Milestones, in order, each with its test. None is done unless its test passes.

| # | Milestone | Needs | Test |
|---|---|---|---|
| M0 | Past the stop at 1.5 M☉ | diagnose why `_c1` and `_c2` both fail at zone 124, near the shock | reproduces, cause named, passes without changing the star elsewhere |
| M1 | 10 M☉, adiabatic and swelling phases | steps growing with M/Ṁ; mesh check | R\*(M), Ṁ against Hosokawa & Omukai 2009 at similar Ṁ; two mesh spacings agree; energy budget closes |
| M2 | Contraction and ignition, about 30–100 M☉ | nuclear abundances carried on the fixed mesh; star shrinking through the mesh | L, T_c, R at arrival against Omukai & Palla 2003 and our 15 and 100 M☉ stars |
| M3 | Ionizing light and the bubble | **the only new solver**: a photon-conserving ray from the star, coupled to the existing chemistry and heating, since one gray field cannot carry it | Strömgren radius and R-type front (Iliev et al. 2006 test 1) |
| M4 | 300 M☉ or halt | M0–M3 together; Eddington-driven outflow, Lyman-α pressure checked at least by estimate (McKee & Tan 2008) | end A or B reached; mesh and restart repeatability |

Biggest risks, as I judge them today: (1) steps stay near the fixed-mesh shock crossing time rather than growing; (2) the contracting star falls below the mesh resolution or the shock front cannot follow it; (3) the coupled ionization front is stiff across a thin zone.

Notes for M3 from reading the code (4 Nov): the network's photo-rates already integrate σ_ν J_ν through `env.J` (`ChemNet.photorates`), but every zone gets the same `env`, J = B(T_CMB). So the chemistry does not yet see the local radiation field at all, including the star's light in the precursor. M3 can therefore be mostly a per-zone J_ν (the star's attenuated spectrum), rather than a new chemistry solver. The H⁻ and H rates in the precursor also need this, and that may matter before M3.

M2 note: nuclear abundances now move with the same mass fluxes as the chemistry (`advect_species!`), untested until hot zones appear.

M3 state (4 Nov). `julia/IonRay.jl`: multi-group radial ray, photon-conserving zone by zone, solved implicitly inside `ChemNet.step`, absorbers H and H₂ (AAZN97 cross section). What has passed and what has not:

| Test | Result | Status |
|---|---|---|
| Iliev et al. 2006 test 1 (`ionfront_test.jl`; static, T fixed) | r_I/analytic with α_B = 2.59e-13: 1.00 at 10 Myr, 1.06 at 500 Myr; the network's own α_B is 2.48e-13 (n ≤ 15), giving 1.04. Unchanged for Δt = 0.5–10 Myr. The paper's codes stay within 5 %. Ionizations + recombinations exceed the photons by 1–2 %, not traced (collisional ionization is the candidate). | passed |
| Iliev et al. 2009 test 5 (`ionfront_hydro.jl`; hydro on TITAN's fixed mesh, blackbody, no gray field) | r_I/r_S0 = 0.87 at 100 Myr, 1.12 at 200, 1.67 at 500, against 0.88, 1.15, about 1.55 (range about 1.5–1.65) read off their Fig. 18. r_I unchanged at 256 zones and half the Δt (9.02 against 9.01 kpc at 500 Myr), so the excess is not mesh error. Mass conserved to 1e-13. Shell peak n 1.5e-3 at r/L = 0.80 against ZEUS-MP's 2.5e-3 at 0.78 at the same Δr (128 zones). | close, slightly fast at 500 Myr, converged in mesh |
| Coupled to the gray field (`RAY=1`: `RAYQ` removes the energy from the gray rows) | runs at 1e44 and 1e46 photons/s on the real state at 3 M☉, Δt falling with the rate. 1e47 and above stop at the first step: test light of 1e5 K on gas that is still 6000 K and opaque, which is not a state the star reaches. | not validated |
| HII temperature with cooling through the gray field | not done | open |

Not in the ray yet: helium (fixed neutral at y_He, so no absorption by He and no He⁺ heating), Lyman-Werner dissociation of H₂, and the diffuse field. Gray emission uses Planck means at LTE composition, so it includes recombination to the ground state, while the network is case B in lit zones.

Bug found by test 5: when `ChemNet.step` halved a step, the second half restarted from the original u, dropping what cooling or heating did in the first half. Fixed. It touched any network zone with non-zero cooling (H₂ lines) whose Newton failed at full Δt, so runs before `_h5` have this energy error at those zones. `_h4` (old) against `_h5` (fixed) from one checkpoint measures its size; `_h6`, `_h7` restart the mesh check with the fix.

Size of that bug: `_h4` against `_h5` from step 35500 (3.5 M☉) to 5.2 M☉ differ by under 0.2 % in M, up to about 1 % in L, one zone in R.

## Numerics found on the way to M1 (4 Nov)

Runs were not gaining from longer accretion times: Δt followed a sawtooth (×2.5 per step up to 2e6 s, then ÷50–100), mean 5e5 s, which projected to months. Causes, each measured from the checkpoint at step 40000 (5.2 M☉, 400–800 step tests):

| Cause | Measured | Change |
|---|---|---|
| A failed attempt cut Δt by 20 | NITER=25 alone no gain; FAILTFAC=0.5 doubles mean Δt | `FAILTFAC=0.5`, `NITER=25` |
| Round-off from the gas-radiation exchange: c κ ρ Δt reaches 1e20 in the opaque core, so Newton stalled near 1e-4 in zones 3–12, growing with Δt. The upwind flip was the first suspect; freezing the branches changed nothing (95 failures against 94) | capping that product at 1e8: 0 failures in 400 steps (from 95). 1e4 is too low, stalling at the shock | `EXCAP=1e8` (same factor in both energy rows, so energy is conserved; E = aT⁴ holds to about 1e-8) |
| `CDYN=1000`, 1000 × central dynamical time, capped Δt at 2.8e6 s | no longer binds at 1e4 or 1e5 | `CDYN=1e5` |
| H₂ line transfer every step, 75 % of the cost | every 0.1 t_cool: M and infall ρ agree within 6e-4, cost halved at this stage, and falling as Δt grows | `LINEFRAC=0.1` |

All together against the old numerics, `_h8` against `_h5` from 5.2 to 7 M☉: M within 0.4 %, about what the snapshot time offset gives, and L within a few %, inside the star's fluctuations. Δt then reaches 1e7–6e7 s by 10–18 M☉ (`_h9`, `_h8`), set by E_r in the zones under the photosphere. `ERTHICK=10` removes that repeat of the T criterion, for ×1.4 net at 5 M☉; left off for now.

**The energy budget did not close, and this is a physical matter, not bookkeeping.** `popiii_track.jl` books a fixed sphere around the star against what crosses its face. With the gas energy row, energy appeared at the accretion shock: +7 % of the radiated energy at 6–7 M☉ on the 0.03 mesh, +24 % on the 0.06 mesh (growing to +33 % at 14 M☉), all between 150 and 240 R☉. So L above GMṀ/R in those runs, and the swelling they show, were at least partly numerical; `_h9` reached L/L_acc = 1.8 at 15 M☉, and that number should not be used.

Fix, with TITAN's own total-energy row (P4.2, `ETOTAL=1`, gas row kept in the hypersonic infall, `ETASW=0.1`). It had no gravity work, as it had never run with self-gravity. Work computed from the momentum row's force g m_k still left +5 %, because m_k, the mean of the neighbouring zones, is not the mass crossing the face where ρ is steep. Computing it from the continuity row's mass flux (`ETOT_GRAV`) conserves gas + kinetic + potential energy by construction: imbalance now −3 % (0.03 mesh) and −1 % (0.06 mesh) at 1000 R☉, from an 800-step test after a restart. Longer windows are to come from the runs below. The whole-run test restarts from 0.5 M☉, because earlier stars carry the shock error in their history: `_i3` (0.03), `_i6` (0.06), `_i1` (0.015).

**The dual-energy switch and the surface (evening, 4 Nov).** With the total-energy row, η = e/(e + K) was singular for a material whose e has a negative zero point: a zone reached 18 K at 0.55 M☉ (`_i6`, first attempt). The test now uses e + UOFF_TH. At ETASW = 0.1 one surface zone at 15 M☉ then fell to e + U < 0 (1700 K, x_H2 0.45) and back, and L jumped between 5e4 and 1.1e5 L☉. At ETASW = 0.3 L stays near 5.5e4 L☉, so 0.3 replaces 0.1. R and M up to 10 M☉ agree between the two within about 3 %.

From about 15 M☉ the outer layers oscillate at ±2–5 km/s with a period near 1.5e7 s, and E_r in the surface zone then limits Δt to about 1e5 s. Whether this is physical has not been decided; finer meshes will show. `STOL=0.3` (from 0.1): mean Δt about 8e6 s over 15.4–25 M☉; R within 1 %, L within 2 % of the STOL = 0.1 run at 15.5 and 15.9 M☉ (`_l1` against `_j6`); budget closes to +2–4 % at 500–3000 R☉ over 15.4–25 M☉. STOL = 1 fails to converge. `NTRY=25` with `FAILTFAC=0.5` allows Δt to fall by 3e-8 before a stop; ten tries stopped `_k3` at 12.2 M☉.

**The switch pins the zone it sits on (night of 4 Nov).** At ETASW = 0.3 `_p6` stalled at 24.6 M☉ with Δt about 3e3 s. One zone held T ≈ 2300 K (radiation temperature 6800 K) with e flat on average. A term-by-term look showed de/dt taking only two values, +7e6 erg/g/s on the gas row and −2.3e7 on the total row, flipping as η crossed ETASW, so the switch holds the zone at η = ETASW (the same pinning gave 1700 K at 0.1). The cold state then sits on the H₂ dissociation plateau, where the Newton's frozen-composition T overshoots by 20–30 % each step and chem! brings it back, and that is what limited Δt. ETASW = 0.6 puts the switch where gas near the radiation temperature keeps η below it: from the stalled state 24.6→26.4 M☉ in 20 min, budget within 3 %. ETASW = 0.8 is slower, budget 7–17 %. Rerun from 0.5 M☉ at 0.6 (`_f6`): 25 M☉ in 45 min of wall time, budget within 1 % at 500–3000 R☉ over the whole run. That the row budget at low mass is not worse at 0.6 than at 0.1 is measured at this mesh only.

## Rules for this work

- The infall stays on a fixed logarithmic mesh, and no part of it is dropped.
- Native TITAN solvers first; new machinery only where they cannot do the job, and M3 is the only such case known.
- Do not quote a number as a result until its test above passes; say what is measured and what is inferred.
- **Budget gate:** if the projected wall time of the full-mesh run exceeds about two weeks after M1, stop and bring the choice to Tom: a star plus a separate infall calculation with a radiating inner boundary (option B of 3 Nov), which would keep the fixed mesh for the infall.
- Other decisions to bring back rather than make alone: dropping any physics from the end state, changing the goal, anything public.

## Progress

| Date | State |
|---|---|
| 2026-11-04 | Goal set. Best run 1.49 M☉ (`_c1`), 1.47 M☉ (`_c2`, with convection). M0 in progress. |
| 2026-11-04 | **M0 passed.** Cause: when a zone's pressure crosses PSWITCH it changes from network to table at fixed e. At 7000–9500 K the network's H ionization is ×5 below Saha, so T drops up to 16 % (9520 → 7977 K) and κ falls ×9. E_r must relax in ~7 ms, Newton fails down to that Δt, and DTSTOP = 1e-2 stopped the run. Reproduced at the same step. Moved the switch to PSWITCH = 300 with T ≥ 4500 K (`PSWTMIN`), where the mismatch is ~1 %; M, T_c and L agree with the old switch to 0.04 % at t = 4.5 yr (restart comparison). Long run `_h1` started. Whether the network's low ionization or Saha is right in the switching gas has not been decided: the zones there are optically thick, which favours Saha. |
| 2026-11-04 | `_h1` reached 3.2 M☉ with no stop, about 5800 steps per e-fold, 1.8 s per step, 78 % of it H₂ lines. M3 ray passes Iliev test 1; hydro test 5 close. |
| 2026-11-04 | Numerics above. On the gas row, R(shock) against the adiabatic scaling (Hosokawa & Omukai 2009 eq. 9): ratio 0.66 at 0.7 M☉, 0.81 at 1.5, 0.92 at 3, 0.98 at 5 M☉, at Ṁ about 0.035–0.05 M☉/yr; to be repeated on `_i3`. Interior T peaks off centre (3.5e5 K at m = 0.55 M☉, 5 M☉), near the virial estimate. Deuterium is not in the network; Hosokawa & Omukai find it unimportant at such rates. |
| 2026-11-04 | Flags now `NITER=25 FAILTFAC=0.5 LINEFRAC=0.1 EXCAP=1e8 CDYN=1e5 ETOTAL=1 ETASW=0.3 STOL=0.3 NTRY=25`. R(shock) against eq. 9 at Ṁ 0.05→0.028: 0.06 mesh (`_j6`) 0.72 at 0.8 M☉, 0.99 at 4, 1.08 at 10; 0.03 mesh (`_i3`, ETASW 0.1) 0.67, 0.93, 0.97. The meshes differ by 2–11 %, up to two zones of the coarser one; 0.015 mesh still early. Budget −3 to −5 % of radiated to 12 M☉ on the 0.03 mesh. Swelling starts near 12 M☉: 0.06 mesh R/eq. 9 = 1.30 at 15 M☉, 1.45 at 20.6, L/(GMṀ/R) 2.7 at 20.6; R ≈ 400 R☉, L ≈ 2e5 L☉ at 26.6 M☉ (scout `_l1`, mesh not yet checked). Production `_p6`; 0.03 and 0.015 restarts in progress. |
| 2026-11-05 | ETASW = 0.6 in all runs. `_f6` (0.06, from 0.5 M☉): R/eq. 9 = 0.63 at 1 M☉, 0.93 at 4, 1.07 at 9.3, 1.33 at 15, 1.61 at 19.5 (333 R☉); L/(GMṀ/R) 2.6 at 19.5. `_e6` stopped at 27.8 M☉ after one 2.2e8 s step left a zone unsolvable; `CDYN=1e4` caps Δt near 3e7 s (`_g6`). 0.03 and 0.015 meshes: `_r3` (17 M☉), `_r1` (2.8 M☉). |
| 2026-11-05 | `_f6`, `_g6`, `_e6` all stopped at zone 118 (28 M☉) with corrections that did not shrink with Δt. Cause: a table zone at 8509 K fell below `PSWITCH`; the network holds x(H⁺) below Saha at fixed e, so T jumped to 13938 K. `PSWTEX=5000` keeps table zones on the table until T < 5000 K (E_ion/E_th ≤ 4 % there). `_s6` passes 28 M☉. New limit: from 29 M☉, Δt ~ 4e4 s against 2e7 s before, about 0.05 M☉/hr. The surface zone runs a 1e6 s sawtooth (T 11000–27000 K at constant ρ), and Newton converges linearly (ratio ~0.7) and fails when Δt grows. Ruled out: Δt limiter (`ERTHICK=10` made no difference), κ derivatives, Jacobian columns of zone 126 (`JACCHK`), linear solve. Mechanism still open. `_r3` runs at 0.3–0.4 M☉/hr at 17.7 M☉, past the two-week gate; brought to Tom. |
| 2026-11-05 | The slowdown at 29 M☉ (Δt ~ 3e4 s, ~0.05 M☉/hr) had three causes, each tested from one checkpoint. (1) The E_r floor in `maxdel!` makes Newton fail: in thin zones the F row sets E and the E row sets F, so forcing E = F/c fights both (44 zones on the floor per iteration). `ERTAUCONV=0`: 0 divergences, was 896. (2) One zone holds the whole photosphere (dτ = 120, T 25585→8544 K, Δln r = 0.06 is ~9 scale heights); its E_r swings ±40 % and sets Δt. `ERTAUMIN=3000` takes it out of the rule (T stays in) and Δt returns to 3e7–1e8 s. (3) At p < PSWITCH the outer edge of the opaque layer stays on the network, which has no radiation field and sits at coronal x(H⁺) ≈ 0.05 against Saha 0.99, so it is not a lag; Newton stopped at 39.5 M☉. `PSWTAU=1` uses the table at τ ≥ 1. New flags: add `ERTAUCONV=0 ERTAUMIN=3000 PSWTAU=1`. `_j7` (0.06 mesh) reaches 46.6 M☉, R = 84.5 R☉, Tmax 6.5e6 K, about 35 M☉/hr, 0 divergences. **Open:** with the floor off, F/cE reaches 2.9 in thin zones (f ≈ 0.65 at 1000–25000 R☉ violates f ≥ (F/cE)²), so E_r and T ahead of the photosphere are not trustworthy; enforcing the bound (`EDDREAL=1`) collapses Δt, not understood. The fast contraction is one mesh only; `_k3`, `_k1` run the same flags. |
| 2026-11-05 | **Correction to the row above for M > 75 M☉.** `_j7` stopped at 117 M☉ with T_c = 4e8 K. Not ignition: the equilibrium table began at log ρ = −12 and `bicubic` extrapolated a cubic in log p below that; `PSWTEX`/`PSWTAU` put zones there. p was wrong by 13 % at 80 M☉, 7× at 90, up to about 75× by 113 M☉, next to the photosphere, then ran away in zone 133 (5672 dyn/cm² against 1), and Δt fell to 0.03 s. `_j7` envelope and photosphere above ~80 M☉ are invalid; up to ~75 M☉ no table zone was below the edge. Table rebuilt to log ρ = −20 (`eqeos_table_lr20.jls`; the old file is left in place), and `eos` now errors below it. Step-size test from 36 M☉: Δt ~ 40× smaller gives the same Tmax to 0.03 % at 37 M☉. For M < 75 M☉ the structure is a contracting core inside an opaque accreting envelope with the photosphere in the flow, as Hosokawa et al. 2012 describe for 6e-3 M☉/yr. `_a8`, `_a3`, `_a1` rerun with the fixed table. |
