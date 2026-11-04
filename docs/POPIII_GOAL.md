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
| Iliev et al. 2009 test 5 (`ionfront_hydro.jl`; hydro on TITAN's fixed mesh, blackbody, no gray field) | r_I/r_S0 = 0.87 at 100 Myr, 1.12 at 200, 1.67 at 500, against 0.88, 1.15, about 1.55 (range about 1.5–1.65) read off their Fig. 18. Mass conserved to 1e-14. Shell peak n 1.5e-3 at r/L = 0.80 against ZEUS-MP's 2.5e-3 at 0.78 at the same Δr. | close, slightly fast at 500 Myr; finer mesh pending |
| Coupled to the gray field (`RAY=1`: `RAYQ` removes the energy from the gray rows) | runs at 1e44 and 1e46 photons/s on the real state at 3 M☉, Δt falling with the rate. 1e47 and above stop at the first step: test light of 1e5 K on gas that is still 6000 K and opaque, which is not a state the star reaches. | not validated |
| HII temperature with cooling through the gray field | not done | open |

Not in the ray yet: helium (fixed neutral at y_He, so no absorption by He and no He⁺ heating), Lyman-Werner dissociation of H₂, and the diffuse field. Gray emission uses Planck means at LTE composition, so it includes recombination to the ground state, while the network is case B in lit zones.

Bug found by test 5: when `ChemNet.step` halved a step, the second half restarted from the original u, dropping what cooling or heating did in the first half. Fixed. It touched any network zone with non-zero cooling (H₂ lines) whose Newton failed at full Δt, so runs before `_h5` have this energy error at those zones. `_h4` (old) against `_h5` (fixed) from one checkpoint measures its size; `_h6`, `_h7` restart the mesh check with the fix.

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
