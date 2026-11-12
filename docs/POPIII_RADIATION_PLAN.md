# One radiation physics for the run and its spectrum (plan, 12 Nov 2026)

Decision (Tom, 12 Nov): all shortcomings of the emission physics are to be fixed now, before resolved runs, and the emergent
spectrum, including helium, must be the one the simulation itself uses. Same atoms, same rates, same opacities and
emissivities in the evolution and in the snapshot analysis.

## What is wrong today

| | Run | Snapshot spectrum (`snapspec.jl`) | Older line tools (`NLTEH.jl` …) |
|---|---|---|---|
| H atom | 3-level (n = 1, 2, continuum), Sobolev Lyα escape, in network zones; Saha in table zones | n = 2 approximated (min of Saha, Boltzmann) | 12-level non-LTE, Sobolev lines; only infall below 2e4 K |
| Helium | inert, neutral in network and in low-p EOS; ionized in the opacity table | no opacity | LTE (He I 10830: 3-state post-process) |
| Lines | H₂ (63 levels, ALI); Lyα local | H₂ only | Balmer, Hα/Hβ, He I 10830, O I Bowen |
| Ionizing light | none (`IonRay` off) | continuum only | n = 1 photoionization, X-ray secondaries |
| Velocity | — | static | Sobolev in free fall |
| Consistency | grey means from table/three-level opacities | ∫L_ν = 0.8 × run's L | separate input format |

At 180 M☉ the envelope is 24 M☉ of H II/He II/He III at 35–160 kK under a 3400 K photosphere: exactly where all three columns are
incomplete.

## Architecture: one spectral solver, refreshed in the run, rerun on snapshots

The run keeps its implicit grey two-moment transport for the hydrodynamics. Everything frequency- and level-dependent comes from
one module (working name `RadSpec`). It is called in the run every few steps, the same way the H₂ ALI is refreshed now, and in
the analysis on the recorded snapshot.

1. **Atoms.** H I n ≤ 12 (n = 1, 2 resolved where needed); He I singlets and triplets (1¹S, 2³S, 2¹S, 2³P, 2¹P, n = 3–5 collapsed);
   He II hydrogenic n ≤ 8; H⁻; H₂ through the existing `H2Levels`/`H2Lines`. Level-resolved recombination, photoionization,
   collisional excitation and ionization, two-photon decay. Data from refereed sources only, verified on ADS.
2. **Chemistry and EOS.** He⁰/He⁺/He²⁺ become network species (charge, electrons, energy). The table EOS gets helium Saha in its
   low-pressure part, matching the opacity table. One energy zero across the network/table switch.
3. **Transfer.** Continuum plus lines on one frequency grid, spherical tangent rays, ALI. Lines are Sobolev in the run's own
   velocity field. Lyα gets escape and destruction probabilities in the run, calibrated against a Monte Carlo transfer with
   partial redistribution in the analysis. The ionizing continuum (`IonRay`) becomes the radial part of the same grid,
   with He I and He II absorption.
4. **Feedback into the run, per zone, at each refresh.** Net heating and cooling by every atomic process; photo-rates for the
   network (H, He, H⁻, H₂, LW from level populations); κ_P, κ_E, κ_R and the **flux mean** from the actual J_ν (this also gives
   the non-grey radiation force); line force (Lyα and the rest); Eddington factors as now. Between refreshes the departure
   coefficients b_n are held and the rates rescaled with ρ, T, as is done for the H₂ lines.
5. **Records.** Snapshots store b_n of every level, the ionization state and the refresh's J_ν moments. The analysis rebuilds the
   same state, then adds what only post-processing needs (fine frequency grid, line profiles, observer frame).

**Consistency tests:** (i) the analysis solver on a snapshot, at run resolution, reproduces the run's per-zone heating, cooling and
rates to round-off; (ii) ∫L_ν dν equals the run's L within 2 % (today 0.8). If (ii) cannot be met with grey hydrodynamics plus
refreshed means, the fallback is multigroup transport in the run. **That fallback would be a method change and comes back to Tom.**

## Order of work

| Step | Content | Test before it counts |
|---|---|---|
| R0 | Inventory of the existing pieces; atomic data with refereed sources | data reproduce published case-B emissivities (Hβ, He II 4686, He I 5876) and recombination coefficients |
| R1 | Helium in the network and the low-p EOS (= P1) | Saha limit; energy conservation through He ionization; restart from `c128` at 180 M☉ against the unchanged run |
| R2 | `RadSpec` non-LTE H + He + lines with velocity, on snapshots | case-B limits in a static nebula; Iliev test 1 with He; earlier `NLTEH` results reproduced where they overlap |
| R3 | Coupling into the run (default off), refresh like the H₂ ALI | bit-identical when off; consistency test (i); restarts at 4, 40 and 180 M☉ |
| R4 | Lyα: run escape/destruction + force; Monte Carlo for the spectrum | MC against analytic slab solutions (Neufeld/Harrington); run vs. MC destruction rates |
| R5 | Ionizing ray in the same grid, He photoionization, momentum | Iliev tests 1 and 5 with the coupling on; H II temperature |
| R6 | Records and observables export (spectra panel in the explorer) | consistency test (ii) |

Cost target: a refresh no more expensive than today's H₂ ALI call (GPU kernel path through KernelAbstractions as for the H₂
lines, CUDA and Metal).

## Model split

- Opus: design, solver and coupling code, interpretation, decisions for Tom.
- Sonnet subagents: atomic data compilation and ADS verification, unit tests against published tables, running the analysis
  on snapshots, comparisons against written criteria, first-pass crash triage.
- Haiku / shell watchdog: routine monitoring.

The running chains (`c128`, `f128`, `b224`) continue unchanged as the scaffold; each new piece is tested by restarts from
their frozen checkpoints.
