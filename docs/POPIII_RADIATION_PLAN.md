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

## What exists to build on (inventory, 12 Nov)

| Component | Build on | Gap |
|---|---|---|
| H atom | `NLTEH.jl` (n ≤ 12, Sobolev β, `solve_point`); exact f-values, bound-free and free-free Gaunt factors (van Hoof 2014), occupation probabilities and Stark profiles in `z0/atomic.jl` | crude collision rates (van Regemorter, Seaton); dense Λ-iteration; old snapshot format; no case-B check |
| He I / He II | `z0/` (Verner ground-state cross sections, NIST He I multiplets, hydrogenic He II n ≤ 16), `he10830.jl` (3-state, approximate) | no He I level model, no He collision or recombination data, no He II statistical equilibrium |
| Chemistry | `ChemNet` machinery (species, `Photo` with Milne inverse, element elimination already done for H and Li) | He hard-wired inert in `full`/`fullvec`; `fill_opac!` takes electrons from H⁺ only |
| EOS | `Z0.saha` (full H + He Saha) | `PopIIIEqEOS.node` has neutral He at all T, p; energy zero across the switch unchecked; **Y = 0.25 in the run (yHe = 1/12) against 0.24 in Z0** |
| Continuum transfer | `Spectrum.jl` (rays, Feautrier, ALI + Ng), `snapspec.jl` | static, serial in ν, no moments handed back, ∫L_ν = 0.8 L |
| Lines with velocity | H₂ line transfer (comoving Doppler shifts on tangent rays from the run's velocity, GPU path); `NLTEH.βsob` (free fall only) | atomic lines need the run's own du/dr incl. non-monotone flow |
| Lyα | `LyaTransfer.jl` (Voigt, CRD core / coherent wings, ALI eliminated into the three-level atom; off) | static, no destruction by H₂/H⁻, no force, no MC calibration |
| Ionizing ray | `IonRay.jl` (photon-conserving, Iliev 1 and 5 pass) | H and H₂ only; no He, no diffuse field, no momentum |
| Opacity means | `nlte_means`, `ContinuumOpacity.means`, StellarOpacity table (has He ionization) | **no flux mean at all** (χ is Rosseland); no mean from the actual J_ν |
| Refresh pattern | H₂ ALI refresh (`ALIINT`, `ALIFAC`, rescaled escape probabilities between calls), `RAYQ`/`GASQ` energy hooks | refresh logic H₂-specific; no line-force term in the momentum row |
| Records | `recsnap!`, `lineage` | three-level parameters, starlight rates, H₂ populations, J_ν moments, b_n not recorded |

Design decisions from it: the 12-level H atom (extended to He) is the one reference; the run's three-level atom stays as
the cheap form between refreshes, with its n = 2 population and Lyα escape set from the reference at each refresh. Constants
and the Verner fits, duplicated in five modules with small differences, move to one place. Y is fixed to one value in all
modules (choice of 0.25 or 0.24 to be confirmed).

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

## Progress

| Date | Step | State |
|---|---|---|
| 12 Nov | R0 | inventory done (above); atomic-data sources compiled with bibcodes verified refereed on ADS (list in the session scratchpad `atomdata/bibcodes_verified.txt`; test values still to be read from the tables themselves) |
| 12 Nov | R1a | `HECHEM=1`: He, He⁺, He²⁺ in `ChemNet`/`ChemRates` (8a7fbd8): Voronov 1997 collisional ionisation with detailed-balance three-body recombination, Verner 1996 ground-state + hydrogenic excited-level photo-ionisation with Milne inverses. Off: bit-identical. On: reaches its own Saha limit in J = B; energy conserved to 4e-16. Missing: dielectronic recombination, He⁺ + H charge exchange, He line cooling (model atom) |
| 12 Nov | R1b | helium in the equilibrium EOS table and the matter module (34cba0a); restart padding `he_pad!`. Off: bit-identical over 300 steps (`t_heb0`/`t_heb1`). On, in the `c128` envelope (4e4 K, 1.2e-11 g cm⁻³): u +69 %, p +8 %. Partition-function cut-off differs from `Z0` (n ≤ 15 vs. occupation probabilities): to unify with the model atom |
| 12 Nov | R1 test | `c128` restarted at 183.5 M☉ with `HECHEM=1` (`t_hec1`) and without (`t_hec0`), 6000 steps (170 yr): no failed solves in either, same cost. Between 7.6 R☉ and 4300 R☉ helium is fully He²⁺ (x = 0.0833 = y_He), so the 69 % in u is a constant offset there and the structure is unchanged to ≤ 3 % (ρ), ≤ 1 % (T). He²⁺ → He⁺ recombines exactly in the two shell zones (5600–7536 R☉: He²⁺ 0.076 → 0.0013). After 170 yr: 0.19 M☉ less in the inner shell zone, L −0.4 %, R★ and Ṁ unchanged. The helium ionization front sits on the unresolved shell, so its effect cannot be judged on this mesh |
| 12 Nov | R2 | `ModelAtom.jl`: rate-matrix solver (LTE exact for J = B), hydrogenic H I / He II (n ≤ 40–200) and He I (99 terms: CHIANTI n ≤ 5 + hydrogenic shells); collisions CHIANTI n = 1 → n′ (Anderson et al. 2002; Ballance; Sawey & Berrington 1993, Bray et al. 2000) + Johnson 1972 (from the ADS scan). Tests (`modelatom_test.jl`) against Storey & Hummer 1995 and Porter et al. 2012/13: α_A,B 0.2–0.6 %; case-B Hβ 0–10 % at n_e 1e8–1e14; He II 4686/1640 13–16 % low (l-collapse); He I lines within ~25 % except up to 2× low at n_e 1e10–1e12 (approximate excited-term data). Voronov 1997 in `ChemRates` verified against two copies |
| 12 Nov | R2 snapshot | `snapnlte.jl`: the three atoms per zone on a snapshot, run's own radiation field (B(T_r) inside, diluted T_eff outside), Sobolev β in the snapshot's velocity field, n_e iterated. `c128` at 186.3 M☉ (step 1156100): envelope 15–6600 R☉ in LTE (b₂ = 0.998), H⁺, He²⁺ to the shell; lines produced there (Lyα 3e7, He II 1640 1.1e7, He I 584 8e6 L☉ locally) are buried under grey τ ≈ 1e4: **no recombination line escapes in this phase**; the visible lines come from the 5400 K infall outside the 3400 K photosphere (Hα ≈ 290 L☉, 4e-5 L; n = 2 overpopulated, b₂ 1.7e3–1e5). With the frequency-dependent LTE continuum at each line's wavelength (Z0Physics, τ_eff = √(3τ_a(τ_a+τ_s))) the result is the same: Hα 278 L☉, Hβ 13 L☉, He II and He I lines < 1e-6 L☉. Ten snapshots from step 640000 to 1200000 (≈142 → 190 M☉): Hα 206 → 278 L☉, Hβ 13–15 L☉, He II 1640 and He I 10830 buried throughout (< 1e-5 L☉) |

## R3 design: coupling the atoms into the evolution (draft, 12 Nov, for Tom)

**What each zone type already carries** (read from the code):

- *Table zones* (LTE; p ≥ 300 dyn cm⁻² and T ≥ 4500 K, or T ≥ 5000 K and τ ≥ 3): gas–radiation exchange through κ_P and κ_E
  from `StellarOpacity`, whose **Planck mean above 1e4 K already includes the H and He lines** (Z0Physics: Saha–Boltzmann levels,
  occupation probabilities), and OPAL's Rosseland mean. In LTE and at large optical depth this is the right physics, and the
  non-LTE solve on `c128` confirms it (b₂ = 0.998 through the envelope). Adding a separate line cooling there would count the
  lines twice. Missing: the flux mean (χ is the Rosseland mean) and helium ionisation in the EOS (now `HECHEM`).
- *Network zones* (non-LTE): three-level H (n = 1, 2, continuum) with Sobolev Lyα escape, H⁻, H₂ (63 levels, ALI), HD,
  continuum-only κ_P/κ_E of the network composition (`NLOPAC`): **no atomic line other than Lyα, no helium.**

**Coupling, as the H₂ lines are coupled now:** every refresh (same scheduling as the ALI: interval adapts, rescaled in between),
for the network zones and the table zones with τ_grey < τ_cut:

1. Solve H I, He I, He II (`ModelAtom`) with the run's own J (as `snapnlte.jl`; later the frequency-dependent J of item 4) and
   Sobolev β in the run's velocity field.
2. Feed back per zone:
   - **network zones**: (a) effective ionisation and recombination rates of H and He into `ChemNet` (replacing the three-level
     atom's n = 2 closure and adding He photo-ionisation from excited levels), as rate multipliers held between refreshes;
     (b) net atomic line cooling Λ_atom = Σ_lines hν (n_u A β − β (n_l B_lu − n_u B_ul) J_c) into the same `Λext` channel the
     H₂ lines use (implicit in T with its slope); (c) He line and recombination cooling likewise;
   - **all solved zones**: the departure coefficients b_i, recorded in the snapshots (records item 5 of the plan).
3. Table zones keep LTE; at the photosphere, where a table zone's b_i departs from 1 by more than a threshold, it is switched to
   the network (the existing `PSW*` hysteresis, with b_i as an added criterion) instead of correcting the table.
4. Frequency-dependent J and the opacity means: a continuum + line formal solution (`Spectrum.jl` rays, GPU path later) on the
   refresh gives J_ν per zone, hence κ_P, κ_E, κ_R and **the flux mean κ_F** (the radiation force), and the line force
   Σ (hν/c) κ_line F_ν β. The grey two-moment system keeps the hydrodynamics; these means replace the table's where they differ.
   Acceptance: ∫L_ν dν = L of the run within 2 % (else multigroup transport: a method change, back to Tom).

**Tests:** off ⇒ bit-identical; the analysis solver on a snapshot reproduces the run's per-zone terms to round-off; restarts
from frozen checkpoints at 4 (`b224`/`fz_eh10`), 42 (`f128`) and 182 M☉ (`fz_c128`) against the unchanged runs; energy budget.

**Open questions for Tom:** (i) is it acceptable that table zones stay LTE (lines through κ_P) and only the photosphere and outer
zones are non-LTE, with b_i as an extra switch criterion; (ii) the order: (2a+2b) first in the network zones of the infall
(effect at all masses), then (4) the flux mean (decisive for Γ ≈ 1 in the envelope).

Pengelly & Seaton 1964 (MNRAS 127, 165; read from the ADS scan) l-mixing rate for the l-resolved low levels (to do):
q_nl(nl → nl±1) = 9.93e-6 (μ/m)^{1/2} D_nl T^{-1/2} {11.54 + log10(T m/(D_nl μ)) + 2 log10 R̄_c} cm³ s⁻¹, D_nl = (Z/z)² 6n²(n² − l² − l − 1)
(projectile charge Z, ion charge z), 2 log10 R̄_c = min(10.95 + log10(m T τ_nl²/μ), 1.68 + log10(T/N_e)); split between l − 1 and
l + 1 in proportion to l(n² − l²) and (l + 1)(n² − (l + 1)²). Needs l-resolved photo-ionisation (Burgess 1965 recursion) as well.

**First quantitative motivation for R3 (12 Nov, `b224` step 84 560, 4.22 M☉, `snapnlte.jl` with `ModelAtom.radiative_cooling`).**
In the network zones of the infall (76–1200 R☉, 2000–5600 K) the full H atom (n ≤ 20, run's J: diluted T_eff blackbody) gives
net *heating* of 1e-5 (80–200 R☉) to 1e-8 erg cm⁻³ s⁻¹ (1200 R☉): photo-ionisation of n ≥ 2 by the star's light exceeds
line + recombination emission (lines alone 1e-5–1e-8 cooling, bound–free 2–6e-5 heating). That is 10–50× the run's grey gas–
radiation exchange there and 1e3–1e4× its H₂ line cooling. The run's three-level atom has only the n = 2 Balmer-continuum part:
its `chemq` in `t_dz1` (mesh C, 2.6 M☉) is 2e-6 → 2e-7 erg cm⁻³ s⁻¹ over 85–130 R☉, 5–50× below the full atom at similar
radii (different snapshots: an indication, not yet a like-for-like comparison). Energy balance of the atom verified: J = B and
LTE give lines + bf = 0 to 1e-14 of the gross bound–free emission (`modelatom_test.jl` §6).
