# Strategy to the end of accretion: 300 M☉ or a conclusive halt

Written 12 Nov 2026 (draft for Tom's decision). Goal and end criteria A/B/C as in `POPIII_GOAL.md`; the physics gaps as
ranked in `POPIII_PHYSICS_AUDIT.md`. This file says how to get there from today's state and in what order.

## Where we are (12 Nov, 02:30)

| Run | Physics | Mesh | State | Progress |
|---|---|---|---|---|
| `c128` | full production set¹, no ionizing light | 128 fixed | **181.7 M☉**, H burning (T_c 1.37e8 K), R★ 7.6 R☉, photosphere ≈ 7500 R☉ at T_eff ≈ 3400 K, L ≈ 7.0e6 L☉ ≈ L_Edd(es), Ṁ 2.8e-3 M☉/yr (steady over the last 4000 yr) | 2.5 M☉ per wall hour → 300 M☉ in ≈ 2 days if Ṁ holds |
| `f128` | same | 128 fixed | 42 M☉, divergence burst | 1.1 M☉/h |
| `b224` | same + HD | 224 fixed (mesh C) | 4.1 M☉, quiet phase | 0.8 M☉/h; envelope under-delivers mass by 9–21 % beyond 10⁴ R☉ (too coarse there) |
| `t_mm2`, `t_mm3` | same | moving mesh (`GRIDMAP`), mm3 with frozen envelope beyond 3000 R☉ | mesh tests | — |

¹ `STARLIGHT SL3 NLOPAC NLTE=hyb HDCHEM LINEINC EDDGMU0 EDDSCAT ETOTAL HOT`.

Earlier runs without radiative feedback (`_a8` 352 M☉, `_q9` 312 M☉) already passed 300 M☉, so **reaching 300 M☉ is not the
open question; whether a physically complete star still does is.** The one process expected to end accretion, the
H II region and its pressure, is not in any running calculation (audit gap 1).

`c128` is already in the regime where the end-B outcome is decided: L at the Eddington value, an inflated envelope, and a
23.5 M☉ shell piling up near 5600 R☉ (`POPIII_RUNS.md`). Every one of the top three audit gaps acts in exactly this layer.

## The plan in four stages

### Stage 1 — finish the coarse solution and map where each missing process matters (now, ≈ 2 days, no new physics)

- Let `c128` run on to 300 M☉ or to a halt, unchanged. It is the scaffold: one complete trajectory M★(t), R★, L, Ṁ.
- Freeze a checkpoint every 10 M☉ from 180 M☉ on (plus the frozen checkpoints of earlier masses that already exist) as restart points for stage 2.
- **Physics-trigger map** from the `c128` snapshots (post-processing only): at each mass, the ionizing photon rate of
  the hot surface below the envelope (blackbody at the local T where τ_R = 1 from inside, and at the stellar surface)
  against the recombinations the envelope column can absorb (∫ n_e n_p α_B dV); Γ = κ_F L / 4πGMc zone by zone; the
  mass and pressure of the He/He⁺ ionization layer; the H₂ fraction and LW optical depth outside the photosphere. Output:
  the mass at which an H II region could break out, and which layers carry Γ ≥ 1. This tells us which of the physics
  below changes the answer and from which mass, instead of adding everything everywhere.
- Energy-budget check of any slowdown before calling it a halt (rule from `POPIII_GOAL.md`).

### First trigger map (12 Nov, `c128` 141.5–181.8 M☉; script `trig.jl`, checked by hand on step 1101250)

- **No H II breakout possible in this range.** The whole envelope between R★ (7.6 R☉) and 5600 R☉ is at 35–160 kK, i.e.
  already ionized in the LTE table. Its recombinations (≈1e59 s⁻¹) exceed any ionizing photon rate of the star by 10⁵–10⁹.
  The H II region is *trapped* inside the envelope.
- **A 24 M☉ radiation-supported shell.** Two zones, 4300–7536 R☉ (Δln r 0.27 and 0.29), hold 182.8 → 205.7 M☉ of enclosed
  mass, i.e. 23.9 M☉, nearly at rest (u ≈ −0.01 to −6 km/s), with Γ up to 1.55. Outside it the density drops by 2500 in
  one interface and the infall resumes at 23–60 km/s. The photosphere sits on this shell. The star still gains
  2.8e-3 M☉/yr, the same as the inflow at 10³–10⁵ R☉, so this is not a halt.
- **Everything that decides end B here lives in two zones.** The shell and the jump outside it are not resolved on the
  128-zone mesh, and the shell is in the layer where the EOS lacks helium ionization (all 17 envelope zones, 23.9 M☉,
  are table zones with p < 1e4 dyn cm⁻² and 1e4 < T < 1.5e5 K). L/L_Edd(es) 0.98 → 1.027.
- Consequence for the plan: **a resolution restart of the shell region at ≈180 M☉ joins stage 1** (refine 10³–3×10⁴ R☉,
  restart from a frozen `c128` checkpoint), and P1 (helium) stays first in stage 2. The ionizing ray (P2) matters only
  once the envelope is no longer opaque, which the map does not show before 182 M☉.

### Stage 2 — complete the physics that can change the end, one switch at a time (≈ 1–2 weeks)

Order set by impact on end B. Each item: default-off switch, a unit/benchmark test, then a restart test from the stage-1
checkpoints at the masses where the trigger map says it matters, compared with the unchanged `c128`. Adoption in the
production set is Tom's decision for each.

| # | Process | Why it can change the end | Work | Test |
|---|---|---|---|---|
| P1 | **Helium ionization in the low-pressure table EOS** (audit gap 3) | the inflated layer at L_Edd is set by the opacity bumps; today the opacity has He ionization and the EOS does not (79 zones, 0.31 M☉ of `fx3`). Inconsistent exactly where Γ ≈ 1 | rebuild the equilibrium table with He⁰/He⁺/He²⁺ Saha in the low-p part (the opacity table already has it) | table vs. direct Saha; restart from 180 M☉: does the envelope radius / the 5600 R☉ shell change? |
| P2 | **Ionizing radiation: the H II region** (gap 1, milestone M3) | the expected end of accretion | `IonRay.jl` exists and passes Iliev tests 1 and 5; still missing: validated coupling to the grey field, H II temperature (photoheating vs. cooling), helium absorption, on-the-spot diffuse field, and photoionization **momentum** (the force on the infall) | Iliev test 1/5 with the coupling on; Strömgren radius in a frozen `c128` envelope; then restarts from the mass where the trigger map predicts breakout |
| P3 | **Non-grey radiation force** + Lyα pressure estimate (gap 2, 6) | Γ ≈ 1–2.4 in the ionized layer: a 20 % error in κ_F decides whether the layer is bound | flux-mean weighted with the local spectrum in thin zones (from the formal solution already computed for the Eddington factor), Lyα pressure by the McKee & Tan 2008 estimate first, transfer only if it matters | flux-mean vs. Rosseland in optically thick zones (must agree); restart comparison |
| P4 | Lyman–Werner with level populations, H⁻ heating at the radiation temperature (`NLOPACJ=1`, exists), Compton exchange, shock radiation (`SHOCKRAD`, exists) (gaps 4, 5, 7, 8) | set the H₂ fraction and the temperature of the infall, so the line emission and the infall rate at large radii; unlikely to stop accretion on their own | mostly switches that already exist: test them | restart comparisons at 40 and 180 M☉ |

Not planned (stated as limits of the 1D model in `POPIII_GOAL.md`): rotation, disk, fragmentation, magnetic fields.

### Stage 3 — make the final run affordable and resolved (in parallel with stage 2)

- **Mesh.** Decide between the fixed mesh C, a fine-envelope fixed mesh, and `GRIDMAP` with a frozen envelope
  (`t_mm3`, running): criterion = enclosed mass at fixed radii against the fine fixed mesh (`t_fe2`), R★(M★), and wall
  time per simulated year. The infall must be resolved beyond 10⁴ R☉ (b224 is not).
- **The repeating failed solve at the dissociation front** (moving mesh ×2.5 steps per year): at Δt ≈ 7e6 s the gas
  in the zones just outside the front is replaced ≈ 450 times per step, yet the chemistry step integrates each zone as
  a closed box for the whole Δt. A chemistry step that follows the flow through the front (named in `POPIII_RUNS.md`,
  not built) would remove the split-step swing on all meshes. Check first: does the H₂ profile there change with Δt?
- **Envelope cycle** (unresolved on the coarse meshes): resolution study from restarts at 104 and 180 M☉ only.

### Stage 4 — the final calculation

One run from the initial conditions with the adopted physics and mesh, through to end A or B, plus resolution and
restart repeatability checks at three or four epochs (adiabatic, swelling, ignition, end). Budget gate unchanged: a
projected wall time above two weeks goes back to Tom.

## What counts as conclusive

- **End A:** 300 M☉ inside the shock in the final run, with P1–P3 on, mesh-converged R★(M★) and Ṁ to ~10 %.
- **End B:** the halt criterion of `POPIII_GOAL.md` met, the energy budget closed, the process that set it named (H II
  pressure, radiation force on the ionized layer, or envelope inflation), and the halt mass agreeing within ~10 % between
  two meshes and between restarts at nearby epochs.
- Either way: M★(t), R★, L, Ṁ, and the state of the infall at each stage (the input to the line emission).

## How the work is split between models

| Task | Model |
|---|---|
| Physics design, solver code, deciding what a result means, anything brought to Tom | Opus |
| Post-processing scripts over snapshots (trigger map, enclosed-mass comparisons, mesh census), test-run comparisons against a written criterion, literature/ADS checks, first-pass crash triage (read the log tail, find the zone and equation, report) | Sonnet subagents |
| Routine watch turns (status lines, "is anything stopped", panel refresh) | Haiku subagents or the existing shell watchdog |

Subagents get the run rules in their brief (never touch Tom's jobs or kill by pattern, never print environment tokens,
no `rm`, no in-place edits of `fullrun.sh`/env files, scratch writes only under the ZFS scratch) and report; they do not
restart production or change flags.

## Decisions needed from Tom

1. Keep `c128` running unchanged to 300 M☉ / halt as the stage-1 scaffold (recommended). Plus a resolution restart of its envelope shell at ≈180 M☉ (recommended).
2. `f128` (42 M☉, slow, same physics as `c128`): stop it to free cores for stage-2 restarts, or keep it as a second realization?
3. Order of stage 2: P1 (He EOS) before P2 (H II region), as above, or the H II region first?
4. `b224`: keep as the HD from-IC baseline despite its coarse envelope, or replace it once the mesh decision of stage 3 is made?
