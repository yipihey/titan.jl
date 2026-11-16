# Moving mesh for the evolved zdiagB state (branch `mm-evolved`, 16 Nov 2026)

Question: can the evolved state of the from-IC run `zdiagB` (abn02red, fixed log-Eulerian mesh of 290 zones, full stack),
which stalled at Δt ~ 1e2 s from step ~800 and died at step 3032 (Δt < DTSTOP) because the photosphere / H ionization / H₂
dissociation front sits across one interface, be continued on a moving, adaptive mesh that resolves the front?
(The sub-zone front treatment on the fixed mesh is another agent's work and is not touched here.)

Branch `mm-evolved` = popiii-s6 0781136 + `grideq` b50ec58 (GRIDEQ) + the changes below; everything default off.

**Short answer.** Yes, GRIDEQ restarts directly from zdiagB's fixed-mesh checkpoint and carries the run past the step-3032
failure (where the fixed mesh dies within ~30 steps, Δt → 1e-5 s), with the front spread over 6–11 interfaces at first,
Δt rising to 1e2–7e3 s and only 1–10 failed solves per 100 steps; the physical time advanced by 2.2e5 s in 300 steps
(zdiagB needed 2200 steps for 1.8e6 s). **But it does not reach a production Δt (≥ 1e4 s) and it does not keep the front
resolved:** within ~200–400 steps the cliff re-sharpens to 1–2 interfaces even where Δln r = 0.0023, and Δt falls back to
1e1–1e2 s. The cliff is pinned to the **table → network switch** (eqzone 1 → 0): in every snapshot of every run the steepest
interface (τ 800 → 3, T 6800 → 3300 K, x_H₂ 2e-4 → 0.14, T_rad/T 1.8 on the network side) is exactly the interface between the
last table zone and the first network zone. That jump is a per-zone discrete model choice decided once per step (a lagged
input), and no radial resolution can spread it. Moving mesh = necessary-ish (it removes the one-zone geometry and gets the run
going again) but not sufficient (also with an iterate-only monitor, t_mm22); the switch has to become a smooth (iterate-evaluated) blend across the front first.

## 1. The merge

`grideq` (b50ec58) was built on iceul/g4 (ce50412), 34 commits behind popiii-s6. Conflicts were in `iceulbuild!`
(popiii-s6: abn02fin/abn02red ICs and ENVLIN central zones; grideq: the `rIset` starting mesh and IC=p100) and in
docs/POPIII_RUNS.md. Resolved by keeping popiii-s6's `iceulbuild!` and adding only `rIset` (GRIDEQ's starting mesh at ICEUL);
p100 is not taken.

**Not taken from grideq's history, on purpose:** `PSWDEEP`, `NLOPACJTAU`, `FLUXTIECORE` (PopIIIMaterial.jl, TitanCore.jl) and
the `CONV` env switch of popiii_m2c.jl. They are not in popiii-s6, **but zdiagB's run.env sets all four (`PSWDEEP=1e3
NLOPACJTAU=1,10 FLUXTIECORE=1 CONV=1e-3`) and they are inert there** (zdiagB ran popiii-s6 94f75f3, which has none of them).
Merging them would silently switch them on for every run that restarts with zdiagB's env (the bit-identity test below would
fail). If they are wanted they need their own decision (and a renamed env for zdiagB-like runs). scripts/popiii_run is
unchanged (grideq's iceul.env and the iceulrun.sh comment are dropped).

## 2. What it takes (question 1)

GRIDEQ (mode 1, the implicit grid equation, TitanResidFast.grideq_rows!) is already the right structure for this restart:

- **Mesh requirements** — `ENVEUL=0,…` (regular centre), `lgrid = 2`: zdiagB's checkpoint qualifies as is. The first interface
  above the centre (0.1 R☉, ENVLIN=1 central zone) and the outer edge (1e9 R☉ = 22.5 pc) are the fixed ends of the grid
  equation; every interface in between is an unknown of TITAN's r row, replaced by n̂_j = n̂_{j−1}.
- **State** — none beyond the mesh: a fixed-mesh checkpoint restarts directly with `GRIDEQ=1`. No remap is needed: TITAN's
  moving-grid (ALE) form of the conservation rows moves mass, momentum, energy and E_r through the moving faces (the same
  residual as GRIDWIN/GRIDMAP), i.e. the "remap" is done conservatively by the solver itself, continuously; the composition
  follows with `advect_species!` (mass crossing each interface = change of x_m,e there, which includes the grid motion; XADV2
  van Leer). Zones keep their identity, so the per-zone material state (eqzone table/network flag, COUPLED species, ALI
  factors, LINEIMP models, LOCALJ, NLOPAC/NLOIT) stays attached to its gas.
- **VEFIMP** — the formal solution is evaluated at every Newton iterate on the iterate's radii: nothing to change.
- **COUPLED rows** — local per zone; species advection is operator-split after the step (as on the fixed mesh).
- **DTERR** — TR-BDF2 with lag = GRIDMOVE (radii combined as for lgrid 3); its error norm includes r.
- **Jacobian** — the grid row of interface k couples r_{k−2…k+2}, ρ, P, T, E_r of zones k−2…k+1 and x_m,e: inside the 5-block
  band of the coloured ForwardDiff sweep (5 colours, ±2 blocks). The τ and composition terms are frozen per step (no
  Jacobian entries, no extra unknowns).

What was missing / wrong, and is now on the branch (all default off):

| change | why |
|---|---|
| fix: `tr_ = max(min(κ t_c, κ_f t_fe), …)` with κ_f = 0 gave 0·Inf = NaN in every grid row | the default GRIDEQ (no feature tracking) failed in its first Newton iteration (`t_mm14`: NaN in zone 1, 60 attempts) |
| `GEQX=wx[,xfl]`: composition term wx Σ_s (∂ln(x_s + xfl)/∂ln r)² over every network species, frozen per step (and in the explicit target) | the front is an ionization + H₂ dissociation front; the rule had no composition |
| `GEQRAMP=κ0,n[,f]`: κ (and the speed floor) rise geometrically κ0 → κ over n steps, only if the mesh is far from the rule at the start (max M Δln r / median > f, default 10) | κ 30 adapts the mesh over ~30 crossing times of the front (~1e7 s); the fixed-mesh state fails within 30 steps (`t_mm13`, `t_mm15`) |
| `GEQREMAP=f`: one-time rebuild of the restarted state on the rule's mesh (envbuild! from an lgrid 2 source, ENVCONS mass, XMZERO again), with a mass/energy report | the alternative to the ramp; tested and rejected (below) |
| `GEQW` entry 8 `wχ`: (∂ln(χρ)/∂ln r)² of the Newton iterate in the grid row (and in the explicit rule) | a local, non-lagged substitute for the frozen ln τ term (Tom: no lagged inputs in the implicit solve) |
| grid-eq line: interfaces with 1 ≤ τ ≤ 1e5, largest \|Δln T\| per interface, current κ | front diagnostics |

## 3. Monitor function (question 2)

GRIDEQ's rule (GridEq.jl): every zone between the fixed ends covers the same ∫M d ln r, M = √(1 + G²),
G² = Σ_q w_q (∂ln q/∂ln r)², q = ρ, P, T, E_r, m, τ (clamped to [0.1, 1e4]), + (GEQX) ln(x_s + 1e-3) of every species;
G² capped smoothly at gmax²; neighbours within (1+α)/α = 1.67. All terms are logarithmic derivatives: no scale; a power law
gives a constant M (the logarithmic mesh stays). The only non-logarithmic choice is the abundance floor x_fl = 1e-3 (species
below it do not drive the mesh; it is a share of nuclei, not a tuning to this front). Recommended general set (iterate-only
monitor, section 4 `t_mm22`; zone budget stays ~290):

    GRIDEQ=1 GEQC=100,1.5,30,0.5,0.5,0.1,1e4,0.3,0.3 GEQW=1,1,1,1,1,0,0,3 [GEQRAMP=1e-3,300 on a restart from a fixed mesh]

(the grideq agent's t_gse7 set, gmax 300, w_τ 3, κ_f 0.3, + GEQX resolves more zones at the front, t_mm17, but costs Δt and
uses two frozen terms; floor/ceiling: the finest zone is ~(∫M d ln r / N)/gmax ≈ 0.006 (gmax 100) or 0.0023 (gmax 300) in
ln r, the coarsest ~0.19 in the free-falling cloud.)

Where the rule puts 290 zones for the zdiagB S5 state (GridEq.target on the snapshot of step 3000; `scratchpad/mm/target_scan.jl`):

| monitor | interfaces 1 ≤ τ ≤ 1e5 | 0.1 ≤ τ ≤ 1e4 | Δln r min | zones < 10 / 10–20 / 20–1e4 / > 1e4 R☉ | Δln r max |
|---|---|---|---|---|---|
| fixed mesh (zdiagB) | 8 | 25 | 0.0797 | 58 / 8 / 79 / 144 | 0.080 |
| GRIDEQ default (w 1, gmax 100) | 27 | 38 | 0.0070 | 60 / 39 / 64 / 126 | 0.166 |
| + GEQX 1 | 31 | 41 | 0.0060 | 60 / 46 / 61 / 122 | 0.182 |
| gmax 300, w_τ 3 | 38 | 50 | 0.0043 | 56 / 53 / 62 / 118 | 0.184 |
| gmax 300, w_τ 3, GEQX 1 | 41 | 53 | 0.0039 | 57 / 57 / 61 / 114 | 0.195 |
| gmax 1000, GEQX 1 | 36 | 46 | 0.0047 | 59 / 51 / 60 / 119 | 0.188 |

Arc length ∫M d ln r = 105 (111 with GEQX): star < 10 R☉ 21, front 10–20 R☉ 15 (18), envelope 23, cloud > 1e4 R☉ 45. So ~15 %
of the zones go to the front and the cloud coarsens from Δln r 0.08 to ≤ 0.19 (it is a ρ ∝ r⁻² free fall: harmless).
(τ counts here come from ln τ interpolated across the old one-interface jump, i.e. what the starting mesh offers; the
resolved front is then set by the physics, section 4.)

Which quantity: ρ, P, T, E_r and m are at the Newton iterate in the grid row; τ (non-local: a sum to the edge) and the
composition (operator-split species) are frozen per step in `sx`. Following Tom's rule against lagged inputs, `GEQW` entry 8
(wχ) adds ln(χρ) of the iterate (local; the opacity cliff is what τ measures), so `GEQW=1,1,1,1,1,0,0,3` with GEQX off is a
monitor without any frozen term (t_mm22). The remaining per-step coefficients are the relaxation time τ_v = κ t_cross (and
κ_f's feature speed from the last two steps' M) and the Δt cap fmove: they set how fast the mesh may move, not where it goes.

Accuracy can be raised later at low cost by the zone count alone: the rule is resolution-free (no budgets), so N → 2N halves every
zone's ∫M d ln r; the front's share stays ~15 % of the zones.

## 4. Tests from zdiagB S5 (question 3)

All runs: copies of `zdiagB/ck/m2c_z19zdiagBS5.*` (step 3000, t = 7.540741e12 s), zdiagB's run.env + phaseB.env (+ the
settings shown, `MATTHREADS=4`), 4 threads, `numactl --physcpubind=96-99|100-103 --membind=12`, via the test tree's
scripts/popiii_run/iceulrun.sh (`git archive` trees under /zpool/nvme/data/tabel_scratch/titan_popiii/mm_trees).

| run | settings | result |
|---|---|---|
| `t_mm11` | popiii-s6 0781136 tree, fixed mesh (reference) | Δt 2.6e-5 s at step 3019 (as zdiagB S6, which died at 3032) |
| `t_mm12`, `t_mm23` | mm-evolved trees 1ccd23b / b9ad619, GRIDEQ off | bit-identical to t_mm11 (section 7) |
| `t_mm14` | GRIDEQ=1 GEQX=1 GEQRAMP (default κ_f = 0) | NaN in every grid row at the first iterate (the κ_f = 0 bug, fixed in 9294424) |
| `t_mm13`, `t_mm15` | GRIDEQ=1 + ramp from κ0 = 0.03 with the speed floor kept | the mesh barely moves (floor Δln r t_c/f_v ≈ 6e4 s at the front vs Δt ~ 1e2 s): Δt 1e-8 s by step ~3030, as the fixed mesh |
| `t_mm16`, `t_mm18` | GEQREMAP=3 (one-time rebuild on the rule's mesh), then moving (16) or frozen (18, κ 1e12) | Δt → 5e-9 s by step 3045: the interpolated front (T, x linear in ln r across the old jump: 6000 K gas with x_H₂ 0.07, T_rad ≠ T) is far from equilibrium; repeated Newton divergence in the radiation row of the first coarse zone after the refined block (zone 112). Rejected. Remap energy: E_int −8e-6, E_rad −3.6e-4 (the CMB bath of the coarsened cloud), mass exact (ENVCONS). |
| `t_mm17` | GRIDEQ=1 GEQRAMP=1e-3,300 GEQC=300,1.5,30,0.5,0.5,0.1,1e4,0.3,0.3 GEQW=1,1,1,1,1,3 GEQX=1 | **passes 3032**; 3050–3250: 10–11 interfaces in the cliff (3 < τ < 1e5), largest per-interface jumps Δln T 0.16–0.18 (fixed mesh 1.04), Δln τ 1.8–2.0 (12.0), Δln(x_H⁺+1e-3) 0.87–1.0 (3.7), Δln(x_H₂+1e-3) 1.1–1.3 (5.2); Δln r 0.0023–0.0034 at the cliff. Δt 0.9 → 19 → 285 s (steps 3100/3200/3300), peaks 1.9e3 s; failed solves per 100 steps 6, 1, 9, 21. From ~3300 the cliff re-sharpens: 8 → 4 → 2 interfaces (3300/3350/3400), Δln τ 6 across one interface of Δln r 0.0025, always at the table→network interface; Δt 10–170 s. Stopped at ~3450. 7.2 s/step (4 threads). |
| `t_mm19` | as t_mm17 but GRIDEQ's default gmax 100, w_τ 1 | passes 3032; 6 cliff interfaces to 3150 (Δln τ 3.0–3.3), Δt 20 → 555 s, peaks 7.2e3 s, fails 5/6/18 per 100 steps; at 3200 the front jumped from 14.1 to 12.35 R☉ and re-sharpened to 2 interfaces (Δln τ 8.5, Δln T 0.7–0.9, again at the switch); Δt 1–9 s at 3250–3300. +2.2e5 s of physical time and M★ 0.0611 → 0.0624 M☉ in 300 steps. Stopped at 3330. |
| `t_mm22` | GRIDEQ=1 GEQRAMP=1e-3,300 GEQC=100,…,0.3 **GEQW=1,1,1,1,1,0,0,3**, no GEQX: no frozen monitor term (ln χρ of the iterate replaces ln τ) | passes 3032 fastest: +8.2e4 s of physical time in the first 100 steps (Δt peaks 1.8e3 s at 3059, 10 failed solves per 100 steps, 9.4 s/step); at 3050 6 cliff interfaces (Δln τ 3.2, Δln T 0.26); by 3100 the front has moved in from 14.1 to 12.2 R☉ and is again one interface, the table→network one (zone 65 table 7840 K, τ 1.3e4 / zone 66 network 3450 K, τ 15; Δln τ 6.8, Δln T 0.82) at Δln r 0.007; Δt 2 s. Stopped at ~3110. Same outcome as with the frozen τ/composition terms, so the lagged monitor terms are not what pins the front. |

Overlap with zdiagB (it has no time beyond step 3032, so the comparison is at the restart and against zdiagB's own trend):
T_eff 2900 K (zdiagB S5) → 2180–2300 K in t_mm17/19 within ~100 steps and 2500–3040 K later (zdiagB swung 2460–2940 K
in its last 600 steps); R★ (perf.csv) 11.0 R☉ → 11.2–11.4 (t_mm17), 10.9–12.1 (t_mm19); Ṁ at the surface 0.095 M☉/yr →
0.05–0.097 (noisy, as the surface radius jumps), at 1e4 R☉ unchanged 0.0228 M☉/yr. M★ 0.0611 → 0.0613 (t_mm17, 3400),
0.0624 (t_mm19, 3300); zdiagB grew 0.0558 → 0.0611 in 2200 steps.

**Conservation.** Mass: |ΔM/M| ≤ 4e-14 over every run (outer wall, regular centre; the grid motion moves no mass across the
ends). Energy, from TITAN's own ledger (`tetot` of totnrg, erg sr⁻¹, read from the checkpoints): zdiagB 2400 → 3000
(1.5e6 s): +1.75e42. The ramp phase (steps 3000 → 3199, fast mesh motion) adds +3.2e42 (t_mm17, gmax 300) / +1.3e42 (t_mm19)
— a one-time 1–3 % of the protostar's binding energy per steradian (G M★²/R★/4π ≈ 1e44), the cost of moving ~100 interfaces
by up to several zone widths through the star and front with the van Leer donor-cell fluxes; afterwards the drift is
−2.4e41 (t_mm17) and +5.6e40 (t_mm19) per 100 steps, as small as the fixed mesh's. During the ramp the star's mass moves
measurably between fixed radii (t_mm17 3050→3100: −3.1e-5 M☉ inside 1 R☉, +2.4e-5 in 1–10 R☉ in 19 s), i.e. numerical
diffusion from the fast motion (the grid speed reached 5000 (|u| + c_s) with the floor ramped by 3e-5): the ramp's floor
should not be relaxed this far (keep e.g. f_v ≥ 10) — a refinement for the next attempt.

## 5. Cost and risks (question 4)

- **Cost per step**: 5.4–7.3 s/step on 4 threads (zdiagB: 5.5 s on 8); the grid rows add no unknowns (they replace the fixed
  r rows) and stay within the ±2-block colouring; the monitor's frozen terms and GEQX cost one ChemNet.full per zone per step.
- **Δt from the mesh**: fmove (≤ 0.5 of the narrower zone per step) caps Δt while the mesh adapts (t_mm17: Δt 0.2–3 s for the
  first ~100 steps); afterwards the DTERR limiter is E_r / x_H⁺ / x_H₂ in the cliff zones, as on the fixed mesh. Finer zones at
  the cliff cost Δt: gmax 300 (Δln r 0.0023) ran at 3–10× smaller Δt than gmax 100 (0.006), the grideq agent's
  "factor 5–10 per factor 2–3 in zone width".
- **Tangling**: none seen (neighbour ratio ≤ 1.40, min Δln r 0.0023); the Dorfi–Drury bound and ctol (0.3 relative width
  change per iteration) hold; GEQREMAP is the risky path (inconsistent state), not tangling.
- **Species advection**: x moves with the mass through moving faces (advect_species!, XADV2); fast mesh motion (ramp) smears
  the composition and the star (above). With the ramp done the mesh moves at ≤ f_v (|u| + c_s).
- **VEFIMP / COUPLED**: no change needed and no failure attributable to them; the limiting rows in failed solves were the
  radiation rows of the cliff zones (moving) or of the first coarse zone after a refined block (remap).
- **Lagged inputs in the front zones** (Tom's point): the moving mesh's own frozen terms (ln τ, GEQX composition, the
  relaxation time/feature speed) can be avoided (t_mm22: ln χρ at the iterate). The physics in the cliff zones has several
  per-step (lagged) inputs, and the cliff sits exactly on one of them:
  1. **the table/network switch** (eqzone, decided in chem_core! after each step from the step's T, P, τ): the cliff is always
     the interface between the last table and the first network zone (t_mm17 3400/3450, t_mm19 3200/3300, zdiagB 3000); the
     network zone just outside has T 3300–3900 K under T_rad 3800–6100 K. Not removable by the mesh.
  2. starlight rates of the three-level H / H⁻ / LW network (STARLIGHT, SL3, SLEB, SLBLEND: "rates lagged over the step,
     start-of-step T", T_c of the light), which set x_H⁺ in the first network zones;
  3. NLOPACJ's κ_E weighting with the lagged T_c (NLOIT evaluates the opacity at the iterate's populations but with the
     lagged radiation colour);
  4. LINEIMP's radiative rates of the last transfer call and the ALI factors (H₂ line cooling of the cold side).
  The flux tie's h is not lagged (VEFIMP=1 evaluates it at every iterate). The DTERR limiting variables in the moving runs were
  E_r, x_H⁺ and x_H₂ of the two zones at the switch, consistent with 1–2.

## 6. GRIDMAP, GRIDSEG

- **GRIDMAP** (in popiii-s6): one mesh from arc length in ln(r + r_c), ln(m + f_c M★) and asinh((m − M★)/(ε M★)) around the
  star's anchor mass, prescribed interface velocities (GRIDVEL), optional GMTAU (ln τ coordinate) and GMWIN (radius window).
  It restarts from any lgrid 2 mesh, but it is built around markers (M★ from the velocity profile, windows in R★) and its
  τ coordinate needs a Δln r floor < 0.02; with floor 0.005 two tests crashed (`t_gt0/1`, Δt → 5e-10 s). Not general.
- **GRIDSEG** (branch, e5d0bd4): three segments with fixed zone budgets and markers (m_cl, r_b, T < 2500 K & x_H₂ ≥ 0.05):
  case-specific by design, ruled out by Tom; every mesh fine at the stellar surface stalled in phase B within ~1000 steps
  (table/network flip outside the shock).
- **GRIDEQ** is the only one without markers or budgets; it went to 54.4 M☉ from the IC (`t_gse3`).

## 7. Bit identity

With GRIDEQ unset, 30 steps from a copy of zdiagB S5 (4 threads, cores 96–99/100–103): popiii-s6 0781136 (`t_mm11`) against
mm-evolved 1ccd23b (`t_mm12`) and against the final code b9ad619 (`t_mm23`): all 591 arrays of the final state byte-identical
(md5 of the decompressed tstate 28ef61e416837aaf in all three) and the .aux identical (08a3ecc0171446a8). The runs used
zdiagB's env, which names PSWDEEP/NLOPACJTAU/FLUXTIECORE/CONV: this is why those grideq-history switches are not merged.

## 8. Recommendation

1. Keep GRIDEQ (mode 1) as the moving-mesh route: it restarts from fixed-mesh checkpoints, needs no remap and carried the dying
   zdiagB state 300–450 steps further with the front on 6–11 interfaces for ~200 steps. Use the iterate-only monitor
   (`GEQW=1,1,1,1,1,0,0,3`, no GEQX) and `GEQRAMP` only for the restart from a fixed mesh (do not use GEQREMAP).
2. It will not give production Δt while the photosphere is pinned to the table → network switch. The next step is not more
   mesh but a smooth, iterate-evaluated table/network blend across the front (or the network everywhere above the switch
   pressure in the front), together with iterate (not lagged) starlight/T_c inputs in the cliff zones; then the same GRIDEQ
   mesh can be tested again from S5 (or better from an earlier, healthy checkpoint, e.g. S1 at step 600, Δt 1.5e4 s, before the
   front locked into one zone). This is complementary to the sub-zone front work on the fixed mesh.
3. Accuracy later at low cost: double the zone count (no budgets in the rule) and/or raise gmax; expect Δt to drop with the
   cliff zone width as long as the switch is a step.

## 9. Second round (16 Nov, on popiii-s6 ae589ab: SUBZONE + SLIMP=2 + g2-material, network carries the photosphere)

`mm-evolved` merged popiii-s6 ae589ab cleanly (9d09fbe). Bit identity with GRIDEQ off: 30 steps from `t_g2e5/ck/m2c_g2e5`
(step 3300) with z300g2's env, ae589ab (`t_mm26`) against mm-evolved f49c878 (`t_mm27`): final state md5 2ff42784585537fb and
aux 7a27a1c69f224d26 in both (591 arrays). New since round 1: `GEQRAMP` entry 4 (fvr, cap of the ramped speed floor, used 20),
`GEQDBG=1|2` (diagnostics: first non-finite grid row / residual or Jacobian entry).

All runs: copy of m2c_g2e5 (step 3300), z300g2's run.env + phaseB.env (+ settings below, MATTHREADS=4, EQEOS_CACHE copied
into the run dir), 4 threads on 96–99 / 100–103. Monitor everywhere the iterate-only one: `GEQW=1,1,1,1,1,0,0,3` (ln ρ, P, T,
E_r, m and ln χρ of the Newton iterate; no frozen τ or composition term), `GEQC=100,1.5,30,0.5,0.5,0.1,1e4,0.3,0.3`.

| run | settings | result (status at report time, 16 Nov ~11:45) |
|---|---|---|
| `t_mm24` | + `GEQRAMP=1e-3,300,10,20`, SUBZONE=1,2 (as z300g2) | ramp phase 3300–3600: Δt 0.5–10 s (held by the fmove motion cap, error 0.01), 50 failed solves in 500 steps; after it Δt 1.5e3–1.6e5 s (perf.csv 4.4e4 at 3700, 4.5e4 at 4000), 3–6 failed solves per 100 steps, 4.4–6 s/step; segment S2 did not ramp again (quality 2.25). At step ~4050, running. |
| `t_mm31` | no ramp (κ 30 from the start), SUBZONE=1,2 | the mesh adapts over ~200 steps at Δt 30–1e3 s (fmove), 20 failed solves per 100 steps; then Δt 3e4–8e4 s; at ~3700, running. |
| `t_mm25`, `t_mm28`, `t_mm29`, `t_mm30` | as t_mm24 with **SUBZONE=0** (28: GEQDBG=1, 29: NANSCAN=1, 30: GEQDBG=2) | **fails at step 3319** (19 steps after the restart, deterministic): NaN in every Newton correction at the first iterate down to Δt 1e-36. The state after the hydro and chemistry is finite (NANSCAN silent) and no grid row is non-finite (GEQDBG=1 silent); the first non-finite Jacobian entries (row 5, block 227 in the cloud) are there from step 3300 on, so they are not the trigger. Not pinned down further; with SUBZONE on the same steps pass. |

**Does the front stay resolved? No.** Per snapshot (t_mm24, 15 snapshots 3350–4000; t_mm31 7 snapshots): after a short
spell at the end of the first ramp steps (t_mm24 step 3400: 6 interfaces with 3 < τ < 1e5, Δln T 0.34 per interface) the
cliff is across **one interface** in every later snapshot: Δln T 0.78–0.94, Δln τ 6.3–9.1, Δln(x_H⁺+1e-3) 2.5–4.1,
Δln(x_H₂+1e-3) 2.9–5.2 per interface, at Δln r 0.0076–0.0092 — ten times finer than the fixed mesh (0.08), where the same
jump is Δln T 0.95 (Δlog T 0.41). Both zones are network zones now (no switch behind it). **The jump does not shrink with
the zone width**, i.e. the front is thinner than any zone the mesh can afford.

**What drives it (diagnosis).** The opaque side of the cliff is optically very thick per zone even on the refined mesh:
t_mm24 step 3500, zone 66: Δτ_zone 4.4e5 (T 11 400 K), zone 67: 1.6e4 (8 740 K), zone 68: 1.5 (3 510 K, T_rad 4 170 K),
zone 69: 0.035. Zone 67 is 1.6e4 photon mean free paths thick over Δr = 8.6e9 cm, i.e. ℓ ≈ 5e5 cm = 4e-7 in ln r. The H⁻ /
H bound-free opacity rises steeply with T (and x_e through Saha), so the gas cools until it becomes transparent and the
transition from opaque to transparent happens within about one mean free path of the opaque side: a radiative
transparency (cooling) front, a discontinuity on every affordable mesh. A gradient monitor then always saturates at its
cap (gmax 100 → Δln r 0.0075 here; gmax 300 → 0.0023 in round 1, still 1–2 interfaces): resolving the τ 1e5 → 1 transition
needs Δln r ~ 1e-6, and at the front's speed (it moves 20 → 17 → 22 R☉, ~10 km/s at times) the motion cap (fmove) would
then hold Δt at seconds. So the sub-zone closure is the physics, not a crutch: **SUBZONE stays necessary** (and its
switch-off fails, above).

**Physics vs z300g2 over the overlap** (z300g2 is the fixed mesh; t at the same physical time; z300g2's R★ moves in steps of
the zone width 0.08 in ln r): t = 238 933.94–.96 yr: M★ 0.0999 (z300g2) / 0.1004 (t_mm24), R★ 17.75 / 17.06 R☉; t = .12 yr:
M★ ≈ 0.112 / 0.1130, R★ 19.2 / 19.37 R☉. Ṁ at 1e4 R☉ 0.0229 / 0.0229–0.0230 M☉/yr; at the surface 0.075–0.094 / 0.079–0.091.
T_eff 3850–4940 K / 3230–5900 K (noisy in both; it rises with time in both). The moving mesh gives a continuous R★; otherwise
the same evolution.

**Conservation.** Mass |ΔM/M| ≤ 6e-14. TITAN's energy ledger tetot (erg sr⁻¹): t_mm24 3300 → 3999 +1.39e43 in 1.0e7 s
(1.39e36 per s, ramp included); z300g2 3300 → 4399 +2.59e43 in 2.3e7 s (1.13e36 per s): the moving mesh drifts about as much
as the fixed mesh.

**Cost.** After the start-up the moving mesh is not cheaper: z300g2 (8 threads) runs 2.3–3.7 s/step with 0–1 failed solves
per 100 steps and Δt 1e4–4e4 s; t_mm24 (4 threads) 4.4–6 s/step, 3–6 failures per 100, Δt 1.5e3–1.6e5 s (limited by x_H₂ /
x_H⁺ of the zones next to the cliff and by the motion cap when the front jumps).

**Monitor: what is general.** Chosen: ln ρ, ln P, ln T, ln E_r, ln m and ln χρ of the Newton iterate (`GEQW=1,1,1,1,1,0,0,3`),
smooth cap gmax 100, Dorfi–Drury α 1.5, κ 30 crossing times + κ_f 0.3 feature speed; no markers, no budgets, no scales (all
log-derivatives; a power law gives a uniform log mesh), nothing lagged in the residual. It moves zones to whatever structure
appears (accretion shock, photosphere, ionization fronts, the cloud's own gradients) at any accretion rate, and with ~300
zones it puts ~15 % of them near the stellar surface. Composition is left out (operator-split: it would be a lagged input);
its fronts coincide with the χρ and T fronts here. τ is left out (non-local). What it cannot do is resolve a discontinuity:
then gmax is the resolution and the cost knob, and the sub-zone closure must carry the unresolved structure. A natural
general refinement (not implemented): stop refining a zone once its own optical depth on both sides of the front is below
~τ_f (nothing is gained below the mean free path of the transparent side), which bounds the cost where gradients are
unbounded.

**Recommendation (round 2).** The moving mesh works on the new physics (production Δt after a 300-step start-up, smooth R★,
conservation as good as the fixed mesh) but does not resolve the front, and the fixed mesh + SUBZONE (z300g2) is cheaper and
more robust. Keep z300g2 as the main run. GRIDEQ is worth it where the mesh must follow the star over decades in radius (the
swelling to ~100–200 R☉ and the contraction later, where a fixed log mesh of 0.08 puts few zones in the star) — there it
replaces re-meshing; use the iterate-only monitor, GEQRAMP only for restarts from a fixed mesh, SUBZONE on.
