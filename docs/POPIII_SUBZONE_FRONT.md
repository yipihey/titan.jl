# Sub-zone treatment of an unresolved photosphere / recombination front (`SUBZONE`)

Branch `subzone-front` (from `popiii-s6` 0781136), 16 Nov 2026. Test runs: `D/t_sz<N>`, D = /zpool/nvme/data/tabel_scratch/titan_popiii/runs,
each from a `git archive` tree under /zpool/nvme/data/tabel_scratch/titan_popiii/sz_trees, restarted from zdiagB's checkpoint
`z19zdiagBS5` (step 3000) with zdiagB's own `phaseB.env` + `run.env`, `numactl --physcpubind=56-59|60-63 --membind=7`, 4 threads
(launcher: scratchpad `sz/szrun.sh`, one segment as `iceulrun.sh` runs it). Zone numbers below are TITAN indices k (snapshot
index = k − 4; the brief's "zones 62→63" are k = 67 → 68).

## Summary

* **What breaks.** Restarting zdiagB's S5 checkpoint (step 3000) with the unmodified code reproduces segment S6 exactly: Newton
  failures from the first step on, Δt 149 → 0.58 s → 1e-2 s, DTSTOP at step 3032 (t_sz4: 26 Newton failures in 32 steps).
  The front sits between zone 67 (9045 K, χ 21.7 cm²/g, Δτ 4.7e5, table zone) and zone 68 (3187 K, χ 3e-5–5e-4, Δτ ≈ 1,
  irradiated network zone). The flux crossing interface 68 is ∝ 1/χ₆₈ (harmonic mean, §1), and χ₆₈ is steep and lagged:
  in one step at fixed T its x_H⁺ rose ×3.7 and χ ×16 (NLOITDIAG, t_sz16), the flux into zone 68 was choked, E₆₈ fell ×10 in the
  next Newton (t_sz9) and the step failed.
* **SUBZONE** (interface optical depth across an unresolved front, §2) with τ_f = 2 (default) or 1 removes the collapse: Δt
  rises without failures from 2.3 s to 1e4 s within 25 steps after the restart (t_sz18, t_sz23), past step 3032; τ_f = 2/3
  (the bare Eddington–Barbier depth) does not (t_sz17: as the baseline).
* **What remains** with SUBZONE alone: a relaxation cycle of ~20–25 steps (Δt climbs to 1e4–3.6e4 s, then a Newton failure
  drops it to ~20–200 s). In the cycle the irradiated zone's quasi-steady state depends on Δt (T_rad 4840 → 4400 K, χ ×1.4, L
  into zone 69 147 → 95 L☉ at large Δt; back within seconds of model time at small Δt): the signature of an input that is
  held at the start-of-step state. **Freezing the star light** (`SLFREEZE`, diagnostic) **removes the cycle**: SUBZONE + SLFREEZE
  runs from Δt 2.3 s to 2.6e5 s and then holds 4e4–9e4 s (t_sz27, 80 steps, 4 Newton failures all in the first step, 2 error
  rejections). The star-light dilution W and colour temperature T_c of the photo-rates (three-level pumping, H⁻ detachment, LW,
  LOCALJ) are evaluated once per step from the start-of-step E_r, F_r (PopIIIMaterial.starlight_setup via cpl_prepare!) — the
  lagged input that matters. At the front zone itself (τ ≈ 3, SLBLEND weight wq = 1) W ≡ 1 and T_c = T_r(E_r) of the start of step, so the lag that matters
  there is T_c.

* **SLIMP** (§6): the star light of the coupled zones at the Newton iterate. `SLIMP=1` (W only) leaves the cycle unchanged
  (t_sz30: collapse at step 3033 as t_sz23 at 3031). `SLIMP=2` (W and the local part of T_c) removes it: **SUBZONE=1,2 + SLIMP=2**
  (t_sz31) runs Δt 0.58 s → 1e4 s in 30 steps and then 2.6e4–2.1e5 s, 6 Newton failures in steps 3000–3100, while the front
  advances two zones (63|64 → 65|66) — the hand-over that cost SUBZONE alone its Δt is now a factor ≲ 4 dip.

## 1. How the interface optical depth and the exchange are formed now

* **Flux (radiation-momentum) row** at interface k between zones k−1 and k (TitanResidFast, `RADIATION MOMENTUM`): the drag term is
  ½(ρ_k dV_k + ρ_{k−1} dV_{k−1}) · avchi_k · F_k/c with **avchi_k = 2/(1/χ_{k−1} + 1/χ_k)**, the harmonic mean of the two zones'
  flux-mean χ (θ-level, χ = θχ_n + (1−θ)χ_o). The interface optical depth is therefore Δτ_T = (M₋ + M₊)·avchi, M± = ½ρ dV/r^μ.
  With a = χ₋M₋, b = χ₊M₊ the two half-zone depths, Δτ_T = 4ab/(a+b) for equal M and **Δτ_T → 2b(1 + M₋/M₊) when a ≫ b**:
  the opaque half-zone is given the transparent zone's opacity. The same avchi enters the radiation-force work of the gas/total
  energy rows. With VEF_SPH the static terms are r²(q f E)-differences; FLUXTIE replaces the row by F = h c E (upwinded) where
  τ (from the outside) ≤ 0.1, blended to the moment row by τ = 1 (zone 68 sits at τ ≈ 3: weight 0).
* **Exchange** (both energy rows): excap · (c κ_E E − 4π κ_P B) ρ dV with zone-centred θ-level κ_E, κ_P, E, B (one value per zone).
* **VEFIMP formal solution** (TitanCore.eddfac!): piecewise-constant dtau_k = χ_k ρ_k Δr_k per zone, source from κ_P (EDDSCAT: χ − κ_P
  scattering), linear-in-τ source with a harmonic slope limiter; f, g, h of every Newton iterate, held fixed in the Jacobian.
* **Opacity** (PopIIIMaterial `material_opac(::NetGas)`): table zones OPAC(ρ, T); network zones NLOPAC/NLOIT: OPAC(ρ, T) × the
  continuum code's non-LTE/LTE ratio for the zone's composition, the composition from the CK expansion of the coupled species
  with this step's prepared terms (photo-rates of the start-of-step radiation field), steep in T (q = d ln R/d ln T ~ 10–20,
  noted in nlo_it) and in x_H⁺.

At the zdiagB front (step 3000): zone 67 T 9045 K, χ 21.7 cm²/g, Δτ_zone 4.7e5 (table zone); zone 68 T 3187 K, χ 3.1e-5 (table) –
6.1e-5 (NLOIT), Δτ_zone 0.85, T_rad 5235 K (network zone, irradiated: it absorbs ~340 of the 588 L☉ crossing interface 68).
Δτ_T(68) ≈ 1.1 is set entirely by zone 68.

## 2. Candidates and the choice

* (a) interface κ from a reconstructed sub-zone profile between the centres: with κ ∝ T^q (q ≈ 13 across 67→68) any monotone
  profile through the two centre values (log-linear in r: the logarithmic mean; MUSCL/minmod per zone; the series a + b of
  Larsen, Morel & Miller 1987) gives Δτ(67→68) ≈ 3e4 – 1e5 instead of 1.1. That is the depth of a uniform 9045 K half-zone; a
  T(τ) grey atmosphere puts T = 9045 K only τ ≈ 4.6 below the photosphere, i.e. the photosphere sits next to zone 67's centre and
  the profile reconstruction is inconsistent with the zone value. Blocking the flux by 3e4 would heat zone 67 to ~1e5 K: rejected.
* (b) sub-zone reconstruction of T, x within the front zone with κ, emission and exchange integrated over sub-cells: needs the
  NLOIT opacity at sub-zone states, i.e. CK compositions at states that are not iterates (cache thrash, ~3–6 nlo_comp per zone per
  Jacobian refresh already dominate the step); and the lever rule on T (mixture fraction φ = (T − T_c)/(T_h − T_c)) mistakes a
  smooth irradiated gradient for a mixed zone (φ = 0.07 in zone 68 would make it opaque, Δτ 2e4). Rejected for now.
* (c) two-phase zone with the front position as Newton unknown: a tenth unknown per zone (neqn 9 → 10, every row, the CK
  coupling, the error norm). Not attempted.

Chosen: **(a) in the form that is consistent with an unresolved photosphere**: the opaque half-zone contributes its own depth only
up to the depth τ_f ≈ 2/3 from which the emergent radiation originates (Eddington–Barbier relation; Mihalas 1978, *Stellar
Atmospheres*, 1978stat.book.....M), the transparent half-zone its own depth b in series:

    D₋₊ = τ_f a/(a + τ_f) + b,   D₊₋ = τ_f b/(b + τ_f) + a,
    s(x) = x⁴/(x⁴ + τ_f⁴),  ω₋₊ = s(a)(1 − s(b)),  ω₊₋ = s(b)(1 − s(a)),
    Δτ_SZ = (1 − ω₋₊ − ω₊₋) Δτ_T + ω₋₊ D₋₊ + ω₊₋ D₊₋,   avchi_SZ = Δτ_SZ/(M₋ + M₊).

ω is ≈ 1 only at a thick–thin interface; away from fronts the harmonic mean is kept (relative change ≲ (τ_f/min(a, b))⁴ in thick
interiors, ≲ (max(a, b)/τ_f)⁴ in thin envelopes). Rational in χ±, ρ±, r (ForwardDiff Duals carry it, exact Jacobian); the interface
still couples only its two zones. At the zdiagB front d ln Δτ/d ln χ₆₈ drops from 1.0 (TITAN) to 0.27; as zone 68 thickens
(b through τ_f) ω → 0 and the next interface takes over. The front literature for unresolved ionization fronts (photon
conservation: Abel, Norman & Madau 1999, 1999ApJ...523...66A; Mellema et al. 2006, 2006NewA...11..374M) and sub-cell
discontinuity location (Harten 1989, 1989JCoPh..83..148H) were considered for (b)/(c); TITAN: Gehmeyr & Mihalas 1994,
1994PhyD...77..320G. All bibcodes checked on SciX (refereed; Mihalas 1978 is a monograph).

Switches (all default off; off is bit-identical):
* `SUBZONE=1[,τ_f]` — the interface form above (TitanResidFast `SUBZONE`, `sz_tau`, `sz_avchi!`); every PRINTEVERY steps the log
  lists the interfaces with ω > 0.05 (`SUBZONE: n front interface(s) | k: w … dtau/dtau_T …`, k counted from the inner boundary).
* `SZJCHECK=z1,z2` — AD vs central-FD Jacobian of the gas-energy, E_r and F_r rows of zones z−2…z+2 w.r.t. ρ, T, E_r, F_r, x_H⁺,
  x_H₂ of zone z, at the first Newton iteration (TitanResidFast `rad_jcheck`).
* `OLDOPAC=1` — after the COUPLED prepare pass the start-of-step opacities are re-evaluated with the material the Newton uses, so
  the θ-scheme's old level (chifo, xkeo, xkpo) and the new level are one opacity law (see §4).
* `SLFREEZE=1` — diagnostic only (not physics): every zone keeps its first star-light W, T_c.

## 3. Diagnosis (step by step, baseline code restarted from S5)

| run | what | result |
|---|---|---|
| t_sz4 | baseline restart, 60 steps | = S6: DTSTOP at step 3032, 26 Newton failures / 32 steps, all E_r-limited |
| t_sz2, t_sz9 | NEWTDBG/NEWTKX of steps 3001, 3002 | step 3001 (Δt 149) converges (TR 18 + BDF 8 iterations), E₆₈ 5.7 → 13.1; step 3002: zone 68 starts at χ 4.9e-4 (was 3.1e-5), Newton drives E₆₈ 13 → 1.3, diverges; final contraction of the converging TR stage ≈ 0.11 per iteration (linear, not quadratic) |
| t_sz16 | NLOITDIAG, CMP of zones 67–69 | x_H⁺(68) 8.4e-9 → 3.1e-8, χ_NLOIT(68) 4.9e-4, κ_P 1.1e-3 after step 3001 at unchanged T |
| t_sz25 | ZDIAG zones 66–70 | zone 68 per gram over a 149 s step: exchange 1.06e9, RAYQ 1.2e6, line cooling Λ_ext 5.5e5, GASQ 80 erg/g: the lagged RAYQ/GASQ/line sources are ≤ 0.1 % |
| t_sz3 | LOCALJ=0 | same collapse |
| t_sz5 | NLOIT=0 | collapse to Δt ~ 1–3 s (instead of 1e-2) |
| t_sz19 | OLDOPAC=1 | same collapse |
| t_sz14 | SLFREEZE=1 (alone) | same collapse |
| t_sz7, t_sz17 | SUBZONE=1, τ_f = 2/3 | same collapse |
| t_sz20 | SUBZONE=1,1 | no collapse; Δt cycles 2–1.4e3 s (12 Newton failures / 60 steps) |
| t_sz18, t_sz23 | SUBZONE=1,2 | no collapse; Δt cycles 20 s – 3.6e4 s (13 / 60 steps) |
| t_sz27 | SUBZONE=1,2 + SLFREEZE=1 | Δt 2.3 s → 2.6e5 s, then 4e4–9e4 s; 4 Newton failures (first step) + 2 error rejections / 80 steps |

### Lagged inputs of the front zones' residual (PI question, 16 Nov)

Held over the whole step (start-of-step state, constant in the Newton):
1. **Star-light photo-rates** of the coupled chemistry and of the NLOIT opacity composition — W (from c(E_r − E_cmb), F_r: SLWFLUX),
   T_c (SLBLEND/SLEB: T_eff at the τ = 2/3 zone blended with T_r at τ + 2/3), kHm, kH, kLW, the three-level parameters β, P12, RL
   (Lyα/two-photon pumping ∝ W), LOCALJ's J = W B(T_c): `starlight_setup` → `chem_core!` prepare pass → `CPLTERMS[k]` (env, phc,
   extra), PopIIIMaterial `cpl_prepare!` (called by popiii_m2c's STEPFN before every step). **Matters** (t_sz27 vs t_sz23).
2. NLOPACJ: κ_E weighted with B_ν(T_c) of the same lagged T_c (PopIIIMaterial nlte_means `Tw`). Part of 1.
3. GASQ (LW heat), RAYQ (three-level photon energy of non-coupled zones), CHEMB: last chemistry (PopIIIMaterial chem!), energy rows
   TitanResidFast lines `isempty(GASQ[])`/`RAYQ[]`. ≤ 0.1 % of zone 68's exchange: **does not matter here**.
4. LINEIMP: H₂ line cooling of the iterate's populations with the radiative (escape/ALI) rates of the last transfer call (popiii_m2c
   `li_store!`, EP/ALI). Λ_ext ≤ 5e-4 of the exchange in zone 68: **does not matter here**.
5. The dual-energy row choice (TitanResidFast `dualswitch`, old state; ETASWTIE uses TIEW). Discrete, not a source; not tested.
6. The θ-scheme old level (chifo, xkeo, xkpo, plfo, Compton old) — by construction. After a restart the old opacities are the
   table's (NLOST not yet rebuilt); OLDOPAC re-evaluates them with the Newton's material: **no effect** (t_sz19).

Evaluated at every iterate but without derivatives (the Newton sees them as constants):
7. VEFIMP Eddington factors f, g and the flux tie's h and weights (TitanCore.eddfac!, TD.HEDD, TD.TIEW): f₆₈ moves 0.48 → 0.39
   between iterates/steps. The TR stage converges linearly (×0.11 per iteration) — consistent with this. An exact ∂f/∂E is a
   dense (formal-solution) dependence that the banded Newton cannot hold; not tested.
8. EXCAP factor (`excap`, value only) — the one Jacobian mismatch SZJCHECK finds (zone 67, ρ column, row E: 0.7 % of the
   1e-6-floored row scale), identical with SUBZONE off.
9. NLOIT's T, x, ρ slopes (finite-difference ratios refreshed when the iterate moves by more than NLOITJ) and CPLJLAG's reused local
   Jacobians: approximate derivatives, exact values.

## 4. Tests

1. **Bit-identity, switches off** (SUBZONE, OLDOPAC, SLFREEZE, SZJCHECK unset), against `popiii-s6` 0781136 (`sz_trees/base`):
   * restart of z19zdiagBS5, 32 steps to the DTSTOP at 3032, 4 threads: t_sz24 (96c895e) vs t_sz4 — **all 591 state arrays and
     the .aux identical** (also t_sz8, 6bb247e);
   * ICEUL start (abn02red, 128 IC shells), 6 steps: t_sz15 (09b7414) vs t_sz11 — **all 591 arrays and the .aux identical**
     (code after 09b7414 changes only default-off paths and the τ_f default).
2. **Jacobian vs finite differences.** `julia/subzone_test.jl` (TITAN_NOLIB=1): ∂ln avchi/∂ln(χ, ρ) of a 4-zone strip by AD vs
   central FD, worst absolute difference 6e-10 – 1.2e-9 for front / hand-over / thick interior / thin envelope; limits away from
   fronts |Δτ_SZ/Δτ_T − 1| ≤ 4e-7. In the run (`SZJCHECK=62,63,64,65`, t_sz21, SUBZONE=1,2 at S5): rows E, F, gas energy of zones
   64–71 w.r.t. ρ, T, E, F, x_H⁺, x_H₂ of zones 66–69: worst 1.1e-5 (F row vs ρ: FD truncation), except the EXCAP entry above (same
   with SUBZONE off, t_sz22).
3. **From S5 with the switch on** (t_sz23, iceulrun.sh, 4 threads, still running at hand-back): past step 3032 without collapse;
   steps 3000–3213: Δt median 1.8e3 s, 63 of 213 steps at Δt ≥ 1e4 s (maxima 1.2e4–4.3e4 s in each ~20–25-step cycle, minima
   ~20–200 s); Newton failures (TR stage) 19 in steps 3000–3100 and 7 in 3100–3200 (baseline S5 segment 11/100 at Δt ~150 s;
   baseline restart 26 in 32 steps before DTSTOP); 7.8–8.8 s/step at 4 threads (perf.csv). M★ 0.0611 → 0.0634 M☉ by step 3185 (t 2.38933476e5 → 2.38933505e5 yr), R★ 11.003 R☉ (the zone edge, as
   before), L_ph 40–115 L☉ and T_eff 2600–3270 K over the cycle against 57–94 L☉ and 2520–3100 K in the baseline steps 2700–3000;
   the steps both runs cover (3000–3032) the baseline spends at Δt < 1e2 s, so only the state at 3000 is common.
   Front weight in the run: two interfaces with ω > 0.05 (68: ω 0.82–0.86, Δτ/Δτ_T 0.61 at τ_f = 1; 69: ω ≈ 0.1).
4. **No front, no effect.** ICEUL start with SUBZONE=1 (t_sz28) vs off (t_sz11), 6 steps: no interface with ω > 0.05; all
   physical state arrays identical, only the Jacobian work array e00 differs (3.9e-16 relative). Unit test: thick interiors
   |Δτ_SZ/Δτ_T − 1| ≤ 1e-11, thin envelopes ≤ 4e-7 (τ_f = 2).

## 5. Open issues

* Star light: done as SLIMP (§6). Still lagged: the non-local part of T_c (T_eff at the τ = 2/3 zone, T_r at τ + 2/3; enters
  through T_src^(1−wq)), the LW shielding columns, and the GASQ/RAYQ bookkeeping of non-coupled zones.
* τ_f is a closure constant (depth of the opaque zone's centre below its unresolved surface): 2/3 fails here, 1 and 2 work; a
  self-consistent τ_f (from T_h/T_eff of the grey T(τ) relation) was not tried. It also sets the front window; separating the
  two roles may allow the physical 2/3.
* The interface treatment does not address the exchange of the irradiated zone (its κ_E, κ_P remain zone-centred) nor the
  VEF formal solution's piecewise-constant dtau.
* The 556 M☉ warm front of the production run (zones ~98–106) was not tested.

## 6. Star light at the Newton iterate (`SLIMP`, 16 Nov, PI request)

`SLIMP=1|2` (with COUPLED; default off, off is bit-identical). Per coupled zone k, with the start-of-step star light recorded by the
COUPLED prepare pass (`starlight_setup` → `SLW0`):

* q_k = W_k(iterate)/W_k(start) and t_k = T_c,k(iterate)/T_c,k(start) from zone k's E_r and its two interface fluxes (`sl_qt`):
  SLBLEND T_c = T_src^(1−wq) T_r^wq with T_r = (max(E_r, E_cmb)/a)^¼ (T_src, wq held), W = 1 (wq ≥ 1) or min(c u/(4σT_c⁴), 1),
  c u = max(c(E_r − E_cmb), ½(F_k + F_k+1)) (SLWFLUX). Local (zone k and interface k+1): the block structure holds. At the start
  state q = t = 1 exactly. `SLIMP=1`: t ≡ 1.
* Every star-light rate is linear in W at fixed T_c, so q scales the prepared terms exactly (`cpl_tq`): three-level β, P12 (reactions,
  Peebles C, heats, cooling closure), the network photo-rates of J = W B_ν(T_c) (R → qR, α → α₀ + q(α − α₀), α₀ at W = 0), the
  star-light `extra` (kHm, kH, kLW). With t ≠ 1 the same terms are rebuilt at T_c = t T_c,start (β₂ₛ(T_c), n̄_α, ⟨hν⟩, photo-rates,
  kHm, kH, kLW, NLOPACJ's weighting T_c).
* q, t are local variables of the coupled zone (`cpl_zone_q`: z = ρ, T, x_H⁺, x_H₂, q, t; 7×6 local Jacobian, AD in q, central
  differences in ln t), of the NLOIT opacity (value at (q, t) exact, slopes in ln q, ln t by differences as NLOIT's other slopes)
  and of the post-step species (`cpl_finish!` at the solved E_r, F_r).

Tests:
* Bit-identity off (d5d63fb): restart of S5, 32 steps to DTSTOP (t_sz34 vs t_sz4) and ICEUL 6 steps (t_sz35 vs t_sz11): all 591
  arrays and the .aux identical.
* SZJCHECK (t_sz29 SLIMP=1, t_sz32 SLIMP=2, t_sz33 off; zones 67–69, rows gas energy, E_r, F_r, x_H⁺, x_H₂): the E_r and F_r columns
  of the species rows agree to ≤ 9e-7 (SLIMP=1) and 6e-5 / 9e-5 (SLIMP=2, zone 68/69, the T_c column by differences). The
  species-row entries of the ρ, T, x columns differ by 1e-4–2e-3 identically with SLIMP off (pre-existing: NLOIT slopes / CPLJLAG).
* From S5 with SUBZONE=1,2:

| run | steps | Δt | Newton failures | s/step (4 thr) |
|---|---|---|---|---|
| t_sz23, SLIMP off | 3000–3396 | median 3.1e2 s, 63 of 396 ≥ 1e4 s; cycles 20 s ↔ 1e4–6e4 s; after the front moved (3195) ≤ 3.6e3 s | 19, 7, 15 per 100 | 7.4–7.8 |
| t_sz30, SLIMP=1 | 3000–3068 | same cycle (collapse at 3033, 3049) | – | – |
| t_sz31, SLIMP=2 | 3000–3105+ | 1e4 s from step 3030, then 2.6e4–2.1e5 s | 6 in 3000–3100 | ~9.8 |
| t_sz27, SLFREEZE (diagnostic) | 3000–3080 | 4e4–2.6e5 s | 4 + 2 rejections / 80 | ~6.8 |

* Physics (runrec scalars): at equal model time the two runs agree — t ≈ 238933.54 yr: M★ 0.0667 (t_sz31, interpolated) vs 0.0669
  (t_sz23), R★ 11.9 R☉ in both; T_eff 2900–3700 K and L_ph 60–140 L☉ in both. t_sz31 covers 0.17 yr in 105 steps (M★ 0.0611 →
  0.0755, R★ 11.0 → 12.9 R☉), t_sz23 0.065 yr in 396 steps.

Table/network boundary (moving-mesh agent's question). The front is the last equilibrium-table zone | first network zone in every
snapshot of t_sz23 and t_sz31 (zdiagB: zone 67 table, 68 network). t_sz23: 1 switch in 396 steps (3195); the cycle runs without
switches. t_sz31: 2 switches in 105 steps (3055, 3080), each followed within ~5 steps by the largest Δt dips of the run (1.7e5 →
4.3e4, 1.0e5 → 2.6e4 s): with the star light at the iterate the remaining Δt limitation coincides with the per-step table↔network
switch of the zone the front enters. (Sampling: snapshots every 5 steps.)

VEF f, g, h derivatives (linear TR-stage contraction ×0.11). f_k, g, h_k come from the formal solution (TitanCore.eddfac!, short
characteristics through every zone): ∂f_k/∂(T_j, E_j, ρ_j) is dense, which the block-pentadiagonal coloured-AD Jacobian (ncol 6)
cannot hold. Feasible variants: (i) a truncated local derivative — coloured finite differences of eddfac (perturb every 6th zone's
E_r and T together, 6 colours × 2–3 variables ≈ 12–18 extra formal solutions per Jacobian) entering only the k−2…k+2 blocks;
cheap, approximate, most of the dependence in the thin zones where f ≠ 1/3 is local; (ii) exact: Jacobian-free Newton–Krylov
(GMRES on the residual with f at the iterate, the present banded Jacobian as preconditioner), one formal solution per Krylov
iteration — a solver change. Not attempted.
