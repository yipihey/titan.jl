# Interface radiation flux from the steady slab solution between zone centres (`SLABFLUX`)

Branch `slabflux` (from `popiii-s6` ac7c73c), 16 Nov 2026. Code: julia/TitanResidFast.jl (`SLABFLUX`, `sf_E`, `sf_mean`,
`sf_avchi!`), julia/popiii_m2c.jl (switch, log line), unit test julia/slabflux_test.jl. Test runs `D/t_sf<N>`,
D = /zpool/nvme/data/tabel_scratch/titan_popiii/runs, from `git archive` trees in /zpool/nvme/data/tabel_scratch/titan_popiii/sf_trees,
4 threads each, `numactl --physcpubind=32-35|36-39 --membind=4`.

## Summary

* **Formula.** The steady first moment across the dual cell between two zone centres (no source inside: the exchange is in the zone
  rows) has F constant, F = c Δ(q f E)/Δτ, so the slab problem is the optical depth Δτ between the centres. For a grey layer with
  χ ∝ (τ + τ₀)^n (exact for a grey LTE layer at constant flux with χ ∝ T^q, n = q/4) through the two centre values,
  **Δτ = ΔM χ₊ ℰ(x/n)/ℰ(−x(n−1)/n)**, x = ln(χ₋/χ₊), ℰ(y) = (e^y − 1)/y; n = 1 is log-linear in r (logarithmic mean), n = 2 the
  geometric mean, **n → ∞ (default) the exponential in optical depth, χ̄ = ln(χ₋/χ₊)/(1/χ₊ − 1/χ₋)** — the Scharfetter–Gummel
  Bernoulli form. Symmetric, between the harmonic and the geometric mean, = χ for χ₋ = χ₊ (diffusion limit / series a + b), ≤ max χ
  (thin limit untouched: free streaming stays with VEF f and FLUXTIE h), C∞, two-zone local, E-independent (diffusive). `SLABFLUX=1`.
* **What it does and does not fix.** The photospheric fluctuation of the 100-zone runs is a switch between a bright state
  (L ≈ 2e3 L☉, T_eff ≈ 5000 K) and a dark state (L ≈ 0.3–20 L☉, T_eff ≈ 550 K) in which the hot front zone (9350 K, Δτ 3e4) sits under
  an interface that carries < 1 L☉. The choke is the **FLUXTIE donor cell** (h = H/J at the donor's centre, 1e4 optical depths
  below the interface in the formal solution's piecewise-constant picture: h ≈ 3e-5), not the interface opacity. SLABFLUX alone
  reduces the dark time (89 % → 36 % moving mesh, 66 % → 30 % fixed) but not the switching. **`SLABTIE=1`** applies the same slab
  profile to the tie's donor (§2b: the donor centre lies τ_i ≈ 0.1 below the interface; h = h_FS + f/(τ_i + 2f) − f/(a + 2f),
  exactly h_FS without a front) and removes the dark state: L = 1.7–3.9e3 L☉, T_eff 5310 K ± 0.7–1 % (baseline 650 K, ± 2 % moving /
  26 % fixed), on both meshes, with the M★ growth unchanged.
* **Cost.** SLABFLUX+SLABTIE: 20–22 Newton failures/100 steps (baseline 0.1–1.5), 2.4 s/step (baseline 1.4) at 4 threads; the failed
  solves converge slowly in E_r and the species of the zones at and just above the photosphere, which the dark baseline never
  exercises. Median Δt 3e3–4.5e3 s (baseline 6e3).
* **Bit-identity off:** z300g2 (290 zones) and z100mm (100 zones), 30 steps each: all 591 state arrays and the .aux identical.
* **Recommendation:** for the emergent radiation, `SLABFLUX=1 SLABTIE=1 SUBZONE=0` (SUBZONE's τ_f cap is not needed and raises the
  front Δτ 90–190× above TITAN's in z100mm); before production, the Newton convergence in the bright photosphere needs work (§4.5).

## 1. The problem

The photosphere of the accreting protostar is a transparency front thinner than a photon mean free path: across one interface χ
drops by 1e3–1e5 (H⁻/H bound-free, χ ∝ T^q with q ≈ 10–20) and the optical depth of the hot half-zone is 1e3–1e5 while the cold
neighbour's is ≲ 1. TITAN's radiation-momentum row (TitanResidFast, `RADIATION MOMENTUM`) at interface k is, in the static part,

    r_k^μ [q f E]_{k−1}^{k} + ΔM_k χ̄_k F_k/c = 0,      ΔM_k = M₋ + M₊,  M± = ½ ρ± dV±/r_k^μ   (column of the two half-zones)

with χ̄_k = avchi_k = 2/(1/χ₋ + 1/χ₊), the harmonic mean (the same avchi enters the radiation force of the energy rows). Across the
front χ̄ → 2χ₊: the opaque half-zone is given the transparent zone's opacity, and the flux leaving the star follows the steep,
irradiated, lagged χ of the cold zone one to one (d ln Δτ/d ln χ₊ = 1). SUBZONE (docs/POPIII_SUBZONE_FRONT.md) caps the opaque
half-zone's depth at τ_f; in the live z100mm run its Δτ is 88–194× TITAN's at the front interface (z300g2: 7.6–331×; the runs' SUBZONE log lines, 16 Nov).

## 2. Derivation

**Finite-volume statement.** Take the dual cell between the two zone centres. The unknown F_k is the interface flux; the zone values
(E, T, χ, f) are point values at the centres. Write the grey moment equations in the column mass m (dm = ρ dr) and τ (dτ = χ dm):

    d(q f E)/dτ = −q F/c,                    (first moment; q = Auer's sphericality factor, VEF_SPH)
    dF/dτ = c (κ_E/χ)(E − 4πB/c) + geometry  (zeroth moment)

(Mihalas & Mihalas 1984, *Foundations of Radiation Hydrodynamics*, 1984frh..book.....M; Auer's sphericality factor as
used by TITAN's VEF form, Gehmeyr & Mihalas 1994, 1994PhyD...77..320G). In the finite-volume split the absorption/emission of each
half of the dual cell belongs to its zone's energy rows (they hold the exchange with zone-centred κ_E, κ_P, E, B). The steady
first moment across the dual cell is then solved exactly with F = const:

    F_k = c [ (q f E)₋ − (q f E)₊ ] / Δτ_k,      Δτ_k = ∫ χ dm   (centre to centre).

The sphericality part (q) is already exponentially fitted in TITAN's VEF_SPH form (vefq: the exact ratio of q between the centres for
constant f, an exact solution of the χ = 0 problem). What is left of "the analytic solution across the slab" is the optical depth
between the centres for an opacity profile through the two centre values χ₋, χ₊ over the column ΔM. The error of F = const is the
second moment of the exchange inside the dual cell, ∫(τ − τ_k) dF/dτ dτ; at the zdiagB front (the cold zone absorbs ~60 % of the flux
over its whole depth, its half-depth is ≲ 0.4 of a total Δτ of ≳ 7) it changes the relation by ~1 %, below the closure uncertainty.

**Opacity profile.** χ must be modelled between the centres. The natural one-parameter family is the grey LTE layer at constant flux:
with E = (3F/c)(τ + τ₀) (Eddington) and χ ∝ T^q = E^{q/4}, χ ∝ (τ + τ₀)^n with n = q/4. Given χ₋, χ₊ (ratio R = χ₋/χ₊, x = ln R)
and ΔM, τ₀ and the constant follow, and

    Δτ = ΔM χ_n,      χ_n = χ₊ ℰ(x/n) / ℰ(−x (n−1)/n),      ℰ(y) = (e^y − 1)/y.

Derivation: s± = τ± + τ₀, s₋/s₊ = R^{1/n}; ΔM = ∫dτ/χ = (s₊/χ₊ − s₋/χ₋)/(n − 1); Δτ = s₋ − s₊ = s₊(R^{1/n} − 1); eliminate s₊:
Δτ = ΔM χ₊ (n − 1)(R^{1/n} − 1)/(1 − R^{1/n − 1}), which is the ℰ form above (written so that it is smooth through R = 1). Special
members:

| n | profile | χ_n |
|---|---|---|
| 1 | χ linear in τ, i.e. **exponential in column/radius** (log-linear in r) | logarithmic mean (χ₋ − χ₊)/ln R |
| 2 | χ ∝ (τ + τ₀)² | geometric mean √(χ₋χ₊) |
| q/4 | grey LTE layer with χ ∝ T^q | exact for that layer |
| ∞ | **exponential in optical depth** | χ₊ x/(1 − e^{−x}) = ln R/(1/χ₊ − 1/χ₋) |

The n → ∞ member is the Scharfetter–Gummel form: the Bernoulli function B(x) = x/(e^x − 1) appears exactly as in the exponentially
fitted drift-diffusion flux (Scharfetter & Gummel 1969, IEEE Trans. Electron Devices 16, 64, 1969ITED...16...64S), here because the
diffusion coefficient c/(3χ) varies exponentially across the cell. All χ_n are symmetric in χ±, lie in [min χ, max χ], decrease with
n, and H ≤ χ_∞ ≤ G = χ_2 ≤ χ_1 (H harmonic, G geometric; unit test 1).

**Why not log-linear in r (n = 1), as first proposed.** For the zdiagB front (χ 21.7 | 3.1e-5, cold half-zone depth 0.42) n = 1 gives
Δτ ≈ 3e5 and n = 2 ≈ 700: the front would sit in the middle of the slab, the hot zone's centre ~1e5 optical depths below it. A grey
T(τ) atmosphere puts the 9045 K centre only ~5–12 optical depths below the photosphere (T⁴ = ¾T_eff⁴(τ + ⅔) with T_eff 3000–5200 K),
which is the SUBZONE doc's reason (§2a) for rejecting profile reconstructions. The family resolves this: the implied front position
(as a fraction φ of the slab from the opaque centre, for a sharp step Δτ = ΔM(φχ₋ + (1−φ)χ₊)) is φ = 1.4e-6 (harmonic), 1.8e-5 (n = ∞),
2.5e-4 (n = 3), 1.2e-3 (n = 2), ~0.1 (n = 1) (unit test 2b): only large n puts the photosphere next to the hot zone's centre, i.e.
treats the hot zone value as a point a few optical depths below an unresolved surface. n = ∞ gives Δτ = 11.3 there, the grey-atmosphere
estimate is 5–12. Default: n = ∞.

**Kirchhoff transform, rejected.** For the LTE layer the exact flux is the Kirchhoff form F = (c/3ΔM)∫dE/χ(E). Fitting χ(E) as a power
law through the two centre values (exponent from the E's) makes F_k depend on E₊ with ∂F_k/∂E₊ > 0 at a front (flux into the cold
zone grows with its own energy: anti-diffusive, unstable in the Newton). χ_n uses χ± only; n is a fixed closure constant.

**Limits (requirements).**
* Thick, smooth χ (both half-zones ≫ 1, R → 1): χ_n = χ(1 + O(ln R)) → F = −(c/3χρ)∂E/∂r with f = 1/3: the diffusion limit, and the
  scheme reduces to TITAN's cell-edge diffusion discretisation (the asymptotic diffusion limit in the sense of Larsen, Morel & Miller
  1987, 1987JCoPh..69..283L; equal χ gives exactly the series a + b of the two half-zones).
* Thin (Δτ ≪ 1): χ_n ≤ max χ, so the drag term stays negligible against r^μ Δ(q f E); free streaming is then carried, as before, by the
  VEF factors f, g from the formal solution (f → 1) and FLUXTIE's h (F = h c E). Nothing changes there.
* Front: Δτ = ΔM χ₊ ln R/(1 − 1/R) ≈ (ln R) × (twice the cold half-zone depth): 4.15–23 in the tests against TITAN's 0.2–4.
  The sensitivity d ln Δτ/d ln χ₊ = 1 − 1/ln R ≈ 0.90–0.93 (n = ∞), 1 − 1/n ≈ 0.66 (n = 3): only slightly below TITAN's 1.0.
* Smooth: C∞ in χ±; ℰ by its Taylor series (to y⁶) for |y| < 1e-2 (value, first and second derivative continuous to 2e-16, 7e-15,
  3e-12); ForwardDiff Duals carry it, so the Newton gets exact derivatives through χ(ρ, T, x, E…) as before.
* Local: χ_k−1, χ_k only: the block structure is unchanged.
* Spherical geometry / VEF: the row is unchanged except avchi; q, f, g, h and FLUXTIE as before.

**Short characteristics.** The Olson & Kunasz (1987, 1987JQSRT..38..325O) short-characteristic solution with a source linear (or
exponential) in τ across the cell gives the intensity at the interface; its moment F = 2π∫Iμ dμ is the transport analogue of the above.
In the moment system the corresponding statement is the first-moment relation with F constant: the source enters only through the zone
exchange, so the linear-source part of OK87 contributes the same ∫(τ − τ_k)dF/dτ term discussed above. The VEF formal solution itself
(TitanCore.eddfac!) still uses piecewise-constant dtau per zone (unchanged).

**Multigroup.** `sf_mean(χ₋, χ₊, n)` depends only on the two zones' opacities of the quantity being transported, so per group g it is
`sf_mean(χ_g,k−1, χ_g,k, n)` with the group's flux-mean opacity; the source/absorption of the group stays in its zone energy rows
(F_g constant across the dual cell). `sf_avchi!(avchi, chif, ngrs, ngre)` takes the opacity array as an argument.

Sources (all checked on SciX, 16 Nov 2026; refereed unless a book): Olson & Kunasz 1987 JQSRT 38, 325 (1987JQSRT..38..325O);
Scharfetter & Gummel 1969 IEEE TED 16, 64 (1969ITED...16...64S); Larsen, Morel & Miller 1987 JCoPh 69, 283 (1987JCoPh..69..283L);
Mihalas & Mihalas 1984 (1984frh..book.....M, book); Gehmeyr & Mihalas 1994 (1994PhyD...77..320G).

## 2b. The same slab argument at the flux tie (`SLABTIE`) — where the photospheric fluctuation comes from

The first SLABFLUX runs (t_sf8, t_sf9, §4) and the baseline both alternate between a **bright** state (L ≈ 1.7e3 L☉ through
100 R☉, T_eff ≈ 4900 K, the cold zone above the front at T_rad ≈ 4200 K) and a **dark** state (L ≈ 0.3–20 L☉, T_eff ≈ 540 K) in which
the hot front zone sits at T = T_rad ≈ 9350 K and its outer interface carries 0.3 L☉ (t_sf8 step 5805: zone 22 Δτ 2.7e4, χ 15.7;
zone 23 Δτ 0.008, χ 7.3e-6, T_rad 900 K). A 9350 K surface at 58 R☉ should radiate ~2e4 L☉. The interface is in the FLUXTIE region
(τ from the outer boundary 0.034 ≤ 0.1, weight 1), where the row is F_k = h c E upwinded with h = H/J of the formal solution **at the
donor zone's centre**. The formal solution (TitanCore.eddfac!, piecewise-constant dtau per zone) places that centre 1.4e4 optical depths
below the interface, so h₂₂ ≈ 1/(3 × 1.4e4) ≈ 2.5e-5 and F = h c E₂₂ gives ~1 L☉: the flux is choked by the donor-cell tie, not by the
interface opacity — SLABFLUX alone cannot remove the dark state. The bright state is the same configuration with a front zone of
Δτ ≈ 1 (h ≈ 0.3).

SLABTIE applies the slab solution to the donor: in SLABFLUX's exponential-in-τ profile the donor centre lies only
τ_i = Δτ_SF ln(1 + (R − 1)M₋/ΔM)/ln R below the interface (0.12 at the dark front, against the piecewise-constant 9e3). The steady first
moment from the donor centre to the interface (fE − f_s E_s = F τ_i/c) with a Marshak/Eddington free surface (F = c E_s/2; Mihalas &
Mihalas 1984) gives h_slab(τ) = f/(τ + 2f); the tie's outward donor factor becomes

    h = max(h_FS + h_slab(τ_i) − h_slab(a), 0),     a = χ₋M₋ (the formal solution's piecewise-constant half-zone depth),

i.e. the formal solution's h corrected by the change of the slab estimate when the donor centre is moved from depth a to τ_i. No front
(χ₋ = χ₊): τ_i = a, h = h_FS exactly; thick slab: 1/(3τ_i + 2) (diffusion); τ_i → 0: ½. At the dark front h = 0.43 instead of 3.7e-5.
Local (zones k−1, k), smooth in χ±, ρ±, r (the max only at h = 0). The inward term of the tie is unchanged.

## 3. Implementation

* `SLABFLUX=1[,n]` (popiii_m2c; default off; n omitted or 0: n = ∞; n ≥ 1). After TITAN's harmonic avchi is formed, `sf_avchi!`
  replaces avchi_k of the interior interfaces ngrs+1..ngre by χ_n(χ_{k−1}, χ_k) (θ-level χ). The boundary interfaces keep TITAN's.
* With `SUBZONE` also on, SUBZONE's front form is applied on top (its Δτ_T is then the slab Δτ; at a front ω ≈ 1 and SUBZONE's
  D = τ_f a/(a + τ_f) + b wins). With `SUBZONE=0` (the stack default is `SUBZONE=1,2`) SLABFLUX alone.
* `SLABTIE=1` (needs `SLABFLUX=1`, n = ∞; default off): the FLUXTIE row's outward donor factor of §2b (TitanResidFast `st_tauint`,
  `st_hdonor`); everything else in the tie (weights TIEW, inward term, r² dilution) unchanged.
* Log (every PRINTEVERY steps): `SLABFLUX: n interface(s) with dtau/dtau_T outside [2/3, 3/2] | k: ratio …` (k from the inner
  boundary).

## 4. Tests

All runs restart z100mm's checkpoint `z19z100mmS2_auto0` (step 5000, t = 237828.078 yr, M★ 0.281 M☉) with z100mm's own phaseB.env
+ run.env, 1000 steps, SNAPEVERY=5, 4 threads; "fixed" = the same checkpoint with `GRIDEQ=off` (PI, 16 Nov: 100-zone tests only).
Analysis scripts (scratchpad sf/): ana.jl (photosphere = zone where τ from outside reaches 2/3, L_ph = 4πr²F at its outer interface;
L(100 R☉) through the interface nearest 100 R☉; T_eff of the star-light record; rms ln about an 11-snapshot running median; perf.csv),
series.jl, prof.jl, failz.py; bit compare audit/ckcmp.jl.

### 4.1 Unit tests (julia/slabflux_test.jl, TITAN_NOLIB=1): all passed

1. Limits/ordering: χ_n = χ for χ₋ = χ₊ (≤ 4e-16); n = 1 = logarithmic mean, n = 2 = geometric mean (1.9e-15); H ≤ χ_∞ ≤ G ≤ χ_1,
   symmetric, decreasing in n, within [min χ, max χ] on a 21 × 21 grid χ± ∈ 1e-6…1e4; ℰ series branch at |y| = 1e-2: jumps of value,
   first and second derivative 1.5e-16, 6.7e-15, 3.2e-12.
2. (a) Fine-grid numerical solution of the two-zone slab (grey LTE layer at constant flux, χ ∝ T^q across a transparency front,
   χ ratio R up to 1e5, shooting on 2e5 graded RK4 cells): fine grid = Kirchhoff closed form to 1e-14; SLABFLUX with n = q/4 exact
   (≤ 9e-15). Flux predicted from the centre values alone, relative to the fine grid:

   | q, R, cold half-depth b | Δτ fine | harmonic (TITAN) | SLABFLUX n = ∞ | geometric | log mean (n = 1) |
   |---|---|---|---|---|---|
   | 13, 1e3, 0.3 | 10 | F × 8.4 | × 2.4 | × 0.53 | × 0.12 |
   | 13, 1e5, 0.3 | 45 | × 38 | × 6.6 | × 0.24 | × 0.009 |
   | 10, 1e5, 1.0 | 297 | × 74 | × 13 | × 0.47 | × 0.017 |
   | 8, 1e4, 0.05 | 10 | × 50 | × 11 | × 1.0 | × 0.09 |

   TITAN's harmonic mean overestimates the flux out of such a layer by 8–74×, n = ∞ by 2.4–13×; the exact member needs q.
   (b) Sharp step (zdiagB 67|68 values): the front position each closure implies (fraction of the slab from the opaque centre):
   harmonic 1.4e-6, n = ∞ 1.8e-5 (Δτ 11.3), n = 3 2.5e-4 (Δτ 148), geometric 1.2e-3 (Δτ 703).
3. AD vs central FD of avchi on a 4-zone strip (front, hand-over, thick, series branch, thin; n = ∞ and 3): worst 1.4e-9 in the
   elasticities. Front sensitivity d ln χ̄/d ln χ₊: 0.90–0.93 (n = ∞), 0.66 (n = 3), TITAN 1.0.
4. SLABTIE: τ_i = χM₋ for χ₋ = χ₊ exactly, series branch continuous (3e-9); dark front of t_sf8 step 5805: τ_i 0.117 (piecewise-constant
   half-depth 9.1e3), h 3.7e-5 → 0.43; AD vs FD of h 2.8e-11; no front: h = h_FS exactly; thin envelope gradient (χ ratio 1.3): 4e-6.

### 4.2 Bit-identity with the switches off

* z300g2 checkpoint `z19z300g2S5_auto1` (290 zones, step 12900), 30 steps: t_sf3 (popiii-s6 ac7c73c) vs t_sf4 (5854839, SLABFLUX
  code): **591 / 591 state arrays and the .aux identical** (raw bytes).
* z100mm checkpoint (100 zones, moving mesh), 30 steps: t_sf13 (ac7c73c) vs t_sf14 (7cbe6e8, SLABFLUX + SLABTIE code): **identical**.

### 4.3 Jacobian in the run (SZJCHECK=19,20,21,22,23 at the z100mm front, first Newton iteration)

t_sf15 (SLABFLUX + SLABTIE) vs t_sf16 (off, SUBZONE=0): AD vs central FD of the gas-energy, E_r, F_r, x_H⁺, x_H₂ rows of zones z ± 2
w.r.t. ρ, T, E_r, F_r, x_H⁺, x_H₂ of zone z. The F_r column (the one the switches change): worst 1.9e-7 – 5.8e-6 on (2.2e-7 – 7.9e-6
off). All other columns are identical on and off: the pre-existing mismatches (ρ column of the E_r row = the EXCAP factor, value only,
docs/POPIII_SUBZONE_FRONT.md; T column of the species rows, 1e-3 – 3e-2 = NLOIT slopes) are unchanged.

### 4.4 Restarts, 1000 steps

| run | mesh | switches | L(100 R☉) median [L☉] | dark (L < 100) | T_eff mean (min–max) [K] | rms ln T_eff | rms ln L_ph | failures /100 | s/step | Δt median [s] | t covered [yr] | M★ end [M☉] | R★ end [R☉] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| t_sf10 | moving | stack (SUBZONE=1,2) | 22 | 89 % | 651 (526–1731) | 0.023 | 0.09 | 0.1 | 1.37 | 6.0e3 | 0.96 | 0.329 | 59.5 |
| t_sf8 | moving | SLABFLUX, SUBZONE=0 | 1697 | 36 % | 2357 (536–5095) | 0.091 | 0.35 | 9.9 | 2.03 | 7.3e3 | 1.00 | 0.330 | 52.0 |
| t_sf11 | moving | SLABFLUX + SLABTIE, SUBZONE=0 | 2346 | 0 % | 5308 (4923–5666) | 0.007 | 0.03 | 22.1 | 2.43 | 3.0e3 | 1.25 | 0.343 | 49.3 |
| t_sf7 | fixed | stack (SUBZONE=1,2) | 55 | 66 % | 647 (429–5951) | 0.26 | 1.04 | 1.5 | 1.40 | 5.8e3 | 1.65 | 0.361 | 62.0 |
| t_sf9 | fixed | SLABFLUX, SUBZONE=0 | 584 | 30 % | 1731 (505–6330) | 0.16 | 0.61 | 7.6 | 2.11 | 5.9e3 | 1.43 | 0.351 | 53.4 |
| t_sf12 | fixed | SLABFLUX + SLABTIE, SUBZONE=0 | 2525 | 0 % | 5310 (4827–5929) | 0.010 | 0.04 | 20.5 | 2.39 | 4.5e3 | 1.52 | 0.355 | 50.1 |

(rms ln about the 11-snapshot running median; failures = TR + BDF Newton failures from perf.csv, first 100 steps excluded; the first
snapshot after each restart is the checkpoint's dark state, L(100 R☉) 66–68 L☉, also in t_sf11/12. The live run z100mm, 8 threads,
over the same interval: L_ph median 0.32 L☉, 4.0 failures/100, as t_sf10.)

* **Photosphere.** The baseline spends most of the time in the dark state on both meshes (moving: L_ph 0.3 L☉ with T_rad(ph) 8400 K;
  the fixed mesh flips between the states, rms ln L_ph 1.0). With SLABFLUX + SLABTIE the bright state is reached within 10 steps of the
  restart and kept for the whole run: L_ph 1.7e3–3.9e3 L☉, T_eff 5310 K with 0.7 % (moving) / 1.0 % (fixed) scatter, photosphere in
  zone 21–22 at 55–62 R☉. For comparison the accretion luminosity G M★ Ṁ/R★ ≈ 7.5e3 L☉ (M★ 0.33 M☉, Ṁ ≈ 0.04 M☉/yr from the
  baseline's M★ growth, R★ 55 R☉): the bright state radiates ~30 % of it, the dark state ~0.3 %.
* **Physics over the same time.** M★ grows at the same rate: moving baseline 0.291 → 0.329 M☉ by t = 237829.03 yr; t_sf11 0.337 M☉
  at 237829.21 (+0.17 yr at Ṁ ≈ 0.04 → 0.336 expected); fixed baseline 0.361 at 237829.73, t_sf12 0.355 at 237829.60. R★ (perf.csv)
  is 49–50 R☉ with SLABTIE against 59–62 R☉ in the dark baseline.
* **Step cost.** The bright photosphere costs Newton convergence: 20–22 failures per 100 steps; the failed solves are slow, not
  divergent (largest corrections in E_r near the photosphere and in the front zone's species; the E_r correction sits at zones 18–26
  in 38 of 117 failures in t_sf11 and 106 of 128 in t_sf12, failz.py); Δt is limited by the failures (stol/smax 300–3000), 1.7× the
  wall time per step. In t_sf11 Δt reached 2.6e6 s right after the restart (1.1 yr in the first 100 steps) and then settled at ~3e3 s.
* **SUBZONE on top** (SLABFLUX + SUBZONE=1,2) was not run: at a front SUBZONE's D = τ_f + b (≈ 2, against Δτ_SF ≈ 0.1) replaces the
  slab Δτ (ω ≈ 1), i.e. it undoes SLABFLUX there; at a tie interface it enters only the (1 − w) part of the row.

### 4.5 Open

* The Newton convergence of the bright photosphere: candidates are the per-iterate VEF f, g, h without derivatives
  (docs/POPIII_SUBZONE_FRONT.md §6: linear contraction) and the max(·, 0) of the tie.
* n = ∞ is the smallest optical depth of the family; the grey-LTE value n = q/4 ≈ 3 (q from the opacity table at the front) would
  give Δτ ≈ 150 at the zdiagB front; not tried in a run.
* The formal solution (eddfac!) still uses piecewise-constant dtau per zone; SLABTIE corrects only the tie's donor factor. A formal
  solution on the same sub-zone profile would give h at the interface consistently.
* Multigroup: `sf_mean`, `st_tauint`, `st_hdonor` take the group's opacities and Eddington factors as arguments; per group g the same
  calls with χ_g and the group's h, f (the source of each group stays in its zone rows).
