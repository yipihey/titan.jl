# Interface radiation flux from the steady slab solution between zone centres (`SLABFLUX`)

Branch `slabflux` (from `popiii-s6` ac7c73c), 16 Nov 2026. Code: julia/TitanResidFast.jl (`SLABFLUX`, `sf_E`, `sf_mean`,
`sf_avchi!`), julia/popiii_m2c.jl (switch, log line), unit test julia/slabflux_test.jl. Test runs `D/t_sf<N>`,
D = /zpool/nvme/data/tabel_scratch/titan_popiii/runs, from `git archive` trees in /zpool/nvme/data/tabel_scratch/titan_popiii/sf_trees,
4 threads each, `numactl --physcpubind=32-35|36-39 --membind=4`.

## Summary

RESULTS_SUMMARY

## 1. The problem

The photosphere of the accreting protostar is a transparency front thinner than a photon mean free path: across one interface χ
drops by 1e3–1e5 (H⁻/H bound-free, χ ∝ T^q with q ≈ 10–20) and the optical depth of the hot half-zone is 1e3–1e5 while the cold
neighbour's is ≲ 1. TITAN's radiation-momentum row (TitanResidFast, `RADIATION MOMENTUM`) at interface k is, in the static part,

    r_k^μ [q f E]_{k−1}^{k} + ΔM_k χ̄_k F_k/c = 0,      ΔM_k = M₋ + M₊,  M± = ½ ρ± dV±/r_k^μ   (column of the two half-zones)

with χ̄_k = avchi_k = 2/(1/χ₋ + 1/χ₊), the harmonic mean (the same avchi enters the radiation force of the energy rows). Across the
front χ̄ → 2χ₊: the opaque half-zone is given the transparent zone's opacity, and the flux leaving the star follows the steep,
irradiated, lagged χ of the cold zone one to one (d ln Δτ/d ln χ₊ = 1). SUBZONE (docs/POPIII_SUBZONE_FRONT.md) caps the opaque
half-zone's depth at τ_f; in the z100mm run its Δτ is 50–300× TITAN's at the front interface (log, 16 Nov).

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

## 3. Implementation

* `SLABFLUX=1[,n]` (popiii_m2c; default off; n omitted or 0: n = ∞; n ≥ 1). After TITAN's harmonic avchi is formed, `sf_avchi!`
  replaces avchi_k of the interior interfaces ngrs+1..ngre by χ_n(χ_{k−1}, χ_k) (θ-level χ). The boundary interfaces keep TITAN's.
* With `SUBZONE` also on, SUBZONE's front form is applied on top (its Δτ_T is then the slab Δτ; at a front ω ≈ 1 and SUBZONE's
  D = τ_f a/(a + τ_f) + b wins). With `SUBZONE=0` (the stack default is `SUBZONE=1,2`) SLABFLUX alone.
* Log (every PRINTEVERY steps): `SLABFLUX: n interface(s) with dtau/dtau_T outside [2/3, 3/2] | k: ratio …` (k from the inner
  boundary).

## 4. Tests

RESULTS_TESTS
