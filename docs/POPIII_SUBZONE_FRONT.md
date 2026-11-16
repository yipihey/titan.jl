# Sub-zone treatment of an unresolved photosphere / recombination front (`SUBZONE`)

Branch `subzone-front` (from `popiii-s6` 0781136), 16 Nov 2026. Test runs: `D/t_sz<N>`, D = /zpool/nvme/data/tabel_scratch/titan_popiii/runs,
each from a `git archive` tree under /zpool/nvme/data/tabel_scratch/titan_popiii/sz_trees, restarted from zdiagB's checkpoint
`z19zdiagBS5` (step 3000) with zdiagB's own `phaseB.env` + `run.env`, `numactl --physcpubind=56-59|60-63 --membind=7`, 4 threads
(launcher: scratchpad `sz/szrun.sh`, one segment as `iceulrun.sh` runs it). Zone numbers below are TITAN indices k (snapshot
index = k − 4; the brief's "zones 62→63" are k = 67 → 68).

__RESULTS__

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
