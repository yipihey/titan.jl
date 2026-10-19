# Population III protostellar collapse with titan.jl: notes

Project: extend titan.jl so that one 1D spherical radiation-hydrodynamics calculation follows a metal-free cloud from
collapse through core formation and accretion towards the ZAMS. First milestone: reproduce the 1D collapse of
Ripamonti, Haardt, Ferrara & Colpi (2002; hereafter RHFC02) as summarised in Ripamonti & Abel (2005, §1.3, §1.5,
Fig. 1.19; hereafter RA05), compared with Omukai & Nishi (1998; ON98). Branch `popiii` (from `titan-port`).

Evidence tags (as in `PORTING_NOTES.md`, plus the ones this project needs):
- **[run]** verified by running code here; **[read]** verified by reading code in this repository.
- **[read-src]** taken from the primary source, with table/equation number. *None yet; see M0.0.*
- **[snippet]** from a web-search result summary, not from the source itself. A lead, not a fact; never used in code.
- **[memory]** my recollection, unverified. Never used in code or in the book; listed only so Tom can correct it.
- **[prompt]** stated in the project brief; to be checked against the sources.
- **[design]** a proposal, not yet implemented.

---

## M0.0 Environment (19 Oct 2026) [run]

- **Julia.** `julialang-s3.julialang.org` and `pkg.julialang.org` are refused by the session's network policy (HTTP 403
  at the proxy). Julia 1.10.4 was installed from conda-forge (`micromamba create -p /opt/jl -c conda-forge "julia>=1.10,<1.12"`,
  linked to `/usr/local/bin/julia`); this is the version pinned in `julia/env/Manifest.toml`. The General registry was
  cloned from GitHub and the packages fetched by git (`JULIA_PKG_SERVER=""`). `Pkg.instantiate()` and precompilation
  succeed. To redo after a container restart: see the commands in this paragraph.
- **Literature: not reachable.** `arxiv.org`, `export.arxiv.org`, `ar5iv`, `api.crossref.org`, `doi.org`,
  `ui.adsabs.harvard.edu`, `academic.oup.com`, `iopscience.iop.org`, `hitran.org`, `zenodo.org`, `web.archive.org`,
  Semantic Scholar and the institutional mirrors returned by search (`ricerca.sns.it`, `boa.unimib.it`) all fail
  (curl and WebFetch). Only web-search *summaries* come through. **Consequence:** no specification below is taken from
  a primary source yet, and no `.bib` entry can be verified. The extraction (M0.3, M0.4) is a checklist of what to
  extract, with [snippet]/[memory] leads marked as such. Hosts needed: `arxiv.org` (+ `export.arxiv.org`),
  `api.crossref.org`, `ui.adsabs.harvard.edu`; for the CIA data `hitran.org`; optionally the publishers.
- Machine: 4 cores reported by `nproc` (the brief says 2; runs will be sized for 2).

---

## M0.1 What TITAN does at r = 0 and with self-gravity [read]

The verified paths never had a zone edge at r = 0: SUBC/SUPC and the quasi-star runs all have r(ngrs) > 0. Every place
below breaks or is untested with r(ngrs) = 0, or with gravity from the gas alone.

| # | where | what it does | at r(ngrs) = 0 |
|---|---|---|---|
| 1 | `TitanResidFast.residual!` Newton variables (l. 336–345) | r, x_me, ρ, T, E are corrected multiplicatively, `base·exp(x)` | r = 0 stays 0, but ∂R/∂x = 0 for the grid row of interface ngrs: **singular row**. Needs an additive correction for that interface or a fixed-value row `R = x`. |
| 2 | viscosity `dudr`, energy rows `aur` (l. 394, 596, 652) | u_k/r_k terms | **0/0 → NaN** at k = ngrs (u = 0, r = 0). Needs the regular limit (u/r → ∂u/∂r at the centre, or drop the term: u(ngrs) = 0 makes it vanish analytically). |
| 3 | `TitanCore.timstp!` (l. 108ff) | relative change (rn − ro)/ro, (x_me − x_meo)/x_meo | **0/0 → NaN** for r at ngrs. Exclude the centre. |
| 4 | `TitanCore.eddfac!` spherical (l. 536ff) | `ncor` core rays with impact parameters p_j = (j−1) r(ngrs)/ncor and angle √(1 − p²/r(ngrs)²) | **0/0** in the core-ray angles. Needs a centre-regular variant: no core, rays tangent to each shell pass through the centre (symmetric). |
| 5 | `TitanCore` totnrg (l. 629, 639) | gravitational energy terms x_m/r | 0/0 at the centre (diagnostic only). |
| 6 | gravity (`gasmom`, l. 542; MLT l. 564) | g = 4πG x_m / r² with x_m = x_mtot − x_me (per steradian, θ-centred) | Used only at interior interfaces ngrs+1…ngre, so r = 0 itself is not evaluated. **Never exercised with self-gravity**: every benchmark has cgrav = 0; the quasi-star runs had a dominant point mass. |
| 7 | enclosed mass x_m = x_mtot − x_me | v1.2's exterior-mass variable, chosen for steep *outer* envelopes | **Precision at the centre.** Newton converges on δ ln x_me ≤ conv = 1e-5, i.e. x_m near the centre is known to ~1e-5 x_mtot. For a ~10³ M☉ cloud that is ~1e-2 M☉, more than the initial core (~3e-3 M☉ [prompt]). Round-off alone loses log₁₀(x_mtot/x_m) digits. **A regular-centre mode needs the interior mass as the variable** (or a Lagrangian grid, where x_me is fixed per interface). |
| 8 | inner hydro boundary | `leibc = 1` wall: R_r = r_n − r_o, u(ngrs) fixed (row `R = x`), ghosts by reflection; `llibc = 1` piston | Zero velocity and zero mass flux at r = 0 are available (phil0 = 0). Mass row at ngrs: x_me(ngrs) − x_meo(ngrs) = 0 → x_m(ngrs) = 0 if set initially. |
| 9 | inner radiation boundary | `lribc = 2` reflecting: F(ngrs) = 0, ghost F mirrored (l. 691) | Zero flux at r = 0 is available. The residual code refuses lribc = 1. |
| 10 | adaptive grid (`lgrid = 1`) | Dorfi–Drury equation; the rows next to each boundary set ν(ngrs+1) = ν(ngrs) | `ladx = 2` uses ν = (r_{k+1}+r_k)/(2Δr): at the centre zone ν = ½ for any Δr, so ν carries no information about the central zone width. `ladx = 1` (ν = x_scale/Δr) is regular. Resolution weights in ln ρ over 15–20 decades: untested. P7.6–P7.7 show the grid works on a shock but needs τ_grid tuning. |
| 11 | `peq!` extrapolation (l. 360–385) | r, x_me extrapolated with old dr/dt, dm/dt; monotonicity check with retries | fine at r = 0 (dr/dt = 0). |
| 12 | array size | `TITAN_MGR` (P7.7), default 330 | fine. |

**Conclusion [read].** A regular centre is new code. Minimal set: (1) additive/fixed centre radius row; (2) regular
limits of u/r; (3) timstp exclusion; (4) centre-regular tangent-ray Eddington factors (or f = 1/3 inside the
optically thick core, as a first step); (5) diagnostics; (7) interior-mass variable or a Lagrangian grid. All behind
a switch (`REGULAR_CENTRE[]`), so the verified paths stay unchanged; SOD/SUBC/SUPC and the quasi-star hydrostatic test
are rerun after the change.

**Grid strategy for the collapse [design].** A Lagrangian grid (`lgrid = 3`) in mass, with shell masses decreasing
geometrically towards the centre, makes x_me exact per interface (item 7 disappears) and needs no species advection.
[memory] RHFC02 and ON98 were Lagrangian codes. The adaptive grid becomes necessary at the accretion shock (M3–M4); by
then the interior-mass variable must exist. Plan: M1 Lagrangian first, then the adaptive grid with the interior-mass
variable as a second, cross-checked configuration.

## M0.2 What the material interface must provide for a multi-species, non-equilibrium gas [read + design]

Today (`TitanMaterials.jl`, `TitanResidFast` l. 234–265) a material is a function of (ρ, T) only:
`material_gas → (p_g, e_g, a_s)`, `material_gasne → (…, n_e)`, `material_opac → (χ_F, κ_E, κ_P)`,
`material_thermo → (∇_ad, c_P, δ)`. All are called with ForwardDiff Duals inside the residual, and with Float64 in
the state fills. There is no composition slot, and no source term other than gas–radiation exchange.

Needed for primordial chemistry:

1. **Composition per zone.** x = (H, H⁺, H⁻, H₂, H₂⁺, e⁻, He [, He⁺, He²⁺, D, D⁺, HD]) per unit mass. Interface:
   `material_gas(m, ρ, T, k)` (zone index, so the material reads its composition array), or a `Composition` argument.
2. **EOS with composition.** p = n k T with n(x); e(ρ, T, x) = translational + H₂ rotation/vibration (ortho/para
   choice explicit) + chemical binding energies (H₂ dissociation 4.48 eV, H ionisation 13.6 eV, He). Γ₁, c_P, ∇_ad
   from the same e(T). T is TITAN's variable, so no e → T inversion inside the residual.
3. **Opacity with composition.** (χ_F, κ_E, κ_P)(ρ, T, x) from one κ_ν model (CIA + H⁻ + lines + Rayleigh + e⁻
   scattering + Z0 atomic at high T); see M2b.
4. **A cooling/heating source term** in the gas-energy row that does not feed E_r: optically thin or escape-corrected
   H₂ line cooling (the line photons are not in TITAN's grey E_r), chemical heating/cooling. Where a process's photons
   *are* in E_r (CIE continuum, H⁻), it must enter through κ_P/κ_E and nowhere else (double counting, M2b.4). Natural
   home: the `GASROW` gas-only energy row (P7.5), which already has an exchange term.
5. **Species evolution.** Two options:
   - **(A) operator split [design, recommended to start].** Species stored per zone, attached to mass shells
     (Lagrangian: no advection; adaptive: advect mass fractions with TITAN's mass fluxes and limiter). After each
     hydro step, integrate the stiff network implicitly over Δt at the (ρ, e) path of the step, then re-solve T at
     fixed e (Strang if needed). Cheap; testable on its own; but the 3-body heating couples chemistry and T on
     timescales ≪ Δt during the chemo-thermal instability, where splitting may bias the track.
   - **(B) fully coupled.** Species as extra Newton unknowns (neqn = 7 + n_s). Exact coupling and conservation by
     construction; the block solver is static 7×7 (`BlockSolver.Penta`), the AD sweep cost grows with ncol·neqn.
   - Test that decides: the one-zone and centre tracks of M2 with (A) at decreasing Δt against (B) (or against a
     tightly coupled one-zone ODE solution). Document in `decisions.md`.
6. **Conservation to round-off.** H nuclei and charge (and D, He) enforced by construction (eliminate one species per
   conservation law, or project), checked every step.

## M0.3 Specification of the RHFC02 calculation (to extract) — status: NOT EXTRACTED

Source unreachable (M0.0). The table lists what has to come out of the paper; the right-hand column is a lead only.

| item | needed | lead (not usable as a fact) |
|---|---|---|
| bibliographic | authors, title, journal, vol, pages, DOI, arXiv | [snippet] "Radiation from the first forming stars", MNRAS 334, 401 (2002), arXiv:astro-ph/0107095; DOI 10.1046/j.1365-8711.2002.05516.x [prompt] |
| method | Lagrangian/Eulerian, zone count, time integration, artificial viscosity | [snippet] "spherically symmetric hydrodynamical scheme including radiative transfer and primordial chemistry" |
| initial conditions | profile (ρ, T, u, x_H₂, x_e), total mass, radius, central density; whether several ICs (incl. from 3D cosmological runs) | [snippet] several initial configurations, including ones taken from 3D cosmological simulations |
| chemistry | species, reactions, rate sources; three-body rate; collisional dissociation; heating per reaction | — |
| H₂ thin cooling | which cooling function / level populations | — |
| H₂ line transfer | method (frequency-dependent? escape probability?), number of lines, line profile, velocity field | — |
| CIE | included? which fit or data | — |
| continuum opacities | sources (H⁻, CIA, Rayleigh, e⁻), mean or frequency-dependent | — |
| EOS | γ(T), H₂ internal degrees of freedom, ortho/para, dissociation energy | — |
| boundaries | inner (r = 0), outer (pressure, fixed radius, free) | — |
| stopping criterion, validity | why the run ends; reliability limits | [prompt] RA05: reliable only for M_core ≲ 0.1 M☉ |
| outputs | 10 stages of Fig. 1.19 (ρ, u, T, x_H₂ vs enclosed mass); times/central densities of the stages | — |
| core | n at core formation, initial core mass, accretion rate | [prompt] n ≳ 1e21 cm⁻³, ~3e-3 M☉, 0.01–0.1 M☉/yr; [snippet] core grows from ~3e-3 to ~4.5e-2 M☉ "in about 0.30 years" (text garbled or exact; check) |
| luminosity | L(t) | [snippet] L rises slowly to ~1e36 erg/s during the self-similar phase, mostly H₂ IR lines; collapse of the centre takes 1e5–1e6 yr |

## M0.4 Specification of the Ripamonti & Abel (2004) CIE treatment (to extract) — status: NOT EXTRACTED

| item | needed | lead |
|---|---|---|
| bibliographic | | [snippet] arXiv:astro-ph/0311355 (Nov 2003); MNRAS 348, 1019, DOI 10.1111/j.1365-2966.2004.07422.x [prompt] |
| H₂ line cooling thin→thick approximation | functional form, the optical-depth variable, validity range | [snippet] "new analytic approximations for H₂ cooling bridging the optically thin and thick regimes" |
| CIE cooling | fit, data source, range | [prompt] review Eq. 1.29: L_CIE ≈ 0.072 ρ T⁴ X f_H₂ erg g⁻¹ s⁻¹, 400–7000 K |
| CIE optical depth | escape factor and its τ estimate | — |
| instabilities | linear growth rates (3-body at 1e8–1e11 cm⁻³; CIE at ≳1e13–1e14) | [snippet] neither leads to fragmentation; the CIE one is "as weak as" the low-density one |
| 1D test of simplified transfer | set-up, comparison quantities, outcome | [snippet, citing review] the 1D model becomes optically thick to CIE at n_H ≈ 1e16 cm⁻³ |

## M0.5 ChemistryKernels.jl (yipihey, local clone) [read]

Built for cosmological cells (reduced Abel–Anninos–Zhang–Norman 1997 network, operator-split sub-cycling). What it
offers here, and what it does not:

- **Reusable:** rate-coefficient functions k1–k22 (+ D, HeH⁺), scalar, generic in the number type (`R = typeof(T)`),
  pure: likely AD-compatible (to test with Duals in M2). Chemical heating `_h2_chem_heat` (`src/fast.jl:91`) with
  per-channel energies (3.53 / 1.83 / 4.48 eV) attributed to Palla–Salpeter–Stahler 1983 / Omukai 2000 (to verify).
  H₂ cooling (Glover & Abel 2008 thin + LTE, two-level interpolation) is in EmissionKernels.jl (`cooling_h2.jl`).
- **Species storage:** "mass-equivalent ×N" (`yH2I = 2 n(H₂)`, `yHeX = 4 n(HeX)`, …). H⁻, H₂⁺ algebraic equilibrium;
  HII, H₂, HD advected; He in collisional–radiative equilibrium; n_e from charge.
- **Three-body H₂ (k22):** one fixed fit, 1.3e-32 (T/300)^−0.38 (T ≤ 300 K) and (T/300)^−1 above (`src/rates_h2.jl:209`),
  comment "three_body_rate = 0". Source not stated in the code. Must become a selectable parameter with each group's
  rate copied from its source (Glover & Abel 2008 compare them).
- **Collisional dissociation (k13):** Abel et al. 1997 fit, **zero below T_eV = 0.3 (≈ 3480 K)**, low-density limit
  only (no n-dependence towards LTE). The "three-body equilibrium" shortcut (`_threebody_h2_equilibrium`,
  `network_step.jl:46`) balances k22 against this k13, which are not detailed-balance partners: below 3480 K it
  predicts full conversion. For n ≳ 1e14 cm⁻³ and T ≳ 2000 K (the CIE regime) the H₂ ⇌ 2H balance must come from
  detailed balance (partition functions) or a consistent pair of rates. **Not usable as is for this problem.**
- **Chemical heating suppression:** `fcr = n_H/(n_H + 1e8)` (hard-wired critical density) on formation heating. To
  check against the sources.
- **Deuterium:** commit `faaf324` clips a negative k50 fit to zero (non-negativity). D nuclei are updated as three
  independent equations, with no projection onto total-D conservation (H nuclei are projected, incl. the H in HD).
  So total D is not conserved to round-off by construction [read; to measure in M2].
- **Solver:** clamps and branch switches (non-smooth, not suited to Newton/AD); sub-cycling with a 10 % limiter.
  For TITAN a new implicit species solve is needed; the rates are what to reuse.

**Recommendation [design].** Reuse the CK rate functions where the sources match RHFC02/modern compilations (each
checked against its source table); write a new, small, implicit, conservation-by-construction network inside a new
material `PrimordialGas`; keep CK as a cross-check for one-zone tracks.

## M0.6 Gap list: titan.jl vs the RHFC02 calculation

| physics | titan.jl now | needed | milestone |
|---|---|---|---|
| self-gravity from gas | code path exists, never tested | test: Bonnor–Ebert static, free fall, Larson–Penston | M1 |
| regular centre | absent (M0.1 items 1–5, 7) | switch + tests | M1 |
| dynamic range | adaptive grid tested over a shock; mgr ≤ ~500 | Lagrangian mass grid over 15–20 decades in ρ_c; later adaptive with interior mass | M1, M3 |
| species, chemistry | absent | non-equilibrium H, H⁺, H⁻, H₂, H₂⁺, e⁻, He (+D?) | M2 |
| EOS | ideal (γ, μ) or Z0 Saha (atomic) | variable μ, H₂ rot/vib, binding energies, non-ideal at n ≳ 1e18–1e21 | M2–M3 |
| H₂ thin cooling | absent | cooling function (source per RHFC02 and modern) | M2 |
| H₂ thick lines | absent (NLTEH has Sobolev for atomic H) | escape probability → RHFC02's line transfer | M3 |
| CIA κ_ν, CIE | absent | κ_ν from tabulated pair coefficients (H₂–H₂, H₂–He, H₂–H), Kirchhoff emission | M2b |
| low-T continuum opacities | Z0 table, log T ≥ 3.0, no molecules (κ a lower bound below 2000 K) | H⁻ bf/ff, Rayleigh, CIA; means on (ρ, T, x) | M2b |
| frequency-dependent transfer | grey two-moment in the hydro; `Spectrum.jl` post-processing | grey/multigroup with the computed means; multigroup optional | M2b, M3 |
| energy accounting | total-energy row exists (P4.2), gas row (P7.5) | chemical energy in the ledger | M2 |
| nuclear burning | absent | D, pp, CN | M5 |

## M0.7 Milestone plan and test matrix

Each milestone ends with a passing test, a section here, a commit on `popiii`. Logs in `lib/popiii/` (committed).

**M1 Self-gravity and regular centre (simple ideal gas, γ = 5/3 and isothermal-like γ → 1).**

| test | pass criterion | reference |
|---|---|---|
| regular centre, static uniform sphere without gravity | stays static to round-off; no NaN; Jacobian = FD | analytic |
| Bonnor–Ebert sphere (ξ_max below critical), isothermal-like EOS | ρ(r) unchanged to < 1e-3 over 10 sound-crossing times; u/c_s < 1e-3 | Lane–Emden isothermal solution, integrated here [computed] |
| pressure-free collapse (uniform sphere, p ≈ 0) | central collapse time = t_ff = √(3π/32Gρ₀) to < 1 %; r(m, t) follows the parametric solution | analytic (derivable) |
| Larson–Penston | approach to the self-similar profile (ρ r² → const, u → −3.3 c_s [memory, verify]) over ≥ 8 decades of ρ_c | self-similar ODE solved here [computed]; constants checked against Larson 1969, Penston 1969, Hunter 1977 when reachable |
| energy | total energy (incl. gravitational) conserved to the scheme's level; report the number | energy ledger |
| regression | SOD, SUBC, SUPC, quasi-star hydrostatic: unchanged | `metrics_reference_fixed.csv`, P7.3 |
| resolution | zone count 2× / 4× on LP | convergence |

**M2 Chemistry and thin cooling.** One-zone free-fall track (ρ ∝ free fall) with the new network: T(n) near ~200 K at
~1e4 cm⁻³, three-body conversion near 1e8–1e10, full molecularity; numbers checked against the sources (RA05 Fig.,
ON98, Glover & Abel 2008). Conservation of H, He, charge to round-off. Split (A) vs coupled (B) at decreasing Δt.
Three-body rate as parameter, bracketed. Then the same at the centre of the 1D code.

**M2b CIA opacity and the continuum transition.** κ_ν^CIA from tables (units checked per table); emission by Kirchhoff;
thin CIE rate = 4π∫κ_ν B_ν dν vs RA04 and RA05 Eq. 1.29 (table + figure); Planck/Rosseland means on (ρ, T) for the
non-equilibrium composition and, for equilibrium composition, vs Lenzuni et al. 1991 and Mayer & Duschl 2005; RA04
thin/thick comparison along their track; reproduction of RA04's 1D simplified-transfer test.

**M3 Thick H₂ lines, CIE, core.** Run to core formation; the 10 stages of RHFC02 / RA05 Fig. 1.19; ON98 comparison
(T–n track, core mass, onset Ṁ); three-body bracket; resolution and time-step studies.

**M4 Accretion phase**, **M5 nuclear burning**: planned in writing when M3 is done.

## M0.8 Questions for Tom

1. **Network access** for the literature (arXiv, Crossref, ADS; HITRAN for CIA) — or PDFs uploaded to the session.
   Without it M0.3/M0.4 cannot be done, and no rate can be entered.
2. Helium and deuterium chemistry in the first version? (Proposal: He as neutral spectator + CIA partner first;
   D/HD later — HD cooling matters little in this regime [memory].)
3. Faithfulness: copy RHFC02's line transfer and CIE fits (to reproduce), or modern line and CIA data (to be right)?
   (Proposal: both, as switches; reproduce first, then modernise and show the difference.)
4. The third book: `popiii/` at the repository top level, published under `/popiii/` like `spectra/`?
5. Grid: Lagrangian mass grid for the collapse (M1–M3), adaptive with interior mass for accretion (M4)?

## M0.9 Baseline of the verified paths in this environment [run]

Before any change, on branch `popiii` = `titan-port` @ 67d89d3, Julia 1.10.4 (conda-forge build), 1 thread:

```
cd julia; for P in SUPC SUBC SOD; do Z0_TABLE=ext TITAN_NOLIB=1 julia -t 1 --project=env phase6_run.jl pure $P; done
julia --project=env popiii_regress.jl          # new: tier-(c) observables vs reference_runs/metrics_reference_fixed.csv
```
Logs: `lib/popiii/m0_baseline_{SOD,SUBC,SUPC}.log`, `lib/popiii/m0_baseline_regress.log`.

| problem | steps here | P6 record | fixed Fortran reference | tier (c) at tierc.jl's times | ms/step |
|---|---|---|---|---|---|
| SUPC | 713 | 713 | 713 | ACCEPT | 5.2 |
| SUBC | 587 | 587 | 587 | ACCEPT | 6.2 |
| SOD | **464** | 455 | 452 | ACCEPT (+2.7 % steps, inside the 3 % criterion) | 2.9 |

- SOD's step count differs from P6's 455. P3.3/P6 record that SOD's count moves with last-bit differences; this is
  a different Julia build and CPU. Not investigated further. [run]
- `popiii_regress.jl` also prints the reference's SOD sample at t = 0.5 (after the wall reflection), which `tierc.jl`
  never used for acceptance: positions agree to ≤ 1.1e-3 L, but p* differs by 1.1 % and u* by 4.4 % from the
  reference. Kept as information; the P6 run's values at t = 0.5 were not stored, so it is unknown whether this is new.
- The quasi-star hydrostatic test (P7.3) was not rerun: M0 changes no code. It is part of the M1 regression set.
