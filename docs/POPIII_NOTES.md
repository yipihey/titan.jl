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

## M0.3 Specification of the RHFC02 calculation — status: EXTRACTED (19 Oct 2026) [read-src]

Source: Ripamonti E., Haardt F., Ferrara A., Colpi M., 2002, MNRAS 334, 401–418, "Radiation from the first forming
stars", DOI 10.1046/j.1365-8711.2002.05516.x (Crossref, checked), arXiv:astro-ph/0107095v2. Page numbers below are
those of the arXiv v2 preprint (39 pp.); equation, table and figure numbers are the paper's. Local copy (not in the
repository): `~/Projects/popiii_lit/RHFC02/`.

| item | RHFC02 says | where |
|---|---|---|
| hydro scheme | 1D **Lagrangian**, explicit finite differences; Eqs. (A1)–(A3) in M_r with radiative force f_rad = (1/c)∫κ_ν F_ν dν; discretised as Thoul & Weinberg 1995 (TW95) Eqs. 9–13, 15; "second order in space and time" | §2.1 Eqs. 2–3; App. A1.1, A2.1 |
| artificial viscosity | q ∝ ρ Δv² where shells would interpenetrate (Bowers & Wilson 1991; TW95), c_q = 4. Shocks spread over several shells; "inaccurate when dealing with strong shocks" | Eq. A4; A2.1; §3.2.2 |
| energy equation | from TW95 Eq. 14; the chemical term makes u(T) implicit, solved by **bisection**, so **first order in time** | A2.1 |
| time step | min of: 2Δt_prev; C_dyn (G r³/M_r)^½; C_Cour Δr/c_s; C_cool uρ/(Γ−Λ); C_cross Δr/Δv; with C_dyn = 1, C_Cour = 0.2, C_cool = 0.1, C_cross = 0.05; a posteriori ΔT/T ≤ 0.5 % per step | A2.2 |
| boundaries | centre: r = 0, M_0 = 0, v_in = 0; outside: constant pressure beyond the outermost shell | A2.3 |
| zoning | **100 Lagrangian shells** (P100/2: 140), massive outside, light inside: in P100 outer shells ≳ 1 M☉, central shell **3e-7 M☉**. P100/2's extra 40 shells refine −3 ≲ log M_r ≲ −1 (where the shock forms). The mass distribution law is not given | §3.1 p. 7 |
| initial conditions | at rest; constant f_H2; hydrostatic with polytropic index 5/3 (P) or isothermal (I); N1500 from ABN00. f_e = 1e-10. He: n_He/n_H = 1/12 (X = 0.75, Y = 0.25) | §3.1, Table 1; A1.2 |
| runs (Table 1) | ρ_0,c in m_H cm⁻³, T_0,c, f_H2: P100 2.6e6, 270 K, 5e-4 (cf. ON98 A,B); P100/2 same, 140 shells; P1000 2.3e4, 246, 5e-4; P10 3.5e8, 249, 5e-3 (cf. ON98 C); I100 7.8e7, 300, 5e-4; I1000 1.1e7, 300, 5e-4; N1500 2.0e6, 440, 5e-4. Cloud radius is **not stated** (follows from hydrostatic + M_tot + outer pressure, also not stated) | Table 1 p. 7 |
| chemistry | species H⁰, H⁺, H⁻, H₂, H₂⁺, He⁰, He⁺, He⁺⁺, e⁻ (9); **no D**. Reactions (Table A1): H⁺+e→H+γ (AANZ98 k2, **case A**), H+e→H⁻+γ (k7), H+H⁻→H₂+e (k8), H+e→H⁺+2e (k1) from "AANZ98"; **3-body H+H+H→H₂+H (PSS83 k4)**, H₂+H→3H (k5), H+H+H₂→2H₂ (k6), H₂+H₂→2H+H₂ (k7), H+H→H⁺+H+e (k9) from Palla, Salpeter & Stahler 1983. No H₂⁺ channel listed although H₂⁺ is a species | Table A1 p. 29; A1.2 |
| chemistry switch | for ρ ≳ ρ_eq = **5e-11 g cm⁻³** switch to **equilibrium** chemistry: Saha for p ≤ 1e4 dyn cm⁻², SCVH95 tables for p ≥ 1e6, interpolated in between | A1.2 pp. 29–30 |
| EOS | ideal (p ≤ 1e4): p = nkT; u = u_th + u_chem, u_th = N_tot kT/(γ_ad − 1), u_chem = −N_tot Σ y_i χ_i; γ_ad a weighted mean: 5/3 monatomic, 7/5 for H₂⁺, and for H₂ **1/(γ−1) = 5/2 + x² eˣ/(eˣ−1)², x = 6100 K/T** (from ON98: rotation fully excited + vibration). Non-ideal (p ≥ 1e6): p(ρ,T), u(ρ,T) from **Saumon, Chabrier & Van Horn 1995** (SCVH95), interpolated between | A1.3 Eqs. A5–A9 |
| H₂ lines | levels **v = 0–2, J = 0–20**, all quadrupole-allowed transitions; A-values **Turner, Kirby-Docken & Dalgarno 1977**; **LTE populations**; α_ν from Eq. A11; **Gaussian profile, 30 frequency bins**, thermal broadening and bulk-velocity **Doppler shifts**; optically thin approximation where line absorption depth < 0.1. No separate thin cooling function | A1.4; A2.4 Eq. A11 |
| transfer method | tangent rays (Hummer & Rybicki 1971; ON98), 1 + 2i angles at the outer radius r_i of shell i; cooling Λ_rad,i = [L(r_i) − L(r_{i−1})]/[(4π/3)(r_i³ − r_{i−1}³)] | A2.4 Eq. A10 |
| when transfer starts | radiative transfer from stage 1 (Table 2 "Start Rad. Trans."); continuum from **T > 1000 K** | Table 2; A1.4; §3.3.1 |
| continuum | **grey**, LTE source function; **Planck mean of Lenzuni, Chernoff & Salpeter 1991** (Z = 0) for 1000 < T < 7000 K; **Rosseland mean of Rogers & Iglesias 1992** for T > 7000 K; **diffusion approximation** (Eq. A12) when τ > 4 in a single shell. **No separate CIE cooling term**: CIA enters only through Lenzuni's Planck mean (to check in Lenzuni 1991 that their mean includes H₂–H₂ CIA) | A1.4; A2.4 Eq. A12 |
| freezing (accretion) | after core formation, a shell is frozen (v = 0, all else fixed) if v/v_ff < 1e-3, T > 5e4 K and ≥ 5 shells inside the shell of maximum infall. Core radius shrinkage then neglected (innermost 10 shells: 1.6e10 → 1.0e10 cm between stages 8 and 9 in the full run); shock position differs by ≃ 4 % | §3.2.4 p. 13 |
| stopping | **T_c ≃ 1e5 K** (stage 9): ρ_c ≃ **0.3 g cm⁻³**, Courant time ≲ 300 s | §3.2.3 p. 12 |
| stages (Table 2, "run P100"; Fig. 1 caption says **P100/2**) | stage: T_c [K], Δt [yr] (since previous), t_f − t_i [yr] (to the end): 0: 270, —, 740000 (IC); 1: 450, 720000, 20000 (start RT); 2: 650, 19000, 930; 3: 1500, 920, 9.6; 4: 3000, 9.2, 0.40 (H₂ dissociating); 5: 6500, 0.060, 0.34; 6: 15000, 0.005, 0.34 (core forming); 7: 30000, 0.002, 0.34 (core formed); 8: 65000, 0.011, 0.33; 9: 100000, 0.33, — (final). Central densities per stage are **not tabulated** (only in Figs. 1, 2, 7) | Table 2 p. 9 |
| collapse | Larson–Penston-like, ρ ∝ r^−α, **α ≃ 2.2**; self-similarity ends between stages 6 and 7. 3-body H₂ formation from n ≳ 1e8 cm⁻³ (stage 1); H₂ formation heating dominates Γ_ad for T ≲ 1550 K; dissociation (4.48 eV/molecule) above 1550 K; line cooling dominates for T ≲ 1800 K, continuum above; centre nearly adiabatic for T ≳ 2350 K; f_H2 ≃ 5 % at stage 6 (T_c ~ 15000 K) | §3.2, §3.2.1 |
| optical depth | point A (τ = 1, Figs. 2, 7) and B (core formation); at T_c ≃ 2000 K only ≃ 0.02 of the line flux escapes locally; stage 2: 14 lines with τ > 0.1; T_c = 1000 K: 42 lines (half τ > 1) | §3.3.1 |
| core | core = mass inside the radius where infall < **0.1 v_ff**, v_ff = (2GM_r/r)^½. At formation (P100, Table 3): M = 6e30 g (≈ **3e-3 M☉**), R = 1e11 cm, T_c = 2e4 K, ρ_c = **1e-2 g cm⁻³**, f_H2 ≲ 0.03. Initial core mass ≈ M_Jeans(T_B, n_B) and ≈ the minimum of M_Jeans(T_c, n_c) (Fig. 12). Pressure ionisation: at stage 9 n_H⁺/n_H ~ 2 | §3.2.2, §3.4.1, Table 3 |
| accretion | core grows from ~0.3e-2 to ~4.5e-2 M☉ in ~0.30 yr: **Ṁ ~ 0.14 M☉/yr**; fit M_*(t) = 9.2e-2 M☉ (t/yr + 0.0063)^0.657 (Eq. 5) vs ON98 M_*(t) = 0.11 M☉ ((t − t_cf)/yr)^0.73 (Eq. 4); P100 reaches ~0.3 M☉ in ~6 yr (frozen-shell run, "indicative") | §3.2.4 Eqs. 4–5, Fig. 6; §3.4.2 |
| luminosity | lines dominate; at τ_c = 1: lines 1.06e36, continuum 7.5e35, total 1.81e36 erg/s (Table 4, run P100/2); continuum photosphere at ≃ 2300–2400 K; full line list Table B1 at T_c = 650, 1000, 1500 K and τ_c = 1 | §3.3, Table 4, Table B1 |
| radiation force | negligible: κ < 0.05 κ_es everywhere before core formation; removing f_rad changes nothing | §3.3.3 |
| validity | shock treatment approximate (artificial viscosity; transition region of H⁰ possibly numerical); results after the first 0.34 yr of the frozen run "indicative"; stated comparison: "good agreement with ON98" for P100, P100/2, P10 | §3.2, §3.2.4 |

**Reference-list entries to verify in M0.3b** (as printed in RHFC02): AANZ98 is listed as "Abel T., Anninos P.,
Norman M. L., Zhang Y. 1998, ApJ, 518" — incomplete; the rates k1, k2, k7, k8 are those of Abel, Anninos, Zhang &
Norman 1997 (New Astron. 2, 181) [memory] and/or Anninos et al. 1997; to resolve. PSS83 = ApJ 271, 632; Lenzuni,
Chernoff & Salpeter 1991, ApJS 76, 759; Rogers & Iglesias 1992, ApJS 79, 507; Turner, Kirby-Docken & Dalgarno 1977,
ApJS 35, 281; TW95 = ApJ 442, 480; ON98 = ApJ 508, 141; Omukai 2000, ApJ 534, 809; L69 = MNRAS 145, 271; Penston
1969, MNRAS 144, 425; Hunter 1977, ApJ 218, 834. **SCVH95 is cited in the text but missing from the reference list**
(Saumon, Chabrier & Van Horn 1995, ApJS 99, 713 [memory]).

**Consequences for the plan [design].**
1. **Two gaps M0.6 did not list:** (a) a non-ideal EOS + equilibrium composition table (SCVH95) for p ≥ 1e6 dyn cm⁻²
   (stage ≳ 6), with Saha below; (b) frequency-resolved H₂ line transfer (≈ 300 lines × 30 bins, Doppler-shifted,
   tangent rays) feeding a **local cooling rate** Λ_rad,i, not a grey E_r. In TITAN terms the lines are a
   *nonlocal source term* in the gas-energy row, computed by a formal solution per step (operator split, as RHFC02
   effectively do with an explicit scheme); only the continuum is TITAN's grey E_r.
2. **"Faithful" mode = RHFC02's switches**: equilibrium chemistry above 5e-11 g cm⁻³; continuum only above 1000 K;
   Lenzuni Planck / Rogers–Iglesias Rosseland; LTE line populations; PSS83 3-body rates; case-A recombination.
   **No separate CIE term** — the CIE physics of RA04 is the "modern" switch (M2b), and Lenzuni's mean is the
   faithful one (Lenzuni 1991 tables needed).
3. **Lagrangian with 100–140 shells** reproduces their resolution; the shell-mass law must be chosen by us (central
   shell 3e-7 M☉ of 100 M☉; outer shells ≳ 1 M☉) and documented. Their result has ~2 shells per decade of M_r: our
   resolution study (M1 test matrix) will show what that resolution misses.
4. **Comparison data available:** Table 2 (T_c and times per stage), Table 3 (core at formation), Eqs. 4–5 (core
   growth), Table 4 / B1 (luminosities), Figs. 1, 2, 6, 7 (to digitise for ρ, v, T, f_H2 profiles and the n_c–T_c
   track). The cloud's outer pressure and radius are not given: the P-runs are defined by (M_tot, ρ_0,c, T_0,c,
   polytrope 5/3 in hydrostatic equilibrium), which fixes the profile up to the truncation radius — to be chosen and
   stated.
5. The Table 2 / Fig. 1 run label conflict (P100 vs P100/2) stays open; Table 4 is P100/2.

## M0.4 Specification of the Ripamonti & Abel (2004) treatment — status: EXTRACTED (19 Oct 2026) [read-src]

Source: Ripamonti E., Abel T., 2004, MNRAS 348, 1019–1034, "Fragmentation and the formation of primordial
protostars: the possible role of collision-induced emission", DOI 10.1111/j.1365-2966.2004.07422.x (Crossref,
checked), arXiv:astro-ph/0311355v1. Pages are the preprint's (1–17). Full extraction with every formula verbatim:
`~/Projects/popiii_lit/RA04/extraction.md` (extracted by a subagent from text + PDF; the items below marked ✓ I
checked against the text myself; the others are to be re-read before entering code).

| item | RA04 says | where |
|---|---|---|
| thin H₂ lines | Hollenbach & McKee 1979 LTE fit (as GP98's high-density limit): L = (X f_H2/m_H) H(T₃), T₃ = T/1000 K, H = 9.5e-22 T₃^3.76/(1+0.12 T₃^2.1) e^{−(0.13/T₃)³} + 3e-24 e^{−0.51/T₃} + 6.7e-19 e^{−5.86/T₃} + 1.6e-18 e^{−11.7/T₃} erg s⁻¹ (per molecule). X = 0.75 | Eqs. 6–7 p. 4 ✓ |
| **factor 2** | Eq. 6 has X f_H2/m_H, Eq. 8 (line sum) has X f_H2/(2 m_H) for the same quantity; with f_H2 = fraction of H nuclei in H₂, molecules per gram = X f_H2/(2 m_H), so Eq. 6 is 2× too large as printed. Which one their code used: to decide against their Fig. 1 | Eqs. 6, 8 p. 4 ✓ |
| thick lines (used) | **L_thick = L_thin · min(1, (n/n₀)^−β), β = 0.45, n₀ = 8e9 cm⁻³** (ρ₀ = 1.34e-14 g cm⁻³), n = ρ/m_H; "usable to n ~ 1e16" if the collapse stays roughly spherical. (A full per-line model, Eqs. 9–21 with ξ = 0.20, is derived but the 1D test uses Eq. 22.) Eqs. 36, A1 print max[…] where min is meant | Eq. 22 p. 9 ✓ (min) |
| CIE data | H₂–H₂ Borysow 2002 (400–1000 K), Borysow et al. 2001 (1000–7000 K); H₂–He Jørgensen et al. 2000; H₂–H Gustafsson & Frommhold 2003; H–He Gustafsson & Frommhold 2001. Complete only for 1500–2500 K; extended to 400–7000 K by keeping each pair's ratio to H₂–H₂ fixed | Table 1 p. 3 |
| CIE fit | L_CIE ≈ A_CIE ρ T^α X f_H2, **A = 0.072 erg cm³ g⁻² s⁻¹ K⁻⁴, α = 4**; valid f_H2 ≈ 1 (roughly to 0.5); He and H pairs folded into A; used for 400–3000 K | Eq. 5 p. 4 ✓ (form) |
| CIE optical depth | L_CIE (1 − e^{−τ_c})/τ_c, **τ_c = (ρ/3.3e-8 g cm⁻³)^2.8**, an empirical fit to RHFC02's continuum optical depth (not from an opacity), with A lowered to 0.054 to match RHFC02; "good fit up to n ≳ 1e19". τ_c = 1 at n ≈ 2e16 cm⁻³ [computed] | §3.3 p. 14; Eq. A1 |
| chemistry (analytic model) | dF/dt = n_f[2 n_p k₄ (1−f)² − k₅ f], n_f = n_p(1 − 15f/16), n_p = ρX/m_H; k₄ = 5.5e-29/T cm⁶ s⁻¹ (PSS83); k₅ = 2.2e-9 T^0.2 e^{−51800/T}(1 − e^{−6000/T}) cm³ s⁻¹ (PSS83 modified to Martin et al. 1996); χ_H2 = 4.48 eV | Eqs. 29–33 p. 11 |
| instability | criterion C < 0 (from Omukai & Yoshii 2003; A, B coefficients only there). Along RHFC02's path: 3-body instability at n ≲ 2e10, t_growth/t_ff ≳ 7; CIE instability ~1.3e14–9e14 cm⁻³, t_growth/t_ff ≳ 20 (read from Fig. 9). Conclusion: fragmentation unlikely | §4, Eq. 28, Figs. 7–9 |
| 1D test | RHFC02's Lagrangian code with the transfer replaced by the local cooling of Eq. 36 (thin-corrected lines Eq. 22 + CIE with escape factor); "virtually indistinguishable" from the full code in the n–T and n–v tracks of four shells and profiles at T_c = 700, 1000, 1600, 2200, 6000 K. ICs = RHFC02's | p. 9, Figs. 4–5 |

Other printing problems flagged in the extraction (affect only the unused per-line model or the appendix): Boltzmann
factors printed e^{+E/kT}; n vs n_p in Eq. 29 vs A9/33; a sign in Eq. A12; missing s⁻¹ in A1.

**Consequence [design].** RA04's 1D test gives a **cheap local-cooling configuration that reproduces RHFC02's
collapse to T_c ≈ 6000 K**: thin HM79 lines × min(1, (n/8e9)^−0.45) + CIE × escape factor. This is the natural first
target for M2/M3 before full line transfer: it needs no transfer at all, and RA04 Figs. 4–5 provide the comparison.
The factor-2 question must be settled first (reproduce RA04 Fig. 1).

### M0.4b Ripamonti & Abel (2005) review [read-src]

arXiv:astro-ph/0507130v1, chapter 1 of a SIGRAV school volume; the book is not named in the text (believed: *Joint
Evolution of Black Holes and Galaxies*, eds. Colpi, Gorini, Haardt & Moschella 2006 [memory], to check on ADS).
Extraction: `~/Projects/popiii_lit/RA05/extraction.md`.
- H₂ cooling Eqs. 1.22–1.24 (Λ = Λ_LTE/[1 + Λ_LTE/(n_H Λ_n→0)], HM79 LTE, GP98 low-density polynomial). CIE Eq. 1.29
  = RA04 Eq. 5. Three-body: k = 5.5e-29 T⁻¹ (PSS83) and k/8 for H₂ as third body, the only source cited. Table 1.1
  (19 reactions) has unit misprints; rates are to come from the primary sources.
- **Fig. 1.19 is RHFC02's Fig. 1.** No stage times; central densities read off the figure [read-src, approximate]:
  log n_c ≈ 6.4, 8.4, 10.2, 13.6, 17.9, 20.2, 21.2, 22.2, 23.0, 23.3 for stages 0–9; peak infall ≈ 18 km/s; core edge
  ≈ 3e-3 M☉ (stage 7) and ≈ 0.06 M☉ (stage 9). Outer edge of the plot ≈ 100 M☉.
- Core at n ≳ 1e21 cm⁻³, M_core,0 ~ 3e-3 M☉, Ṁ 0.01–0.1 M☉/yr declining ∝ t^−0.3. Validity, verbatim (p. 40): "they
  do not properly trace the internal structure of the core, so their results cannot be trusted except in the very
  initial stages of accretion (say, when Mcore ≲ 0.1M⊙)."

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

### Answers (Tom, 19 Oct 2026)
- Literature hosts opened in the environment settings (not yet effective in the cloud container at the time of
  writing: still refused); the work moves to Tom's Mac, where they are reachable.
- Agreed: neutral He (as CIA partner) first, D/HD later; RHFC02's line transfer and CIE fits AND modern data, as
  switches (reproduce first, then modernise); third book `popiii/` at top level, published under `/popiii/` like
  `spectra/`; Lagrangian mass grid for the collapse (M1–M3), adaptive grid with interior mass for accretion (M4).
- Next: M0.3/M0.4 extraction from the papers, then M1.
- On the Mac: Julia 1.10.4 as pinned in `julia/env/Manifest.toml` (juliaup: `juliaup add 1.10.4`), then
  `cd julia && julia +1.10.4 --project=env -e 'using Pkg; Pkg.instantiate()'`.
- **Superseded (Tom, 19 Oct 2026): on the Mac the work uses the current Julia, 1.12.6** (juliaup default), not
  1.10.4. `julia/env/Manifest.toml` is git-ignored (`.gitignore` l. 19), so "pinned" only ever meant the cloud
  container's local Manifest. `Pkg.instantiate()` resolves ForwardDiff 1.4.6, StaticArrays 1.9.23. Baseline: M0.10.

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

## M0.10 Baseline on the Mac, Julia 1.12.6 [run]

Same commands as M0.9 (`-t 1`, the three runs in parallel as separate processes), branch `popiii` @ a40acbd, Apple
silicon (18 cores), ForwardDiff 1.4.6. Logs: `lib/popiii/m0_mac_{SOD,SUBC,SUPC}.log`, `lib/popiii/m0_mac_regress.log`.

| problem | steps (Mac, 1.12.6) | steps (cloud, 1.10.4) | tier (c) | ms/step (steady) |
|---|---|---|---|---|
| SUPC | 713 | 713 | ACCEPT | 1.8 |
| SUBC | 587 | 587 | ACCEPT | 2.0 |
| SOD | 464 | 464 | ACCEPT | 1.0 |

`popiii_regress.jl`: ACCEPT. Against the cloud log, the only differences are in the last printed digits (largest:
SOD info sample at t = 0.5, p* 0.81090162 → 0.81090289, relative 1.6e-6). The change of Julia version (and of
whatever package versions the cloud resolved; not recorded there) does not alter the verified paths.

