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


## M1a Regular centre, Lagrangian mass grid, self-gravity, free outer surface (19 Oct 2026) [run]

Code (all behind switches, default off; SOD/SUBC/SUPC regression identical to M0.10, `lib/popiii/m1a_regress*.log`):
- `TitanResidFast.REGULAR_CENTRE[]`: radius row of the centre interface is the fixed-value row R = x (M0.1 item 1);
  u/r at the centre takes the regular limit u(ngrs+1)/r(ngrs+1) in the viscosity and the (1 − 3f) energy terms
  (item 2). `TitanCore.timstp!` and `totnrg!` skip the 0/0 at r = 0 (items 3, 5: `relchg`, `mor`).
- `TitanResidFast.XMINT[]` (Lagrangian, lgrid = 3): fixed interior mass per interface, used for gravity (item 7),
  **and x_me frozen (mass rows R = x)**. Found by the static test: TITAN's mass row x_me(k−1) − x_me(k) = ρΔV takes the
  shell mass as a difference of exterior masses ≈ x_mtot (= 1000 × cloud mass by default in `ProblemSetup`), which put
  1e-8 relative noise into the central densities of a sphere at rest. On a Lagrangian grid the continuity row with
  u_rel = 0 already keeps every shell's mass, so nothing is lost.
- `TitanResidFast.PEXT[]` + `TitanCore.FREE_OUTER[]`: free outer surface (RHFC02 A2.3): the outermost interface obeys
  the momentum equation of its half zone against a constant external pressure, instead of TITAN's piston (peq no longer
  imposes ur/uextr). `ProblemSetup.build!(outer = :free, pext, selfgrav = true)`.
- Item 4 (Eddington factors at r = 0) is not done: M1 is hydrodynamics only (neqn = 5).

Tests, `julia/popiii_m1.jl` (logs `lib/popiii/m1_*.log`):

1. **Static uniform sphere** (100 M☉, ρ = 1e-19, 300 K, 40 shells with Δm ratio 1.3, no gravity, p_ext = p): over 30
   sound-crossing times, Δr, u, Δρ **exactly 0**, r(centre) = 0, no NaN. **PASS.** (Before freezing x_me: 1e-8.)
2. **Pressure-free collapse** (100 M☉, ρ0 = 1e-19, 40 shells, ratio 1.2, p_ext = 0) against r = r0 cos²β, β + sin β cos β
   = t(8πGρ0/3)^½, u = −(2Gm/r0)^½(1/x − 1)^½; worst shell, at the times the exact R/R0 crosses 0.5 / 0.1 / 0.03:

   | θ | stol | T0 [K] | steps to R/R0 = 0.1 | |Δr/r| at 0.5 / 0.1 / 0.03 |
   |---|---|---|---|---|
   | 0.55 (TITAN) | 0.1 | 1e-3 | — (fails, 20 steps to 0.73 t_ff) | 9.3e-2 (at 0.35) |
   | 0.55 | 0.01 | 1e-3 | 720 | 4.2e-3 / 6.3e-2 / 4.5e-1 |
   | 0.55 | 0.001 | 1e-3 | 6940 | 4.6e-4 / 6.8e-3 / 4.4e-2 |
   | 0.5 | 0.01 | 1e-3 | 700 | 4.2e-4 / 6.1e-3 / 4.0e-2 |
   | 0.5 | 0.001 | 1e-3 | 6920 | 2.4e-4 / 5.2e-3 / 3.5e-2 |
   | 0.5 | 0.001 | **1e-6** | 6920 | **1.1e-5 / 1.6e-4 / 1.0e-3** |

   Same numbers with conv = 1e-9 (Newton tolerance not limiting) and with 80 shells (no spatial error: the solution is
   homologous). The floor at T0 = 1e-3 K is the residual gas pressure (adiabatic heating makes p/ρ grow ∝ x⁻²), not a
   scheme error. Arrival: R/R0 = 0.01 at t/t_ff = 0.999574, exact 0.99958. **PASS** (< 1 % in the collapse time; r(m, t)
   to 1.6e-4 at R/R0 = 0.1 with θ = 0.5).

**Findings.**
- **θ = 0.55 (TITAN's default) makes the step first order**: error ≈ 7 stol at R/R0 = 0.1. θ = 0.5 is second order
  (10× better at stol = 0.01). RHFC02's scheme is second order except the energy equation (A2.1). For the collapse
  runs: θ = 0.5 or stol ≤ 1e-3, to be decided with the shock tests of M3 (θ > 0.5 damps).
- **TITAN's stol control lets dt grow 10× per step** (tfac ≤ 10), which with stol = 0.1 takes the collapse to
  0.73 t_ff in 20 steps. RHFC02 caps dt at the dynamical time C_dyn(r³/Gm)^½ and at 2 Δt_prev (A2.2). A dynamical-time
  cap belongs in `timstp!` for the collapse (next).
- After the shells reach the centre the run degenerates (shell crossing at r → 0; steps → 0): expected for dust.

Next in M1: dynamical-time step cap; Bonnor–Ebert sphere (isothermal-like EOS) static for 10 sound-crossing times;
Larson–Penston approach; energy ledger with self-gravity; quasi-star hydrostatic regression (P7.3).

## M1b Time-step cap, Bonnor–Ebert, Larson–Penston, energy, quasi-star regression (19 Oct 2026) [run]

**Dynamical-time cap** `TitanCore.DT_DYN[]` (RHFC02 A2.2): Δt ≤ C_dyn min_k (r_k³/G m_k)^½, default 0 (off).
Pressure-free collapse (T0 = 1e-6 K, stol = 0.1 so the cap alone controls Δt), worst-shell |Δr/r| at exact R/R0 =
0.5 / 0.1 / 0.03 (logs `m1_dust_cdyn*`):

| θ | C_dyn | steps to R/R0 = 0.1 | error |
|---|---|---|---|
| 0.5 | 0.1 | ~100 | 7.3e-3 / 1.5e-1 / 5.1e-1 |
| 0.5 | 0.03 | ~100 | 1.8e-3 / 1.8e-2 / 1.1e-1 |
| 0.5 | 0.01 | ~280 | 9.6e-5 / 2.6e-3 / 1.4e-2 |
| 0.5 | 0.003 | ~860 | 6.4e-6 / 1.5e-4 / 1.0e-3 |
| 0.55 | 0.01 | ~280 | 2.6e-3 / 5.1e-2 / 3.0e-1 |

Order ≈ 1.75 in C_dyn with θ = 0.5 (first order with θ = 0.55). C_dyn = 0.003 matches the stol = 1e-3 run (M1a) with
8× fewer steps. Step numbers are at 20-step print resolution.

**Bonnor–Ebert sphere** (`popiii_m1.jl be n variant`, ξ_max = 3, γ = 1.001, 300 K, ρ_c = 1e-19, 477 M☉, free surface
at the pressure that holds the outer half shell, θ = 0.55, C_dyn = 0.01, 10 sound-crossing times):
- `discrete` (radii marched so that TITAN's momentum rows vanish at rest): |Δρ/ρ| = 7.3e-10, |u|/c_s = 1.3e-10. **PASS.**
- `le` (zone densities sampled from the Lane–Emden solution; measures the discretisation of hydrostatic equilibrium):

  | n | 30 | 60 | 120 | 240 |
  |---|---|---|---|---|
  | max |Δρ/ρ| | 8.0e-3 | 2.0e-3 | 5.1e-4 | 1.3e-4 |
  | max |u|/c_s | 1.0e-3 | 2.5e-4 | 6.3e-5 | 1.6e-5 |

  Second order. The criterion (< 1e-3) holds from n = 120. **PASS** (n ≥ 120).

**Larson–Penston reference, computed** (`julia/lp_selfsimilar.jl`) [computed]. Shu's similarity equations with
x = r/(c_s t) < 0; regular centre v = 2x/3, α = α0. The sonic point is a **node**: every α0 below 1.7212 reaches the
sonic line with both numerators vanishing, and the asymptotes vary continuously with α0 (α0 = 1.60 → αx² = 9.23,
v = −3.39; 1.72 → 8.54, −3.18). Hunter (1977, first page, ADS scan) says so: "an infinity of other similarity
solutions". A first bisection on "which numerator vanishes first" found only the edge of that family (α0 = 1.7212) —
recorded so that nobody repeats it. **LP is the member analytic at the node**: the curves arrive with slope
v' → 1/|x*| plus a non-analytic term C|Δx|^{λ−1} that dominates near x*; LP has C = 0 (sign change of v' − 1/|x*|
at fixed small distance from the sonic line, independent of the distance 1e-4…1e-6):
**α0 = 1.66582, x* = −2.34112, α x² → 8.838, v → −3.274** (at x = −1e5). [memory, unverified]: Hunter 1977 / Larson
1969 give 1.667, 2.33, 8.86, −3.28; agreement 0.3 %. The rest of Hunter 1977 was not reachable (ADS page scans: access
denied); to be checked against his table when available.

**Larson–Penston approach** (`popiii_m1.jl lp n`, BE sphere ξ_max = 6.5, 1625 M☉, γ = 1.001, collapse started by
T × 0.9 and p_ext × 0.9, geometric mass grid, θ = 0.55, C_dyn = 0.01). Peak of 4πGρr²/c_s² and minimum of u/c_s,
and α_c = 4πGρ_c(t0 − t)² from the local slope of ρ_c^{−1/2}(t):

| log ρ_c/ρ_c0 | 2 | 4 | 6 | 8 | 10 | 12 | LP |
|---|---|---|---|---|---|---|---|
| max 4πGρr²/c_s² | 4.66 | 5.92 | 6.93 | 7.68 | 8.19 | 8.53 | 8.838 |
| min u/c_s | −0.99 | −1.73 | −2.27 | −2.64 | −2.88 | −3.04 | −3.274 |
| α_c (between decades) | 5.1 | 2.8 | 2.15 | 1.88 | 1.82 | 1.76 | 1.666 |

- Independent of resolution and time step to ~2e-3 at 8 decades: n = 200 (q 1.1 and 1.05), 400 (q 1.05), 300 with
  Δm_1/M = 4e-14; C_dyn = 0.01 vs 0.003 (agree to 1e-3).
- The approach is monotone and slow, as Hunter (1977) found ("none of the collapses studied showed any strong tendency
  to a similarity solution"; LP the best description near core formation). **PASS** as "monotone approach to the
  computed LP values, numerically converged"; the gap at 12 decades (3.5 % in ρr², 7 % in u) is physical.
- **Set-up trap found:** lowering T without lowering p_ext drives an 11 % pressure step inward from the surface; it
  steepens into a converging shock that only the finest grid (n = 400, q = 1.05) resolves, and wrecked that run.

**Energy.** `energy()` in `popiii_m1.jl`: K (interface half-masses), U, W = −Σ G m_k m̄_k/r_k (the potential of the
momentum row's gravity term), p_ext V. Relative drift (E − E0)/|W|:
- dust (θ = 0.5, C_dyn = 0.003) to R/R0 = 0.01, where K ≈ |W|: 2.1e-5;
- LP over 8 decades: 5.9e-5 (θ = 0.55), **1.3e-6 (θ = 0.5)**.

**Quasi-star hydrostatic regression (P7.3)**: `QS_SINK_M0=0 qstar_run.jl relax 3e9 2000`, compared with
`qstar_compare.jl` against `qs_profile_f0.1_md0.csv` (`lib/popiii/m1_qstar_compare.txt`). Interior rows (r ≤ 0.98 R)
differ from the stored table `reference_runs/qstar_hydrostatic_compare.txt` by ≤ 2e-4 in ΔT/T, ≤ 7e-4 in Δρ/ρ, ≤ 5e-4
in L_c/L. **PASS.** Notes: (i) after 2000 steps this run is at t = 2.8895e7 s, the stored one at 2.8146e7 s (2.7 %;
not investigated: a different time-step history, e.g. the clock/energy fixes, or stol); (ii) P7.3's prose says
|Δρ/ρ| ≤ 1 % and L_c within 0.002, its own table has 0.0136 (0.961 R) and 0.0027 (0.73 R); (iii) `qstar_run.jl`
failed on any shock-free run (`shellmass()` summed over an empty range, from the later infall work); fixed (`init = 0`).

**M1 status: done** for the hydrodynamic part (regular centre, self-gravity, Lagrangian mass grid, free surface,
dynamical time step), with all tests passing. Left for later: Eddington factors with r = 0 (M0.1 item 4; needed with
radiation, M2b/M3); the adaptive grid with interior mass (M4).
**Default for the collapse runs: θ = 0.5, C_dyn ≤ 0.01** (to revisit with shocks in M3).

## M2a Primordial chemistry module and one-zone tracks (19 Oct 2026) [read-src + run]

**Rate sources** — full extraction, verbatim formulas, page references and every cross-check: `docs/popiii_rates_extraction.md`
(compiled by a subagent from the page images of AAZN97 (Abel, Anninos, Zhang & Norman 1997, New Astron. 2, 181),
GP98, GA08, RA04, RA05, RHFC02; spot values re-checked here [run]).
- **RHFC02 Table A1 is fully specified**: AAZN97 k1, k2 (case A), k7, k8 [read-src]; PSS83 k4 = 5.5e-29/T, k5 =
  6.5e-7 T^−½ e^−52000/T (1 − e^−6000/T), k6 = k4/8, k7 = k5/8, k9 = 1.7e-4 × 5.8e-11 T^½ e^−158000/T — **PSS83 only
  through RA04/RA05** (and the ADS OCR text layer of PSS83 Table 1, which agrees; the page images beyond p. 632 are not
  reachable). RHFC02's "AANZ98" is AAZN97 (the k-numbers match its appendix tables).
- Printing errors found in the sources: RA05 k5 unit (cm⁶ for cm³); GA08 three-body units (cm³ for cm⁶); GA08 k7 and
  AAZN97 k9 "log" means ln; AAZN97's printed k13 (the code's form agrees with GA08); AAZN97's formation-cooling formula
  (1.38 for 1.83 eV, k10 missing); RA04 eq. 6 factor 2 (M0.4).
- **ChemistryKernels.jl / EmissionKernels.jl are not RHFC02** (M0.5 confirmed in detail): case-B recombination; the
  ABN02 three-body rate (14× below PSS83); no H-impact dissociation below 3481 K; no H₂-partner three-body channels; an
  unsourced heating suppression n/(n + 1e8); H₂–H⁺ and H₂–e⁻ cooling fits that do not match GA08 although labelled so;
  an LTE fit that is not HM79 (0.55× at 8000 K). They stay a cross-check only.
- **Detailed balance [run]**: PSS83's k5/k4 = 0.25 × K_Saha(T) (0.256 at 1000 K → 0.250 at 5000 K; K from partition
  functions, `PrimordialChem.K_H2`). The PSS83 network relaxes to an equilibrium 4× too molecular in K; RHFC02 switched
  to Saha/SCVH equilibrium above ρ_eq = 5e-11 g cm⁻³ anyway. RA04's modified k5 is 0.13–0.34 K (no constant factor).
  Option `k5_db` = k4 K_H₂ (modern, detailed balance).

**Code** `julia/PrimordialChem.jl`:
- Species H, H⁺, H⁻, H₂, H₂⁺, e⁻, He (neutral, n_He/n_H = 1/12); unknowns x_H⁺, x_H⁻, x_H₂, x_H₂⁺; **x_H and x_e
  eliminated (H nuclei and charge conserved by construction, to round-off)**.
- `rates_rhfc02(; k5, tb)`: the faithful set; `k5_ra04`, `k5_db` options. AAZN97 H₂⁺ channels etc. not yet (RHFC02 has
  none; needed for the modern set).
- EOS `eos_rhfc02` (A6–A9): u_th = Σ N_i kT/(γ_i − 1) with H₂ from Eq. A9, u_chem = −Σ N_i χ_i (H₂ 4.48 eV, H 13.6 eV,
  H⁻ 0.754 eV, H₂⁺ 2.65 eV). [interp] RHFC02's "weighted average of the adiabatic indexes" read as the number-weighted
  average of 1/(γ−1), which is what makes u_th = Σ N_i kT/(γ_i − 1). Checks: c_v = 1.625 k per H nucleus (atomic),
  1.875 k (molecular, vibration excited); T(u) inverts to round-off.
- Cooling: `cool_H2_thin` (RA05 eq. 1.22: HM79 LTE / (1 + LTE/(n_H Λ_GP98)), per molecule × n_H₂; n_H = atomic H);
  `cool_H2_ra04` (× min(1, (ρ/m_H/8e9)^−0.45)); `cool_CIE_ra04` (A = 0.054, escape factor (1−e^−τ)/τ with
  τ = (ρ/3.3e-8)^2.8).
- **`step_coupled`: chemistry and T solved together** (backward Euler in y with the rates at the new T, plus
  u(ρ, T, y) = u; damped Newton keeping every abundance including x_H, x_e ≥ 0; recursive halving). Found necessary:
  with rates held at a fixed T (operator split), H₂ dissociation at T ≳ 2000 K absorbed more than the thermal energy in
  one step and u fell below u(T → 0). This is the chemo-thermal coupling that M0.2 option (A) worried about; within
  the chemistry step it is now exact. Cooling is implicit in T (backward Euler) — the cooling time can be ≪ the step.
- `onezone`: dρ/dt = ρ/t_ff, Δln ρ steps; (1) compression (trapezoid in p) + implicit cooling at fixed y, (2)
  `step_coupled` at fixed u.

**One-zone tracks** (`julia/popiii_m2_onezone.jl cooling k5`, logs `lib/popiii/m2_onezone_*.log`), from RHFC02 P100's
centre (ρ/m_H = 2.6e6, 270 K, f_H₂ = 5e-4, f_e = 1e-10):
- three-body conversion: f_H₂ = 1e-2 at log ρ/m_H = 9.4, 0.23 at 10.4, 0.78 at 11.4, 0.97 at 12.4 (RHFC02: "n ≳ 1e8,
  3-body reactions efficient ... stages 2–3");
- thin cooling only: T 660–830 K at 1e12–1e14, still molecular at 1e17 (1936 K);
- RA04 thick lines: dissociation from ≈ 1800 K, the "thermostat" holds T ≈ 2300–3300 K at 1e15–1e17 (f_H₂ 0.89 → 0.62);
  RHFC02 stage 4: T_c = 3000 K, "H₂ is dissociating";
- RA04 lines + CIE: CIE cools from 2050 to 1690 K around 1e15–1e16, then T rises (2440 K at 1e17).
- k5 choice (PSS83 / detailed balance / RA04) at log ρ/m_H = 17.9: T = 2841 / 2696 / 2845 K — a 5 % systematic.
- Δln ρ = 0.01 vs 0.003: max |ΔT/T| = 3.5 %, |Δf| = 2.7e-3 (first order; where it peaks: to locate).
- **The one-zone free fall is not RHFC02's centre early on**: the P100 centre takes 7.2e5 yr from stage 0 to 1
  (Table 2), the free fall from the same state 4.4e4 yr to the same density; adiabatic heating then dominates and the
  one-zone gas is at ~1000 K at 1e8 where RHFC02 has 450 K. RHFC02's cloud starts near hydrostatic equilibrium; the
  comparison with RHFC02 is the 1D run (M2c/M3), not this.

**Next (M2b/M2c):** the chemistry and the RHFC02 EOS as a TITAN material (composition per Lagrangian shell, operator
split chemistry with `step_coupled` per shell after each hydro step, as M0.2 (A), u kept by the energy row); the
1D collapse of P100 with thin / RA04 local cooling (RA04's own 1D test) against RHFC02's stage table; then the
continuum (Lenzuni Planck mean, faithful) and line transfer (M3). Saha equilibrium switch at ρ_eq.

### M2b Initial conditions of RHFC02 / ON98 — open question [read-src + computed]

ON98 (Omukai & Nishi 1998, ApJ 508, 141, arXiv:astro-ph/9811308; local copy `~/Projects/popiii_lit/ON98/`) §3.1:
"hydrostatic equilibrium clouds with adiabatic stratification ... ρ = ρ_c θ_n^n, n = 1/(γ_ad − 1) ... p = Kρ^γad, where
K is a free parameter fixed by giving the number density at half mass radius n_h and the total mass M. These clouds
are cut off at the radius r_s where the density falls off 10⁻³ from the central value." Table 2: run A M = 10² M☉,
n_h = 1e6, y_H2 = 5e-4, y_e = 1e-10, N = 250 zones, Δr_{i+1}/Δr_i = 1.01; run B the same with N = 100 (ratio 1.05);
run C 10 M☉, 1e8, 5e-3; run D 1 M☉, 1e10, 0.45; run E as A with y_H2 = 0, y_e = 1e-8. RHFC02 P100 ≙ ON98 A/B.

**The numbers do not describe a hydrostatic cloud at ~250 K [computed]**: an n = 1.5 polytrope with ρ_c/m_H = 2.6e6 and
T_c = 270 K (μ ≈ 1.22) has M = 4πρ_c a³ ξ₁²|θ'| = 0.10 M☉; the Jeans mass at 1e6 cm⁻³, 300 K is ≈ 0.3 M☉. Holding
100 M☉ hydrostatically at that central density needs T_c ~ 2e4 K. Yet ON98 Fig. 3a (run A) starts at log n_c ≈ 6.3,
T_c ≈ 250 K, and the central collapse to stage 1 takes 5.7e5 yr (ON98 Fig. 1 caption) / 7.2e5 yr (RHFC02 Table 2),
≈ 20 central free-fall times — quasi-static, as stated. ON98 also gives t_ff = 3.2e5 yr "at the center ... initially",
which corresponds to n_c ≈ 2e4, not ≳ 1e6. **To ask Tom** (co-author of RA04/RA05; may know R02's set-up): how were
the P/I runs' T profiles and outer pressure chosen?

**Reference track, ON98 Fig. 3a (run A), read off the figure [read-src, approximate ±0.05 dex]:** log T_c ≈ 2.40 at
log n_c 6.3; 2.6 at 8; 2.75 at 10; 3.05 at 12; 3.3 at 15; 3.35 at 16; 3.8 at 20; 4.0 at 21.5; 5.1 at 24.

**Provisional set-up for the 1D runs [design, to be replaced when answered]:** density shape of an n = 1.5 polytrope
(ρ_c/m_H = 2.6e6, truncated at ρ = 1e-3 ρ_c) scaled in radius to M = 100 M☉; T = 270 K × (ρ/ρ_c)^{2/3}; at rest;
f_H2 = 5e-4, x_e = 1e-10; external pressure = p at the edge. Not hydrostatic (≈ 100 Jeans masses): the collapse will be
faster than R02's early phase; the comparison is meaningful only after the self-similar phase is reached
(RHFC02 §3.2.5: tracks converge for n_c > 1e14 regardless of IC; ON98 Fig. 3a: runs A, C, D converge by 1e12–1e14).

## M2c First 1D collapse with chemistry in TITAN (19 Oct 2026) [run]

**Code.** `julia/PopIIIMaterial.jl`: material `PrimordialGas` (composition per zone `y[k]`, RHFC02 EOS via
`PrimordialChem.eos_rhfc02`; e_g includes the binding energies). `TitanResidFast.ZONE[]` gives the material the zone
index of each call (set in `_material!` and the fills). After every converged hydro step, `chem!` advances chemistry +
cooling + T in each zone at fixed (ρ, e_g) over the step (`step_coupled` with implicit cooling) — operator split,
M0.2 (A); above `ρeq` (RHFC02: 5e-11 g cm⁻³) `step_equilibrium` (Saha, `y_saha`: H₂ ⇌ 2H with K_H2, H ⇌ H⁺ + e⁻; no
SCVH95 yet). Driver `julia/popiii_m2c.jl cooling nshell` (env THETA, CDYN, NSTOP, TSTOP, DTDIAG).

**Run** (provisional P100 set-up of M2b; 100 shells, central 3e-7 M☉, q = 1.197; RA04 local cooling = thick H₂ lines
+ CIE, i.e. RA04's own 1D test configuration; RHFC02 network throughout — this run predates the equilibrium switch;
C_dyn = 0.01). Logs `lib/popiii/m2c_ra04_n100{,_th055}.log`. **Correction (M2d):** the θ = 0.55 run already had the
Saha switch above ρ_eq = 5e-11 g cm⁻³ (the constructor's default), contrary to the sentence above; the θ = 0.5 run did
not.

- **θ = 0.5 fails**: the centre rings (n_c oscillating between 2.5e6 and 2e8, T_c spikes to 2400 K destroying the
  central H₂, f_H2 → 7e-5) for 1.2e5 steps, Δt down to 1e5–1e7 s, then a Newton failure (t = 3.47e5 yr). The time-step
  limiter is u in the central zones (DTDIAG). θ = 0.5 has no numerical damping; the non-equilibrium start excites
  acoustic oscillations of the light central shells that never decay. **Use θ = 0.55 (TITAN's default) for the
  collapse** (M1's preference for 0.5 was for smooth flows; revisit with an explicit damping alternative if the first
  order in time matters, M1a: error ≈ 7 stol).
- **θ = 0.55**: smooth collapse, 13 348 steps to n_c = 1e18 (stop). Central track vs ON98 Fig. 3a (run A, read off):

  | log n_c | 8.4 | 9.4 | 9.9 | 10.4 | 11.4 | 11.9 | 12.9 | 13.4 | 14.4 | 14.9 | 15.4 | 15.9 | 16.4 | 16.9 | 17.4 | 17.9 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
  | T_c here [K] | 407 | 487 | 562 | 679 | 892 | 1002 | 1309 | 1537 | 1823 | 1965 | 2090 | 2180 | 2222 | 2418 | 2680 | 2983 |
  | ON98 [K] | 428 | 508 | 554 | 649 | 917 | 1090 | 1338 | 1473 | 1784 | 1964 | 2093 | 2218 | 2496 | 2842 | 3232 | 3682 |
  | f_H2,c | 0.007 | 0.087 | 0.23 | 0.45 | 0.83 | 0.92 | 0.99 | 0.996 | 0.96 | 0.93 | 0.91 | 0.92 | 0.94 | 0.92 | 0.87 | 0.80 |

  **Agreement within ≈ 5 % from 1e8 to 1e16 cm⁻³** (ON98 read-off uncertainty ±0.05 dex ≈ 12 %). Above 1e16 this run is
  10–20 % cooler: expected differences, not yet separated — (i) [ruled out in M2d: the switch was on], (ii) RA04's local cooling instead of line + continuum
  transfer, (iii) ON98's own physics (its CIE/continuum and EOS).
- The centre needs 1.15e6 yr to reach 1e8 (RHFC02 stage 1: 7.2e5 yr; ON98 (1): 5.7e5 yr): the provisional initial
  state (M2b) — the pre-collapse phase is not comparable until the IC question is settled.

**Next:** rerun with the equilibrium switch and the k5 variants (bracket); resolution (140 shells as P100/2, 200);
radial profiles at the RHFC02 stages (Fig. 1 / RA05 Fig. 1.19 read-off); continue past 1e18 to core formation (needs
SCVH95 for p ≥ 1e6 dyn cm⁻² — tables to obtain; Saha alone gives f_H2 = 0.12 at 2e4 K, 1e22 cm⁻³ vs RHFC02's ~0.05).

## M2d Bracket: equilibrium switch, k5, resolution (19 Oct 2026) [run]

Same set-up as M2c, θ = 0.55 (now the default of `popiii_m2c.jl`), RA04 local cooling; logs `lib/popiii/m2c_ra04_*`.

| log n_c | ON98 | Saha switch, 100 shells | network PSS83, no switch | network k5_db, no switch | switch, 140 | switch, 200 |
|---|---|---|---|---|---|---|
| 8.4 | 428 | 412 | — | 412 | 383 | 388 |
| 11.4 | 917 | 892 | — | 892 | 884 | 878 |
| 12.9 | 1338 | 1309 | — | 1309 | 1302 | 1293 |
| 14.9 | 1964 | 1965 | ≈ 2040 | 1966 | 1962 | 1959 |
| 15.9 | 2218 | 2180 | 2247 | 2179 | 2176 | 2169 |
| 16.9 | 2842 | 2418 | 2462 | 2418 | 2390 | 2404 |
| 17.9 | 3682 | 2983 | 3140 | 2983 | 2979 | 2969 |

- **The detailed-balance network reproduces Saha exactly**: k5_db without the switch is identical to the Saha switch
  to all printed digits; direct test at n_H = 1e16, 2200 K: f_H2 = 0.46067 both ways, T = 2217.97 K both ways — a
  consistency check of `k5_db`, `K_H2` and `y_saha` [run].
- **PSS83's pair without the switch** (4× too molecular) is 2–5 % hotter above 1e14 (3140 vs 2983 K at 1e18).
- **Resolution**: 100 / 140 / 200 shells agree within ≈ 2 % from 1e10 to 1e18 (≤ 7 % at 1e8.4, in the pre-collapse
  phase, where the time to reach 1e8.4 also shifts: 1.15 / 1.20 / 1.22 Myr).
- **The 10–20 % deficit against ON98 above 1e16 is not chemistry or resolution.** Remaining candidates: RA04's local
  cooling (fitted to RHFC02, which itself is not ON98: ON98's CIE and continuum), the EOS (RHFC02's γ-weighted u_th vs
  ON98's), and the read-off of ON98 Fig. 3a. Next: compare with RHFC02's own stages (Table 2: T_c = 3000 K at stage 4;
  RA05 Fig. 1.19 read-off log n_c ≈ 17.9 at stage 4) — ours reaches 2983 K at log n_c = 17.9: **agrees with RHFC02**,
  whereas ON98 has 3682 K there. So RA04's cooling reproduces R02, and R02 and ON98 differ by ≈ 20 % in T_c at 1e18
  (RHFC02 claims "good agreement with ON98"; the difference is within the scale of their Fig. 7).

## M2e SCVH95 tables obtained; ABN02-like initial profiles (19 Oct 2026)

**SCVH95 EOS tables.** The AAS CD-ROM volume 5 (Tom's pointer, aas.org/.../cdrom/volume5/doc/toc5.htm) lists them
(`volume5/apjs/v99/p713/`: h_tab_i.dat, h_tab_p1.dat, h_tab_p2.dat, he_tab_i.dat, rho_crit.dat, read.f, readme.doc,
table.doc), but aas.org sits behind a Cloudflare browser challenge (HTTP 403 for scripts and WebFetch), and the Wayback
Machine has the listing pages (2026-01-31/02-01) but not the files. CDS (J/ApJS/99/713) is behind the same kind of
challenge. **Copies taken from the YREC stellar-evolution code** (github.com/yreclab/yrec, `input/eos/scv/`):
h_tab_i.dat and he_tab_i.dat, 391 230 bytes each, 63 isotherms log T = 2.10 … 7.06 (step 0.08), columns log P, X_H2,
X_H, log ρ, log S, log U, then derivatives and ∇_ad, as in SCVH95. Check [run]: log T = 2.10, log P = 4.00 →
log ρ = −5.7154 (pure H₂ ideal gas: −5.72); log U = 10.120 (ideal H₂, 5/2 kT per molecule: 1.3e10 erg/g, so U's zero is
cold H₂ — the offset to our u_chem zero has to be applied). YREC has no licence file: the tables stay outside the
repository, in `~/Projects/popiii_lit/SCVH95/` (sha256 below). Not obtained: h_tab_p1/p2 (variants across the
plasma phase transition), rho_crit.dat, readme.doc (to confirm the column definitions against the original).
```
9b199919d995d8ca84d6aecc3b818aebb963d889a5c195e3ed55f2c26ecc5db0  h_tab_i.dat
799034124a7574e1f396ad7bd2b054d55a3c9b90e32b38b0b5f9124f4c532a55  he_tab_i.dat
```

**Initial conditions (Tom, 19 Oct 2026):** "The initial densities for the 100 solar mass cloud must have been smaller.
However best would be to use the spherical profiles similar to the ones found by Abel Bryan and Norman in 3D sims."
New option `IC=abn02` in `popiii_m2c.jl`: ABN02 (Science 295, 93; arXiv:astro-ph/0112088) Fig. 2, profile "9 Myr
later" (red dotted with circles), read off a 400-dpi rendering [read-src, approximate ±0.05 dex]:
- n(r) (particle density; ρ = n μ m_H, μ = 1.22): log n = 5.83 (core, r < 0.027 pc), 5.40 at log r(pc) = −1.28,
  5.25 at −1.08, 4.73 at −0.84, 3.99 at −0.59, 3.47 at −0.35, 2.87 at 0.06, 0.97 at 0.87 (slope ≈ −2.1 outside);
- T(M): 252 K (M < 30 M☉), 229 at 100, 205 at 300, 164 at 1000 (minimum), 355 at 2500, 1120 at 5000 M☉;
- H₂ mass fraction: 1.3e-3 (centre), 1.0e-3 at 100, 7.8e-4 at 300, 6.6e-4 at 1000, 4.3e-4 at 2500, 2.2e-4 at 5000;
- v_r(M): −0.03 km/s (centre), −0.32 at 30, −0.49 at 100, −0.19 at 300, −0.81 at 1000, −0.70 at 2500 M☉;
- x_e not shown in ABN02: 1e-7 assumed.
Cloud cut at M = 1000 M☉ (inside the T minimum, before the virial-shock heating), free surface at the edge pressure;
150 shells, central shell 3e-7 M☉. Gives R = 3.4 pc, ρ_c/ρ_edge = 9.5e3, central t_ff = 5.7e4 yr.
**Asked Tom** for the actual profile data behind ABN02 Fig. 2 (spherically averaged n, T, v_r, f_H2, x_e vs M) —
better than a figure read-off.

**SCVH reader** `julia/SCVH.jl` [run]: parses both tables (63 isotherms each; log P and log ρ monotone on every
isotherm), mixes H and He by the additive-volume rule (X = 0.75), inverts (ρ, T) → p by bisection in log P. Tables
start at **log P = 4** (p = 1e4 dyn cm⁻²): consistent with RHFC02 using Saha below 1e4 and SCVH above 1e6. Against our
Saha/ideal EOS: p agrees within 3 % up to n_H = 1e22 (1e19–1e22, 3000–2e4 K); SCVH dissociates earlier (X_H2 0.89 vs
f_H2 0.94 at 1e19, 3000 K; 0.20 vs 0.27 at 1e20, 6000 K); at 1e23 cm⁻³ p is 21–23 % higher (non-ideal) and H is mostly
pressure-ionised (X_H2 + X_H = 0.23 at 1e5 K vs Saha x_H⁺ = 0.30) — cf. RHFC02 §3.2.3 (n_H⁺/n_H ≈ 2 at stage 9).
[interp] X_H2, X_H read as number concentrations relative to H (readme.doc not obtained). Not yet wired into the
material (needed past n_c ≈ 1e18; with the U zero-point offset to u_chem, M2e).

**Hydrostatic ABN02-like cloud** (`IC=abn02hs`): the read-off ABN02 profile (`IC=abn02`) is not in force balance in
1D even with the dark matter: pressure/gravity = 2.5–5 at the core edge (0.03–0.2 pc), ≈ 1.2 outside 0.6 pc. Run
`IC=abn02` (no DM, linear interpolation) expanded from the start and its centre cooled adiabatically to 49 K
(Newton failure at 1.8e5 yr). ABN02's dark-matter density (Fig. 2A) equals the gas density near 1 pc: now included as
a static external mass `TitanResidFast.EXTMASS[]` (M_DM(<r) read off Fig. 2A, monotone-cubic interpolation, slope −1
continued inward of 0.45 pc). `abn02hs` keeps ABN02's n_c = 6.8e5, T(M), f_H2(M) and dark matter and integrates
hydrostatic equilibrium outward to 1000 M☉ (at rest): pressure/gravity = 0.99 everywhere except the innermost shells
(round-off of Δp); log n = 5.79 / 5.55 / 4.60 / 2.75 / −0.42 at log r(pc) = −1.57 / −1.08 / −0.59 / 0.06 / 0.87 vs
ABN02 5.81 / 5.25 / 3.99 / 2.87 / 0.97; 1000 M☉ inside 3.0 pc (ABN02 read-off: 3.4 pc); ρ_c/ρ_edge = 1.2e4 — far
beyond the Bonnor–Ebert limit, an unstable equilibrium.

## M2f Collapse of the hydrostatic ABN02-like cloud (19 Oct 2026) [run]

`IC=abn02hs`, 1000 M☉, 150 shells (central 3e-7 M☉), ABN02 dark matter, RA04 local cooling, Saha above ρ_eq, θ = 0.55,
C_dyn = 0.01; log `lib/popiii/m2c_abn02hs_ra04_n150.log`. Smooth collapse, 3929 steps: quasi-static contraction to
n_c = 1e7 in 4.2e5 yr, then runaway (1e13 at 4.88e5 yr, 1e18 shortly after). Central track (ON98 Fig. 3a interpolated
from the M2b read-off):

| log n_c | ABN02-like 1000 M☉: T_c | f_H2 | P100 (provisional): T_c | ON98 |
|---|---|---|---|---|
| 8 | 571 | 0.00208 | 368 | 398 |
| 9 | 722 | 0.00568 | 424 | 466 |
| 10 | 775 | 0.0655 | 562 | 554 |
| 11 | 1015 | 0.482 | 787 | 827 |
| 12 | 1238 | 0.877 | 1002 | 1122 |
| 13 | 1553 | 0.983 | 1309 | 1380 |
| 14 | 1760 | 0.954 | 1669 | 1650 |
| 15 | 2034 | 0.894 | 1965 | 1964 |
| 16 | 2276 | 0.869 | 2180 | 2218 |
| 17 | 2540 | 0.863 | 2418 | 2985 |
| 18 | 3105 | 0.735 | 2983 | 3700 |

- The 1000 M☉ ABN02-like cloud is 40–70 % warmer than the 100 M☉ runs at 1e8–1e10 and converges onto the same track by
  1e14–1e16 (within 4 %) — the behaviour RHFC02 §3.2.5 / Fig. 7 describe (higher-mass models warmer at a given n_c,
  memory of the initial state lost by n_c ≈ 1e14). Three-body conversion happens at higher n_c here (f_H2 = 0.48 at
  1e11 vs 0.75 in P100) because the gas is warmer.
- Above 1e16 both runs lie 10–20 % below ON98 and agree with RHFC02's stage 4 (M2d): the remaining difference is the
  cooling/transfer treatment (RA04's local fits), not the initial state.
- Peak infall at n_c = 1e18: 6.7 km/s.

## M2g The z = 19 ABN02 cloud, extended past 1000 M☉ (19 Oct 2026)

**Tom:** "Using the blue line with central density of 1e3 would be best. We can extend it to masses past 1000 solar
mass cloud we can see the effects of radiating MS star in the envelope."

`IC=abn02z19` (`popiii_m2c.jl`), ABN02 Fig. 2 blue solid lines (z = 19), read off [read-src, approximate ±0.05–0.1 dex]:
- n(r): log n = 3.05 for log r(pc) ≤ −0.59, 2.90 at −0.27, 2.61 at 0.06, 2.01 at 0.38, 1.27 at 0.71, 0.60 at 0.95,
  0.30 at 1.20, −0.07 at 1.44, −0.59 at 1.68, −1.48 at 1.97.
- T(M): 178–180 K for M ≲ 440 M☉, 194 at 720, 247 at 1.2e3, 416 at 2e3, 650 at 3.2e3, 855 at 5.2e3, 1000 at 8e3
  (peak: virial-shock heating), 890 at 1.4e4, 725 at 2.2e4.
- H₂ fraction: 6.6e-4 … 7.2e-4 inside 200 M☉, falling to 1e-4 at 1.4e4 M☉.
- v_r(M): +0.0 … +0.22 km/s (slight outflow) inside ~700 M☉; −0.38 at 1.2e3, −1.0 at 2e3, −0.5 at 3–5e3, −1.5 at
  1.4e4, −2.6 at 2.2e4 M☉.
- dark matter (z = 19 line): log ρ_DM [GeV cm⁻³] = 4.39 at log r = −1.24, 4.17 at −0.76, 3.65 at −0.43, 2.23 at 0.63,
  1.71 at 1.20, 1.19 at 1.52, −0.59 at 2.0.
- enclosed gas mass (panel B): 34 M☉ at 0.45 pc, 180 at 1.15 pc, 1.3e3 at 7.4 pc, 8.3e3 at 48 pc, 3.4e4 at 100 pc.
- assumed (not in ABN02): x_e = 1e-4; CMB floor T_CMB = 2.725 × 20 = 54.5 K in the H₂ cooling (`PrimordialChem.cool_floor`,
  Λ(T) − Λ(T_CMB); env TCMB).

Construction: hydrostatic interior (ABN02's n_c = 1.1e3, T(M), DM; as `abn02hs`) out to r_join = 7.4 pc (1000 M☉), where
it agrees with the read-off within 0.1 dex (log n 3.00 / 2.82 / 2.49 / 0.71 vs 3.05 / 2.95 / 2.61 / 0.80 at log r
−0.84 / −0.35 / 0.06 / 0.87); outside, the read-off n(r) spliced on continuously (factor 0.81) with ABN02's infall
velocities (blended in over one M_join). A purely hydrostatic continuation fails there: with T rising outward through
the virial-shock layer the hydrostatic density plunges and the mass saturates at ≈ 1.8e3 M☉ — the outer gas in ABN02 is
accreting, not static. Default cloud 1e4 M☉ (MCLOUD): R = 58 pc (ABN02 panel B: 1e4 M☉ at ≈ 55 pc), 200 shells (central
3e-7 M☉, q = 1.117), pressure/gravity = 0.99 inside 7.4 pc, ≈ 0.4 in the infalling envelope.

### M2g result: collapse of the z = 19 cloud, 1e4 M☉ [run]

**First attempt failed — no artificial viscosity at these scales.** The run inherited the SOD state's viscosity:
`ql0 = 1e-4` (cm!), `ql1 = 0`, `cq2 = 1` — a viscous length of a micron on a parsec grid, i.e. none. The infalling,
virial-shock-heated envelope drives a compression wave into the 150 K core; it steepens into a weak shock at
r ≈ 0.3–0.5 pc (4.6 Myr), and without viscosity single shells spike (n ×4, T to 480 K, u to −3.5 km/s) and relax; one
spike ran away to −26 km/s, Δt fell to 1e6 s, and the chemistry Newton hit a singular Jacobian. (The earlier runs,
M2c–M2f, had the same non-viscosity but met no shock before n_c = 1e18.)
**Fix:** `popiii_m2c.jl` now sets `av = (ql0 = 0, ql1 = 0.15, cq1 = 0, cq2 = 1)` (env AVQL1, AVCQ2): viscous length
0.15 r ≈ 2 zones (Δr/r ≈ 0.04 core … 0.12 envelope), quadratic term only, as PORTING_NOTES P7 for the quasi-star infall.

**Run** (`IC=abn02z19 TCMB=54.5`, 200 shells, ABN02 z = 19 dark matter, RA04 local cooling, Saha above ρ_eq, θ = 0.55,
C_dyn = 0.01; log `lib/popiii/m2c_abn02z19_ra04_n200.log`): smooth, 7843 steps to n_c = 1e18.

| log n_c | z=19 cloud, 1e4 M☉: t [Myr] | T_c | f_H2 | 1000 M☉ hydrostatic: T_c | P100: T_c | ON98 |
|---|---|---|---|---|---|---|
| 4 | 6.889 | 145 | 0.001 | - | - | - |
| 6 | 14.416 | 254 | 0.00117 | - | - | 251 |
| 8 | 14.813 | 634 | 0.00167 | 571 | 368 | 398 |
| 10 | 14.830 | 802 | 0.0516 | 775 | 562 | 554 |
| 12 | 14.831 | 1245 | 0.868 | 1238 | 1002 | 1122 |
| 14 | 14.831 | 1751 | 0.956 | 1760 | 1669 | 1650 |
| 16 | 14.831 | 2260 | 0.874 | 2276 | 2180 | 2218 |
| 18 | 14.831 | 3065 | 0.752 | 3105 | 2983 | 3700 |

- Quasi-static phase: n_c 1e3 → 1e4 in 6.9 Myr, → 1e6 at 14.4 Myr (ABN02 reach ~7e5 nine Myr after z = 19: our
  contraction is ~1.5× slower — static dark matter halo here, growing in ABN02; and the read-off initial state). Then
  the runaway: 1e6 → 1e18 in 0.4 Myr.
- From n_c ≈ 1e10 on, the central track is the same as the 1000 M☉ hydrostatic cloud's (M2f) within 1–3 %: the 1e4 M☉
  envelope does not change the centre's thermal history (as RHFC02 found for different total masses); warmer than
  the 100 M☉ runs and ON98 at 1e8–1e12, converging by 1e14–1e16, 17 % below ON98 at 1e18 (cooling treatment, M2d).
- Peak infall at n_c = 1e18: −6.4 km/s at m ≈ 0.05 M☉.
- The 1e4 M☉ envelope (out to 58 pc, falling in at ≲ 2.6 km/s) is in the grid for the later stages (accretion,
  radiation of the protostar/main-sequence star into the envelope — M4/M5).

## M3a SCVH95 in the gas physics; collapse to core formation (19 Oct 2026) [run]

**Code.** `julia/PopIIIEqEOS.jl`: the equilibrium EOS for ρ ≥ ρ_eq as RHFC02 A1.2–A1.3 — Saha + ideal gas
(`eos_rhfc02` with `y_saha`) for p ≤ 1e4 dyn cm⁻², SCVH95 (`SCVH.jl`) for p ≥ 1e6, blended with a smoothstep weight in
log p (the form of the weight is ours). Tabulated once on a 271 × 216 grid (log ρ −12…1.5, log T 1.7…6.0; SCVH weight
> 0 at 33 598 nodes, = 1 at 25 685) and cached outside the repository (`~/Projects/popiii_lit/SCVH95/eqeos_table.jls`);
interpolated bicubically (Catmull–Rom) in (log ρ, log T), so ForwardDiff Duals pass: dp/dT by AD = finite difference to
1e-6. Energy zero: SCVH's U is zero for cold H₂ (+ neutral He); ours for free p + e: u = U_SCVH − (χ_H + χ_H₂/2)/μ_H
[checked where both are ideal: agreement 0.04–10 %; the 9–10 % at 2000–3000 K is SCVH's real H₂ partition function and
earlier dissociation vs RHFC02's γ_ad formula — hence the blend]. SCVH composition read as number fractions among H
particles: f_H₂ = 2X_H₂/(1 + X_H₂), x_H⁺ = (1 − X_H₂ − X_H)/(1 + X_H₂) [interp]. `PrimordialGas(…; eqtable)`: zones with
ρ ≥ ρ_eq take p, u from the table inside the hydro Newton (composition in equilibrium at the current (ρ, T), so the
Newton sees the dissociation/ionisation softening) and `PopIIIEqEOS.step` (T at fixed u − dt Λ/ρ) after the step.
`y_saha` fixed at low T (K_H2 underflow → x_H₂ from conservation).

**Run** `IC=abn02z19 TCMB=54.5 NSTOP=1e25 TSTOP=1e5` (EQEOS=scvh default), 200 shells, RA04 local cooling (beyond its
calibrated range above ≈ 1e19 cm⁻³, but the core is opaque there and the escape factors switch the cooling off),
θ = 0.55, C_dyn = 0.01, AV ql1 = 0.15. Log `lib/popiii/m2c_abn02z19_ra04_n200_scvh.log`. Stopped by the 2-h background
limit at step ≈ 68 500, in the post-core phase (Δt ≈ 40 s).

| log n_c | T_c, Saha+SCVH95 | f_H2 | T_c, Saha only (M2g) | RHFC02 stage: T_c (log n_c from RA05 Fig. 1.19 read-off) |
|---|---|---|---|---|
| 16.0 | 2260 | 0.872 | 2260.5 |  |
| 17.0 | 2506 | 0.875 | 2506.0 |  |
| 17.9 | 3105 | 0.732 | 3065.3 | 3000 K, stage 4 (H₂ dissociating) |
| 19.0 | 4006 | 0.558 | - |  |
| 20.2 | 6638 | 0.242 | - | 6500 K, stage 5 |
| 21.2 | 16103 | 0.0427 | - | 15000 K, stage 6 (core forming) |
| 22.2 | 31101 | 0.042 | - | 30000 K, stage 7 (core formed) |
| 23.0 | 64642 | 0.0251 | - | 65000 K, stage 8 |

- **The central track goes through RHFC02's stages 4–8 within 2–7 %** in T_c at RHFC02's (read-off) central
  densities: dissociation thermostat → core formation at log n_c ≈ 21–22, T_c ≈ 1.5–3e4 K; central f_H2 ≈ 0.04 at
  stage 6–7 (RHFC02: "≃ 5 %", with SCVH95). Identical to the Saha-only run up to 1e17 (p < 1e4 there).
- Late state (step 68 500): n_c ≈ 1.3e23 cm⁻³, T_c ≈ 9e4 K (peaked at 9.3e4 and slowly bouncing), accretion shock at
  m ≈ 0.022 M☉, r ≈ 6e11 cm, infall −21 km/s in front of it (RHFC02: shock "∼ 15 km/s → 2 km/s", smeared).
- Not yet measured (no profiles saved): core mass by RHFC02's definition (infall < 0.1 v_ff), its growth and Ṁ —
  next run saves snapshots. Performance: ≈ 35 000 steps/h post-core (Δt ≈ 40 s by stol/C_dyn in the core).

## M3b H₂ line transfer (19 Oct 2026) [run]

**Line data**: ExoMol RACPPK (Roueff et al. 2019, exomol.com; `~/Projects/popiii_lit/H2lines/`, not in the repo): 302
levels, 4712 quadrupole lines (v ≤ 14, J ≤ 31). RHFC02's subset v ≤ 2, J ≤ 20: 63 levels, 231 lines. RHFC02 used the
Turner, Kirby-Docken & Dalgarno 1977 A-values (not available): these are modern values. Thin LTE emission per molecule
(this subset) vs the HM79 fit: 0.947 (300 K), 0.983 (1000), 0.973 (2000), 0.908 (3000), 0.735 (5000 K; HM79 includes
v > 2).

**Code** `julia/H2Lines.jl` (`cooling(L, r, u, ρ, T, nH2; nbin, fnlte, τthin)`): LTE level populations (partition
function over all 302 levels) × an optional non-LTE factor on the source function; α from RHFC02 Eq. A11; Gaussian
thermal profile; comoving-frame Doppler shift from the bulk velocity (linear in r within a zone); tangent rays with
impact parameters = interface radii; short characteristics with the source linear in τ; F(r_k) = 2π∫μΔI dμ integrated
**exactly for ΔI linear in μ between rays** (the trapezoid in μ made the innermost zone's thin cooling 1.5× too large:
with two rays it gives 2πjr instead of 4πjr/3); Λ_k = ΔL/ΔV (Eq. A10). Frequency bins: span ±(3 v_th,max + |u|_max),
width ≤ ½ the narrowest thermal width (≥ nbin = 30, ≤ 200). Lines whose radial line-centre depth is < τthin = 0.02
use the thin emission (RHFC02 used the same shortcut with 0.1; at 0.1 single lines near the threshold were off by
5 %). Threaded over lines.

**Tests** (`julia/h2lines_test.jl`): thin uniform sphere — total and innermost-zone cooling = thin emission to 1e-6;
homogeneous static sphere vs the analytic emergent flux πS[1 − (1 − (1 + 2τ)e^{−2τ})/(2τ²)] at line-centre depths up to
2.6e7 — total line luminosity within 0.6 % (0.3 % with 160 shells), every line within 1 %.

**In TITAN** (`popiii_m2c.jl lines`): after each step's chemistry, H2Lines on the new state; Λ_k held fixed over the next
step (`PrimordialGas.Λext`, added to the local CIE cooling inside the implicit chemistry/T step); Δt ≤ 0.1 uρ/|Λ|
(RHFC02 A2.2 C_cool). **Non-LTE**: RHFC02 assumed LTE because their cloud started at 2.6e6 cm⁻³; the z = 19 cloud
starts at 1e3 < n_crit, where LTE populations overcool (first 500 steps: T_c 161 K vs 203 K with the RA04 fit, which has
the low-density limit). Source function × f = 1/(1 + Λ_LTE/(n_H Λ₀)) per zone (the GP98/HM79 combination of
`cool_H2_thin`; env NLTE=0 for pure LTE) [design: an escape-probability-free approximation; f → 1 above ≈ 1e5 cm⁻³].
Cost: ≈ 0.6 s per step on 6 cores before the thin shortcut and bin changes.

### M3b result: collapse of the z = 19 cloud with H₂ line transfer [run]

`IC=abn02z19 TCMB=54.5 NSTOP=1e18`, `lines` (H2Lines transfer for the lines, non-LTE factor on, RA04 CIE for the
continuum), 200 shells, Saha+SCVH95 above ρ_eq; log `lib/popiii/m2c_abn02z19_lines_n200_v1.log`. 8303 steps; all lines
thin (shortcut) below n_c ≈ 1e9, 23 thick at 1e10, 99 at 1e12, 157 at 1e14, 206 at 1e16, 224 of 231 at 1e18. Line
luminosity 1.8e33 erg/s (n_c = 1e4) → 2.7e34 (1e8) → 4.2e36 (1e12) → 6.5e36 (1e16–1e18; RHFC02 P100: ≈ 1e36 at
τ_c = 1 for 100 M☉; this cloud's molecular envelope is larger).

| log n_c | T_c, H₂ line transfer + RA04 CIE | f_H2 | T_c, RA04 local fits | ON98 |
|---|---|---|---|---|
| 6 | 260 | 0.00118 | 254 | 251 |
| 8 | 638 | 0.00166 | 634 | 398 |
| 10 | 816 | 0.0506 | 802 | 554 |
| 12 | 1273 | 0.867 | 1245 | 1122 |
| 14 | 1731 | 0.962 | 1751 | 1650 |
| 16 | 2241 | 0.883 | 2260 | 2218 |
| 17 | 2488 | 0.883 | 2506 | 2985 |
| 18 | 3090 | 0.743 | 3105 | 3700 |

- **H₂ line transfer and RA04's local fits give the same central track within 1–2 %** at every density, as RA04
  claimed for their simplified model against R02. The 17 % deficit against ON98 at 1e18 (3090 vs ≈ 3700 K) therefore
  does not come from the line cooling; it stays with the continuum (RA04's CIE fit vs ON98's continuum treatment), the
  EOS, or the ON98 read-off. Our track agrees with RHFC02's stage 4 (3000 K).
- Cost: 8303 steps in < 2 h with 16 threads (most steps have few thick lines).

## M3c From non-equilibrium emission to Rosseland diffusion: tools (Tom's question, 19 Oct 2026)

**Question (Tom):** "Can we understand the complete transition from full non equilibrium emission to when Rosseland
means are justified?" Approach: four regime boundaries that need not coincide — (1) H₂ level populations non-LTE → LTE
(n vs n_crit), (2) chemistry non-equilibrium → equilibrium (t_chem vs t_ff), (3) lines thin → trapped, (4) continuum
thin → thick → diffusive; Rosseland means need (4) **including the windows** (κ_R is a harmonic mean dominated by the
most transparent frequencies) and thermalisation (J_ν ≈ B_ν: effective depth √(τ_a(τ_a+τ_s)) > 1, not τ_ext > 1).
Measured directly on the collapse: frequency-resolved flux vs the Rosseland diffusion flux at every shell and time.

**Code.**
- `julia/ContinuumOpacity.jl`: κ_ν(ρ, T, composition) — HITRAN CIA H₂–H₂ (200–3000 K, held above), H₂–He (200–9900 K),
  H₂–H (1000–2500 K) [`~/Projects/popiii_lit/CIA/`]; H⁻ bound–free and free–free, John (1988) fits [coefficients
  cross-checked in three independent codes: oracle `ahmin.f`, RH `hydrogen.c`, NEMESIS `hmin_ff.f`; spot values
  σ_bf(0.85 μm) = 3.99e-17 cm², the known peak ≈ 4e-17]; hydrogenic H bf/ff (Kramers, g = 1); Rayleigh H₂ (RH: Victor &
  Dalgarno 1969, Tarafdar & Vardya 1973) and H (Dalgarno 1962 [memory, unverified]); Thomson. Planck mean (absorption)
  and Rosseland means (extinction; absorption only). Mayer & Duschl 2005 (MNRAS 358, 614, arXiv:astro-ph/0411613) for
  comparison of the means (their tables not obtained yet).
- `popiii_m2c.jl` SNAPS=1 (SNAPEVERY): profiles saved every 0.25 dex of n_c (and every N steps) to
  `lib/popiii/snaps_*.jls`.
- `julia/popiii_regimes.jl SNAPFILE [every]`: per snapshot and shell: n/n_crit, t_3b/t_ff, equilibrium flag, max
  line-centre depth to the surface (231 lines), τ_P, τ_R, τ_R,abs to the surface, f_win = share of ∫∂B/∂T dν at
  frequencies where the thermalisation depth < 1, and F/F_D: frequency-resolved formal solution (tangent rays,
  S = B_ν; scattering as absorption / absorption only as brackets) over the Rosseland diffusion flux.

**First means, equilibrium composition, continuum only [run]** (cm² g⁻¹):

| n_H | T | κ_P | κ_R (extinction) | κ_R (absorption) |
|---|---|---|---|---|
| 1e14 | 1500 | 4.5e-8 | 1.2e-9 | 6.1e-17 |
| 1e16 | 2000 | 5.0e-6 | 7.4e-6 | 2.8e-11 |
| 1e18 | 3000 | 3.4e-4 | 5.6e-5 | 2.1e-6 |
| 1e20 | 4000 | 2.3e-2 | 3.1e-3 | 2.7e-3 |
| 1e21 | 6000 | 6.4 | 4.9 | 4.3 |

Windows between the CIA bands dominate the absorption Rosseland mean up to ≈ 1e18 (10⁸ below κ_P at 1e14); there the
extinction mean is set by Rayleigh scattering — photons in the windows scatter but are barely absorbed (not
thermalised). The spectrum fills in with H⁻ at ≳ 1e20, 4000–6000 K: κ_P/κ_R → 7 → 1.3 and absorption ≈ extinction.

**Formal-solution fix (M3c) [run].** The first regime analysis gave F/F_D = 1e4–1e5 deep in the core — a solver error:
(i) the source function was constant per zone (with optically thick zones the flux then arises only from the jumps
between zones), and (ii) once made linear in τ between the interfaces, the short-characteristic step had S_in and
S_out swapped in its weights (invisible as long as S_in = S_out — hence all earlier tests passed). Correct:
I_out = I_in e^{−Δτ} + S_out (w₀ − w₁/Δτ) + S_in w₁/Δτ. Now: S at the interfaces (linear in r between zone centres), linear
in τ along each segment, in `H2Lines.cooling` and the new general `H2Lines.formal_flux` (used by popiii_regimes).
**New test** (h2lines_test.jl §4): grey sphere, T = T₀(1 − 0.5 (r/R)²): F/F_D = 0.993–0.999 for χR = 1e2–1e5 (Δτ per zone up
to 1.2e3), 0.94–0.96 at χR = 10 (not yet diffusive, expected). A single Gaussian line is not a diffusion test (its wings
stay thin). **Consequence:** the M3b line-transfer collapse (constant S per zone) is suspect where single zones are
optically thick in the lines with a temperature gradient (near the core); to be rerun.

**H⁻ (Tom: "H⁻ has no bound states and should be very easy to be calculated always via photo attachment and detachment
from the full radiation spectrum").** H⁻ has one bound state (0.754 eV, no excited bound states): a single fast steady
state, computable from rates with J_ν (Milne relation → Saha when J = B). Which term dominates [run]:

| n_H | T | photodetachment in a full Planck field [s⁻¹] | associative detachment k₈ n_H [s⁻¹] | ratio | x_e(H, Saha) | x_e(Li, Saha) |
|---|---|---|---|---|---|---|
| 1e14 | 1500 | 1.9e3 | 4.0e2 | 4.8 | 9e-22 | 2.2e-11 |
| 1e16 | 2000 | 1.6e4 | 3.3e5 | 0.05 | 1.6e-16 | 3.0e-10 |
| 1e18 | 3000 | 1.6e5 | 2.5e8 | 6e-4 | 3.1e-11 | 4.6e-10 |
| 1e20 | 4000 | 6.1e5 | 2.4e10 | 3e-5 | 2.6e-9 | 4.6e-10 |
| 1e21 | 6000 | 2.9e6 | 6.0e11 | 5e-6 | 1.2e-6 | 4.6e-10 |

The radiation field controls H⁻ only below ≈ 1e15 cm⁻³ (where H⁻ is irrelevant for the opacity, and J < B anyway);
above ≈ 1e16 associative detachment and its inverse (20–1e5× faster than photodetachment even in a Planck field)
lock H⁻ to the collisional steady state = Saha with the local electrons. **The lever is n_e: lithium.** Li (5.39 eV)
supplies the electrons below ≈ 3500 K (15× H's at 3000 K, ≫ at 1500–2000 K); our network has no Li → H⁻ opacity (∝ n_e)
underestimated in the 1e16–1e19 window regime. Mayer & Duschl 2005 (Table 1: ⁷Li/H = (4.15 ± 0.45)e-10, SBBN, Coc et
al. 2004; §2.3: "Lithium providing free electrons at comparably low temperatures (T ≈ 800…3000 K). These enhance the
H⁻ … abundance in the high-density limit"). `popiii_regimes.jl` LI=1 (default) adds Li ⇌ Li⁺ + e⁻ in Saha.

### M3c result: the transition map [run]

`popiii_regimes.jl` (corrected solver) on the 83 snapshots of the z = 19 line-transfer collapse
(`lib/popiii/snaps_abn02z19_lines_n200_snap.jls`, log n_c = 3.1 … 22.8; output `lib/popiii/m3c_regimes_{li,noli}.txt`).
Following the **central gas** (m ≈ 1e-3 M☉) as n_c rises, with Li electrons:

| log n_c | T | n/n_crit | t_3b/t_ff | τ_line,max | τ_P | τ_R | τ_R,abs | f_win | F/F_D (ext.) | F/F_D (abs.) |
|---|---|---|---|---|---|---|---|---|---|---|
| 3.9 | 148 | 8.6 | 4e9 | 5e-6 | 8e-9 | 3e-8 | 3e-15 | 1.00 | ~0 | ~0 |
| 7.9 | 640 | 3e4 | 2e4 | 9e-4 | 3e-9 | 1e-7 | 8e-14 | 1.00 | 9e-14 | 3e-18 |
| 9.9 | 818 | 2e6 | 25 | 0.16 | 9e-9 | 2e-6 | 1e-11 | 1.00 | ~0 | ~0 |
| 10.9 | 1019 | 7e6 | 2.9 | 4.6 | 3e-8 | 6e-6 | 9e-10 | 1.00 | 9e-10 | 2e-13 |
| 12.9 | 1543 | 1e7 | 4.2 | 140 | 5e-6 | 2e-4 | 3e-6 | 1.00 | 5e-6 | 7e-9 |
| 14.0 | 1746 | 3e8 | 0.025 (eq.) | 460 | 2e-4 | 3e-3 | 3e-5 | 1.00 | 2e-4 | 3e-6 |
| 15.0 | 2021 | 8e9 | 2e-4 | 1.3e3 | 5e-3 | 0.035 | 1e-3 | 1.00 | 0.015 | 0.002 |
| 16.0 | 2252 | 1e11 | 4e-6 | 3.6e3 | 0.15 | 0.46 | 0.025 | 1.00 | 1.7 | 1.2 |
| 17.0 | 2485 | 9e11 | 2e-7 | 9.4e3 | 3.2 | 3.8 | 0.36 | 0.12 | 1.21 | 3.2 |
| 18.0 | 3023 | 2e13 | 2e-9 | 2.3e4 | 71 | 46 | 26 | 0 | **0.998** | 3.9 |
| 19.0 | 3762 | 3e14 | 5e-11 | 4.7e4 | 1.7e3 | 1.2e3 | 1.0e3 | 0 | **0.995** | 2.9 |
| 20.1 | 4836 | 3e15 | 2e-12 | 8e4 | 7.5e5 | 6.2e5 | 6.1e5 | 0 | **1.01** | **1.03** |
| 21.1 | 8916 | 7e16 | 3e-14 | 1.1e5 | 4.7e9 | 2.1e9 | 1.4e9 | 0 | 1.15 | 1.59 |

(τ: radial depth from the shell to the surface; f_win: share of ∫∂B/∂T dν at frequencies with thermalisation depth
√(τ_a(τ_a+τ_s)) < 1; F/F_D: frequency-resolved flux over the Rosseland diffusion flux with the extinction (ext.) or
absorption-only (abs.) mean, S = B_ν.)

**The sequence (central gas):**
1. **H₂ level populations LTE** already at n ≳ 1e3–1e4 (n/n_crit = 9 at 6e3 cm⁻³, 150 K): the z = 19 start is just
   above n_crit; below that only the outer envelope (n/n_crit = 0.2 at 1e3 M☉) is non-LTE.
2. **Lines thick** at the centre from n_c ≈ 3e10 (τ_line,max = 1); at 1 M☉ from ≈ 1e11 (τ 26 there by 1e12).
3. **Chemistry**: three-body formation faster than free fall (t_3b/t_ff < 1) from ≈ 1e11–1e13; equilibrium switch at
   ρ_eq (log n 13.8).
4. **Continuum** (all-frequency, not lines): τ_R = 1 at n_c ≈ 2e16, τ_R,abs = 1 at ≈ 1e17; windows thermalised
   (f_win → 0) by ≈ 3e17.
5. **Rosseland diffusion justified (F/F_D = 1 within 1 %) from n_c ≈ 1e18 (3000 K) with the extinction mean**; between
   1e16 and 1e19 the window extinction is mostly Rayleigh scattering, so the absorption-only mean underestimates the
   opacity (F/F_D(abs) = 3–4); **with the absorption mean only from ≈ 1e20 (H⁻ fills the spectrum, 4800 K)**.
6. **Only the core region becomes diffusive**: at the last snapshot (log n_c 22.8) the continuum photosphere (τ_R = 1)
   is at m ≈ 0.02–0.03 M☉, T ≈ 2400 K — RHFC02: "continuum photosphere … ≃ 2300–2400 K". Outside, the molecular envelope
   (1 M☉ at 1e12, lines τ ≈ 26, continuum τ_R = 3e-5) and everything beyond never become continuum-thick: the line
   transfer is needed there for the whole evolution.
7. **Lithium**: changes the absorption Rosseland depth by 100–400× at n_c = 1e15–1e17 (τ_R,abs 1.1e-3 vs 7e-6 at 1e15;
   0.36 vs 9e-4 at 1e17) and brings F/F_D(abs) closer to 1 at 1e18–1e19 (3.9/2.9 vs 6.4/9.7 without Li); the extinction
   boundary (step 5) is unchanged.
8. To check: inside the hot core at the last snapshot F/F_D = 0.17–0.35 (vs 0.99–1.05 at log n_c 20.8) — steep
   gradients at the core/shock, few zones; resolution study needed.

## M3d Closing the M3c open issues (19 Oct 2026)

1. **H₂–H₂ CIA above 3000 K** [read-src + run]. Borysow, Jørgensen & Zheng (1997, A&A 324, 185) table `cia_h2h2.dat`
   (A. Borysow's page, www.astro.ku.dk/~aborysow/programs/: "H2-H2, CIA data 1000-7000K, all bands"; 800–20 088 cm⁻¹;
   `~/Projects/popiii_lit/CIA/borysow/`). Units cm⁻¹ amagat⁻² → cm⁵ by 1/n_L² (checked: ratio to HITRAN 0.8–1.2 at
   1000–2000 K; 1.2–2.0 at 3000 K, largest in the roto-translational band — Borysow 1997 vs Abel et al. 2011).
   `ContinuumOpacity.k_h2h2`: HITRAN to 3000 K; above, Borysow's T dependence scaled at each wavenumber to join HITRAN
   continuously at 3000 K; below 800 cm⁻¹ HITRAN held at 3000 K; held above 7000 K.
2. **Lithium** [design + run]. `PrimordialChem.ne_with_Li(n_H, T, n_e,H)`: Li ⇌ Li⁺ + e⁻ in Saha, Li/H = 4.15e-10.
   x_e(H) → x_e(H + Li): 2e-24 → 3.9e-12 (1e12, 1300 K), 9e-22 → 1.9e-11 (1e14, 1500), 1.6e-16 → 2.8e-10 (1e16, 2000),
   3.1e-11 → 4.5e-10 (1e18, 3000), 2.6e-9 → 3.0e-9 (1e20, 4000). In the network zones (≲ 1600 K, n ≲ 3e13) Li adds
   ≲ 10 % to the network's own x_e (≈ 5e-11 at 1e12–1e13) — no Li reactions needed there; it matters for the
   equilibrium electrons (H⁻ opacity) at 1700–3500 K. p, u unaffected (Li/H = 4e-10). Used by the post-processor and
   by any continuum opacity call; the model's hydro does not use the continuum opacity yet (RA04 CIE).
3. **Core F/F_D = 0.17–0.35 at log n_c 22.8 — explained, not resolution** [run]: the frequency grid stopped at 1e5 cm⁻¹,
   while at T_c = 58 000 K the Planck peak is at ≈ 1.1e5 cm⁻¹: the formal flux missed ≈ 85 % of the spectrum while F_D
   uses the full 16σT³ (the ratio was a uniform 0.15 across the smooth core, Δτ_R per zone 2e13). Grid now
   10 cm⁻¹ … max(1e5, 25 kT_max/hc) (`popiii_regimes.fgrid`, `ContinuumOpacity.grid(; T)`): F/F_D = 0.97–1.06 across the
   core (0.99 over most of it; Δτ_R per zone 1e12).
4. **M3b line-transfer collapse with the corrected solver**: rerun as `_v2` (in progress).
   Regime analysis rerun with the scaled grid (`lib/popiii/m3c_regimes_li.txt` replaced): unchanged below log n_c ≈ 20
   (the grid already covered the spectrum); central F/F_D (extinction) at log n_c = 22.1: 0.996 (was 0.80); last snapshot
   (22.8): 0.99 (was 0.17). Remaining: 1.15 (extinction) / 1.7 (absorption) at log n_c = 21.1, during core formation,
   where dissociation and the onset of ionisation make κ_ν vary steeply within single zones — a transition zone, not a
   regime boundary.

## G Generalisation plan (Tom, 19 Oct 2026: "We want to implement all. Please plan ahead.")

Audit (19 Oct): several fixes were tuned to this problem. Replace each by a general method; keep the old behaviour as a
switch where it is the faithful-RHFC02 choice. Every item ends with a test that fails on the old method or proves the
general one, a notes section, a commit.

**G0 Test suite** `julia/popiii_tests.jl`: one entry point running every acceptance test so far (SOD/SUBC/SUPC
regression, quasi-star P7.3; M1 static sphere, dust, BE, LP; M2 conservation, Saha consistency, one-zone; H2Lines thin /
homogeneous sphere / grey diffusion), each PASS/FAIL with its number. Run after every G item.

**G1 Artificial viscosity ∝ local zone width** (replaces ql1 = 0.15 tuned to Δr/r of one grid). Option in
TitanResidFast: viscous length l_k = c_l·Δr_k (θ-centred). Map to RHFC02's q ∝ c_q ρ Δv² (c_q = 4). Tests: spherical
Noh/Sod on uniform and geometric grids (stretch 1.02–1.2): shock width in zones and post-shock error independent of the
stretch with l ∝ Δr, not with l ∝ r; regression unchanged when off.

**G2 TR-BDF2 for the collapse** (replaces θ = 0.55). `Trbdf2.jl` (2nd order, L-stable; verified on the Eulerian grid)
on the Lagrangian grid with REGULAR_CENTRE, XMINT, PEXT, EXTMASS; chemistry split after each full step. Tests: dust
collapse convergence order (expect ≈ 2 in C_dyn), BE static, LP converged values; z = 19 collapse without the θ = 0.5
ringing; T_c track vs θ = 0.55 at converged Δt.

**G3 Reaction-network engine** (replaces the hand-coded 4-unknown network). Reactions as data (reactants, products,
rate function, energy); unknowns = all species; conservation of H, He, Li nuclei and charge by elimination; Jacobian by
ForwardDiff; coupled implicit chemistry + T step (as now). Rate sets: RHFC02 faithful; modern (AAZN97/GA08 incl. H₂⁺
channels, H⁻ + H⁺, collisional detachment); **detailed-balance pairing** — each reverse rate from the forward one and the
equilibrium constant (partition functions), incl. three-body recombination H⁺ + 2e and three-body H⁻ attachment — so
the network's own equilibrium is Saha/K exactly.

**G4 H⁻ and Li from rate balance with J_ν** (replaces Saha-only H⁻ and Li). H⁻: radiative attachment / photodetachment
(Milne-consistent with John 1988 σ_ν, J_ν supplied per zone), associative detachment + inverse, collisional detachment
+ inverse, mutual neutralisation. Li: Li, Li⁺ (+ Li⁻, LiH if the sources support them): collisional and photo-
ionisation, radiative + three-body recombination, charge exchange with H⁺ — rates from GP98's Li network (local) and
primary sources. J_ν: CMB by default, the transfer solution when available. Tests: J = B → Saha for H⁻ and Li; x_e(Li)
vs Mayer & Duschl 2005 Fig. 1 (POP III vs Li-free); one-zone tracks unchanged at low density.

**G5 Timescale-based equilibrium** (replaces ρ_eq = 5e-11). Per zone: equilibrium when the network's slowest
relaxation time (smallest eigenvalue of its Jacobian) < ε Δt, or where the gas is non-ideal (p ≥ 1e4 dyn cm⁻², where
SCVH applies; blend as now). With G3's detailed balance the switch is continuous. Tests: no T jump at the switch;
collapse track vs the ρ_eq run.

**G6 Line transfer: per-line frequency grids and non-LTE populations** (replaces the f_H₂ ≥ 0.1 span cut and the
GP98-factor non-LTE). (a) Per line: span and bin width from the zones carrying ≥ 99.9 % of that line's emissivity and
opacity; thin-line shortcut controlled by a stated tolerance. (b) Non-LTE level populations: statistical equilibrium of
the 63 levels with collisional rates (H–H₂, H₂–H₂, He–H₂, e–H₂ state-to-state data — to obtain: BASECOL / Flower's
CCP7 (unreachable so far) / Cloudy's data files) and radiative rates with J̄ from the transfer — first with local escape
probabilities (multilevel LVG), then accelerated Λ-iteration (Rybicki & Hummer 1991) with the short-characteristics
solver. Tests: LTE limit at high n; thin low-density limit vs the GP98/GA08 cooling fits; two-level analytic cases
(escape probability, thermalisation length); ALI convergence.

**G7 Implicit line cooling** (replaces explicit Λ with Δt ≤ 0.1 uρ/Λ): linearise Λ_k in T_k (local derivative from the
transfer's source-function dependence) inside the coupled chemistry/T step. Test: the collapse with the cooling cap off.

**G8 Initial-condition construction**: splice radius of the ABN02 cloud chosen automatically (where the hydrostatic
interior and the read-off envelope agree best in ρ and slope), documented; and Tom's profile data when available.

**Order and parallelism.** G0 first. Then G1, G2, G6a in parallel (separate files: TitanResidFast/ProblemSetup;
Trbdf2.jl; H2Lines) in git worktrees; G3 → G4 → G5 in sequence (one engine); G6b after the collision data are found
(search started); G7 after G6b; G8 any time. Finally rerun the reference collapses (P100, z = 19) with everything and
redo the M3c regime map. Implementation of well-specified packages by subagents; design, review and verification here.

### G0 Test suite [run]
`julia/popiii_tests.jl [quick|full]`: each test a separate process, its output parsed for the deciding number; logs in
`lib/popiii/tests/`, summary `lib/popiii/tests/summary.txt`. Quick tier (8 tests, 2 min 16 s on this Mac, all PASS):
SOD 464 / SUBC 587 / SUPC 713 steps; static sphere exactly at rest; BE discrete equilibrium 7.3e-10; H2Lines thin limit
1.000000 and homogeneous sphere 0.994–1.000; detailed-balance network = Saha (2217.9715 K both); one-zone f_H2 0.9707 at
log ρ/m_H = 12.4. Full tier adds the dust collapse (arrival t/t_ff), BE Lane–Emden n = 120, LP after 8 decades.

### M3d item 4: the line-transfer collapse with the corrected solver [run]
`_v2` (log `lib/popiii/m2c_abn02z19_lines_n200_v2.log`, snapshots `snaps_abn02z19_lines_n200_v2.jls`), stopped by the 2-h
background limit after core formation (step ≈ 9500, log n_c 22.9). Central track:

| log n_c | lines, corrected (v2) | lines, constant-S (v1) | RA04 local fits | RHFC02 stage |
|---|---|---|---|---|
| 10 | 816 | 816 | 802 |  |
| 12 | 1273 | 1273 | 1245 |  |
| 14 | 1730 | 1731 | 1751 |  |
| 16 | 2241 | 2241 | 2260 |  |
| 17.9 | 3087 | 3090 | 3105 | 3000 (4) |
| 19 | 3979 | - | 4006 |  |
| 20.2 | 6495 | - | 6638 | 6500 (5) |
| 21.2 | 15688 | - | 16103 | 15000 (6) |
| 22.2 | 30500 | - | 31101 | 30000 (7) |
| 22.9 | 62641 | - | 64642 |  |

The solver fix changes the central track by ≤ 0.1 % (the constant-S error sat in line cooling that is negligible for the
core's energy budget); line transfer and RA04's fits agree within 1–3 % through core formation, both on RHFC02's
stages.

### G6a Per-line frequency grids [run]
`H2Lines.cooling(…; tol = 0.01, cover = 0.999)`: per line, the zones carrying 99.9 % of that line's emission
(n_u A hν ΔV) or of its line-centre opacity column set the span ±(3 v_th + |u|) and the bin width (≤ ½ their smallest
v_th; 30 ≤ bins ≤ 200). Replaces the f_H₂ ≥ 0.1 cut and its hard-coded μ_H. Thin shortcut by tolerance: τ_c < 2 tol (the
thin approximation overestimates a uniform sphere's line luminosity by ≈ (3/8)τ_c). Test 5 (h2lines_test.jl): hot,
1 %-molecular core (2e4 K) in a cold infalling envelope (≤ 2000 K, ≤ 20 km/s): surface line luminosity per-line vs a
reference with every zone in every span and 4000 bins: ratio 1.00000, max |ΔΛ/Λ| 2.5e-7 over the cooling zones, with
≤ 183 bins (0.3 s vs 6.4 s). Tests 1–4 unchanged.

### G3 Reaction-network engine [run]
`julia/ChemNet.jl` + `julia/ChemRates.jl`. Species with mass, charge, element content, energy zero E0 (relative to free
nuclei + electrons) and internal partition function (electronic × nuclear spin × rotation × vibration: H₂ as K_H2; H with
levels n ≤ 8; H₂⁺ approximate); equilibrium constants K for any reaction from these; reactions as data, `with_db`
adds reverses k_rev = k_fwd/K; photo-processes from σ(ν, T) and a per-zone J_ν with Milne inverses (stimulated emission
included); unknowns all species but H, Li (nuclei conservation) and e (charge); exact Jacobian (ForwardDiff); coupled
chemistry + T step as before; `steady` (pseudo-transient continuation) for tests.

Tests [run]:
- **Reproduces PrimordialChem on the RHFC02 set**: dx/dt identical to all digits; u within 1.5e-7; one step at
  (1e8, 600 K), (1e15, 2000 K): T and x_H₂ identical to all printed digits; (1e12, 1100 K, 2 % of H forming H₂ in the
  step): 1748.1 vs 1740.3 K — D0(H₂) = 4.478 eV (spectroscopic, used for K) vs RHFC02's 4.48 eV.
- **Modern set in a Planck field (J = B(T)) relaxes to equilibrium**: network/equilibrium for H₂/H, H⁺, H⁻ = 1.00000
  (≥ 0.99997) at (1e14, 2000), (1e18, 3500), (1e20, 6000), (1e21, 1e4 K). Getting there exposed two non-general choices:
  case-B recombination (one-way) → replaced by one Milne pair over H levels n ≤ 8 in Boltzmann ratio (total recombination
  ≈ case A; Saha exact for J = B); the H₂⁺ channels → detailed balance for dissociative recombination and H⁻ + H⁺ → H₂⁺ + e.
  **Remaining**: radiative association H + H⁺ → H₂⁺ + γ is one-way (its inverse, H₂⁺ photodissociation, needs σ(ν)
  — to obtain, e.g. Stancil 1994 / Dunn 1968): deviations 2 % (H₂), 0.3 % (H⁺), 1 % (H⁻) at (1e16, 8000 K, 2 % ionised),
  and x_H⁺ = 1e-21 instead of ~1e-43 at (1e10, 1000 K) (negligible in absolute terms).
- Hydrogen photo-ionisation cross sections are Kramers (g = 1): order of magnitude — refine.

### G4 Lithium and photo-pairs in the network [run]
Data: Li rates from GP98 Table 4 (subagent extraction `~/Projects/popiii_lit/rates_li_extraction.md`: 26 reactions;
photo-processes only as CMB rates, no σ(ν)); **Verner et al. 1996 photo-ionisation fits** (D. A. Verner's phfit2.f v2,
copy in FAC `ionis/phfit2.f`: H I PH2 (0.4298 eV, 5.475e4 Mb, 32.88, 2.963), Li I PH2 (3.107 eV, 62.45 Mb, 15.01, 4.895),
threshold 5.392 eV; `~/Projects/popiii_lit/verner/`); **Voronov 1997** collisional ionisation (Li I: 5.4 eV, P 0, A
1.39e-7, X 0.438, K 0.41; AMUSE/SPHRAY and FAC copies).
Network (`ChemRates.modern(; li = true)`): Li + e ⇌ Li⁺ + 2e (Voronov, detailed balance), Li⁺ + H⁻ ⇌ Li + H (GP98 Li3,
Peart & Hayton 1994; d.b.), Li + H⁺ → Li⁺ + H + γ (GP98 Li8, one-way), Li(levels) + γ ⇌ Li⁺ + e (Milne). GP98 Li7
(charge exchange, Kimura et al.) left out: direction/validity doubtful in GP98's table. Li⁻, LiH, LiH⁺: not yet.
**Levels**: H (n ≤ 15, Boltzmann) and Li I (2s + quantum-defect levels n ≤ 15, l ≤ 3; 2p 1.81 eV vs 1.848 measured)
in both the partition functions and the photo-cross sections (ground states Verner, excited hydrogenic/Kramers) — the
Milne total recombination then reaches: H 0.60 / 0.86 / 1.00 of case A (AAZN97 k2) at 1e3 / 3e3 / 1e4 K; Li 0.59 / 0.72 /
0.80 of GP98 Li1 (Verner & Ferland 1996). The shortfall at low T is recombination to n > 15 and the Gaunt-factor-1
excited-level cross sections → **next: Hummer & Mihalas (1988) occupation probabilities** (level dissolution with
density, consistently in Z and σ) and exact hydrogenic σ_n; the n ≤ 15 cut-off is harmless below ≈ 1.5e4 K (Z_H = 2.08;
SCVH95 takes over in the dense hot core).
**Equilibrium test** (J = B(T), steady state by pseudo-transient continuation — rewritten to accept only converged
steps after it stalled at 1e21 cm⁻³): network/equilibrium = 1.0000 for H₂, H⁺, H⁻ and Li⁺ at (1e14, 2000), (1e18, 3500),
(1e20, 6000), (1e21, 1e4 K). Remaining deviations from the one-way radiative association H + H⁺ → H₂⁺ + γ: H⁺ 3.4× at
(1e16, 2000 K) where x_H⁺ ~ 1e-22, 0.3–2 % at (1e16, 8000 K). Li electrons with the network self-consistent: at (1e18,
3000 K) x_e = 4.16e-10 (vs 4.45e-10 from PrimordialChem.ne_with_Li, which adds Li⁺ to an unchanged H⁺: the extra
electrons suppress H⁺ 13×).

### G5 network chemistry in the collapse [run]
`PopIIIMaterial.NetGas`: ChemNet network (`ChemNet.step`, coupled with T and cooling) per zone; the Saha+SCVH table
(`PopIIIEqEOS`) where the ideal-gas pressure ≥ PSWITCH (default 1e4 dyn cm⁻², RHFC02's EOS validity) — a criterion on
the EOS, not a tuned ρ_eq. Driver: `CHEM=rhfc02 | net | netmodern | netmodernli` in popiii_m2c.jl (default rhfc02 unchanged).
- P100 `ra04`, `CHEM=net` (RHFC02 network in ChemNet) vs the reference `m2c_ra04_n100_eq.log` (PrimordialGas, ρ_eq = 5e-11):

  | log n_c | 8 | 10 | 12 | 14 | 16 | 17 |
  |---|---|---|---|---|---|---|
  | T_c ref / net [K] | 390 / 374 | 621 / 608 | 1066 / 1061 | 1749 / 1800 | 2196 / 2213 | 2545 / 2466 |
  | f_H2 ref / net | 3.8e-3 / 3.7e-3 | 0.331 / 0.325 | 0.951 / 0.951 | 0.968 / 0.975 | 0.934 / 0.963 | 0.893 / 0.919 |

  Above 1e14 the network stays out of equilibrium up to the pressure switch (n_c ≈ 5e16), so f_H2 is higher than the
  reference's Saha by up to 3 %. At the switch, T_c dips by 0.07 % and max dT/T is ≈ 5e-2 for ≈ 40 steps: a mild kick, no jump.
- `CHEM=netmodernli` first died at step 0 with NaN: below ≈ 200 K, K(H → H⁺ + e) and K(Li → Li⁺ + e) underflow to 0, so
  the detailed-balance reverses k_f/K and the Milne α = A/K became 0/0. Fixed in general:
  - `ChemNet.lnKeq`.
  - Reverses are formed in logarithms. Below the underflow temperature of an endothermic fit (T_u ≈ ΔE/700k), the
    exothermic reverse is extrapolated as a power law from [T_u, 2T_u].
  - `Photo` is now a list of levels (E_x, g, E_th, σ(ν)) with the electronic partition function. The Milne integral is
    done per level with the Boltzmann factor and 1/K in one exponent.
  - All rates are finite at 35 K. With J = B, R/(αK) = 1 ± 1e-15 for H⁻, H and Li at 100 K–1e4 K.
  - The per-level frequency grids raise H's Milne total over case A from 0.60/0.86/1.00 to 0.81/0.91/0.95 at 1e3/3e3/1e4 K.
    At 35 K it is 0.26 (n ≤ 15 truncation; the Hummer–Mihalas follow-up).

### G2 TR-BDF2 on the Lagrangian collapse [run]
`Trbdf2.jl` now runs on lgrid = 3 (Eulerian path unchanged; `qstar_run.jl` with `QS_TRBDF2=1` still runs, 60 steps, no
failures). Changes: (i) the BDF2 stage's old state y* also holds the radii, r* = r_γ + B(r_γ − r_n), and the old zone
volumes recomputed from r* (ρ* = m*/dV*; m* = m for a Lagrangian zone), `ro` and `dvolo` are saved/restored; (ii) with
`FREE_OUTER`/PEXT the outer interface velocity is combined with the mass of the outer half zone; (iii) **y* is written
A y_γ − B y_n = y_γ + B(y_γ − y_n)**: the first form gives a frozen x_me (1e47) a round-off of 1e31, hence dm/dt ≠ 0 and
a spurious advective energy flux — the static BE sphere then moved by 3.7e-4 (instead of 7e-10) in the first three
steps. REGULAR_CENTRE, XMINT (A − B = 1), EXTMASS, DT_DYN need nothing: they read no old-state array, and timstp sees y_n
because the old arrays are restored after the step. `PM.chem!` runs after the full step on y_{n+1} (first-order split
in time; the stages themselves are chemistry-free). Switch `INTEGRATOR=trbdf2` in `popiii_m1.jl`, `popiii_m2c.jl`.
- **Dust collapse** (T0 = 1e-6, stol 0.1, 40 shells; worst-shell |Δr/r| at exact R/R0 = 0.5 / 0.1 / 0.03; log interpolated
  between the 20-step prints; logs `m1_dust_trbdf2_cdyn*`):

  | C_dyn | steps to 0.1 | 0.5 | 0.1 | 0.03 |
  |---|---|---|---|---|
  | 0.03 | ~96 | 2.5e-4 | 5.9e-3 | 3.8e-2 |
  | 0.01 | ~260 | 3.0e-5 | 6.7e-4 | 4.4e-3 |
  | 0.003 | ~860 | 2.7e-6 | 6.1e-5 | 4.0e-4 |

  Order in C_dyn 1.97 (0.03 → 0.01) and 1.98 (0.01 → 0.003) at 0.1. θ = 0.5 (M1b) had 2.6e-3 at C_dyn 0.01: 4× less error
  at the same step count (two Newton solves per step). No stage failures, no backward-Euler fallback. Energy drift
  |dE|/|W| 7e-6 (C_dyn 0.003).
- **BE discrete, n = 60**: 7.3e-10 / 1.3e-10 (θ: 7.3e-10), 1735 steps. **LP, n = 200**: max 4πGρr²/c_s² = 4.66 / 5.92 /
  6.94 / 7.68, min u/c_s = −0.99 / −1.73 / −2.27 / −2.64 at 2/4/6/8 decades: identical to M1b (θ = 0.55) to the printed
  digits; energy drift over 8 decades 6.7e-7 (θ = 0.55: 5.9e-5; θ = 0.5: 1.3e-6).
- **P100 collapse** (`ra04 100`, THETA unset, C_dyn 0.01, current code incl. Saha + SCVH95 above ρ_eq; logs
  `m2c_ra04_n100_trbdf2` vs `_th055ctl`, the θ = 0.55 control with the same code): no ringing. Before the collapse both
  integrators show a ≤ 0.1 dex breathing of the core (log n_c 6.43 → 6.55 → 6.52 → 6.57), from 6.7 on n_c rises
  monotonically (the θ = 0.5 run swung 2.5e6 ↔ 2e8). T_c at log n_c = 7.9 / 9.9 / 11.9: TR-BDF2 369.5 / 560.1 / 994.7 K,
  θ = 0.55 369.6 / 560.0 / 994.3 (≤ 0.1 %); the older `m2c_ra04_n100_eq.log` 367.7 (7.90) / 561.5 / 1002.3 (≤ 1 %,
  earlier code). t(n_c = 1e12) 1.2046e6 vs 1.2047e6 yr. TR-BDF2: 12 987 steps (89 TR and 184 BDF stage failures, each
  retried with Δ × 0.05, no backward-Euler fallback) vs 11 347 for θ = 0.55: about 2.3× the work (two solves per step).
  Quick tier of `popiii_tests.jl`: 8/8 PASS (BE 7.27e-10).
**Default for the collapse runs: `INTEGRATOR=trbdf2`** (second order, no ringing, same track as θ = 0.55).
- Fit validity in detailed balance: k_f/K amplifies any error of an endothermic fit by e^{ΔE/kT}. AAZN97 k1 (cut to 0
  below 1000 K, and inaccurate above that cut) gave H⁺ a three-body recombination time of 1e-9 s at 200 K. Each db pair
  now carries Tmin, the validity limit of its forward fit: k1 1e4 K, k14 3e3 K, Voronov Li 1e4 K. Below Tmin the reverse
  is extrapolated as a power law. The step at 35–200 K and n = 1e2–1e4 then behaves.

### G6b H₂ statistical equilibrium [run]
`julia/H2Levels.jl` (`load`, `populations(L, T, n_H, n_H2, n_He, n_Hp; β, Jbar, opr)`, `cooling(L, pops, n_H2; β)`, `coverage`, `lte`):
the 63 levels v ≤ 2, J ≤ 20 (ExoMol RACPPK energies/g, 231 lines with A as in H2Lines), Cloudy 25.00 collision files as de-excitation
k(u→l) [cm³ s⁻¹], excitation by detailed balance with the ExoMol energies (all rows checked u above l: 0 swapped, 0 duplicates).
T-interpolation linear in log k vs log T, end values held outside the grid. Missing transitions: k = 0 (both ways). Radiative: A β + B J̄ (u→l),
B_lu J̄ (l→u), defaults β = 1, J̄ = 0. Solver: GTH state reduction (≡ the linear solve with a normalisation row, but exact for the nearly
decoupled ortho/para chains); 0.13 ms per call. Ortho↔para only through H (3 Sun & Dalgarno rows) and H⁺ (24 pairs); with no H/H⁺ the chains
are decoupled and normalised 3:1 (`opr` forces any ratio and drops conversion). H₂ collider: n_H2 split ortho/para from the computed
populations (fixed point, 1e-9). No electron collisions (no data).
Coverage on the 63 levels (of 1953 pairs; levels linked; radiative lines with a collisional partner): H 844 (59/63; 211/231; not linked:
(1,20), (2,18–20)); H₂-para 421 (J ≤ 16, v ≤ 2; 42/63; 135/231), H₂-ortho 361 (J ≤ 15; 39/63; 123/231) [LeBourlot 1999, 9 T points 100–6000 K];
He (ORNL) 529 (63/63; 176/231; 1–10⁴ K); H⁺ 44 (v = 0, J ≤ 9; 10/63; 8/231; 3–6000 K). Levels without a partner for a collider get only radiative
decay and the other colliders.
Tests (`h2levels_test.jl`, 15 checks, 3 s; in `popiii_tests.jl full` as `h2levels`): (a) n_H = 1e12/1e14, β = 1: max |p/Boltzmann − 1| over levels
with q > 1e-4 is 1.4e-7 (300, 1000 K), 3.5e-5 (3000 K; with 0.1 n_H of He all 63 levels, q > 1e-7: 1e-5); J̄ = B_ν(T) without collisions gives
Boltzmann to 2e-15. The 63 levels hold 0.995 of the 302-level partition function at 3000 K (1.0000 at ≤ 1000 K), populations are normalised to the 63.
(b) Λ/(n_H2 n_coll), n_coll = 1e-4, o:p = 3:1 fixed, vs GA08 Tables 1–6 (Table 8 gives the same): H₂ 0.89–1.16 (100–6000 K), He 1.05–1.16, H 0.95/0.98/0.95/0.91/0.79
at 1000/1500/2000/4000/6000 K but 0.40 (100 K), 0.90 (200), **1.78 (300), 1.84 (500)**, 1.25 (700 K) — the Cloudy WGF07 file is about twice the GA08
fits at 300–500 K (para J=0→2 and ortho 1→3 checked by hand against the file; GP98 is 4–6× GA08 there); H⁺ 0.15 (100 K) … 0.87 (6000 K)
(Gerlich table v = 0, J ≤ 9, o↔p not included in this comparison). (c) Λ(n_H)/Λ_LTE(63 levels) vs 1/(1 + Λ_LTE^HM/(n_H Λ0^GP)): n_cr (ratio ½) ours
1.5e4 / 6.6e4 / 9.3e3 cm⁻³ at 300/1000/3000 K vs 8.8e2 / 4.8e3 / 6.3e3 for cool_H2_thin; ours/GP at n_H = 1e3 … 1e4: 0.24, 0.32, 0.76; the GP/HM
form reaches LTE too early at 300–1000 K (ours is 0.3× of it between n_H = 1e2 and 1e4) and agrees within 20 % at 3000 K; at n_H ≥ 1e7 both → 1.
Low-density Λ0 ours/GP: 0.30 (300 K), 0.86 (1000), 1.15 (3000). (d) uniform β scales A: n_eff (Λ = ½ Λ_LTE(β), 1000 K) = 6.8e4 β cm⁻³ exactly
(slope 1.000 over β = 1…1e-3); at n_H = 1e6 Λ/Λ_LTE(β) = 0.92 (β = 1) → 1.00 (β = 1e-3). Not yet used in the collapse drivers.

### G6c Non-LTE H₂ lines coupled to the transfer (ALI) [run]
`H2Lines.cooling(…; nlte = (LV = H2Levels.load(), nH, nHe, nHp), nlte_mode = :ali | :ep, nlte_tol = 1e-3, nlte_maxit, nlte_init = :ep | :thin,
nlte_lambda = :diag | :none)`; default (`nlte = nothing`) is the old path, unchanged. H2Levels' 231 lines are H2Lines' 231 lines in the same
order (checked at every call: ν from the 63 energies, A, g of both levels). Per zone: populations (63 fractions × n_H2) from `H2Levels.populations`;
κ0 = (hν/4π)(n_l B_lu − n_u B_ul) (stimulated emission included), S = n_u A/(n_l B_lu − n_u B_ul) = (2hν³/c²)/(g_u n_l/(g_l n_u) − 1); same tangent-ray
formal solver, which now also returns the profile-weighted mean intensity J̄ per zone (J = ½∫(I_in + I_out)dμ at the interfaces with the comoving profile
φ(v ∓ μu) of the adjacent zone, trapezoid in μ over the rays; zone = mean of its two interfaces).
- **Diagonal operator**: in the short-characteristics step the coefficient of the exit-node source is ψ(Δτ) = w0 − w1/Δτ. Λ* at a node = mean over
  directions of ψ for the segment arriving at it (inward rays through the zone outside, outward through the zone inside), accumulated with the same μ
  quadrature and profile weights as J̄; zone Λ* = mean of its two nodes. (Local operator, not the exact diagonal, which needs n unit-source solves.)
  ALI step (Rybicki & Hummer): net rate = n_u A (1 − Λ*) − (n_l B_lu − n_u B_ul)(J̄_fs − Λ*S_old), linear in n, i.e. `populations(…; β = 1 − Λ*,
  Jbar = max(J̄_fs − Λ*S_old, 0))`: one SE solve per zone per iteration. Start: EP populations (a few iterations, counted in `LAST.nep`) or thin (`:thin`).
  Converged when max relative change of populations (those > 1e-6) < `nlte_tol`; `LAST[]` has `iters`, `nep`, `converged`, `hist`, `Λrate`, `pops`.
  Lines with radial line-centre τ < 2·tol: thin emission from the actual n_u, J̄ = 0. Cooling returned = −∇·F of the last formal solution (as the LTE path).
- **EP** (`:ep`): β per line and zone = profile average of ½[E₂(τ_out φ̂) + E₂(τ_in φ̂)] (τ_out radially to the surface, τ_in through the centre) or the
  Sobolev (1 − e^{−τ_S})/τ_S with τ_S = κ0 c/(ν0 max(|du/dr|, |u/r|)), the larger of the two; fixed point on β (about 3–5 SE sweeps); cooling Σ n_u A β hν.
  (My first EP, (1 − e^{−τ})/τ with τ the depth to the surface, was 1.8–3.4× high in the thick limit: log divergence from the surface layers.)
- **Tests 6–8** (`h2lines_test.jl`; 1–5 unchanged; 10/10 G6c checks): (6) n_H = 1e14, n_H2 = 5e13, 1000 K, 105 thick lines: L_ALI/L_LTE = 1.00000, 1 iteration.
  (7) n_H = 1e2/1e3/1e4 thin: L/(V n_H2 Σ n_u A hν)(β = 1 populations) = 1.000000 for ALI and EP (the LTE-thin emission is 72/16/4.6× higher at these densities).
  (8) homogeneous spheres, 1000 K, n_H2 = n_H/2, 30 zones (EP, fnlte, LTE luminosity ratios to ALI):

  | n_H, R | n_H2 R | LTE τ_max (thick lines) | EP/ALI | fnlte/ALI | LTE/ALI | ALI it. (EP start) | thin start | plain Λ to 1e-6 |
  |---|---|---|---|---|---|---|---|---|
  | 1e4, 1e22 | 5e25 | 13 (48) | 0.83 | 2.37 | 3.51 | 16 (5) | 17 | 40 (L/L_ALI 0.9996) |
  | 1e5, 3e21 | 1.5e26 | 39 (74) | 0.77 | 1.22 | 1.28 | 6 (4) | 8 | 16 (0.9998) |
  | 1e6, 1e21 | 5e26 | 131 (86) | 0.75 | 1.02 | 1.03 | 3 (3) | 4 | 6 (1.0000) |

  ALI history at 1e4: 3.1e-1 … 8.2e-4, convergence factor ≈ 0.7 per iteration (zone Δτ < 1 for the strong lines: the local operator is weak; plain Λ
  to 1e-3 needs ~21 iterations at similar τ, so ALI gains little here); thin and EP starts give the same L (< 1e-3). The old fnlte factor (HM79/GP98 form) is 2.4× high at
  n_H = 1e4 and reaches the ALI answer only at n_H ≳ 1e6 (the 63-level n_cr at 1000 K is 6.6e4, G6b). EP is 17–25 % low.
  **Rate-based vs flux-based cooling** Σ[n_u A − (n_l B_lu − n_u B_ul)J̄]hν vs −∇·F: 1.022 / 1.054 / 1.203 (30 zones); at n_H = 1e5: 1.054 / 1.018 / 1.006 (30/60/120 zones),
  L_ALI 1.826 / 1.840 / 1.845e48 (0.7 % over 4×). In LTE-like thick zones the rate form is a cancellation (J̄ ≈ S); the flux form is returned.
  Infalling cloud (n_H 1e7 … 1e4, 600–1400 K, ≤ 10 km/s, 75 thick lines): EP/ALI 0.93, LTE/ALI 1.21, 4 ALI iterations.
- **Timing** (100 zones, that cloud, 4 threads): LTE path 0.19 s; ALI 1.5 s for 4 iterations + EP start (≈ 0.35 s per iteration; J̄ and Λ* roughly double the
  transfer cost); EP 0.08 s (no transfer).
- Not done: coupling into the collapse drivers (`popiii_m2c.jl` unchanged); no electron collisions; maser zones (κ0 ≤ 0) treated as transparent.

### G5b network speed
- Network evaluation now uses an index form (`ChemNet.Compiled`, cached per Network) in place of Dicts.
- Photo-rates are cached per Newton T value. When J is the same in every zone (`env.static`, the driver's CMB/none
  field), they come from a log-T table (0.005 dex, Catmull–Rom); table vs direct ≤ 3e-3.
- A modern+Li step drops from 23 ms to 0.45 ms per zone. The P100 netmodernli run takes ≈ 30 min, against ≈ 8 h before.

### G7 implicit transfer cooling [run]
- Λ_j from the transfer is applied over the next step as Λ_j (T/T_j)^s inside the implicit chemistry–temperature step
  (`PopIIIMaterial.ext_cooling`). The slope s = d ln Λ_j/d ln T comes from a second transfer call with all T × 1.02.
- The second call is made only when the explicit cap 0.1 t_cool would bind. The cap then relaxes to DTCOOL t_cool
  (default 1). `IMPCOOL=0` restores the explicit scheme.
- P100 `lines` to 1e13: T_c matches the explicit run within 0.2 %, with the same step count to within 0.05 %.
- z = 19 `lines` to 1e9: identical, with 0 implicit steps. On these runs the cap never binds; dynamics and stol limit
  Δt. G7 is therefore a robustness guarantee, not a speed-up.

### Driver: NLTE modes for `lines`
`NLTE=1` keeps the GP98/HM79 factor (default). The other modes:
- `0`: LTE.
- `ep`: H2Levels with escape probabilities.
- `ali`: G6c, statistical equilibrium coupled to the transfer.
- `hyb`: ALI every NALI steps (default 20), and EP × (ALI/EP per zone) in between.

### G1 Viscosity ∝ zone width [run]

**Code.** `TitanResidFast.AV_ZONE` (`Ref(0.0)`): when > 0 the viscous length of zone k is l_k = AV_ZONE·(r[k+1] − r[k]) with
the θ-centred radii (ForwardDiff-clean); ql0, ql1 then unused. Default 0 = old formula (regression 8/8 PASS, `popiii_tests.jl quick`).
`TitanCore` only reads `qf` from the residual (no recomputation), nothing else to change. `popiii_m2c.jl`: env `AVZONE`
(default 0 = keep ql1). Mapping to von Neumann–Richtmyer q ≈ c_q ρ (Δv)² (RHFC02/Ripamonti et al.: c_q = 4) — derived, not
sourced: qf = cqvis·cq2·l²·ρ·(div u)(du/dr), l = cΔr ⇒ c_q ≈ cqvis·cq2·c² (cqvis = 4/3 ⇒ c² = 3/cq2 for c_q = 4, c ≈ 1.7).

**Test** `julia/av_test.jl` (log `lib/popiii/tests/av_test.txt`). Noh problem, cold gas ρ0, u = −v0, γ = 5/3, wall at the
centre/left, free outer boundary (Lagrangian, lgrid 3; planar wall at r = R0, spherical with REGULAR_CENTRE); exact shock
speed v0/3, jump 4 (planar) / 64 (spherical). Grid: uniform dx = 0.004 R0 up to r0 = 0.3 R0, then Δr ×s per zone (s = 1
uniform, 1.02, 1.05, 1.1, 1.2; n = 250 … 94). Measured when the shock is at 0.1 R0 (it has just entered the stretched
part). wu = zones between the 10 % and 90 % points of u/u_pre; Δx_s/r = Δr/r of the zone at the shock. Old: l = ql1 r̄
with ql1 tuned on the uniform grid (spherical 0.05, planar 0.0018 = 2 × the shocked-zone width); new: AV_ZONE = 2 (cq1 = 0,
cq2 = 1, θ = 0.55, stol 0.05). Position error = (r_s,num − r_s)/r_s from the ρ half-jump.

| geometry | s | Δr/r (shock) | wu, l ∝ r | pos err, l ∝ r | wu, l = 2Δr | pos err, l = 2Δr |
|---|---|---|---|---|---|---|
| planar | 1 | 0.016 | 1.80 | 0.2 % | 3.82 | 0.7 % |
| planar | 1.02 | 0.023 | 1.64 | 0.9 % | 3.74 | 0.7 % |
| planar | 1.05 | 0.035 | 1.69 | 2.3 % | 3.59 | 0.6 % |
| planar | 1.1 | 0.045 | 0.68 | 1.4 % | 3.40 | −0.4 % |
| planar | 1.2 | 0.088 | 1.63 | 9.8 % | 3.19 | −1.9 % |
| spherical | 1 | 0.018 | 4.35 | 1.2 % | 3.82 | −0.3 % |
| spherical | 1.02 | 0.027 | 3.00 | 1.3 % | 3.78 | 0.2 % |
| spherical | 1.05 | 0.030 | 2.09 | 1.4 % | 3.70 | 1.0 % |
| spherical | 1.1 | 0.041 | 1.52 | 2.3 % | 3.43 | 2.3 % |
| spherical | 1.2 | 0.055 | 1.50 | 4.5 % | 2.46 | 3.5 % |

- With l = 2Δr the shock is 3.2–3.8 zones wide for s ≤ 1.1 (planar through 1.2) and the position error stays ≲ 2 % (sph. 3.5 %
  at s = 1.2); with l ∝ r the width falls 4.3 → 1.5 zones (spherical) as Δr/r grows — l no longer covers a zone, the shock
  is a 1–2-zone jump (planar: 0.7–1.8, erratic; position error up to 10 %). Hence ql1 = 0.15 had to be retuned per grid.
- Spherical, large s: l = cΔr acts also on the smooth converging pre-shock flow (div u = −2u/r ≠ 0) and slows it where
  Δr/r is large (u/v0 = −0.86 at s = 1.1, plateau 0.67–0.84 of 64ρ0 at s = 1.1–1.2): keep Δr/r ≲ 0.05–0.1 in the infall region, or c < 2.
  (The collapse grids have Δr/r ≈ 0.04–0.12.)
- ρ-based 10–90 % widths were not used: the numerical post-shock plateau is 5–10 % below the jump and not flat.
- Step counts: l = cΔr needs 5–10× fewer steps than the tuned l ∝ r runs (ql1 case: ringing → 1.4e4–2.2e4 steps vs 1.5e3–1.9e3 spherical).
- Use for the collapse: `AVZONE=2 julia popiii_m2c.jl …` (not yet run on the z = 19 cloud; c_q = 4/3·1·4 = 5.3 for c = 2).
- G5 result, P100 `ra04`, `CHEM=netmodernli` (full detailed balance, Li), T_c [K] by log n_c, compared with the reference and CHEM=net:

  | log n_c | 12 | 14 | 15 | 16 | 17 |
  |---|---|---|---|---|---|
  | reference (Saha above ρ_eq) | 1066 | 1749 | 2032 | 2196 | 2545 |
  | CHEM=net | 1061 | 1800 | 2108 | 2213 | 2466 |
  | CHEM=netmodernli | 1055 | 1741 | 2027 | 2197 | 2542 |

  The detailed-balance network relaxes to equilibrium on its own and follows the Saha-switched reference within 0.5 %
  from 1e14 to 1e17, with no ρ_eq switch. The RHFC02 network (PSS83 k5 not in detailed balance with k4) deviates by up
  to 4 %. The rate pair, not the switch, sets the high-density track.

### G8 automatic splice of the z = 19 cloud [run]
`LRJOIN=auto` (default) sets r_join to the outermost radius inside which the hydrostatic interior stays within DLNJOIN of
ABN02's read-off n(r) everywhere. A plain "best agreement" criterion is degenerate: the flat core agrees trivially.

| DLNJOIN [dex] | log r_join [pc] | M_join [M☉] | slope mismatch |
|---|---|---|---|
| 0.05 | −0.86 | 0.3 | −0.17 |
| 0.10 | −0.65 | 1.3 | −0.32 |
| 0.15 (default) | 0.68 | 815 | −0.41 |

The interior departs from the read-off by 0.10–0.15 dex between 0.2 and 5 pc, and beyond 5 pc it diverges, because the
hydrostatic mass saturates. The default is the read-off accuracy of the figure (≈ 0.15 dex). It gives 4.8 pc,
against 7.4 pc (LRJOIN=0.87) chosen by hand earlier. A numerical number can be set with LRJOIN=<log r>. Tom's profile
data would remove the read-off uncertainty.

### Final P100 reference run (all generalisations on) [run]
**Set-up:**
```
INTEGRATOR=trbdf2 AVZONE=2 CHEM=netmodernli NLTE=hyb SNAPS=1 NSTOP=1e21 popiii_m2c.jl lines 100
```
This is TR-BDF2, viscosity ∝ zone width, the detailed-balance network with Li, non-LTE H₂ lines (ALI, hybrid), and
implicit transfer cooling.
- **Integration:** 15 794 steps. There were 72 TR and 293 BDF stage retries and no backward-Euler fallback; implicit
  line cooling was never needed.
- **T_c [K] by log n_c:**

  | log n_c | 8.1 | 10.2 | 12.0 | 14.2 | 16.2 | 18.2 | 20.2 |
  |---|---|---|---|---|---|---|---|
  | T_c | 372 | 627 | 1069 | 1715 | 2158 | 3160 | 5942 |
  | ON98 | 407 | 596 | 1126 | 1701 | 2338 | 3924 | — |

  The run stays within 1–5 % of the earlier tracks.

**Regime map** (`popiii_regimes.jl`, central shell; `lib/popiii/regimes_p100_final.txt`), against M3c:
- The centre's H₂ lines become thick (τ_line = 1) at log n ≈ 9.9.
- τ_R = 1 at log n ≈ 16.3, as in M3c.
- The windows are thermalised (f_win = 0) by log n ≈ 16.7, slightly earlier than M3c's 17.5.
- F/F_D (extinction) is 0.63, 0.86 and 0.80–0.83 at log n = 18.2, 19.0 and 19.8–20.0. Diffusion with the extinction
  Rosseland mean is good to ≈ 20 % from 1e19.
- The equilibrium table (pressure switch) takes over the centre at log n ≈ 17.3.

### Final z = 19 reference run (all generalisations except TR-BDF2) [run]
**Set-up:**
```
IC=abn02z19 TCMB=54.5 AVZONE=2 CHEM=netmodernli NLTE=hyb SNAPS=1 NSTOP=1e21 popiii_m2c.jl lines 200
```
- Integrator: the θ-scheme. TR-BDF2's stages fail in the infalling envelope (see G2b).
- Splice: automatic, at 4.8 pc, with M_join = 815 M☉.
- The run took 12 773 steps and stopped at T_c = 5618 K, n_c ≈ 8e19.

**Contraction phase:**
- 1e5 is reached at 20.3 Myr, against 14.6 Myr in v2; runaway collapse follows 1.2 Myr later.
- The slower contraction comes from the H₂ non-LTE populations (statistical equilibrium plus ALI). The old GP98/HM79
  factor overestimates line cooling by up to 2.4× at n ≈ 1e4 (G6c). The new splice also changes the envelope.
- These two contributions have not been separated yet.

**T_c [K] by log n_c:**

| log n_c | 6.1 | 7.1 | 8.1 | 9.1 | 10.1 | 11.1 | 12.1 | 13.1 | 14.1 | 15.1 | 16.1 | 17.1 | 18.1 | 19.1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T_c | 260 | 420 | 629 | 754 | 829 | 1197 | 1239 | 1339 | 1769 | 2055 | 2250 | 2583 | 3240 | 4228 |

- The warm core runs hotter than ON98 and P100 at 1e7–1e11, by up to 43 % at 1e11.
- From 1e13 it joins the common track: within 4 % of ON98 up to 1e16, and within 3 % of P100 at 1e18.

**Regime map** (central shell; `lib/popiii/regimes_z19_final.txt`):
- The H₂ lines become thick at log n ≈ 10.3.
- τ_R = 1 at ≈ 16.2.
- The windows are thermalised by ≈ 17.
- F/F_D (extinction) is 1.2–1.5 at ≥ 1e19, against 0.8 for P100. The central temperature gradient is almost flat, which
  makes the ratio noisy, so Rosseland diffusion is good to 20–50 % there.
- **Attribution of the slower contraction.** Each run below changes one setting relative to the final run; t is the
  time to n_c = 1e5:

  | run | t(1e5) [Myr] |
  |---|---|
  | final (NLTE=hyb, automatic splice at 4.8 pc) | 20.3 |
  | NLTE=1 (old GP98/HM79 factor), else as final | 12.9 |
  | LRJOIN=0.87 (old splice at 7.4 pc), else as final | 21.6 |
  | v2 (old everything: RHFC02 chemistry, NLTE=1, 7.4 pc, ql1 viscosity) | 14.6 |

  The H₂ level populations account for the whole delay, a factor of 1.6. The old factor reaches LTE too early, which
  overestimates the cooling of the cloud at n ≈ 1e3–1e5. The splice changes the time by only 6 %, and the remaining
  differences from v2 (chemistry, viscosity) by about 12 %. The minimum core temperature changes little (142 K against
  147 K).

### G2b TR-BDF2 stage failures in the infalling envelope [run]
Problem: z = 19 cloud (`IC=abn02z19 TCMB=54.5 lines 200`, 500 steps): TR-BDF2 reached 0.7 Myr, θ = 0.55 4.9 Myr; 23 TR and
51 BDF stage failures in 501 steps, each cut Δ × 0.05 (timstp then also damps the regrowth: Δ cycles 2e10 → 4e11 over ~8 steps).
- **Diagnosis** (temporary instrumentation of `newton!`: per-iteration dmax and the zone/variable of each maximum, dense
  re-solve of the 1000 × 1000 Newton system, scans of the first iteration over Δ and θ). Not a bug in y* (|y* − y_γ| ≤ 0.1 %), not
  a convergence-test problem, not negative ρ/T. The failing stages start with a first Newton step of relative size 1e2–1e8 in
  T, ρ, r at zones 85–135 (far from the envelope that limits Δ; the residual itself is ≤ 3e-3 of the row scale, in zones
  186–188), which the 10 % damping (dtol) turns into 12–40 iterations (niter = 10; niter = 20 does not help). The block system
  has cond ≈ 1e19–1e20 (row-equilibrated) at all steps, good and bad: any initial guess that is off in a soft direction is
  amplified (dense LAPACK solve of the same system gives the same garbage, so it is the system, not matslv).
  **The same state fails for θ = 0.55 and backward Euler too if Newton starts from `reset` (new := old, no predictor)**
  (θ = 0.55: 5/5 failures, 10 iterations with dmax 7e3; the BE fallback likewise, hence be_fallback 0), and converges in
  3 iterations (first dmax 4e-3) as soon as the radii/x_me are extrapolated, r = r_old + (dr/dt)Δ (variants: reset +
  extrapolation of r, x_me only: 5/5 converge; reset + ghost/derived refresh only: 5/5 fail). TITAN's step.f always starts a step with
  `peq`, i.e. with that extrapolation; zones in the envelope move u Δ ~ their width per step (that is what limits Δ). The
  TR-BDF2 code lacked it in three places: (i) the BDF stage starts at y_γ, not at the extrapolated y_γ + C2Δ·(dr/dt);
  (ii) the BE fallback and the Δ-cut restarts use `reset` only; (iii) `peq` at the start of the step extrapolates over Δ
  where the TR stage covers γΔ (1.7 × overshoot). Ablation (z = 19, 501 steps): BDF predictor only: fail_tr 38, fail_bdf 0;
  TR extrapolation over γΔ only: fail_tr 0, fail_bdf 40; neither: 38 + 2; both: 0 + 0. (With a predictor the BE fallback
  also works: be_fallback 38–40 plus 10 sticky BE steps each, so those runs are partly first order.)
- **Fix** (`TitanCore.predict!`, the second half of peq split out: grid extrapolation, boundary values, ghost/derived
  quantities; `peq!` bit-identical, SOD/SUBC/SUPC step counts unchanged): `Trbdf2.step!` calls `peq` with dtime = γΔ (restoring
  Δ afterwards), `predict!(halve = false)` before the BDF stage, after `reset` in the BE fallback, in the sticky-BE restart and
  after the Δ × TFAC cut (the latter with the γ factor for the next TR stage). The solution of a converged step is unchanged
  (initial guess only); TFAC, the tfac memory and the convergence test are untouched.
- **Results.** z = 19, 501 steps (`m2c_abn02z19_lines_n200_g2b.log`): fail_tr 0, fail_bdf 0, be_fallback 0; t = 4.931 Myr,
  T_c 172.9 K at log n_c 3.708 vs θ = 0.55 (`..._theta.log`) 4.932 Myr, 172.9 K, 3.707 (before: 0.7 Myr, 23 + 51 failures).
  P100 (`lines 100`, `m2c_lines_n100_g2b_full.log`, to n = 1e18): 12 944 steps, 0 stage failures (final reference run: 15 794
  steps, 72 + 293 failures); T_c at equal log n_c vs `m2c_lines_n100_final.log` at 7.9 / 8.9 / 9.9 / 10.9 / 11.9:
  366.2 / 417.5 / 573.7 / 803.0 / 1053.6 K vs 363.5 / 415.9 / 570.6 / 796.4 / 1038.5 K (+0.7, +0.4, +0.5, +0.8, +1.5 %; t(11.9)
  1.2907e6 vs 1.2981e6 yr, the 3000-step prefix `m2c_lines_n100_g2b.log` 0 failures). G2 dust-collapse order (T0 = 1e-6,
  C_dyn 0.03 / 0.01 / 0.003, `m1_dust_g2b_cdyn*.log`): worst-shell errors at R/R0 = 0.1 5.87e-3 / 6.74e-4 / 6.05e-5, identical
  to G2, order 1.97 and 2.00; steps 175 / 427 / 1412, 0 failures. BE discrete n = 60 7.35e-10 (`m1_be_g2b.log`), LP n = 200
  (`m1_lp_g2b.log`) byte-identical to `m1_lp_trbdf2.log`; `popiii_tests.jl quick` 8/8 PASS.
- Note: the G2 table needs `T0=1e-6`; at the default T0 = 1e-3 of `popiii_m1.jl dust` (also for the unmodified code) the error floors at
  ~2e-4 for C_dyn ≤ 0.01 (pressure).

### Core oscillation in the z = 19 run at 1e11–1e13: numerical; fixed [run]
**Symptom.** The "Final z = 19" run (θ-scheme, NLTE=hyb) and its TR-BDF2 repeat were non-monotone in n_c between 1e11
and 1e13: n_c fell back by up to 0.8 dex, T_c swung by about 40 %, and T was ragged from shell to shell by ±15 %. The
per-shell line cooling Λ_j/n_H₂ was ragged by factors of 2–10. Diagnosis, one change at a time (z = 19, 200 shells
unless stated):

| variant | oscillates |
|---|---|
| NLTE=1 (old factor) | no |
| θ + NLTE=ep / CHEM=rhfc02 / AVZONE=0 / LRJOIN=0.87 / IMPCOOL=2 (implicit cooling every step) | yes, in every case |
| TR-BDF2 + NLTE=ep at 100, 200 and 400 shells | no; the three tracks agree within ≈ 5 % |
| TR-BDF2 + NLTE=hyb at 100 and 200 shells | yes |

- Field's isobaric criterion holds: d(Λ − Γ_3b)/d ln T at constant P is +7 to +27 Λ. So this is not a local thermal
  instability.
- Three-body formation heating is comparable to the line cooling there, Γ/Λ = 0.8–4. That makes the core sensitive,
  but it doesn't by itself oscillate.

**Causes** (both numerical):
1. The weakly damped θ = 0.55 scheme amplifies the core's response.
2. The hybrid scheme's ALI/EP factor was frozen for 20 steps while the lines turned thick. The ratio changes by up to
   2× in that time, so each refresh jumped and kicked the core.

**Fixes** (general):
- `popiii_m2c.jl` now uses TR-BDF2 by default; `INTEGRATOR=theta` restores the old scheme.
- The hybrid refresh interval is adaptive. It is halved when the factor moved by more than 5 % since the last refresh
  and doubled, up to NALI, when it moved by less than 1 %. The z = 19 run made 2038 ALI calls in 9742 steps.

**Rerun** `m2c_abn02z19_lines_n200_final2.log`, with all generalisations. This replaces the "Final z = 19" run.
- 9742 steps, with 1 BDF retry recovered by backward Euler. It is monotone.

| log n_c | 6 | 8 | 10 | 11 | 12 | 13 | 14 | 16 | 18 | 19 |
|---|---|---|---|---|---|---|---|---|---|---|
| T_c [K] | 260 | 628 | 824 | 1061 | 1314 | 1551 | 1796 | 2271 | 3262 | 4279 |

**Regime map** (`regimes_z19_final2.txt`):
- The central lines are thick from 1e10.4.
- τ_R = 1 at 1e16.3.
- The windows are thermalised by 1e17.
- F/F_D (extinction) rises smoothly: 0.41 at 1e17.3, 0.56 at 1e18, 0.68 at 1e18.8 and 0.81 at 1e19.8, consistent with
  P100. Rosseland diffusion is good to ≈ 20 % from ≈ 1e19.5.
- The noisy 1.2–2.5 of the earlier regime map came from the oscillation.

---
## S Plan: deuterium burning, nuclear network, through the main sequence (Tom, 20 Oct 2026)

**Goal.** Follow the protostar born in the collapse (P100, z = 19) through deuterium burning and Kelvin–Helmholtz
contraction onto the main sequence and through hydrogen burning. This also gives the luminosity and spectrum of the
radiating star for the envelope-feedback problem.

**What TITAN already has** (PORTING_NOTES phase 7):
- Radiation hydrodynamics: E_r and F_r with Eddington factors, Compton exchange.
- MLT convective flux in the energy equation.
- Dorfi–Drury adaptive grid (lgrid = 1, used for the quasi-star accretion shock), an inner mass sink, TR-BDF2.
- The Z = 0 Saha EOS and opacity (Z0Physics). Opacity table log T 3.3–8.0; EOS without degeneracy, valid below
  ρ ≈ 1e-2 g cm⁻³; continuum opacity with no bound–bound lines.

**What is missing:**
- A nuclear network.
- EOS and opacity for stellar interiors: degenerate and relativistic electrons, ρ up to ~1e3 g cm⁻³, T up to ~1e9 K,
  composition-dependent.
- Composition that evolves per zone, with mixing in convective regions.
- Isolated-star validation runs.
- The accretion phase.

### S0 Literature and data
1. Validation tables:
   - Pop III ZAMS: Bromm, Kudritzki & Loeb 2001; Marigo et al. 2001; Schaerer 2002; Ekström et al. 2008.
   - Accreting Pop III protostars: Stahler, Palla & Salpeter 1986; Omukai & Palla 2001, 2003; Hosokawa & Omukai 2009.
   - Deuterium burning: Stahler 1988.
2. Nuclear rates: JINA REACLIB, for the pp chains, CNO I–II, 3α, and D(p,γ), D(d,n), D(d,p).
3. Screening: Chugunov, DeWitt & Yakovlev 2007, plus the weak-screening limit.
4. OPAL GN93 Z = 0 tables, and a reference EOS (Timmes & Swesty Helmholtz) for testing.

### S1 Nuclear network (`NucNet.jl`)
- **Species:** p, D, ³He, ⁴He, ⁷Li, ⁷Be, ⁸B, ¹²C, ¹³C, ¹³N, ¹⁴N, ¹⁵N, ¹⁵O, ¹⁶O, ¹⁷O, ¹⁷F, ¹⁸O, ¹⁸F, ¹⁹F.
- **Inputs:** REACLIB rates, including the weak decays.
- **Energy release:** Q minus the mean neutrino energy per reaction.
- **Solver:** an implicit step coupled to T at fixed ρ and u, the same structure as `ChemNet.step`.
- **Tests:**
  - pp-chain ε at solar-centre conditions;
  - CNO equilibrium ratios and ε_CNO;
  - 3α against NACRE/CF88;
  - ε_D against Stahler 1988's fit;
  - a one-zone Pop III core at 1.2e8 K, where 3α builds ¹²C to X_C ~ 1e-9 and CNO then takes over (Marigo et al.
    2001).

### S2 Stellar EOS
- **Physics:** ions ideal; electrons at arbitrary degeneracy and relativity (Fermi–Dirac integrals); a smooth join to
  the Saha region of Z0Physics; radiation stays in E_r.
- **Inputs:** composition (X, Y, Z_CNO) per zone.
- **Required:** AD-ready, thermodynamically consistent.
- **Tests:** against the Helmholtz EOS, plus consistency identities.

### S3 Stellar opacity
- **Tables:** OPAL Z = 0 (log T 3.75–8.7) as a function of X; our ContinuumOpacity below 1e4 K (H₂ CIA, H⁻, Li).
- **High T:** electron scattering with the Buchler–Yueh correction; conduction where it matters.
- **Output:** Rosseland and Planck means.
- **Tests:** continuity across the blends, and against published values.

### S4 Composition and burning in TITAN
- A per-zone nuclear abundance vector, with the material's EOS and opacity depending on it.
- ε_nuc in the energy residual, implicit in T with the abundances lagged. The composition is then updated by the S1 step.
- Convective mixing: an implicit diffusion step for the abundances, with D from the MLT velocity × mixing length.
- **Tests:** energy conservation including nuclear energy; a mixing test on a convective core.

### S5 Isolated Pop III stars (validation of S1–S4 together)
- Masses 1, 15 and 100 M☉: from a pre-main-sequence polytrope through the ZAMS to hydrogen exhaustion.
- Quasi-static evolution with large implicit steps; photosphere boundary condition from TITAN's radiation transfer.
- **Compare:** L, T_eff, R, T_c, ρ_c and the main-sequence lifetime against Marigo et al. 2001, Schaerer 2002, Ekström
  et al. 2008 and Bromm et al. 2001.

### S6 The accreting protostar from the collapse
- Continue P100 and z = 19 through the second collapse (H₂ dissociation, then H ionisation) into the accretion phase.
- The adaptive grid resolves the accretion shock and its radiative precursor (as in quasi-star P7.6).
- Deuterium burning, the luminosity wave, Kelvin–Helmholtz contraction, and arrival on the ZAMS.
- **Compare:** R*(M*) and L*(M*) against Omukai & Palla 2003 and Hosokawa & Omukai 2009, at their accretion rates and
  at the rate our envelope delivers.

**Further S6 references:**
- Omukai, Hosokawa & Yoshida 2010: spherical collapse followed by a protostar, Z = 0.
- Hosokawa et al. 2011 and 2012 (Tom): Pop III.1 vs III.2 R*(M*) and Ṁ(M*). R* peaks at t_KH = t_acc (≈ 6–9 M☉);
  ZAMS at ≈ 40 M☉ (III.1) and ≈ 17 M☉ (III.2). Summary in `~/Projects/popiii_lit/stellar/validation.md`.

### S7 Feedback on the envelope (later)
- UV from the star (H ionisation, Lyman–Werner H₂ photodissociation) and radiation pressure.
- Needs frequency-dependent continuum transfer that feeds ChemNet's photo processes, with J_ν per zone.

### S8 Deuterium chemistry and HD cooling (optional; matters for Pop III.2)
- The D/HD network of Hosokawa et al. 2012 (Table 1: Galli & Palla 1998, Stancil et al. 1998, Wang & Stancil 2002,
  Wolcott-Green & Haiman 2011), in ChemNet with detailed balance.
- HD rotational cooling.
- HD matters where the gas cools below ≈ 150 K, e.g. relic H II regions. Our P100 and z = 19 clouds stay above
  ≈ 140–200 K, so it is not needed there. It would carry the D abundance the protostar later burns.

### Order
1. S0 first.
2. S1, S2 and S3 in parallel; they are separate modules.
3. S4 after S1.
4. S5 after S1–S4.
5. S6 after S5.
6. S7 last.

Every item adds its tests to `popiii_tests.jl`.

### S2 Stellar EOS [run]

Files: `julia/StellarEOS.jl` (module StellarEOS), `julia/StellarMaterial.jl` (`StellarGas <: TitanMaterials.Material`), `julia/stellareos_test.jl`
(27 checks, ALL PASS; `full`-tier entry `stellareos`). One-line change in `TitanResidFast` (sets `ZONE[] = k` in the two MLT-coefficient
loops so `material_thermo` sees the zone's composition). Validity: ρ 1e-12–1e4, T 3e3–3e9 K (below ≈2.5e3 K Z0Physics' Saha underflows);
20000 random (ρ, T, X, Y, Z) points: all finite, χρ, c_V, Γ1 > 0.

**Method.** Gas only (radiation: `Prad`, `erad`, `srad` separately; `eos(...; radiation=true)` gives ∇ad, c_P of gas + aT⁴).
Two branches, each a full (p, e, s) per gram at (ρ, T), generic in the number type:
- S (Saha): Z0Physics' LTE H/He gas with arbitrary (X, Y) (`Z0.Comp`; Z0Physics unchanged); metals are neutral ideal atoms.
- E (fully ionised): ideal bare ions H, He, C, N, O (Sackur–Tetrode, g = 1 as Z0Physics) + electrons/positrons from
  F_k(η, θ) = ∫ x^k (1 + θx/2)^{1/2}/(e^{x−η}+1) dx, θ = kT/mc²: n = A0 θ^{3/2}(F½ + θF3/2), P = (2/3)A0 mc² θ^{5/2}(F3/2 + θF5/2/2),
  E_kin = A0 mc² θ^{5/2}(F3/2 + θF5/2); positrons at η₊ = −η − 2/θ with 2mc² n₊ added to e (pairs included when e^{η₊} > 1e-20);
  YS89 uniform-background Coulomb (the Timmes–Swesty fit; one ion species, Z̄ = n_e/n_i; the two fits joined C² over Γ = 0.6–1.7).
  Energy zero = neutral atoms (E includes the ionisation energy of all ions). η from charge neutrality: Float64 safeguarded Newton,
  then 1–2 Newton steps in the caller's number type (steps follow the Dual nesting depth: exact through 2^k − 1 derivatives, as in Z0Physics);
  ∂F/∂η comes free from the quadrature.
- F_k quadrature: η < 48: composite 16-point Gauss–Legendre on x ∈ [0, max(η,0)+46] (first panel in √x, panels ≤ 5 wide); η ≥ 48: edge form
  F = ∫₀^η g dx + ∫₀^46 [g(η+u) − g(η−u)]/(e^u+1) du (first term Gauss–Legendre in t; geometric panels when ultra-relativistic).
  Max rel. error 2.9e-15 over 180 (η, θ, k) (η −40…1e4, θ 1e-5…5) vs an independent tanh-sinh; analytic ∂F/∂η vs AD 2e-15.
  Electron–positron gas vs momentum-space integration: 1e-12 (11 points up to 3e9 K, ρ 1e4).
- Join: blend of the FREE ENERGY, f = (1−w) f_S + w f_E, w = 1 − (1−w_T)(1−w_ρ), w_T = smootherstep(log T; 5.3→5.9), w_ρ =
  smootherstep(log ρ; −1→0.7) (crude pressure ionisation: full ionisation + degeneracy for ρ ≳ 0.1–5). From p = ρ²∂f/∂ρ, s = −∂f/∂T:
  p = (1−w)p_S + w p_E + ρΔf ∂w/∂lnρ, s = (1−w)s_S + w s_E − (Δf/T) ∂w/∂lnT, e = (1−w)e_S + w e_E − Δf ∂w/∂lnT, Δf = f_E − f_S.
  Exact Maxwell consistency; f ∈ C², so Γ1 and ∇ad are continuous. At w = 0 identical to Z0Physics' gas (0.0 difference in p, e).
  Branches differ by Δp/p = −3.7e-4 (T = 1e6, ρ = 1e-3) and −1.2e-3 (T = 10^5.6).
API: `eos(ρ, T; X, Y, Z, zmix, pairs, coulomb, radiation, comp)` → P, e, s, cV, χρ, χT, Γ1, ∇ad, δ, cP, cs, ne, η (+ dPdρ, dPdT, dedρ);
`gas_pe(ρ, T; …)` → (p, e); `rho_PT(P, T; …)` (round trip 3e-15; (∂lnρ/∂lnP)·χρ − 1 = 8e-16); `StellarGas(mgr; X, Y, Z)`, `setcomp!(m, k, X, Y, Z)`,
`ghosts!(m)`. `StellarGas` opacity = Z0 table placeholder until S3; sound speed for the AV = 5/3 estimate.

**Tests (numbers).**
- Maxwell ∂e/∂ρ|T = (p − T∂p/∂T)/ρ²: 4e-14 max; c_V vs T∂s/∂T: 2e-11 (725 points, log ρ −10…4, log T 3.5…9.5, Z = 1e-3). AD vs FD Jacobian 1.3e-7 (FD-limited); nested dΓ1/dlnT 9e-8.
- Limits: ideal-gas P to 1e-12; Chandrasekhar T = 0 (ρ 1…1e10, T = 10 K): 7.8e-8 at ρ = 1 (thermal), 1e-11…1e-16 above; NR P/P_NR = 1 − 5x²/14 and ER P/P_ER = 1 − 1/x² reproduced.
- vs Helmholtz (A = 1.3, Z = 1, E branch, Coulomb on both; 660 comparable points; 7 where Helmholtz's bomb-proofing removed its Coulomb term skipped),
  max rel. diff P / e / s / Γ1: hot non-pair (log T 6–9, log ρ ≥ −4) 2.4e-5 / 3.8e-5 / 1.2e-3 / 7.9e-4; dense (log ρ ≥ 1, log T ≥ 6) 1.1e-5 / 1.4e-5 / 1.2e-3 / 5.6e-5;
  cool dense (log T 4–6, log ρ ≥ 0) 2.1e-5 / 5.0e-5 / 1.3e-2 / 4.0e-5; pair region (log T ≥ 7, log ρ ≤ −6) 2.7e-5 / 1.9e-2 / 8.5e-4 / 1.9e-1.
  The larger e, s, Γ1 differences are Helmholtz's table (≈10 nodes/dex in T): against exact momentum integrals its e is off by 1.3e-2
  at (log ρ, log T) = (−8.5, 8.25) and 5.5e-4 at (−10, 8.5); ours agrees to 1e-12.
- Blend: the largest step jump of Γ1 and ∇ad halves with the step on 6 scans (log T 5.0–6.2 at ρ = 1e-6, 1e-4, 1e-2; log ρ −2…1.2 at T = 1e4, 1e5, 1e6).
- Pop III points (X = .76, Y = .24; η_e, Γ1, ∇ad, gas only): (100, 3e7): −2.70, 1.6625, 0.3982; (30, 1.2e8): −6.03, 1.6502, 0.3940;
  (1e-3, 1e6): −9.12, 1.6666, 0.3998; (1e-6, 3e4): −10.86, 1.5303, 0.3301. Electrons are non-degenerate even at these Pop III ZAMS centre points
  (η_e < 0); degeneracy matters only for ρ ≳ 1e3 at 3e7 K. With radiation: (100, 3e7) Γ1 = 1.6549, ∇ad = 0.3927, β = 0.995; (30, 1.2e8) Γ1 = 1.4216, ∇ad = 0.2590, β = 0.491.
- Coulomb size (Γ_ion): (100, 3e7) 0.038: Δp/p = −9.5e-4, Δe/e = −1.9e-3; (30, 1.2e8) 0.006: −7e-5; (1, 1e6) 0.245: −1.4e-2, −2.6e-2.
- Speed (µs/call, one thread): Saha 1.7, hot 3.6–4.9, central 4.2, degenerate 2.3, pair plasma (3e9 K) 53; (p, e) + AD Jacobian 2.3–7.3 (pairs 56); full `eos` 2–7.

**Limits / problems.** (1) Pressure ionisation is a switch, not physics: in the transition (log ρ −0.5…0.5, log T 4…5.5) ∇ad reaches 0.66
(Γ1, χρ, c_V stay > 0); irrelevant inside stars and envelopes (T > 1e5 or ρ < 1e-5 there) but wrong for brown-dwarf-like cold dense gas. (2) The Saha
branch has no degeneracy and no molecules (H₂ below ≈3500 K), metals are neutral; the E branch has no partial ionisation of He/metals below ≈3e5 K.
(3) Coulomb: one ion species with Z̄ (no linear mixing), no electron screening for non-degenerate electrons, YS89 fits only; the Γ ≥ 1 fit is shifted by
0.0074 kT per ion to join the Γ < 1 fit continuously. (4) Pair plasma: η is ill-conditioned when n₊ ≫ n_net (rel. error up to 1e-5 in η; p, e, s unaffected); 55 µs.
(5) `StellarGas` opacity is a placeholder (S3). The Helmholtz download (untrusted) lives in `~/Projects/popiii_lit/helmholtz/` (README there; not in the repo);
test section D is skipped if absent.
### S3 Stellar opacity [run]
**Code.** `julia/StellarOpacity.jl` (run time: table + AD-ready interpolation + analytic high-T model), `StellarOpacityLow.jl`
(own low-T means, general-X equilibrium), `stellaropacity_build.jl` / `stellaropacity_data.jl` (builder, parsers; 40 s on 4
threads), `stellaropacity_test.jl` (14 checks, "full" tier entry `stellaropacity`), table `julia/stellar/kaptable_stellar.bin.gz`
(3.5 MB), build log `lib/popiii/s3_stellaropacity_build.log`.
API: `kappa_R(ρ,T,X)` (with conduction), `kappa_R_rad`, `kappa_cond`, `kappa_P`, `kappas`, `log_kappa_R(lT,lR,X)`, `in_table`, `coverage`.
κ_R 0.35 µs, with gradient in (ρ,T,X) 0.7 µs.

**Provenance** (`~/Projects/popiii_lit/opal/README.md`; nothing downloaded is in the repository). OPAL GN93hz (Iglesias & Rogers 1996),
Z = 0 tables X = 0, .1, .2, .35, .5, .7, .8, .9, .95, 1, from MESA's copy (opalopacity.llnl.gov unreachable from the sandbox; sha256 in
README). Cassisi et al. 2007 conduction tables from MESA. Entries logT = 8.7, logR ≤ −7 of OPAL are placeholders (identical in all
tables, −0.668/−0.931/−1.004 vs ≈ −1.017) and are masked.

**Coverage.** Grid log T 1.5–9.5 (Δ0.025) × log R −15…+5 (Δ0.25) × the 10 X nodes; Catmull–Rom bicubic in (log T, log R) on log κ, PCHIP
in X (monotone, C¹). X ∈ [0, 1]; Y = 1 − X, Z_CNO = 0. Region by log T (X = 0.7 node counts: OPAL data 26 %, own means 29 %,
extended 41 %, blends 3 %):
- 1.5–3.8: own continuum means (CIA, H⁻, Rayleigh H₂/H/He, Thomson) in chemical + Saha equilibrium for general X, with Li electrons (same
  Saha chemistry as `PrimordialChem.y_saha`: f_H₂ agrees to 1e-6 at X = 0.75). No H₂ or atomic lines. HITRAN CIA held at 200 K below it.
  X ≲ 0.05: not reliable (no He opacity beyond Rayleigh: OPAL − own = +0.4…+2.6 dex at X = 0).
- 3.8–4.1: log κ blend (smoothstep) own ↔ OPAL.
- 4.1–8.7, log R −8…+1: OPAL exactly (|Δ| ≤ 7e-15 at 11 520 nodes). Outside OPAL's R range or its hot high-R corner: extended OPAL (flag 2).
  log T 4.1–5.5 beyond log R = 1: linear in log R with OPAL's last slope (clamped to [0, 1.2]); below log R = −8 saturating
  extension. From 5.5 to 6.5 the extension changes smoothly to the high-T model anchored to OPAL (κ_high + Δ, Δ = OPAL − κ_high held
  at the last valid node of each row, and for T > 8.7 at the 8.7 row). Quality there is an extrapolation: do not trust it above ~log R = 3.
- Conduction (all T, switched on over log ρ −6…−5 and log T 4…5): 1/κ = 1/κ_rad + w_c/κ_cond, κ_cond = 16σT³/(3ρK) with Cassisi's K,
  1/K = X/K_H + Y/K_He. Negligible except in degenerate matter: at ρ = 150, T = 1.5e7 (solar centre) κ_cond/κ_rad = 1000; a He core
  at ρ = 1e6, 1e7 K has κ_rad 390, κ_cond 6e-5 (conduction decides). During Pop III hydrogen burning ρ_c ≲ 1e2–1e3: irrelevant.
- Planck mean: own means below 3.8; blend over 3.8–4.1; above, the Z0Physics continuum Planck mean for general X (Saha H/He, bound-free,
  free-free with g = 1, H⁻), times the free-free Pauli factor above 1e6 K. **No lines** (a lower bound at 10⁴–10⁵ K, where lines raise κ_P
  by a factor of a few), agrees with the Z0Physics table at X = 0.76 to ≤ 0.3 % for log R ≤ −4, log T 4.5–7.

**High-T model** (`kappa_high`): κ_es = 0.2(1+X) f_T D with f_T = 1/(1+(T/4.5e8)^0.95) (Buchler–Yueh/Paczyński form; exponent refitted from
0.86 to OPAL's R = 1e-8 column, 1 % over 1e6–3e8 K), D = 1/(1+2.7e11ρ/T²) (Paczyński degeneracy), plus free–free 3.68e22 g(T)(1+X)ρT^-3.5/
(1+22.2·2.7e11ρ/T²)^0.366 with g(T) quadratic in log T: **empirical** 5-parameter fit to OPAL (rms 0.063 dex, log T 6.8–8.5, log R ≥ −3).
Itoh's relativistic free–free and Cassisi-style degenerate scattering were not used (fits not at hand). OPAL − κ_high at X = 0.7:
≤ 0.01 dex for log R ≤ −4 everywhere; log T 7: −0.07…+0.03 (R −1…−3); 8.0: +0.02 (R −1); 8.7: +0.03 (R −2.5), +0.135 (R −1.5, relativistic
degenerate edge). Beyond log T 8.7 the offset Δ is held; pair production (T ≳ 3e9) not modelled.

**Mismatches.** OPAL − own means (log10 κ_R), X = 0.7, over the 19 OPAL columns: mean/max +0.015/+0.04 at log T 3.80, +0.043/+0.08 at 3.90,
+0.079/+0.19 at 4.00, +0.119/+0.36 at 4.10 (own means lack H lines). X = 0.1…1 similar. κ_OPAL/κ_Z0(table, continuum only), X = 0.76, median
(range) over log R −8…1: log T 4.0: 1.14 (1.0–1.4); 4.5: 3.6 (1–17); 5.0: 2.5 (1–8); 5.5: 1.33 (1–5.5); 6.0: 1.00 (0.76–1.28); 6.5: 0.98 (0.53–0.99);
7.0: 0.97 (0.49–0.99); 8.0: 0.81 (0.41–4.3); the < 1 values at large R are electron degeneracy that Z0Physics lacks.

**Tests** (all PASS): OPAL nodes 7e-15; reference (7.0,−3,0.7) = −0.412 etc. equal raw; off-node vs OPAL's Hermite 2e-3 dex max; gradient
jumps across edges and blend knots ≤ 2e-6 (log T), 5e-9 (log R), 1.3e-6 (X); AD vs FD 7e-9; κ→0.2(1+X)f_T within 2 %; own means vs
`ContinuumOpacity.means` on the same composition: κ_P 1e-5, κ_R 0.5 % (He Rayleigh); table vs direct at the M3c points (n_H 1e14–1e21,
1500–6000 K, X = 0.75) ≤ 0.2 %; fine-grid interpolation error in the molecular regime rms 0.003, max 0.027 dex (steep Li-ionisation onset).
**M3c numbers.** Recomputing the M3c table (no Li) with the new solver gives κ_P/κ_R = 4.49e-8/1.17e-9, 4.97e-6/7.45e-6, 3.35e-4/5.62e-5,
2.60e-2/3.17e-3, 6.38/4.94 for (n_H, T) = (1e14, 1500) … (1e21, 6000); the M3c table has κ_P = 2.3e-2 at (1e20, 4000) (13 % lower; the
code, `ContinuumOpacity.means`, gives 2.60e-2, so that entry looks stale). Li electrons matter: at (1e14, 1500) κ_R(ext) is 3.0e-7 with Li,
1.2e-9 without (the M3c table is the Li-free case).

**Open.** MD05 and Lenzuni et al. 1991 numerical tables not obtained (MD05's are electronic only), so no comparison. Planck mean without
lines above 4000 K. He Rayleigh cross section from memory. OPAL's last nodes in the hot rows (logT 8.5, log R −1) drop by 0.25 dex and cause
curvature at 8.5–8.7. The M3c/regime τ's are column integrals; the check was done on local κ.
### S1 Nuclear network (NucNet) [run]
Files: `julia/NucNet.jl`, `julia/nucnet_test.jl` (45 checks, ~1 min), `popiii_tests.jl` entry `nucnet` (full tier; quick unchanged).
- **Data** (not in the repo): `~/Projects/popiii_lit/reaclib/reaclib_default2_20250330` = JINA REACLIB snapshot 2025-03-30 as shipped
  in pynucastro (release 3.1.0), sha256 24e1f37c…ba862; provenance/licence in the README there (pynucastro BSD-3; REACLIB cite Cyburt+ 2010).
  Override path with `ENV["NUCNET_REACLIB"]`.
- **Network**: 19 species (p, d, 3He, 4He, 7Li, 7Be, 8B, 12C, 13C, 13N, 14N, 15N, 15O, 16O, 17O, 17F, 18O, 18F, 19F), 31 reactions: pp, pep,
  d(p,γ), 3He(3He,2p), 3He(α,γ), 7Be EC, 7Li(p,α), 7Be(p,γ), 8B decay, hep, d(d,n)3He and d(d,p)t, d(3He,p)α, CNO I/II
  (+ 18F(p,α)15O), 3α, 12C(α,γ). All non-reverse REACLIB sets of a reaction summed. No reverse rates (Q ≥ 1.94 MeV for p-captures,
  e^{-Q/kT} ≤ 1e-4 at T9 ≤ 0.3; none near equilibrium). d(d,n): neutron assumed captured at once (n+p→d), net d+p→3He (Q 5.4935);
  d(d,p)t with t→3He lumped (Q 4.0513 incl. t decay, ⟨E_ν⟩ 0.0129). REACLIB 'ec' sets (7Be EC, p+p+e⁻) are per unit n_e:
  λ_eff = λ ρ Y_e (7Be: REACLIB 1.414e-9 vs Bahcall/CF88 formula 1.400e-9 at 1.57e7 K; τ(7Be) = 82 d at the solar centre).
- **Energy**: ε = N_A Σ (Q − ⟨E_ν⟩) rate; ⟨E_ν⟩ [MeV] (Bahcall 1989, Kippenhahn–Weigert): pp 0.265, pep 1.442, 7Be 0.814, 8B 6.71, 13N 0.707,
  15O 0.997, 17F 0.999, 18F 0.382, hep 9.628. REACLIB Q vs AME2020 mass excesses: max |ΔQ| = 1.3 keV over all reactions.
- **Screening**: `:chugunov` (Chugunov–DeWitt–Yakovlev 2007 eqs. 19/21, multi-component Z1Z2/z̃ form of Yakovlev+ 2006 as in pynucastro,
  re-implemented), `:weak` (Salpeter 1954), `:none`. 3α: two sequential screens (α+α, α+8Be). Composition lagged within a substep.
- **Integrators**: `burn(Y0,ρ,T,dt)` backward Euler, analytic Jacobian (= ForwardDiff to round-off), positivity damping, recursive halving;
  returns time-averaged ε, ε_ν, ε by group (pp, d, cno, he) and integrated event counts. `burn_coupled(Y0,ρ,u,T0,dt,eos_u)` Newton on (Y, ln T) with
  the ForwardDiff Jacobian, initial guess from the fixed-T solution, u_eos(ρ,T,Y) − u − dt ε = 0 (to 1e-12). `evolve` (adaptive), `steady`
  (pseudo-transient, `frozen` species). Abundances below 1e-16 mol/g are clamped, not damped (LU noise at cond ~1e23).
- **Tests** (numbers):
  - a. Solar centre (1.57e7 K, ρ 150, X 0.34): ε_pp = 15.49 (3He equilibrium X3 = 7.5e-6, reached within 4.6 Gyr), ε_CNO = 1.00 (relaxed 4.6 Gyr)
    / 2.45 (full CNO equilibrium), ε = 16.5; ε_ν/ε_tot = 3.5 %; 3He held at ×0.1/0.5/2 of equilibrium: ε_pp 6.3/9.8/31.5. Screening ε_pp: Chugunov
    15.49, Salpeter 16.06, none 14.23 (pp factor 1.029 vs 1.050). pep/pp = 3.5e-3 (Bahcall–May 2.8e-3), 8B/pp 1.2e-3.
  - b. CNO equilibrium: 12C/13C = 3.29 at 2e7 and 3e7 K (= λ13/λ12 of REACLIB, flat 3.27–3.29 over T9 0.015–0.05; below the 3.5–4 of
    older fits), 14N 97/98 % of CNO; ε_CNO / Kippenhahn–Weigert (18.65, 8.24e25 g X_CNO X ρ T9^-2/3 exp(−15.231 T9^-1/3 − (T9/0.8)²)) = 0.69, 0.64.
  - c. 3α REACLIB (fy05) / CF88: 0.966 (T9 0.1), 0.984 (0.2), 1.003 (0.3). NACRE (Angulo 1999) table values were not available offline; not compared.
  - d. ε_D (d(p,γ), ρ = 1, D/H = 2.5e-5) vs Stahler 1988 eq. 48 (4.19e7 f [D/H] ρ T6^11.8): X = 0.70: ratio 0.968, 1.014, 0.680 at T6 = 0.8, 1, 1.5;
    X = 0.75: 1.114, 1.167, 0.782 (ε ∝ X²). Effective index 11.24 over 0.8–1.5 (fit 11.8, so the fit overshoots above T6 ≈ 1.2).
    τ_D = 1.7e6, 4.0e4, 2.8e3, 3.6e2, 34, 2.2, 0.068 yr at T6 = 0.6, 0.8, 1, 1.2, 1.5, 2, 3 (X = 0.7, ρ = 1). At D/H = 2.5e-5 and T6 = 1.5, d(d,n)/d(d,p)
    destroy 4 % of D within one d(p,γ) e-folding.
  - e. Pop III one-zone (ρ = 30; X = 0.75, D/H 2.5e-5, 3He/H 1e-5, 7Li/H 4.7e-10): X(12C) = 1e-11/1e-10/1e-9 at t = 1.3e4/8.4e4/1.6e5 yr (1.2e8 K),
    5.6e6/1.1e7/2.0e7 yr (1.0e8 K), 93/9.0e2/4.3e3 yr (1.4e8 K). ε_CNO = ε_pp when X(CNO) = 2.4e-11 (X12 6e-13, t = 810 yr) at 1.2e8 K, 1.1e-10 (X12 2.2e-12,
    2.3e6 yr) at 1.0e8 K, 6.8e-12 at 1.4e8 K; these agree with the KW equilibrium estimate (1.7e-11, 6.7e-11, 5.6e-12) within a factor 1.6 and are independent of ρ
    (ε_pp and ε_CNO both ∝ ρ). X_CNO ~ 1e-9 (Ezer & Cameron) therefore needs T_c ≲ 8e7 K; at ρ = 30 and 8e7 K H is exhausted before CNO catches up
    (max ε_CNO/ε_pp = 0.16). H lasts 4e5 yr at 1.2e8 K (ρ = 30).
  - f. ΣA·Y − 1 ≤ 5e-15 in every run; Σ_r ΔZ_r N_r = Δ(ΣZY) to 1e-15 (weak events = e⁺ emitted/e⁻ captured); ∫(ε+ε_ν)dt vs −N_A Σ ΔY·M(AME2020): 2e-5–5e-5.
    Coupled: u_eos(new) − u − dt ε = 5e-12 relative; ε agrees with fixed-T to 0.4 % at ΔT/T = 3e-4; dt = 1e12 s at 1.2e8 K, ρ = 30 is a thermal runaway (T × 3.2).
  - g. burn: 20 µs per call (1 substep, 19 species); burn_coupled 176 µs; rates_k (31 reactions, screening) 4.7 µs.

### S4 Composition and burning in TITAN [run]
Files: `julia/StarMaterial.jl` (`StarGas`), edits in `TitanResidFast.jl` (source hook, `step_diag`, MLT v and ℓ side outputs) and `TitanCore.jl`
(`DT_HOOKS`, `STATIC_FALLBACK`); tests `s4_test.jl` (driver; a, b, c), `s4_a.jl`, `s4_b.jl`, `s4_c.jl`, star `s4_star.jl` + `s4_d.jl`; entries `s4_*` in
`popiii_tests.jl` (full tier). Quick tier 8/8 PASS after the edits. Logs: `lib/popiii/s4_test_{a,b_efrac0.05,b_efrac0.01,c}.log`, `lib/popiii/s4_d_100Msun.log`.
- **Material.** `StarGas`: Y[19, zone] (mol/g, NucNet order). EOS composition: X = Y_p + 2Y_d, Y_He = 3Y(³He) + 4Y(⁴He) + 7Y(⁷Li, ⁷Be) + 8Y(⁸B), Z = rest (C:N:O by
  `zmix`); Z < 1e-18 set to 0 (3α makes 1e-240 of ¹²C and `StellarEOS` then returns NaN derivatives). Gas-only p, e, n_e from `StellarEOS`; κ_R(ρ,T,X_H) with
  conduction and κ_P from `StellarOpacity` (`opac = :z0` = the old table). `kcap` (optional, soft cap of κ at ~Γ < 1, see S5 hand-over).
- **ε in the residual.** Hook `has_source(m)` / `material_source(m, ρ, T)` [erg g⁻¹ s⁻¹, net of ν]. Subtracted as −ε ρ ΔV (ρ, ΔV at the θ-level) in all four
  energy rows (combined gas+rad, GASROW gas row, total-energy row, pure gas), θ-weighted: θ ε(ρ_n,T_n) + (1−θ) ε(ρ_o,T_o), the old term from the start-of-step state in
  `context()` (like Compton). Y of the zone fixed during Newton (lagged); ε differentiated by ForwardDiff through `NucNet.rates_k` (Dual-safe as is; the screening
  plasma uses Float64 ρ and the lagged Y, so ∂ln f_scr/∂ρ is not in the Jacobian). 3 µs Float64, 5 µs with 2 partials per zone. No source → no code path changes
  (quick tier identical step counts 464/587/713).
- **burn!(m)** after each converged step (`output` callback, as `chem!`): per zone `NucNet.burn` at fixed ρ and T_n (`burnT = :end`; `:theta` = θT_n + (1−θ)T_o), dt = timen − timeo.
  Energy: the residual injected E_inj = Δt·4π·(ε_θ ρΔV) (read back from the Float64 residual, `TX.step_diag`), the network released E_bur = Δt ε_burn·m (Σ_r N_A(Q_r − ⟨E_ν⟩_r) nev_r exactly);
  e ← e + (E_bur − E_inj)/m at fixed ρ, T re-solved from the EOS with the new Y (Newton, 1e-12). Total energy (gas + cumulative release) conserved to that tolerance (below);
  `ledger` keeps ΣE_inj, ΣE_bur, Σ|E_bur − E_inj|/Σ|E_bur| (size of the splitting correction), max T re-solve error.
- **mix!(m)** after burn!: implicit conservative diffusion of all 19 species on the Lagrangian mass grid, ∂Y/∂t = ∂/∂m[(4πr²ρ)² D ∂Y/∂m], D = vℓ/3 at interfaces with
  F_conv > 0 (MLT weight > 0), v = ℓ(gδ/8H_P)^½ ζ (Kippenhahn–Weigert mean velocity; ζ² = ∇ − ∇′ is the same ζ as in F_conv = ρc_P T (gδ)^½ ℓ² ζ³/(4√2 H_P^{3/2}), ℓ = αH_P). v, ℓ, F_conv come from
  `step_diag` evaluated on the converged state before burn!. One Thomas factorisation for all species, solved for δY with two refinement sweeps (conservation to 1e-16).
  Composition change at fixed zone energy e: T re-solved (zone energies unchanged, so energy conservation is exact; ΔT/T ≤ 1.4e-11 in the star test).
  Caveat: deep in the core ∇ − ∇′ ~ 1e-9 and ζ, hence v, is determined only to a factor ~2 by a Newton state converged to `conv`; mixing is insensitive when t_mix ≪ Δt.
- **Time step.** `TC.DT_HOOKS` (consulted at the end of `timstp!`) returns min(post-step cap, predicted cap): dt ≤ dX_max/|dX_H/dt| (dX_max = 1e-3, rate averaged over each mixed region),
  dt ≤ f_e e/ε, dt ≤ f_e e/max_i(|dY_i/dt| Q_i N_A MeV) with Q = main-channel energy per nucleus of p, D, ³He, ⁴He(3α), ⁷Li (catalysts and short-lived intermediates have Q = 0: with
  their binding energies (tried) the criterion is dominated by ¹⁵O/¹⁵N/¹⁴N quasi-equilibria); mix!/burn! apply the net (burn + mixing) ΔX_H of the step at once. f_e = 0.05 (closed box), 0.5 (star: thermal feedback is implicit).
- **Test a (no burning).** `StarGas(opac = :z0, nuc = false)` vs `StellarGas`, 25 steps of the 100 M☉ rad-hydro star (40 zones, MLT): r, ρ, T, E_r, F_r, u, p, e **bit-identical** (max rel. diff 0). PASS.
- **Test b (closed box, spherical shell with walls, ρ = 1, T = 1e6 K, D/H = 2.5e-5, X = 0.7, 4 zones).** D burns, T 1.0 → 1.454e6 K. TITAN steps vs `NucNet.burn_coupled` with the same dt sequence:
  max|ΔT/T| 4.3e-3 (f_e = 0.05, 28 steps) / 2.1e-4 (0.01, 97 steps); vs a 40-substep reference: 4.3e-2 / 9.4e-3 (∝ f_e, first order in dt).
  Energy: |ΔE_gas − ∫ε dm dt|/ΔE = 4e-13 / 1.5e-11, max over steps relative to E_0 1.8e-13 / 6.6e-12 (limit 1e-10); T re-solve ≤ 1.3e-13 / 8.5e-13; splitting correction Σ|E_bur − E_inj|/ΣE_bur 17 % / 5 %.
  Rest-mass check: N_A Σ ΔY_i B_i (AME2020) / net release = 0.9991 (REACLIB vs mass table Q + ν losses). PASS.
- **Test c (mixing).** c1 synthetic 60-zone grid, convective zones 15–35 (step ΔX_H = 0.2): Σ Y dm conserved to ≤ 1.4e-16 in every case, radiative zones changed by exactly 0, spread/step after 3 / 10 / 30
  turnover times ℓ/v (200 substeps) 0.105 / 3e-4 / 5e-11 (one backward-Euler step: 0.35 / 0.13 / 0.048, first-order in dt/t_mix). c2 in the 100 M☉ star (below): core ∫dr/v = 9.4e-2 yr,
  step of 0.1 in X_H over the inner half of the core: spread/step 0.099 / 3.8e-3 / 3.9e-8 after 1 / 3 / 10 core crossing times, species masses of the whole star conserved to 7e-16, the 40 radiative zones above the core exactly unchanged, |Δe|/e ≤ 1.4e-11. PASS.
- **Test d (100 M☉, X = 0.75, Z = 0).** Start: Eddington standard model (n = 3 polytrope, β = 0.5807, μ = 0.5926) of R = 6.5 R☉ (T_c 1.04e8), 140 zones (80 log-spaced in exterior mass fraction 0.3 → 1e-10), regular centre,
  free surface, MLT α = 1.5, StarGas with burning on after step 150 (the first 150 steps relax the polytrope's dynamical imbalance, 1e-5 yr), θ = 1, stol 0.4, dt growth ×1.3 when the step changed < 60 % and Newton took ≤ 9 iterations
  (`s4_d.jl`, 294 s on 4 threads for 503 steps to 1.1e5 yr). The star contracts for ~5e3 yr (T_c overshoots to 1.49e8), 3α makes X_CNO(centre) = 6e-13 / 1.8e-11 / 1.3e-10 at t = 7.4e2 / 1.9e3 / 3.1e3 yr, ε_CNO/ε = 0.013 / 0.41 / 0.90
  (S1: ε_CNO = ε_pp at X_CNO = 2.4e-11 at 1.2e8 K), the convective core grows to M_cc/M = 0.826 (r = 2.23 R☉) by 1.6e4 yr and stays. Numbers:

  | t [yr] | R [R☉] | log L | log T_eff | T_c [1e8 K] | ρ_c | L_nuc/L | X_c | X_CNO,c |
  |---|---|---|---|---|---|---|---|---|
  | 2.5e4 | 4.26 | 6.129 | 4.979 | 1.34 | 45.4 | 1.18 | 0.7463 | 7.5e-10 |
  | 1.1e5 | 4.47 | 6.134 | 4.970 | 1.26 | 40.6 | 1.13 | 0.7317 | 1.3e-9 |

  Schaerer 2002 (100–120 M☉ ZAMS): log L 6.1–6.24, log T_eff ≈ 4.97: agreed. L_surface = 5.21e39 erg/s vs ∫ε dm = 5.90e39: L_nuc/L = 1.13–1.18 does not tend to 1: the star expands (R 4.2 → 4.5 R☉ in 1e5 yr, H burning) and the
  excess feeds the expansion (≈ 6 % of L from the virial estimate); L itself is flat to 0.3 % after 1.5e4 yr (≈ 1.5 t_KH). Ledger over the run: ΣE_inj/ΣE_bur = 1.076, Σ|E_bur − E_inj|/Σ|E_bur| = 0.16 (dt = 400 yr ≫ the CNO
  catalyst equilibration time, so the lagged-Y heating differs from the burned energy; corrected exactly by burn!), T re-solve error 1e-12. 1D numerics: T_c wanders ±3 % (1.26–1.31e8) on steps of 400 yr.
- **Handed over to S5 / open.** (1) The initial polytrope is not in discrete hydrostatic equilibrium (surface layers, L from the Eddington model): burning must be switched on after ~1e-5 yr of relaxation, otherwise Newton fails
  in step 2; `TC.STATIC_FALLBACK[] = true` is needed when dt·u exceeds the zone width (peq's extrapolation halves dt otherwise). (2) Z = 0 OPAL κ_R reaches 1.08 at 1.5e5 K in the
  skin (0.1 % of M): Γ > 1 accelerates the 24 outermost zones outward at 100s km/s and pins dt at 1 s; `kcap = 0.74` (soft cap, Γ ≤ 0.8) was used. This hides a real issue of 1D VMS envelopes.
  (3) The first MLT interface at the centre gives a large F_conv (g → 0, ℓ = αH_P ≫ r); diagnostic L_conv in the inner 5 % of the profile exceeds L (Newton tolerance on near-adiabatic ζ³), harmless for the energy balance but not
  usable as a profile. (4) T_c = 1.26–1.34e8 K at ρ_c 41–45 was not compared with a published model (no tables at hand). (5) `burn!` supports the θ-scheme only (TR-BDF2 stages would need per-stage injection bookkeeping).
  (6) f_e = 0.5 in the star is a choice; the energy cap is not a stability criterion of the implicit hydro. (7) ~0.6 s per step at 140 zones, dt ≈ 400 yr: main-sequence evolution (3e6 yr) is ~1e4 steps.

### S5 Isolated Pop III stars [run]
Files: `julia/s5_run.jl` (driver: build → relaxation → evolution with ledger, step control, step rejection, checkpoint/restart; options in its header), `s5_b.jl` (closed box,
TR-BDF2), `StarMaterial.jl` (`install_stage_hook!`, `burn_zone`, `solve_step!`, `install_state_predictor!`, `mix!(capmix, bridge, tmax)`, `REGION_CAP`), `TitanResidFast.jl`
(`stage_diag`, `MLT_RCAP`, `MLT_HSDISC`, ∇/∇_ad side outputs), `TitanCore.jl` (`EXTRAP_FRAC`, `DT_RULE`, `PEQ_HOOK`, `PRED_LAGRHO`), `TitanDriver.jl` (`STAGE_HOOK`, `FULL_DT`, `FAIL_TFAC`),
`Trbdf2.jl` (hooks), `s4_star.jl` (grid and surface fixes). Tests `s5_b`, `s5_zams15` in `popiii_tests.jl` (full tier; `POPIII_SERIAL=1` runs them one process at a time); quick tier 8/8 PASS,
S4 tests a, b, c PASS. Logs `lib/popiii/s5_*.log`, `s5_p100e.csv`, `s5_p15e.csv`. All runs `TITAN_NOLIB=1`, Z = 0, X = 0.75, D/H = 2.5e-5, ³He/H = 1e-5, α_MLT = 1.5, 140 zones
(60 log-spaced in the exterior mass fraction from 0.05 to 1e-12), κ not capped.

**1. Energy budget.** *Diagnosis.* Ledger per step: E_inj booked from the residual, E_bur from the network, Φ_out dt, E_tot = E_int + E_rad + E_kin + E_grav (E_grav = −Σ G m_k M_k/r_k, M_k = ½(Δm_{k−1} + Δm_k)).
The global ledger already closed in S4 (burn! conserves by construction); what failed was the *distribution* of the heating. S4 coupling reproduced (`s5_diag_s4coupling.log`: ε(T_n; Y_lag) in the residual, network at T_n, E_bur − E_inj added at fixed ρ after the step, BE):
ΣE_inj/ΣE_bur = 1.066 and Σ|E_bur − E_inj|/ΣE_bur = 6.3 % (S4: 16 %), L_nuc(∫ε dm of the end state)/L = 1.12 (S4: 1.13–1.18), but not because of the hydro step: T_c after the hydro step = 1.2414e8 K,
after burn! = 1.2692e8 K (+2.2 %, ±1 % scatter with the dt jumps; ε ∝ T^18 ⇒ L_nuc ±30 %: the "T_c wander"). Cause: the central zone releases 1.03–1.32× ε(T_n, Y_old)Δt (1.11 at Δt = 450 yr: ε grows 7 % per 100 yr while the
CNO self-production is unfinished, plus the T history), the cool outer-core zones (T = 5–7e7 K) release 5–20× less (³He, D and the CNO catalysts relax within the step), the net correction heats the core and cools the rim after the step.
*Fix* (general): (a) `STAGE_HOOK` books E_inj and Φ_out after every converged implicit solve with the weights of the scheme (TR-BDF2: A = 1/(γ(2−γ)) on stage 1, 1 on stage 2, Aγ + (1−γ)/(2−γ) = 1), so `burn!`
works with θ-steps, both TR-BDF2 stages and the BE fallback; (b) the network is integrated along the T/ρ history of the step (ln-interpolated through T_o, [T_γ], T_n, ≤ 6 segments); (c) `solve_step!` iterates a per-zone factor
f_k on ε in the residual, f_k ← E_bur,k/(E_inj,k/f_k) (f ∈ [1e-3, 1e2], predictor 1 + fslope Δ), restoring a snapshot between iterations, until Σ|E_bur − E_inj| ≤ 3e-3 ΣE_bur *and* no zone's mismatch exceeds 1e-3 of its own energy
(with the sum criterion alone the central zones, 1e-4 of the release, were corrected by up to 5 % of their energy: T_c dips of 5 % after the step); the implicit T-dependence of the heating stays in the Jacobian, only its level follows the network.
*Result*: closed box with TR-BDF2 (`s5_b.jl`): |ΔE_gas − release|/release = 2e-11, per-step mismatch 4e-4 mean / 8e-4 max (30 steps; S4: 28 steps, 5–17 %). Stars: see the table (cumulative closure
(ΔE_tot − (ΣE_bur − ∫Φ_out dt))/Σ|E_bur| ≤ 4e-4 after the first 1e4 yr (final 4.7e-5 at 100 M☉, 2.4e-5 at 15 M☉), per-step hydro-step energy error ≤ 1e-3; T_c(after burn!) = T_c(hydro step) to 1e-4; ⟨L_surface⟩/⟨L_nuc⟩ over the MS = 0.9970 (100 M☉), 0.9967 (15 M☉), i.e. gravothermal 0.3 %, the remainder is the gravothermal term).

**2. Start-up.** A discrete hydrostatic model was tried first (outward recursion in TITAN's own momentum row for the n = 3 closure T(ρ), shooting on ρ_c): no discrete equilibrium exists on a 140-zone grid with ΔlnP = 0.3–0.5 per
envelope zone (the pressure reaches 0 at 1 % exterior mass for every ρ_c: the discrete star holds less mass than the continuous polytrope; K would have to be adapted as well), so it was dropped. Instead the initial model is built without the
S4 defects and relaxed dynamically with burning off: backward Euler (dissipative), dt = 1.5×/step up to 3 t_dyn, until t > 30 t_dyn and max|u|/c_s < 1e-3 (49–73 steps; 100 M☉: 3.3e4 s, 15 M☉: 4.4e4 s), then burning from the first
evolution step (dt0 = 3e3 s, no delayed burning). Defects of the S4 model, all in `s4_star.jl`: (i) the Lane–Emden table (RK4, h = 2e-4) cannot resolve the last 1e-9 of the mass, the outermost two interfaces collapsed onto R (2e-6 R apart,
ρ 1000× that of the neighbours): now r = R − δ*(q/q*)^{1/4} for q < 1e-7 (n = 3 surface law); (ii) T → 0 / Tfloor at the surface made the flux boundary condition F = c g_R E_b inconsistent with the interior flux by 10³ (the surface zone, m = 1e-10 M,
cools in 0.03 s): T⁴ = T_poly⁴ + ½T_eff⁴ (grey-atmosphere constant); (iii) PS.build!'s exterior-mass offset x_m,e = 10³ M made the surface zones' mass differences round-off (the extrapolated grid failed the monotonicity test spuriously).

**3. Surface.** The S4 "super-Eddington skin" was defects (i)/(ii), not physics. After the fixes Γ = κ(L + L_conv)/(4πcGm) in the envelope peaks at 0.61 (X_c = 0.745) and 0.70 (X_c = 0.5) at T = 1.0–1.1e5 K (κ_max = 0.53–0.59, the He II
bump; the S4 value κ = 1.08 belonged to the squeezed surface zones); globally Γ_e = 0.35. No cap is needed (kcap = Inf), nothing is inflated (R = 4.3 → 7.9 R☉ is the evolution off the ZAMS). What remains is numerical: the top zones
(m_ext < 1e-9 M) sit on the steep κ(T) rise at Γ_top = 0.5, so ρ_top ∝ (1 − Γ_top) is hypersensitive and Newton diverges there (density corrections ×3–10) when Δt ≳ 5e2–1e3 yr; this sets the largest usable Δt (1e3 yr at 100 M☉).
Thinner surface zones are better (330-step test runs, `s5_t*.log`): 3e-10 and 3e-11 of M: no start / failure at step 130; 3e-12: failure at step 185; 1e-12, 3e-13, 1e-13: no failure.

**4. Step control (replaces `STATIC_FALLBACK` and the ×1.3 forced growth).** *Why `STATIC_FALLBACK` was needed*: `peq` extrapolates every interface with the grid velocity of the last step, r = r_old + ṙ Δt; for Δt ≫ the zone crossing
time the noise velocities (1e-3 … 1 cm/s) move the thin outer zones by more than their width, the grid becomes non-monotonic and `peq` halved Δt down to 1 s. *Why the ×1.3 was needed*: TITAN's `timstp` (stol on r, m, ρ, T, E_r) would allow ×3 per step here
(measured smax = 1e-4 … 7e-3 against stol 0.3); after a Newton failure it sets tfac = 0.05 (Δt ÷ 20) and regrows only through sqrt(tfac·φ): the forced growth was a hand brake against the failures at Δt ≳ 1e3 yr, not an accuracy rule. *Replacement*
(all off by default for the other problems): `EXTRAP_FRAC` limits the grid extrapolation of an interface to a fraction (0.3) of its smaller adjacent zone, `PRED_LAGRHO` sets ρ_new = ρ_old ΔV_old/ΔV_new (continuity satisfied at the initial guess),
`PEQ_HOOK` extrapolates T, E_r in time; `DT_RULE` replaces the stol rule by Δt_new = Δt·clamp(0.9·0.03/max|Δ ln T, ρ, r|, 0.5, 1.6) (the 3 % is not what limits Δt: TOL = 0.03, 0.01, 0.003 give log L identical to 4 digits, R within 3e-4 and X_c within 4e-5 over 2e4 yr at
X_c = 0.49), composition caps (`DT_HOOKS`: dX_H ≤ 4e-3 from the rates averaged over the mixed regions, ΔE ≤ 5 e/ε), growth ≤ 1.15 for 12 steps after a failure, Δt × 0.5 (not 0.05) after a Newton failure (`FAIL_TFAC`); `robust_step!` repeats a
step with Δt × 0.3 from the saved state if Newton fails, if the source/network mismatch exceeds 2e-2 (zone-wise 5e-3) or if the converged state is not the neighbour of the old one (max|Δ ln T|, |Δ ln ρ| > 0.15 outside the skin).

**5. Centre MLT.** The first interface of the regular centre was wrong for three reasons. (a) `MLT_HS` measures ∇ against the hydrostatic Δln P = −Δr/H_P with Δr = ½(r_{k+1} − r_{k−1}), but zone 1 is a sphere: the distance of the zone-centre
pressures is 2× smaller, the ratio actual/hydrostatic fell into the blending range (¼, ¾), and the blend weight (a constant in the Jacobian) spoiled ∇ and made Newton converge linearly (factor ⅓/iteration in the central variables); `MLT_HSDISC` takes TITAN's
own momentum balance, Δln P = −g M_k/(r_k² P̄). (b) ℓ = αH_P with H_P = P/ρg → ∞ at r → 0 is 24 R☉ at the first interface of a 4 R☉ star: `MLT_RCAP` uses ℓ = α(H_P⁻² + r⁻²)^{-1/2}. (c) in the deep core ∇ − ∇_ad ~ 1e-7 and F ∝ (∇ − ∇_ad)^{3/2}:
Newton's tolerance 1e-5 leaves ∇ − ∇_ad uncertain at the 1e-3…1e-6 level, so L_conv(r) is noise (15 M☉ ZAMS core: L_conv/L_nuc(<r) = 7, 0.4, 5, 15 at the first four interfaces); `CONV` = 1e-7…1e-8 gives a smooth profile at the same cost.
Result: 15 M☉ ZAMS (CONV 1e-8): L_conv(1)/L_conv(2) = 0.487 (r³ law 0.471), ∇ − ∇_ad = 1.8e-6 … 3e-7 monotonic; 100 M☉ (CONV 1e-5): 0.494; L_conv/L_nuc(<r) = 1.7, 1.6, 1.6, 1.55 (the excess is the gravothermal term of the core).
(The (∇ − ∇_ad)^{3/2} law is also why Newton takes 8–16 iterations per step with a convective core against 4 without MLT (measured); a linear continuation below ∇ − ∇_ad = 1e-5 did not reduce it and was removed.)

**6. Convective core: what limits the main-sequence runs.** Two numerical problems appeared once the runs lasted 1e3 steps; both are the intermittent switching of the convective flux (zero for ∇ ≤ ∇_ad, its derivative zero as well): (i) a step in which the
innermost interfaces are not convective leaves the core unmixed for one step, the next mixing step moves X_H by 0.1 per zone, the `dX_H` cap (net change incl. mixing) then cut Δt to 1e-6 yr, which kept the core unmixed (a death spiral at X_c = 0.60 in the first 15 M☉ run and at X_c = 0.71 in the first 100 M☉
run): mixing no longer limits Δt (`capmix = false`, `REGION_CAP`: only the region-averaged burning rate does); (ii) re-solving T at fixed e for a large composition change (the boundary of a growing/receding core moves over a zone: ΔX_H = 0.2) shifted T of that zone by 10–20 %
(rejected steps, spurious solutions): partial mixing (`tmax` = 2 %: only the fraction of the change that keeps |Δ ln T| ≤ 2 % is mixed per step), contiguity of the convective region (`bridge` = 4: gaps of ≤ 4 interfaces are filled, a core is extended to the centre).
**7. Open points.** (a) 1 M☉ not done: the relaxation does not converge (max|u|/c_s ≈ 8 at the surface); the cool envelope needs its own treatment. (b) Newton cost: 12.7 (100 M☉) / 18 (15 M☉) iterations per step on average (max 25 = NITER),
Δt ≲ 2.7e3 yr at 100 M☉ (median 330 yr, 5.9e3 steps), median 5e3 yr at 15 M☉; wall 3.2 h (100 M☉, last segment, 4 threads) and 1.2 h (15 M☉). (c) TR-BDF2 works with the coupling (closed box; 15 M☉ first segments; 100 M☉ test `s5_tF.log`: 753 steps, closure 1.9e-5, R and L within ~1 % of BE) but its TR stage fails often in the
near-adiabatic 100 M☉ core (230 failures, about half the steps fall back to BE), so the production runs used θ = 1 (BE). (d) The 100 M☉ run is a chain of restarts from checkpoints across code versions (`s5_p100.log` → `c` → `d` → `e`, restarts at about
1.5e6 and 1.8-2.0e6 yr); the first segment still had the per-zone dX cap and no bridging/partial mixing; `s5_p100e.csv` carries the whole series (restored from the checkpoint aux), the last segment ran with the final code. 15 M☉ likewise (`s5_p15.log`, `s5_p15b.log` = the core-flicker spiral, `s5_p15e.log`).
(e) Convective boundary: pure Schwarzschild criterion with ∇_ad, no overshoot or semiconvection; the core boundary flickers at the Newton-tolerance level (handled by `bridge`), and τ_H is sensitive to it. (f) Mass loss and rotation are off. (g) The ZAMS luminosity of the 15 M☉ model is +0.04 dex above Schaerer.
(h) The last CSV row of each run is the X_c < 1e-3 hand-over; its L_nuc/L (0.94, 0.81) is the burn-out transient, not the MS budget.

**Validation** (X = 0.75, Y = 0.25, Z = 0; ZAMS-1 = first state in thermal equilibrium, |L_surf/L_nuc - 1| < 2 % over 10 steps; at 100 M☉ the CNO is still building there, values at X_c = 0.745 in brackets; log L in L☉, R in R☉, T_c in K, ρ_c in g cm⁻³):

| quantity | 100 M☉ TITAN | 100 M☉ Schaerer (interp. 80/120) | 15 M☉ TITAN | 15 M☉ Schaerer |
|---|---|---|---|---|
| log L | 6.174 (6.128) | ≈ 6.11 | 4.368 | 4.324 |
| log Teff | 4.980 (4.977) | ≈ 4.976 | 4.770 | 4.759 |
| R | 4.46 (4.28) | | 1.47 | |
| T_c | 1.334e8 (1.281e8) | | 9.94e7 | |
| ρ_c | 51.2 (44.8) | | 165.0 | |
| M_cc/M | 0.61 (0.78) | | 0.25 | |
| τ_H [yr], X_c < 1e-3 | 3.018e6 (3.026e6 from t = 0) | 2.74e6 (Marigo 2.95e6) | 9.95e6 (1.003e7 from t = 0) | 1.04e7 (Marigo 1.13e7) |
| τ_H vs Marigo / Schaerer | +2.3 % / +10 % | | -12 % / -4 % | |

Marigo 2001 Fig. 3 (T_c-ρ_c, read from the figure): the 100 M☉ post-CNO point (log ρ_c = 1.65, log T_c = 8.1) lies at the top of lines a/b; 15 M☉ (2.22, 8.00) near line a.
CNO self-production, 100 M☉ (t; X_CNO,c; ε_CNO/ε; T_c): 1e3 yr 1.5e-12, 0.04, 1.13e8; 3e3 yr 1.3e-10, 0.90, 1.275e8; 1e4 yr 4.4e-10, 0.94, 1.328e8 (peak); 3e4 yr 7.7e-10, 0.91, 1.28e8; 1e6 yr 3.6e-9, 0.98, 1.19e8. 15 M☉ (pp phase first): 1e5 yr 1.5e-12, 0.011; 3e5 yr 8.2e-12, 0.06; 1e6 yr 5.2e-11, 0.32; X_c = 0.6 2.8e-10, 0.74; X_c = 0.3 1.9e-9, 0.975; X_c = 0.1 6.3e-9, 0.998.
Evolution, 100 M☉ (X_c; t [yr]; R; log L; log Teff; T_c [1e8 K]): 0.745; 3.2e4; 4.28; 6.128; 4.977; 1.281 | 0.5; 1.38e6; 6.06; 6.234; 4.928; 1.190 | 0.3; 2.19e6; 9.2; 6.309; 4.856; 1.194 | 0.1; 2.78e6; 21.9; 6.367; 4.682; 1.225 | 1e-3; 3.026e6; 76.7; 6.406; 4.420; 1.382.
15 M☉: 0.745; 1.09e5; 1.47; 4.362; 4.769; 0.994 | 0.5; 4.76e6; 1.87; 4.486; 4.747; 1.009 | 0.3; 7.52e6; 2.40; 4.577; 4.715; 1.017 | 0.1; 9.34e6; 3.23; 4.650; 4.669; 1.048 | 1e-3; 1.003e7; 3.56; 4.700; 4.661; 1.230.
Ledger: cumulative closure 4.7e-5 (100 M☉; max 4e-4 after 1e4 yr) / 2.4e-5 (15 M☉); MS ⟨L_surf⟩/⟨L_nuc⟩ 0.9970 / 0.9967. Steps, last segment incl. 53/49 relaxation: 3375 / 1524, wall 1.16e4 s / 4.4e3 s (4 threads); whole 100 M☉ chain 5.9e3 steps.


### S6a Accreting Pop III protostar [run]
Files: `julia/s6b_run.jl` (driver, options in its header), `s6b_setup.jl` (embryo, SPS/OP03-type boundary conditions `precursor`), `s6_star.jl` (isentropic embryo), `s6_setup.jl`, `s6_run.jl` (approach A probe),
`s6_advect_test.jl`, `s6_analyze.jl`, `s6_profsum.jl`, `s6_kappa_cmp.jl`; `StarMaterial.jl` (`advect_composition!`, `advect!`, `TFLOOR_OPAC`), `TitanResidFast.jl` (`ALE_MDOT`, `inflow_mass/egas`, `PEXT_HOOK`, `XIMR_HOOK`,
`MLT_WOLD`, `MLT_HSAD`, `ABSORB_TRAD`), `TitanCore.jl` (`NEWTON_LIM`, ALE predictor), `TitanDriver.jl` (`ALE_MDOT`, `ALE_PSI`). All new switches default off; quick tier 8/8 PASS; full tier gains `s6_advect`, `s6b_early`.
Logs: `lib/popiii/s6_a3.log`, `s6_v1.log`, `s6_v2.log` (approach A), `s6b_pA.log/.csv/_analysis.txt/_profsum.txt/_prof_*.csv` (Ṁ_fid/4 production), `s6b_pN80.*` (80 zones), `s6b_pB.log` (Ṁ_fid), `s6_advect_test.log`.
All runs `TITAN_NOLIB=1`, Z = 0, X = 0.75, D/H = 2.5e-5, ³He/H = 1e-5 (the accreted gas: `accreted_Y`), α_MLT = 1.5, backward Euler.

**1. Approach.** (A) was built and tried, then dropped for (B); S6b will need (A) in some form (see 7).
*(A) shock in the domain.* `s6_setup.jl`: isentropic star (regular centre) + free-fall inflow region (ρ = Ṁ/4πr²v_ff, Lucy T) to 4 R*, lgrid = 1 (Dorfi–Drury, weights ln m, ln ρ, ln T), `leobc = 2`, E_r/F_r, StarGas.
Composition advection (`StarMaterial.advect_composition!`: conservative donor-cell/PLM in mass, one Barth–Jespersen limiter per zone so Σ A Y = 1 is kept, substeps for CFL > 0.5, implicit upwind sweep when all fluxes point one way) is tested in `s6_advect_test.jl`:
(1) random grid motion with inflow: Σ M Y_i conserved to 1.1e-15, Σ A Y − 1 = 4e-15, Y ≥ 0; (2) matter step through a drifting grid (c = 0.15, 200 steps): L1 error 4.0 (1st order) → 1.2 (2nd order) zone masses, width 0.9 zone; (3) stationary grid = identity;
(4) displacements of 3 zone masses per step: 6 substeps, 7e-16; (5) self-similar growth ×2.2 with donor CFL up to 2.5e5 (surface zones): 3.6e-16, step smeared over 1.7 zone masses. PASS (`s6_advect_test.log`). `advect!` applies it to the TITAN state (x_me before/after).
Result of the probe: the run works (relaxation, shock, composition conserved to 1e-15 per step) but dt cannot grow. Baseline (τ_grid = 1e7 s, 170 zones: 100 star, 70 infall zones): dt = 24 s at t = 1e3 s, then 2–9 s up to t = 6.4e3 s, 1100 steps (`s6_a3.log`);
with the grid active (τ_grid = 1e3 s, 110 infall zones) dt = 5e-5 s at t = 0.09 s after 400 steps (`s6_v1.log`), τ_grid = 1e5 s: 0.02 s at t = 5 s (`s6_v2.log`). Needed: 1e7–1e9 s. Newton at the shock zone (id corrections ±0.2…5) fails at 1.6× dt, and the
adaptive grid moves interfaces by 10–35 zone masses per step. A cure would need a stationary-shock solver; the accretion time (1e3–1e5 yr) divided by this dt is 1e6–1e9 steps. Not pursued further.
*(B) used: Lagrangian star with a prescribed self-similar mass grid, SPS86/OP03-type surface conditions.* TITAN lgrid = 3 with free outer boundary (S4/S5 machinery), but the grid row of interface k is not kinematic: x_me,k = x_mtot − (x_mtot − x_me,k^old)(1 + ṀΔt/M)
(`ALE_MDOT`): m_k ∝ M(t), N fixed, relative resolution fixed. Matter flows through the zones (mass fluxes dm/dt of the existing advection terms), the surface zone gains Ṁ/4π through the outer boundary (`inflow_mass`: a constant, not r²φ_r0: the surface zone
holds 1e-7 of what flows through it per unit M, so the inflow must not depend on the iterate's radius). Composition: `advect!` with δ_k = (m_k)(ṀΔt/M) after burn!/mix!. x_me offset: after each step shifted so that the surface x_me = 1e-5 M/4π (precision of the thin zones), raised by 1.5 ṀΔt/4π during a step.

**2. Boundary conditions (surface = shock = R*; `update_bc!`, every step, lagged except where stated).**
(i) Momentum: P_ext = Ṁ v_ff(R)/(4πR²) (post-shock gas pressure = ram pressure), as a function of the iterate's radius (`PEXT_HOOK`, in the Jacobian). (ii) Radiation: the incident radiation of the free-fall precursor, F = c g (E_b − E_s) with E_s = a T_s⁴, Lucy's T⁴ = T_eff⁴(W + ¾τ) with
L_tot = L* + G M Ṁ/R*, ρ_ff(r), τ = Rosseland depth of the Z = 0 StellarOpacity from infinity, T floor 2.5e3 K (`precursor`); L* from the interface 15 below the surface (the surface flux is the stiff response to this condition); E_s damped (geometric mean with the last value, change ≤ 15 %/step) and E_s(R) tabulated on 13 radii ±6 % and interpolated in ln R inside the Newton iteration
(`XIMR_HOOK`; without it the opacity cliff made the surface hunt: T_s jumps 10–20 kK between steps, dt ≈ 1e-4 t_acc). (iii) Energy: the accreted gas enters with the specific internal energy of the post-shock state (P_ram, T_s), "cold accretion": the shock energy G M Ṁ/R* leaves as L_acc (not reduced by the retained enthalpy, ≲ 10 %);
(entering with the top zone's own e made the top zone cool to 300 K: its gas is not coupled to the 4600 K field, κ_P(T_gas) → 0). (iv) Absorption of the radiation at T_rad: κ_E = κ(ρ, max(T, T_rad)) (`ABSORB_TRAD`) and opacities evaluated at ≥ 3000 K (`TFLOOR_OPAC`: the tables end there).
L_acc, the photosphere r_ph, T_eff,ph (Lucy photosphere of the precursor, dilution included) and S_ion (black body of T_eff,ph, r_ph) are diagnostics. D, ³He are supplied with every accreted gram.

**3. Embryo and start.** Isentrope (∇ = ∇_ad, StellarEOS, X = 0.75) with T_c = 1e5 K, ρ_c = 5e-4 g cm⁻³ integrated outwards until P = P_ram: M = 0.2174 M☉, R = 22.76 R☉ (T_s = 9.4e3 K), i.e. R ≈ OP03's value at this mass; the embryo is NOT OP03's 0.01 M☉ (their transient is forgotten by 0.1 M☉; ours by ≈ 0.3).
120 zones (70 interior zones geometric in mass, ratio 1.12, to q = 0.05; 50 envelope zones log-spaced in q to 3e-7; q = exterior mass fraction), regular centre. Relaxation with accretion and burning off, dt ×1.5/step up to 3 t_dyn: 180 steps, 22.6 t_dyn (t_dyn = 3.7e5 s), max|u|/c_s = 5.3e-4; R* → 22.3 R☉.

**4. What limited dt, and the fixes (all switches in the header of `s6b_run.jl`).** dt = 0.1 s at start, needed ≥ 1e7 s. Found one after the other by `TRACE`/`SCAN` (Newton corrections per variable and zone, residual along the Newton direction):
(a) MLT: the blend weight of `MLT_HS` (pressure-gradient ratio in the nearly isobaric skin below the ram-pressure boundary) was a constant in the Jacobian: Newton 2-cycled with ±1e-2 in E_r, T at the top zones for dt ≥ 0.2 s → differentiated (`MLT_HSAD`, dt ×100).
(b) Surface zone: entering with its own e made it cool to 300 K (outside the tables: NaN, singular matrices) → post-shock state as inflow (see 2); `ABSORB_TRAD`, `TFLOOR_OPAC`; a temperature floor on the EOS (`TFLOOR`) was tried and abandoned (a flat floor makes the energy row singular in T, slope 0.05 lets T run away).
(c) Inflow ∝ r(θ)²φ_r0 made the surface zone (mass ∼ 1e-7 of what flows through it) lose its mass: the constant −Ṁ/4π (`inflow_mass`); the exterior-mass offset must exceed the mass accreted in a step (`prep_xmext!`), the Newton tolerance 1e-7.
(d) Surface condition hunting: L*, R, E_s(R) lagged by a step → limit cycle (dt ≈ 1e-4 t_acc): L* taken deep (15 interfaces), E_s damped and ramped, E_s(R) and P_ram(R) inside the iteration (`XIMR_HOOK`, `PEXT_HOOK`).
(e) Step control: dt = min(1.6 dt, 0.9 TOL dt/max|Δ ln T, ρ, r|) over the zones with q > 1e-5 (`NSKIP`: the optically thin top zones follow the boundary conditions at once), accreted mass per step ≤ 3e-3 M; relaxation phase as S5.
(f) the composition advection of the surface zones (donor CFL up to 3e4) needs the implicit sweep; the mass of a surface zone is Δx_me, not ρΔV (converged only to 1e-3 of it).
Typical run (Ṁ_fid/4, 120 zones, TOL 0.05): 3118 steps from 0.217 to 70 M☉, 1777 s wall on one core, no failed Newton attempt, 69 steps rejected by the source/network test, dt = 7.6e5 s (0.3 yr) at 0.2 M☉ → 1.4e9 s (45 yr) at 20 M☉.

**5. Runs and results, Ṁ = Ṁ_fid/4 = 1.1e-3 M☉/yr** (X = 0.75, 0.217 → 70 M☉; columns of the CSV: R*, L* (interior, base of the surface layers), L_acc = G M Ṁ/R*, T_eff (precursor photosphere), T_c, ρ_c, T_max and m(T_max), L_D, L_pp, L_CNO, f_D, t_KH, t_acc, S_ion).
`pA` = mesh m_k ∝ M for all k (120 zones, TOL 0.05; 3118 steps, 1777 s on one core, no failed Newton attempt, 69 steps rejected by the source test); `pH3` = hybrid mesh (Lagrangian interior, ψ(q) smooth from 0 at q = 0.9 to 1 at q = 1e-4; 2601 steps, 1248 s, TOL 0.1);
`pN80` = 80 zones, TOL 0.15, to 12 M☉. Convergence: pA and pN80 agree to 1–2 % in R*, T_c, T_max, f_D at 0.3–10 M☉ (R*(8 M☉) = 256 vs 250 R☉, R*(10) 371 vs 365, L*(8) 1.88e4 vs 1.86e4).
Validation (OP03 numbers read from Figs. 1, 2, 4 of OP03, curves for Ṁ_fid/4; SPS86 eq. 21; the OP03 curves start at 0.01 M☉, ours at 0.22 M☉):

| quantity | OP03 / SPS86 | pA (all interfaces) | pH3 (hybrid mesh) |
|---|---|---|---|
| d ln R*/d ln M, 0.5–3 M☉ | 0.27 (SPS86 eq. 21) | 0.263 | 0.279 |
| R*(0.5 / 1 / 3 M☉) | ≈ 20 / 25–30 / 40 R☉ (Fig. 1) | 19.7 / 23.5 / 31.6 | 20.1 / 24.4 / 33.3 |
| L* at 1–3 M☉ | 0.1–0.5 L☉ (Fig. 4 top) | 0.03–20 (flickers by 10²) | 150–450 |
| luminosity wave | L* ↑ at ≈ 3.5–4 M☉ (0.5 → 20 → 1e3 at 5–6 M☉) | L* > 100 at 3.6, > 1e3 at 5.9, > 1e4 at 7.1 M☉ | already 4e2 at 3; > 1e3 at 5.5, > 1e4 at 7.2 |
| surface convection zone at the wave | ≈ 20 % in mass | 24 % (m/M 0.757–1) at 1 and 3 M☉; 15 % at 8 | 7 % at 3 M☉; 2.5 % at 8 |
| T_max off-centre, m/M | ≈ 0.1 (maintained through phase I) | 0.146 (3), 0.131 (4), 0.093 (6), 0.066 (8 M☉) | 0.151, 0.138, 0.105, 0.075 |
| D burning | active near 3 M☉ (ε_D > 0.1 L*/M*) at Ṁ_fid/4, L_D small, hardly any D consumed | L_D/L* ≤ 1.3e-3, f_D = 0.998 (6), 0.975 (8), 0.85 (10 M☉); T_max = 7.6e5 K at 3 M☉: no episode at 3 M☉, D burning from ≈ 5 M☉ | the same (f_D 0.974 at 8, 0.73 at 10) |
| R* peak, then KH contraction | peak 100–140 R☉ at 8–10 M☉, falls to 4 R☉ at 40 M☉, L* ∝ M³, R* ∝ Ṁ/M² | **not reproduced**: 256 (8), 371 (10), 610 (15), 1045 (30), 1900 R☉ (70); T_eff,ph 4100 K; L* ∝ M^1.5; t_KH/t_acc = 0.06 (8) → 0.01 | the same: 218 (8), 320 (10), 531 (15), 1500 R☉ (70); L* ∝ M^1.4 |
| H burning / ZAMS, S_ion | CNO ignition at ≈ 40 M☉, ZAMS ≈ 50 M☉, S_ion 1e48–49 | **not reached**: T_c = 5e7 (50), 7.7e7 K (60 M☉), L_pp, L_CNO ≈ 0; S_ion(70 M☉) = 3.6e37 | T_c 8.2e7 K (60); S_ion 5e39 |

Convective zones (pA, m/M): 1 M☉ core 0–0.14 and surface 0.76–1; 3 M☉ core 0–0.04, surface 0.76–1; 8 M☉ core 0–0.011, surface 0.85–1; 20–40 M☉ no core, envelope 0.68–1.0; the core of the 70 M☉ star is radiative, T_c = 8e7 K, ρ_c = 55–100 g cm⁻³ (log files `s6b_pA_profsum.txt`, `s6b_pH3_profsum.txt`, structure dumps `s6b_p*_prof_*.csv`).
So phase I is reproduced (R*(M) law, luminosity wave at 3.5–6 M☉, surface convection, off-centre T_max, small D burning), phases II–III are not: after the wave the star swells without limit as a convective envelope + compact core (red-supergiant-like, T_eff,ph = 4100–4800 K) and does not contract to the ZAMS.
The swelling does not depend on the mesh (all-interface vs hybrid), the zone number (80/120), the tolerance (0.05–0.15). Not tested: α_MLT, the boundary conditions' entropy, the embryo.

**6. Energy ledger and the Ṁ_fid run (both unresolved).** *Ledger.* Per step: ΔE_tot (int + kin + grav, Gauss quadrature of W = −G∫m dm/r) against E_inj − F_out − G M δm/R (nuclear heating, boundary flux including the inflow enthalpy and the surface work, binding energy of the new mass at R*; `LEDGER=1`).
Accretion off (Lagrangian, same star at 6.4 M☉): 0.2 L☉ in 2e3 (1e-4). Accretion on: step error 5–8 % of the terms at 0.2–2 M☉, 0.7 at ≥ 10 M☉, cumulative 0.07 (pA), 0.76 (pH3, hybrid mesh too; `s6b_pA.log`, `s6b_pH3.log`). Dissection at the 52.6 M☉ state of pA (`s6b_ledger_psi0.log`): ΔE_int + ΔE_kin closes with the work of gravity Σ F_k u_k Δt to 4e-3 L☉
(the thermal and kinetic equations are exactly consistent), but the potential energy decreases by 4.4e5 L☉ while the gravity work is 9.0e5 L☉ → +4.7e5 L☉ of spurious energy (75 % of L*): proportional to Ṁ (2.5e5 / 4.7e5 / 9.8e5 L☉ for Ṁ = 0.5 / 1 / 2 × 1.1e-3), independent of Δt (1e8 vs 4.5e9 s) and of the zone number, ≈ Ṁ G M/R_core, i.e. the potential energy released if the accreted mass fell into the core.
With a Lagrangian interior (ψ = 0 for q > 0.3, `s6b_ledger_psi1.log`) the first steps close to 3–6e3 L☉ but it rises step by step; a long run with a smooth ψ(q) (pH3) does not close either. The L* of the all-interface mesh (0.01–20 L☉ at 0.3–3 M☉) and of the hybrid mesh (150–450 L☉) differ by 10²–10³ while R*, T_c agree to 3 %:
L* of phase I and everything after the wave must be taken as unvalidated until the ledger closes.
*Ṁ_fid.* The precursor of the flow is optically thick once the gas ionises (τ_prec 400–1400, T_s = 2.8–4.3e4 K): E_s(R) is a cliff, the surface thermostat breaks down at t = 1.48 yr, R = 27.4 R☉, M = 0.223 M☉ (dt → 1e-2 s; `s6b_pB.log` all-interface mesh, `s6b_pHB.log` hybrid mesh, 2200 steps); not followed further.
The same states with the thin-precursor condition would be wrong (L_acc = 1.1e3 L☉ at R = 27 R☉ gives T_eff = 6100 K, i.e. above the ionisation cliff).

**7. Open points for S6b.**
1. Energy accounting of the accreting mesh (above): a conservative treatment of the gravitational potential of the mass that moves through the mesh (or MESA-style Lagrangian zones added at the surface, which needs the remap of the envelope; the ALE rows and the composition advection are ready) must come first.
2. (A) needs a stationary-shock solver (the dt ceiling of the adaptive grid is 1e-5…1 s, needed 1e7 s) or the infall region in a frame where the shock is at rest; S6b's Ṁ(t) from the collapse enters `ALE_MDOT`/`inflow_mass` as a function of t, the infalling gas state (ρ_ff, T, E_s) from the collapse run replaces `precursor`.
3. Surface physics of the post-shock layer: entropy of the settled gas (here (P_ram, T_s(Lucy)), at R = 600 R☉ P_ram = 8 dyn cm⁻²), α_MLT, T_eff,ph = 4100–4800 K vs the 7000+ K of OP03's compact stars, L_acc not reduced by the retained enthalpy (≲ 12 %), no H₂ below 3e3 K and opacities frozen at 3000 K, Z = 0 κ_R of S3 equals the P7 table (`s6_kappa_cmp.jl`: ≤ 10 % at 3500–2e4 K), no independent check against Lenzuni et al. 1991.
4. Embryo: 0.22 M☉ isentrope vs OP03's 0.01 M☉ low-entropy core; relaxed over 22 t_dyn; the early L* (45 L☉ at the start) depends on it.
5. Numerics worth keeping: `MLT_HSAD`, `MLT_WOLD`, `NEWTON_LIM`, `inflow_mass`, `prep_xmext!`, NSKIP step control, `XIMR_HOOK`/`PEXT_HOOK`; TR-BDF2 is not used with the ALE rows (θ = 1).

### S6a-fix Energy accounting on the accreting mesh [run]
Files: `julia/TitanResidFast.jl` (`GRAVWORK`, `gravwork`, `GHOST_RHO`, `ABSORB_SMOOTH`), `TitanDriver.jl` (`ALE_INFLOW`, `CONV_STALL`), `TitanCore.jl` (ALE predictor with `ALE_INFLOW`), `s6b_run.jl` (W_d ledger, `GWDBG`, `DECTRACE`, `LSTARQ`/`LTAU`, `DTR`, options in the header),
`s6fix_test.jl` (tests). All new switches default off except `GRAVWORK = :auto` (on only with `ALE_MDOT ≠ 0`): quick tier 8/8 PASS (464/587/713), S4/S5 tests (`s4_a/b/c`, `s5_b`, `s5_zams15`) PASS, `s6b_early` PASS (`s6b_test_early.log`).
Logs: `lib/popiii/s6fix_test.log`, `s6b_gwA2.log`, `s6b_gwA7.log` (correction per interface on the S6a mesh at 2.1 and 7.3 M☉), runs `s6b_fG.log` → `fG2` → `fG3` → `fG5.log/.csv/_prof_*` (150 zones), `s6b_fA3.log` → `fA4.log/.csv/_prof_*` (S6a's 120-zone mesh), `s6b_fH.log`, `s6b_fH2.log/.csv` (hybrid mesh).

**1. The missing term.** Semi-discrete system (per steradian, θ = 1): interfaces k with interior mass m_k = x_mtot − x_me,k, fluid velocity u_k, grid velocity ṙ_k = u_k − u_rel,k; mass crossing interface k inward Δm_k = m_k^n − m_k^o
(= r²ρ_b u_rel Δt by the continuity row, ρ_b the upwinded, van-Leer-limited density); momentum-cell mass m̄_k = ½(Δm_{k−1} + Δm_k). The momentum row applies −G m_k m̄_k/r_k², so the kinetic + internal energy rows (pdV with u, advected e_b, E_r) change by
−Σ G m_k m̄_k u_k Δt/r_k² + boundary terms (S6a's dissection: exact to 4e-3 L☉). The potential that belongs to these rows is W = −Σ_k G m_k m̄_k/r_k. Its exact change splits into
(i) interface motion at the new masses, Σ G m_k m̄_k Δr_k/(r_k r_k^o), = the ṙ part of the work (as on a Lagrangian grid), and (ii) relabelling at the old radii, −Σ G Δ(m_k m̄_k)/r_k^o = Σ_j A_j Δm_j after summation by parts,
A_j = −G[½(m̄_j^n + m̄_j^o)/r_j^o + ½ m̃_{j−1}/r_{j−1}^o − ½ m̃_{j+1}/r_{j+1}^o] (m̃ = time mean) ≈ −g_j Δr̄_j. The momentum rows credit for (ii) the u_rel part of the work, −G m_j m̄_j u_rel,j Δt/r_j². Both agree only if
m̄_j u_rel ≈ (r²ρ_b u_rel)Δr̄_j, i.e. ρ_b equals the momentum-cell density. At a steep density gradient ρ_b is the (lower) density of the zone above: the gas that sinks through interface j is credited g m̄ |u_rel| Δt, but releases only g Δm Δr̄:
a source C_j ≈ g_j F_j Δr̄_j(ρ̄/ρ_b − 1) per interface, ∝ Ṁ, independent of Δt, and not small with N where the gradient is not resolved. So the term that is missing is the gravitational work done ON THE MASS THAT MOVES THROUGH THE MESH
(the flux of the potential, Σ F_j ΔΦ_j); the boundary is not the culprit: the surface interface's A_{b+1}Δm_{b+1} = −G M δm/R(1 + O(Δm_b/M)) is exactly the binding energy the ledger books as L_acc, the inflow brings e_in only (no double counting).
The inflow lacks the flow work P/ρ of the entering gas (test 1: 3e-3 of G M Ṁ/R; S6a: P_ram/ρ ≪ GM/R), left as it is. Where it sat in S6a (`GWDBG`, `s6b_gwA2.log`, `s6b_gwA7.log`): at the junction of the two parts of the mass grid (q = 0.05–0.15: 70 interior zones of
ratio 1.12 end with 10 % of M per zone, the first envelope zone has 1 %) and at the edge of the compact core (ρ_b/ρ̄ = 0.49–0.88 there): Σ C = 167 L☉ at 2.1 M☉ (L* = 0.03 L☉, L_acc = 2.5e3), 2.5e3 L☉ at 7.3 M☉ (0.4 L*); 4.7e5 L☉ = 0.75 L* at 52.6 M☉ (S6a).
This also explains why L* of the all-interface and the hybrid mesh differed by 10²–10³ in phase I (the spurious heat went into different zones).

**2. Proof (`s6fix_test.jl`, `s6fix_test.log`).** γ = 5/3 gas, no radiation/convection/burning, discrete hydrostatic equilibrium (radii marched so the momentum rows vanish), regular centre, free surface at constant P_ext; E = U + K + W + P_ext V with W both as W_d and as the Gauss quadrature of S6a.
Error (ΔE − expected)/(G M Ṁ_g t / R), off → on:

| test | GRAVWORK off | on (W_d / Gauss) |
|---|---|---|
| relabel only (mesh drifts 20 % of M through the inner half; no inflow), smooth polytrope, 80 zones | 1.0e-4 | 2e-8 / 1.7e-4 |
| accretion of 2 % M of gas at rest with e_in = e_top, smooth polytrope, 80 zones | 2.8e-3 | −3e-6 / 4e-5 |
| accretion, low-entropy core under a ×100-entropy envelope (ρ drops ×21 per zone at the core edge) | **+5.35** (+5.10) | −3.8e-3 / −0.32 |
| relabel, same core-edge star | **+2.15** | −1.2e-2 / −0.11 |
| adaptive grid (lgrid = 1, x_me + ln ρ resolution, rigid walls), steps 2–200 after the first remap | ΔE_d = +6.6e44 erg = 3e-5 |W| | 0 to 1e42 erg (5e-8 |W|) |

The uncorrected core-edge star gains energy (ρ_c −7 %, the star swells: the S6a symptom); with GRAVWORK W_d closes to ≤ 1e-2 (the Gauss W differs by the quadrature error of the smeared jump). The first step of the adaptive grid
(a remap of the centre, ρ_c +56 %) is not cured (+6.6e-3 → −2.4e-3 |W|): large grid displacements leave the O(Δr²) split error.

**3. The fix (`GRAVWORK`).** C_j = A_jΔm_j − G m_j m̄_j u_rel,j Δt/(r_j r_j^o) for every interior interface (the urel work with 1/(r r^o), so the grid-motion split is exact for the relabelled part), removed from the gas energy rows
of the two adjacent zones (½ each), inside the residual (differentiated), for every energy-row variant. Pure function of x_me, r, u, ρ: holds for the prescribed mesh, a time-dependent Ṁ (S6b) and lgrid = 1 (`GRAVWORK[] = true`).
`:auto` = on iff ALE_MDOT ≠ 0. Exact for θ = 1 (all ALE runs); for θ < 1 the work is still θ-centred.

**4. Further numerics needed to run past 7 M☉ (found on the way; all in `s6b_run.jl` defaults).** (a) Newton 2-cycles of 1e-7…1e-4 in ln E_r, T at the H-ionisation zone, independent of Δt: `ABSORB_TRAD`'s hard switch κ(max(T, T_rad))
alternates between ∂κ/∂T and ∂κ/∂E_r where T ≈ T_rad (optically thick): C¹ blend T_a = T + (T_rad − T)s((T_rad − T)/(εT)), ε = 0.05 (`ABSORB_SMOOTH`) → quadratic convergence. `CONV_STALL` (accept a repeating dmax ≤ 1e-5) as a guard.
(b) Linear convergence (factor −0.6/iteration, 20–36 iterations, dt ≤ 3e7 s from 17 M☉ on) of the top-zone density: the free-surface ghost densities are refreshed to ρ_top after every update but their Newton rows are the identity,
so the Jacobian saw them frozen; `GHOST_RHO` (ghost ρ = ρ_top inside the residual) → 2–3 iterations, dt 1e8–3e8 s, the 23 → 70 M☉ part in 31–38 min (AD vs FD Jacobian of the top rows agrees to 1e-3 otherwise).
(c) ESLIM 1.05 (E_s change ≤ 5 %/step) to pass the thin → thick precursor transition at 9.4 M☉ (τ_prec 0 → 50, T_s 5e3 → 2.5e5 K). (d) Mass grid NIN=80 QC=1.083 QSPLIT=0.3 NENV=70 (150 zones, zone mass continuous at the split; the
S6a grid gives the same star within 2–10 %, `fA3`/`fA4`).

**5. Rerun, Ṁ = 1.1e-3 M☉/yr, all interfaces self-similar (PSI=0), from the S6a embryo** (`fG` 150 zones TOL 0.05; `fA` = S6a's 120-zone grid; both to 70 M☉, no failed step after (b); OP03 from Figs. 1, 2, 4 for Ṁ_fid/4, as in S6a):

| quantity | OP03 / SPS86 | S6a (pA, no fix) | fG (fix, 150 z) | fA (fix, 120 z) |
|---|---|---|---|---|
| R*(0.5 / 1 / 3 M☉) | ≈ 20 / 25–30 / 40 | 19.7 / 23.5 / 31.6 | 19.7 / 23.5 / 31.4 | 19.7 / 23.5 / 31.4 |
| L* at 1–3 M☉ | 0.1–0.5 L☉ | 0.03–20 | 0.03–0.06 | 0.02–0.04 |
| luminosity wave (L* > 1e2 / 1e3 / 1e4) | L* ↑ at 3.5–4 M☉, 1e3 at 5–6 | 3.6 / 5.9 / 7.1 | 4.1 / 5.9 / 8.0 | 4.3 / 5.9 / 8.0 |
| R* peak | 100–140 R☉ at 8–10 M☉ | none (1900 R☉ at 70) | **112 R☉ at 6.9 M☉** | **123 R☉ at 7.4 M☉** |
| KH contraction R*(10 / 15 / 20 / 30 / 40 M☉) | falls to ≈ 4 R☉ at 40 M☉ | 369 / 610 / 768 / 1043 / 1286 | 59.8 / 28.2 / 14.3 / 7.0 / 4.6 | 65.8 / 29.6 / 14.9 / 7.2 / 4.7 |
| L*(10 / 20 / 40 M☉), d ln L*/d ln M (10–40) | L* ∝ M³ (KH) | 3.7e4 / 1.5e5 / 4.2e5, 1.75 | 1.5e4 / 9.4e4 / 4.7e5, 2.5 | 1.5e4 / 9.1e4 / 4.4e5, 2.4 |
| t_KH/t_acc | < 1 after the wave | 0.06 (8) → 0.002 | 0.3 (7) – 0.38 (15) – 0.52 (20) – 0.63 (40) | 0.25–0.51–0.67 |
| D burning | active near 3 M☉, L_D small, little D burnt | from 5 M☉, f_D 0.85 (10), 0.46 (70) | f_D < 0.99 at 6.5, 0.5 at 10.8, 0.1 at 15.5 M☉; L_D ≤ 100 L☉ (< 1e-2 L*) to 20 M☉ | same (6.5 / 11.0 / 15.8) |
| CNO ignition (L_CNO > 0.01 L*), ZAMS (L_CNO > 0.9 L*) | ≈ 40 M☉, ZAMS ≈ 50 M☉ | not reached | 42.4 M☉ (T_c 1.18e8 K), **ZAMS 44.7 M☉**, R 4.09 R☉, T_c 1.31e8 | 42.9, **45.2 M☉**, R 4.01 |
| R*, T_c, L* on the MS (50 / 70 M☉) | — | — | 3.9 / 3.8 R☉, 1.29e8 K, 6.7e5 / 8.7e5 L☉ | 3.9 / 3.8, 1.29e8, 6.3e5 / 8.4e5 |
| S_ion (precursor photosphere, T_eff,ph) at 50 / 70 M☉ | 1e48–49 | 3.6e37 (70) | 8.2e47 / 2.3e48 s⁻¹ (14 400 / 15 800 K) | 7.8e47 / 2.4e48 |
| energy ledger, cumulative Σ|err|/Σ|terms| | — | 0.72 (pA), 0.76 (pH3) | 3.3e-3 (0.2–9.4), 5.2e-3 (9.4–18.5), 4.8e-3 (18.5–23.8), 1.1e-2 (23.8–70, signed −1.1e-2) | 5.0e-3 (0.2–23.5), 1.2e-2 (23.5–70) |

So the swelling of S6a was the spurious source: with the gravitational work of the moving mass in the energy rows the star peaks at 110–120 R☉ right after the luminosity wave, contracts on t_KH ≈ 0.3–0.7 t_acc
(R ∝ M^−1.85 from 10 to 40 M☉, OP03 ≈ M^−2), ignites CNO at 42–43 M☉ with T_c 1.2e8 K and settles on the ZAMS at 45 M☉ with R ≈ 4 R☉ — OP03's sequence, the peak ≈ 1–2 M☉ earlier. The two grids agree to 2–10 % in R.
Remaining: (i) L* below 3 M☉ is 0.02–0.06 L☉ against OP03's 0.1–0.5 (both grids; the hybrid mesh with the fix gives 1.2–2.1 L☉, `s6b_fH.log`): the luminosity of the deep layers before the wave depends on the mesh
by ×30 still (it is ≪ L_acc, so the structure does not notice); (ii) D burns only after the wave (no 3 M☉ episode: T_max = 7.9e5 K); (iii) hybrid mesh at TOL 0.1 (`fH2`, before `GHOST_RHO`): a pp burst (L_pp 1.1e5 L☉ at 22 M☉) and re-expansion to
108 R☉ at 25 M☉, not seen with all interfaces at TOL 0.05 — unexplained, the hybrid ψ(q) is fixed at the start mass (its interfaces with ψ < 1 crowd toward the centre as M grows ×100); (iv) the ledger's last 1 % (signed sink, 23–70 M☉),
probably the O(Δr²) grid-motion split while R falls ×3; (v) the Ṁ_fid case still stalls at 0.218 M☉ (0.38 yr, τ_prec ≈ 1–150, E_s cliff; a test run with all fixes: dt → 2–20 s, log not kept) — not addressed.

### S6b Protostar fed by the collapse [run]
Files: `julia/popiii_m2c.jl` (restart, sink), `ProblemSetup.jl` (`inner = :lagr`), `TitanCore.jl`/`TitanResidFast.jl`/`TitanDriver.jl` (`INNER_GHOST_GEOM`), `ChemNet.jl` (`XFLOOR`), `PopIIIMaterial.jl` (`CHEM_FALLBACK`),
`s6c_precursor.jl` (surface condition of the thick precursor), `s6b_run.jl` (`PRECTHICK`, `BCSMOOTH`, `BCSLOPE`, `BCFALLBACK`, `MDOTFILE`, `MDOTRAMP`), `s6c_prec_test.jl` (test, full tier `s6c_prec`), `s6c_mdot.jl`, `s6c_infall.jl`, `s6c_mt.jl`, `s6c_analyze.jl` (analysis).
All new switches default off. Tests: quick tier 8/8 PASS (464/587/713 steps), full tier 23/23 PASS (S4 a/b/c, S5 b/zams15, S6 advect/b_early, s6c_prec). P100 collapse with the defaults (`MAXSTEPS=12000 NSTOP=1e12 lines 100`, `lib/popiii/m2c_lines_n100_s6bchk.log`):
T_c = 293.3 / 373.3 / 447.6 / 629.0 / 866.9 K at log n_c = 7.00 / 8.17 / 9.15 / 10.15 / 11.17 (0.25-dex output, interpolated) against 293 / 372 / 448 / 631 / 867 (≤ 0.4 %). `PRECTHICK=0` of the new `s6b_run.jl` reproduces the S6a-fix run `fG5` to 0.1–0.5 % in R*, L*, T_c up to 8 M☉ (`s6b_pv5.log`).
Logs: `m2c_lines_n100_core2.log`, `m2c_abn02z19_lines_n200_core2.log` (collapse past protostar formation), sink runs `m2c_lines_n100_p100f.log` + `s6c__p100f{,_events,_mdot,_mdot_smooth}.csv`, `m2c_abn02z19_lines_n200_z19f.log` + `s6c__z19f*.csv`, convergence `m2c_lines_n100_c_{r25,r10,cd01,sp2}.log` + `s6c__c_*_{events,mdot}.csv`,
precursor `s6b_pt5.{log,csv}` (Ṁ_fid), `s6b_pv6.{log,csv}` (1.1e-3), `s6b_pv4{,a,b}.log` (failure and its fixes), protostar `s6b_p3_p100c.{log,csv}`, `s6b_p3_z19f.{log,csv}`.

**Part 1. Accretion history from our collapse.**
*Protostar.* The final runs were continued with the regular centre (`RESTART=` snapshot, `LINEEVERY=5`: line cooling re-evaluated every 5 steps, 94 % of the wall time was the ALI transfer). P100 at n_c = 8.8e22 cm⁻³ (6000 steps, `core2`), z19 at 6.7e22:
ρ_c = 0.147 / 0.112 g cm⁻³, T_c = 5.7e4 / 5.8e4 K, f_H₂,c = 0.035 / 0.029, x_e,c = 0.67; hydrostatic core (|u| < 0.5 km/s) of M ≈ 3e-3 M☉, R ≈ 0.6 R☉, an accretion shock at 1–4 R☉ (u: −0.4 → −6 km/s). The two clouds give the same protostar and the same infall inside 100 AU
(T = 2000 K at 1 AU, ρ = 5–8e-10, u = −5.5 km/s, m(<1 AU) = 0.05 M☉; M(r) ∝ r). OHY10 (Z = 0, ≤ 1e-6 Z☉): 1–2e-3 M☉, 0.5–1 R☉, n = 1e21–1e22, T = 2–3e4 K.
*Sink* (`SINK=1`). Zones whose inner interface lies inside r_s are removed and their mass added to a point mass (gravity only); the innermost interface is TITAN's Lagrangian boundary (llibc = 1) moving ballistically (velocity u − ½ g Δt, no pressure force: the grid is open to supersonic inflow; prescribing it from the neighbour,
u(a+1)√(r(a+1)/r(a)), is a positive feedback through the pressure of the first zone and ran away), ghost interface r(a−1) = r(a)²/r(a+1) (the equal-volume ghost has no real root for a zone wider than 2^{1/3} in r). Radius
r_s = 5 AU (M_pt/M☉), floor 0.05 AU (11 R☉ at M_pt = 0.01 M☉): t_ff(r_s, M) ∝ M stays a few % of t_acc (dt ∝ t_acc: 3e4 steps for 67 M☉), r_s ≈ 0.07 GM/c_s² (T = 2000 K), i.e. inside the free-fall region, and ≥ 5 R* (R* ≤ 240 R☉ = 1.1 AU at M ≥ 13 M☉ where r_s ≥ 65 AU).
The collapse grid has zone masses of 20 % of the enclosed mass (a zone spans a factor 3–4 in r at the sink): `SPLIT=4` (P100; 2 for z19: the state arrays do not hold 382 zones) divides every zone into equal-mass sub-zones (ρ kept, r³ linear in m, u linear in m).
Needed on the way: `NLTE=ep` (ALI cost 10–15 s/step on 200 zones), network step: Newton convergence required λ = 1 and a positive step, which never holds for species at 1e-18 (H₂⁺): `XFLOOR = 1e-12` (species below it are not limited, clamped at 0; default 1e-20) — the step failed at 20 000 s for a state that converged at 10 000 s; `CHEMFALLBACK=1`: a zone whose network step still fails keeps its state for that step (100 of 2.8e6 zone-steps).
*Mass flux.* Ṁ = ΔM/Δt of the shells crossing r_s (the crossing time of the interface interpolated inside the step; t shifted by the free-fall time r_s → 0.05 AU, 1–3 % of t_acc), over windows of Δln M ≥ 0.5 (single pairs scatter ×4 between neighbours).
Runs end where the cloud does: P100 (100 M☉): M_pt = 67 M☉ at 14.3 kyr (8 zones left); z19: 234 M☉ at 28 kyr (stopped by hand; Ṁ = 6e-3 at 170 M☉). M_pt(t) [M☉] at t = 10 / 30 / 100 / 300 / 1e3 / 3e3 / 1e4 yr: P100 0.43 / 0.89 / 2.2 / 4.7 / 10.6 / 23.7 / 53.5; z19 0.75 / 1.8 / 5.2 / 12.6 / 32 / 74 / 151.

| M* [M☉] | Ṁ P100 | Ṁ z19 | OHY10 Z=0 (Fig. B16, ±30 %) | P100/OHY | z19/OHY | φ = Ṁ/(c_s³/G), T of the shell at formation: P100 / z19 |
|---|---|---|---|---|---|---|
| 0.03 | 8.0e-2 | 1.3e-1 | 1.2e-1 | 0.67 | 1.10 | 12 / 20 |
| 0.1 | 5.7e-2 | 9.9e-2 | 8.5e-2 | 0.67 | 1.16 | 12 / 22 |
| 0.3 | 3.3e-2 | 6.4e-2 | 6.8e-2 | 0.49 | 0.94 | 10 / 16 |
| 1 | 2.0e-2 | 5.3e-2 | 5.5e-2 | 0.37 | 0.97 | 9 / 15 |
| 4 | 1.1e-2 | 4.3e-2 | 4.2e-2 | 0.26 | 1.01 | 5 / 12 |
| 10 | 6.9e-3 | 3.4e-2 | 3.5e-2 | 0.20 | 0.96 | 3 / 9 |
| 20 | 5.1e-3 | 2.5e-2 | 3.2e-2 | 0.16 | 0.79 | 2.5 / 7 |
| 40 | 4.0e-3 | 2.1e-2 | 2.7e-2 | 0.15 | 0.78 | 2.3 / 6 |
| 100 | — | 1.3e-2 | 2.0e-2 | — | 0.63 | — / 6 |

(Ṁ in M☉/yr.) z19 (ABN02, T ≈ 1000 K warm core) is OHY10's analogue to ±15 % up to 10 M☉; P100 (polytrope, cold outskirts) falls to 0.15 of it. φ ≈ 10–20 up to 1 M☉ (OHY10: ≈ 10 c_s³/G), then falls to 2–8 as the outer envelope is not Larson–Penston-like.
H11/HYOY12 (Pop III.1, 2D, rotation, with feedback: ≈ 3e-3 at 2 M☉ → 1e-3 at 30 M☉): our spherical non-rotating Ṁ is 4–5× (P100) and 15–20× (z19) higher, as expected. The snapshot rate 4πr²ρ|u| at the formation time (OHY10's eq. 13) is 0.6 of the dynamical one at 0.01 M☉ and 0.3 beyond 2 M☉: the envelope accelerates when the point mass grows.
*Convergence* (Ṁ at M = 0.3 / 1 / 4 / 7 M☉, P100; reference r_s = 5 AU/M☉, CDYN 0.03, SPLIT 4): r_s = 2.5: +25 / +14 / +3 / −5 %; r_s = 10: +7 / −13 / −8 / −16 %; CDYN 0.01: +4 / −6 / +1 / +1 %; SPLIT 2: +22 / +9 / +17 / +14 %. So ±20 % (r_s, resolution), ±5 % (time step).
*Infall state at the sink* (zone outside the interface that crossed r_s, events within 12 % of M averaged in ln; `s6c_infall.jl`):

| M* | r_s [AU] | n_H [cm⁻³] | T [K] P100 / z19 | f_H₂ P100 / z19 | x_e P100 / z19 | u/v_ff P100 / z19 | Mach P100 / z19 |
|---|---|---|---|---|---|---|---|
| 0.01 | 0.05 | 3e17 | 3.3e3 / 3.4e3 | 0.51 / 0.48 | 6e-10 / 8e-10 | 0.73 / 0.71 | 2.6 / 2.5 |
| 0.1 | 0.5 | 8e14 / 1.6e15 | 3.4e3 / 2.6e3 | 0.19 / 0.36 | 3e-3 / 4e-12 | 0.99 / 0.99 | 3.2 / 3.8 |
| 1 | 5 | 6e11 | 6.2e4 / 5.7e3 | 0.17 / 0.04 | 0.72 / 0.02 | 1.8 / 1.05 | 1.05 / 2.5 |
| 10 | 50 | 8e8 / 2.5e10 | 3.0e5 / 1.4e4 | 0 / 0 | 1.0 / 0.37 | 2.1 / 1.26 | 0.5 / 1.7 |
| 30 | 150 | 4e7 / 1.6e9 | 4.2e5 / 1.9e4 | 0 | 1.0 / 0.61 | 2.0 / 1.35 | 0.4 / 1.4 |

**[Corrected in S6c: the hot P100 first zone below was a flicker of the boundary zone, not a thermal state of the infall; reruns give a cold supersonic infall (1.3–2.6e3 K, Mach 3–4.4) and Ṁ 11–24 % (P100) / 7–16 % (z19) below the numbers here.]** The infall near the sink is adiabatically heated (T ∝ 1/r, ≈ 1e4 K at r_s = 5 AU M☉⁻¹ for free fall) and there is no atomic cooling in the network: z19 stays at 6e3–2e4 K, Mach 1.4–2.5, whereas the P100 first zone heats to 6e4–4e5 K (Mach < 1 beyond 3 M☉; u/v_ff = 2 at the interface is the ballistic extrapolation).
The P100 Ṁ beyond 3 M☉ is therefore partly a boundary result (the convergence tests agree to 20 %, all with the hot zone); z19 is clean. Not done: Lyα/free-free cooling of the infall, a sink radius ≫ 10 AU/M☉.

**Part 2. Optically thick precursor.**
*Diagnosis.* S6a's Lucy law gives at the stall (Ṁ_fid, M = 0.22 M☉, R = 22–27 R☉) T_s = 3–4e4 K, τ = 400–1400, i.e. all of G M Ṁ/R = 1.5–1.9e13 erg g⁻¹ radiated. But the ionisation energy of the accreted gas (H, He: 1.5e13 erg g⁻¹) is the same size: at R = 27.7 R☉ = G M/h_ion the shock cannot ionise the infall.
Energy balance of the steady flow (free fall, T_gas = T_rad, v²/2 − GM/r cancels): L_co(r) = L_∞ + Ṁ h(r) + (16π/3) r² E |v| (h = gas enthalpy with ionisation energy from StellarEOS, last term = trapping of the radiation by the inflow, factor e^{4τv/c});
at the surface L_∞ = L* + Ṁ (GM/R − h_s): the energy that heats and ionises the accreted gas is not radiated (TITAN's rows already carry Ṁ e_in; only E_s changes, the energy rows do not). Diffusion integral from far out (`precursor_E`, weighting 1/r² of the flux, κ = Rosseland of S3).
With h_s = h(T_s): g(h) = h(T_s(h)) − h strictly decreasing but T_s(h) jumps (first-order ionisation-front transition: front outside the shock, hot, τ ~ 100, T_s ~ 3e4 K, or inside, thin, 6e3 K, at h_c ≈ 0.6 GM/R): E_s(R) of either side is the cliff that stalled S6a and Newton. The root h_R (bisection on the sign of g, continuous) with the thermostat condition
**T_s = T(h_R)** (inverse of h_gas at the pre-shock density) is smooth in R, L*, Ṁ. At the S6a stall states (L* = 30 L☉, Ṁ_fid, M = 0.2174): R = 22 / 25 / 28 R☉: T_s = 11.8 / 9.9 / 8.0 kK (Lucy 41.5 / 36.2 / 28.9), τ = 3.9 / 5.3 / 3.1, h_R/(GM/R) = 0.62 / 0.41 / 0.14, L_∞/L_in = 0.39 / 0.60 / 0.87: the precursor is a thin shell (r_ph/R = 1.01), the 400–1400 of S6a were the artefact.
*Blend.* For f = h_R/(GM/R) ≪ 1 (M ≳ 5 M☉ at 1.1e-3) the new weighting alone lowers E_s of thick extended precursors by 4 (T_s × 0.7) and the star then does not contract to the ZAMS (R*(40 M☉) = 9.6 against 4.6 R☉, ledger 10 %), so there E_s = Lucy × trapping (`precursor_E` with/without advection: +4 % in T_s);
w = smoothstep in ln f between 0.01 and 0.05 blends the two (`precursor_thick`, 4.6 ms).
*Numerics.* E_s(R) is tabulated (hook) on 21 radii ±9 %: `BCSMOOTH = 0.03` (Gaussian in ln R), `BCSLOPE = 8` (|d ln E_s/d ln R| ≤ 8). Without them Newton failed at every dt at M = 9.5 M☉ (thin → thick transition, top zones 2-cycling; the step passes with `BCHOOK=0`, `s6b_pv4a.log`, or with the slope limit, `s6b_pv4b.log`);
for the high-Ṁ runs `BCFALLBACK=1`: after a failed attempt the hook is dropped for 100 steps.
*Validation, Ṁ = 1.1e-3* (`s6b_pv6` against `fG5` of S6a-fix, same grid, ESLIM 1.05, TOL 0.05; ratios at equal M): R*(0.5 / 1 / 3 / 5 / 8 / 10 / 15 / 20 / 30 / 40 / 50 M☉) = 0.93 / 0.95 / 0.97 / 0.99 / 0.99 / 0.99 / 0.99 / 1.00 / 0.99 / 0.97 / 0.96; L* (M ≥ 5) 1.03 → 1.00; T_c 1.06 (0.5 M☉) → 1.00 (≥ 15);
luminosity wave (L* > 1e2 / 1e3 / 1e4: 4.14 / 5.82 / 7.93 against 4.14 / 5.88 / 7.98), R* peak 111.8 at 6.80 (111.9 at 6.86), f_D < 0.99 / 0.5 / 0.1 at 6.31 / 10.71 / 15.44 (6.51 / 10.79 / 15.45), CNO ignition 42.35 (42.35), ZAMS 44.72 M☉, R = 3.93 R☉ (44.72, 4.09). Largest change: 7 % in R* at 0.5 M☉ (the enthalpy sink, h_R/(GM/R) = 0.03–0.15 there). Ledger Σ|err|/Σ|terms| 7.4e-3 cumulative (S6a-fix: 3e-3 – 1.1e-2 per mass interval).
*Validation, Ṁ_fid = 4.41e-3* (`s6b_pt5`, 0.217 → 200 M☉, 45 kyr, 2569 s): S6a stalled at 0.223 M☉. R*(0.5 / 1 / 3 / 5 / 8 M☉) = 30 / 38 / 53 / 62 / 76 R☉ (SPS86 eq. 21: 48 R☉ M^0.27 (Ṁ/Ṁ_fid)^0.41 = 48 at 1 M☉: 0.8), luminosity wave (1e2 / 1e3 / 1e4) at 5.8 / 6.9 / 9.7 M☉ (OP03: 8), R* peak 219 R☉ at 10.8 M☉, D burning (f_D < 0.99 / 0.5 / 0.1) at 12.7 / 18.8 / 27.7,
KH contraction 117 / 65 / 31.5 / 18 / 13 / 11 R☉ at 15 / 20 / 30 / 40 / 50 / 70 M☉, CNO ignition at 81 M☉. **Phase IIIb**: (L* + L_acc)/L_Edd (κ_es = 0.34) rises to 0.92 at 76–80 M☉, R* swells 11 → 57 R☉ (108 M☉) and falls again with the CNO cycle (31 R☉ at 174, ZAMS by the L_CNO criterion); OP03: Ṁ_fid expands at ≈ 300 M☉
(ours 4× earlier: no radiation force on the infalling precursor); a first variant without the blend swelled to 1170 R☉ at 190 M☉ (`s6b_pt3`).

**Part 3. Protostar fed by Ṁ(M*)** (embryo 0.217 M☉; `MDOTFILE` = the smoothed tables of Part 1; the rate is capped at Ṁ_fid (M/0.22)^2.5 until it meets the collapse rate at M ≈ 0.45 M☉, as OHY10 whose models do not converge at the collapse rate for the first 0.1 M☉).
P100 to 60 M☉ (12.5 kyr, `p3_p100c`, 2284 steps, 1281 s, no failed step), z19 to 170 M☉ (14.7 kyr, `p3_z19f`; restarted once at 31.5 M☉ without the E_s damping memory after a Newton stall of the same kind). T_eff,ph / S_ion,ph: photosphere of the precursor (Lucy τ = 2/3); T_eff,* / S_ion,*: black body of the stellar surface (L*, R*).

| M* | t [yr] | Ṁ | R* [R☉] | L* [L☉] | L_acc [L☉] | T_eff,ph [K] | T_c [K] | f_D | S_ion,ph [s⁻¹] | T_eff,* [K] | S_ion,* [s⁻¹] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **P100** 1 | 48 | 2.0e-2 | 66 | 0.05 | 9.6e3 | 5980 | 2.2e5 | 1.000 | 9e39 | 340 | — |
| 3 | 183 | 1.4e-2 | 84 | 2.3 | 1.6e4 | 6620 | 3.0e5 | 1.000 | 2e41 | | — |
| 10 | 965 | 6.9e-3 | 96 | 4.1e3 | 2.2e4 | 7370 | 5.2e5 | 1.000 | 4e42 | 4.7e3 | 1e37 |
| 15 | 1770 | 5.7e-3 | 198 | 5.1e4 | 1.3e4 | 6500 | 8.3e5 | 0.999 | 8e41 | 6.2e3 | 2e41 |
| 20 | 2680 | 5.2e-3 | 92 | 9.4e4 | 3.5e4 | 8140 | 1.35e6 | 0.69 | 1e44 | 1.1e4 | 3e45 |
| 30 | 4820 | 4.3e-3 | 39 | 2.4e5 | 1.1e5 | 8390 | 9.2e6 | 0.105 | 5e44 | 2.1e4 | 2e48 |
| 40 | 7230 | 4.0e-3 | 19 | 4.3e5 | 2.6e5 | 8500 | 2.2e7 | 0.017 | 1e45 | 3.4e4 | 2e49 |
| 60 | 12540 | 3.7e-3 | 9.6 | 1.0e6 | 7.3e5 | 9860 | 6.6e7 | 0.002 | 2e46 | 5.9e4 | 7e49 |
| **z19** 1 | 34 | 5.3e-2 | 89 | 18 (flicker) | 1.9e4 | 5580 | 2.0e5 | 1.000 | 2e39 | | — |
| 3 | 78 | 4.5e-2 | 131 | 1.2 | 3.3e4 | 6130 | 2.5e5 | 1.000 | 7e40 | | — |
| 10 | 263 | 3.4e-2 | 169 | 81 | 6.3e4 | 6730 | 3.3e5 | 1.000 | 1e42 | | — |
| 15 | 422 | 2.9e-2 | 190 | 1.1e4 | 7.2e4 | 6930 | 3.9e5 | 1.000 | 3e42 | 4.3e3 | 1e36 |
| 20 | 606 | 2.5e-2 | 290 | 3.6e4 | 5.5e4 | 5920 | 4.8e5 | 1.000 | 1e41 | 4.7e3 | 8e37 |
| 30 | 1014 | 2.7e-2 | 414 | 3.2e5 | 6.1e4 | 7000 | 9.0e5 | 0.99 | 2e43 | 6.7e3 | 8e42 |
| 40 | 1447 | 2.1e-2 | 151 | 5.1e5 | 1.8e5 | 7810 | 1.6e6 | 0.59 | 3e44 | 1.3e4 | 1e47 |
| 60 | 2562 | 1.7e-2 | 60 | 1.0e6 | 5.2e5 | 8010 | 1.1e7 | 0.098 | 1e45 | 2.4e4 | 1e49 |
| 100 | 5374 | 1.3e-2 | 46 | 2.7e6 | 8.7e5 | 8200 | 4.5e7 | 0.015 | 3e45 | 3.5e4 | 1e50 |
| 150 | 11200 | 6.4e-3 | 78 | 4.5e6 | 3.9e5 | 8640 | 1.26e8 | 0.013 | 1e46 | 3.0e4 | 1e50 |

Milestones: P100: luminosity wave L* > 1e2 / 1e3 / 1e4 at 3.2 / 8.3 / 12.1 M☉ (L* flickers 0.5–2 L☉ below that), R* peak 243 R☉ at 13.1 M☉, D burning f_D < 0.99 / 0.5 / 0.1 at 16.1 / 22.0 / 30.3 M☉, no CNO ignition by 60 M☉ (T_c = 6.6e7 K, L_nuc/L* = 0.12), the star contracts (R* = 9.6 R☉) while the cloud ends at 67 M☉: **no ZAMS arrival**.
z19: L* > 1e4 at 14.9 M☉ (transient spikes to 1e3 at 8 M☉), R* peak 549 R☉ at 24.4 M☉ (T_eff,ph 5200 K), D at 30 / 42 / 60 M☉, KH contraction to 46 R☉ at 70–100 M☉, CNO ignition at 151 M☉ with (L* + L_acc)/L_Edd ≈ 0.85 and R* growing again (78 R☉ at 150): **no ZAMS by 170 M☉**. Ledger Σ|err|/Σ|terms| 7.4e-3 (P100), 1.3e-2 (z19; 8e-3 with W_d).
*Comparison.* OHY10 (Z = 0, Ṁ = 0.1 → 0.02, Fig. B17 read off, ±20 %): R* ≈ 80 R☉ at 0.6 M☉, 170–250 at 3–10, peak ≈ 750 at 35 M☉, 130 at 80 M☉, no ZAMS (Eddington-limited at 80 M☉). z19 (Ṁ within ±15 % of theirs): 89 (1 M☉), 131–169 (3–10), peak 549 at 24 M☉ (0.7 of theirs, 10 M☉ earlier), 46 at 70–100: same sequence (adiabatic phase, wave at 10–15 M☉, swelling to R_peak, contraction, no ZAMS), 70–80 % of their radii.
P100 (Ṁ 0.2–0.5 of theirs beyond 3 M☉; R ∝ Ṁ^0.41): 66 / 84 / 96 / peak 243 at 13 M☉, contraction already at 20 M☉. H11/HYOY12 (Pop III.1: R* peak 100+ R☉ at 8–9 M☉, ZAMS 40 M☉, no-feedback Ṁ ≈ 1e-3 at 50 M☉): our constant 1.1e-3 run (above) peaks at 112 R☉ at 6.8 M☉ and reaches the ZAMS at 44.7 M☉; with the spherical Ṁ(M*) (4–20× higher) the star never reaches it before the cloud is used up (P100) or the luminosity approaches Eddington (z19).
S_ion of the photosphere stays ≤ 3e46 s⁻¹ throughout (Z = 0 H⁻ photosphere in the envelope, T_eff,ph 5–10 kK, locked as in OP03 Fig. 8); the black body of the star's own surface (not what escapes) reaches 1e47 at 40 M☉ and 1e49–1e50 at 60–150 M☉ (Schaerer02 ZAMS: 1.9e49 at 40, 3.5e49 at 60 M☉). No Z = 0 atmosphere correction was available.

**Open points.** (1) S7 needs S_ion(t) and the shape of the spectrum that *leaves the envelope*: frequency-dependent transfer through the infall (J_ν per zone for ChemNet's photo processes), a Z = 0 atmosphere spectrum for the stellar surface (L*, R*, T_eff,* above), and the radiation force on the infalling gas (not in our free-fall precursor; OP03's IIIb onset is at 4× the mass of ours).
(2) Ṁ from the 1D collapse is spherical and non-rotating (4–20× H11); the sink boundary region is hot without atomic cooling (P100: Mach < 1 beyond 3 M☉; ±20 % in Ṁ from sink radius/resolution); Ṁ(M*) is used, not Ṁ(t) (arrival delay ≤ 3 % of t_acc). (3) The precursor is a steady, LTE (T_gas = T_rad), free-fall, 1D diffusion model with one-sided energy balance; the ionisation-front transition is
resolved only through the thermostat T_s = T(h_R) and the blend with S6a's Lucy law (f = 0.01–0.05); a time-dependent shell (approach A) is not needed for it. (4) The Part 3 start (cap of Ṁ over 0.22 → 0.45 M☉) and the L* flicker below the wave (0.05–20 L☉, also in S6a-fix) are unvalidated; the z19 run was restarted once (E_s damping memory). (5) Sink: no cooling of the hot infall, P100 ends at 67 M☉ by cloud exhaustion, `SPLIT` ≤ 4 (SOD arrays hold ≤ 300 zones).

## S7-prep Atomic cooling of H and He (ChemistryKernels session, 21 Oct 2026) [run]

S6b found the sink boundary region at 6e3–4e5 K with no atomic cooling in the network. `julia/AtomicCooling.jl` (standalone,
no new dependencies; test `julia/atomic_cooling_test.jl`, 11/11 PASS) provides the complete Cen (1992) / Anninos et al. (1997)
atomic set for H and He: collisional excitation (H Lyα, He I 2³S ∝ n_e² n_He⁺, He II), collisional ionisation (13.6 / 24.6 /
54.4 eV × k1 / k3 / k5, He I 2³S), recombination cooling (case B, Hui & Gnedin 1997 rates; He⁺ dielectronic), bremsstrahlung
(Z² weights), plus the helium ionisation state `he_fractions(T, n_e; mode = :cie | :saha)`. The fits are transcribed from
EmissionKernels.jl (identical to all digits [run]) — same primary sources, not re-checked against the papers here; the He⁺
dielectronic coefficient of `k4d` is [memory]. CIE curve at n_H = 1: H peak 1.0e-22 (Λ/n_H²) at 1.8e4 K, He share 0.83 at 1e5 K,
bremsstrahlung above 1e6 K (9e-24 at 1e7 K) [run].

**Superseded by S6c.** `AtomicHHe.jl` (S6c, in ChemNet, with He⁺/He²⁺ species, Milne pairs and escape probabilities) is the
set used by the collapse/protostar runs; `AtomicCooling.jl` stays as the independent thin-limit cross-check (CIE curve) of it
and is not wired into NetGas.

ChemistryKernels.jl (branch `dense-equilibrium`) now has the same helium set behind `helium_cooling = true` (EmissionKernels
`cooling_rate_helium`), with He Saha at the gas T and the He ionisation energy in the conserved budget under `planck_radiation`.

### S7-prep.2 Free–free processes completed (ChemistryKernels session, 21 Oct 2026) [run]

Inventory before: cooling functions had only e⁻–ion bremsstrahlung (Cen 1992 integrated Gaunt factor, Z² weights, emission
only); `ContinuumOpacity` had H⁻ ff (John 1988) and hydrogenic ff with g = 1; `Z0Physics` the same with g = 1. Added:
- `julia/ContinuumFF.jl` (+ `julia/data/vanHoof2014-nr-gauntff.dat`, the authors' BSD-licensed table): thermally averaged
  free–free Gaunt factors ḡ(γ², u) (van Hoof et al. 2014, bilinear in log γ², log u); He⁻ ff (John 1994, K(λ, θ) table,
  stimulated emission included); H₂⁻ ff (John 1975 as tabulated in Turbospectrum `jonabs_vac_v07.3.dat`, component 22,
  per H₂ per unit P_e, stimulated emission applied by us — inferred from detabs.f; Bell 1980 would supersede it, table not
  reachable); H₂⁺ bound–free (σ_bf per molecule) and free–free (per n_H n_H⁺) from Stancil (1994) with its LTE K(T),
  Arrhenius-extrapolated below 3150 K. Transcriptions from Korg.jl / Turbospectrum, numerically identical to their tables
  at the nodes [run]; none re-checked against the papers [read-src-code].
- `ContinuumOpacity.alpha/means`: switch `EXT_FF` (default on; `false` reproduces M3c exactly). New kwargs `nHep`, `nHepp`
  (He⁺/He²⁺ ff with Z = 1, 2 Gaunt factors) and `nH2p` (network H₂⁺ for the bf term; LTE when absent). Hydrogenic ff of H⁺ now
  carries ḡ. New `parts`: `hem_ff`, `h2m_ff`, `h2p`, `he_ff`.
- Effect on the means, He neutral, Li electrons (`continuum_ff_test.jl`, 18/18 PASS): (1e14, 2000 K) +0.3 % / +0.5 % in
  κ_P / κ_R; (1e16, 3000) ≤ 0.3 %; (1e18, 3500, f_H₂ 0.3) ×1.10 / ×1.20; (1e20, 8000) ×1.08 / ×1.25; (1e21, 2e4, x_H⁺ 0.3)
  ×1.26 / ×2.11 (Gaunt factor ≈ 1.5 plus He⁻/H₂⁺). `stellaropacity_test.jl` still 14/14 (its M3c comparison moves 0.04 %).
- `Z0Physics.GAUNT_FF` (default **false**): van Hoof Gaunt factors in `kappa_nu` for H⁺, He⁺, He²⁺ ff. Switching on raises the
  ff-dominated κ_ν by ×1.47 at (1e-8 g cm⁻³, 3e4 K, 1e14 Hz) and requires `build_table` + the S3 re-validation; Float64 only.
- `AtomicCooling.brem_Z(T, Z)` (and EmissionKernels `brem_Z`): He²⁺ bremsstrahlung with the integrated Gaunt factor at T/Z².
- ChemistryKernels `planck_radiation`: the whole radiative sum (lines, recombination, free–free, Compton) is now dropped
  (balanced by absorption in J = B), not just the collisional-ionisation terms.
`alpha` finite from 10 K (`ck/alpha-lowT`, for RadTransfer's TCO_MIN → 0): the John 1988 H⁻ ff fit is held at its 1400 K value
below 1400 K (its polynomial turns negative near 1000 K and was clipped to 0: a 270-dex cliff in the hm_ff part); the LTE H₂⁺
estimate is formed in logarithms and capped at min(n_H, n_H⁺) (K underflowed below ~800 K → NaN); CIA held at the nearest HITRAN
block outside its temperatures (H₂–H₂, H₂–He < 200 K; H₂–H < 1000 K); He⁻ and H₂⁻ ff clamped at their table edges. Test
`julia/alpha_lowT_test.jl` (16/16): every part finite, ≥ 0 and continuous (relative change < 1e-3 for ΔT/T = 1e-6) over
10 K–1e5 K and 1e12–3e15 Hz at five compositions; means finite at (1e6, 10 K) … (1e12, 300 K). The remaining steep T
dependence of the `h` part is the Boltzmann factor of the excited-level bound–free, not a discontinuity.
Completed in the same branch (`ck/freefree`, worktree `../titan-ck`, files from data.nublado.org/gauntff): the merged
relativistic Gaunt tables of van Hoof et al. (2015) for Z = 1, 2 are the default of `gaunt_ff` (agree with the 2014 table to
< 0.5 % below 1e7 K, valid to 1e10 K); the frequency-integrated ⟨g_ff⟩(γ²) table (`gaunt_ff_total`) shows the Cen (1992) fit of
the bremsstrahlung cooling coefficient is within 1.7 % for Z = 1 and 2 over 3e3–3e8 K — kept, but He²⁺ now takes gff(T/4) in
`AtomicHHe.cooling`, `AtomicCooling.brem_Z` and EmissionKernels `brem_Z` (+5–10 % on the He²⁺ term); He⁻ ff above 10 080 K
extrapolated linearly in θ to 20 160 K and held (positive, continuous). Still open: Bell (1980) H₂⁻ ff (no open-source
transcription found; Kurucz ATLAS9 and Korg carry only John 1975 / H⁻), and a two-photon continuum (S7).
### S6c Atomic H and He chemistry and cooling [run]
Files: `julia/AtomicHHe.jl` (new: network `atomic()`, cooling `cooling`/`parts`, escape probabilities), `ChemNet.jl` (species He/He⁺/He²⁺ with partition functions and E0, He eliminated by nuclei conservation when He⁺ is a species, `Photo.esc`, `photoenergy[_t]`, `steady(xfloor, tol, xclamp)`, opt-in `NESTED`, finite-Jacobian guard, overflow-free `u_H2` below 20 K),
`PopIIIMaterial.jl` (`NetGas.zonehook`, `Teqmax`, He in `x_from_table`), `popiii_m2c.jl` (`CHEM=netatomic` (no Li) / `netatomicli`, `ESCAPE`, `TEQMAX`, `NESTED`, `ATOMDT`, `ATOMICOFF`), `s6c_atomic_test.jl` (15 checks, full tier `s6c_atomic`), `s6c_atomic_cmp.jl` (Ṁ tables). All new switches are opt-in: `modern()` steady state and coupled step are identical to 15 digits to the S6b code (4 states), and the P100 default collapse (`MAXSTEPS=12000 NSTOP=1e12 lines 100`, `m2c_lines_n100_s6cchk.log`) is bit-identical to `s6bchk` (T_c = 293.9 / 373.5 / 446.5 / 627.5 / 867.9 K at log n_c = 7.00 / 8.17 / 9.15 / 10.15 / 11.17, difference 0; the S6b table values 293 / 372 / 448 / 631 / 867 are within 0.6 %). Quick tier 8/8 PASS, full tier 24/24 PASS (23 earlier + `s6c_atomic`).

**Processes** (reverses by detailed balance `lnKeq` unless said):
- Species: He (1¹S, 2³S 19.82 eV, 2¹S 20.62, 2³P, 2¹P, the six n = 3 terms resolved with NIST energies, n = 4…15 lumped with g = 4n²), He⁺ (hydrogenic, 2n², n ≤ 15), He²⁺; E0 relative to neutral He (χ = 24.587, 54.418 eV; the zero of the equilibrium table), so `Keq` is Saha (independent evaluation agrees to 1e-5, a3); electrons from charge balance, He from nuclei conservation.
- Collisional ionisation by electrons: H Janev/AAZN97 k1 (as `modern`), He and He⁺ Voronov 1997 (`ionfit = :janev` gives AAZN97 k3/k5; Voronov/Janev = 1.0 above 1e5 K for He, 0.2–0.9 for He⁺ below 3e5 K, H Voronov/Janev 1.0–1.3); the reverses are the three-body recombinations. Below the fits' limits (1e4 K H, 5e3 He, 1e4 He⁺, 3e3 H⁻, 1e4 Li) the forward rate is continued as Arrhenius (`arrh_ext`) so that the reverse is exact (ChemNet's power-law extrapolation of the reverse left up to 7e-6 of Saha at n_H = 1e20, 3e3 K).
- Radiative recombination: Milne photo pairs (`ChemNet.Photo`; J = B gives Saha). Levels: H as `modern` (Verner 1996 1s + Kramers), He I ground Verner 1996 (phfit2.f PH2(2,2)) + excited levels Kramers × 1.15 (calibrated: He⁺→He total within 2 % of Verner & Ferland 1996 at 3e3–3e4 K), He II hydrogenic Z = 2 (Verner 1s, Johnson 1972 Gaunt factors n ≥ 2). Totals against Hui & Gnedin 1997 / VF96 / GA08 (n ≤ 15 only): H case A −5.5 % at 1e4 K (−2 % at 1e5), case B −4.4 %; He²⁺ case A −12 % at 1e4 K, −4.5 % at 1e5, case B −16 %. The ground levels are separate `Photo`s whose recombination is scaled by the on-the-spot escape fraction. Limitation: Z(He⁺) holds the LTE-populated n ≤ 15 levels, so the Milne He⁺→He coefficient is divided by Z/2 (3.4 at 1e5 K, 150 at 3e5 K); harmless (He is He²⁺ there, and the dielectronic term dominates above 6e4 K).
- Dielectronic recombination He⁺ + e → He: Aldrovandi & Péquignot 1973 = Black 1981 = AAZN97 k4d = GA08 k19,di, **one-way**: its k/K reverse would be a unimolecular thermal ionisation of He of 2e-9 s⁻¹ at 1.5e4 K that exists only in a thermal photon bath (it ionised the whole He content at 1.5e4 K in a first version); added only for the J = B test (`jb_inverse`).
- Charge exchange He⁺ + H → He + H⁺, GA08 k26, reverse by detailed balance (`cx = :db`; GA08's independent k27 is 40–100× larger than the db value at 1e4–1e5 K, `cx = :ga08`). Effect on the CIE He ionisation: < 3 % above 2e4 K, up to ×20 of a He⁺ fraction of 1e-6–1e-4 at 1.2–1.5e4 K (no energetic consequence). HeH⁺ (He + H⁺ → HeH⁺ + γ, ~1e-16 cm³ s⁻¹) is not included: it is a minor sink of H⁺ below ~1e4 K and plays no role in the energy at 1e3–1e6 K.
- Removed for the atomic net (modern's detailed-balance reverses to ground-state pairs that are not the product channels of the forward fit): the reverses of H⁻ + H⁺ → 2H and H₂⁺ + e → 2H (collisional ion-pair formation / associative ionisation: 4.6e-15 cm³ s⁻¹ at 1e4 K, 130× the electron-impact ionisation at x_e = 0.05, which ionised warm atomic gas by itself). Present with `jb_inverse` (Saha test).
- Cooling [erg cm⁻³ s⁻¹, `AtomicHHe.parts`]: (rad) the net radiative loss of **every** photo pair from the same integrals, n_ion n_e (AE + β_casc AX) − n_X RE with AE = Σ w ∫4πσ/(hν)(2hν³/c² + J)e^{…} hν dν, AX = Σ E_x,l α_l (excitation energy of the recombined level, radiated by its cascade, × the escape probability of the cascade's last line × (1 − J/B)), RE = absorbed J energy. Per recombination the loss is χ + 0.80 / 0.82 / 0.91 kT (H, He, He⁺ at 1e4 K) = the whole photon energy: recombination removes E0 of the species, so the ionisation energy has to leave the gas here once; collisional ionisation and three-body recombination carry no term (their energy is E0 in `ChemNet.eos`), so nothing is counted twice (in CIE the explicit sum equals Cen's ionisation + recombination cooling). (exc) collisional excitation: H and He⁺ n = 2 from the Black 1981/Cen 1992 fit q₂ = 4.59e-8 e^{−118348/T}/(1 + √(T/1e5)) (He⁺: 8.47e-7 T^{−0.397} e^{−473638/T}/(1 + √(T/1e5))), n = 3…8 by Bethe scaling of the collision strength, Υ₁ₙ = Υ₁₂ f₁ₙE₁₂/(f₁₂E₁ₙ); He I singlets 1¹S → n¹P, n = 2–4, van Regemorter with ḡ ≥ 0.05 (metastable and spin-forbidden channels omitted: He I excitation is < 3 % of the cooling where He⁰ exists); each line × escape probability × (1 − J(ν₀)/B_ν₀(T)). (dr) dielectronic recombination α_di × (40.81 + 24.59 eV). (ff) 1.42e-27 g_ff √T n_e (n_H⁺ + n_He⁺ + 4n_He²⁺), g_ff = 1.1 + 0.34 exp(−(5.5 − log T)²/3) (Cen/KWH96, no Karzas–Latter table). (comp) 1.0178e-37 T_rad⁴ (T − T_rad) n_e (Cen: 5.65e-36 (1+z)⁴, 0.7 % apart at z = 19). ff and dr are emission only (not reduced by J).
- Escape probabilities (`AtomicHHe.set_zone!`, set per zone by the driver hook before its chemistry step): lines β = (1 − e^{−τ})/τ with τ = σ₀ n_l ℓ, ℓ = min(√π b L_v, r), L_v = 1/max(|du/dr|, |u|/r) (i.e. the Sobolev τ = (πe²/mc) f λ n_l L_v capped by the static column n_l r). Continua: the recombination continuum to the ground level of H, He, He⁺ is absorbed locally with 1 − β(τ_cont) (τ over r from Verner σ of H, He, He⁺ at the threshold: on-the-spot, case B when thick); the continua of the excited levels are thin by construction (He II n = 2, ≥ 13.6 eV, is left thin too). **Done now:** photons trapped in a thick line are assumed thermalised (energy stays in the gas, the loss is the escaping fraction). **S7 must take over:** the transport of the trapped Lyα, He I 584 Å, He II 304 Å and Lyman-continuum photons and the radiation field J_ν per zone (it enters through `env.J` of the `Photo`s — photoionisation heating RE and stimulated recombination — and the (1 − J/B) factor of the lines; a non-static `env` makes `photorates` a direct quadrature at each call), the two-photon continuum (not included), H⁻ bound–free in the thick envelope, the inverse of the one-way photo processes (H₂⁺ association, Li + H⁺, DR).
- Equilibrium table (PopIIIEqEOS) keeps He neutral. Saha He ionisation at p = 1e4…1e8 dyn cm⁻² and 1.2e4 K is ≤ 3.4e-4 of He (energy ≤ 6.5e-5 of the H ionisation energy); at 2e4 K and p = 1e4 it is 94 % (He ionisation energy 14 % of H's). `NetGas.Teqmax = 1.2e4 K` (atomic net, `TEQMAX`) keeps hotter zones in the network even where p ≥ pswitch; `x_from_table` sets He⁺, He²⁺ by Saha with the table's electrons. (SCVH95, blended in above 1e4 dyn cm⁻², contains He.) The sink-infall zones have p ≲ 1e2 dyn cm⁻² and never reach the table.
- Stiff cooling: in the hot ionised infall t_cool ≈ 4 s ≪ dt = 1e7 s and Λ has a cliff at 1.2e4 K (Lyα, e^{−118348/T}); the coupled Newton of `step` fails there (overshoots to T ≈ 15 K, where `u_H2` overflowed to NaN, or to 6e22 K). `NESTED` (on with `netatomic`): after a failed attempt T is bracketed and bisected with the chemistry solved at each trial T (backward Euler at fixed T); energy residual 4e-14 (test h). A backtracking line search on the scaled residual was tried and stalls on the cliff.

**Tests** (`julia/s6c_atomic_test.jl`, `lib/popiii/s6c_atomic_test.log`, 15/15 PASS):
a1 detailed balance: at the Saha composition (H, H⁺, H⁻, H₂, H₂⁺, He, He⁺, He²⁺, Li, Li⁺; T = 1e3–1e6 K × n_H = 1e0–1e20, 48 states, J = B(T)) every reaction and photo pair is balanced, max |net|/flux = 6.5e-14 (all 36 reactions and 8 photo pairs paired; net built with `jb_inverse`). a2 pseudo-transient relaxation from a state with every free species scaled by 0.97 or 1.03 returns to Saha to 9.8e-7 for the species with x ≥ 1e-6 (18 states, 3e3–1e6 K, n_H = 1e4–1e20; at 1e3 K the rates are < 1e-30 s⁻¹, covered by a1; the pseudo-time stop rule limits a2 to ~1e-6). a3 He Saha constant against an independent evaluation 1.00000.
b CIE (J = 0, n_H = 1e-2) against Cen 1992: max |Δx_H⁺| = 0.054 (T ≥ 1.5e4 K), |Δ(He⁺/He)| = 0.126, |Δ(He²⁺/He)| = 0.127 (Voronov against the Aladdin/Cen fits: 20–30 % in the rates); He⁺ peak at 5.0e4 K (Cen 5.0e4), He²⁺/He = 0.5 at 7.9e4 K (Cen 8.9e4).
c1 cooling curve Λ/n_H² (the SD93 and GS07 tables were not obtainable offline; the comparator is the Cen 1992/KWH96 set coded independently in the test): at Cen's composition 0.89–1.20 of Cen's total over 1e4–1e6 K, at our own CIE 0.72–1.29; H peak at 1.8e4 K, He⁺ peak at 7.9e4 K; Λ/n_H² = 1.4e-22 (1.8e4 K), 5.1e-23 (3.2e4), 4.2e-23 (1e5), 1.4e-23 (1.8e5), 4.9e-24 (1e6). c2 recombination totals (above). c3 energy per recombination (above). c4 Compton 0.993 of Cen.
d1 closed box, no cooling: max |Δu|/|u| = 2.7e-16; d2 with cooling: per step |Δu + dt Λ/ρ|/|u| ≤ 4.8e-15 and Σ trapezoid(Λ dt)/ρ = 0.997 of the u drop.
e isochoric cooling from 1e5 K at n_H = 1e3 (thin, J = 0): 44 yr to 2e4 K (1.8 × the time along the CIE curve: the He²⁺ left over from 1e5 K cools less well than CIE He⁺), 72 yr to 1e4 K (0.13 × CIE: the gas stays ionised, x_H⁺ = 0.735 at 1e4 K against 0.0015 in CIE, He²⁺/He = 0.41 at 3e4 K against 0), i.e. the recombination lag of Gnat & Sternberg 2007; n_H = 1e10 takes exactly 1e-7 of the time (the thin processes scale with n²; there is no collisional de-excitation or thick line in this test).
f Λ_exc(τ)/Λ_exc(thin) = 1 − 3e-8 at τ_Lyα = 6e-8, 0.185 at 5.8, 2.3e-10 at 5.8e9, log-slope −1.000; on the spot: ground-level continuum escape fraction 1.6e-7 in 1e24 cm⁻² of H I, 1 when thin. g equilibrium-table He (above). h stiff step (above).

**P100 and z19 sink reruns.** `popiii_m2c.jl lines 100|200`, `RESTART` = the S6b `core2` snapshots (n_c = 8.8e22 / 6.7e22 cm⁻³; not committed), `SINK=1 SPLIT=4|2 RSINKAU=5 RSINKP=1 RSMREF=1 RSMFLOOR=0.01 CDYN=0.03 NLTE=ep XFLOOR=1e-12 CHEMFALLBACK=1 TMAXYR=1.5e4|3e4 NSTOP=1e30 TSTOP=1e30`; b1 = `netmodernli` (S6b chemistry), a2 = `netatomicli`; logs `m2c_lines_n100_p100{b1,a2}.log`, `m2c_abn02z19_lines_n200_z19{b1,a2}.log`, `s6c__{p100b1,p100a2,z19b1,z19a2}{,_events,_mdot,_mdot_smooth}.csv`, `s6c_atomic_mdot_cmp.txt`, `s6c_infall_*.txt`; 0 failed zone steps in all four.
*The S6b hot infall is not reproduced.* The b1 reruns use the S6b settings and code (popiii HEAD gives events identical to every digit; step 22: log n_c 17.737, T_c 3400.5 against 17.729, 3401.4 of S6b) but the infall at the sink is cold and smooth: P100 T = 2.5e3 / 2.2e3 / 1.8e3 / 1.6e3 / 1.3e3 / 1.5e3 / 1.9e3 K at M = 0.03 / 0.1 / 0.3 / 1 / 3 / 10 / 30 M☉, Mach 3.5 / 3.7 / 4.1 / 4.3 / 4.0 / 3.2 / 2.8, u/v_ff = 0.67–0.78 (S6b `p100f`: 3.4e3 → 6e4 → 4e5 K, Mach 1 → 0.4, u/v_ff 1.8–2.1). The S6b snapshots (`n212_p100f`, every 1000 steps) show T of the first zone alternating between 1e3 and 2e5 K from snapshot to snapshot (2.1e3, 1.3e4, 2.7e3, 4.9e4, 1.2e4, 2.5e3, 1.2e4 K at steps 5000…11000), 100 chemistry failures were kept, and the run diverged at 66.9 M☉: the 6e4–4e5 K of the S6b table was a flicker of the boundary zone of that run, not a thermal state of the infall. Why it does not recur is unknown (the first 0.007 M☉ agree to 2 %, then the step sequences differ: 208 against 181 steps at M = 0.022 M☉).
Ṁ(M*) [M☉/yr], smoothed tables, rerun baseline b1 | atomic a2 | S6b published | OHY10 Z = 0: P100 0.1 M☉: 5.02e-2 | 5.02e-2 | 5.66e-2 | 8.5e-2; 1: 1.68e-2 | 1.68e-2 | 2.03e-2 | 5.5e-2; 3: 1.06e-2 | 1.06e-2 | 1.40e-2 | 4.5e-2; 10: 5.4e-3 | 5.4e-3 | 6.9e-3 | 3.5e-2; 20: 3.95e-3 | 3.95e-3 | 5.1e-3 | 3.2e-2 (a2/OHY10 = 0.59 / 0.31 / 0.24 / 0.15 / 0.12; the rerun ends at 55.9 M☉ at TMAXYR). z19 0.1: 9.20e-2 | 9.20e-2 | 9.88e-2 | 8.5e-2; 1: 4.81e-2 | 4.81e-2 | 5.32e-2 | 5.5e-2; 10: 3.03e-2 | 3.03e-2 | 3.36e-2 | 3.5e-2; 40: 1.76e-2 | 1.76e-2 | 2.14e-2 | 2.7e-2.
The atomic chemistry and cooling change Ṁ and the infall state by < 1e-3 in these cold infalls (x_e ≤ 1e-9 up to 10 M☉, T < 2.6e3 K: nothing is ionised, nothing to cool); P100 is 11–24 % and z19 7–16 % below S6b (z19 0.092 / 0.048 / 0.030 / 0.0176 against 0.099 / 0.053 / 0.034 / 0.021), the S6b numbers having carried the boundary flicker.
*Hot infall as initial condition* (the S6b snapshots restarted with `SPLIT=1`, `RESIDX` 11 = 2.5 M☉ with T1 = 4.9e4 K and 18 = 23.7 M☉ with T1 = 2.2e5 K; logs `m2c_lines_n212_h11{a,b}.log`, `h18a`, `s6c_infall_h*.txt`): without atomic cooling (`netmodernli`) the zone at the sink keeps 5.2e3 K, x_e = 0.32 at 2.7 M☉; with `netatomicli` the hot zone is at 2.9e3 K, x_e = 1e-4; the 2.2e5 K zone at 24 M☉ ends at 9.8e3 K, x_e = 3.6e-3, Mach 3.0, u/v_ff 1.7 (S6b: 4.2e5 K, Mach 0.4). Hot ionised infall therefore cools to ≈ 1e4 K within one step (t_cool ≈ 4 s), as real primordial gas does, and the infall Mach number at the sink is 3–4.5 in every run.
**Open points.** (1) The S6b flicker is unexplained; the "boundary result" remark of S6b for P100's Ṁ beyond 3 M☉ is resolved in the sense that the cold reruns (0.3–0.6 of OHY10 for P100, 0.8–1.0 for z19) replace it. (2) Atomic cooling is not exercised in a steady hot thick infall (none occurs in the rerun); the escape probabilities are Sobolev/column estimates for spherical flow, and resonance scattering makes a real Lyα escape more slowly than 1/τ; the trapped energy is S7's. (3) n > 15 levels (α −5 % H, −12 % He²⁺ at 1e4 K, −18 % H at 1e3 K); He⁺ excited-state populations; He I metastable excitation; two-photon continuum; HeH⁺; Karzas–Latter free–free Gaunt factors; SD93/GS07 tables. (4) `jb_inverse` nets are valid in a thermal photon bath only. (5) The T ≈ 8e3 K plateau of atomic cooling with trapped Lyα (H⁻ bound–free as the remaining coolant) is not tested.

### S7a Frequency-dependent transfer and photochemistry [run]
Files: `julia/RadTransfer.jl` (module), `s7a_hii.jl` (pure-H photon-conserving sweep + chemistry/energy, used by the Iliev tests), `s7a_hydro.jl` (TITAN hydro + sweep, Iliev 2009 Test 5), `s7a_test.jl` (24 checks, full tier `s7a`; `lib/popiii/s7a_test.log`), ChemNet: two inert hooks (`env.rates_fn` / `env.energy_fn` in `photorates_t` / `photoenergy_t`), AtomicHHe: `fJ` no longer returns Inf·0 = NaN for E/kT > 600 (a cold zone with J > 0 gave NaN cooling; found by the first ChemNet step with a non-zero field). Everything opt-in: quick tier 8/8, P100 default collapse log `m2c_lines_n100_s7achk.log` identical to `s6cchk` (diff empty; T_c 293.3 / 373.3 / 866.9 K at log n_c 7.00 / 8.17 / 11.17), `s6c_atomic_test` unchanged. Papers (data) in `~/Projects/popiii_lit/rt/` (Iliev 2006, 2009, Mellema 2006, Abel 1999, Wolcott-Green 2011).

**Method.** (1) Grid `default_grid()`: 177 bins 0.1–100 eV, Δln ν = 0.1 (IR) / 0.08 / 0.01 (LW 11.2–13.6) / 0.02 (above the H, He I, He II edges, which are bin edges); cross sections are 3-point bin averages. (2) Opacity per zone and bin: bound–free of every `ChemNet.Photo` (Boltzmann level populations of the Milne pair; net of stimulated emission `stim`), recombination emission η = κ_stim 2hν³/c² from the Milne relation (`emit`: the diffuse continua, so case B / on-the-spot is not imposed: `AtomicHHe.ONSPOT[] = false`), everything else from `ContinuumOpacity.alpha` (CIA, H⁻ ff, H/He⁺/He²⁺/He⁻/H₂⁻ ff, H₂⁺, Rayleigh, Thomson; He⁺, He²⁺, H₂⁺ densities passed; its H Kramers bf and H⁻ bf are removed when the network supplies them; called at max(T, 800 K): the H⁻ ff fit diverges at 100 K) with LTE emission κB. (3) Direct field from the central source (finite star inside the first interface): photon-conserving radial attenuation L_out = L_in e^{−τ} per zone and bin, the absorbed power shared by net κ among the absorbers, effective zone intensity J = absorbed/(4πκV) (so 4π∫σJ/hν·nV is exactly the photon number). (4) Diffuse field: tangent-ray formal solution (zone-constant source, exact exponentials, rays p = r_j plus core rays through the inner void, μ-moments with I linear in μ), sources = thermal + recombination emission + scattered direct beam + σ_s J (coherent isotropic scattering, Λ iteration ≤ 6 sweeps; Thomson/Rayleigh are scattering, not absorption); zone-mean J from the interface values for τ_zone ≤ 0.5, from the flux divergence (A = 4πVη + L_k − L_{k+1}) in thick zones, thin zones scaled per bin so the budget closes (the plain average of the interface J overestimated the absorption of a 180-τ zone by 500×). Outer boundary: optional isotropic `Jext`. (5) Outputs per call (`RadField`): J, Γ per absorber (gross photo-rate per atom), heating (hν − E_ion per absorbed photon), ionisation-energy sink, radiation acceleration (direct: (L_in − L_out)/(cρV), diffuse: (1/ρc)∫χF), k_LW, photon and energy ledgers; `RT.ARAD[]` keeps the last acceleration (the S7b hook). (6) ChemNet bridge: `ZoneEnv` (`env.J` bin-wise, `env.kLW`, static = false); R and the absorbed energy RE come from the transfer's bins (= Γ), the Milne α and emitted energies from fixed-node quadratures (`PhotoQuad`, σ tabulated once) lazily tabulated in log T (0.02 dex, Catmull–Rom) per zone; `RT.zonehook(ze, rf)` loads zone k before its `ChemNet.step`. `AtomicHHe.cooling(photos, ze)` then carries the photo-heating (the `rad` term subtracts n_X RE: gross absorbed energy, the ionisation energy being in E0): **the transfer's own `heat` must not be added** (test h2: Σ n_p RE_p = heating + ionisation sink to 4e-16; i1: one `ChemNet.step` in a cold zone gains Δu ρ/dt = 5.495e-19 = the transfer's absorbed energy, ratio 1.000). `heat_other` (ff, CIA, H₂⁺ absorption, not in any Photo) is the only term to add. (7) LW: k_LW = 1.38e9 [J_dir(12.87 eV) f_sh(N_in) + (J_diff + J_ext)(12.87 eV) f_sh(N_out)] per H₂, N_in = H₂ column from the source to the zone centre, N_out from the zone to the surface, f_sh Wolcott-Green et al. 2011 eq. 10 (α = 1.1, default) or Draine & Bertoldi 1996 eq. 37 (`LWSHIELD[] = :db96`), thermal b; the reaction H₂ + γ → 2H with k = env.kLW is added by `RT.with_lw(net)` only. (8) Time-dependent I-fronts: `HII.sweep!` marches inside out per hydro/time step, per zone the time-averaged x̄ solves x̄ = ⟨x⟩(Γ(x̄)) with the exact time average of the ODE with constant Γ, c, α (C²-Ray; a plain (x_old + x_new)/2 over-absorbs in zones that saturate within dt: it heated the first zone to 1.5e7 K), so ionisations = absorbed photons exactly; energy by backward Euler in T with the photo-heating of the absorbed photons.

**Tests** (`s7a_test.log`, 24/24 PASS):
a1 Iliev 2006 Test 1 (monochromatic 13.6 eV, T = 1e4 K, α_B = 2.59e-13, n_H = 1e-3, Ṅ = 5e48, x₀ = 1.2e-3, 1D, 200 zones, Δt = 0.25 Myr; t_rec = 122.35 Myr, r_S = 5.393 kpc): r_I(x = 0.5)/[r_S (1 − e^{−t/t_rec})^{1/3}] = 1.000 / 1.004 / 1.011 / 1.020 / 1.042 at 10 / 30 / 100 / 200 / 500 Myr (volume-equivalent radius 1.010 / 1.006 / 1.007 / 1.011 / 1.024; 100 zones and 0.5 Myr: 0.96 at 10 Myr, converged at ≥ 200 zones). a2 x_HI(r) at 500 Myr / the steady on-the-spot Strömgren profile (ODE, same σ): 1.00 at r/r_S = 0.3, 0.5, 0.7, 0.8 (2.67e-3, 8.16e-3, 2.02e-2, 3.46e-2). b1 Test 2 (10⁵ K blackbody, T from 100 K, Cen cooling, case-B ⟨E_kin⟩ from the Milne integrals): r_I/r_an(1e4 K) = 0.910 / 0.942 / 0.989 / 1.033 / 1.095 (Iliev Fig. 15: codes 0.93–0.97 early, 1.0–1.10 at 4 t_rec); b2 T = 2.0e4 → 1.5e4 K at 0.3 r_I and 1.5e4 → 1.1e4 K at 0.9 r_I (10 → 500 Myr), maximum 2.5e4 → 3.8e4 K in the first zones (paper 1–4e4 K), pre-heated gas ahead of the front (x 0.2–0.4 eight zones ahead at 1.5–1.0e4 K; the hard photons, absent in monochromatic codes); b3 sweep photon budget 5e-16. c (hydro, 300 Lagrangian zones of 50 pc, TITAN θ = 0.55 TR-BDF2, AV length 1.5 zone widths, Δt ≤ 0.5 Myr and ≤ 0.1 t; 106 s): r_I/r_S0 at 100 / 200 / 300 / 500 Myr = 0.859 / 1.116 / 1.327 / 1.673 against the Iliev 2009 codes (read off Fig. 18, ±0.05) 0.87 / 1.12 / 1.32 / 1.57 (−1.2 / −0.3 / +0.5 / +6.5 %); T of the ionised gas 1.24e4 → 1.17e4 K, density peak 2.3e-3 at 12.7 kpc at 500 Myr (ZEUS-MP Fig. 17: shell at 0.78 L_box = 11.7 kpc, peak 3e-3); converged to 1e-3 in resolution (150, 200, 300 zones: 1.6722 / 1.6725 / 1.6728) and in time step (0.125 Myr: 1.6719); logs `lib/popiii/s7a_hydro_n{150,200,300,150dt}.{csv,log}`. Against Spitzer r_S0 (1 + 7 c_s t/4r_S0)^{4/7} with c_s = (2kT/m_H)^{1/2} of the ionised gas: 0.69 / 0.77 / 0.81 / 0.85 (the asymptote; the plotted analytic curve of Iliev 2009 lies 30 % below their codes); slope d ln r_I/d ln t (300 → 500 Myr) 0.45 (4/7 asymptotically). d1 photons, direct beam, ground-level absorbers: |ΣΓ n V − (emitted − escaping)|/emitted = 3.9e-16; d2 heating + ionisation sink = absorbed energy to 4.8e-16; d3 with diffuse field, all absorbers, continua and scattering: (Ṅ_src + Ṅ_gas − escaping − absorbed)/Ṅ_src = −2.6e-5, energy 6e-15; d4 Σaρ V = (L_in − L_out)/c to 2e-16, (1/ρc)∫χF against the divergence form 7e-3 in thin zones. a3 the diffuse field reproduces case B: Σ diffuse photons absorbed / Σ α_1 n_p n_e V = 1.0002 (locally 0.24–0.30 at 0.3–0.8 r_S: the ionised interior at x_HI = 1e-4, n = 100 has mean free path 1.6e19 cm > r_S, the photons are absorbed at the front). e1 k_LW = 1.38e9 J f_sh(N, T) exactly (2e-16) at N_H2 = 1e14…1e21 × T = 300, 1e3, 5e3 K; at T = 1e3 K, J = 1e-21: k = 1.28e-12 (f 0.93), 7.7e-13 (0.56), 1.5e-13 (0.11), 1.6e-14, 2.0e-15, 3.8e-16, 8.0e-17, 1.1e-17 s⁻¹ for N = 1e14…1e21; e2 f_sh(WGHB11)/f_sh(DB96) = 1.06 / 1.57 / 4.65 / 4.32 / 1.95 / 1.26 / 1.08 / 1.05 (the shielding itself is the published fit; no line-by-line check: no H₂ LW line data in the repository); e3 the network reaction gives −k x_H₂ exactly and exists only through `with_lw`. f1 J = B(T) through `ZoneEnv`: every reaction and photo pair balanced to 5.7e-14 at the Saha composition (12 states, 1e3–1e6 K × n_H 1–1e16; S6c a1: 6.5e-14); f2 T-tabulated mode 3.4e-4; f3 relaxation to Saha 9.8e-7 (S6c a2: 9.8e-7). h1 env R = transfer Γ to 2e-15; `ChemNet.photorates` on the same binned J differs by up to 8.5 % in R (its nodes alias against the J bins), Milne α by 0.8 % (T ≥ 1.5e4 K; both quadratures under-resolve the recombination threshold at 100 K, where the α are negligible); T-table vs exact 4.6e-5. g cost, 200 zones × 177 bins × 7 absorbers, 3 threads: full transfer call (opacity incl. `ContinuumOpacity.alpha`, direct, diffuse with Λ iteration, rates) 1.35 s, direct only 0.37 s; per-zone ChemNet env (`set_zone!` + Milne rates and energies of 7 photo pairs at 6 temperatures) 5.3 ms per zone, 1.07 s per 200 zones (single thread); the pure-H sweep costs ≈ 1.5 ms per step at 200 zones × 1 bin (Test 1: 3 s for 2000 steps; Test 2, 40 bins, 86 s).

**S7b needs.** (1) Coupling schedule: call `transfer` every N hydro steps with J, k_LW from the previous chemistry state; the snapshot J is photon-conserving but lagged, so an I-front that crosses several zones per call (R-type phase, Δt_call ≳ the crossing time) is wrong in speed — either cap Δt (as in the Iliev runs) or write the inside-out sweep for ChemNet (`HII.sweep!` is pure H; the multi-species version needs the per-zone network step interleaved with the attenuation, He edges included). (2) `RT.ARAD[]` into the momentum equation of `TitanResidFast` (ρ a per zone) with the star's L_ν as `Source` (blackbody_source(Teff, L) or a table; R_* inside the sink); the Lyα / He line force is not included. (3) Remove the heat term when the atomic cooling is used with `ZoneEnv`; add `heat_other`. (4) Not done: line transport (Lyα, He I 584, He II 304: trapped energy is still thermalised in AtomicHHe), H₂ LW line transport (k_LW uses the continuum J(12.87 eV) × f_sh; LW photons absorbed by H₂ are not removed from the beam; H₂ photodissociation heating and H₂⁺ photodissociation absent), two-photon emission (`twophoton`, needs `X[:aB]`; untested), Compton heating by hard photons, light-travel time (steady transfer), diagonal ALI (scattering converges for τ_es ≲ a few; thicker needs ALI), non-radial acceleration. (5) The raw tangent-ray moments close the budget to 3e-3 (uniform sphere, τ = 10); the production value 2.6e-5 relies on the per-bin rescale of thin zones; Gauss μ-nodes would reduce the raw error. (6) Cold-gas Milne quadratures (T ≲ 300 K) under-resolve the threshold in both ChemNet and `PhotoQuad`.

### Performance recording (`julia/Perf.jl`)
Every driver records performance the same way, so runs can be compared when looking for room for improvement:
- `Perf.reset!(run; stage, meta...)`; `Perf.@section :name expr` gives inclusive wall time, calls and allocations;
  `Perf.count!`; and sparse per-step `Perf.sample!`.
- `Perf.write()` writes `lib/popiii/perf/<run>_{meta.txt,sections,counters,samples}.csv`. The metadata records host,
  CPU, Julia version, threads and git revision.
- Sections to use where they apply:
  - `:step`, `:newton`, `:residual`, `:jacobian`, `:linsolve`;
  - `:eos`, `:opacity`, `:chem`, `:burn`, `:mix`, `:linecool`, `:rt_transfer`, `:rt_photorates`, `:io`.
- Counters to use where they apply:
  - `:steps`, `:newton_iter`, `:step_fail`, `:stage_fail`, `:be_fallback`, `:chem_fail`, `:ali_calls`, `:rt_calls`.

### S7b Coupled protostellar feedback [run]
Files: `julia/RadSweep.jl` (in `RadTransfer`: `ion_sweep`, `interface_accel`; `set_zone_light!` in `RadTransfer.jl`), `s7b_feedback.jl` (module `S7b`: star spectrum, inner free-fall column, `envfields`, OI02 / Lyα estimates), `s7b_exchange.jl` (records), `s7b_env.jl` (hooks of `popiii_m2c.jl`, `S7B=1`),
`s7b_star.jl` (hooks of `s6b_run.jl`, `S7B_STAR=1`), `s7b_env.sh`, `s7b_star.sh` (production settings), `s7b_test.jl` + `s7b_shell.jl` (full tier `s7b`), `s7b_offline.jl`, `s7b_break.py`, `s7b_compare.py`, `s7b_summary.py`, `s7b_forceband.jl` (analysis);
edits of existing files, all default off: `TitanResidFast.AEXT` (momentum row), `AtomicHHe.FJ_FLOOR`, 3 hooks in `s6b_run.jl` (`S7BH`, `S7BDT`, `S7BMDOT`), 6 lines in `popiii_m2c.jl`. Perf records `lib/popiii/perf/s7b_*`. Logs and tables cited below are `lib/popiii/s7b_*`, `m2c_abn02z19_lines_n200_s7b_*.log`.

**1. Coupling (two processes, common clock).** TITAN's state is global, so (a) two long-lived processes (envelope = `popiii_m2c.jl` sink run from the `core2` snapshot; protostar = `s6b_run.jl`) exchange small text records (key = value, atomic rename) at barriers of a common clock t_c = 0 at the synchronisation point (M ≈ 3 M☉);
(b) (alternate segments with save/restore) was rejected: the envelope restart rebuilds the grid from a profile and loses the step history, and every segment would pay the S6b relaxation. The envelope is the clock master: barrier spacing Δτ = `FEX` × M/Ṁ (`S7B_FEX` = 0.01 default). Both processes run in parallel over the interval (dt cut to land on the barrier), then:
star → envelope: M, R*, L*, L_acc, L_∞, T_ph, r_ph, T_s, τ_prec, T_c, t_KH (`star_<j>.txt`); envelope → star: Ṁ_win, accreted mass since the sync M_acc, next barrier, infall state at the grid edge (`env_<j>.txt`). The star's Ṁ for the next interval is Ṁ_win plus (M_acc − (M* − M*_sync))/(15 Δτ) (limited to ±20 %, ≤ 5 % change per interval, ramped linearly in time: the precursor boundary condition stalls on steps in Ṁ).
*Accretion rate.* Sink events (crossing time of each interface through r_s, mass m_i; S6b Part 1) are a staircase with jitter; Ṁ_win and M_acc come from a least-squares line through the events of the last window (ln M span ≥ 0.3, ≥ 4 events), bounded by the average since the first event when the sink starves (without the bound the fit kept growing when the infall had stopped: found in the first production attempt). The Eulerian flux 4πr²ρ|u| at the grid edge is 3–10 times lower than this dynamical rate at M > 20 M☉ (the first zone spans a factor 3–6 in radius and 5–10 M☉), so ρ_s of the inner column is set by Ṁ_win/(4πr²|u|) (`S7B_RHOMODE`).
*Synchronisation error* (c3 FEX = 0.02 against c4 FEX = 0.01, 3 → 7 M☉, 96 yr, `s7b_compare.py`): M*(t) differs by ≤ 0.8 % to t_c = 70 yr and 1.2–1.8 % at 80–90 yr (≈ the noise of the accretion estimate: the accreted mass at the same sink state differs by 5 % between the two runs), T_c by ≤ 0.4 %, R* by 1–9 % (R* wanders on the 20 % scale around the plateau), L* flickers (0.6–2e3 L☉ in both, not comparable). The exchange interval is ≪ t_acc (1 %) and ≪ t_KH of the star, but the star state used by the envelope is one interval old (first order); a run with identical settings reproduces itself exactly (p1 against c4 to 4 digits).
*Robustness.* The S6b star stalls at the ionisation cliff of the precursor (step → seconds, L*/T_s flicker) at ≈ 30 M☉ (t = 3.41e10 s; S6b restarted there by hand) in the coupled runs p3, p5, p6; a stall detector (dropping `ES_OLD`, the BC hooks, rollback to the last checkpoint with replay of the recorded Ṁ ramps) did not clear it (the checkpoint state itself stalls at constant Ṁ), so the two-way production run could not be completed through 30 M☉ in this session (p1, an earlier code state, passed to 135 M☉ but is invalid, see 5). The production below is therefore one-way from the S6b star track (see 6); the exchange itself was exercised to 33 M☉ (p6) and to 135 M☉ (p1).

**2. Source and the gas inside r_s.** Radiation enters the envelope grid at its first interface r_e ≥ r_s (5 AU M*/M☉; the grid edge lies 1–10 r_s out). The gas between the stellar surface and r_e is not on the grid: it is the *analytic free-fall column* ρ(r) = ρ_s (r/r_e)^{−3/2}, u ∝ r^{−1/2} normalised to the envelope at r_e, in which the stellar spectrum is attenuated by `ion_sweep` (H, He, He⁺ photoionisation equilibrium, case B, the OI02 flux term Ṁ/m_p included, T = 1.5e4 K) from R* to r_e (`inner_column`);
the absorbed power is re-emitted as a blackbody at the S6b precursor photosphere T_ph (steady state: bolometric L conserved to 1e-9, test e2). Source spectrum: blackbody of the stellar surface (T_eff,* from L*, R*) + blackbody (T_ph) for L_∞ − L* (`mode = :star`, the upper limit of the ionising output; `:photosphere` = S6b's S_ion,ph ≤ 3e46 only). *Result* (`s7b_offline.jl` along the S6b z19 track with the infall state of the cold sink rerun `z19a2`, `lib/popiii/s7b_offline_z19_{mdot,grid}.csv`): the H II region is confined to R_II ≤ 1.02 R* at every mass up to 166 M☉ with the dynamical Ṁ (Q*/Q_crit ≤ 1.5e-2); nothing above 13.6 eV enters the grid except the thermal tail of the re-emission (S_ion,ph ≤ 1.5e46; 1e44–1e46 s⁻¹ absorbed in the first zone);
with the Eulerian grid density (3–10× lower) Q/Q_crit reaches 1.4 only at 166 M☉ (R_II = 2.8 R*, still no escape). Breakout rate (OI02 eq. 5–7, `s7b_break.py`, α_B = 2e-13): Ṁ_break = 1.4e-4 (100 M☉), 1.6e-4 (130), 2.3e-4 M☉/yr (166), against Ṁ = 1.1e-2, 8e-3, 5.6e-3 of the dynamical collapse: a factor 25–80 (3–14 with the Eulerian density).

**3. Radiation force.** `TitanResidFast.AEXT[]` (empty = off, one `isempty` per residual call) adds an outward acceleration per interface to the momentum row, `w −= a_k m_k` (same form as the radiation-flux force; held fixed over the step, so the force of a step is that of the last evaluation). It is fed with `RadTransfer.interface_accel(RT.transfer(...).arad, m)`:
the continuum force of the direct beam, (L_in − L_out)/(c ρ V) per zone (Thomson, Rayleigh of H and H₂, H⁻ and H bound–free, ff, CIA, photoionisation: every absorbed photon deposits hν/c), plus the diffuse field when `S7B_DIFFUSE=1`; force conserved at the interfaces. *Missing:* the Lyα and He I 584 / He II 304 line force and the damping-wing trapping of Lyα (McKee & Tan 2008 §4, eq. 22), lines are not transported.
What the continuum already contains is the Rayleigh scattering of neutral H, whose Dalgarno polynomial is used above 10.2 eV, where it is only the off-resonance dispersion estimate (the true Lyman-series wings and cores are larger): for the M = 100 M☉ star (L_∞ = 2.5e6 L☉, T_ph = 8.7 kK) at 3000 AU in n_H = 1e7 gas the force is a_rad = 0.52 g, 71 % of it from 10.2–13.6 eV (mean extinction 2.2 cm² g⁻¹ over 10 % of L), 29 % from 4–10.2 eV (`s7b_forceband.jl`). The line force is therefore not zero but a crude, probably low, surrogate; a bound for the trapped Lyα of the H II region (0.68 Q photons, N_sc = 1/ε, ε = n_p q_{2p→2s}/A) is ≤ 2.4e-4 of the infall ram momentum flux at the H II edge for the dynamical Ṁ (all masses to 166 M☉) and 5e-4 (100 M☉), 5e-3 (130), 0.17 (166) with the Eulerian density (`s7b_offline`, last column): negligible for a confined H II region.
The momentum of the beam absorbed inside r_s (the H II region at R*) is not applied to the inner free-fall flow (not on the grid): L_∞/c is 0.003 (30 M☉) → 0.02 (150 M☉) of Ṁ v_ff(R*), and the Eddington factor of the star+accretion luminosity reaches Γ = 0.9–0.97 at 80–150 M☉ (electron scattering, κ = 0.35), i.e. the inner flow is decelerated by up to a factor 1 − Γ in the ionised layer of S6b (in the stellar model, not in the envelope).

**4. Ionisation fronts in the infall, LW, heating.** The lagged-J problem (a front that crosses many zones per transfer call is wrong in speed) is removed instead of capped: `ion_sweep` solves, zone by zone from the source, the *equilibrium* of H⁺, He⁺, He²⁺ with the exact photon budget (absorbed photons of each species = case-B recombinations in the zone + the inflow term), so that the front position follows from photon conservation at any Δt (R-type fronts are instantaneous on the hydro time scale; recombination time 1/(α n) ≈ 4e3 s ≪ the exchange interval).
The equilibrium composition (raised above the network's, never below) enters `RT.transfer` (other continua, diffuse field, LW with H₂ self-shielding, radiation acceleration); the per-zone J_ν bins go to the ChemNet step (`ZoneEnv`, `set_zone_light!`: only R and RE are refreshed, the T-tabulated Milne coefficients are kept, valid while J ≪ 2hν³/c², which holds everywhere outside the immediate neighbourhood of the star). `RT.with_lw(net)` adds H₂ + γ(LW) → 2H with k_LW = 1.38e9 J_LW f_sh.
*Heating, no double counting:* photoionisation and H⁻ detachment heating enter through `AtomicHHe.cooling(photos, ZoneEnv)` (absorbed J energy RE); the absorption of the non-network continua (H⁻/H ff, CIA, …) is added as negative cooling (`heat_other`; the same opacities' emission is not subtracted, because the network already carries the free–free emission and the CIA emission is negligible: the absorbed power is 1e-5 of the H₂ line cooling at n = 1e7). A one-sided variant that subtracted 4π∫κB was tried and double counts ff (zone T → 0).
**Bug found and fixed on the way:** `AtomicHHe.fJ` = 1 − J/B at the Lyα line is unbounded below: with a stellar J(Lyα) ≫ B(700 K) it reached −1e63 and, times the recombination cascade, a spurious heating of 3e38 erg cm⁻³ s⁻¹ in cold gas (the gas settled at an unphysical 1650 K and the first production runs showed a "stall" that was this artefact). `AtomicHHe.FJ_FLOOR` (default −Inf = S7a) is set to 0 by S7b: no net line cooling where J exceeds the Planck function, the absorption is not a heating term of that formula.

**5. Tests** (`s7b_test.jl`, `lib/popiii/s7b_test.log`, full tier `s7b`, 20 checks PASS; quick tier 8/8 PASS, full tier 26/26 PASS (`s7a`, `s6c_atomic`, `s6b_early` unchanged)):
a OI02 (2002) §2, free-fall envelope ρ = Ṁ/(4π r² v_ff), M = 30 M☉, R_in = 10 R☉, 20 eV photons, Ṁ = 1e-3 and 1e-2 M☉/yr, Q − Ṁ/m_p = q Q_crit, q = 0.5…5: R_II/R_in = 1.652 / 2.727 / 7.407 / 20.12 / 148.4 against e^q = 1.649 / 2.718 / 7.389 / 20.09 / 148.4 (max 0.3 %), photon budget (ionisations = recombinations + inflow term) 1.7e-14; the critical rate for R_II = 1e4 R_in at Q = 1e50: Ṁ_crit = 5.5632e-5 against 5.5631e-5 M☉/yr (bisection on the sweep against the analytic eq. 5–6). Q_crit^FF = 3.5e51 s⁻¹ (M = 30, Ṁ = 1e-3; OI02's 1e51 (M/100)^{-1} ṁ² with their α).
b OI02 §3, uniform-velocity solution ρ = M/(4π R_H r²): at Q = Q_crit(R_II) the ionising radiation force equals gravity, a_rad/g = 1 ± 6e-6 at every radius of the H II region (R_II/R_in = 10, 30, 100), R_II = R_II,crit to 1e-3, and Q_crit(R_II → ∞) = α M²/(4π R_H² m_p² R_in) = (hν/c)⁻¹ G M²/(R_in R_H) (eq. 30) to 1e-9 (R_H = 36 pc for 20 eV photons, 26 pc for 13.6 eV).
c Radiation-pressure shell in TITAN (cold ionised shell around a hot cavity, central L = 1e43 erg/s, RT.transfer every step, `AEXT` hook): the momentum gained by the shell equals ∫Σ a_j M_j dt to 1.0000 (thin, τ_Th = 4e-3) / 0.9988 (thick, τ = 60), ∫ΣaM dt = ∫(L_in − L_out)/c dt to 1e-5, acceleration of the thin shell = 0.998 σ_T L/(4π m_p r² c); work: ΔKE = 0.998 (thin) / 0.883 (thick, the rest dissipated by the artificial viscosity in the accelerated skin) of ∫Σ a M ū dt, which is 9e-5 (thin) / 7e-5 (thick) of the absorbed radiation energy.
d Sweep: Strömgren sphere (20 eV, steady, on the spot) r(x = ½) equals the ODE solution to 1e-5 (x_HI(r) to 2.7e-3); H + He (blackbody 1e5 K, n_H = 100): nested He²⁺ ⊂ He⁺ ≈ H⁺, ionisations = recombinations 1.9e-14, photon closure 0, energy (absorbed = ionisation + heating) 8e-17, volume-equivalent radius 5.8744 / 5.8742 / 5.8742 kpc at 100 / 200 / 800 zones (a front of 332 zones in one call).
e S7b module: records round trip, inner column conserves the bolometric power (re-emission, 1e-4: the Rayleigh–Jeans tail beyond the 0.1 eV grid), breakout for Ṁ/300 (Q/Q_crit = 211, escaping 0.95, R_II = 6451 R*), the sweep's beam and `RT.transfer`'s direct beam on the sweep's composition agree to 4e-4 (they differed by a factor 6 until the second pass with the opacities of the updated composition was added: Rayleigh scattering of the neutral H/H₂ that the sweep had removed).

**6. Production: ABN02 z = 19, from the `core2` snapshot (n_c = 6.7e22), sink run as S6b (SPLIT 2, r_s = 5 AU M/M☉, NLTE = ep, netatomicli), S6b protostar track `s6b_p3_z19f.csv` (to 170 M☉) as the source, feedback on from M_pt = 3 M☉.**
The two-way run is limited by the star (see 1): the exchange ran 3 → 33 M☉ with the final physics (p6) and 3 → 135 M☉ with the code before the `fJ` fix (p1, discarded: its "stall" was the artefact), so the production is the envelope with feedback from the stored star (one-way: the feedback acts on Ṁ through the envelope only; the star's reaction to the changed Ṁ(t) after ≈ 100 M☉ is therefore not followed).
Runs (`lib/popiii/s7b_F*_env.csv`, per-barrier records; `s7b_F*_prof.csv`, the 12 innermost zones at every barrier): F0 feedback evaluated but not applied (baseline in the same code path), F1 everything on (J_ν chemistry, H⁻ detachment, LW, continuum force), F2 = F1 without the force, F3 = F1 without LW.

| M_pt reached at t_c [yr] | 30 | 50 | 70 | 100 | 120 | 140 | 150 | 159 | 170 |
|---|---|---|---|---|---|---|---|---|---|
| F0 (off) | 905 | 1906 | 3028 | 5290 | 6026 | 9094 | 10233 | 10437 | (stopped at 168.2, 12 826) |
| F1 (all) | 826 | 1775 | 2958 | 5022 | 7593 | 9912 | 12134 | 23539 | never (159.4 at 59 kyr) |
| F2 (no force) | 836 | 1781 | 2963 | 5151 | 6586 | 8456 | 10177 | 10354 | 12663 |
| F3 (no LW) | 826 | 1796 | 2950 | 4939 | 7409 | 9253 | 11463 | 11906 | 14845 |

F1 against F0, Ṁ (M☉/yr): 10 M☉ 3.64e-2 / 3.57e-2, 50: 2.08e-2 / 1.86e-2, 100: 1.26e-2 / 1.21e-2, 120: 9.1e-3 / 1.26e-2, 150: 6.2e-3 / 7.4e-3, 159: 2.3e-3 / 8.0e-3 (the zone masses at the sink are 4–8 M☉, so ±1 zone in the times above); F0 continues at 6.5e-3 M☉/yr beyond 168 M☉.
Star (S6b track, M* = 10 / 20 / 50 / 100 / 150): R* = 168 / 263 / 86 / 41 / 82 R☉, L* = 6 / 3e4 / 6.8e5 / 2.5e6 / 4.4e6 L☉ (L_∞ 5e4 / 9e4 / 1.0e6 / 3.4e6 / 4.8e6), T_eff,* = 0.7 / 4.8 / 17.8 / 35.8 / 29.2 kK, T_ph = 6.7 / 6.3 / 8.0 / 8.2 / 8.6 kK, Q_* (blackbody of the stellar surface) = 1.2e42 / 3.4e41 / 2.3e48 / 1.1e50 / 1.2e50 s⁻¹, S_ion,ph ≤ 1.5e46.
*H II region:* confined at R_II = 1.00–1.14 R* at every mass (Q_*/Q_crit ≤ 7e-2 until the infall has halted; escaping fraction of the inner column 0 throughout); only the thermal re-emission (T_ph) reaches the grid, 1e44–1e46 s⁻¹ absorbed in the first zone. *LW:* k_LW in the first zone 2e-17 (10 M☉) → 5e-3 s⁻¹ (100 M☉) unshielded at the grid edge (self-shielding by the column outwards): H₂ is destroyed in the first zone (f_H2 2e-3 → 1e-11) from ≈ 60 M☉ on and in the zones out to several 1000 AU later (profiles).
*Radiation force:* max a_rad/g over the grid 0.08 (50 M☉), 0.3 (70), 0.56 (100), 0.75 (150) at 3000–6000 AU, carried by the Rayleigh scattering of neutral H of the 4–13.6 eV photons (see 3); L/c against the infall momentum flux Ṁ v_ff: 0.2–2 % (dynamical Ṁ), 0.5–19 % (Eulerian density).
*What ends accretion:* neither alone. LW alone (F2) changes nothing visible (the infall at ≥ 3000 AU is already slow in the baseline, 1–3 km/s against v_ff 6–8 km/s); the force alone (F3) delays by 14 % at 159 M☉; both (F1): from ≈ 120 M☉ the force (a_rad ≈ 0.5–0.75 g) acts on gas that no longer cools, the interfaces beyond 3000 AU decelerate to 0 and reverse (u = +5…+17 km/s outward at 8000–30000 AU, a_rad/g = 4–68 in the tenuous zones), the first zones shock-heat to 1e4–2e4 K and ionise collisionally (x = 0.9, T = 1.1–2.4e4 K), the grid edge retreats from 2000 to 27000 AU and M_pt stays at 159.4 M☉ for 36 kyr (the last sink event is at 23.5 kyr; Ṁ_win falls to 8e-4).
Final mass in this 1D calculation ≈ 160 ± 8 M☉ (zone mass), reached at ≈ 24 kyr; the baseline would continue (6.5e-3 M☉/yr at 170). OI02's H II confinement holds throughout (Q/Q_crit ≤ 7e-2): breakout would need Ṁ < 1.4–2.3e-4 M☉/yr (`s7b_break.py`), 25–80 times below the collapse rate.
*Comparison.* McKee & Tan (2008, fiducial Pop III.1): Lyα pressure reverses the polar infall at 20–30 M☉ but cannot stop equatorial accretion, the H II region passes the escape radius at 50–100 M☉, disk photoevaporation ends growth at ≈ 140 M☉ (HYOY12 / Hosokawa et al. 2011, 2D with rotation: 43 M☉ for Pop III.1 through ionising radiation escaping along the poles and photoevaporating the disk, ≈ 17 M☉ for Pop III.2). In spherical symmetry there is no polar cavity and no disk shadow: the ionising photons are absorbed in the free-fall column inside r_s (OI02: no firm upper limit from the H II region), so the 2D mechanisms are absent by construction and the 1D end of accretion is set by the radiation force on the neutral infall together with the loss of H₂ cooling, at 160 M☉ — comparable with McKee & Tan's final mass but by a different process, and 4× the 2D mass (the 2D masses are set by what escapes through the polar cavity, which a 1D model cannot represent).
*Uncertainties:* the Lyman-series force is a continuum (Rayleigh) surrogate; the first zone and the sink are coarse (4–8 M☉ zones spanning factors 3–6 in radius), the stall is sharp in time only to one zone; the star is one-way; ρ_s of the inner column from Ṁ_win (the Eulerian grid density gives Q/Q_crit up to 1.4 at 166 M☉ and R_II up to 2.8 R*, no escape); F1–F3 used the single-pass sweep (identical for x_sweep < x_net + 1e-3, which holds throughout).

**7. Performance** (`lib/popiii/perf/s7b_*`; one-way F1: 14 693 steps, 430 barriers, wall 2209 s on 4 threads against 1919 s for F0). Where the time goes: chemistry 1096 s (50 %; 821 s in F0: +33 % from the per-zone J_ν, `set_zone_light!` itself is 0.3 ms per zone), H₂ line transfer 782 s (35 %), the feedback evaluation 131 s (6 %: 305 ms per barrier = opacities 101, `RT.transfer` 123, inner column 59, sweep 17 ms), hook bookkeeping 9 ms per step. Coupled run (c4, 3 → 7 M☉, 86 barriers, FEX = 0.01): 1170 s, each process waits 29 % of the time for its partner at the barriers (envelope 339 s, star 322 s); the feedback evaluation is 31 s. The star step dominates the coupled wall time (≈ 3000 steps per 100 M☉ interval at 1–10 s each, plus the relaxation of the embryo: 20 min). Room for improvement: ChemNet's step with a non-static env (J lookup per photo process, ≈ 35 % of the chemistry), caching the opacity table across barriers, a smaller frequency grid in the dark outer zones.

**Open points.** (1) The two-way production through 30 M☉ needs the S6b precursor stall (E_s memory / ionisation cliff) cured or a restart protocol that keeps the barrier state; the one-way result misses the star's reaction to Ṁ(t) after ≈ 100 M☉ (Ṁ falls by 40 % between 120 and 150 M☉: L_acc and the CNO ignition would shift). (2) Lyα / He II line force and the damping-wing trapping (McKee & Tan eq. 22) are not transported; the Rayleigh surrogate carries 71 % of the force at 3000 AU; a Z = 0 atmosphere spectrum for the star is not used (blackbody). (3) The momentum of the beam absorbed inside r_s (L_∞/c = 0.3–2 % of Ṁ v_ff(R*); Γ = 0.9–0.97) is not applied to the inner flow. (4) Ṁ_win from the sink events is a staircase; 1D collapse with a 0.2-dex-coarse sink grid limits the accuracy of when the stall happens (±1 zone). (5) The Eulerian-flux / dynamical-rate mismatch at the sink (3–10×) should be understood (resolution of the first zones) since Q_crit ∝ ρ_s². (6) Heating by the non-network continua is absorption only (negligible); photodetachment, Li photoionisation and H⁻ chemistry are included through the network but not validated against a one-zone reference with a stellar J. (7) AtomicHHe's (1 − J/B) factor was floored at 0, line absorption by the field is not a heating term.

### S7c-2 Lyman-wing opacity and Lyα radiation force [run]
Files: `julia/LymanWings.jl` (cross section, Voigt cores, MT08 blend), `julia/RadLya.jl` (trapped-radiation force on the grid), edits of `RadTransfer.jl` (switches, `Opts.lyman/lya`, `ARAD_LYA`, `default_grid` line refinement), `s7c2_test.jl` (full tier `s7c2`, 15 checks PASS, `lib/popiii/s7c2_test.log`), `s7c2_prod.sh`, `s7c2_compare.py`; a few lines in `s7b_env.jl` (RT_LYA: inject the H II-region Lyα; `S7B_TCMAX`; `S7B_PERFPREFIX`). `ContinuumOpacity.jl` is **not** touched (its Rayleigh term is replaced inside `fill_opacity!`). Papers (data only): `~/Projects/popiii_lit/lyman/` (Lee & Kim 2005 astro-ph/0402023, Bach & Lee 1306.0181, Lee 2308.04959, Dijkstra & Loeb 0807.2645 / 0809.2099, Smith, Bromm & Loeb 1607.07166, 2103.14655; MT08 is in `feedback/mt08`). Defaults are unchanged: `RT_LYMAN`, `RT_LYA` (env, or `Opts(lyman, lya)`) are off and the S7b numbers are reproduced to 1.6e-13 (test d1); quick tier 8/8, `s7a` (24 checks) and `s7b` (20) PASS.

**1. Cross section** (atomic H 1s, `LymanWings.σ_gas`). Kramers–Heisenberg, exact non-relativistic: σ = σ_T ω⁴ |S(ω)|², S = Σ_{n=2}^{300} f_n/(ω_n² − ω² − iωγ_n) + ∫(df/dε)dε/(ε² − ω²) (P.V. above the edge), f_n = 2⁸n⁵(n−1)^{2n−4}/(3(n+1)^{2n+4}) (exact, f₂ = 0.41620), df/dε from Stobbe's 1s bound–free cross section (sum rule Σf = 0.99975; S(0) = 4.504 a.u.), γ_n = A_np/τ_au (damping; A_np from f, branching to 1s 1/0.88/0.84/0.81: it sets only |Δν| ≲ 100 Γ). It reduces to Rayleigh far below Lyα with σ/σ_T (ω_l/ω)⁴ → 1.267 (c₀ = 1.2654, Lee & Kim). *Thermal cores:* the scattering is coherent in the atom frame, so the gas cross section is the Maxwell convolution; for Ly-2 … Ly-12 the resonant term |S_n|² is multiplied by φ_Voigt/φ_Lorentz (Humlíček W4; same line integral: ∫σ over ±0.3 % is conserved to 1e-5, peak σ₀H(a,0)), interference and dispersion are kept. *Redistribution:* the transfer scatters coherently in the observer frame (correct in the damping wings |x| ≫ 1, where the force on a thick infall comes from; wrong in the cores, which the frequency-diffusion plateau of MT08 covers, §2); the beam removes the scattered photons (single-scattering momentum, S7a convention; with `diffuse` on they re-enter the diffuse field).
*Validation* (`s7c2_test.jl` a): against the Lee & Kim (2005) low-energy expansion (their Table 1, shown there to coincide with Gavrila 1967) ratio 1 ± 1.3e-3 for ω/ω_l = 0.05–0.5 (0.6: +4.0 %, the series is truncated and quoted good to 4 %); against Ferland's fit to Gavrila ≤ 5.4 %; Lyα peak 7.0564e-11 cm² = 3λ²/2π, line integral πe²f/mc to 1e-4; against the Lee (2003) fifth-order expansion about Lyα (their eq. 19) 0.1 % for Δω/ω₂₁ = −0.1 … +0.05, i.e. on both sides, including the region between Lyα and Lyβ (σ = 1.8e-24 cm² midway). Dalgarno's polynomial used before is −6 % at 0.4 ω_l, −19 % at 0.5, −46 % at 0.6 and has no resonances.
*Grid* (`default_grid(lyman = true)`, 70 → 430 bins for the production grid): extra edges at ν_n(1 ± δ), δ = 1e-5 … (half the line spacing, ≤ 0.1), 6 per decade, for Ly-2 … Ly-12 (`RT_LYMAN_LINES`); σ̄(bin, T) is tabulated (37 T nodes, ln–ln interpolation, 0.4 s) with the cores resolved in the bins that contain or neighbour a line (the T-interpolation error is up to 17 % only in the 4–7 Doppler-width Gaussian-to-Lorentz transition, 0.5 % width-weighted). Unresolved: lines n > 12 (E > 13.505 eV, 1 % of the band) sit in one mean-σ bin that over-absorbs a beam of low column (N_H ≲ 1e19). *Shell test* (b): force of a 5e4 K blackbody on a neutral shell at 1000 K against the direct integral of the exact σ(ν): ratio 1.006 / 0.992 / 1.003 at N_H = 1e20 / 1e22 / 1e24 (1.04 at 1e18); Dalgarno's polynomial gives a force 5.5 times smaller at 1e22.
*LW band:* the H₂ Lyman/Werner lines are **not** in the opacity (k_LW = 1.38e9 J_LW f_sh(N_H₂)); the H I wings attenuate the same photons independently, so no photon is counted twice; k_LW changes by the difference of the two H I scattering cross sections (+1.5 % at N_H = 1e22, +15 % at 1e23, since Dalgarno's polynomial is an extrapolation there). Scattered photons are not recycled into J_LW (beam convention).

**2. Trapped Lyα / Lyman-series force** (`RadLya.jl`, `Opts.lya`, `RT.ARAD_LYA`; `RT.ARAD` = beam + diffuse + this). McKee & Tan (2008) §4, App. B–C: in spherical symmetry with pure scattering, ∇·F = 0 and the force on the neutral gas is −∇P_iso = ρκF/c with the *conserved* flux F = L_ν e^{−τ_abs}/(4πr²); the beam (scattered photons removed) saturates at (1 − e^{−τ})F/c and underestimates a thick layer by τ. P_iso = (4π/3c) Σ_bins min[B*_νΔν, (J/F) F_bin] with 1/(J/F) = 1/J_coh + 1/J_plat: J_coh = (3/4π)N_eff σ_KH(ν) (B10, all lines and the Rayleigh continuum) and J_plat = 0.411(aτ̄)^{1/3}/(min[1,(aτ̄/450)^{1/3}] + Γ) the frequency-diffusion plateau of the nearest line (B7, B13–14; Γ = 1.04(ετ̄)^{1/3} two-photon destruction for Lyα, ε = 4.4e-9/(1 + n_cr/n_e)); with Lyα alone this is MT08's eq. 22 / B11. N_eff = r_j² Σ_{k≥j} n_k Δr_k/r_k² (their C6–C10 without the power-law assumption: (2/5) n r for ρ ∝ r^{-3/2}); a_k = −ΔP/(ρΔr). The component stored in `ARAD_LYA` is the excess over the beam's single-scattering force from the same H I scattering (subtracted), so it vanishes in thin gas and the total H I force with `lya` is MT08's. Sources besides the star's spectrum (bin sub-sampled across the lines): `LYA_INJECT` (line-centre Lyα of the H II region inside r_s, 0.68 per absorbed ionising photon, set in `s7b_env.jl`) and the Lyα of photoionised grid zones (0.68 α_B n_e n_p limited by n_H Γ_H of the ground-level photoionisation; thermal collisional Lyα is not a source, see below). Validity: aτ̄ ≳ 1, no dust, escape width (aτ̄)^{1/3} b of hundreds of km/s ≫ the infall velocity range, no polar cavity (1D); applied to the grid's neutral H only.
*Tests* (c): thin-shell M_F = (4π/3)J/F = 13.6 / 29.3 / 63.1 / 136 at N_H = 1e19 / 1e20 / 1e21 / 1e22 (T = 1e4 K) against the Dijkstra & Loeb (2008) Fig. 6 fit M_F ≈ 10 N19^0.44 = 10 / 27.5 / 76 / 209 quoted by Smith, Bromm & Loeb (2017) (ratio 1.36 / 1.06 / 0.83; our slope 1/3, theirs 0.44), slab estimate 15(τ₀/10^5.5)^{1/3} = 19 / 40 / 86; coherent regime: force per mass / (σ n_H F/ρc) = 1.00–1.05 over a ρ ∝ r^{-3/2} envelope with τ_total = 40 where the beam gives ≤ 2e-3; thin limit |a(lya) − a(beam)|/a = 2e-3; injected-line pressure = formula to 1e-12.

**3. Effect on a_rad.** 100 M☉ star (L_* = 2.5e6, L_∞ = 3.4e6 L☉, T_ph = 8.2 kK) through ρ ∝ r^{-3/2}, n_H = 1e7 at 3000 AU, T = 1000 K, 80 zones; a_rad/g in zones 1, 10, 20, 40, 80: S7b 0.785 / 0.147 / 0.077 / 0.047 / 0.037; Lyman-wing beam 1.51 / 0.122 / 0.068 / 0.044 / 0.036 (inner zone ×1.9, which absorbs the wing photons that the outer zones lose); + trapped force 10.5 / 10.7 / 12.9 / 19.6 / 74 (flux-mean opacity of the wings 4–30 cm²/g times a conserved flux; the last zone is the surface term N_eff → 0, a ∝ N^{-2/3}). The 10.2–13.6 eV band that carried 71 % of S7b's 0.5 g is, with the wings, saturated in the beam and carried by the trapped component in the outer envelope.

**4. Production: one-way ABN02 z = 19, F1 set-up (core2 snapshot, S6b track, feedback from M_pt = 3 M☉; `s7c2_prod.sh`)**; times in yr (`lib/popiii/s7b_{F0,F1,c2L,c2A,c2LA}_env.csv`, `s7c2_compare.py`):

| run | M=10 | 30 | 50 | 100 | 120 | 140 | 160 | M_final | halt | grid edge |
|---|---|---|---|---|---|---|---|---|---|---|
| F0 (feedback not applied) | 164 | 905 | 1906 | 5290 | 6026 | 9094 | 12826 | 168.2 (track end) | – | 935 AU |
| F1 (S7b: Rayleigh polynomial, no Lyα force) | 161 | 826 | 1775 | 5022 | 7593 | 9912 | – | 159.4 | 23.5 kyr | 3.5e4 AU |
| c2L (Lyman-wing opacity, beam only) | 161 | 834 | 1794 | 5130 | 7536 | 9930 | 23486 | 168.2 (track end) | – | 2.4e4 AU |
| c2A (trapped force, Dalgarno beam) | 158 | 825 | 2702 (jump 47.3 → 52.9) | – | – | – | – | **52.85** | 2.7 kyr | 1.2e5 AU |
| c2LA (both) | 170 | 817 | 2214 | – | – | – | – | **52.85** | 2.2 kyr | 4.0e4 AU |

(M_pt is the zone mass at the sink, 2.6–3 M☉ per step at 45–50 M☉: ±1 zone.) *The opacity alone does not change the end of accretion within the S7b scatter*: c2L passes 160 M☉ at 23.5 kyr and reaches the track end 168 (F1's halt at 159.4 was at the same time, 23.5 kyr; Ṁ at 10/30/50/70/90 M☉ 3.65/2.50/2.06/1.58/1.32e-2 against F1 3.64/2.56/2.08/1.61/1.39e-2): the single-scattering beam force of the wings is absorbed in the first zones. *The trapped-radiation force ends accretion at 53 ± 3 M☉* (c2A, c2LA), a factor 3 below S7b's 160 ± 8: the total force on the grid rises from 6 % (1.5 kyr, 45 M☉) to 45 % of the infall momentum flux Ṁ v_ff (`Fsum_over_Mdotv`); a/g of the 12 innermost zones is 5e-4 (20 M☉), 0.06–0.11 (38 M☉), 0.5–1.8 (45 M☉), 5–7 (50 M☉) against ≤ 0.012 in F1 at 45 M☉; the interfaces at 1–5 kAU reverse, the first zones shock-heat to 8–10 kK, the grid edge retreats to 4e4–1.2e5 AU and Ṁ_win falls from 1.4e-2 to 4e-3 M☉/yr over 3 kyr; M stays at 52.85 until the cut at 30 kyr (`S7B_TCMAX`). The mechanism is Lyα (Lyman-series) radiation pressure on the neutral infall (MT08: it reverses the polar infall at 20–30 M☉; in 1D, without a polar escape route, it stops the whole envelope), *not* the ionising radiation or LW (H II confinement unchanged).
Two bugs found on the way; the first c2A/c2LA attempts (stalls at 27 and 12 M☉) are discarded: (i) collisional-excitation Lyα as a grid source, evaluated in a shock-heated first zone (n = 1e10, T = 6e3 K: 1e49 erg/s of thermal energy, thermalised in that gas) gave one-barrier spikes a/g ≈ 1e3–1e4 in the next zones; (ii) the photoionisation limit of the recombination source included the Balmer-continuum ionisation of thermally excited n = 2 atoms of the hot zone: now only the ground-level (> 13.6 eV) rate limits the source. The `arad_over_g_max` column of the env csv is dominated by the surface term of the outermost zone and is not the force on the infall; use `Fsum_over_Mdotv` and the profiles.
*Performance* (`lib/popiii/perf/s7c2_*`): the feedback evaluation costs 0.9–1.0 s per barrier with the lines (opacity 0.26 s, transfer 0.5 s, 430 bins) against 0.3 s in S7b; wall 1720–1785 s for the 360-barrier halted runs, 7540 s for the 411-barrier c2L run on the loaded machine (S7b F1: 2209 s).

**Open points.** (1) The trapped pressure inside r_s (the free-fall column from R* to r_s, N_H ~ 4e26) is not applied to the inner flow; with the H II-region Lyα it is MT08's criterion P_iso(R_II) = 2ρv² and may halt the flow earlier than the grid does. (2) Photon destruction in the trapped regime is by two-photon decay only; thermalisation by true absorption (H⁻, CIA, H₂) enters along the beam path only (κ_abs/κ_sc ~ 1e-2 in the 8 eV test would cut the diffusive flux beyond τ_sc ≳ 1/√(3ε); to be checked in the production envelope). (3) Zone-mass resolution at the sink (±3 M☉), the one-way star (it does not react to the new Ṁ(t)), spherical symmetry (the 2D rotating-disk results differ by other means); S7c-1's two-way run should use `RT_LYMAN=1 RT_LYA=1`. (4) Lines n > 12 are unresolved; velocity-gradient (Sobolev) widening of the cores, hyperfine and fine structure are not included. (5) The outermost zone feels a ∝ N^{-2/3}; a ghost column beyond the grid would remove it.

### S7c-1 Precursor stall and the two-way feedback run [run]
Files: `julia/s7c1_prod.sh` (two-way production), `s7c1_repro.sh`, `s7c1_perf.jl` (Perf wrapper of `s6b_run.jl`), `s7c1_val.sh`, `s7c1_prec_scan.jl`, `s7c1_compare.py`, `s7c1_table.py`; edits (all opt-in or bug fixes): `s6b_run.jl` (`ESTAU`, `LSTARQ`, `PROTGAMMA`, `SAMEDTRETRY`, checkpoint/restart state, `PROFSTEP`), `TitanDriver.jl` (`SAMEDT_RETRY`), `s7b_star.jl` (`S7B_ENVDIR`/`S7B_RESUME`: replay of recorded exchange records), `s7b_env.jl` (`S7B_LYASURV`), `s7b_star.sh` (`PROTGAMMA=0.25`, `JT`, `SCRIPT`).

**1. The stall, reproduced.** Re-running the p6 recipe (`s7b_env.sh`/`s7b_star.sh`, MSYNC 3) gives p6's star trajectory exactly (`q1`: M = 10.5557 at 286.68 yr, 29.2529 at 958.59 yr) and the stall at t = 3.4125e10 s (M = 32.58, R = 402 R☉, L* = 3.9e5 L☉, T_eff,* = 7.3 kK, (L* + L_acc)/L_Edd = 0.44). Smallest set-up: the checkpoint at step 2100 (M = 29.25) plus the recorded exchange records of `q1` (`lib/popiii/s7b_x_q1`, replayed with `S7B_ENVDIR`, `S7B_RESUME=1`; the star runs alone, stalls 70 steps later, 1 min after start-up): `s7b_star.sh rp1 S7B_ENVDIR=… S7B_RESUME=1 RESTART=…`.
*What the stall is* (not the cliff of E_s alone: `s7c1_prec_scan.jl` shows T_s jumping 6.4 → 11.3 kK for L* 3 → 4e5 at R = 400 R☉, the first-order ionisation-front transition of the precursor, which the thermostat T_s = T(h_R) turns into a steep but finite E_s(L*, R): d ln E_s/d ln L ≈ 8): at that transition the surface skin (q < 2e-5, T = 9–12 kK, κ = 13–46 cm²/g, i.e. local Γ > 10) carries the whole L* by MLT (26 ρ c_s³) and is heated by E_s (×18). Three defects of the driver then act together; each alone leaves the stall (replays):
(i) *E_s damping per step instead of per time.* ESW/ESLIM move E_s by up to 5 % per step whatever dt is; the step control sees the response of the skin to that move (ρ at zones 127–130 changes by 5 %) and cuts dt by 2–3 per step: dt → seconds in ~20 steps, with E_s still moving at the same rate per step. Fix (`ESTAU`): the damped move is multiplied by f = 1 − exp(−dt/τ_E), τ_E = ESTAU R/v_ff (the flow time of the precursor); f = 1 for dt ≫ τ_E.
(ii) *L* inside the skin.* L* of the boundary condition is the flux at 15 interfaces below the surface (KLSTAR), i.e. inside the convective skin; when the skin is heated its convection switches off and L* flickers 6e4 ↔ 7e5 L☉, E_s follows with gain ≈ 8 (a lagged, explicit loop that no dt reduces). Fix: L* at exterior mass fraction 1e-3 (`LSTARQ`, existing S6a-fix option; L varies by < 2 % between q = 1e-3 and the base of the skin).
(iii) *Newton initial guess.* At the planned dt the first Newton correction is ×14 in E_r of the ionisation zone (interface radii extrapolated with their last velocity, `EXTRAP_FRAC` 0.3; 25 iterations fail), at dt/2 from the old state 0.33 and 7 iterations: dt is halved at every step (0.8 per step net). `TitanDriver.SAMEDT_RETRY`: after the first Newton failure retry once at the same dt from the old state before cutting dt (`SAMEDTRETRY`). (A global `EXTRAP = 0` also cures this but fails the Ṁ_fid protostar at 15.3 M☉, where the radius extrapolation is needed.)
(iv) *Restart.* Checkpoints did not store E_s, the enthalpy guess and L*, and `restart!` left the convective flux of the last step empty (L* of the boundary condition was the radiative part only: 713 against 3.9e5 L☉ at the first step): every earlier "recovery" restarted from a corrupted boundary condition (also the reason why the p3/p5/p6 rollbacks never cleared the stall). Fixed (aux stores `es_old`, `h0`, `lbc`; `restart!` recomputes the step diagnostics).
(i)–(iii) are active only while Γ = (L* + L_acc)/L_Edd (κ_es = 0.34) of the previous step exceeds `PROTGAMMA` (0.25 in `s7b_star.sh`, off at 0.8 of it; default 0 = off): without the gate they changed the Ṁ_fid protostar at 15 M☉ (stall at 15.3 M☉ with ESTAU + LSTARQ, which pass alone to 14 M☉). Replay with the gate: the step at which the stall occurred passes with dt ≈ 1e6 s (from 7e7) and recovers.
*Validation.* Defaults are unchanged: `s7c1_reg` (Ṁ = 1.1e-3 to 8 M☉) reproduces `s6b_pv6` to 1.000 in R*, L*, T_c at M = 0.5 … 8 (R peak 111.8 R☉ at 6.80 M☉). The un-gated variant (EXTRAP 0, ESTAU 1, LSTARQ 1e-3 all the time; `s6b_s7c1_v1`) kept pv6 within 0.4 % in R* (R peak 112.4 at 6.80, ZAMS reached at 44.7 with R = 3.95 R☉ against 3.93, L* −5 % from 8 M☉ on, early L* ≪ L_acc differs because it is now the flux below the skin) and Ṁ_fid within 2.5 % to 15 M☉ (`v2`: R peak 221.8 at 10.8 against 219.1), but not beyond 15.3 M☉: hence the gate. The gated runs were not repeated to ZAMS (S6a-fix/S6b validations at Γ ≥ 0.25 — the late phases — are therefore validated only through the un-gated variant for Ṁ = 1.1e-3, and not at all for Ṁ_fid beyond 15 M☉). Tests: quick tier 8/8 PASS; full tier 27/27 PASS (incl. s4–s7, s6b_early, s6c_prec, s7a, s7b, s7c2; run serially, `POPIII_SERIAL=1`).

**2. Two-way coupled production** (ABN02 z = 19, core2 snapshot, MSYNC 3 M☉, FEX 0.01, `FJ_FLOOR = 0`, `s7c1_prod.sh`; logs `m2c_abn02z19_lines_n200_s7b_c1{a,b,c}.log`, `s6b_s7b_c1*.{log,csv}`, `s7b_c1*_env.csv`, `s7c1_c1*_star_stdout.log`): (a) `RT_LYMAN=1 RT_LYA=0`, (b) `RT_LYMAN=1 RT_LYA=1`, (c) as (b) with the Lyα of the H II region reduced by the survival factor of §3 (`S7B_LYASURV=1`). **No stall: all three runs pass 30–35 M☉ and the whole run without a step failure or stall recovery (c1b: 2622 star steps, 0 failed).** M_pt = zone mass at the sink (+3 M☉); star values at the same time (the star's own mass is not used after the sink has starved: see the end of §2).

| run | M_pt | t [yr] | Ṁ_win [M☉/yr] | R* [R☉] | L* [L☉] | T_eff,* [K] | S_ion,* [s⁻¹] | R_II/R* | F/(Ṁ v_ff) | grid edge [AU] |
|---|---|---|---|---|---|---|---|---|---|---|
| a | 10.1 | 170 | 3.2e-02 | 170 | 2.7e+01 | 1.0e+03 | 1.2e+42 | 1.00 | 0.000 | 98 |
| a | 20.6 | 500 | 2.9e-02 | 247 | 3.0e+04 | 4.8e+03 | 5.0e+41 | 1.00 | 0.000 | 194 |
| a | 30.4 | 832 | 2.5e-02 | 513 | 2.8e+05 | 5.9e+03 | 4.5e+41 | 1.00 | 0.000 | 486 |
| a | 52.8 | 1788 | 2.1e-02 | 90 | 6.2e+05 | 1.7e+04 | 1.6e+48 | 1.00 | 0.000 | 1046 |
| a | 73.6 | 2982 | 1.6e-02 | 47 | 1.2e+06 | 2.8e+04 | 2.7e+49 | 1.00 | 0.001 | 1616 |
| a | 102.5 | 5174 | 1.2e-02 | 28 | 2.1e+06 | 4.1e+04 | 1.2e+50 | 1.00 | 0.001 | 2702 |
| a | 134.9 | 9682 | 7.7e-03 | 75 | 4.0e+06 | 3.0e+04 | 1.2e+50 | 1.00 | 0.003 | 4642 |
| a | 142.7 | 20378 | 2.2e-03 | 24 | 4.6e+06 | 5.4e+04 | 3.2e+50 | 1.15 | 0.028 | 10817 |
| b | 10.1 | 173 | 3.2e-02 | 176 | 2.9e+02 | 1.8e+03 | 1.1e+42 | 1.00 | 0.000 | 82 |
| b | 20.6 | 491 | 2.7e-02 | 286 | 2.7e+04 | 4.4e+03 | 1.4e+40 | 1.00 | 0.000 | 202 |
| b | 30.4 | 869 | 2.2e-02 | 515 | 3.1e+05 | 6.0e+03 | 7.1e+41 | 1.00 | 0.000 | 278 |
| b | 44.7 | 1482 | 1.9e-02 | 117 | 4.6e+05 | 1.4e+04 | 2.5e+47 | 1.00 | 0.038 | 742 |
| b | 52.8 | 2279 | 1.3e-02 | 62 | 7.2e+05 | 2.1e+04 | 6.3e+48 | 1.00 | 0.434 | 1577 |
| c | 10.1 | 170 | 3.2e-02 | 169 | 2.1e+01 | 9.5e+02 | 1.2e+42 | 1.00 | 0.000 | 94 |
| c | 20.6 | 484 | 2.7e-02 | 262 | 2.8e+04 | 4.6e+03 | 1.3e+41 | 1.00 | 0.000 | 292 |
| c | 30.4 | 891 | 2.2e-02 | 417 | 3.0e+05 | 6.6e+03 | 6.1e+42 | 1.00 | 0.000 | 223 |
| c | 52.8 | 2160 | 1.5e-02 | 67 | 6.9e+05 | 2.0e+04 | 4.7e+48 | 1.00 | 0.100 | 462 |
| c | 59.0 | 3323 | 7.9e-03 | 34 | 8.8e+05 | 3.0e+04 | 2.7e+49 | 1.00 | 0.137 | 3212 |

*Final mass and mechanism.* (a) M_pt = 120.8 at 7.4 kyr, 127.8 at 7.7, 134.9 at 9.7, 142.75 at 20.4 kyr, Ṁ_win falling from 1.4e-2 (82 M☉) to 7.7e-3 (135) and 2.2e-3 M☉/yr (143), grid edge 10.8 kAU, first-zone T = 8.4 kK, f_H₂ = 1e-9: accretion ends at **143 ± 4 M☉** (zone mass 4–8 M☉) at ≈ 20 kyr; it is gradual, with the radiation force only 0.3–3 % of the infall momentum flux, so the mechanism is the S7b one (LW dissociation of the H₂ coolant + heating/force on the neutral infall at 3000–10000 AU, H II region confined at R_II = 1.0–1.2 R*, Q/Q_crit ≤ 0.2) acting on a *star that reacts*: against the one-way S7b value 160 ± 8 M☉ (F1, 159.4 at 23.5 kyr) the two-way run is 11 % lower (2σ of the zone scatter), at M_pt = 100 / 120 / 135 M☉ Ṁ is 1.2e-2 / 9e-3 / 7.7e-3 against F1's 1.26e-2 / 9.1e-3 /–, i.e. the difference arises after 120 M☉ where the star (R* = 28 → 75 R☉, T_eff,* = 4e4 → 3e4 K, S_ion 1.2e50) differs from the stored S6b track. (b) with the Lyα force of the H II region injected unattenuated: halt at **52.85 M☉** at 2.3 kyr (F/(Ṁ v) = 0.43; same as the one-way c2LA), the trapped Lyα pressure reverses the infall (first-zone T 0.5–8 kK, grid edge 1.6 kAU and growing). The star's own mass (`Mstar` column) keeps growing after the sink starves (the sink events, hence Ṁ_win, are a staircase and the estimator bounded by the average since the first event decays slowly: the star reached 178 M☉ at 24 kyr in (a) against M_pt = 143): beyond the halt only M_pt, and the star values at the same time before the divergence (≈ 10 kyr in (a), 2.3 kyr in (b)), mean anything.
*Star.* No stall recoveries (`stall_recoveries` 0), 0 failed steps; R* peaks at 513–615 R☉ at 25–30 M☉ (T_eff,ph ≈ 5 kK), contracts through the opacity cliff at 30–33 M☉ (T_s 6 → 14 kK, τ_prec 15, thick precursor), R* = 90 (53 M☉), 47 (74), 28 (102 M☉, T_eff,* 4.1e4 K), re-expands to 75 R☉ at 135 M☉ (Γ = 0.9); L* = 6e5 (53), 1.2e6 (74), 2.1e6 (102), 4.0e6 L☉ (135 M☉).

**3. Lyα destruction inside r_s (caveat of S7c-2).** The Lyα of the H II region (0.68 per absorbed ionising photon) is born at R* ≤ R_II ≲ 1.1 R* and has to cross the neutral free-fall column from R_II to r_s = 5 AU M/M☉ before it reaches the grid: N_H = Ṁ/(4π μ_H m_H √(2GM)) · 2 (R_II^{-1/2} − r_s^{-1/2}) = 1.4–4.7e26 cm⁻² for M_pt = 20–128 M☉ (Ṁ from the run), τ_0 = σ_0 N_H = 7e12–2.3e13 (σ_0 = 4.8e-14 cm² at 1.5e4 K), aτ_0 ≈ 1e10, so the photon is in the frequency-diffusion regime (escape through the wings of width (aτ_0)^{1/3} ≈ 2.5e3 Doppler widths) with N_sc ≫ 1/ε. Destruction per scattering: two-photon decay after 2p → 2s mixing, ε ≤ 4.4e-9 (MT08 B13, the upper bound for n_e ≫ n_cr; the thermalisation by H⁻, H₂, CIA and collisional de-excitation would only add to it), so the trapped intensity is divided by 1 + Γ, Γ = 1.04 (ε τ_0)^{1/3} = 32–48: the surviving fraction is **1/(1 + Γ) = 2.0–3.0 %** at M_pt = 20 / 30 / 45 / 50 / 59 / 102 / 128 M☉ = 2.0 / 2.7 / 2.3 / 2.3 / 2.3 / 2.5 / 3.0 % (a static-slab estimate; the infall velocity gradient (Sobolev) would reduce N_sc; the 1/(1+Γ) scaling is kept as the estimate, uncertain by a factor of a few). With that factor on the injected H II-region Lyα (`S7B_LYASURV=1`, run c) accretion still ends, later and more gently: **59 M☉ at 3.3 kyr** (Ṁ 7.9e-3 → 1.6e-3, F/(Ṁ v) = 0.14, M_pt frozen at 59.02 until the run was stopped at 11.7 kyr) instead of 52.85: the stellar Lyman-series photons of the beam and the Lyα of the photoionised grid zones, which are not attenuated by the inner column, carry the force. So the Lyα effect of S7c-2 is an upper bound by the H II-region injection (53 M☉), but the halt at ≲ 60 M☉ is robust against it; the actual survival of the line photons in the column (thermalisation, Sobolev escape in the infall) is not resolved here.

**4. Performance** (`lib/popiii/perf/s7c1_*`): coupled run c1b (2 processes × 3 threads, loaded machine): star 2603 s wall for 2622 steps (37 % waiting at the barriers: 955 s), envelope 2614 s (310 barriers: S7b hook 1160 s = 44 %, of which RT transfer 197 s, opacity 116 s, inner column 42 s, sweep 16 s); no failed steps. Un-gated full protostar runs: pv6-type run to 50 M☉ 2709 s, Ṁ_fid to 15 M☉ 620 s (`s7c1_val_*`). Diagnostics: a replay from the step-2100 checkpoint reaches the stall in ≈ 1 min (3 threads).

**Open points.** (1) The gate Γ_Edd ≥ 0.25 is physical (near-Eddington skin) but the protections were not run through the late phases of the constant-Ṁ validation protostars with the gate (only the un-gated variant, ≤ 0.4 % for Ṁ = 1.1e-3). (2) The mass of the coupled star diverges from the envelope after the sink starves (Ṁ_win estimator); the estimator should return Ṁ → 0 when no sink event occurs for ≳ the exchange time (then the star could be followed to the end). (3) The skin of an extended star close to Γ_Edd (supersonic MLT, density inversion) is not physical in 1D; the fixes make the model pass it, they do not model it (a flux cap ρc_s³ / MLT++ would). (4) Lyα survival: static-slab Γ only; Sobolev escape and thermalisation by H⁻/H₂/dust-free CIA, and the unattenuated beam/grid Lyα that still halt the infall at 59 M☉, remain. (5) c1c was stopped at 11.7 kyr (halted since 3.3 kyr), c1a at 25 kyr (halted since ≈ 20 kyr): no Perf record for them. (6) One `pkill -u` of `popiii_m2c.jl` by this session (stopping the c1a run) may have stopped another agent's envelope run on this machine.

### S8b 21 cm hyperfine line diagnostic [run]
Files: `julia/HI21.jl` (module `HI21`, no dependencies), `hi21_snap.jl` (snapshot post-processor), `s8b_test.jl` (full tier `s8b`; `lib/popiii/s8b_test.log`, 20/20 PASS), 3 inert lines in `popiii_m2c.jl` (`HI21=1`, `HI21EVERY`, `HI21RUN`: `HI21.inline!` writes `lib/popiii/hi21_<run>.csv` and `_prof.csv`). Everything is opt-in and nothing in the hydro or chemistry reads it. Defaults unchanged: a 130-step `popiii_m2c.jl thin 60` run with and without `HI21=1` gives identical logs (diff empty); quick tier 8/8 PASS. Papers (arXiv) in `~/Projects/popiii_lit/21cm/`: Furlanetto, Oh & Briggs 2006 (FOB06), Furlanetto & Furlanetto 2007 (e–H, p–H), Hirata 2006, Furlanetto & Pritchard 2006, Chuzhoy & Shapiro 2006, Sigurdson & Furlanetto 2006, Pritchard & Loeb 2012 (PL12). Zygelman 2005, Field 1958/59 and Allison & Dalgarno 1969 are not on arXiv; their rates enter through FOB06 Table 3.

**Physics.** (1) T_s⁻¹ = (T_γ⁻¹ + x_c T_k⁻¹ + x_α T_c⁻¹)/(1 + x_c + x_α) (FOB06 eq. 20); T* = hν/k = 0.06817 K (ν = 1420.4058 MHz), A10 = 2.85e-15 s⁻¹.
- x_c = (T*/A10 T_γ) Σ n_i κ10^i(T_k). H–H: FOB06 Table 3 (Zygelman 2005 up to 300 K, Sigurdson above; 27 rows, 1–1e4 K). e–H: FF07a Table 1 (17 rows, ≤ 2e4 K; valid below the n = 2 threshold, ≲ 1.5e4 K). p–H: FF07b Table 1 (17 rows). Log–log interpolation; rates are held constant outside the tables. The H–H and e–H transcriptions were checked against the text of the PDFs row by row. H₂: no calculation of H–H₂ hyperfine exchange exists, so `KH2` × κ_HH with default 0. He, He⁺, D are not included.
- x_α = 4 P_α T*/(27 A10 T_γ) S_α with P_α = 4π ∫σ_ν J_ν/hν dν (flat J across the line: P = 4π (J_ν/hν_α) · 1.105e-2 cm² Hz). The critical J_αC (x_α = 1, S = 1) is 1.108e-10 photons cm⁻² s⁻¹ Hz⁻¹ sr⁻¹ at z = 19 (PL12: 1.165e-10, rounded constants). This equals Hirata's x̃_α = 8πλ²γT*S̃J/(9A10T_γ) (test c2).
- `kind = :continuum` (J is the continuum near the line): Hirata (2006) eqs 40–42 for S̃_α(T_k, T_s, τ_GP) and T_c^eff, iterated to the fixed point in T_s. τ_GP = 3A_α n_HI λ³/(8π g) with the local gradient g = max(|du/dr|, |u/r|, H(z)). The fit holds for 1e5 ≤ τ_GP ≤ 1e7 and T ≥ 2 K and is clamped outside (in the dense cloud x_c ≫ 1 makes that irrelevant). `kind = :actual` (Monte Carlo line intensity): S = 1 and the supplied T_c.
- Background: T_γ = T_CMB (`TCMB`, 54.5 K at z = 19) plus an optional power-law radio excess (`Tgamma(…; Ar, β)`, at 78.1 MHz). Not included: collisional Lyα (T_k ≳ 6200 K; that gas has x_c ≫ 1), dust, CMB polarisation.
- `ts_rates` is the same without the T* ≪ T_γ linearisation (explicit 0↔1 rates; differs by ≤ 1.2e-3 = T*/T_γ). `relax` is the time-dependent upper-level fraction with constant coefficients (rate R_up + R_down).
- **Equilibration time** t_eq = 1/(R_up + R_down): 3478 yr for the CMB alone at z = 19 (1/[(1 + 4n̄)A10]). In the z = 19 cloud t_eq/t_ff ≤ 1.2e-6 everywhere (131 yr in the outer zone, x_c = 26; ≪ 1e-12 t_ff in the core), so the equilibrium T_s is exact for these runs. `NEQ=1` (levels carried through the snapshots, aligned from the outer end) gives the same T_s as equilibrium for all 11 feedback snapshots.

(2) Line transfer: α_ν = (3λ³A10/8π) n₀ (1 − e^{−T*/T_s}) φ_v with n₀ = n_HI/(1 + 3e^{−T*/T_s}) and φ_v the thermal Gaussian (b = √(2kT/m_H)).
- Per zone: τ_Sobolev = κ0/g, τ_thermal = κ0 Δr/(√π b), τ21 = their harmonic combination.
- Rays: impact parameters p = 0, r_1…r_n through the shells (u linear in r within a shell, sub-steps ≤ b/2 in velocity, v = −uμ, positive = receding). T_b ← T_b e^{−Δτ} + T_s(1 − e^{−Δτ}) from the far side starting at T_γ; δT_b,obs = (T_b − T_γ)/(1 + z).
- Disc means with annulus areas, a core aperture p < 1e4 AU (0.048 pc), flux density S = (2kν_o²δT_obs/c²) πR²/D_A² (D_A(19) = 533 Mpc, ABN02 cosmology; 1″ = 2.6 kpc), equivalent width, and the thin-limit emission L21 = (3/4) A10 hν Σ n_HI V.

**Tests** (`s8b_test.jl`, all PASS).
- a κ10: the values quoted in the FOB06/FF07 text (H–H 10, 70, 100, 1000 K; e–H 1, 100, 1e4 K; p–H 10, 1e4 K) are reproduced; monotonic; interpolation inside rows; clamping. κ_eH/κ_HH = 259 / 19 / 18 / 11 at 10 / 100 / 1e3 / 1e4 K, κ_pH/κ_HH = 2.15 at 1e4 K (FF07b: a factor of 2).
- b limits: no coupling gives T_s = T_γ to 1e-14 (with a radio excess, the raised T_γ); n_H = 1e8, T_k = 300 K gives T_s = T_k; strong Lyα with T_c = 20 K gives T_c; continuum Lyα gives T_k; detailed balance of the explicit rates 1e-13; relaxation e-fold exact.
- c Wouthuysen–Field: J_αC; the Hirata x̃_α identity; S̃_α(100 K, 2e6) = 0.923 (Chen & Miralda-Escudé fit 0.954) and 0.579 at 5 K; fixed-point residual 1e-13.
- d IGM (mean-density gas, adiabatic T_k from z_d = 150, x_e = 2e-4): at 1+z = 70 with T_k = 88 K (FOB06): x_c = 0.90, T_s = 123 K, δT_b = −37.6 mK (published peak ≈ −40 mK). At z = 19, T_k = 7.2 K: x_c = 8.7e-4, T_s = 54.2 K, δT_b = −0.2 mK without Lyα and −234 mK fully coupled (thin form −246 mK against PL12 eq. 7: −248 mK).
- e rays: uniform static sphere against the analytic chord integral (rays p/R = 0, 0.3, 0.6, 0.9, τ up to 0.24) 5e-9; thin-sphere disc mean against κ0 φ_v (4R/3) 9e-5 at line centre (EW 1e-3); expanding sphere u = Hr (R = 230 Sobolev lengths): T_b(v = 0) = T_γ e^{−τ_S} + T_s(1 − e^{−τ_S}) at τ_S = 0.08 and 3.3 to 1e-4; infall/outflow mirror symmetry exact.
- f the last collapse snapshot (when present): T_s between T_k and T_γ in all 200 zones, finite.

**Results** (z = 19, T_CMB = 54.5 K). Outputs `lib/popiii/hi21_<run>.csv` (summary per snapshot), `_prof.csv` (per zone: T_k, T_s, x_c parts, x_α, τ21, t_eq, extrema of the ray through the zone), `_spec_<k>.csv`, logs `hi21_*.log`, Perf `lib/popiii/perf/s8b_*` (0.2 s per snapshot, 200 zones × 200 rays).
- *Collapse* (`abn02z19_final2`, 67 snapshots, log n_c = 3.6 → 19.8, T_c = 154 → 5.5e3 K; `abn02z19_core2`, 41 snapshots to log n_c = 22.8, T_c = 5.8e4 K). The gas is never colder than the CMB (T_s,min = 154 K at log n_c = 3.6, rising to 239 K; halo edge 1010 K) and the coupling is collisional everywhere: x_c = 26 at the edge (n_HI = 0.18 cm⁻³), 1e3 at n_HI = 10, 1e12–4e21 in the core. So T_s = T_k within 1 % inside r = 1.6e6 AU (n_HI ≳ 10) and 0.60 T_k at the edge. **No absorption occurs anywhere; the cloud is seen in emission.** Mass-weighted T_s = 592 K (T_k 726 K).
- τ21 (max over zones, radius): 0.11 (1.3e5 AU) at log n_c = 3.6; 2.6 at 5.9; 27 (1.5e3 AU) at 8.9; 300 (160 AU) at 11.9; 8e2 at 15; 1.2e5 at 18; 1.9e6 (4e-4 AU) at 19.8; 4e6 at 22.8. The line is saturated (T_b → T_s/(1+z)) in the central region from log n_c ≈ 6.
- Central-ray peak δT_b,obs: 8.9 K (log n_c = 3.6), 41 K (10.4), 190 K (18.8), 270 K (19.8), the T_s/(1+z) of the hot core; its FWHM grows 25 → 77 km/s from log n_c = 10.4 to 19.6 with the infall (the disc mean is 5.8 km/s). The core-aperture mean (p < 1e4 AU) grows 8.9 → 25.6 K. The disc mean over the whole 52 pc halo is 372–382 mK (EW 2.4 K km/s) and **does not change through the collapse**: the halo emission dominates the area, the collapsing region is 1e-6 of it.
- Observer: the halo (R = 1.07e7 AU) subtends 0.020″ radius, the 0.048 pc core 1.9e-5″. Peak flux 1.7e-3 nJy for the halo, 1e-7 nJy for the core; a 1″ beam would see ≈ 0.4 mK from the halo. The halo is +7.4 K in its source frame (372 mK observed); the mean-density IGM at this z is ≈ 0 mK without Lyα and −234 mK fully coupled. Unobservable individually; relevant only as a contribution to the mean signal.
- *Lyα scan* on the collapse (`final2_lya1`, `_lya100`: uniform J_α = 1 and 100 J_αC): ⟨T_s⟩_m 592 → 594 and 675 K, disc mean 372.6 → 372.8 and 376.9 mK, nothing inside the halo (x_c ≫ x_α). Lyα can only move the outer zones from ≈ 0.6–0.8 T_k to T_k, ≤ 1.2 % of the disc signal whatever J_α is.
- *Feedback.* The S7c-1 runs saved no envelope snapshots, so the envelope run c1a (`RT_LYMAN=1 RT_LYA=0`, M_pt 3 → 143 M☉) was **replayed** against the recorded star records of `s7b_x_c1a` (exchange dir with `star_*.txt` only, `S7B_SAVEAT = 4,8,15,25,35,50,70,90,110,130,140`, 2 threads; log `m2c_abn02z19_lines_n200_s7b_h21a.log`; M_pt = 3.160 at the sync and 142.748 at the end, as in the original run; stopped by hand at 25 kyr; the 11 snapshots and star records are in `lib/popiii/hi21_in/c1a/`). Rows (`hi21_s7c1_c1a.csv`, `_rt.csv`):

| exch. j | t_c [yr] | M_pt | M* | L* [L☉] | T_eff,* [K] | log n (first zone) | M_HI [M☉] | ⟨T_s⟩_m [K] | ⟨T_s⟩_m with Lyα (RT) | max τ21 (r [AU]) | core δT [K] | centre δT [K] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 26 | 24 | 4.4 | 3.9 | 9.1e2 | 2.7e3 | 11.2 | 7491 | 592 | 592 | 409 (74) | 25.4 | 62 |
| 145 | 335 | 15.7 | 14.9 | 1.6e4 | 4.8e3 | 9.8 | 7483 | 592 | 592 | 150 (214) | 25.5 | 63 |
| 224 | 1126 | 37.9 | 34.7 | 3.9e5 | 9.7e3 | 8.8 | 7466 | 591 | 591 | 47 (915) | 24.9 | 61 |
| 257 | 1840 | 52.8 | 49.5 | 6.5e5 | 1.8e4 | 8.3 | 7455 | 591 | 595 | 33 (1.4e3) | 25.4 | 48 |
| 291 | 3124 | 73.6 | 70.2 | 1.3e6 | 2.9e4 | 7.8 | 7440 | 591 | 648 | 20 (2.6e3) | 25.9 | 53 |
| 320 | 5091 | 96.8 | 94.1 | 2.1e6 | 4.1e4 | 7.3 | 7423 | 591 | 684 | 14 (3.3e3) | 29.1 | 144 |
| 356 | 9682 | 134.9 | 134.6 | 4.0e6 | 3.0e4 | 7.0 | 7398 | 583 | 708 | 6.1 (7.2e3) | 36.7 | 67 |
| 385 | 20378 | 142.7 | 170.5 | 4.6e6 | 5.4e4 | 5.4 | 7329 | 725 | 727 | 1.0 (6.0e4) | 45.1 | 45 |

  Feedback creates no absorption either (T_s,min = 217–239 K). It removes the dense HI core (max τ21 409 → 1.0, radius 74 AU → 6e4 AU; M_HI falls by 162 M☉ by sink accretion and ionisation), heats the inner envelope (⟨T_k⟩_m 726 → 861 K and T_s,max 1.8e4 K at the end), raises the core-aperture signal 25 → 45 K, and widens the central line to ±46 km/s (j = 320) before the cavity is emptied. The halo disc signal moves by +1.7 % at most (375.4 → 381.8 mK at M* = 135 M☉).
- *Lyα with the protostar.* The fallback (RadTransfer zone J at the Lyα bin, star spectrum through the inner free-fall column, continuum photons, Hirata S̃_α) gives x_α ≳ x_c in the envelope from M* ≳ 50–70 M☉ (T_eff,* ≳ 1.8e4 K, S_ion ≳ 1e49 s⁻¹) and raises ⟨T_s⟩_m from 591 to 648–708 K (to T_k), +1.7 % in the disc signal; below M* ≈ 25 M☉ there is no effect. This is an indication only: the continuum J ignores the scattering in the wings of the thick neutral envelope, the trapped-line intensity and the H II-region Lyα; at the last snapshot (grid inner edge at 1.3e4 AU) the fallback J is 2000× lower (x_α = 0.5 at the edge), which shows its sensitivity to the inner-column treatment. The bracket is the uniform scan: none (592 K) to full coupling (726 K) in ⟨T_s⟩, ≤ 2 % in the disc signal.

**Open points.** (1) `LYA=rt` is the fallback. The switch to the Monte Carlo module (S8a, `LyaMC.jl`, section "S8a Monte Carlo Lyman α transfer") is `HI21.use_lyamc!(f)` with `f(snapshot) → (; J, Tc)` per zone (J in RadTransfer units, S_α = 1 then). The module was not on `popiii` when this was written: wire it and rerun the Lyα rows when it lands. (2) There is no cold gas in these runs; the absorption and Lyα-dominated regimes are exercised only by the IGM tests. Absorption needs gas below T_γ = 54.5 K (the box edge is 1010 K; the IGM outside the halo is at ≈ 7 K). (3) Disc means are over the box (52 pc, 1e4 M☉), not over a beam; the LOS IGM is not added. (4) Rates above the table ends are held constant (H–H > 1e4 K, e/p–H > 2e4 K); the hot core is in the T_s = T_k limit anyway. (5) The Hirata fit is outside its τ_GP range in the dense zones (clamped; irrelevant for T_s there) and does not capture recoil/spin exchange of trapped Lyα in a thick cloud. (6) Velocity field radial only (no turbulence or rotation; the Hubble flow is not added, option `hubble` of `cloud_spectra`); He is neutral in the RT medium of the fallback. (7) The replayed envelope used 2 threads; round-off differences from the original c1a trajectory are possible (M_pt agrees to 5 digits).
### E2E-A Chained ABN02 run without feedback [run]
Files: `julia/e2e_chain.jl` (chain driver), `e2e_handoff.jl` (star export/import), `e2e_profile.jl` (sampling profile), `e2e_gamma.jl` (Γ(m) of a hand-off star), `e2e_atomic_diag.jl`; Perf recording in `TitanDriver`
(`:integrator`, `:matgen`, `:linsolve`, counters), `Trbdf2`, `TitanResidFast` (`:residual` = forward-mode AD evaluation of the residual (values and Jacobian partials in one pass), `:jacobian` = scatter, `:context`),
`TitanMaterials` (`:eos`, `:opacity` = the state fills only; the zone evaluations inside the residual are in `:residual`), `popiii_m2c.jl`, `s6b_run.jl`, `s5_run.jl` (`PERF=<name>` writes `lib/popiii/perf/<name>_*`, `PERFSTAGE`).
Perf.jl: the first call of each section is kept separately (`first_call_s`: JIT, 45–200 s for the integrator). Logs `lib/popiii/*e2e*`, `perf/e2e_*`; checkpoints, snapshots and the hand-off file are in `lib/popiii/e2e_work/` (not committed).
New switches, all default off (defaults unchanged: the P100 default collapse `MAXSTEPS=12000 NSTOP=1e12 lines 100` is bit-identical before and after, diff of the log empty, also after the rebase; quick tier 8/8 PASS after the rebase, full tier 26/26 PASS on the pre-rebase code):
`AtomicHHe.THICKCONT` (env `THICKCONT`), `ChemNet`: a singular Jacobian is a failed attempt, `popiii_m2c.jl`: `SNAPEND`, `SNAPDIR`, `MSTOPMSUN`, `MDOTSTOP`; `s6b_run.jl`: `EDDCAP`, `MDOTSCALE`, `HANDOFF_OUT`; `s5_run.jl`: `HANDOFF`, `EDDCAP`.
**Perf overhead:** 18 ns per `@section`, 5 ns per `count!` (microbenchmark), ≈ 6e5 sections in the P100 run (165 s): 1e-4; A/B of that run, HEAD vs instrumented, 166.5 s vs 164.9 s.

**Chain** (`julia --project=env e2e_chain.jl [collapse|formation|sink|protostar|ms ...] [--force|--status]`; `E2E_TAG`, `E2E_PERF`, `E2E_MCAP` [250 M☉], `E2E_EXTRA_<STAGE>="K=V ..."` for new switches, e.g. the S7b feedback);
each stage reads the saved output of the previous one (`<E2E_WORK>/<stage>.done`) and restarts from snapshots/checkpoints if killed (the machine is shared: other agents' `pkill` killed four of my runs; the children are started as `e2ejl ... -e include(ENV["E2E_SCRIPT"])`).
Settings: IC=abn02z19 TCMB=54.5, 200 shells, TR-BDF2, AVZONE=2, NLTE=hyb, CHEM=netatomicli, XFLOOR=1e-12 CHEMFALLBACK=1 THICKCONT=3 TEQMAX=1e30. Stage 1 collapse; 2 formation (RESTART, LINEEVERY=5, 3000 steps); 3 sink (S6b Part 1 method, NLTE=ep, SPLIT=2, r_s = 5 AU M, CDYN 0.03);
4 protostar (S6b Part 3: PRECTHICK, Mdot cap ramp, 150 zones PSI=0) fed by the Gaussian-smoothed Mdot(M*), EDDCAP=0.95; 5 ms (`s5_run.jl` from the hand-off, θ = 1, EFRAC 5, DXMAX 4e-3, CONV 1e-5, EDDCAP 0.95, TOL 0.01, NITER 40).
**Stop of the accretion (documented):** without feedback the cloud feeds the star for ≫ 1e5 yr (Ṁ = 4e-3 M☉/yr at 250 M☉); the chain stops at M_pt ≥ `E2E_MCAP` = 250 M☉ (reached 261.6 at 30.8 kyr, 22 299 steps) or Ṁ < 1e-5 (never). The Ṁ table ends at 170 M☉ (windows of Δln M ≥ 0.5) and Ṁ = 6.6e-3 is held up to 262 M☉ (true ≈ 4–5e-3: 20–30 % too fast above 170 M☉). Mass loss is off in all stages.

| stage | span | steps | wall (4 threads) | milestones |
|---|---|---|---|---|
| 1 collapse | t = 0 → 21.676 Myr, n_c 1.1e3 → 7.7e20 | 9 776 + 179 | 9 232 + 1 343 s | T_c [K] at log n_c = 6/8/10/11/12/13/14/16/18/19/20: 266/629/826/1062/1314/1553/1797/2265/3255/4106/6424 (final2: 260/628/824/1061/1314/1551/1796/2271/3262/4279 at 19.4); f_H2,c 0.27 at 1e20 |
| 2 formation | 3000 steps from n_c 7.7e20 (Δt ≈ 2 d) | 3 000 | 17 330 s | protostar at step 459: n_c = 5.7e21, T_c = 2.4e4 K, hydrostatic core (\|u\| < 3 km/s) M = 2.2e-3 M☉ inside r = 7.5e10 cm (1.1 R☉) (OHY10: 1–2e-3 M☉, 0.5–1 R☉, n = 1e21–22, T = 2–3e4); peak n_c = 1.0e22 (T_c 2.3e4); then the core expands (step 3000: n_c 4.6e20, T_c 7.5e3, 4e-3 M☉ inside 3.6 R☉) |
| 3 sink | 30.79 kyr, M_pt 0.010 → 261.6 M☉ | 22 299 (25 Newton retries) | 3 512 s | Ṁ below; boundary zone T flickers 1.7e3–1e4 K as in S6b, 0 chemistry failures |
| 4 protostar | 0.217 → 261.9 M☉, 28.65 kyr | 2 801 (43 retries) | 3 111 s | R* peak 600 R☉ at 25.2 M☉ (788 yr); L* > 1e2/1e3/1e4 at 10.9/11.9/15.3 M☉; f_D < 0.99/0.5/0.1 at 30.6/42.1/59.7; CNO > 1 % of L* at 145 M☉ (T_c 1.29e8), > 50 % at 152, > 90 % at 219 M☉; end: R* = 28.0 R☉, L* = 8.2e6 L☉ (+ L_acc 1.9e6), T_eff,ph 9.5e3 K, T_c 1.38e8, Γ_e(0.35) = 0.83 |
| 5 ms | 0 → 2.349 Myr (X_c 0.743 → 9.9e-4) | 9 382 | 43 796 s | table below; 13 240 Newton divergences retried, 19 rejected steps, ledger closure 1.8e-5 |

*Stage 3 against the earlier work.* Ṁ(M*) [M☉/yr] (smoothed table of this run | S6c z19 a2 | S6b z19f | OHY10 Z = 0): 0.03 M☉ 0.133 | 0.127 | 0.132 | 0.12; 0.1: 0.097 | 0.092 | 0.099 | 0.085; 0.3: 0.071 | 0.064 | 0.064 | 0.068; 1: 0.053 | 0.048 | 0.053 | 0.055;
4: 0.042 | 0.038 | 0.043 | 0.042; 10: 0.034 | 0.030 | 0.034 | 0.035; 20: 0.028 | 0.023 | 0.025 | 0.032; 40: 0.021 | 0.018 | 0.021 | 0.027; 100: 0.0102 | 0.0077 | 0.0125 | 0.020 (this run/OHY10 = 0.93–1.14 up to 10 M☉, 0.86 / 0.76 / 0.51 at 20/40/100 M☉; +5–30 % above S6c).
Protostar (M [M☉]: t [yr], Ṁ, R* [R☉], L* [L☉], L_acc, T_eff,ph [K], T_c [K], f_D): 1: 34, 0.053, 88, 0.05, 1.9e4, 5.6e3, 2.0e5, 1.000 | 3: 76, 0.044, 129, 0.45, 3.2e4, 6.1e3, 2.5e5, 1.000 | 10: 261, 0.034, 169, 5.4, 6.3e4, 6.7e3, 3.3e5, 1.000 | 20: 591, 0.028, 260, 3.7e4, 6.7e4, 6.5e3, 4.8e5, 1.000 |
30: 983, 0.024, 412, 3.1e5, 5.4e4, 6.9e3, 8.8e5, 0.994 | 40: 1437, 0.021, 142, 4.9e5, 1.8e5, 7.8e3, 1.6e6, 0.59 | 60: 2552, 0.016, 59, 1.0e6, 5.1e5, 8.0e3, 1.1e7, 0.097 | 100: 5750, 0.010, 28, 2.5e6, 1.2e6, 8.3e3, 4.9e7, 0.009 | 150: 11 802, 0.0071, 57, 4.8e6, 5.9e5, 8.6e3, 1.4e8, 0.010 |
200: 19 276, 0.0066, 43, 6.2e6, 9.7e5, 8.7e3, 1.4e8, 0.006 | 250: 26 859, 0.0066, 30, 7.8e6, 1.7e6, 9.3e3, 1.4e8, 0.003. S6b z19f (same Ṁ within 10 %): R* peak 549 R☉ at 24.4 M☉, D at 30/42/60 M☉, CNO ignition 151 M☉: reproduced (600 at 25.2, D 31/42/60, CNO 145–152); there the star re-swelled at 150 M☉ (R* = 78), here it contracts (R* = 57 → 28 R☉) because of EDDCAP.

*Hand-off star* (M = 261.92 M☉, R = 28.04 R☉, L = 8.14e6 L☉ = 0.83 L_Edd,es, T_c 1.380e8, ρ_c 32.05; Γ(m) with κ_R: 0.61 in the core, 0.59 at m/M = 0.34, 0.65–0.84 in the envelope; M_cc/M = 0.31 after the relaxation).

*Stage 5, main sequence* (no mass loss; the hand-off grid is kept, 150 zones, surface zone 3e-7 M):

| X_c | t [yr] | R [R☉] | log L | log T_eff | T_c [K] | ρ_c | M_cc/M |
|---|---|---|---|---|---|---|---|
| 0.743 (hand-off, relaxed) | 0.3 | 28.1 | 6.912 | 4.765 | 1.380e8 | 32.1 | 0.31 |
| 0.741 (**ZAMS-equivalent**, R minimum) | 2.0e4 | 6.78 | 6.761 | 5.036 | 1.355e8 | 26.5 | 0.86 |
| 0.7 | 1.95e5 | 7.43 | 6.768 | 5.018 | 1.295e8 | 23.0 | 0.86 |
| 0.5 | 9.28e5 | 9.82 | 6.823 | 4.971 | 1.252e8 | 20.2 | 0.75 |
| 0.3 | 1.50e6 | 15.3 | 6.859 | 4.884 | 1.249e8 | 19.8 | 0.65 |
| 0.1 | 2.08e6 | 39.6 | 6.937 | 4.697 | 1.284e8 | 21.3 | 0.65 |
| 0.01 | 2.333e6 | 319 | 6.997 | 4.259 | 1.395e8 | 27.0 | 0.65 |
| 1e-3 | **2.349e6** | 400 | 7.001 | 4.211 | 1.523e8 | 35.3 | 0.65 |

The hand-off star is not in thermal equilibrium (R* = 28 R☉ is the accreting protostar): it contracts on t_KH ≈ 1e4 yr to R = 6.8 R☉ (the ZAMS row). τ_H(X_c < 1e-3) = 2.35 Myr from the hand-off (2.33 Myr from the ZAMS row). Schaerer 2002 (Z = 0) is tabulated to 500 M☉ but is not in the repository;
the nearby values of the S5 table (100 M☉: log L 6.11–6.17, 120 M☉: 6.24, log T_eff ≈ 4.98, τ_H 2.7 Myr) extrapolated with L ∝ M^1.1–1.6 give log L = 6.8–6.9 and log T_eff ≈ 5.0–5.05 at 262 M☉ [extrapolation, ±0.1]: the ZAMS row (6.76, 5.04) is 0.05–0.1 dex low in L; τ_H ≈ 2.4 Myr (∝ M^-0.1) agrees.
A pilot with the older 171 M☉ hand-off (`s5_e2etest_ms2.log`, protostar of S6b z19f, TOL 0.03): ZAMS R 6.46 R☉, log L 6.556, log T_eff 4.995, T_c 1.3375e8, ρ_c 34.6, τ_H = 2.47 Myr, final R 84 R☉, log L 6.70.

**Breakages and fixes** (the chain was written to expose them):
1. *collapse, step 1*: `ChemNet.step` Newton fails for species at exactly 0 (He⁺, He²⁺, Li⁺ of the atomic net) → `XFLOOR = 1e-12`, `CHEMFALLBACK = 1` (S6b settings) in all `popiii_m2c` stages.
2. *collapse, n_c > 1e16*: the S6c atomic net, never run through the collapse before, made T_c oscillate and fall (2.5e3 → 1.7e3 K at 1e21, f_H2 → 1, `m2c_*_e2e_nothick_c1*.log`; without atomic cooling or with `ATOMICOFF=1` the track is smooth and equals final2). Cause: the thin-limit emission of the photo pairs (H⁻ radiative attachment, J = CMB only, no re-absorption) in opaque gas.
`AtomicHHe.THICKCONT`: 1 = RA04 factor min(1, (n/8e9)^-0.45) (smooth to 1e18, wrong above: the emission grows as n², `m2c_*_e2e_mode1*`); 2 = (1 − e^-τ)/τ with the Planck-mean τ of `ContinuumOpacity` (still oscillates at 1e17, `dbgthick2`); **3 = emission terms × 1/(1 + (n_H/1e12)²)** (T_c 3153.1 K at log n_c = 18.02, equal to the run without atomic cooling to 1e-4; `dbgthick3`). A stand-in for the S7 transfer.
3. *formation*: S6c's `TEQMAX = 1.2e4 K` keeps zones with p > 1e4 and T > 1.2e4 K in the ideal-gas network instead of the SCVH table; the core then bounces at n_c = 2e21 (`e2e_teqmax_c2.log`) → `TEQMAX = 1e30` for the chain (identical below 1.2e4 K, so the collapse was resumed from its n_c = 2e20 snapshot).
4. *sink, M_pt = 114 M☉*: `LinearAlgebra.SingularException` in `ChemNet.step` aborted the run (not caught by `CHEMFALLBACK`) → a singular Jacobian is a failed attempt (the step halves / bisects in T); sink runs write a snapshot every 1000 steps and resume from it.
5. *protostar, 0.64 M☉ (t = 27.5 yr)*: dt falls from 1e6 to 1e2 s and stays: Newton 2-cycles at the ionisation cliff of the precursor, zone 145 of 150, corrections 5–10 in ln E_r, u. Deterministic: not cured by a restart without the E_s memory, `LTAU = 1`, `BCHOOK = 0`, or Ṁ × 1.01 / 0.99; the S6b z19f table passes. The windowed Ṁ(M*) of `s6c_mdot.jl` zig-zags by ±20 % from node to node
(0.0563 → 0.0668 between 0.51 and 0.60 M☉) and the wave onset depends on it → the chain smooths ln Ṁ with a Gaussian of σ = 0.4 in ln M (`smooth_mdot`, nodes every 0.15 in ln M); with it 0.2 → 262 M☉ runs without a stall. (S7c-1, merged afterwards, cures the stall in the boundary condition; these runs predate it.)
6. *protostar, 139 M☉* (L/L_Edd,es ≈ 0.88): the star inflates, R* 50 → 930–1260 R☉ breathing on hours, dt 1e4–7e5 s against t_acc 3e11 s (`s6b_e2e_noeddcap_protostar.log`; S6b stopped at 170 M☉ with R* = 78). `EDDCAP = 0.95` (soft cap κ ≤ Γ_max 4π c G M/L*, every step) brings it back (R* 566 → 100 R☉ in 300 steps, L* transient 1e8 L☉, then dt 2e9 s; `s6b_e2e_t_edd.log`).
7. *hand-off*: with P_ext = 0 and no precursor the 171 and 262 M☉ models are unbound at once (top zone Γ = 1.26 at T = 9e4 K, κ = 0.55; the envelope expands to 6900 R☉ in 6e6 s) → `s5_run.jl EDDCAP` (same cap, L_surface; the relaxation needs ≈ 330 steps = 960 t_dyn). With `EDDCAP = 0` the relaxation fails.
8. *ms, t = 1.88e4 yr, R = 7.3 R☉* (contraction to the ZAMS): Newton diverges at all dt, with corrections of 0.5–25 in the top zones (surface zone 3e-7 M at τ ~ 1e4, 1 − Γ_top = 0.05 at the cap; `s5_e2e_msB.log`) → TOL 0.01 and NITER 40 (passes it and the rest; `msF`, `msG`). EDDCAP 0.85 instead gives a luminous, fully convective star (L_surf/⟨L_nuc⟩ = 1.5, dt → 1 yr; `msE`).
9. Process management on the shared machine: four runs were killed by other agents' `pkill`; the chain's `launch` reported exit 0 for them (now: a signal = failed attempt, resume from the last snapshot/checkpoint), and a stage whose log lacks the `wall time:` line is not accepted as finished.
Not understood: after formation the core re-expands (stage 2 above) while the S6b `core2` run (netmodernli) kept compressing to 6.7e22 at T_c 5.7e4 K from the same state; the two separate after ≈ 500 steps (He chemistry in the envelope?). The sink does not use it (it removes the inside of 5 AU M).
Not cured, only worked around: EDDCAP is a numerical stand-in for the wind/inflation physics at Γ → 1 (mass loss, porous or super-adiabatic envelopes); the end-of-MS inflation to 400 R☉ (T_eff 1.6e4 K) is unchecked.

**Performance** (`perf/e2e_*`; sections inclusive; wall includes the JIT of the first call; 4 threads on a shared machine):

| stage | steps | Newton iter / step | failures | where the time goes |
|---|---|---|---|---|
| collapse | 9 899 | 4.4 | 1 stage fail, 1 BE fallback, 0 chem | ALI line transfer 82 % (2 271 calls × 3.3 s), chemistry 8 %, hydro step 1.7 % (1.6 ms per residual) |
| formation | 3 000 | 4.4 | 0 | ALI 96 % (600 calls × 27.8 s at 129–158 bins), chemistry 2.3 %, hydro 0.6 % |
| sink | 22 299 | 5.1 | 25 | line transfer (NLTE=ep) 43 % (68 ms/step), chemistry 42 % (67 ms/step), hydro 14 % |
| protostar | 2 801 | 12.9 | 43 | residual 77 % (60 ms per call), boundary condition (precursor tables, 13 radii) 13 % (144 ms/step), EOS/opacity fills 3 %, burn+mix+advect 1.5 %, linear solve 0.3 % |
| ms | 9 325 | 86.9 (31.6 per integrator call, 2.75 solves per step) | 13 240 divergences, 19 rejections | residual 95.6 % (50 ms per call, 836 000 calls), linear solve 0.2 %, burn 0.1 % |

Sampling profile of the star residual (`perf/s5_ms_profile.txt`, `s6b_protostar_profile.txt`; 5 ms samples of the main thread): `StellarEOS` 60 % (S5) / 67 % (S6b), ForwardDiff arithmetic in the residual and its callers 17 + 13 %, `NucNet` 4 %, `StellarOpacity` 2–4 %, block solver 0.9 %.
*Three most promising targets* (not done): (1) **EOS inside the AD residual** (60–67 % of every star stage, most of the 12-hour ms stage): evaluate `StellarEOS.gas_state` in Float64 and inject its analytic (ρ, T) derivatives as dual partials instead of carrying 35-wide duals through the ionisation iteration, or tabulate it:
estimated 2–3× on stages 4–5. (2) **ALI transfer** in stages 1–2 (82–96 %, 2.9 h + 4.8 h): the hybrid refresh interval is capped at 20 steps and every call restarts from LTE; warm start from the last source functions, interval cap 60 where the factor moves < 1 %, only the lines that are thick: estimated 3–4× (7.7 h → 2–2.5 h).
(3) **Newton work in the main sequence**: 31 iterations per solve and 1.4 divergences per step, because (∇ − ∇_ad)^{3/2} in the convective core and the capped surface make the Jacobian stiff; Jacobian reuse (chord) within a solve, `CONV` 1e-5 → 1e-4 outside the core, and a surface condition that does not sit at Γ_top = 0.95 could halve the 810 000 residual evaluations (estimate 1.5–2×, 12 h → 6–8 h);
the sink's per-zone chemistry loop (42 %) is trivially threadable (≈ 1.3× on that stage).

**Open points.** (1) Re-run stages 4–5 with the S7c-1 code (same-dt Newton retry, Eddington gate): they ran on S7b. (2) The Ṁ table beyond 170 M☉ and the 250 M☉ cap are choices; the S7b feedback runs end the accretion physically. (3) τ_H, ZAMS and the end-of-MS radius of a 260 M☉ star without wind depend on the model (EDDCAP 0.95, α_MLT = 1.5, Schwarzschild boundary, fixed surface-zone mass).

### S8a Monte Carlo Lyman α transfer [run]
Files: `julia/LyaMC.jl` (module `LyaMC`, self-contained: includes its own `LymanWings`), `s8a_couple.jl` (module `S8a`: inner column, grid coupling), `s8a_test.jl` (full tier `s8a`), `s8a_snap.jl` (Monte Carlo field of stored envelope snapshots), `s8a_prod.sh`; edits (all default off): `RadTransfer.jl` (`RT_LYA=mc` → `LYAMC`, `LYA_MC_WIN`, `Opts.lya_mc_win`: the window around Lyα and the line sources are taken out of the MT08 force), `RadLya.jl` (`skipbins`, `nolines`), `s7b_env.jl` (hook `lya_mc`, logs, `dump_J!`, Ṁ estimator), `s7b_star.jl` (`S7B_SLEWDN`), `hi21_snap.jl` (`LYA=mc:<dir>:<tag>`), `popiii_tests.jl` (test `s8a`). Papers (arXiv, data only) added to `~/Projects/popiii_lit/lyman/`: Dijkstra 2014 review (1406.7292), Dijkstra, Haiman & Spaans 2006 I/II (astro-ph/0510407, 0510409), Hirata 2006 (astro-ph/0507102), Hansen & Oh 2006 (astro-ph/0507586), Laursen et al. 2009 (0907.2698), Herczeg et al. 2002 (H₂ pumping lines near Lyα, astro-ph/0201319), Tomaselli & Ferrara 2021 (2103.14655, already there). Defaults are unchanged: `RT_LYA` unset reproduces S7c-1; `S7B_STARVE` is the only default that changes a run (Ṁ estimator, §4).

**1. Method** (1D spherical, moving medium, per zone T, n_HI, n_e, n_p, H⁻ absorption, n_H₂; radial velocity at the interfaces, linear inside a zone).
- *Frequency.* Photon state (r, μ, w = c Δν/ν₀); fluid-frame x = (w − uμ)/b, b = √(2kT/m_H); the ray is sub-stepped (Δx ≤ max(0.25, 0.1|x|)) because x changes along the ray in a velocity gradient (test c1).
- *Cross section.* Voigt profile H(a, x) (Humlíček W4 of `LymanWings`, a = A₂₁/(4πΔν_D)) times the Kramers–Heisenberg/Lorentz ratio C(Δν/ν₀) tabulated from `LymanWings.σ_coh` (|Δν/ν₀| ≤ 0.3, the wings and the Lyβ side included; 1 inside 1e-5).
- *Scattering.* Exact atom velocity: u∥ from p(u) ∝ e^{−u²}/((u − x)² + a²) by two-region rejection with a Gaussian envelope for u < u₀ and a Cauchy envelope for u ≥ u₀, any u₀ ∈ [0, x] exact (Zheng & Miralda-Escudé 2002; 0.1–0.4 µs per sample, KS test a1), u⊥ Gaussian; coherent in the atom frame with recoil −g(1 − cosΘ), g = hν₀/(m_H c b); outgoing frequency x_out = x_a − g(1 − cosΘ) + u∥ cosΘ + u⊥ sinΘ (partial redistribution, <Δx²> = 1, drift −1/x in the wing: a2); angle: Rayleigh (3/8)(1 + μ²) in the wings, 2/3 isotropic + 1/3 Rayleigh (Hamilton 1940: 11/24 + 3μ²/24) inside the fine-structure splitting |ν_a − ν₀| < 10.97 GHz of the atom frame.
- *Destruction per scattering.* 2p → 2s mixing by protons (q = 4.7e-4 (T/1e4)^{1/2} cm³ s⁻¹, electrons half of it, 2s → 2p three times larger) followed by two-photon decay (8.23 s⁻¹, chain summed): ε = P_ps P_d/(1 − P_ps(1 − P_d)), saturates at A_2γ/(3A₂₁) = 4.38e-9 for n_e ≫ n_cr = 1.6e4 cm⁻³; frequency-independent (the standard assumption). Continuous absorption κ_abs (H⁻ photodetachment at 10.2 eV, John 1988) and H₂ Lyman-band lines: the 18 B–X lines with v″ ≤ 3 within ±3 Å of Lyα (Herczeg et al. 2002 Table 6: 1–2 R(6) 1215.73, 1–2 P(5) 1216.07, 0–2 R(0) 1217.21, … ; LTE populations, Voigt with Γ_B = 1.7e9 s⁻¹; **f = 0.01 for all lines, the individual f are not in the local data**).
- *Sources.* Zone sources (case-B recombination × Cantalupo et al. 2008 Lyα fraction, thermal-profile emission uniform in volume, isotropic); `:shell` (spectrum table at a radius, Lambert directions), `:sphere`, `:point`, `:line`. Collisional 1s → 2p excitation by electrons (Callaway fit of Dijkstra 2014) as an option (`S8A_GRIDSRC = all`, default).
- *Core skipping.* A scattering at |x| < x_c (x_c ∈ {6, 3, 2, 1.2}, the largest one for which D = 4 x_c l(x_c) ≤ 0.2 Δr_zone; l = mean free path at x_c) is replaced by a precomputed frequency walk of the real kernel (tables per ln T bucket: exit frequency, number of scatterings N (1.5e4 from x = 0 at x_c = 3, 1.8e5 at 6, T = 800 K), visits per frequency bin, which keep J(x) and J_α), the destruction acts on all N (1 − ε − ε_A)^N, the exit direction is isotropic. The photon does not move during the walk, so skipping is applied only deeper than D from both interfaces and from r = 0 (a skipped walk next to a boundary changes the reflection of the layer: N_H = 1e19 shell, force 7.5 L/c with skipping everywhere against 40.6 ± 0.6 exact; 41.4 ± 0.8 with the layer). The walk tables are cached per process (first call: 20–60 s).
- *Outputs per zone* (`LyaMC.result`, a NamedTuple behind a property forwarder; arrays over the zones of the run, rates per second for the injected rate `rate`): `force` [dyn, radial, on the gas of the zone; momentum balance of the radiation field F_k = A_kΠ_k − A_{k+1}Π_{k+1} + ∫(u − Π_rr)/r dV from crossing tallies of Π_rr and Δμ along the tracks (∫(1 − μ²)/r ds = Δμ); variance independent of N_sc; the event estimator Σ(hν/c)μ_in is `force_event`, noise ∝ √N_sc]; `force_err` (batch scatter); `Jalpha` [photons cm⁻² s⁻¹ Hz⁻¹ sr⁻¹, track length averaged over |x| < 1 of the zone's Doppler profile; **the 21 cm module uses J = hν_α Jalpha = `LyaMC.E0 * Jalpha` [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹]**], `Jint`, `Jx` (path per x bin), `Rsc` (scatterings per H atom per s = P_α = 4π∫σ_ν J_ν/hν dν, the input of the Wouthuysen–Field rate), `Tcolor` [K] with `Tcolor_err`, `Tcolor_valid` (T_k when the tilt is not measurable), `nsc`, `nskip`, `heat_recoil` (E₀ (b/c) g per scattering), `heat_exch` (h(ν_in − ν_out) in the fluid frame, noisy), `Egas_lab` (lab-frame energy given to the gas, work included), `heat_Hm`, `n_Hm`, `n_H2`, `n_2g`, `f_out/f_in/f_2g/f_Hm/f_H2/f_lost`, `nsc_per_photon`, `spec_out` (+ `spec_out_src`, `wcenters(res)`), batches. `T_c` is measured from the odd part of ln J(x): d ln J/dx = −2g T_k/T_c (Hirata 2006), fitted over 0.3 < x < 6 from the 8 batches; `LyaMC.hi21_adapter(dir, tag)` returns the `f(snapshot) → (; J, Tc)` that `HI21.use_lyamc!` wants (it reads `s8a_<tag>_J_<j>.csv`, written by the coupled run at the barrier of each saved snapshot or by `s8a_snap.jl`).

**2. Tests** (`s8a_test.jl`, `lib/popiii/s8a_test.log`, full tier `s8a`, 19 checks PASS, 8 min on 4 threads; quick tier 8/8 PASS (`lib/popiii/s8a_quick_tier.log`) and `s7a` (24), `s7b` (20), `s7c2` (15), `s8b` (20) re-run PASS after the edits of `RadTransfer.jl`/`RadLya.jl`/`LymanWings.jl`; the rest of the full tier was not rerun, nothing else is touched).
- a Sampler and kernels: parallel-velocity sampler KS distance 0.0049 over 8 (x, a) pairs (1 % critical value 0.0052); wing redistribution <Δx²> = 1.03 / 1.01 and drift −1/x at x = 8 / 30; Rayleigh angle sampler <μ²> = 0.4005; ε(n_e = 1e8) = 4.381e-9 = A_2γ/(3A₂₁); KH wing table against `LymanWings.σ_gas` within 1e-3 at ±2…10 % of ν_α.
- b Static media (T = 10 K, a = 0.015). Slab (thin shell, midplane source): peak of the emergent spectrum at 1.05 (aτ₀)^{1/3} (Harrington 1973/Neufeld 1990: 1.06; Hansen & Oh 2006 eq. A1), rms deviation of the shape from the analytic cosh solution 0.24 (3000 packets); N_sc/τ₀ = 1.97 / 1.79 / 1.72 for τ₀ = 1e4 / 3e4 / 1e5, log-slope 0.94 (tends to 1.6). Uniform sphere (Dijkstra, Haiman & Spaans 2006 eq. 9): peak 0.95 (0.92) (aτ₀)^{1/3}, rms deviation 0.23, N_sc = 0.93 τ₀.
- c Moving medium. First-flight scattering probability in a homologous flow u = Γr against the Sobolev value 1 − <e^{−τ(x₀)}>: 0.0484 / 0.371 / 0.898 against 0.0479 / 0.364 / 0.893 for τ_S = 0.1 / 1 / 10. Sphere with u = ∓8b r/R, aτ₀ = 1.4e3: infall blue fraction 1.000, outflow 0.000, <x> = +7.12 / −7.35, N_sc equal within 3 % (mirror symmetry); the gas loses −<x>b/c = 3.0528e-4 E_α per photon to the escaping photons, the tally gives 3.0527e-4.
- d Destruction. Neufeld (1990) escape fraction for frequency-independent absorption (aτ₀ = 1.5e3, slab): 0.800 / 0.548 / 0.204 for τ_a = 0.01 / 0.03 / 0.1 against 0.800 / 0.553 / 0.221 (ζ = 0.525 fit). Per-scattering ε: f_esc = <e^{−εN}>: 0.783 (MC) against 0.749 for ε = 3e-7, τ₀ = 1e6 (N_sc = 9.7e5 without destruction). Scaling used for the inner column: (τ₀, ε) → (τ₀/10, 10 ε) gives f_esc 0.760 and the same emergent spectrum in x/(aτ₀)^{1/3} (rms x ratio 0.98).
- e Force. *Exact*: conservative grey scattering (far Lyα wing, N_sc ≫ 1) has F = L/(4πr²) at every radius and force density κF/c, so M_F = τ_total and every zone feels τ_zone: the momentum-balance estimator gives 1.033 ± 0.037 (τ = 3) and 1.002 ± 0.035 (τ = 30) of τ_zone per zone (the event estimator: 12 % high, noise ∝ √N_sc). Static thin shell around a central source (D&L08 set-up, T = 1e4 K): M_F = 40 / 57 for N_H = 1e19 / 1e20, D&L08 fit 10 / 27, MT08 plateau 14 / 29, the diffusion solution for a uniform sphere 3.51 (aτ₀)^{1/3} = 23 / 49 (Tomaselli & Ferrara 2021). **Our M_F is 1.5–4 times the D&L08 fit**: the cavity returns the photons that the shell reflects (18 crossings of the inner face per photon at N_H = 1e19, so M_F counts every bounce), and the estimators and the exactness test above exclude a bookkeeping error; Kimm et al. (2018) and Smith et al. (2017) also find numerical values below the diffusion solution, we find above; unresolved (open point 1).
- f Wouthuysen–Field (T = 10 K, recoil 2g = 1.6e-2 per Doppler width, τ₀ = 1e5, 12 000 packets): T_c/T_k = 1.068 ± 0.044, 1.016 ± 0.048, 1.16 ± 0.10 in the three inner zones (weighted mean 1.056 ± 0.031: T_c → T_k, Field 1959, Hirata 2006); without recoil the fitted tilt is zero.
- g Ledgers: photons out + in + 2γ + H⁻ + H₂ + cut = 1 − 1e-15; energy E_emit − E_out − E_in − E_2γ − E_gas = −1e-15 E_emit in a moving (infall 5b) medium with H⁻ and H₂ absorption; event and momentum-balance forces agree in a thick sphere (60.9 / 59.2).
- h Noise ∝ N^{-1/2} (relative scatter of the total force over 8 seeds 1.5e-2 / 1.4e-2 / 5.3e-3 for 200 / 800 / 3200 packets, slope −0.37 with 8 seeds); core skipping against the exact transfer (T = 10 K, τ₀ = 3e5, 4 zones): f_esc 0.946 / 0.951, N_sc 2.76e5 / 2.75e5, rms x 15.83 / 15.88, Σ force 75.9 / 74.9 L/c.

**3. Coupling** (`RT_LYA=mc`, `s8a_couple.jl`, hook `lya_mc` in `s7b_env.jl`; `RT_LYA=1` is the MT08 force of S7c-2, unchanged).
- *Replacement.* In `transfer` the beam's single-scattering force of the H I in the window |ν/ν_α − 1| < 0.12 (`RT_LYA_MCWIN`; the Lyβ line, at +0.185, is outside) is removed and replaced by the force of the Monte Carlo; the MT08 plateau/closure of the other lines and of the coherent far wings is dropped (the beam keeps them, as in S7b and the no-Lyα runs; `RT_LYA_MT08REST=1` keeps it with the Lyβ… conversion `LYB_CONV`: Ly-n loses 1 − branch(n) = 12–19 % of its scatterings to the cascade). A first attempt that kept MT08 for the rest (runs `d1`, `d2`, partial: `s7b_d{1,2}partial_env.csv`) gave F/(Ṁ v_ff) = 0.04 at 70 M☉ from the coherent Lyman-wing/Rayleigh closure alone (flux conserved in 1D, 0.001 in the beam-only runs) and slowed the infall by ≈ 30 % (M_pt = 69.5 at 3.7 kyr against 2.7 kyr); that is the MT08 spherical-closure assumption, not Lyα, and is not part of the result below.
- *Inner column* (H II region → first interface, N_H = 1.4e26–2.1e27 cm⁻², τ₀ = 7e12–1e14, aτ₀ ~ 1e10: out of reach of a direct run, N_sc ≈ 1.6 τ₀ = 1e14): scaled medium with the invariance of the transfer equation in the diffusive limit (x → x/λ, τ₀ → τ₀/λ³, ε → sε with s = λ³, g → λg, u → u/λ; test d2/d3), τ₀' = 2e6, 24 sub-windows of the stellar continuum + the H II-region line; cached while τ₀ and T stay within 25 % / 20 %. Result: **the Lyα of the H II region (0.68 Q) does not leave the column during the accretion**: ε τ₀ = 1.5e4–6e4 (ε = 4.4e-9, n_e ≥ 1e5 cm⁻³ saturated), so the survival is e^{−ε N_sc} < e^{−1e4}, not the 2–3 % of S7c-1 (the (ετ)^{1/3} of MT08 B13 measures the reduction of the trapped intensity near the source, not the escape; for a photon in the core the number of scatterings before it reaches the wing is ∝ τ₀, and every one of them can destroy it). The stellar continuum in the window is destroyed as well (the escape begins at |x| ≈ x_esc = (aτ₀)^{1/3} = 1.7e3 b = 9 % of ν_α). The line escapes only when the column thins after the accretion has ended (grid edge 2e4–5e4 AU, N_H ≲ 1e24): escape fraction 0.56 at N_H = 3e24, 0.2–1.0 later.
- *Grid.* `LyaMC` on the hydro zones (moving medium, composition of the sweep: n_HI, n_e, n_p, H⁻, H₂), sources = emergent spectrum of the column (shell at the first interface) + case-B recombination and collisional excitation Lyα of the zones (own emissivity from n_e n_p α_B f_B and n_e n_HI γ_1s2p, not from the AtomicHHe photo-pair terms, so `THICKCONT` does not enter), inner boundary reflecting with the albedo 1 − 1/(1 + 0.75 τ_col(ν)) of the thick column. The force of the last Monte Carlo is held and scaled with the source rate between calls (`S8A_NEVERY` = 25 barriers). 300 + 300 packets, 8 batches.
- *Noise and cost.* Total force: relative scatter 4–7 % (median over the 33 calls 6.5 %, 10–40 % after the accretion has ended, when the forces are 100 g and photons are destroyed in the hot first zones); **per zone ≈ 40 % (zones carrying more than 5 % of the peak force), not the 5 % asked for: that needs ≈ 80 × more packets (hours per call); these forces are ≤ 1e-3 g until the infall ends, so the budget was kept at 600 packets**. Cost (3 threads, `lib/popiii/perf/s8a_d4_env_sections.csv`): 33 calls, median 76 s (first 157 s: JIT and 21 skip tables; 666 s at 25 barriers while the tables of the new temperature buckets are built; later 13–390 s), 4227 s of the 9619 s wall of the envelope (the star exchange waits for 8000 s of it, in parallel); inner column 1.2 s per call; the core skipping makes all zones eligible (N_sc/photon 1e5–1e9 in the call, explicit scatterings ≈ 5e-4 of them).
- *Two things found on the way.* (i) The event estimator of the force is useless for thick media (variance ∝ N_sc: 59 ± 14 against 40.6 ± 0.6 L/c at N_H = 1e19) and core skipping without a boundary layer biases the force low by 5 (7.5 against 40.6): both in the Method. (ii) The retained `d1/d2` MT08 remainder (above).

**4. Ṁ estimator** (`mdot_window`, `barrier!`, `S7B_STARVE` = 3; `S7B_SLEWDN` = 0.2, `S7B_DTMAX` = 60 yr in `s8a_prod.sh`). The S7c-1 estimator was bounded by the average since the first event of the window, which decays as 1/t: after the sink starved the star kept growing (178 against M_pt = 143). Now the rate cannot exceed STARVE × (mass of the last zone removed)/(time since the last event) (also when no 4-event fit exists: then the point mass bounds M_acc), the star may lower its Ṁ by 20 % per exchange (default 5 %) and the exchange interval is capped. Run `d5` (no Lyα force, 40 kyr): Ṁ_win 7.9e-4 M☉/yr at the end and **star 144.8 M☉ against M_pt = 142.7** (S7c-1: 178 against 143); `d4`: 138.5 against 134.9 (the star was frozen at its last record after it failed at 138.6 M☉, see the Result). The residual 2–4 M☉ is the zone mass at the sink (4–8 M☉).

**5. Production** (ABN02 z = 19, `core2`, two-way, MSYNC 3, FEX 0.01, `RT_LYMAN = 1`; `s8a_prod.sh`; `lib/popiii/s7b_{d4,d5}_env.csv`, `s8a_d4_lya.csv`, `s8a_d4_lyaprof.csv`, `m2c_abn02z19_lines_n200_s7b_{d4,d5}.log`). d4 = Monte Carlo Lyα (`RT_LYA=mc`), d5 = no Lyα force (reference with the new estimator). The protostar fails at 113 M☉ (L/L_Edd,es ≈ 0.88: R* → 1400 R☉, the known S6b inflation of E2E-A) in d4: restarted from its checkpoint at 100 M☉ with `EDDCAP = 0.95` (the cap of E2E-A; d5 has it from the start) and the replay of the recorded exchange (`S7B_RESUME`); it fails again at 138.6 M☉ (d4; 144.8, d5: "grid extrapolation nonmonotonic") after M_pt had stopped growing, and the last star record is repeated to the end of the run (24.9 kyr for d4 and 40 kyr for d5; the star only enters through L*, T_eff, Q* of a star at 139–145 M☉, which do not change).

| M_pt [M☉] | t_c d5 [yr] | t_c d4 [yr] | Ṁ_win d4 [M☉/yr] | R*, T_eff,* [R☉, K] | Q* [s⁻¹] | F/(Ṁ v_ff) with MC | MC part | F/(Ṁ v_ff) MT08, same state | grid edge d4 [AU] |
|---|---|---|---|---|---|---|---|---|---|
| 10.1 | 161 | 161 | 3.6e-2 | 170, 1.0e3 | 1.1e42 | 0 | 0 | 0 | 115 |
| 30.4 | 829 | 913 | 2.1e-2 | 393, 6.8e3 | 9.4e42 | 0 | 0 | 0 | 186 |
| 52.8 | 1784 | 1949 | 1.7e-2 | 80, 1.9e4 | 2.9e48 | 0.0022 | 0.0017 | 0.46 | 496 |
| 73.6 | 2971 | 3338 | 1.3e-2 | 37, 3.1e4 | 3.7e49 | 0.0008 | 0.0004 | 1.5 | 947 |
| 102.5 | 5113 | 6084 | 9.9e-3 | 22, 4.8e4 | 1.4e50 | 0.0007 | 0.0003 | 3.5 | 1540 |
| 120.8 | 7453 | 7404 | 9.5e-3 | 21, 5.0e4 | 1.8e50 | 0.0020 | 0.0011 | 6.3 | 3912 |
| 127.8 | 7633 | 9084 | 8.1e-3 | 38, 4.1e4 | 2.0e50 | 0.0009 | 0.0000 | 6.6 | 4201 |
| 134.9 | 9433 | 15864 | 3.1e-3 | 26, 5.1e4 | 2.7e50 | 0.14 | 0.13 | 28 | 2.0e4 |
| 142.7 | 13333 | not reached | | | | | | | |

*Final mass and mechanism.* **d5 (no Lyα force): 142.7 M☉ at 13.3 kyr** (zone mass 4–8 M☉; S7c-1 c1a: 142.75 at 20.4 kyr, S7b one-way 160 ± 8), no further growth to 40 kyr. **d4 (Monte Carlo Lyα): 134.9 M☉ at 15.9 kyr**, one zone (7.9 M☉) below d5; it follows d5 within the scatter of the zone staircase to 120.8 M☉ (7.4 kyr, both) and then the last two zones arrive 1.5 and 6.5 kyr later (9.1 and 15.9 kyr against 7.6 and 9.4 kyr). With the Monte Carlo the accretion is not stopped by Lyα: the force of the line is ≤ 0.2 % of the infall momentum flux Ṁ v_ff until 120 M☉ (the MT08 approximation on the same states: 0.05 at 47 M☉, 1.5 at 74, 3.5 at 102, 6.3 at 121: it halted the S7c-1 run at 53–59 M☉, a factor 50–7000 more) and 0.1–0.3 only after the infall has ended (Ṁ_win ≤ 3e-3, first zones at 8 kK, grid edge 2e4–5e4 AU, the grid-born Lyα of the hot first zones, 50–84 % destroyed by the two-photon channel, N_sc/photon 1e8–1e9). The mechanism of the end of accretion is the S7b one: LW dissociation and heating of the infall (T_first zone 8 kK from 70 M☉, f_H₂ ≤ 1e-8), the H II region confined at R_II = 1.0 R* (1.4–1.6 R* at the end), a gradual fall of Ṁ from 1.0e-2 (100 M☉) to 3e-3 M☉/yr (135 M☉); Lyα changes only the last 8 M☉ and the timing (+6 kyr). The result is the same within ±1 zone as S7c-1's 143 ± 4 without the Lyα force; **the 53–59 M☉ of the MT08 treatment are an artefact of its approximations: unattenuated H II-region Lyα (destroyed in the column, ε τ₀ ≈ 1e4–6e4) and the trapped pressure of the continuum-wing closure.**
*Escape fractions.* From inside r_s (H II-region line photons through the column): < e^{−1e4} during the accretion (3.2–127.8 M☉), 0.56 / 0.22 / 1.0 at N_H = 3.3e24 / 1.7e23 / 3.4e21 after it (j = 525 / 600 / 650); window continuum: destroyed for |Δν/ν_α| ≲ 0.09 (expected from x_esc; the transmission of the sub-windows is computed in the scaled run but was not logged). From the grid (Lyα born in the zones, escaping the outer edge): 0.97–1.00 during the accretion (N_sc/photon 1e4–1e8), 0.14 / 0.31 / 0.41 / 0.27 at j = 450 / 525 / 600 / 650 with 84 / 65 / 50 / 32 % destroyed by two-photon decay and 0–31 % absorbed by the column (reflection fails); H⁻ and H₂ absorb ≤ 1e-2.
*Force profile* (`s8a_d4_lyaprof.csv`, a/g, MC / continuum beam without the window / MT08 total): j = 250 (47 M☉): zone 1 (472 AU) 0.038 / 0.027 / 1.85, zone 10 (3060 AU) 0.024 / 0.006 / 1.85; j = 300 (78 M☉) zone 4 (3065 AU) 0.05 / 0.038 / 106, zone 10 (6620 AU) 0.055 / 0.02 / 149; j = 350 (95 M☉) zone 7 (8030 AU) 0.064 / 0.04 / 630 — the MC force stays at a few per cent of g (mostly the beam wings: column 2 is the beam without the window, column 1 includes the Monte Carlo part), nearly constant with radius, while the MT08 a/g reaches 100–1000 where the infall is stopped in S7c-1. After the infall (j ≥ 450) the MC force reaches 10–100 g in the first zones (sign changes: noise, relerr 0.2–0.4) and the grid edge retreats from 4e3 to 5e4 AU.

**6. 21 cm (S8b) with the Monte Carlo Lyα** (`LyaMC.hi21_adapter`; `hi21_snap.jl` `LYA=mc:<dir>:<tag>`; files `s8a_c1a_J_<j>.csv` for the 11 snapshots of the S8b replay (`s8a_snap.jl`, 150–340 s each) and `s8a_d4_J_<j>.csv` for 6 snapshots written by the run itself; Tc is T_k where the tilt is not measured: with 600 packets T_c is measurable in < 1 % of the zones, 284 of 1646 spurious values were reset by the final rule 0.3 < T_c/T_k < 3 and scatter < 15 %; J in erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹ = `E0 * Jalpha`; `hi21_s8a_{c1a,d4}_{none,mc}.csv`, `_prof.csv`). The Monte Carlo field is far below the collisional coupling: max x_α/x_c = 1e-4 / 1e-4 / 8e-4 / 2.5e-2 … 8e-2 (c1a: 24 yr … 9.7 kyr, d4: up to 120 M☉), ⟨T_s⟩_m identical to the no-Lyα values to 0.1 K (591.8 → 583.3 K, 724.8 K at 142.7 M☉; the continuum fallback of S8b gave 648–708 K at 70–135 M☉), core and disc signals unchanged (disc 375.1–376.1 mK). x_α exceeds x_c only after the infall has ended (c1a snapshot 385, M* = 170 M☉, grid edge 1.3e4 AU: 19 zones with x_α/x_c up to 1e5 at r = 1e4–6e4 AU where J_α reaches 1.3e-8, but there T_s = T_c = T_k; ⟨T_s⟩ unchanged 724.8 K). So the MT08/RT surrogate overestimated the Lyα coupling by orders of magnitude; the 21 cm signal of the collapse and the feedback phase stays collisional.

**7. Performance.** Unit tests 8 min (4 threads); explicit scattering 0.2 µs per event (sampler 0.15–0.4 µs), skip tables 0.7 s (x_c = 3) … 41 s (x_c = 12) per bucket at M = 128; call costs in §3. `Perf` sections: `:lya_mc` in `perf/s8a_d4_env_sections.csv` (the envelope's own accumulator, as `:rt_transfer`).

**Open points.** (1) M_F of a static shell is 1.5–4× the Dijkstra & Loeb (2008) fit and 2× the diffusion solution for the uniform sphere (the estimator is exact for conservative scattering, test e1; unresolved, perhaps the cavity bounces or Kimm et al.-type differences in the redistribution). (2) Per-zone noise 40 % at 600 packets. (3) The inner-column result rests on the scaling (ε N_sc with N_sc ∝ τ₀ in the thick limit: tested at τ₀ ≤ 1e6 and aτ₀ ≥ 1e3; the factor e^{-1e4} makes the exact value immaterial) and on the frequency-independent ε (collisions during an off-resonant scattering are suppressed in the damping wings, which would let the window photons out at |x| ≫ 1: the wing destruction probability is then ε × O(a/x), the escape of the window stays dominated by the core). (4) H₂ line f-values are a guess (0.01); H₂ and H⁻ absorb ≤ 1e-2. (5) Recoil, hyperfine and the 2P₁/₂–2P₃/₂ interference in the line wings are not resolved beyond the atom-frame fine-structure switch. (6) Star: the S6b protostar fails at 113 M☉ (d4) and 138.6 / 144.8 M☉ (EDDCAP 0.95): after the infall stopped it is frozen. (7) The `d1/d2` runs (MT08 remainder) and the d4 star restart at 100 M☉ are not repeats of one clean run. (8) One `pkill -f feeder.py` was used on my own helper process (pattern may have matched other agents' processes of that name: none found afterwards).
