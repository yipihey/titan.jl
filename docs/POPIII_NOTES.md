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
