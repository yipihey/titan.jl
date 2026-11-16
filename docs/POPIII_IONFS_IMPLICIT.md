# IonFS in the Newton solve: H I, He I, He II and Lyman–Werner as implicit ray groups (IONFSIMP)

Branch `ionfs-implicit`, from `popiii-s6` (79d4dc1) with `ionray-fs` (35ba03d) merged. Code: `julia/IonFSImp.jl` (new),
hooks in `IonFS.jl`, `PopIIIMaterial.jl`, `TitanResidFast.jl` (`RAYIMP`), `popiii_m2c.jl`, `TitanFortran.jl` (meqn), `Trbdf2.jl`
(`NEWTONERR` diagnostic only). Everything is behind `IONFSIMP=1` (needs `IONFS=1 COUPLED=1 STARLIGHT=1`), default off.

## What it does

The IonFS ray (photon-conserving, radial, sourced by the formal solution of the table region at the first zone k₀ outside the grey
photosphere) becomes part of the hydro Newton:

- **Unknowns.** y_g(k) = ln N_g(k), the photon rate of group g entering zone k, are Newton unknowns 10–13 (`neqn` 9 → 13, `meqn` ≥ 13
  set automatically for `IONFS=1 IONFSIMP=1 COUPLED=1`). The ray row of zone k > k₀ is y_g(k) − y_g(k−1) + Δτ_g(k−1) = 0, of k₀
  y_g(k₀) = ln N_g,src; elsewhere δ = 0. A row reads only its own and the inner neighbour's unknowns, so TITAN's block structure is
  unchanged (`BlockSolver.jl` with 13×13 blocks; the ForwardDiff colouring with `ncol = 6` is untouched).
- **Coupled zones** of the ray (`cpl_rayzone`, the COUPLED `cpl_zone` with an 8-wide local Jacobian in ρ, T, x_H⁺, x_H₂, y₁…y₄,
  cached and reused as `CPLJLAG`): H photoionization `extra` N_g(1 − e^{−Δτ_g}) Δτ_H/Δτ_g, photo-heating Σ N_g(1 − e^{−Δτ_g}) E_g / V in the
  energy row, LW dissociation, He I/He II photo-rates, and Δτ_g for the ray rows, all as functions of y_g(k) and of the zone's own CK
  composition. The material EOS of a ray zone is taken at the iterate's ray too (`TX.CPLWY`).
- **He photo-rates** per atom are N_g fτ(Δτ_g) σ_g,s Δr/V (fτ = (1 − e^{−τ})/τ) with Δτ_g from the zone's own He stages. That is a fixed
  point u = expansion(rates(u)) in u = (x_H, x_He, x_He⁺). It is solved by Newton in ln u (backtracking on |F|) on the Float values, and
  differentiated by one implicit-function step w = w* − J⁻¹F(w*; z) with the residual's dual numbers (`ifsi_ufix`). Converged means
  |F| < 1e-11 in ln u (fractions floored at 1e-40, so that an iterate with x_H⁺ > 1, x_H < 0, does not stall it); where Newton stalls above
  that (x_H = 1 − x_H⁺ − … at x_H ~ 1e-8 carries round-off ~1e-16/x_H; He⁺ under a dominant He²⁺) the best of 60 iterations is used
  (~0.5 % of the solves). Zones whose He photons are negligible (N_He Δt < 1e-12 of their He) take the start-of-step fractions.
  A plain fixed-point iteration was tried first: in a photon-limited He zone its slope is (x₀ − x_He)/x₀ → 1, it did not converge, and
  the residual then depended on the previous call (CPLJCHECK ratios 0.5–1 from that noise). A looser tolerance (|F| < 1e-6) left that
  noise in the residual and the hydro Newton's median Δt over the first 80 steps fell from 380 s to 5 s (t_ifs24 against t_ifs15).
- **Split zones** inside the ray (coupled zones decoupled by `CPLBADDELTA`, table zones): Δτ from the start-of-step absorbers
  (`CPL_RAYTAU`, LW with the carrier unknown); after the Newton their chemistry uses IonFS's Float terms with the converged photons
  (`ifsi_terms_split`) and IonFS's zone-by-zone transport (`ifs_advance!`).
- **Ray extent.** k₀ … k_end, k_end two zones beyond the last zone where the ionizing photons can change x by `IONFSXMIN` in Δt or the LW
  rate changes x_H₂ by `IONFSLWMIN` (1e-2) in Δt (start-of-step sweep). Beyond k_end the LW rate is the start-of-step sweep's.
- **Replaces** `IONFSCPL` (decoupling the zones the ray changes) and `IONFSITER` (re-sweeps after each Newton update); neither is used
  with `IONFSIMP=1`. The `IONFSTOL` Δt hook acts on the converged-ray bookkeeping error.

## Groups (4)

| g | band | absorbers | why |
|---|---|---|---|
| 1 | 13.598–24.587 eV | H I | below the He I edge only H absorbs |
| 2 | 24.587–54.418 eV | H I, He I | |
| 3 | 54.418–300 eV | H I, He I, He II | |
| 4 | 11.2–13.598 eV (LW) | H₂ lines (H I Lyman lines as shielding) | |

The edges are the absorbers' thresholds, so the absorber set is fixed inside a group. σ_g,s are photon-weighted over the formal
solution's own spectrum (`IONFSNSUB` = 8 frequencies per group), E_g the mean photon energy. What is lost against IonFS's 11 groups is the
hardening inside a band: the last photons of a band get the band-mean σ. A front zone has τ_912 ~ 1e3 and absorbs a whole band anyway, so
the ionizations per zone (photon conservation) are unchanged; only the share of the last e-folds moves by a zone. Each group costs one
row per zone and widens the coupled zones' local Jacobian by one; 4 keeps neqn at 13 (the ForwardDiff Jacobian carries 6 × 13 = 78
partials, against 54 with COUPLED alone). H₂ photoionization (15.4 eV) stays out, as in IonFS (H₂ ≤ 1e-15 where the ray reaches).

## Lyman–Werner: photon-conserving and consistent with SLSHIELD

STARLIGHT's rate per molecule is k = 1.1e8 F_ν(12.87 eV) f_sh(N_H₂, T) f_HI(N_HI) (Abel et al. 1997, k27, 1997NewA....2..181A;
Wolcott-Green, Haiman & Bryan 2011 eq. 10 with α = 1.1, 2011MNRAS.418..838W; Wolcott-Green & Haiman 2011 eq. 13, 2011MNRAS.412.2603W;
`fsh_h2`, `fsh_hi` in PopIIIMaterial.jl). Draine & Bertoldi (1996, 1996ApJ...468..269D) relate the shielding function to the band's
absorption (their eq. 38): the UV pumping (absorptions) per area down to column N is ζ_pump(0)/σ-normalised ∫₀ᴺ f_shield dN′, and
"⟨p_diss⟩ ≈ 0.15 averaged over all UV absorptions" (their §4; Table/Fig. values 0.10–0.21) of them dissociate. So:

- σ_d = 1.1e8 L_ν(12.87 eV)/N_LW from the formal solution (dissociations per molecule per LW photon flux; 4.02e-18 cm² at the
  t_c128c B7 source, so that σ_d N_LW/(4πr²) = 1.1e8 F_ν(12.87 eV) exactly);
- G(N, T) = ∫₀ᴺ f_sh dN′ in closed form (`ifsi_G`; checked against the derivative and inverted to 1e-13, `scratchpad/ifs/gtest.jl`);
- the LW unknown carries the band photons not yet absorbed by H₂; the photons used so far give the zone's equivalent static column N₁
  at its inner face, N_LW(k) = Φ[1 − (σ_d/p) G(N₁)], p = ⟨p_diss⟩ = 0.15 (`IONFSPD`);
- zone k absorbs A = Φ(σ_d/p)[G(N₁ + ΔN) − G(N₁)], ΔN = x_H₂ n_H Δr (with SLSHIELD=sobolev the zone's mean f is at least f_sh(N_sob),
  STARLIGHT's min(N_static, N_sob)), written as Δτ_LW = A/N_LW(k): the row is the same as the ionizing groups' and a zone never takes
  more than enters it; dissociations D = p N_LW(1 − e^{−Δτ_LW}) f_HI(N_HI) (N_HI the start-of-step column, as STARLIGHT).

For a static column this gives STARLIGHT's σ_d φ f_sh(N) f_HI exactly; p only matters where the carrier runs out. With σ_d = 3.75e-18 cm²
(flat F_ν), (σ_d/p) G(1e21 cm⁻², 1000 K) = 0.98: the band is used up near the upper validity limit of the fit, which is DB96's statement that
the saturated H₂ lines take an equivalent width Δln ν ≈ 0.2 ≈ the band's ln(13.6/11.2) = 0.19. H⁻ in the LW band (σ ~ 1e-17 at 12 eV,
x_H⁻ ≤ 1e-8) is left to the diluted blackbody below the edges. The LW heat (4.48 + 0.4 eV per dissociation, GASQ) and HD's LW rate use the
ray's rate (`ifsi_klw` replaces STARLIGHT's kLW in every zone k ≥ k₀; the SLSHIELD value is kept for the comparison).

## Rates above the edges come only from the ray

In every LOCALJ zone k ≥ k₀ (outside the grey photosphere; `IONFSLJ=all`: every LOCALJ zone) the diluted blackbody W B(T_c) keeps only the
levels below 13.598 eV (`ljcut`), with case-B He recombination (`caseB`): the ray's zones as in IonFS, and now also the zones the ray does
not reach. Inside k₀ (τ_grey > 2/3, T_c → T_r, W → 1) the field is near-thermal and LOCALJ is kept (its purpose: He in thermalised network
zones). Hydrogen 1s recombination is on the spot through the three-level atom (C α_B), as before.

## He recombination and He⁺ line cooling (missing terms, also in IonFS)

ChemNet keeps binding energies in u: a photoionization takes χ from the gas, a recombination gives χ back. IonFS gives the gas the whole hν
of an absorbed photon, so a photoionization–recombination cycle heated the gas by hν per absorbed photon instead of hν − χ − ⟨E_kin⟩: there
is no He recombination radiation in the network (H has it in the three-level chemq, −C α_B(χ_H + 0.75 kT)). At the He I front a single
zone went to 9e4 K (t_ifs9, 29 kK and 92 kK in zones 86–87). `ifsi_herec` adds, in the ray zones, the network's own case-B radiative
recombinations (the Milne inverses α of its He I and He II photo-processes in the ray's env) × (χ + 0.75 kT) to the cooling (implicit, in
the energy row); 0.75 kT is the three-level atom's (tl3_heat).

With that, zones heated by the ray's H I-band photons still went to 5–7e4 K (t_ifs10: photoionization equilibrium of H with only its
recombination cooling, ⟨hν⟩ − χ_H ≈ 3.4 eV = 0.75 kT at 5e4 K) and produced transient He²⁺ out to 1e5 R☉. The network has no He line
cooling (ChemRates: "He⁺ two-photon and line cooling — these come with the model atom"). `ifsi_herec` therefore also adds He⁺ collisional
excitation (n = 1 → 2), ceHeII = 5.54e-17 e^{−473638/T} T^−0.397/(1 + (T/1e5)^½) erg cm³ s⁻¹ × n_e n_He⁺ (Cen 1992, 1992ApJS...78..341C,
as coded in EmissionKernels.jl cooling_atomic.jl for the AAZN97 network; `IONFSHECE=0` drops it). After it the hottest ray zone is the
first one (2.3e4 K) except for the first step and one of 42 snapshots (4.9e4 K at the He I front; t_ifs15). H's collisional excitation
is the three-level atom's (LYABB). The energy of
both terms goes to the grey field (RAYQ). This probably contributed to t_ion1's hot, Δt-limiting zones as well.

## Tests

### (a) Bit-identical off
200 steps from t_c128c B7_auto1 (production run.env), popiii-s6 79d4dc1 (`t_ifs1`) against this branch with the switches off
(`t_ifs2` at 24d9e23, `t_ifs13` at 8a79b7b): identical uncompressed final state, identical .aux, all snapshots identical in every field
(`scratchpad/ifs/bitsnap.jl`), identical non-timing perf.csv columns; t_ifs1 also reproduces production B8 over those steps. Later commits change only IONFSIMP paths, diagnostics (NEWTONERR, IONFSHEDBG) and comments.

### (b) Jacobian against finite differences
`CPLJCHECK=74,76,82` with `IONFSIMP=1` at the first Newton iteration from t_c128c B7 (t_ifs8, dd20512); `TitanResidFast.cpl_jcheck` now
includes the ray's rows and columns; `ifsi_gcheck` compares a ray zone's 8-wide local Jacobian with central differences.
- Ray rows (10–13) against every column (T, x_H⁺, x_H₂, y₁…y₄): AD = FD to ≤ 6e-5 (the ray rows themselves: 2.5e-9).
- Zone's p, e, S, Λ against y₁…y₄ and ρ: agree to 1e-3 (worst: Λ against y₃ 0.6 %, Δτ₂ against y₃ 5 % in the He II-front zone 74).
- d/dx_H⁺ in zones 74/76 differ by ×2 only because x_H = 2e-8 < the FD step (1e-6 x_H⁺), so x_H clamps at 0 in one of the two
  evaluations; zone 82 (x_H 2e-5) agrees to 1e-3. d/dx_H₂ (x_H₂ ~ 1e-26) is below FD resolution.
- d/dT of the He-dependent outputs (e, Λ, Δτ₂ of the He-ionized zones) differs: the COUPLED design takes the network's photo-rates and
  their Milne inverses at the Float T (`ck_eval`: `phc(ForwardDiff.value(T))`), so ∂α_rad/∂T is not in any coupled zone's Jacobian,
  with or without the ray. The ray's own T dependence (LW b(T), the He fixed point's ∂F/∂T) is in.

### (c) Runs

Settings of the production run directories plus `IONFS=1 IONFSIMP=1`. The final-code runs were stopped by the PI's request (no runs past
500 M☉) before 2000 steps; t_ifs15 (542c651: differs from the final code only in the He fixed point's iteration details and logging) ran
2100 steps. Production numbers are over the same steps from the same checkpoint (c128c B8; nolte2 B3, which restarted from B2_auto0).

| run (from) | steps | median Δt | mean Δt | divergences/step | chem fallbacks | wall/step |
|---|---|---|---|---|---|---|
| t_ifs15, 542c651 (c128c B7) | 2100 | 1.9e3 s | 2.3e5 s | 1.11 | 10 | 3.8 s |
| t_ifs25, 72fa45b (c128c B7) | 879 | 1.4e3 s | 1.3e5 s | 1.25 | 0 | 3.6 s |
| production c128c B8, same span | 2000 | 1.1e7 s (100-step marks) | 3.4e7 s | 0.85 | 1 | 1.3 s |
| t_ifs26, 72fa45b (nolte2 B2_auto0) | 455 | 3.7e2 s | 1.5e5 s | 1.22 | 1 | 7.5 s |
| production nolte2 B3, same span | 400 | 6.0e7 s | 2.4e7 s | 0.98 | 1 | 0.6 s (8 threads) |

(The first step from each checkpoint is accepted at the full Δt, 8.7e7 / 9e6 s; the collapse comes after it.) Δt is 3–4 orders of
magnitude below production; see Open issues 1. t_ion1 (explicit IonFS) had 0.02–300 s.

**He²⁺ inside 3e3 R☉** (snapshots, network zones outside the grey photosphere):
- c128c: production has He²⁺/He 0.02–0.83 at zone 78 and > 1e-3 out to zones 82–84 (2.3–3.8e3 R☉) in all 41 snapshots. With IONFSIMP
  He²⁺ ≤ 1.3e-10 there in all 42 (t_ifs15) and 17 (t_ifs25) snapshots; He²⁺ > 1e-3 only in zone 73 (157 R☉, the first ray zone, inside
  the snapshot's photosphere: the He II photons of the formal solution, Q_HeII 3e45–6e47 s⁻¹) and transiently (4 of 17, 5 of 42 snapshots)
  in a hot He I front zone at 2–3e4 R☉.
- nolte2: production has He²⁺ 0.999 at zone 76 out to zone 81 (1867 R☉) in all 41 snapshots. t_ifs26 still has He²⁺ 0.72–1 in zones
  76–77 (460–660 R☉): there the ray's own He II photons ionize the first network zones; beyond them He²⁺ is gone except in the transient
  front zones (out to 1.2–4e4 R☉ in 4 of 9 snapshots). So LOCALJ's spurious He²⁺ is gone in both; the He III region of the real He II
  photons remains, and in nolte2 it reaches the first two network zones.

**Photon bookkeeping** (emitted = absorbed + escaped through the last zone, per group, from the Newton's ray): |Q − A − esc|/Q median
4e-11 (H I), 2e-5…1e-4 (He I), 0 (He II, all absorbed), 4e-10 (LW); worst steps 0.2 (He I, c128c), 0.09 (nolte2): the ray rows are
converged to the hydro Newton's tolerance, and the worst steps are those where zones were decoupled after the Newton. All H I, He I and
He II photons are absorbed on the mesh; 4.5e-6 of the LW photons are.

**H₂ photodissociation** (Σ over zones k ≥ k₀, at the same composition): ray 8.9e44 s⁻¹ against 2.0e44 s⁻¹ from STARLIGHT's SLSHIELD
rate in c128c (ratio 4.45, 4.3–5.7), 6.3e44 against 1.6e44 in nolte2 (ratio 3.9). The shielding factors agree by construction; the
difference is the field: ×2 L_run/L_fs and the formal solution's harder LW spectrum (see below). The gas outside the photosphere is atomic
(x_H₂ ≤ 1.2e-11), so the LW band is unshielded and the rates matter only for the trace H₂ and HD.

Earlier runs (in the scratchpad notes): t_ifs9/10 without the He recombination / He⁺ line cooling (hot front zones 5e4–9e4 K, transient
He²⁺ out to 1e5 R☉), t_ifs21/24 with the loose He tolerance (median Δt 5–14 s).

## L_run/L_fs ≈ 2

Over the c128c run L_run/L_fs = 1.99 (1.9–2.08): the formal solution of the table region (continuum opacity, static, thermal inner
boundary at grey τ = 300) carries half the grey luminosity through r_b. `IONFSNORM=1` scales every group, LW included, by it. What it does
and does not change here:
- photon budgets ∝ f: the He I front (He⁺ > 0.5) sits near its Strömgren radius, at 2–4e4 R☉, where the recombinations inside balance
  Q_HeI. For the infall's n ∝ r^−1.5…−2 the recombination integral grows only logarithmically or not at all with radius, so the front
  position is very sensitive to Q: a factor 2 in Q can move it by many zones (unlike the uniform-density r ∝ Q^{1/3}).
- the LW rates ∝ f in every unshielded zone (all of them here: x_H₂ ≤ 1.2e-11 in every zone outside the photosphere out to 1e9 R☉, and
  99.9995 % of the LW photons leave the mesh).
  The ray's dissociations exceed the SLSHIELD path's by ×4.6 (median; 4.2–6.2); ×2 of it is f, the rest the spectrum: the formal solution's
  L_ν(12.87 eV) at r_b is 2.6e25 erg s⁻¹ Hz⁻¹ (before f), a T_eff = 2.3e4 K blackbody of the same L_fs at r_b gives 4.3e24, so the
  escaping LW comes from layers hotter than T_eff. STARLIGHT's W B_ν(T_c) assumes the blackbody.
- not the temperatures of the photo-heated zones (set by ⟨hν⟩ − χ against the cooling, i.e. by the spectral shape), not the absence of
  LOCALJ's He²⁺, not the Newton coupling.
So it matters at the factor-2 level for where the He I front sits and for the LW rates, and should be resolved (which T structure is
right at r_b: the grey run's or the formal solution's) before front positions are quoted. Not changed here.

## Open issues

1. **Δt is not production-like.** The coupling is implicit and converges, but the hydro Newton fails (≈1.1–1.3 divergences per step
   against 0.85–0.9 in production) at the zone that holds the He I front, and Δt sits at 1e2–1e4 s (production 1e6–1e7 s). The limiting
   zone is the He I front zone (dtlim T@89–92 in c128c, and the largest T and x corrections of the diverging Newton solves are there), not
   the photosphere zones 73–79 of t_ion1. (A diagnostic with STOL = 1, t_ifs23, ran with the loose He tolerance below and is
   inconclusive.) The front sits near its Strömgren radius in a density profile for
   which that radius is very sensitive to Q and to the inner zones' recombinations, and it moves by several zones within a step (t_ifs15:
   He⁺ > 0.5 out to zone 92 → 82 → 89 within 100 steps around step 2256000, after one step of 1.4e6 s). A zone-mean He front has no
   sub-zone model; IonFS's two-phase recombination term κ₂ exists for H only. Candidates: the same two-phase term for He I/He II, or
   carrying x_He⁺ as a Newton unknown (like x_H⁺) instead of the CK-algebraic He chain with the inner fixed point.
2. The He fixed point stalls above 1e-11 in ~0.5 % of the solves (round-off where x_H ~ 1e-8 or He⁺ ≪ He²⁺); the best iterate is then
   history-dependent at that level. A looser tolerance was worse (Tests). Better: a He-chain formulation without the cancellation.
3. Cost: ~3.5–5 s per step (production 1.3 s): 13 instead of 9 unknowns, the 8-wide local Jacobians of ~45 ray zones and the He Newton.
4. The k_end criterion makes the ray span the whole mesh for LW (k_end → 128); the LW band is ~unabsorbed (x_H₂ ≤ 1e-11), so it costs
   rows without changing anything. An LW-only cut (switch the LW group off where Δτ_LW < 1e-6 for all zones) would save most of the ray
   zones beyond the He I front.
5. L_run/L_fs ≈ 2 (above). NLOIT's opacity iterate (`nlo_comp`) does not include the implicit ray's rates (NLOIT is off in production).
6. H₂ photoionization and H⁻ in the LW band are left out (trace absorbers here).

# ionfs-v2 (16 Nov): onto popiii-s6 1782928, the z300mm moving-mesh protostar

Branch `ionfs-v2` = popiii-s6 1782928 + ionfs-implicit 4698bbe + the changes below; all behind `IONFS=1 IONFSIMP=1` (default off).

## Merge
- The coupled-zone rates now run in the threaded zone sweep (PERF512); ray zones call `CPL_RAY` there (per-zone caches, atomics only).
- SLIMP=2: ray zones take the star light at the iterate (`sl_qt` → `cpl_tq`, Float, lagged in the local Jacobian); `cpl_tq`'s rebuilt LW
  rate is the ray's (0 in a ray zone, the ray's lagged rate elsewhere outside r_ph; `ifsi_lwrec`), so SLIMP no longer re-adds STARLIGHT's
  W B_ν(T_c) LW on top of the ray.
- NLOIT: a ray zone's opacity composition is the ray's CK expansion (`nlo_comp`).
- GRIDEQ: the ray's geometry (Δr, V) is the Newton iterate's (Float, lagged within an iteration).
- **Bit-identical off:** 40 steps from z300mm's z19t_mm24S3_auto0 with z300mm's settings, 1782928 (`t_ifm1`) against 0640880 (`t_ifm2`)
  and cad6ca4 (`t_ifm5`): identical uncompressed final state (md5 28ead850…) and .aux (ae6d4c19…).

## L_run/L_fs: the cause
IONFSLDIAG solves the same formal solution once more with the run's own grey opacity (χ = chifn ρ, κ_P ρ, B(T)). At z300mm (r_b 25 R☉):
- before: L_fs = 2.4 L_run, and the grey formal solution 2.19 L_run. The photosphere is a transparency front inside one zone
  (τ_grey 732 → < 2/3 across zone 77, SUBZONE); the formal solution had that zone's centre T (7302 K) out to r_b, against
  T_eff = 6000 K from L_run: (7302/6000)⁴ = 2.19. The run itself carries the sub-zone front (SUBZONE), the formal solution did not.
- fix (`IONFSSUB`, default on with IONFSIMP): the outer part of the front zone as sub-points on the zone's grey τ scale (τ = 300 … 0.001,
  0 at r_b) with the grey Eddington T⁴ = ¾ T_eff⁴(τ + 2/3), capped at the zone's T, and the source opacity there scaled to the zone's grey
  opacity per length. The grey check then gives 0.94–0.97 L_run: the method is consistent with the run.
- what remains, L_fs ≈ 2.8–3.6 L_run with the frequency-dependent continuum opacity: the source opacity is continuum only (no lines).
  At zone 75 (12 kK) its Rosseland mean is 0.6 of the run's χ and its Planck mean 1/60 of the run's κ_P; the windows let flux out from the
  12–15 kK layers below. That is missing line blanketing in the optical/IR, not in the Lyman continuum (H bound–free dominates there).
  IONFSNORM=1 still scales all groups by L_run/L_fs (0.27–0.6 at z300mm), so the UV numbers below are conservative by that factor;
  `fnorm` is in every CSV row (unscaled = Q/fnorm).
- The c128c ratio (0.5, the other sign) was not re-examined: that run is above 500 M☉.
- The source now starts in the innermost zone with grey τ > 3e3 (`IONFSTAUC`), at least 6 zones inside r_b (`IONFSNZMIN`), so
  everything that can leak through the front (the precursor and shock layers below τ ~ 3e3) is in the formal solution.

## Diagnostics: <RUNREC>/ionfs_snap.csv (every snapshot)
See the docstring of `ifsi_snaplog` (IonFSImp.jl): emitted (at r_ph) and escaping (outer boundary) Q_LW, Q_H, Q_HeI, Q_HeII; the H, He I,
He II ionization fronts and the H₂ fronts (f_H₂ = 0.1, 0.5) relative to r_ph; dissociations per s of the ray and of the SLSHIELD rate; the
ray's ionizing and LW heating powers; and at 1e2, 1e3, 1e4, 1e5 R☉ the photon rates crossing, T, u, x_e, f_H₂, both LW rates, the ray's
heating per mass and t_flow/t_heat. `ionfsi.csv` (every step) keeps the photon bookkeeping per group.

## z300mm: does the light get ahead of the photosphere? (t_ifm7, cad6ca4)

From z300mm's z19z300mmS1_auto0 (step 5000, M★ ≈ 0.135 M☉, r_ph ≈ 27–29 R☉, T_eff 5500–6450 K), its full settings + `IONFS=1
IONFSIMP=1`, 8 threads on cores 96–103, 2040 steps (5001–7041, 0.076 yr). Snapshot rows in `t_ifm7/ionfs_snap.csv`, every step in
`t_ifm7/ionfsi.csv`; comparison with z300mm at matched times (`scratchpad/ifs/infallcmp.jl`).

- **Ionizing photons: no.** Q_H leaving r_ph 3e40–4.6e41 s⁻¹ (after fnorm 0.2–0.4), Q_HeI ~1e33, Q_HeII ≤ 1e3. All are absorbed
  in the first zone outside r_ph: Q_H = 0 at 100 R☉ in every snapshot, rIF_H = rIF_HeI = rIF_HeII = r_ph, photo-heating power ≤ 1e31 erg/s
  (L★ ~ 3e36 erg/s).
- **LW photons: they get ahead, but cannot dissociate.** Q_LW 1.5e41–3.8e42 s⁻¹ leaves r_ph, 56–68 % of it leaves the mesh (1e5 R☉ and
  beyond). The dissociations are confined to the first ~0.1–1 R☉: the H₂ front (f_H₂ ≥ 0.1/0.5) sits 0.11–0.37 / 0.37–1.07 R☉ outside
  r_ph. Beyond, the equivalent H₂ column of the photons already used exceeds the WGHB fit's range and the H I Lyman-line shielding of
  the atomic precursor (N_HI ≫ 1e24) is complete: k_LW = 0 at 100–1e5 R☉, in the ray and in STARLIGHT's SLSHIELD path alike. The ray's
  dissociations (1e14–2e19 s⁻¹) exceed SLSHIELD's (1e11–1e14 s⁻¹) by 10³–10⁶, all in those first zones. That gas is warm and H₂-poor
  anyway, so the front moves by < 0.5 R☉.
- **The infall: unchanged.** At 1e2, 1e3, 1e4, 1e5 R☉, T agrees with z300mm to ≤ 1 % (2120/2125 K at 100 R☉), u_r to ≤ 3 %
  (−21 km/s at 100 R☉), f_H₂ to 3 digits (0.87, 0.97, 0.70, 0.011). x_e at 100–1e3 R☉ is 1.3–2× z300mm's (2e-8 against 1e-8). r_ph and
  the H₂ front differ by ≤ 0.5 R☉, within the matched-time spacing of z300mm's snapshots. Ray heating at those radii is 0 (t_flow/t_heat 0).
- **Photon bookkeeping:** emitted = absorbed + escaped to ≤ 8e-9 (LW). All ionizing groups are absorbed in the ray's first zone.
- **dt and failures:** mean Δt 1.18e3 s (median 1.1e3 s) against z300mm's ≈ 1.1e4 s over the same evolution (steps 5018–5298:
  0.1 yr in 280 steps). The run had 114 Newton failures in 2040 steps (0.056/step; z300mm 0.13/step) and no error rejections. Δt is set by
  DTERR's error estimate on x_H₂ and x_H⁺ in zones 75–76: the ray zones next to the photosphere, where the ray's LW and Lyman-continuum
  rates make the chemistry stiffer. Wall time ≈ 4.4 s/step (z300mm ≈ 3.5 s/step).

So for this protostar, the precursor sees neither ionizing radiation nor LW dissociation beyond ~1 R☉ of the photosphere. The LW photons
that escape do so through the gaps between the saturated H₂ and H I lines. The SLSHIELD treatment already reproduces this. The infall's
thermodynamics does not change.

Open:
- the continuum-only source (L_fs ≈ 3 L_run) is scaled by fnorm. If the UV excess is real, Q_LW would be ×3–5 higher; that does not
  matter while the H₂/H I shielding is complete.
- the WGHB fits are extrapolated far beyond N ~ 1e21 / 1e24.
- the He I front treatment (two-phase term or x_He⁺ unknown) is untested here: at T_eff ~ 6000 K there is no He front (Q_HeI ~ 1e33 s⁻¹,
  absorbed in one zone).
- the Δt cost above.
