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
  |F| < 1e-6 in ln u: x_H = 1 − x_H⁺ − … carries round-off ~1e-16/x_H, and where He²⁺ dominates the CK expansion's fixed sweeps leave
  ~1e-7 in He⁺; Newton stalls there with a well-conditioned J. Zones whose He photons are negligible (N_He Δt < 1e-12 of their He) take
  the start-of-step fractions. A plain fixed-point iteration was tried first: in a photon-limited He zone its slope is (x₀ − x_He)/x₀ → 1,
  it did not converge, and the residual then depended on the previous call (CPLJCHECK ratios 0.5–1 from that noise).
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
(`scratchpad/ifs/bitsnap.jl`), identical non-timing perf.csv columns. RESULTS_A

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
RESULTS_C

## L_run/L_fs ≈ 2
RESULTS_L

## Open issues
RESULTS_O
