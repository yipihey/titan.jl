# G2: one non-LTE material everywhere — audit of table vs network material (16 Nov 2026)

Branch `g2-material` (from `popiii-s6` 0781136). Question (Tom): can the network material (EOS and opacities from the
network's own populations) reproduce the equilibrium-table material where the table is used, so that the per-zone
table/network switch can be dropped? Motivation: in zdiagB (IC abn02red, M* ≈ 0.057 M☉, R* ≈ 10 R☉) the photosphere cliff
that stalls Δt sits at the last table zone / first network zone.

Everything below is measured [run] unless marked. Scripts: `scripts/g2audit/` (standalone; no TITAN needed except the bit test).
Plots: `docs/g2/g2_gap_maps_g2opac.png`, `docs/g2/g2_gap_maps_hm88.png`.

## 0. Summary

1. **At the cliff the two materials agree in LTE where the network is used.** For the network-side zones (T ≤ 6000 K,
   ρ ≤ 1e-6, zdiagB zones 62–160) network-in-LTE vs table: p within 1 %, c_P within 4 %, κ_R within 6 %, κ_P within 5 %.
   The cliff itself is physical: across zones 61 → 62 T falls 9600 → 3350 K and κ_R falls 44 → 5e-5 cm² g⁻¹ (×10⁶);
   in LTE neither material removes it. What the switch adds is (a) the zones' **non-LTE state** (zone 62: x_e 1.2e-8 vs
   Saha 3.8e-9, κ_R ×1.7, κ_P ×2.4 vs LTE) and (b) a **κ_P discontinuity**: a network zone at the table side's
   state (zone 61, 9600 K) has κ_P ≈ 8× below the table zone's (NLOPAC/NLOIT drop H bound–free from κ_P, and the table's
   κ_P above 6300 K is an LTE line Planck mean).
2. **Where the table is used (T > 6300 K, and ρ ≳ 1e-4), the network's own physics falls short by large factors**, for
   three separable reasons: (i) no atomic lines / Stark broadening (OPAL κ_R, Z0 line κ_P); (ii) an ideal EOS without
   occupation probabilities or pressure ionisation (p 0.66–0.73 × table in the interior, c_P 0.5–0.7×); (iii) bookkeeping
   bugs of the network's continuum code (no Lyman continuum, means grid stops at 12.4 eV, all neutral H treated as ground
   state). Along the zdiagB profile (161 zones, 1e3 < T < 3e5 K) κ_R is off by > ×3 in 24 % of the zones, κ_P in 23 %,
   p by > 30 % in 30 %.
3. **Table bug found: the SCVH95 composition columns are misread** (×1.85–1.94 in x_H⁺ where H is weakly ionised). Fixed
   behind `SCVHXE=1`. It feeds `x_from_table` (the composition a zone gets when it leaves the table), the electrons of the
   scattering term, the equilibrium cooling and the NLOIT calibration base.
4. Prototype `G2OPAC=1` (`julia/G2Opacity.jl`): the network continuum with the table side's continuum ingredients. With
   the network's populations it fixes the cold/low-density side, but above 3e4 K the network's own Saha (partition function to
   n = 15, no occupation probabilities) overpopulates neutral H by Z_H/2 ≈ 10–500 and κ stays high. With HM88
   occupation-probability populations (as G2 would have them) the zones beyond ×3 drop from 24 % to 17 % (κ_R) and 23 %
   to 7 % (κ_P); what remains is lines (6300 K–5e4 K, ρ > 1e-8) and the dense corner. Off = bit-identical (8-step restart of
   zdiagB S5: 591 state arrays, 10 530 280 bytes, and the 478 907-byte .aux identical to 0781136).
5. Plan (§4): (P1) HM88 occupation-probability partition functions + excitation energy in ChemNet (EOS and Keq);
   (P2) the G2 continuum, tabulated per absorber to cost ~0.3 ms; (P3) LTE line opacity as a tabulated per-absorber
   line-to-continuum factor (cannot be made non-LTE cheaply); (P4) a non-ideal excess free energy ΔF(ρ, T) from SCVH added to the
   network EOS, which vanishes below ρ ≈ 1e-5. Below a few 1e-6 g cm⁻³ (the photosphere) P1–P2 reproduce the table to ≤ 10 %
   in LTE. A smooth weight at the Newton iterate (§4.6) is the fallback for lines and the dense corner.

## 1. Where each material comes from

**Table material** (zones with `eqzone[k]`; PopIIIMaterial.jl:1344–1366, rules PSWITCH (p ≥ 300), PSWTMIN 4500 K, PSWHYS 0.3,
PSWTEX 5000 K, PSWTAU 3 (τ ≥ 3 and T ≥ 5000 K → table), PSWTAUDEEP 3, PSWNH (off)):
- EOS `PopIIIEqEOS` (cache `eqeos_table_lr20_he.jls` with HECHEM=1): Saha H₂ ⇌ 2H, H ⇌ H⁺ + e, He ⇌ He⁺ ⇌ He²⁺ on one n_e,
  ideal gas, RHFC02 A9 c_v of H₂ (PrimordialChem.eos_rhfc02) for p ≤ 1e4; SCVH95 (`SCVH.jl`, YREC copies of h/he_tab_i.dat)
  for p ≥ 1e6; smoothstep in log p in between; bicubic in (log ρ, log T). Above T1 = 3e5 K the HOT branch (StellarEOS:
  Z0 Saha / FD electrons + Yakovlev–Shalybkov Coulomb, free-energy blend) — not part of this comparison.
- Opacity `OPAC[] = StellarOpacity` (table `kaptable_stellar.bin.gz`, built by `stellaropacity_build.jl`):
  κ_R: own continuum means (`StellarOpacityLow` = ContinuumOpacity physics in Saha equilibrium incl. Li electrons and He
  Rayleigh) below log T 3.8; blend to OPAL GN93hz Z = 0 at 3.8–4.1; OPAL above; Cassisi conduction above 1e4 K.
  κ_P: own continuum below 3.8; Z0Physics Planck mean **with H I/He I/He II lines** (`kappa_means(lines = true, lines_R = false)`)
  above. κ_E = κ_P.
**Network material** (`NetGas`, ChemNet with netmodernli + SL3 three-level H, HECHEM=1):
- EOS `ChemNet.eos`: u = Σ x_i (c_v,i kT + E0_i)/μ_H, p = Σ n_i kT (ideal; frozen composition inside the hydro Newton,
  COUPLED zones with x_H⁺, x_H₂ as unknowns).
- Opacity: NLOPAC/NLOIT (`nlte_means`, `nlo_ratio`): the continuum code of ContinuumOpacity on the zone's populations
  (H⁻ from the network, n = 2 from the three-level atom); κ_E without H bound–free, κ_P with H⁻ in Saha with the zone's H and
  e, κ_R with H bf. Below 5000 K the network's own means; above 6000 K table × (network means / network means at equilibrium)
  (`NLOITW`), so in LTE it equals the table there by construction (except κ_E/κ_P, see §3.3).

### Ingredient table (✓ has it, ~ partial, – no)

| ingredient | table | network | notes |
|---|---|---|---|
| ideal gas, H₂ c_v (RHFC02 A9, 6100 K) | ✓ | ✓ | same formula |
| H₂ ⇌ 2H | ✓ (Saha; SCVH above 1e4 dyn) | ✓ (ChemNet K from Z_rot Z_vib, D0 4.478 eV) | 1–7 % in x_H₂ (SCVH weight) |
| H ionisation | ✓ Saha Z_H = 2 (ground only); SCVH | ✓ Z_H = Σ_{n≤15} 2n² e^{−E_n/kT} | network: no occupation probabilities, **no excitation energy in u** (u and Keq inconsistent); overpopulates H⁰ ×Z_H/2 above 3e4 K |
| He, He⁺, He²⁺ | ✓ (Saha, same levels as ChemNet) | ✓ (HECHEM) | identical below SCVH weight |
| H⁻, H₂⁺, Li/Li⁺ | – in EOS (opacity LOW has H⁻ and Li e⁻) | ✓ | Li e⁻ set x_e below 3000 K: 4e-10 vs table EOS 3e-12 (irrelevant for p) |
| H₃⁺, HeH⁺ | – | – | |
| excited-level partition functions | ground only (Saha); occupation-prob. (SCVH) | n ≤ 15 hard cut | |
| pressure ionisation, Coulomb, excluded volume (non-ideal) | ✓ SCVH (p ≥ 1e4…1e6); HOT: Coulomb | – | dominant EOS gap ρ ≳ 1e-4 |
| Debye–Hückel | ~ (in SCVH) | – | |
| H⁻ bf / ff (John 1988) | ✓ | ✓ | same code |
| H bf | ✓ Z0/OPAL: Gaunt, HM88, dissolved levels; LOW: Kramers | Kramers g = 1, n ≤ 30, all H⁰ as ground state | both LOW and network **lack the Lyman continuum** (grid ends at 1e5 cm⁻¹) |
| H⁺ ff | ✓ van Hoof g_ff (Z0), OPAL | Kramers g = 1 | |
| He I bf, He II bf, He⁺/He²⁺ ff | ✓ (Z0, OPAL) | – | |
| H₂–H₂, H₂–He, H₂–H CIA | ✓ | ✓ | same tables |
| H₂ lines, H₂⁺, H₂⁻ ff, He⁻ ff, molecular bands | – | – | H₂ lines are H2Lines' cooling/transfer, not the continuum |
| Rayleigh H, H₂ | ✓ (capped at 1300 Å in Z0/LOW) | ✓ (uncapped) | He Rayleigh: LOW only |
| Thomson | ✓ (+ degeneracy factor) | ✓ | |
| H/He lines, Stark/pressure broadening, line blanketing | ✓ κ_R: OPAL above 6300 K (blend); κ_P: Z0 lines (Holtsmark + impact, HM88) | – | the line κ_P is an LTE Planck mean |
| conduction | ✓ (Cassisi) | – | irrelevant below 1e5 K here |
| means supplied | κ_R (extinction), κ_P = κ_E (absorption) | κ_E (absorption, B_ν(T_c) weight with NLOPACJ), κ_P (emission), κ_R | |

## 2. The gap in LTE (network populations set to Saha with ChemNet's own species data)

Grid log ρ −14 … 0 (Δ 0.1), log T 3 … 6 (Δ 0.02), 21 291 points; network composition by `net_lte` (all network species,
charge balance); thermodynamic derivatives of the equilibrium composition by centred differences (h = 1e-4 in ln). Table
x_e with the SCVH weight > 0 shown both as read and with SCVHXE. Maps: `docs/g2/g2_gap_maps_g2opac.png` (rows: EOS;
x_e and κ of the network as is; κ with G2OPAC on the network's own populations; x_e diagnostics) and
`docs/g2/g2_gap_maps_hm88.png` (row 4: G2OPAC on HM88-Saha populations). Contours |Δ| = 5 %, 30 %, ×3; zdiagB step 2650 profile overplotted.

### 2.1 Numbers at the cliff (zdiagB step 2650, `scripts/g2audit/cliff.jl`)

| zone | log ρ | log T | log τ | x_e table / SCVHXE / net LTE | c_P table / net | κ_R table / net / G2 | κ_P table / net / G2 | zone's own (non-LTE): x_e, κ_R, κ_P |
|---|---|---|---|---|---|---|---|---|
| 58 (table) | −5.33 | 4.26 | 9.06 | 0.654 / 0.486 / 0.459 | 1.60e9 / 1.40e9 | 5779 / 2337 / 2116 | 4.97e5 / 8782 / 34810 | (table) |
| 61 (table) | −6.42 | 3.98 | 5.76 | 0.067 / 0.038 / 0.031 | 5.21e8 / 5.15e8 | 44.1 / 35.8 / 34.9 | 6226 / 60 / 207 | (table) |
| 62 (network) | −6.66 | 3.53 | 0.50 | 3.9e-9 / 2.7e-9 / 3.8e-9 | 1.75e9 / 1.56e9 | 4.74e-5 / 4.96e-5 / 4.95e-5 | 3.26e-5 / 3.35e-5 / 3.34e-5 | 1.23e-8, 8.3e-5, 8.0e-5 |
| 66 (network) | −7.05 | 3.42 | 0.23 | 3.4e-12 / 3.0e-12 / 4.1e-10 (Li) | 1.54e9 / 1.54e9 | 1.86e-5 / 1.97e-5 / 1.97e-5 | 1.66e-5 / 1.74e-5 / 1.74e-5 | 4.2e-10, 2.0e-5, 1.7e-5 |

κ_P at zone 61's state with the network's convention (no H bound–free, hbfP = false): 26 (G2) — via NLOIT's ratio that is
table × 26/207 ≈ 790, i.e. 8× below the table zone next to it.

### 2.2 Along the profile (161 zones, 1e3 < T < 3e5 K; fraction of zones beyond 5 % / 30 % / ×3)

| quantity | network as is | G2OPAC, network populations | G2OPAC, HM88 populations |
|---|---|---|---|
| p | 0.35 / 0.30 / 0 | | |
| c_P | 0.37 / 0.36 / 0 | | |
| Γ1 | 0.13 / 0 / 0 | | |
| κ_R | 0.43 / 0.37 / 0.24 | 0.42 / 0.30 / 0.06 | 0.40 / 0.24 / 0.17 |
| κ_P | 0.43 / 0.37 / 0.23 | 0.40 / 0.35 / 0.18 | 0.25 / 0.12 / 0.07 |

(G2OPAC on the network's own populations scores better in κ_R on the profile only because two errors cancel there in the interior:
the neutral overpopulation raises bf opacity where OPAL's lines are missing. The HM88 column is the meaningful one.)

### 2.3 Where and why (attribution)

- **T < 6000 K, any ρ ≤ 1e-6** (the photosphere and everything outside): agree within 1–6 % (EOS and κ). Residuals: c_P
  4–11 % where the SCVH weight is 0.25–0.8 (SCVH's H₂ dissociation, X_H₂ 0.131 vs 0.141 at zone 62); κ ≤ 6 % (He
  Rayleigh, Li electrons in the table's means, CIA/H⁻ grid). Network κ_P lacks the Lyman continuum (also the table's LOW part):
  G2 κ_P is ×2–4 higher at 5000–6300 K and ρ < 1e-9.
- **6300 K < T < 5e4 K, 1e-10 < ρ < 1e-3** (the table side of the photosphere, and the outer interior): κ_R network 0.34–0.82 ×
  table: missing lines (OPAL), plus Gaunt factors / He; Z0 with lines reproduces OPAL to ≲ 30 % except at log R > −1
  (docs/Z0_EOS_OPACITY.md: ×1.6–6 low at log R 0–1). κ_P network 0.005–0.04 × table: the table's κ_P is the LTE line
  Planck mean (Z0 lines), ×30–200 above the continuum. Not reproducible by a continuum and arguably the wrong coupling
  quantity where the lines are thick (§4.4).
- **ρ ≳ 1e-4, T 2e4–1e5 K** (the interior, still table zones below the HOT branch): p network 0.66–0.73 ×, c_P 0.5–0.7 ×,
  Δu −0.4 … −1.5 kT per particle, x_e 0.04–0.5 × (profile zones 0–48): the network's Saha has neither occupation probabilities nor pressure ionisation; its
  n ≤ 15 partition function keeps H neutral (x_e 0.28 vs 0.46 Saha-ground vs 0.48 SCVH at log ρ −4, log T 4.4). κ_R ×2–50 too
  high there (bf from overpopulated H⁰, all of it counted as ground state, plus the means grid ending at 12.4 eV).
- **ρ ≳ 0.1, T < 2e4 K** (pressure-ionisation corner, SCVH's interpolated region): p ×0.1–0.5 — off the current profile.
- **SCVH95 composition misread** (`scripts/g2audit/scvhtest.jl`): on SCVH's own isotherms at log ρ −8 … −2 and x_H⁺ 1e-3 … 0.1,
  the table's x_H⁺ = (1 − X_H₂ − X_H)/(1 + X_H₂) is 1.85–1.94 × Saha; reading X as concentrations among all particles
  including electrons (X_H⁺ = X_e = (1 − X_H₂ − X_H)/2) gives 0.97–1.05 ×. [interp; SCVH95's table.doc not obtained — to verify.]
  Between the isotherms (Δ log T 0.08) the linear interpolation of X in log T leaves stripes of ×1.3–3 in x_e (both readings;
  visible in every x_e panel). p and u are unaffected (SCVH's own columns).

## 3. Prototypes on the branch (default off; off = bit-identical)

### 3.1 `G2OPAC=1` — `julia/G2Opacity.jl`
`G2Opacity.means(ρ, T; n1, n2, nH2, nHe, nHep, nHepp, nHp, ne, nHm, nHmS, hbfP, Tw)` returns (κ_E, κ_P, κ_R) with the
conventions of `nlte_means` and, from Z0Physics' own functions evaluated on the network's populations: grid to 25 kT/hc at 3e5 K
(1500 points; includes the Lyman continuum); n1 = all neutral H, ground state = rest (n ≥ 2 from n2 × Boltzmann × HM88 w_n);
H bf with Gaunt factors (Stobbe n = 1, Seaton n ≥ 2) and dissolved-level pseudo-continua (Hubeny, Hummer & Lanz 1994);
He I ground bf (Verner et al. 1996), He II bf (hydrogenic, Gaunt, w_n); H⁺, He⁺, He²⁺ ff with van Hoof et al. (2014) Gaunt
factors; Rayleigh capped at 1300 Å. `G2Opacity.n2_lte` gives the LTE n = 2 with occupation probabilities. With `G2OPAC=1`
`nlte_means`/`nlte_means_fast` call it, passing the zone's He⁺/He²⁺ (and Saha He ions for the equilibrium base). Not tried in
a run (cost, §4.5).

### 3.2 `SCVHXE=1` — `PopIIIEqEOS`
SCVH composition per particle incl. electrons (§2.3); table cached as `…_he_xe.jls` (or `EQEOS_CACHE=path`). Effect: x_e of the
table zones ×0.5–0.6 at 9000 K–2e4 K, ρ 1e-7–1e-4 (zone 61: 0.067 → 0.038; Saha 0.031); p, u unchanged. Not tried in a run.

### 3.3 Bit identity
`scripts/g2audit/bitrun.sh`: zdiagB's S6 settings (procedure.log), RESTART=z19zdiagBS5, MAXSTEPS=8, 8 threads on cores 96–103;
`t_g2a` (0781136) vs `t_g2b` (this branch, switches off): decompressed state (591 arrays, 10 530 280 bytes) and .aux (478 907
bytes) byte-identical (`bitcmp.sh`).

## 4. Plan: what the network needs to reproduce the table

Tolerance target in LTE: 10 % in p, u (in kT per particle), c_P, κ_R, κ_P for ρ ≤ 1e-5 at all T < 3e5 K; 30 % above; the network
stays non-LTE wherever its rates say so (nothing below forces LTE).

**P1 — H I, He I, He II partition functions with occupation probabilities (Hummer & Mihalas 1988), and the
excitation energy in u.** Z_H(T, n_ion, n_neu) = Σ_n 2n² w_n e^{−E_n/kT} replaces ZH_el (and the He analogues) in ChemNet's
Keq and the energy u += x_H ⟨E_exc⟩; plus HM88's non-ideal pressure term (∂ ln w/∂V), so p, u and Keq stay consistent
(Mihalas, Hummer & Däppen 1988 give the algorithm). Fixes the neutral overpopulation above 3e4 K (κ ×2–50, Δu 1–3 kT) and
removes the n = 15 cut. Effect measured with Z0Physics.opacity_state (the same HM88 Saha): along the zdiagB profile x_e is
1.01–1.25 × SCVH(XE) for ρ < 1e-1.5 (network as is: 0.11–0.94 ×) and 1.6–1.7 × at ρ ≈ 0.3, T ≈ 1e5 K (network: 0.04 ×);
off the profile, ρ > 1e-3 and T < 3e4 K, still ×0.1–0.5 (pressure ionisation is SCVH physics, not HM88's hard-sphere w; P4). Cost: Keq becomes
density-dependent for 3 reactions (H, He, He⁺ ionisation and their photo/Milne inverses): ~16 exp + 16 w_n per T and zone,
~2 µs; breaks the per-T memo of `_spterm` for those species only. Duals: smooth (w_n is an exponential of densities), exact.

**P2 — the G2 continuum (G2Opacity) as the network's opacity everywhere, made cheap.** Now 22–33 ms per call
(`timing.jl`) vs ~0.5 ms for `nlte_means_fast` (NLOIT's 485 341 calls took 224 s in 8 zdiagB steps): not usable as is.
Restructure as per-absorber kernels: κ_P is linear in the populations, κ_P = Σ_s (n_s/ρ) σ̄_P,s(T) (+ ff ∝ n_e n_ion g_ff(T)),
with σ̄_P,s(T) tabulated in ln T (H n = 1 … 16, He I, He II levels, H⁻ bf, H⁻ ff, H⁺/He⁺/He²⁺ ff, CIA pairs) — ~20 table
lookups, ~1 µs, exact derivatives in populations and (via the cubic table) in T. κ_R is not linear: opacity sampling on
~150 fixed log-ν nodes with per-absorber σ_s(ν_j, T) tables (precomputed, bilinear in ln T) and the occupation-probability
factors applied per level: ~150 × (20 multiply-adds) ≈ 10–20 µs; in a 2–4-partial Dual ~5×. Fits the implicit Newton (NLOIT
now refreshes slopes by finite differences; this gives them exactly). Accuracy of 150-node sampling vs the 1500-point
quadrature to be checked (< 3 % expected away from edges; edges as nodes, as Z0Physics does).

**P3 — lines.** Cannot be computed per zone (Z0Physics' Stark-broadened line lists: 0.04–4 s per point) and the only complete
source in this regime is OPAL (LTE). Proposal: an LTE line factor per absorber, f_R(ρ, T) = κ_R(OPAL)/κ_R(G2 continuum, LTE
populations) and the line Planck mean σ̄_P,lines per H atom in level n (Z0 line strengths), tabulated once; in a network zone
κ_R = κ_R,cont(x) × f_R, κ_P += Σ_n n_n σ̄_P,lines,n(T, n_e). This keeps the non-LTE continuum, puts the lines' level populations on
the network's n = 1, 2 (three-level atom) and makes the LTE limit the table to the accuracy of f_R's tabulation — the same
assumption the current NLOIT ratio makes above 6000 K, but continuous in T (no NLOITW band) and per absorber. Whether a line
Planck mean should couple gas and radiation where lines are thick (κ_P ×30–200 over the continuum at 6300 K–5e4 K) is a
separate question (G4: multigroup/escape probabilities); the grey TITAN closure inherits the table's choice.

**P4 — non-ideal excess free energy.** ΔF(ρ, T) = F_SCVH(ρ, T) − F_net,LTE(ρ, T) (network with P1, at its LTE composition)
tabulated, smoothly → 0 where |ΔF| < 1e-3 kT/particle (ρ ≲ 1e-5 at T < 1e5 K by the p map); network EOS p = p_net(x) +
ρ² ∂ΔF/∂ρ, u = u_net(x) + ΔF − T ∂ΔF/∂T. In LTE this reproduces SCVH exactly; in the photosphere it is zero, so the
non-LTE network is untouched. Composition: the ionisation equilibrium needs the matching shift of the chemical potentials,
K_eff = K exp(−Δμ/kT) with Δμ from ∂ΔF/∂x (SCVH gives only the equilibrium mixture: Δμ must be fitted as Δμ_H⁺(ρ, T) so
that the network's Saha reproduces SCVH(XE)'s x_e — one table). Requires the SCVH interpolation stripes (×1.3–3 in x_e between
isotherms) to be smoothed first (interpolate log X, or rebuild from the free energy). Cost: two bicubic lookups.

**P5 — small items:** He Rayleigh, H₂⁻ and He⁻ ff (cool dense gas, ρ > 1e-6, T < 3000 K; sources to verify), Li electrons in
the table EOS (consistency only).

### 4.5 Cost per zone evaluation (measured, one core)
table κ (R, P) 1.0 µs, table EOS 0.2 µs; network EOS 2.2 µs, network Saha solve (bisection) 37 µs; network continuum
means (`nlte_means` arithmetic, direct CIA) 1.7–5.4 ms, NLOIT's `nlte_means_fast` 0.46 ms (run average); G2Opacity.means
22–33 ms. Target with P2: ≈ 20 µs for (κ_E, κ_P, κ_R) with exact first derivatives — 20× cheaper than today's network opacity.

### 4.6 What cannot reasonably be reproduced, and the fallback
- OPAL's line blanketing in κ_R at 6300 K–5e4 K and ρ > 1e-8 (G2 continuum on HM88 populations 0.1–0.5 × OPAL) — only by
  calibration (P3).
- SCVH's pressure-ionisation region (ρ > 1e-3) — only by tabulated ΔF (P4).
Fallback without a switch: a weight evaluated at the Newton iterate, w(ρ, T) = smoothstep in log n_H (e.g. 1e17 → 1e19 cm⁻³)
× smoothstep in ln T, generic in the number type; p, u = (1 − w) network(x) + w table, ln κ = (1 − w) ln κ_net + w ln κ_tab;
the network's chemistry keeps running in every zone (in the w → 1 zones relaxed to LTE by its own rates at high density). This
is continuous in the iterate (no per-step discrete choice, no hysteresis, no `x_from_table` reset), energy-conserving through
u(ρ, T, x), and differentiable. Its cost is the inconsistency of a blend inside the band (as the HOT branch blend, PopIIIMaterial
`HOT`); with P1–P2 in place the band only has to bridge lines and non-ideal terms, both of which are ≤ 30 % at n_H ≈ 1e17.

## 5. Sources (checked on SciX, all refereed journals)
Saumon, Chabrier & Van Horn 1995, ApJS 99, 713 (1995ApJS...99..713S); Hummer & Mihalas 1988, ApJ 331, 794 (1988ApJ...331..794H);
Mihalas, Hummer & Däppen 1988, ApJ 331, 815 (1988ApJ...331..815M); Hubeny, Hummer & Lanz 1994, A&A 282, 151 (1994A&A...282..151H);
van Hoof et al. 2014, MNRAS 444, 420 (2014MNRAS.444..420V); Verner et al. 1996, ApJ 465, 487 (1996ApJ...465..487V);
John 1988, A&A 193, 189 (1988A&A...193..189J); Iglesias & Rogers 1996, ApJ 464, 943 (1996ApJ...464..943I);
Lenzuni, Chernoff & Salpeter 1991, ApJS 76, 759 (1991ApJS...76..759L); Mayer & Duschl 2005, MNRAS 358, 614 (2005MNRAS.358..614M);
Ripamonti et al. 2002, MNRAS 334, 401 (2002MNRAS.334..401R); Omukai & Palla 2003, ApJ 589, 677 (2003ApJ...589..677O);
Stancil 1994, ApJ 430, 360 (1994ApJ...430..360S; H₂⁺/He₂⁺ continua); Vidal, Cooper & Smith 1973, ApJS 25, 37 (1973ApJS...25...37V);
Lemke 1997, A&AS 122, 285 (1997A&AS..122..285L); Tremblay & Bergeron 2009, ApJ 696, 1755 (2009ApJ...696.1755T; Stark profiles with HM88);
Borysow, Jørgensen & Zheng 1997, A&A 324, 185 (1997A&A...324..185B); Gustafsson & Frommhold 2003, A&A 400, 1161 (2003A&A...400.1161G).
To verify: SCVH95 table.doc (column definition, §2.3); He⁻ and H₂⁻ free–free sources (P5). No formula in the prototypes is
new: each is Z0Physics' existing implementation (docs/Z0_EOS_OPACITY.md), applied to the network's populations.
