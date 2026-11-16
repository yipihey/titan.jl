# G2: one non-LTE material everywhere — audit of table vs network material (16 Nov 2026)

Branch `g2-material` (from `popiii-s6` 0781136). Question (Tom): can the network material (EOS and opacities from the
network's own populations) reproduce the equilibrium-table material where the table is used, so that the per-zone
table/network switch can be dropped? Motivation: in zdiagB (IC abn02red, M* ≈ 0.057 M☉, R* ≈ 10 R☉) the photosphere cliff
that stalls Δt sits at the last table zone / first network zone.

Everything below is measured [run] unless marked. Scripts: `scripts/g2audit/` (standalone; no TITAN needed except the bit test).
Plots: `docs/g2/g2_gap_maps_g2opac.png`, `docs/g2/g2_gap_maps_hm88.png`.

## 0. Summary

**Update (§6):** SCVH95 definition verified; HM88 and the fast G2 continuum implemented; with the network through the photosphere front (table only at n_H ≥ 1e19) the cliff leaves the material interface and the run matches t_sz31 for 300 steps — SUBZONE and SLIMP both remain necessary.

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
  including electrons (X_H⁺ = X_e = (1 − X_H₂ − X_H)/2) gives 0.97–1.05 ×. **Verified (§6.0)**: SCVH95 §5 p. 724 eqs. 33–36 define the
  concentrations exactly so.
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

## 6. Follow-up: steps 0–2 of the plan and the extended-network experiment (16 Nov, PI go-ahead)

### 6.0 SCVH95 column definition — verified
Saumon, Chabrier & Van Horn 1995 (1995ApJS...99..713S), §5 "Calculation of the EOS tables", p. 724: "The concentrations are defined as
X_i = N_i/N (33), where N is the total number of particles, including electrons. With this definition, and considering charge conservation,
we have for pure hydrogen X_e = X_H⁺ = ½(1 − X_H₂ − X_H) (34) … n_i = (2ρ/m_H) X_i/(1 + 3X_H₂ + X_H) (36)" (and eq. 35, 37 for helium).
So per H nucleus x_H⁺ = 2X_H⁺/(1 + 3X_H₂ + X_H), f_H₂ = 4X_H₂/(1 + 3X_H₂ + X_H): exactly `SCVHXE=1` (now written in that form, same
values). The default reading doubles x_H⁺ where H is weakly ionised. (Source: the scanned paper via ADS; table.doc not needed.)

### 6.1 `HM88=1` — occupation probabilities in ChemNet (exact AD)
`ChemNet.HM88` (docstring there): H I partition function Z_w = Σ_{n ≤ 80} 2n² (w_n/w_1) e^{−E_n/kT} with Hummer & Mihalas (1988) eq. 4.71
(as Z0Physics.occupation); network EOS u += x_H ⟨E_exc⟩, p += n_H x_H kT Σ P_n (−ln w_n/w_1) (cross terms through the perturbers' chemical
potentials omitted), exact in ρ, T, x (ForwardDiff = finite differences to 1e-9). Chemistry: every H recombination channel (reactions
H⁺ + e (+ …) → H (+ …) incl. detailed-balance reverses and the three-level C α_B; the Milne inverse of H photoionization) × φ = Z_w/Z_el, so
the network's own equilibrium is the HM88 Saha; φ from the zone's perturbers at the start of the step (`env.hm`, set in `chem_core`), exact
in T. Network LTE x_e with HM88 vs Z0Physics' HM88 Saha and SCVH(XE) (`hmtest.jl`):

| log ρ, log T | network n≤15 | network HM88 | Z0 HM88 Saha | SCVH (SCVHXE) |
|---|---|---|---|---|
| −5.33, 4.26 | 0.460 | 0.493 | 0.493 | 0.488 |
| −4.67, 4.38 | 0.486 | 0.657 | 0.657 | 0.608 |
| −4.00, 4.40 | 0.287 | 0.458 | 0.458 | 0.478 |
| −3.08, 4.78 | 0.249 | 0.884 | 0.939 | 0.761 |
| −2.00, 4.90 | 0.099 | 0.703 | 0.790 | 0.678 |
| −1.00, 5.00 | 0.041 | 0.443 | 0.951 | 0.628 |

(differences to Z0 at ρ ≥ 1e-3: Z0 also has He occupation probabilities.) Off = bit-identical.

### 6.2 `G2OPAC=2` — the fast continuum (`G2Opacity.means_fast`)
Same physics as `means` on 484 frequency nodes (400 log nodes 10 cm⁻¹ … 25 kT/hc at 3e5 K + node pairs at the H, He II (n ≤ 6), He I and
H⁻ thresholds), T-dependent kernels (g_ff Z = 1, 2; H⁻ ff; three CIA pairs) tabulated in log T (Δ 0.01, linear), per-node cross sections
precomputed, trapezoid rule; generic in the number type (ForwardDiff through T and populations). Measured (`fasttest.jl`): **33–37 µs
per call, 50–56 µs with a 2-partial Jacobian** (vs 22–33 ms for `means` and 0.46 ms for today's `nlte_means_fast`); vs `means` over
log ρ −14 … −1, log T 3 … 5.45: κ_P median 0.3 %, 95 % within 4.4 %, max 6.8 %; κ_R median 0.2 %, 95 % within 10 %, max 31 % (in the
CIA band region, log T 3.2, where `means` itself holds the CIA kernel over 8 of its points). AD vs central differences 1e-9 … 1e-2
(the latter where the difference straddles a table node of the piecewise-linear T interpolation). `G2OPAC=1` keeps the 1500-point form.

### 6.3 Bit identity (all switches off)
`t_g2c` (280ba90) vs `t_g2a` (0781136), 8 steps from zdiagB S5: 591 state arrays and .aux byte-identical.

### 6.4 The extended-network experiment (zdiagB S5, step 3000)
Restarts of zdiagB's checkpoint S5 with zdiagB's S6 settings (procedure.log), 8 threads on cores 96–103, RUNREC every 5 steps
(`scripts/g2audit/exprun.sh`, `g2track.jl`, `cmptrack.py`; runs /zpool/nvme/data/tabel_scratch/titan_popiii/runs/t_g2e1 … t_g2e8). G2
switches in all of them: `SCVHXE=1 HM88=1 G2OPAC=2`. **Switch thresholds (general, by density):** `PSWITCH=1e30 PSWTAU=0 PSWSTAY=0
PSWNH=1e19` with the existing `PSWTMIN=4500 PSWHYS=0.3` — the table only where n_H ≥ 1e19 cm⁻³ (ρ ≳ 2e-5 g cm⁻³) and T ≥ 4500 K;
everything else, the whole photosphere front included (front at n_H ≈ 1e17), on the network. `PSWRAMP=10`: zones handed from the table
to the network one every 10 steps after the restart (converting the five zones 59–63 at once failed within 13 steps, E2).

| run | settings beyond G2 + thresholds | result |
|---|---|---|
| E1 t_g2e1 | (no ramp) | stopped 3069 (Δt < DTSTOP); E_r of zone 66, just outside the front — the same failure as zdiagB at 3032 |
| E2 t_g2e2 | SUBZONE=1,2 SLIMP=2, no ramp | stopped 3013: zone 60, converted at once, diverges in E_r, F_r, x |
| E3 t_g2e3 | PSWRAMP=10 | stopped 3032 (as zdiagB: E_r zone 66) |
| E4 t_g2e4 | PSWRAMP=10 SUBZONE SLIMP | stopped 3046: the converted zones 60–62 (τ ≈ 1e5–1e7, 1–2.5e4 K) diverge in E_r, T |
| **E5 t_g2e5** | PSWRAMP=10 SUBZONE=1,2 SLIMP=2 **TL3DEEP=1 NLOHBFDEEP=1** | **300 steps to 3300** (below) |
| E6 t_g2e6 | as E5 without SUBZONE and SLIMP | stopped 3026 (E_r zone 62, the front) |
| E7 t_g2e7 | as E5 without SLIMP | stopped 3026 |
| E8 t_g2e8 | as E5 without SUBZONE | stopped 3043 |

E4 → E5: network zones deep in the thick layers need the n = 2 bound–free exchange in the grey coupling. NLOPAC's κ_P/κ_E leave out H
bound–free (the three-level atom carries it as a lagged chemistry heat): right in the thin lit layers, wrong at τ ≫ 1, where the converted
zones could not exchange energy with E_r inside the Newton and ran away (T 2.5e4 K next to 7.7e3 K at step 3030 of E4). `NLOHBFDEEP=1`
(new, PopIIIMaterial, default off): with `TL3DEEP=1` (no three-level heats at τ ≥ 3) deep network zones take H bound–free in κ_E and κ_P
(hbfP = true), so the exchange is c κ_E E − 4π κ_P B, implicit.

**E5 vs t_sz31** (t_sz31 = SUBZONE=1,2 SLIMP=2 with the table up to the front; same checkpoint; snapshots 3005–3300):

| | E5 (G2 network through the front) | t_sz31 (reference) |
|---|---|---|
| cliff on a table\|network interface | 2 of 60 snapshots (the first two, before the ramp) | 55 of 60 |
| table\|network interface | zone 58 (n_H = 1e19), 5–9 zones inside the cliff | at the cliff |
| Δt from step 3040: median / min / fraction ≥ 1e4 s | 3.4e4 / 1.5e4 s / 100 % | 3.3e4 / 1.3e4 s / 100 % |
| Newton divergences per 100 steps | 2 (6 in 300; DTERR rejections 11) | 2 (6 in 300) |
| wall per step | 3.3 s (8 threads) | 9.7 s (4 threads, t_sz31's own) |
| ⟨log L_ph⟩, σ (3040–3300) | 2.055, 0.256 | 2.031, 0.277 |
| ⟨T_eff,L⟩, σ | 3754 K, 642 K | 3763 K, 674 K |
| M*, R* at step 3300 | 0.0893 M☉, 16.39 R☉ | 0.0916 M☉, 16.39 R☉ (3e-2 yr later in t) |
| cliff ⟨Δ log T⟩ between two zones | 0.41 | 0.42 |

With the network through the front the cliff no longer sits on a material interface; Δt, failures and L, T_eff, R*, M* are those of the
reference. The cliff itself (Δ log T ≈ 0.4 across one zone: the H ionisation front, κ_R ×10⁶) is unchanged: physical and unresolved,
not a product of the switch. Plot: `docs/g2/g2_extended_network_runs.png`. The NLOIT opacity cost per step fell with G2OPAC=2
(E5: 21.4 M ratio evaluations in 1107 thread-s over 300 steps ≈ 52 µs each, vs 460 µs with nlte_means_fast in t_g2a).

### 6.5 Recommendation
- **SUBZONE and SLIMP: both still needed.** Without both every configuration stops within 26–69 steps at the front (E1, E3, E6), with
  the network through the front as with the table up to it; the failure (E_r of the zone just outside the front) is not caused by the
  material interface. Each alone fails (E7 SUBZONE only: 3026; E8 SLIMP only: 3043).
- **G2 material:** SCVHXE + HM88 + G2OPAC=2 + TL3DEEP + NLOHBFDEEP with a density threshold (PSWNH = 1e19; PSWITCH, PSWTAU off) carries
  the photosphere and the front on the network with nothing worse than t_sz31 over 300 steps; it removes the T/τ switch rules (PSWTAU,
  PSWTEX, PSWSTAY, PSWAGREE, PSWRAD) from the front. The remaining density switch sits 2 dex inside the photosphere, where table and
  network (with HM88) agree to 3–8 % in x_e and p. Next: run longer, lower PSWNH (needs the non-ideal terms, P4, for the deep zones).
- SCVHXE is a correctness fix of the table's composition (SCVH95 eqs. 33–36); make it the default after a production comparison.
- Bit identity at 1b4d7db (all switches off): `t_g2f` vs `t_g2a`, 591 arrays + .aux identical.

## 7. The network's validity domain and a smooth table blend (branch g2-domain, 16 Nov)

Goal (Tom): the network material correct by construction everywhere a Pop III protostar simulation reaches at low density, the table only
in the stellar interior, wherever that lies in a given run.

### 7.1 Criterion
Three standard non-ideality measures, evaluated with the LTE (table) composition at (ρ, T):
- Γ_e = e²/(a_e kT), a_e = (3/(4π n_e))^{1/3} — electron Coulomb coupling (Debye–Hückel/plasma non-ideality);
- D = n_e λ_e³/2, λ_e = h/√(2π m_e kT) — electron degeneracy (Boltzmann electrons need D ≪ 1);
- P = (4π/3) n_neutral (2a₀)³ — neutral packing: the hard-sphere exponent of the ground-state occupation probability of Hummer & Mihalas
  (1988, ApJ 331, 794, eq. 4.71; the network already carries HM88 for the excited levels), i.e. the onset of pressure
  dissociation/ionisation of H₂, H, He.
s = ((Γ_e/0.3)⁴ + (D/0.1)⁴ + (P/0.01)⁴)^{1/4}; table weight w = smoothstep(log₁₀(3s)/log₁₀ 3): network for s ≤ 1/3, table for s ≥ 1
(the thresholds are where each parameter's non-ideality becomes significant; the blend window below them).
The thresholds come from the network (HM88 + H2PF, LTE) vs SCVH95 comparison on log ρ −14 … +1, log T 3 … 5.5 (`gridnd.jl`,
`ndion.py`, `ndcrit.py`), binned by each parameter: in ionised gas |Δ ln p| ≤ 0.02 and |Δu| ≤ 0.03 kT up to Γ_e = 0.3, 0.76 at Γ_e 1–3;
in neutral gas the error is 1 % for P < 0.01, 6 % at 0.01–0.03, 20 % at 0.03–0.1; degeneracy 2 % for D < 0.1, 5 % at 0.1–0.3. Result:

| region | points | \|Δ ln p\| p50/p95/max | \|Δu\| [kT] p50/p95/max | \|Δ ln Γ1\| p95/max | \|Δ ln c_P\| p95 |
|---|---|---|---|---|---|
| s < 1/3 (network) | 12873 | 0.000/0.012/0.033 | 0.001/0.19/0.56 | 0.035/0.064 | 0.14 |
| 1/3 ≤ s < 1 (blend) | 1580 | 0.016/0.035/0.045 | 0.032/0.20/0.46 | 0.034/0.064 | 0.19 |
| 1 ≤ s < 3 (table) | 1140 | 0.016/0.080/0.13 | 0.087/0.59/0.66 | 0.089/0.12 | 0.27 |
| s ≥ 3 (table) | 3433 | 0.24/2.8/5.4 | 0.23/0.76/0.94 | 0.58/7.9 | 0.54 |

The residual Δu and Δ c_P where s < 1 are ideal-gas differences, not non-ideality: SCVH's H₂ (and the table's Saha part, RHFC02 A9)
against the network's level-sum H₂ (H2PF, below) at 2500–5000 K, and the SCVH tables' isotherm interpolation (×1.3–3 stripes in x_e).

**`G2BLEND=1`** (PopIIIMaterial, default off, bit-identical off): w at the Newton iterate (generic in the number type: the table
composition through its bicubic interpolant); EOS (u, p) = (1 − w) network + w table through `ChemNet.EOSMIX`, so the chemistry's T(u) and
the hydro's EOS are the same function; NLOIT opacities ln κ = (1 − w) ln κ_network + w ln κ_table; a zone is a table zone exactly
where w ≥ 0.999 (stays to 0.99), where both materials are the same function. Overrides PSWITCH/PSWNH/PSWTAU/PSWTEX/PSWSTAY.

### 7.2 Where it lands
| state | table (w = 1) | blend (0 < w < 1) | network |
|---|---|---|---|
| early collapse f128c A1 (first core forming) | — | — | everything (s ≤ 1e-3) |
| f128c A2 (first core) | log ρ ≥ −2.4 (core: Γ_e 0.6–1.7, D 0.005–0.5) | — (the core edge is a jump to log ρ −4.7) | log ρ ≤ −4.7 |
| zdiagB S5 (0.061 M☉, R* 11 R☉; H front at 12.9 R☉) | r < 1.4 R☉ (log ρ ≥ −2.6, Γ_e ≥ 0.34) | 1.4–8 R☉ (log ρ −2.6 … −5, T 2–8e4 K, τ ≳ 1e8; w 0.99 → 0.6) | r ≳ 8 R☉: the photosphere and the H front (zone 61, Γ_e ≤ 0.1, w ≈ 0) and everything outside |
| z300g2 step 3400 (0.095 M☉) | 35 zones, r < 1.4 R☉ | 33 zones, to ≈ 10 R☉ (w ≤ 0.58 beyond 7 R☉) | the front and outside |
| t_nolte2 late (556 M☉, R* 14.8 R☉) | the star: all zones are HOT (T > 3e5 K, stellar EOS) | — | everything outside the star (ρ < 1e-10: s ≈ 1e-3) |

So in the low-mass protostars the table is complete at Γ_e ≈ 0.3 (log ρ ≈ −2.5, 1e5 K) and fades out where Γ_e ≈ 0.1 — the ionised
interior a few R☉ inside the H front (at τ ≳ 1e8); for the massive star the whole star is the HOT branch and the network has
everything else; the photosphere and the front are on the network in every case (s ≈ 0.15 at zdiagB's front, 1e-3 outside t_nolte2).
EOS along the profiles (network HM88 + H2PF in LTE vs table, `profverify.jl`): in the network domain max |Δ ln p| ≤ 0.015 everywhere;
max |Δu| 0.27–0.28 kT at 2500–4400 K (H₂: network level sum vs table RHFC02 A9 — the network is the better one there; with H2PF off the
two agree, at the cost of a wrong H₂); in the blend zones |Δ ln p| ≤ 0.043, |Δu| ≤ 0.09 kT, |Δ ln Γ1| ≤ 0.06.

### 7.3 Gaps filled / quantified
- **H₂ partition function (new, `H2PF=1`, ChemNet):** the rigid-rotor × harmonic-oscillator Z_int and RHFC02 A9 energy (which uses the
  vibrational *heat capacity* x²eˣ/(eˣ−1)² as an energy) replaced above 1500 K (blended from 800 K; the cold gas keeps RHFC02) by the level sum
  over all 301 bound levels of ExoMol RACPPK (Roueff et al. 2019, 2019A&A...630A..58R), u_H2 = 3/2 + ⟨E_int⟩/kT, consistent with K_H2.
  Network vs SCVH at 3800 K, log ρ −5.6: Δu 0.77 → 0.56 kT, Δ ln p 0.039 → 0.028 (the rest is SCVH's H₂).
- **HM88** (§6.1) is required: without it the network's x_e is 0.04–0.5 × SCVH at Γ_e 0.1–1.
- **Lines (P3):** in the network domain (s < 1) the table/G2-continuum ratio (G2OPAC=2 on HM88 populations, `linefac.py`) is, in dex:
  below 5000 K κ_R −0.11 … +0.05, κ_P −0.44 … +0.03 (no lines; the table's continuum lacks the Lyman continuum); 5000–6300 K κ_P −0.8 … 0
  (same); 6300 K–5e4 K κ_R median +0.03, p95 +0.57, max +0.82 (×6.6) and κ_P median +1.3, p95 +2.4, max +3.3 (H/He lines);
  5e4–3e5 K κ_R ≤ +0.45, κ_P median +0.4. These are carried by NLOIT's LTE calibration κ = κ_net(x) × κ_table/κ_net(LTE): exactly
  P3's line factor, evaluated on the fly (two means_fast calls, ~70 µs) and continuous; a pre-tabulated factor would halve that cost
  (not done). The κ_P line factor is an LTE line Planck mean; whether thick lines should couple grey gas and radiation that strongly is G4.
- **ΔF (P4):** not needed: the blend completes (s = 3) before the network-table difference exceeds ~10 % in p.
- **H₃⁺, H₂⁻ and He⁻ free–free, H₂ lines in the continuum: not quantified** (no data in hand; no formulas from memory). Mayer & Duschl
  (2005, 2005MNRAS.358..614M) find H₃⁺ affects the Rosseland mean at intermediate densities (≈1e-12–1e-9 g cm⁻³, 800–3000 K, through the
  electron budget) and H₂ lines dominate the Planck mean at low density; the table has none of these either, and the H₂ lines are in the
  runs' line transfer (H2Lines), not in the grey κ. To do: H₃⁺ chemistry/partition function (Neale & Tennyson 1995) and the H₂⁻/He⁻ ff
  fits (Bell 1980; John 1994) — to verify on SciX.

### 7.4 Run test from z300g2 (G2BLEND replacing PSWNH=1e19)
z300g2's run settings (procedure.log) + `G2BLEND=1 H2PF=1`, restarted from copies of its checkpoint `z19z300g2S1_auto0`, cores 104–111
(`scripts/g2audit/exprun2.sh`; z300g2 itself untouched).
- `t_g2h` (blend window s 1 … 3, the first version): from step 3400, stopped at 3621 (Δt < DTSTOP): zone 37 (log ρ ≈ −2.5, T ≈ 9e4 K, Γ_e ≈ 0.3,
  s ≈ 1, w ≈ 0), just handed to the network by the restart ramp, diverged in T, E_r, x_H⁺, x_H₂ — the network's coupled chemistry does not
  converge in zones that dense and hot.
- `t_g2i` (window s 1/3 … 1, the default now: the table is complete where a parameter reaches its threshold): from step 3800 (the newest
  checkpoint then), **400 steps to 4199**, table|network interface at zone 36 (Γ_e ≈ 0.3), the H front 36 zones further out on the network.

| steps 3840–4200 | t_g2i (G2BLEND) | z300g2 (PSWNH = 1e19) |
|---|---|---|
| Δt median / min | 1.5e4 / 6.1e3 s (421 s at the restart) | 1.9e4 / 9.8e3 s (50-step snapshots) |
| Newton divergences per 100 steps | 1 (4 in 400) | 0.5 (2 in 400) |
| wall per step (8 threads) | 3.9 s | 2.9 s |
| M*, R* at t = 1.356 yr after 238933 yr | 0.1301 M☉, 22.5 R☉ | ≈ 0.1295 M☉, 22.5 R☉ (interpolated) |
| L_ph, T_eff | oscillating 1–1200 L☉ (as z300g2) | oscillating 0.6–1160 L☉ |
| cliff on a table\|network interface | 0 of 80 | 0 of 21 |

### 7.5 Remaining gaps
- H₃⁺, H₂⁻/He⁻ free–free: not in either material, not quantified (§7.3).
- The network's coupled chemistry fails to converge at Γ_e ≳ 0.3 (log ρ ≳ −2.5 at 1e5 K): the reason the table is kept from s = 1 on; with
  a robust dense-gas chemistry (LTE limit of the kinetics) the window could move to s 1 … 3.
- The EOS comparison in the network domain is limited by SCVH95's own H₂ (0.3–0.6 kT at 2500–5000 K vs the level sum) and its isotherm
  interpolation; a direct check against an independent EOS (e.g. Chabrier et al. 2019 H/He tables) would settle which is right there.
- Lines: carried by NLOIT's LTE calibration; the κ_P line factor (×20–2000 at 6300 K–5e4 K) is an LTE line Planck mean (G4 question).
- G2BLEND was run once over 400 steps from one checkpoint; longer runs and another accretion history are still to do.
- Bit identity at 668e3a9 (all switches off): `t_g2j` vs `t_g2a` (0781136), 8 steps from zdiagB S5: 591 arrays + .aux identical.

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
