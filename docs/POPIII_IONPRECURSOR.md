# Ionizing photons ahead of the photosphere: why x_e does not respond (t_c128c, 14 Nov 2026)

Analysis note (post-processing only). Script: `julia/ionprec.jl` (zone dump: `julia/ionprec_zones.jl`). Snapshots t_c128c steps 2202650
(M★ 454 M☉, R★ 11.3 R☉, grey r_ph 782 R☉, T_eff(starlight) 10 300 K) and 2192000 (r_ph 1123 R☉, T_eff 9200 K); outputs in
`/zpool/nvme/data/tabel_scratch/titan_popiii/analysis/ionprecursor/`.

## Where the ionizing photons are (step 2202650)

- Source: the hot, fully ionized, LTE (table) post-shock envelope, 13–380 R☉, 1.4e5 → 1.7e4 K, Thomson τ ≈ 20, ε = κ/χ ≈ 1e-2 at 912 Å.
  Q_H(r) (photons/s above 13.6 eV through the sphere, formal solution) peaks at 2.9e50 s⁻¹ at 157 R☉ (L_EUV 2.4e6 L☉, 18 % of L).
  Neither the 1.1e6 K inner boundary (Q 2e47 through 11.3 R☉) nor the edge drives it.
- Sink: an H II region boundary in table zones 76–79 (380–1360 R☉), *at* the grey photosphere. Case-B recombinations per zone
  α_B n_e n_p V = 5.8e49 (76), 6.5e49 (77), 7.2e49 (78), 3.4e49 (79) s⁻¹. Q_H: 1.7e50 (380 R☉) → 5.2e49 (545) → 3.5e48 (r_ph 782;
  explorer cache 2.6e47, SpecFull opacity 1.3e47) → 2.5e45 (1123) → 1.2e42 (1667 R☉, first network zone).
  He I: 3.3e48 at r_ph in the explorer opacity (no He opacity there), 1.6e40 with He I bound–free; He II ≤ 3e47 / 4e29.
- τ_912 (absorption), r_ph outward: 614 across zone 78, 3.1e3 to 1360 R☉, 2.4e4 to the first network zone, 1e5 at 1.2e4 R☉,
  1.3e6 to the edge. Mean free path at r_ph ≈ 0.5 R☉ (zone width 341 R☉). In the front J_912/B(T) = 0.76–1.3: Saha (the table)
  is self-consistent there.
- Edge: Q_H 6e36–1.3e37 s⁻¹ (1e-7 L☉), the same in both snapshots: recombination emission of the cold outer envelope's residual
  x_e ~ 1e-4, not starlight. The explorer's r(τ_ν = 2/3) below 912 Å is the mesh edge (2e9 R☉), because the infall is opaque.

So the "precursor" in the spectrum is the flux through the r_ph sphere into the absorbing front. It is absorbed within ~1 R☉ of r_ph.
It is not a region outside the photosphere where J_ν(>13.6 eV) matters.

## What the network does with photons above 13.6 eV (STARLIGHT=1 SL3=1 HECHEM=1, no RAY, no LOCALJ)

| process | in the network zones ahead of r_ph |
|---|---|
| H(1s) photoionization | **absent**. The tl3 zones drop the H photo pair (`ph3 = filter(p -> p.R !== :H, m.photos)`, chem_core!). Elsewhere that pair sees only the CMB env (TCMB 54.5 K). IonRay (`RAY=1`, the only Lyman-continuum path) is off, and RAYTMIN 1.5e4 K would disable it at this T_eff anyway. |
| He I / He II photoionization | Photo processes exist (Verner σ, Milne inverse) but see the static CMB env: zero. LOCALJ=0. |
| H recombination | three-level C α_B (RECFAST fit): case B, 1s captures on the spot; C ≈ 2e-4–3e-4 in zones 81–87 (the Balmer continuum re-ionizes most n = 2 atoms). |
| He recombination | Milne inverse over all levels ≈ case A. He recombination photons (≥ 24.6 eV, He II Lyα) are not given to H. |
| Lyman-continuum recombination photons | on the spot (case B), no transport. Lyα escapes via Sobolev R_Lyα in C. |

## Rates from the formal-solution J_ν (zones 81–128, both opacity models agree)

Γ_gross(1s) = 4π∫σJ/hν is 0.78 s⁻¹ (zone 81) → 3e-2 (83) → 1e-5 (89). That is the gas's own Lyman-continuum recombination
radiation, J ≈ S = b B. Zone 81 is over-ionized relative to Saha at T (b ≈ 49, J_912/B = 39). Its absorption cancels the 1s captures,
as case B assumes. The net rate Γ_net = 4π∫σ(J − S_th)/hν is |Γ_net| ≤ 4e-8 s⁻¹. That is ≤ 5e-4 of the network's own ionization
rate, which is three-level (1 − C)(n_e q12 + P12) = 4.5e-4 s⁻¹ in zone 81 (Balmer continuum in the diluted ~1e4 K starlight).
Equilibrium x_p with Γ_net added: Δx_e ≤ 1.5e-6 (zone 82), ≤ 1e-8 beyond 6600 R☉. He I: Γ ≤ 2e-8 s⁻¹ against α_He n_e = 2.5e-3 s⁻¹,
so x_He⁺ ≤ 7e-7 and Δx_e ≤ 7e-7. Step 2192000 gives the same picture (|Γ_net/I| ≤ 1.4e-5, Δx_e ≤ 3e-7). The port of the three-level
rates reproduces the snapshot's x_p to 0–20 % in zones 81–100.

## Artefact checks

- Emission at the wrong T: no. The EUV comes from LTE table zones with T_gas = T_r. The inner boundary contributes < 1e-3.
- Missing H(1s) bound–free: no. Both snapspec and SpecFull include it with the snapshot's n_H0. Pitfall found: a frequency grid
  point exactly at 13.598 eV is transparent in SpecFull.FOpac, because the model atom's threshold is 13.5984 eV. SpecFull's own
  nugrid brackets its thresholds by ±0.1 % and is unaffected. An ad-hoc grid leaked 1e47 s⁻¹ to the edge.
- LTE source functions in non-LTE gas: no for H. snapspec's bound–free emission is Milne with the zone's own n_e n_p (b ≠ 1 is
  carried). The explorer opacity has no He opacity, so its He I photon count at r_ph (3e48) is an overestimate by 1e8.

## If it ever matters

Monitor Q_H through the first network zone against its α_B n_e n_p V. The EUV solve of `ionprec.jl` (a) costs ~3 s per snapshot.
If the ratio exceeds ~1e-2 (a hotter star, or an H II region breaking out past the table zones), the clean path is a photon-conserving
ray of the Lyman continuum: the existing IonRay, sourced with the formal-solution Q_ν at the table/network boundary instead of a
T_eff blackbody. It is case B already; it serializes the chemistry sweep (CHEMTHREADS off with a ray). A 4th FLD group above 13.6 eV
(radgroups-3g: the 3-group step costs +6 %) has a stiff group opacity (x_H from 1e-9 to 1 across the front, τ_912 ~ 1e3 per zone)
and cannot place a front thinner than a zone. Feeding the network J_ν from the formal solution needs Γ_net (J − S), not Γ_gross,
or it double-counts the on-the-spot photons.

## IONFS=1 test from t_c128c B7_auto1 (554 M☉, 14 Nov evening)

Switch: IONFS=1 (julia/IonFS.jl; default off; t_ion0 with it off is bit-identical to production B8, 6/6 and 12/12 snapshots).
t_ion1 (production run.env/phaseB.env + IONFS=1), 2040 steps (2255473–2257514):
- Source at the grey photosphere (r_b 129–185 R☉): Q_H 1.2–3.0e49, Q_HeI ≈ 4.8e48, Q_HeII ≈ 6e45 s⁻¹ (L_run/L_fs 1.2–2.4).
  Q reaching 1e3 R☉: 1.6–2.4e49; 1e4 R☉: 1e-3–1.1e49; nothing beyond 4.6e4 R☉ (zone 93).
- He I front 3.3e3 → 1.9e4 R☉; T +600…+1700 K between 130 and 1.9e4 R☉; He²⁺ → He⁺ inside 3e3 R☉. The He²⁺ there was
  LOCALJ's unattenuated W B(T_c); the ray brings only 6e45 He II photons. H: x_H = 0.5 unchanged at 3.3e5 R☉ (ionized by the three-level
  Balmer process).
- Δt 0.02–300 s (production 1e6–1e7 s), limited by T at zones 73–79 (just outside the photosphere); 1.4 Newton divergences per step
  (production 0.9), 246 chemistry fallbacks, no abort. Photon bookkeeping ε ≤ 1.5e-3. Simulated time 7 h: Ṁ unchanged
  (1.20e-4 / 1.22e-4 / 1.13e-4 / 1.13e-4 M☉/yr at 1e2/1e3/1e4/1e5 R☉).
- Strömgren budget (stromgren.jl, production step 2261350): Q = 2.6e49 from 185 R☉ is used up by case-B recombinations at 5.1e5 R☉
  (0.011 pc) with the present ionization, or 1.25e6 R☉ (0.028 pc) fully ionized. Breaking out of the mesh needs Q ≳ 4.1e50. No
  break-out past the 7000 K front (0.13 pc) at this Q.
