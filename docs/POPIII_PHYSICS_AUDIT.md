# Chemistry and radiation physics of the Pop III runs: audit of 10 Nov 2026

What the production configuration (`scripts/popiii_run/phaseA.env`, `stageA3.env`, `phaseB.env`; `fx3`, `full3`, `full4`,
`full6`) actually computes, process by process, read from the code paths those flags select. Switches that exist but are
off in every production run: `RAY` (the ionizing ray, `IonRay.jl`), `SLH`, `LYART` (`LyaTransfer.jl`), `SHOCKRAD`,
`NLOPACJ`, `PSWRAD`; `TX.COMPTON` is never set by `popiii_m2c.jl`. `Spectrum.jl`, `snapspec.jl` and `NLTEH.jl` are
post-processing and feed nothing back.

Short answer: the energy bookkeeping of what is there is consistent (no sign, unit or statistical-weight error found), but
the model is not complete. Photoionization from the ground state and Lyman-continuum transport are absent, the radiation
force is the grey flux term alone, and the Lyman-Werner rate is one coefficient times the Wien tail of a grey colour
temperature, independent of the H₂ level populations.

## Structure

- Zone types (`PopIIIMaterial.jl:774, 795`): table (LTE) zones where p ≥ 300 dyn cm⁻² and T ≥ 4500 K, or T ≥ 5000 K and
  τ ≥ 3, with hysteresis; network zones otherwise (non-equilibrium `ChemNet`, three-level H, non-LTE opacity factors); hot
  branch above 3e5 K (`StellarEOS`, `NucNet`, MLT).
- Step (`popiii_m2c.jl:1457–1560`): implicit hydrodynamics with grey radiation at frozen composition and lagged sources,
  then the chemistry at fixed (ρ, u), then the H₂ line transfer, whose cooling enters the next step.
- The star's light in the network zones (`STARLIGHT`) is not a ray: each zone sees J_ν = W B_ν(T_c), with T_c and the
  dilution W taken from the grey E_r (`PopIIIMaterial.jl:596–668`). There is no frequency-dependent attenuation.

## Process by process

| Process | In production | How | Main approximation | Where |
|---|---|---|---|---|
| H photoionization from n = 1 (> 13.6 eV), H II region | no | — | absent: `RAY` off, `SLH=0` gives `pH = nothing`, the three-level atom has no ground-state photo term | `popiii_m2c.jl:372–431`, `IonRay.jl` |
| H ionization through n = 2, collisional ionization, recombination | yes (network zones); table zones are Saha at T_gas | three-level atom: collisional and Lyα-pumped excitation, Balmer-continuum ionization by W B_ν(T_c), case-B recombination with the Peebles factor; heat exchanged with E_r (`RAYQ`) | Sobolev Lyα escape; no diffuse field | `PopIIIMaterial.jl:520–582, 663`; `TitanResidFast.jl:884` |
| Helium | partly | network: inert, neutral, y_He = 1/12. Table: neutral in the Saha part (p ≤ 1e4 dyn cm⁻²), ionized only through SCVH (blend 1e4–1e6) and the hot EOS. The opacity table does include He ⇌ He⁺ | EOS and opacity disagree at low pressure; never photoionized; no HeH⁺ | `PrimordialChem.jl:253`, `PopIIIEqEOS.jl:33–36`, `StellarOpacity.jl` |
| H⁻ | yes | attachment, associative detachment, mutual neutralization, collisional detachment; photodetachment (John 1988) with W B_ν(T_c); photon energy through κ_E c E_r − κ_P 4σT⁴ | the rate uses B_ν(T_c), the absorbed energy a mean weighted at T_gas (`NLOPACJ=0`) | `PopIIIMaterial.jl:262–301, 648` |
| H₂ formation and collisional dissociation | yes | H⁻ and H₂⁺ channels, three-body (PSS83 and H₂, He partners), charge exchange; dissociation as detailed-balance reverse | reverse rates are LTE at all densities; one three-body rate; formation heat fully thermalized | `ChemRates.jl`, `ChemNet.jl` |
| Lyman-Werner dissociation | partly, only outside the photosphere | k = 1.1e8 · 4π W B_ν(12.87 eV, T_c) · f_sh; heat 4.88 eV per dissociation taken from E_r | one coefficient (AAZN97 k27) at one frequency; no dependence on the level populations; fitted shielding (WG11 H₂ with a Sobolev cap on the column, static H I); no UV pumping heat or cascade, no direct continuum dissociation, no H₂ photoionization, no H₂⁺ photodissociation | `PopIIIMaterial.jl:476–493, 650–653` |
| H₂ ro-vibrational cooling | yes | 63 levels (v ≤ 2, J ≤ 20), non-LTE populations, escape probability with ALI, collisions with H, H₂, He, H⁺ | no electron collisions; H and H₂ rate tables end at 6000 K; no incident continuum; no UV or formation pumping; the line luminosity leaves the grid and is not added to E_r; the populations are not used by the chemistry | `H2Levels.jl`, `H2Lines.jl`, `popiii_m2c.jl:437–548` |
| Lyman α | partly | Sobolev escape and two-photon decay inside the three-level atom; the escaping energy goes into the grey E_r of the same zone | no frequency diffusion or trapping, no destruction by H₂ or H⁻, no pressure | `PopIIIMaterial.jl:523–582`; `LyaTransfer.jl` off |
| Continuum transport | yes | grey two-moment (E_r, F_r) with a variable Eddington factor from a grey tangent-ray formal solution; CMB incident at the outer edge | grey; `EDDSCAT` counts all opacity above κ_P as scattering; `FLUXTIE` in thin zones; Compton exchange off | `TitanResidFast.jl:822–926`, `TitanCore.jl:611–783` |
| Radiation force | yes, one term | χ_F F_r / c in the momentum equation, the same term in the radiation momentum row; radiation pressure work through f E div u | χ_F is the harmonic mean of the zones' Rosseland (table) or non-LTE means, weighted at T_gas; no Lyα, line or ionizing-photon force | `TitanResidFast.jl:744, 760, 899–923` |
| Accretion-shock radiation | partly | viscous heat into the gas row, grey exchange | no hard post-shock spectrum, no ionizing precursor (`SHOCKRAD` off) | `TitanResidFast.jl:733–746` |
| Other | | Li chemistry yes; D and HD no (D only in the nuclear network); free-free yes; nuclear burning in the hot branch; no cosmic rays, X-rays or dust (none expected); hot-branch opacity at fixed X = 0.75 whatever the burned composition | | `popiii_m2c.jl:694–699`, `NucNet.jl` |

## Gaps, ranked for a star accreting to 100–300 M☉

1. **No Lyman continuum.** Quiet today: a blackbody at the present photospheres gives 4e38 photons s⁻¹ (`full4`, 8 M☉,
   T_eff 5009 K) and 2e30 s⁻¹ (`fx3`, 103.4 M☉, 2842 K in the expanded phase of the oscillation) against 2e47 hydrogen
   atoms s⁻¹ falling in at 6e-3 M☉/yr. It decides the goal, though: the feedback expected to end the accretion once the
   star contracts and becomes hot is an H II region with the radiation force at its front, and the run has nothing that
   can form one. `IonRay.jl` passes the static and the
   hydrodynamic Iliev tests but its coupling to the grey field is not validated, the H II temperature is open, and it has
   no helium, no diffuse field and no momentum (`POPIII_GOAL.md`, M3).
2. **Radiation force is grey.** In the optically thick ionized layer the Rosseland force is the right one. In thin
   zones the flux mean should be weighted with the spectrum that is there, not with dB/dT at the gas temperature, and
   Lyα pressure has no term at all. The envelope already sits at Γ ≈ 1 in the ionized layer near 101 M☉ (2.4 in one
   zone of the last `fx3` snapshot), so the outcome is sensitive to this.
3. **Helium ionization is missing from the low-pressure equation of state.** In the last `fx3` snapshot (step 288750)
   79 zones, 22.5–3050 R☉, 1.9e4–1.36e5 K, 0.31 M☉, are table zones with p < 1e4 dyn cm⁻²: the whole inflated layer
   between the stellar surface and the photosphere. There the EOS has neutral helium at any temperature while the
   opacity has the helium ionization bumps. Γ₁ and c_P of the He and He⁺ ionization zones are therefore absent in exactly
   the layer that oscillates. Whether that changes the oscillation is not known.
4. **Lyman-Werner.** A constant coefficient evaluated at 12.87 eV. For a 6000 K spectrum the photon number rises by
   e^(ΔE/kT) ≈ 25 from 12.87 to 11.2 eV and by 1e4 to 8 eV, where vibrationally excited H₂ absorbs (v = 1 holds 14 % at
   3000 K). The true rate from a cool star into warm H₂ can exceed the coded one by a large factor; it needs the level
   populations (available from `H2Levels.jl`) and the band oscillator strengths. Whether photodissociation then competes
   with collisional dissociation and H⁻ photodetachment in the envelope has not been evaluated.
5. **H⁻ heating.** The destruction rate uses the star's colour, the absorbed energy the gas temperature. Estimated a
   factor 2–3 too little heat per detachment for 6000 K light on 3000 K gas; this is the dominant envelope heating inside
   5e4 R☉ (`POPIII_RUNS.md`, "What heats the infalling envelope"). `NLOPACJ=1` exists and is untested in production.
6. **Lyα is local** (no trapping, no pressure). 7. **Compton exchange off; table zones use one κ_P for absorption and
   emission.** 8. **Shock spectrum grey** (see "The accretion shock" below: at 7–8 M☉ the shock supplies the whole
   luminosity from optical depth below one). 9. **H₂ cooling limits** (v ≤ 2, no electron collisions, rates frozen above
   6000 K): adequate below about 3000 K. 10. **No D/HD, HeH⁺, H₃⁺.** 11. **Hot-branch opacity at X = 0.75**: irrelevant
   until hydrogen burns.

## The accretion shock (10 Nov)

Zones around the shock in the last snapshots (`shockz.jl`, `shockx.jl` in the session scratchpad):

| | `full3`, `full4` (7.3, 8.1 M☉) | `fx3` (103.6 M☉) |
|---|---|---|
| infall speed before the shock | 111–115 km/s | 162 km/s (522 km/s in the fall-back phase of the oscillation, step 285550) |
| strong-shock temperature 3μm_H u²/16k, μ = 0.6 | 1.7e5–1.8e5 K (15 eV) | 3.6e5 K (31 eV); 3.7e6 K (320 eV) in fall-back |
| τ_R outside the fastest infall | 0.17–0.20 | 750 (319 in fall-back) |
| column outside | 1860 g cm⁻² | 3340 g cm⁻² (2070) |
| P_rad / P_gas ahead of the shock | 0.05 | 149 |
| T_gas / T_rad through the shock | 0.97–1.10 | 1.00 |
| T_gas in the shock zones | 4700–7900 K | 1.5e5–2.4e5 K |
| GM★Ṁ/R★ over L at the photosphere | 1.9–2.0 | 0.19 |
| zones the viscosity spreads the shock over | 5 (Δr 6 R☉ each); τ_R rises from 0.2 to 500 across them | 5 (Δr 1.2 R☉ each) |

- **`fx3`.** The shock is 750 to 1600 Rosseland depths inside, radiation pressure is 150 times the gas pressure, and
  c/|u| = 1850 exceeds τ_R, so radiation diffuses ahead of the shock faster than the flow carries it in. Gas and radiation
  have one temperature on both sides. The strong-shock temperature is not reached and the two-moment grey transport holds
  the structure. Free fall onto the present R★ would be 1450 km/s (2.9e7 K, 2.5 keV); the column outside is 2000–3300 g cm⁻²
  (Thomson depth about 10³), and it scales as Ṁ/√(M r). Reasoning from these numbers, not a calculation: in spherical
  symmetry keV photons could reach the infalling gas only after the accretion rate has fallen by about 10³.
- **`full3`, `full4`.** The shock sits at the photosphere and supplies the whole luminosity. Its upstream side is at
  τ_R = 0.2. Physically this is a Mach 20 gas shock with a post-shock layer at 1.7e5 K; its cooling length is of order
  1e4 cm (estimate: n = 1.4e14 cm⁻³, Λ = 1e-22 erg cm³ s⁻¹, post-shock speed 29 km/s) against zones of 4e11 cm. The model
  has no such layer: the viscous heat is spread over five zones and exchanged with the grey radiation field, and the gas
  stays within 10 % of T_rad. The energy is conserved; the spectrum is not represented. The escaping light is treated as
  a blackbody at the photospheric temperature (5400 K), while about half of it should be the reprocessed emission of a
  1.7e5 K layer seen through τ_R = 0.2. Not evaluated: the emergent spectrum, the change of the photo-rates in the
  envelope (H⁻, H₂⁺, Lyman-Werner), the ionization of the gas just ahead of the shock and its H⁻ opacity.
- No mesh resolves the cooling layer (seven orders of magnitude). A finer mesh does narrow the five-zone numerical shock,
  which in `full3`/`full4` is as wide as the whole transition from τ_R = 0.2 to 500, so the position of the photosphere
  relative to the shock is not independent of the zone size there.
- `SHOCKRAD=1` sends the viscous heat of a strong shock to the radiation energy row instead of the gas row
  (`TitanResidFast.jl:181–189, 883`). It is grey and changes no spectrum. Untested in production.

## Checked and found consistent

`RAYQ` and `GASQ` signs in the gas and radiation energy rows; the LW constant, the WG11 Doppler parameter and exponents,
the H I shielding form; q21 = q12/4; the dilution factor; the CMB boundary intensity; no double counting of H bound-free in
network zones (`hbfP = false`, the three-level atom carries it).

## Not checked

`StellarEOS.jl`, `SCVH.jl`, the opacity table build, `NucNet` rates, MLT; rate fits against the literature (the H Rayleigh
fit is marked unverified in `ContinuumOpacity.jl`); whether the energy zero is identical across the network/table switch;
whether lagging the non-LTE opacity factors conserves energy; what `PSWRAD` does when on; `SHOCKRAD` in a run. Impact statements
other than the numbers quoted in gaps 1 and 3 are reasoning, not measurements.
