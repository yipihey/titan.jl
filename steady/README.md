# QuasiStarInfall.jl

Quasi-static, spherically symmetric quasi-star envelopes with a spherical free-fall accretion shock as the outer
boundary. Purpose: test whether the ~6000 K LTE photospheres of the MESA-QUEST quasi-star models
(Santarelli, Bellinger et al. 2026, ApJL 998 L4) survive external infall.

## v2 (AD-ready EOS/opacity)
See `API_EOS_OPACITY.md`: `eos_rhoT(ρ,T)` (exact AD derivatives, 2.8 μs), `kap`/`kapP` (Rosseland/Planck, C¹ bicubic), edge-aware frequency quadrature. Steady-state results change by ≤0.4 % vs v1 (`scripts/regress_v2.jl`). Tests: `scripts/test_ad_eos.jl`, `test_ad_eos2.jl`, `test_ad_opac.jl`.

## Physics (what is in, what is not)
* **EOS** (`src/eos.jl`): Z = 0 (X = 0.76, Y = 0.24 as in the MESA-QUEST inlist), ideal gas with LTE Saha ionization of
  H, He I, He II, plus blackbody radiation. Thermodynamic derivatives come from the Sackur–Tetrode entropy, taken in
  (ln ρ, ln T). Checked against the Eddington standard-model Γ1(β) to 1e-5, and against the ionized, neutral and
  radiation-dominated limits. No molecules, no degeneracy, no pressure ionization (all negligible here above ~3500 K).
* **Opacity** (`src/opacity.jl`): Rosseland mean computed from continuum sources only: Thomson, H Rayleigh,
  hydrogenic H I and He II bound-free (12 and 8 levels), He I ground bound-free, free-free, and H⁻ bound-free and
  free-free (John 1988). **No bound–bound lines**, so κ_R is a lower bound near 1e4–1e5 K. Tabulated on
  log T = 3.3–8.0 and log R = −12 to 1 (`scripts/build_kap.jl`, about 2 min).
* **Envelope** (`src/envelope.jl`): hydrostatic, with constant L = L_Edd(κ_es, M) from the black hole, as in
  MESA-QUEST. Local MLT with α = 0.8 (Kippenhahn & Weigert cubic). Point mass M_BH at the centre. The shooting
  condition is m(ε_in R) = M_BH; R is insensitive to ε_in (< 5 % over 1e-4 to 1e-2).
* **Outer boundary** (`src/boundary.jl`):
  * Ṁ = 0: Eddington grey photosphere, T = T_eff, P = (2/3)(g/κ + F/c). This is MESA's default T–τ boundary.
  * Ṁ > 0: steady free-fall inflow with ρ ∝ r^{-3/2} and a radiative shock that is thin compared with R (checked:
    l_cool/R < 0.04 for Ṁ ≥ 0.1 M⊙/yr). The boundary pressure is P = ρ_in v_ff² + (2/3) g/κ + (F/c)(τ_in + 2/3).
    The boundary temperature is T⁴ = ¾ T_L⁴ (τ_in + 4/3), with τ_in from the LTE opacity of the infall and the Lucy
    (1971) T(r). This reduces exactly to the Ṁ = 0 case.
  * A stream covering a fraction f of the surface at rate Ṁ has, at its footprint, the conditions of the spherical
    model at Ṁ/f.
* Steady state is justified because t_KH ≈ 5×10³ yr ≪ M/Ṁ.

## Known limitations
1D; MLT rather than TDC; no turbulent pressure; LTE throughout, including the infall (it ignores photoionization by
the shock's X-rays, see `scripts/run_grid.jl` diagnostics); no lines in κ; the radiation force on the infall is
ignored, although L ≥ L_Edd(κ_es) means *ionized* infall cannot be in free fall. Treat the high-Ṁ (compact) branch
as indicative only. Its consistency needs a dynamical radiation-hydrodynamics calculation (TITAN port).

## Running
```
julia -t 2 scripts/build_kap.jl     # opacity table
julia scripts/validate2.jl           # Ṁ = 0 models vs. MESA-QUEST
julia scripts/run_grid.jl            # infall grid + diagnostics → output/infall_grid.csv
julia scripts/branches.jl            # all steady roots R(Ṁ) → output/branches.csv
julia scripts/plot_grid.jl           # figures (CairoMakie)
```
