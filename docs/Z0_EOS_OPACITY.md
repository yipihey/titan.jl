# Z = 0 EOS and opacity: AD-ready API (v2)

The intended consumer is the Julia TITAN port. Primary variables are (ρ, T) in cgs units. Composition is
X = 0.76, Y = 0.24, Z = 0.

```julia
push!(LOAD_PATH, "/home/claude/QuasiStarInfall/src"); using QuasiStarInfall
tab = load_table()                 # output/kaptable_Z0_v2.jls (Rosseland + Planck)

st = eos_rhoT(ρ, T)                # NamedTuple, 2.8 μs/call
st.P, st.e, st.s                   # pressure (gas+radiation), specific internal energy, specific entropy
st.dPdρ, st.dPdT                   # ∂P/∂ρ|T, ∂P/∂T|ρ        (exact, AD)
st.dedρ, st.cV                     # ∂e/∂ρ|T, ∂e/∂T|ρ
st.χρ, st.χT, st.Γ1, st.∇ad, st.δ, st.cP, st.cs, st.β

κR = kap(tab, ρ, T)                # Rosseland mean (absorption + scattering), cm²/g, 0.1 μs/call
κP = kapP(tab, ρ, T)               # Planck mean (absorption only), cm²/g
in_table(tab, ρ, T)                # false ⇒ arguments were clamped (κ held fixed, derivative 0)
```

## What is guaranteed
| Property | How it was checked (scripts/test_ad_eos*.jl, test_ad_opac.jl) | Result |
|---|---|---|
| AD derivatives of P, e, s | vs centred finite differences, 10 (ρ,T) points from 6×10³ to 10⁷ K | agree to FD accuracy (≲1e-8) |
| ρ(P,T) derivatives | vs the exact identities ∂lnρ/∂lnP = 1/χρ and ∂lnρ/∂lnT = −χT/χρ | 1e-15 to 5e-10 (the worst at β = 6e-6, limited by conditioning; FD is wrong there) |
| Thermodynamic consistency | ∂e/∂ρ\|T = (P − T∂P/∂T\|ρ)/ρ² and ∂s/∂T\|ρ = c_V/T | 3e-15, 2e-15 |
| Nested AD (2nd derivatives) | Hessians of ln P and e; d(Γ1)/dT and d(c_V)/dρ through `eos_rhoT` | ≲1e-9 vs FD of AD |
| Eddington Γ1(β) | analytic formula, β = 0.5 … 0.007 | ≤4e-6 (the residual is real Saha recombination at β = 0.5, ρ ≈ 0.02) |
| κ derivatives continuous | jump in d ln κ_R/d ln T across cell edges | bicubic ≤ 1e-6 (bilinear v1: up to 58–116 %) |
| Interpolation accuracy | vs direct calculation at random off-grid points | κ_R ≤ 1.4e-5, κ_P ≤ 3.2e-5 (max) |

## How it works
* Charge balance (n_e) and ρ(P,T) are solved in Float64 by bisection or Newton. The result is then polished with
  two full Newton steps carried out in the caller's number type, which makes derivatives exact through third order.
* `eos_rhoT` seeds one 2-partial Dual in (ln ρ, ln T) and runs a single Saha solve for P, e and s. Its tag is
  registered via `ForwardDiff.Tag(...)`, so it nests correctly inside an outer ForwardDiff Jacobian.
* The opacity table covers log T = 3.3–8.0 (Δ = 0.005) and log R = −12 to 1 (Δ = 0.025), with R = ρ/T6³.
  Interpolation is Catmull–Rom bicubic in log κ, which is C¹. The frequency quadrature includes every
  bound-free edge as a breakpoint; without that, κ_P had about 1 % noise.

## Opacity physics (rewritten 2026-11; the tables before that were wrong, see below)
`z0/opacity.jl`, data and profiles in `z0/atomic.jl`:
* **Populations.** Saha carries the ground-state weight only, so its neutral density is the ground-state density; excited levels follow from
  Boltzmann with Hummer & Mihalas 1988 occupation probabilities, and `opacity_state` makes the ionization balance consistent with the resulting
  partition functions (opacities only; the EOS keeps plain Saha).
* **Bug fixed.** `level_pops` used to normalise that ground-state density over 12 (8) levels, which removed most of the ground-state opacity when
  the partition sum was large: κ_R was too low by up to ~8× around log T 4.5–5.5, and κ_P by 2–4× in the same region. The earlier statement that "lines"
  explain the difference from OPAL was wrong for most of that range.
* **Free-free:** van Hoof et al. 2014 merged relativistic Gaunt factors (`gauntff_merged_Z0{1,2}.dat`, their licence header kept; data.nublado.org/gauntff).
* **Bound-free:** n = 1 exact (Stobbe), n ≥ 2 Seaton 1960; He I ground Verner et al. 1996, excited levels hydrogenic in n*; dissolved-level pseudo-continua
  (Hubeny, Hummer & Lanz 1994).
* **Lines:** H I and He II with exact n-averaged f (Gordon 1929; TRK sum to 5e-4), exact linear Stark patterns for n ≤ 8, n' ≤ 16 (a model beyond), Holtsmark ions
  (no Debye screening, no ion dynamics), electrons as impact core and quasi-static wing, Doppler, radiative damping, far-red (ν/ν₀)⁴. He I: NIST ASD multiplets
  (`hei_multiplets.tsv`), Doppler + radiative + impact, no ion Stark.
* **Calibrated:** the electron-impact width carries one factor, `IMPACT_CAL = 0.2`, set against the Gigosos et al. 2003 Hα/Hβ FWHM fits (Hβ to 3 %, Hα to 20 %, n_e = 1e16–1e18). Those fit coefficients
  were taken from memory; the paper was found, the numbers were not checked against its text.
* **Table.** κ_P and the continuum κ_R on every node; the line correction to κ_R on a coarse grid (Δlog T 0.05, Δlog R 0.25, log T 3.4–5.8), interpolated.

Against OPAL (recomputed κ_R with lines, X = 0.76, `stellaropacity_test.jl: z0_vs_OPAL`): ≲ 17 % at log R ≤ −2; ×1.3–1.7 low at log R = −1, log T 4.5–5; ×1.6–6 low at log R = 0–1, log T 4.25–4.75
(dense, strongly coupled gas). OPAL remains the κ_R of the stellar runs. Line κ_P is the LTE definition; with optically thick lines it overstates grey gas–radiation coupling.

## Physics limitations
* No molecules: H₂ matters only below ~3000 K at these densities.
* LTE level populations. The Planck mean is meant for the gas–radiation energy coupling and is not an NLTE quantity.
* No degeneracy or pressure ionization, so the EOS should not be used above ρ ~ 10⁻² g cm⁻³.
