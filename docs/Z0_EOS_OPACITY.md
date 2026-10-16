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

## Physics limitations (unchanged)
* No bound–bound lines: κ_R is a lower bound near 10⁴–10⁵ K, where OPAL-type tables are higher by a factor of a few.
* No molecules: H₂ matters only below ~3000 K at these densities.
* Gaunt factors are set to 1.
* LTE level populations. The Planck mean is meant for the gas–radiation energy coupling and is not an NLTE quantity.
* No degeneracy or pressure ionization, so the EOS should not be used above ρ ~ 10⁻² g cm⁻³.
