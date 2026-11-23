# Sub-shock convection in the Pop III accretion runs

Branch `popiii-subconv`, from `logmass` b482bdc. That is the commit zie1 runs since segment S6; S1–S5 ran on iedef aace2a0,
which is an ancestor. Both contain `popiii-combined`. Labels: [run] [test] [read] [estimate] [interpretation].

## M0 — Inventory (22 Nov 2026, no physics changes)

### M0.1 The MLT path in full

**Flux** [read] `julia/TitanResidFast.jl:1538–1591`. On interior interfaces k = ngrs+1…ngre:
- The Kippenhahn–Weigert flux is F = wτ ρ c_P T √(gδ) ℓ² ζ³/(4√2 H_P^{3/2}), with ζ = √(U²+∇−∇_ad) − U.
- The velocity v = ℓ√(gδ/8H_P) ζ and ℓ are stored in `W.vcv` and `W.lcv`.
- Each quantity is evaluated at the interface:
  - T and ρ are arithmetic means of the two zones.
  - P = p_gas + f E_r.
  - g = G m/r² (with the `xmu` geometry).
- Two interfaces are skipped: those with g ≤ 0, and those with |Δln P| < 1e-12.
- The flux enters the gas-energy row as rmu·fcv (line 1608).

**Switches and gates in the zie1 stack** [read; values from `runs/zie1/run.env` and `popiii_m2c.jl:910–913`]:

| item | value in zie1 | where | effect |
|---|---|---|---|
| `CONVECTION` | on (because HOT=1) | m2c:911 | MLT block active |
| α (`ALPHAMLT`) | 1.5 | m2c:911 | |
| `MLT_RCAP` | on | m2c:911 | ℓ = α(H_P⁻²+r⁻²)^{-½} |
| `MLT_HSDISC` | on | m2c:911 | hydrostatic Δln P taken from the discrete momentum row |
| `MLT_HS` | on (code default) | TRF:858 | ∇ is blended towards the hydrostatic value where Δln P/Δln P_hs < ¾ (see below) |
| `MLT_AD` | on (code default) | TRF:850 | ∇_ad, c_P, δ, κ are evaluated at the iterate, as Duals |
| `TAU_MLT` | 1 | TRF:843 | weight wτ: smooth step in ln τ_R, from 0 at τ = ½ to 1 at τ = 2 |
| `MLT_NOINFALL` | on (because COLDMLT > 0) | m2c:912, TRF:413 | wτ × (1 − smoothstep): 1 for inflow at 10 km/s, falling to 0 at 30 km/s |
| `COLDMLT` | 1e4 K | stack | table (eq) zones with T > 1e4 K get the table's ∇_ad |
| `HOTT1`/`HOTT2` | 3e5 / 9e5 K (defaults) | m2c:840 | hot branch: stellar EOS, real ∇_ad |
| network zones, T < T1 | ∇_ad = 1e30 | `PopIIIMaterial.jl:251–256` | **never convective** |

Origin of the 1e30 [read]: commit b8f0bd9, 1 Nov. With ∇_ad = 10, a T jump inside the protostar crossed it. The placeholder
c_P = 1e9 then gave 2e40 erg/s between two zones, and Δt collapsed at step 14800. Notes are in `docs/POPIII_NOTES.md` S6-C1c.
M2 must reproduce that run state.

**Not consistent with the Jacobian rule** [read]. These are inputs to M2/M3, not to be fixed in M0:
- The weight wτ and the `MLT_NOINFALL` factor are Float64 (`BUFS[:cvw]`), built in a pre-pass (TRF:383–418) from `base`, i.e.
  the iterate. They are held constant in the Jacobian: the τ(ρ, T) and u dependence is not differentiated.
- The `MLT_HS` blend weight `sw` is computed from `ForwardDiff.value` (TRF:1566), so it is constant in the Jacobian. It is
  also built on a `clamp` and a smooth step with thresholds ¼–¾.
- `mix!` (species mixing) is operator-split after the step, with D taken from the converged `step_diag` (see M0.3).

**Which zones are network, table or hot** [run, snapshots of zie1 S1 read with `scripts/subconv/branches.jl`]:
- Table zones are selected by n_H ≥ `PSWNH` = 1e19 cm⁻³ and T ≥ 4500 K (with hysteresis); `PSWITCH` = 1e30 is inactive.
- The post-shock layer has n_H ~ 1e15–1e18, so it is network gas.

| epoch (step) | M* | shock zone | N² < 0 band of the artefact | its zones | branch |
|---|---|---|---|---|---|
| 0 yr (115) | 0.017 | 27 (6.0 R☉) | 4–5 R☉ | 25 (3.6 R☉, 19 600 K) / 26 (4.6 R☉, 14 200 K) | table (MLT possible) / **network** |
| 5 yr (295) | 0.29 | 36 (60 R☉) | 32–53 R☉ | 33–35 (17 000–22 000 K, 12 300 K) | **network** |
| 111 yr (655) | 3.2 | 42 (278 R☉) | 147–245 R☉ | 40–41 (25 300 K, 6840 K) | **network** |
| 1.76 kyr (2395) | 34 | 42 | stable | 34–37 hot, 38 onwards network | |

Preliminary answer to M1(a) [read + run-snapshot]: for M* ≳ 0.1 M☉, MLT is off in the N² < 0 zones because of the thermodynamic
gate (network ∇_ad = 1e30). `COLDMLT` does not reach them. This still has to be confirmed by M1 with the gate weights at the
iterate.

### M0.2 Where the explorer data come from [read]

**Recording:**
- Writer: `julia/popiii_m2c.jl`. Settings at L366–377 (`RUNREC=<dir>`), `recopen!` at L388, `recsnap!` at L406.
- Cadence (L2068): steps 1 and 2, then every `SNAPEVERY`-th step. zie1 records every 5th step.
- Output: one serialised NamedTuple per snapshot, appended to `snaps_<k0>_<pid>.jls`, plus `meta_*.txt` holding the git hash
  and the env.
- Per-zone content: r, u, m (interfaces); ρ, T, e, pg, er, fr, chi, kp, fedd, qf, eq, Λext; species Y; hotY. Scalars: Lnuc,
  Llines, kph, Teff.
- MLT quantities (fcv, ∇, gates) are **not** recorded.

**Export:** `julia/runrec.jl` (`julia runrec.jl <dir> [every|budget=] [skip=f1,f2]`).
- Segments are stitched by lineage.
- `frame()` is at L108; the output is written to `<dir>/export/{frames.f32, frames.json, track.csv, nuc_v1.csv}`.

**Format:**
- `frames.f32` holds little-endian Float32 values, frame after frame. Within a frame the layout is field-major, shape
  (nfields, n), starting at byte offset `off`.
- `frames.json` = `{fields:[[key,desc]…], scalars:[…], nuc:{…}, frames:[{n, off, s:[…]}]}`.
- All quantities are cell-centred: interface quantities are averaged onto cells.
- Full field list (27): logr, logm, logrho, logT, logTrad, u, L, xHp, xH2, logxHm, logtau, Gamma, FcE, logkap, logkapP, Mdot,
  eq, logLam, logr_o, logm_o, logXH, logXD, logXHe3, logXHe4, logXC, logXN, logXO. The zie1 export skips logkap, logkapP, eq
  and logLam, leaving 23.
- Scalars (23): step, t_yr, dt_s, logn_c, T_c, M_star, R_star, r_ph, T_eff_lum, T_eff_sl, L_ph, Mdot_1e4, u_min, L_nuc,
  L_lines, nz, nchemfb, L_D, L_pp, L_He3, L_CNO, L_3a, L_ph_med.

**Reading:**
- `scripts/abnfig/abnfig.py` (L82–124) base64-decodes `frames.b64.txt`, which is `base64 -w0 frames.f32`.
- Missing fields become all-NaN arrays and their panels are skipped; missing scalars become NaN columns.

**Not in the repository:**
- The b64 step, the publishing, and the Explorer HTML itself, which live in an old scratchpad:
  `/tmp/claude-1002/-home-tabel-Projects/3ebf2758-…/scratchpad/` (`panel_all.sh`, `panel_manifest.py`, `pubsplit.sh`,
  `panel_index_v80.html`, …).
- M1 will copy these into `scripts/explorer/` before extending them, so that the viewer is versioned.

### M0.3 Species mixing in the hot core [read]

- `StarMaterial.mix!` (`StarMaterial.jl:520–624`) is called from `hotburn!` (`popiii_m2c.jl:883–888`) after `burn!`. It is
  **operator-split** after each converged step.
- It solves implicit conservative diffusion on the mass grid, ∂Y/∂t = ∂_m[(4πr²ρ)² D ∂_m Y], with D = v ℓ/3 taken from the
  converged step's `fcv`, `vcv` and `lcv`. The factorisation is tridiagonal, with T re-solved at fixed e.
- It carries patches: `bridge = 4` (gap filling), `tmax = 0.02` (partial mixing) and `capmix = false`.

**Could the network species use it?**
- Mechanically yes: the network species are per-zone Y as well.
- But for H₂, H⁺ and e⁻ in the sub-shock layer, the chemical relaxation time and the turnover time are comparable (M1 will
  check this), and the species are coupled implicitly in the Newton system (`COUPLED=1`, CK rows).
- An operator-split diffusion after the step would be a lagged input there.

[interpretation] The consistent path is a turbulent-diffusion flux term in the coupled species rows, D = α_D ℓ√e_t, at the
iterate (M3). Once that exists, the hot-core `mix!` should move to the same path and the `bridge`/`tmax` patches should go.

### M0.4 Checkpoints available for the M1 diagnosis [run]

zie1 keeps only `_auto0` and `_auto1` at segment ends (CKEVERY = 100, so they overwrite each other). Epochs from
`export/track.csv`, with t₀ = 242 378.26 yr (first M* > 0):

| checkpoint | step | t − t₀ | M* |
|---|---|---|---|
| S1 | 4000 | 1.80 kyr | 34.5 |
| S2 | 8000 | 1.84 kyr | 35.0 |
| S3 | 12000 | 2.70 kyr | 46.7 |
| S4 | 16000 | 2.79 kyr | 48.0 |
| S5 | 19333 | 25.7 kyr | 177 |
| S6 | latest | ~30.8 kyr | 193 |

The 0, 5 and 111 yr epochs (steps 115, 295, 655) lie **inside S1, and no checkpoint is kept there**. Options:
1. Replay S1 from the IC in a separate run directory, with `SAVE` checkpoints every 50 steps up to step ~2500, on aace2a0 (the
   commit of S1). That should be bit-identical to the production run. Check by comparing against the recorded snapshots.
   Estimate: ~1 h of one 8-thread job, from S1's 4000 steps in 1 h 50 m.
2. Rebuild the state from snapshots (`rebuild_from_profile`): not bit-identical, so not acceptable for the diagnosis.

zie3 is analogous (S1… segments). Its equal-mass epochs will be mapped from its own track in M1.

### M0.5 Test matrix (M2–M4)

| # | test | milestone | gate |
|---|---|---|---|
| T1 | `NKCHECK=2` on the post-shock zones of replay checkpoints at 5 yr, 111 yr and 1.76 kyr, and on zie1 S6 | M2, M3 | every new row and term passes |
| T2 | `popiii_regress.jl` (SOD 397, SUBC 587, SUPC 713) + `noh_test.jl` | M2, M3 | bit-identical when off; close when on, where no convection is expected |
| T3 | b8f0bd9 state (step 14800 of the Nov 1 run; locate its checkpoint) | M2 | no 1e40 flux, no Δt collapse |
| T4 | static unstable polytropic layer between stable layers | M3 | growth vs N; saturation = MLT flux; overshoot vs α_D |
| T5 | Ledoux-stable, Schwarzschild-unstable layer | M3 | the stated behaviour |
| T6 | advected unstable layer, t_settle·\|N\| = 0.1 to 100 | M3 | amplitude vs exp(∫σ dt) |
| T7 | supersonic inflow through an unstable profile | M2, M3 | no convection, with no gate |
| T8 | `ETOTMON` drift on T4–T7 and on a zie1 checkpoint | M3 | at the level of production |
| T9 | steady local limit: TDC flux = MLT flux (α mapping) | M3 | to round-off/tolerance |
| T10 | 111-yr layer at 100, 300 and ~1000 effective zones (GRIDEQ band around R*) | M3 | convergence of F_conv/L, depth, R* |
| T11 | B0–B5 from replay checkpoints, into phase C (≳ 30 M☉); B0, B2, B5 repeated at M3 resolution and on zie3's mesh | M4 | comparison table |

## Design of M2 (buoyancy and thermodynamics) and M3 (TDC), 22 Nov 2026

Defaults decided while Tom was away (he asked for autonomous work to a working model, 22 Nov):
- replay S1 for checkpoints;
- version the explorer scripts;
- turbulent species diffusion in the coupled species rows;
- base branch `logmass`.

### D1. Buoyancy from the density contrast (replaces ∇_ad = 1e30, `COLDMLT` and `MLT_HS`)

**Blob criterion.** Across interface k (zones k−1 → k, outwards), a blob that moves adiabatically has dln ρ_b = dln P/Γ₁,b. The
background has Δln ρ = α Δln P − δ Δln T + φ Δln μ, with α, δ, φ at fixed composition (frozen EOS of the zone gas). Define
- A_k = α Δln P − δ Δln T + φ Δln μ − Δln P/Γ₁,b,
- N² = −g A/Δr.

The blob is unstable if A > 0.

**Superadiabaticity used by MLT and TDC:**
- sup ≡ A/(δ_b |Δln P_hs|).
- Δln P_hs = −g M_k/(r_k² P̄) is the hydrostatic pressure difference of the discrete momentum row (`MLT_HSDISC`).
- [read/derivation] In hydrostatic layers, where Δln P = Δln P_hs and Γ₁,b is the frozen value, sup = ∇ − ∇_ad + (φ/δ)∇_μ.
  That is the Ledoux form of the production ∇ − ∇_ad, correct to first order in the differences.
- Behind the shock (Δln P → 0), sup = −δΔln T/(δ_b|Δln P_hs|) + …. This is the correct buoyancy of a flat-pressure layer, so
  no blending (MLT_HS) is needed. The hydrostatic factor multiplies sup instead of being subtracted from it, so deep convection
  (sup ~ 1e-6) keeps its accuracy.

**Blob thermodynamics.** The blob's thermodynamics interpolate between equilibrium chemistry (the table EOS: Γ₁,eq, δ_eq,
c_P,eq, latent heat included) and frozen chemistry (the zone's network/CK EOS at fixed x):
- 1/Γ₁,b = f/Γ₁,eq + (1−f)/Γ₁,fr; δ_b and c_P,b likewise linear in f.
- f = t_turn/(t_turn + t_chem) = ℓ/(ℓ + t_chem·v).
- t_chem = 1/max_i |∂S_i/∂x_i| of the coupled CK rates at the iterate (nested ForwardDiff, so it is differentiated).
- v is the MLT velocity (M2), or w (M3).
- Table zones and hot zones have f = 1 (there is only an equilibrium EOS).
- Network zones that are not coupled have t_chem = ∞, so f = 0 (their composition is held fixed during the hydro step).

### D2. The TDC equation (Kuhfuß type, w = √e_t, MESA-TDC-like Riccati form)

**Unknown.** w_k ≥ 0 [cm s⁻¹] on the interfaces (the staggering of u). It is Newton unknown `iw` = the last row, a linear variable
with the scale of u, and w = 0 at the centre and the outer boundary.

**KW quantities at the interface,** as now: a = ℓ√(gδ/(8H_P)), U = 3acT³/(c_P ρ² κ ℓ²)√(8H_P/(gδ)), H_P = P/(ρg),
ℓ = α(H_P⁻² + r⁻²)^(−½). Then:
- Ỹ ≡ sup − 2U w/a, the blob's ∇ − ∇′ with its radiative losses;
- Dw/Dt = (1/ℓ)(A_s a² Ỹ − A_d w²) + (1/ρ)∇·(α_t ρ ℓ w ∇w) + α_t ℓ |∇w|²;
- F_c = A_c ρ c_P T w (ℓ/2H_P) Ỹ, which is negative in overshoot regions;
- P_t = α_p ρ w² (zone-centred, ē = mean of the two interface w²);
- D_mix = α_D ℓ w.

**MLT limit.** With A_s = A_d = A_c = 1, the steady, local root of Dw/Dt = 0 is w² + 2Uaw = a² sup, i.e. w = aζ, the KW velocity.
F_c is then exactly the production MLT flux for the same sup and α [derivation; checked numerically in T9].

**Growth.** At w = 0 with sup > 0, w grows linearly. With constant coefficients and U = 0, w = w_s tanh(w_s t/ℓ),
w_s = a√sup = (ℓ/√8)N. The time to saturation is ~√8/N (T4, T6).

**Complementarity w ≥ 0.** Let F_k = (w_n − w_o) + Δt·[transport − source], in cm s⁻¹. The row is the Fischer–Burmeister function
w + F − √(w² + F²) = 0, which is semismooth Newton. Where sup < 0 and there is no turbulence nearby, F > 0 strictly, so w = 0
exactly. This is MESA TDC's "decay to zero in finite time", made implicit with transport.

**Diffusion.** The diffusion term is the e_t flux −α_t ρ ℓ w∇(w²) divided by 2w. It is non-singular because the Kuhfuß
diffusivity ∝ w, and it seeds overshoot into w = 0 zones.

**Advection.** The interface control volume (mass m_k = ½(ρΔV)_{k−1} + ½(ρΔV)_k) is advected with the same ṁ, face values and
θ-weighting as the momentum row. Inflow through a supersonic layer carries w = 0, and the transit time ≪ ℓ/w_s, so there is no
convection without any gate (T7).

**Energy, exactly conservative.** The realised local exchange per interface CV is
Q_k ≡ (m e)ⁿ − (m e)ᵒ over Δt + div_k(advective e flux) + div_k(diffusive e flux) + P_t dV_k, with e = w². The e fluxes are
discrete and telescoping. The gas energy row of zone j gets −½(Q_j + Q_{j+1}). P_t enters the momentum row as rmu_k(P_t,k − P_t,k−1),
paired with the P_t dV term like p_gas (compatible discretisation). Gas + E_t + kinetic + gravitational energy is then conserved
to round-off. `ETOTMON` adds E_t = Σ m_k w_k².

**Coupling:**
- `fcv` = F_c when TDC is on.
- The τ gate (`TAU_MLT`) and the inflow gate (`MLT_NOINFALL`) are not used under TDC. The KW U → ∞ as κρℓ → 0 suppresses thin
  layers physically (the Henyey thin-blob correction is an open item).
- Species: D_mix at interfaces enters the coupled species rows (8, 9) as −4πr²ρ D_mix Δx/Δr fluxes, at the iterate.
  - The other network species are operator-split in the chemistry step already, so they are mixed in that same step with the
    converged D_mix (as they get the converged hydro state).
  - The hot-core `mix!` moves to the same D_mix (later).
- Eddy viscosity: α_m ρ ℓ w as an additional viscous pressure, with its heating booked into e_t (default α_m = 0 in the first
  version).

**Switch.** `TDC=1` (provisional, default off, bit-identical when off; it adds the row, so neqn is unchanged when off). Listed in
the switch inventory, with the removal plan "replaces the MLT block, `MLT_HS`, `MLT_NOINFALL`, `COLDMLT` and the network
placeholder once M4 is accepted".

**Code layout.** The coefficients, row and flux terms go in a module `julia/TDC.jl` as pure functions on arrays. The residual and
the stand-alone test harness (T4–T7, with prescribed background profiles) call the same functions.

## M1 — Diagnostics (22 Nov 2026, `CONVDIAG=1`; solution bit-identical)

Labels: [run] [test] [read] [estimate] [interpretation].

### M1.1 What was added [read]

One switch, `CONVDIAG=1` (env, default off; in `docs/POPIII_SWITCH_INVENTORY.md` §13: diagnostic, permanent). With `RUNREC` set every snapshot gets an extra field
`cd` (a NamedTuple of Float64 vectors; old readers ignore it).

| piece | where |
|---|---|
| `conv_thermo(m::NetGas, k, ρ, T)`: ∇_ad, c_P, δ at equilibrium chemistry (table EOS, cold_thermo's formula) and at frozen chemistry (coupled zone: CK EOS `cpl_g` at fixed x₁, x₂ of the iterate; otherwise `ChemNet.eos` at `m.x[k]`; hot zone: `HOT.thermo` for both), φ = 1, μ, branch code (0 network, 1 coupled, 2 table, 3 hot, in `material_gas0`'s order). Generic in the number type; `material_thermo` is unchanged | `julia/PopIIIMaterial.jl:261` |
| cache save/restore around the diagnostic evaluation | `PopIIIMaterial.jl:293` |
| `TX.conv_diag_pass(dt)`: one Float64 residual evaluation at the accepted state, exactly the mechanism of `TX.step_diag` (which is how `StarMaterial.mix!` gets `fcv`, `vcv`, `lcv`), with `CVD_ON[]` set | `julia/TitanResidFast.jl:2080` |
| recording blocks (guarded by `cvp = CVD_ON[] && T === Float64`: never active in a solve): raw and blended ∇, HS weight, ∇_ad/c_P/δ/κ as used, gate-free flux and velocity, and the terms of the ECOMPAT gas-energy row | `TitanResidFast.jl:1557–1730` |
| τ and inflow factors of `cvw` in extra buffers (`BUFS[:cvwt]`, `:cvwi`, `:cvtau`; `cvw` itself computed as before) | `TitanResidFast.jl:408` |
| `convdiag_snap`: assembles `cd` (N², ∇_L, ∇_rad, timescales, t_chem, layer scalars) from the pass and `conv_thermo` | `julia/popiii_m2c.jl:425`, call in `recsnap!` |
| export: `cv_*` fields and scalars, only when some snapshot has `cd` | `julia/runrec.jl:82–110` |
| abnfig page 3 "convection" (only if `cv_*` fields exist), `--cvframe` | `scripts/abnfig/abnfig.py:904` |
| explorer v81: convection panel group + layer-history cards, hidden without `cv_*` fields (v80 untouched; scripts copied verbatim to `scripts/explorer/`) | `scripts/explorer/panel_index_v81.html` |

### M1.2 Definitions [read]

Interface vectors have length n+1 with index j = k−a+1 (interface k between zones k−1 and k), zone vectors length n.
- `nab`, `nab_raw`, `hs_w`: ∇ the MLT block used (after the `MLT_HS` blend), before it, and the blend weight sw (1 = untouched). `nad_used`: what `material_thermo` gave (1e30 = network zone).
- `nad_eq`, `nad_fr` (and `del_*`, `cp_fr`): `conv_thermo` at the interface mean (ρ, T) with the composition of zone k (the one MLT uses).
- `nab_rad` = 3κρF_rad H_P/(4acT⁴) ≡ 3κLP/(16πacGmT⁴) with L = 4πr²F_rad, F_rad the radiative flux of the interface actually present (`frn`), κ the Rosseland κ MLT uses.
- N² of a zone-centre difference, without dividing by Δln P: N²_S = g_eff δ_fr (dlnT/dr − ∇_ad,fr dlnP/dr), `N2_eq` the same with the equilibrium ∇_ad, δ, N²_L = N²_S − g_eff φ dlnμ/dr; g_eff = −ΔP/(ρΔr), P = p_gas + fE_r. `nab_mu` = Δlnμ/ΔlnP, `nab_L` = ∇_ad,fr + (φ/δ)∇_μ.
- Gates: `w_tau` (optical depth), `w_inf` (`MLT_NOINFALL`), `w_th` = 1 if `nad_used` < 1e29. Flux as used `fconv`; `Fconv_over_Ftot` with F_tot = F_rad + F_conv. `fpot_th` ("all gates open"): w_τ = w_inf = 1 and, where the thermodynamic gate is closed, the frozen ∇_ad, c_P, δ in its place, with the same ∇, ℓ, H_P, κ.
- Mach = |u|/c_s with c_s = √(5/3 p_gas/ρ); `t_settle` = H_P/|u|, `t_N` = 1/√|N²_L|, `t_turn` = ℓ/v (potential v where `fconv` = 0), `t_chem` = 1/max|∂S_i/∂x_i| of the two coupled species (NaN for non-coupled zones).
- `e_*`: stor, pdv, av, kin, adv, dfe, rad, conv, chem, src and their sum `e_res` for the gas-energy row of each zone [erg s⁻¹ sr⁻¹]. [test] The ten terms sum to `e_res` to 3e-16 of the largest term in every zone (IC step 2 and 40).
- **Sub-shock unstable layer** (`cd.layer` = [r_in/R*, r_out/R*, m_layer/M☉, max(−N²)t_dyn², max F_conv/F_tot, max F_pot/F_tot, n, R*/R☉]): the run of interfaces just inside the shock (R* = first interface from the centre with u < −2 km/s, as in `runrec.jl`), contiguous inwards, with Mach < 0.5 and N²_L < 0 or N²_eq < 0. Bounds = outermost and innermost interface of the run; N² is the smaller of the two; t_dyn² = r³/Gm.

### M1.3 Tests

[test] **Bit-identity**, `scripts/subconv/bitid.sh` + `compare_snaps.jl` (every field of every snapshot `isequal` bitwise, `cd` excluded; the Newton counters of the two runs are equal too):
- IC, zie1 settings, 40 steps, SNAPEVERY=5: 10 snapshots, 300 fields, 0 differences. `docs/logs/subconv/bitid_ic40.txt`.
- Restart from the replay checkpoint `rp_k650` (zie1rp, 111-yr epoch, M* = 3.2 M☉), 20 steps: 4 snapshots, 120 fields, 0 differences. `docs/logs/subconv/bitid_rp650.txt`.
The material caches (NLOC, CPLCACHE, …) are saved and restored around the extra evaluation; whether that is needed was not tested separately.

[run] **Zones around the shock, step 655 and 670** (`docs/logs/subconv/shocktable_rp650.txt`; 100 zones, R* = 278 R☉, shock at interface 42):
- Interfaces 40 and 41 (r/R* = 0.60 and 0.77, T = 37 000 and 16 000–17 000 K, Mach ≤ 0.02) have N²_eq < 0 and N²_L < 0, so the instability is not an artefact of equilibrium versus frozen ∇_ad. They are coupled network zones (branch 1), `w_th` = 0, with `w_tau` = `w_inf` = 1: only the thermodynamic gate (∇_ad = 1e30) closes MLT. ∇ exceeds ∇_ad,fr by 0.04 and 0.14 (step 655).
- F_conv = 0 there. With all gates open the MLT flux would be 58 and 124 F_tot (step 655), 58 and 376 F_tot (step 670), [estimate: same ∇ and ℓ, frozen c_P and δ, no feedback]. Layer: 2 interfaces, 0.38 M☉, max(−N²)t_dyn² = 7.4 (step 655).
- This confirms the M0 reading [run] that for M* ≳ 0.1 M☉ the gate, not τ or infall, switches MLT off in the layer, and that `COLDMLT` (table zones only) does not reach it.
- Not exercised: hot-branch zones (T > 3e5 K is not reached at these epochs), so branch 3 of `conv_thermo` is untested [read only]; `conv_thermo` has not been called under an outer ForwardDiff.

## M2 — Buoyancy and thermodynamics of the network gas (22 Nov 2026, `MLTBUOY=1`, provisional, default off)

Labels: [run] [test] [read] [estimate] [interpretation]. Logs in `docs/logs/subconv/m2/`; run scripts `scripts/subconv/m2run.sh`, `m2snap.jl`.

### M2.1 What was added [read]

| piece | where |
|---|---|
| `TDC.jl`: `gamma1_inv/gamma1` (1/Γ₁ = α − δ∇_ad), `blob_f` (f = ℓ/(ℓ + t_chem a)), `blob_sup` (design D1; derivation in its docstring). Pure, number-type generic. Included in `TitanResidFast` | `julia/TDC.jl:18,27,55`; `TitanResidFast.jl` top |
| `conv_thermo` returns also `alpha_eq/fr`, `gam1_eq/fr`; keywords `x` (the zone's iterate x_H⁺, x_H₂), `tablefr` (table zone: frozen := equilibrium), `fd` (EOS derivatives by central differences) | `PopIIIMaterial.jl:266` |
| `conv_tchem(m, k, ρ, T)`: 1/max_i\|∂S_i/∂x_i\| from `cpl_g` rows 3, 4; 0 for table and hot zones, Inf for non-coupled network zones; does not touch `CPLCACHE` | `PopIIIMaterial.jl:314` |
| `MLT_BUOY` switch, per-zone thermodynamics (`_buoy_zone`), `buoy_blob`, and the interface branch of the MLT block (before the old ∇ computation, so with the switch off the old path is untouched) | `TitanResidFast.jl:884, 945, 962`, branch in the MLT loop (~L1640) |
| `MLTBUOY` env, exit summary | `popiii_m2c.jl:1022` |
| `CONVDIAG` records `sup`, `fchem`, `tchem`, `delb`, `cpb` (`cd.sup…`; exported `cv_sup`, `cv_fchem`, `cv_tchem`) | `TitanResidFast.jl` `CVD_NAMES`, `popiii_m2c.jl convdiag_snap`, `runrec.jl` |

### M2.2 Choices [interpretation, with the reasons]

1. **Zone values, interface means.** α, δ, c_P, Γ₁ (equilibrium and frozen) are evaluated in each of the two zones at the zone's own (ρ, T, x_H⁺, x_H₂) of the Newton iterate, then averaged. So the Jacobian has the species and the state of both zones in every coefficient (no `ForwardDiff.value`).
2. **Δlnμ.** The `Δlnμ` argument of `blob_sup` is the remainder of the actual density contrast, Δlnρ − α_fr ΔlnP + δ_fr ΔlnT (φ_fr = 1). Then A = Δlnρ − ΔlnP/Γ₁,b exactly; α_fr and δ_fr cancel in A. Reasons: (i) the linearisation α ΔlnP − δ ΔlnT has second-order errors of the size of the signal where ΔlnT ~ 1 (ΔlnT = −0.7 and −1.3 at interfaces 40, 41); (ii) a μ from the zones' EOS would count the fast species twice (they follow ρ, T in the CK expansion and are already in α, δ). In a hydrostatic layer it is the same sup = ∇ − ∇_ad + (φ/δ)∇_μ [derivation, `blob_sup` docstring].
3. **f.** Per zone f_j = ℓ/(ℓ + t_chem,j a), a = ℓ√(gδ_fr/(8H_P)); interface f = mean of the two. It is smooth, finite where sup ≤ 0, and needs no velocity iteration (v = aζ would need a fixed point). a ≥ v for sup ≤ 1, so f is a lower bound of the f with the true turnover time. Table zones: t_chem = 0, f = 1. Non-coupled network zones: f = 0 (the equilibrium arguments are then not touched: NaN-safe).
4. **Derivatives of the EOS by central differences, ln step 1e-3 (`THERMO_FD_H`).** `conv_thermo` with inner ForwardDiff under an outer ForwardDiff gives finite values but NaN partials through the CK expansion of `cpl_g` (it supports one level of Duals) [run: first MLTBUOY run, all Newton attempts failed with "non-finite J v"]. The stencil is differentiated exactly by the outer Duals, so the Jacobian is that of the discrete thermodynamics. Not yet quantified against the exact derivative on a state (the M1 `fd=false` form is still the diagnostic default) [open].
5. **Dual thermodynamics only where sup > 0.** First every needed zone is evaluated in Float64 and the interface sup from its value; zones of interfaces with value(sup) > 0 are evaluated again with Duals. An interface with sup ≤ 0 has no flux and no derivative (the old code `continue`s at the same place), so this is exact, not a gate.
6. **Opacity.** κ of the flux is the interface χ of the radiative transport, `avchi[k]` (harmonic mean of the zones' χ at the iterate). `material_opac` at the interface-mean state returns NaN at some iterates (NLOIT opacity of the zone, found at jstep 700 of the k650 restart); the old code evaluated it too, but never used it where ∇_ad = 1e30.
7. **Kept as in production:** τ weight, `MLT_NOINFALL`, ℓ (MLT_RCAP), g, H_P, `MLT_HSDISC` form of ΔlnP_hs (taken from the discrete momentum row whatever `MLT_HS` is). Not applied: `MLT_HS` blend, the |ΔlnP| < 1e-12 skip, the 1e30 placeholder. Where the table cannot be evaluated (`EQ.eos` throws or NaN) f := 0 and the interface counts in `BUOYEQFAIL` (printed at exit).

### M2.3 Tests

[test] **Bit-identity, MLTBUOY unset** (`bitid_ic40.txt`): IC, zie1 settings, 40 steps, SNAPEVERY=5, HEAD bbb60b9 against the reference c523e91 (the commit before M2; worktree built and removed): 10 snapshots, 300 fields, 0 differences. Same TitanCore (the branch's) in both.

[test] **NKCHECK=1** (finite-difference J v against the band J, restart rp_k650, dt 1e5 s, `NWPROBE=1.0e5 NWPROBEAT=0`, `NKZ=36:42`; iteration 1, because iteration 2 is not reached: the probe converges). `nkcheck_rp650_base.txt` (MLTBUOY unset) and `nkcheck_rp650_buoy.txt`:

| equation | max rel, all zones: base / MLTBUOY | zones 36–42: base / MLTBUOY |
|---|---|---|
| r, m, ρ, u | 5.0e-6, 3.9e-6, 4.7e-6, 4.2e-6 / same | 4.9e-6, 3.1e-6, 4.5e-6, 1.8e-6 / same |
| **T (gas energy)** | 1.706e-2 (zone 1) / 1.708e-2 (zone 1) | **3.380e-6 (zone 38) / 3.380e-6 (zone 38)** |
| E_r, F_r | 9.6e-4, 17.6 (zone 43) / same | 9.6e-4, 1.2e-4 / same |
| x_H⁺, x_H₂ | 2.7e-4, 3.8e-6 / same | 8.5e-7, 2.2e-6 / same |

The gas-energy rows of the post-shock zones agree with the FD product at the noise floor of the product (ε = 1e-5), with and without the new flux; the large F_r rows (zone 43, 71–78) and T zone 1 are the same in both (closure rows, known in docs/POPIII_NEWTON_KRYLOV.md). The flux was active in that state: interfaces 40, 41 have sup = 0.147, 0.356, F_conv = 1.9e11, 5e10 erg s⁻¹ cm⁻² (first Newton iteration, same state) [run, debug print of an earlier build; not in a committed log]. All new rows go through the Duals of the iterate. The first version (inner ForwardDiff through `cpl_g`) was non-finite here; this is the version that passes [run]. The T row in zones 40–41 is not listed separately by `check!` [limitation].

[run] **Pilot B1** (`pilots_summary.txt`; TitanCore.jl of aace2a0 for BOTH baseline and MLTBUOY=1, as the restarts are bit-identical to the replay only with it; commit HEAD code otherwise; CONVDIAG=1, SNAPEVERY=5, zie1 settings; 8 threads):

| restart, steps | Newton-failure rejections base / buoy | smallest accepted dt base / buoy [s] | F_conv/F_tot in the layer base / buoy | sup in the layer buoy | max(−N²)t_dyn² (layer) base → buoy |
|---|---|---|---|---|---|
| k650 (3.2 M☉), 200 | 33 / 30 | 1.07e5 / 1.53e5 | 0 / 0.999 (step 655–680), later 0.09–0.98 | 0.38 → 0.045 in 25 steps, then 0.08, 0.17, 0.005 | 7.4 → 0.7 (680); 11–12 later (R* jumps 278 → 359 R☉, the layer moves to interface 41–42), 1.9 at 835 |
| k290 (0.29 M☉), 100 | 13 / 18 | 1.3e4 / 2.2e4 | 0 / 1 (interfaces 35–36) | 0.14 | 1.6 → 0.4 / 0.25 → 0.4 (not reduced) |
| k115 (0.017 M☉), 100 | 12 / 13 | 9.8e4 / 6.9e4 | 0–1 (table zones) / 1 (interface 2, 28–30) | 0.11–0.23 | 0.01–1.2 / 0.07–0.47 |

[interpretation] In the 111-yr layer the flux switches on at once (F_conv/F_tot 0.999 from step 655) and drives sup from 0.38 to 0.045 within 25 steps (ca. 3×10⁷ s), with N²t_dyn² from 7 to 0.7: the layer is driven towards, but not onto, the adiabat while the flux carries the whole luminosity. After the shock moves outward (R* 278 → 359 R☉, steps 690–700) the layer is rebuilt by accretion and the cycle repeats. No Δt collapse, no runaway flux (F_conv/F_tot ≤ 1) with the corrected opacity, and Newton failures are as in the baseline (30 vs 33). At k290 the flux saturates (F_conv/F_tot = 1) at sup = 0.14: the layer stays superadiabatic, so N² is not reduced there. These are 100–200 steps; no claim on the later phases (that is M4).
[run] The first pilot (κ from `material_opac` at the interface mean) collapsed: Newton failures at every Δt from step 700 on, Δt 5e-12 s at step 733, because κ = NaN at interface 42 (see choice 6) [run, not kept as a log].
Cost [run]: thermodynamics 100 s (632 093 Float64 zone evaluations, 0.16 ms) plus 59 s (26 319 with derivatives, 2.3 ms) in the 200 steps of k650; wall time of the pilots 1.1–1.6× the baseline. [read] The Float64 evaluations are made at every needed zone of every residual evaluation, including the Krylov products; the table range is not covered at some interfaces (600 and 300 counts in k650, k290/k115 runs), where the blob is frozen.

[run] **The 1e40 erg/s divergence of b8f0bd9** (T3): no checkpoint of the Nov 1 run (step ~14 700–14 850) exists. The checkpoints under `lib/popiii/ck` are `m2c_z19.aux` (step 9742, Oct 30), `m2c_z19_auto0/1` (≤ 9500) and `m2c_z19c3` (33 200, Nov 3); `runs/*/ck` start on Nov 11. Not run. Open: T3 needs a replay of that run on its commit.

[read] **T2.** `popiii_regress.jl` only reads stored snapshots (lib/p2_runs) and does not run the solver; `noh_test.jl` uses Z0Gas with `CONVECTION` off and cannot be started here (TitanFortran library not loadable in this environment). Neither exercises the changed path (the MLT block is skipped when `c.conv` is false); not run.

### M2.4 Open

- The EOS derivatives by central differences (choice 4) are not compared with the exact derivative on a state.
- The `NKCHECK` iteration was 1 (state converged at iteration 1 at dt 1e5 s); the T-row of the interfaces 40–41 alone is not extracted.
- The cost of the Float64 screening (every residual evaluation) is a candidate for M3 (the TDC row replaces the MLT block); the Dual thermodynamics are only paid where sup > 0.
- `BUOYT`/`BUOYEQFAIL` are global counters (not thread-safe; the MLT loop is serial).
