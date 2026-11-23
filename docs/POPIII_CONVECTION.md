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

## M3a — The TDC rows as pure functions and a stand-alone harness (22 Nov 2026, no solver change)

Files: `julia/TDC.jl` (M3a block after `blob_sup`; M2 functions unchanged), `julia/tdc_tests/{harness,t4_layer,t5_ledoux,t6_transit,t7_supersonic,t8_energy,t9_mlt,mesa_compare}.jl`,
logs and CSV tables in `docs/logs/subconv/m3a/`. Run: `julia -t 4 --project=<jlenv> julia/tdc_tests/<script>.jl` (each prints PASS/FAIL, writes its log, exit code = number of failures). TitanResidFast.jl, popiii_m2c.jl and PopIIIMaterial.jl are untouched.

### M3a.1 What was added [read]

- `TDCPar(; A_s, A_d, A_c, α_p, α_t, α_D, α_m)`; `fb(w,F)` (Fischer–Burmeister, cancellation-free for F ≥ 0, subgradient (1,1) at exactly (0,0)).
- `tdc_ay(ℓ,g,δ,H_P,c_P,ρ,κ,T) → (a,U)`; `tdc_local(w,a,U,sup,ℓ,par)` = (A_s a²(sup − 2Uw/a) − A_d w²)/ℓ; `tdc_flux(...)` = A_c ρ c_P T w (ℓ/2H_P) Ỹ; `tdc_pturb`, `tdc_dmix`; references `tdc_mlt_zeta`, `tdc_mlt_w`, `kw_mlt_flux`.
- Grid objects: `TDCGeom(nz, r, Az, mzn, mzo, dmdt, dVdt)` and `TDCBack(a, U, sup, ell, rhoz, ellz)`. Indices as in TITAN: zone j between interfaces j, j+1; interface k between zones k−1, k; CV mass m_k = ½(mz[k−1] + mz[k]) with mz = ρ·dvol (= the momentum row's `mk`); `dmdt` at the interfaces with TITAN's sign (positive = inward); zone-centre flux Φ_j = ½(dmdt[j] + dmdt[j+1]).
- `tdc_F!` (the F_k of D2, cm/s), `tdc_rows!` (φ_k = fb(w_k, F_k); φ at k = 1, nz+1 is w: Dirichlet w = 0), `tdc_exchange` (Q_k for all nz+1 CVs), `tdc_gas_exchange` (gas zone sources), `tdc_budget_check`.
- Arrays may carry ForwardDiff Duals (no `Float64` annotations in the functions; checked: the dense `ForwardDiff.jacobian` of the rows equals the 3-colour tridiagonal Jacobian of the harness to 0.0, T8).

**Discretisation of F_k (mirrors the momentum row)** [read]:
F_k = [(m_n w_k − m_o wo_k) − Δt(Φ_hi wb_hi − Φ_lo wb_lo) − w_k((m_n − m_o) − Δt(Φ_hi − Φ_lo))]/m_n − Δt·diff_k − Δt·R_loc(w_k).
The bracket is the momentum row's `((mu)ⁿ − (mu)ᵒ)/Δt − ½((ṁ_k + ṁ_{k+1})ub_k − …)` with u → w, with w_k times the discrete continuity residual subtracted so that it is the *advective* form (a uniform w is exactly invariant even when continuity is not closed). Face value wb at a zone centre: first-order donor (Φ ≥ 0: the upper interface) value, optionally blended with the centred value by `cen` (0 = upwind, the default; centred is unlimited; no limiter was implemented — the θ-weighting of the momentum row is `cen`). Diffusion: diff_k = (A_k G_k − A_{k−1} G_{k−1})/m_n + α_t ℓ_k (g_{k−1}² + g_k²)/2 with G_j = α_t ρ_j ℓ_j w̄_j g_j, g_j = (w_{j+1} − w_j)/Δr_j, w̄_j = ½(w_j + w_{j+1}), A_j the zone-centre area (4π r_c² for μ = 2, 1 for μ = 0): (1/ρ)∇·(α_t ρ ℓ w ∇w) + α_t ℓ|∇w|² as in D2 (non-singular at w = 0).

**Exchange Q_k** [read]: Q_k = (m_n e_k − m_o eo_k)/Δt + (Fe_hi − Fe_lo) + (Fd_hi − Fd_lo) + P_dV,k, e = w², Fe_j = −Φ_j wb_j² (the same face value as the row), Fd_j = −α_t A_j ρ_j ℓ_j w̄_j (e_{j+1} − e_j)/Δr_j (outward positive; zero outside the domain), P_dV,k = ½(Pt_{k−1} dVdt_{k−1} + Pt_k dVdt_k), Pt_j = α_p ρ_j ½(w_j² + w_{j+1}²), `dVdt` = Lagrangian volume rate of each zone supplied by the caller (zeros on a static grid). The two boundary CVs (k = 1, nz+1; half a zone, e = 0) are included, so Σ_k Q_k Δt = ΔΣ m e + Δt Σ_j Pt_j dVdt_j with *no* boundary term (the e that diffuses/advects into a Dirichlet boundary CV is booked as a sink Q < 0 there, i.e. it heats the gas). `tdc_gas_exchange` gives zone j −½(Q_j + Q_{j+1}) and the boundary CVs entirely to their only zone, Σ gas = −Σ Q. Q is the *realised* production − dissipation of the row (including its truncation error); the energy-conserving gas source is therefore consistent by construction.

### M3a.2 Defaults and choices [interpretation]

- A_s = A_d = A_c = 1 (steady local root = KW/MLT exactly, T9). α_D = 1, α_m = 0 (not used by the rows).
- α_p = 2/3 as requested (RSP `RSP_alfap` standard, MESA `mlt_Pturb_factor` = 2). **Finding [derived]:** our w = a ζ is the MLT velocity v_c = Kuhfuß's d v = (⟨v'²⟩/3)^{1/2} = (2/3)^{1/2} ϖ^{1/2}, not MESA's w_J = (3/2)^{1/2} v_c (w_J² = e_t). Kuhfuß's isotropic p_t = (2/3)ρϖ = ρ w² therefore corresponds to **α_p = 1** in our variable (= `RSP_alfap` 1 = MESA `mlt_Pturb_factor` 3 = `TDC_alpha_Pt` 1). 2/3 is the (2/3)^2-reduced value. The integration test should run α_p = 0, 2/3 and 1.
- α_t = 0.01 (RSP suggested starting value for `RSP_alfat`; Paxton 2019). In our variable this is `RSP_alfat` = 0.0082 (α_t^ours = √(3/2) RSP_alfat, below). The tests scan 0 … 3.

### M3a.3 Coefficient mapping [derived from docs/popiii_conv_lit.md; no other source used]

Variable map: w_J = (3/2)^{1/2} w, e_t = w_J² = (3/2) w² (RSP `e_t` and MESA e_t; Kuhfuß ϖ). Our local equation in MESA's form is 2Dw_J/Dt = ξ0 + ξ1 w_J + ξ2 w_J² with ξ0 = 2(3/2)^{1/2}A_s a² sup/ℓ, ξ1 = −4A_s U a/ℓ, ξ2 = −2(2/3)^{1/2}A_d/ℓ, to be compared with MESA's ξ2 = −(8√2/(3√3))α_D/(α h) = −2.177 α_D/ℓ.

| titan.jl | meaning | RSP (S&M08 / Paxton19): name = base value × control | MESA TDC (Jermyn23 / controls.defaults) | Kuhfuß 1986 |
|---|---|---|---|---|
| A_s | source S = A_s a² sup/ℓ | A_s = (4/3) `RSP_alfas` (α_s = ½√(2/3)·alfas, default alfas 1) [derived] | A_s = (4/3) `TDC_alpha_S` (α_S 1; the paper fixes α_s = 1/√6) | C = √(2/3)α_ML (43) |
| A_d | dissipation −A_d w²/ℓ | A_d = (4/3) `RSP_alfad` (α_d = (8/3)√(2/3) alfad) | A_d = (4/3) `TDC_alpha_D` (α_D 1) | c_D (24); c_D = α_d of RSP (same definition ϖ^{3/2}/Λ); the MLT-matching c_D = 32/(9C) or 128/(9C) is tied to his α_ML convention [not reconciled] |
| A_s/A_d | steady w² = a² (A_s/A_d) sup | alfas/alfad | 1/`TDC_alpha_D`·`TDC_alpha_S` (MESA (31): v_c² ∝ 1/α_D) | — |
| A_c | F_c = A_c ρ c_P T w (ℓ/2H_P)Ỹ | A_c = `RSP_alfac` (α_c = ½√(2/3)alfac; F_c = α α_c ρ T c_P e_t^{1/2}Y) [derived: ½√(2/3)·√(3/2) = ½] | `TDC_alpha_C` (L_conv = 4πr²αρc_PT w_J Y/√6) | enthalpy flux (41b) |
| radiative loss: Ỹ = sup − 2Uw/a | −2A_s U a w/ℓ in Dw/Dt | D_r (13): γ_r² = 12 A_s [derived]; `RSP_gammar` default 0, W&F standard γ_r = 2√3 ⇒ A_s = 1 | α_r = A_s [derived] (`TDC_alpha_R` default 0, α_r = 1 ⇔ γ_r = 2√3) | none in 1986 (W&F 1998) |
| α_p | P_t = α_p ρ w² | `RSP_alfap` [derived] (P_t = (2/3) alfap ρ e_t = alfap ρ w²); RSP default 0 | `mlt_Pturb_factor` = 3 α_p (default 0); `TDC_alpha_Pt` (in the e_t equation, default 0) | p_t = (2/3)ρϖ ⇒ α_p = 1 (26) |
| α_t | e flux −α_t ρ ℓ w ∇(w²) | α_t = √(3/2) `RSP_alfat` [derived] (F_t = −α_t ρΛ e_t^{1/2} ∇e_t; default alfat 0; 0.01 suggested) | L_t (3) defined, **set to 0** in MESA | j_t (36) |
| α_D | D_mix = α_D ℓ w | — | — | (34),(35): D = α_i Λ ϖ^{1/2}, α_i = √(2/3)α_MLT ⇒ α_D = √(3/2) α_i = α_MLT [derived, assuming Λ = ℓ] |
| α_m | η = α_m ρ ℓ w (not in the rows) | μ_Q = α_m ρ Λ e_t^{1/2} ⇒ α_m = √(3/2) `RSP_alfam` [derived] (default 0.25) | `TDC_alpha_M` (default 0) | μ = α_μ ϖ^{1/2}Λ (31) |
| P_t dV work | in Q (exchange) | + P_t dV/dt in (3) | `TDC_alpha_Pt` | −p_t div v (25d) |
| ℓ | α_ML min(H_P, r)-type cap (`MLT_RCAP`) | Λ = α H_P (`RSP_alfa` 1.5) | α h, h = P/(ρg) (`mixing_length_alpha` 2) | Λ |

Equations side by side (e = w², as in the sources: D(e_t)/Dt):

- Ours: Dw/Dt = (1/ℓ)(A_s a² Ỹ − A_d w²) + (1/ρ)∇·(α_t ρ ℓ w ∇w) + α_t ℓ|∇w|², Ỹ = sup − 2Uw/a, F_c = A_c ρ c_P T w(ℓ/2H_P)Ỹ.
- RSP (S&M08 eqs 3, 11–18): De_t/Dt + P_t DV/Dt = −ρ⁻¹r⁻²∂(r²F_t)/∂r + E_q + C, C = S − D − D_r, S = α α_s e_t^{1/2}(TPQ/H_P)Y/… (14–15), D = α_d e_t^{3/2}/Λ (12), D_r (13), F_t = −α_t ρΛ e_t^{1/2}∇e_t (18), F_c (17).
- MESA TDC (Jermyn23 eqs 2, 4, 8): De_t/Dt + α_Pt P_t Dρ⁻¹/Dt = q + C − ∂L_t/∂m (L_t not used), C = α e_t^{1/2}(TPQ/(h√6))Y − α_D(8√2/(3√3))e_t^{3/2}/(αh) − 48σα_r T³e_t/(α²ρ²c_Pκh²); local 2Dw/Dt = ξ0 + ξ1 w + ξ2 w² with coefficients frozen at the end of step; exact tanh/tan step (12)/(13). Y → YΓ/(1+Γ) (Cox) in MESA only; not in ours.
- Kuhfuß 25d: dϖ/dt = (source) − c_D ϖ^{3/2}/Λ − ρ⁻¹ div j_t − ρ⁻¹p_t div v + viscous, j_t = −α_t Λ ϖ^{1/2} ρ ∇ϖ.

### M3a.4 Tests [test; logs in docs/logs/subconv/m3a]

All eight scripts pass (units ℓ = a = 1, g = 8, δ = H_P = ρ = c_P = T = 1 in the layer tests: N² = gδ sup/H_P = 8a² sup/ℓ², w_s = a√sup = ℓN/√8).

- **T9 MLT limit.** The local Newton root of R_loc = 0 equals aζ to 4.1e-16 (525 cases, U 1e-3…1e3, sup 1e-8…1e2). F_c(aζ) = KW flux ρc_PT√(gδ)ℓ²ζ³/(4√2 H_P^{3/2}): 3.3e-13 relative (104 physical cases with U²/sup ≤ 1e3), and 6.9e-16 over all 240 cases with the steady-root form Ỹ = (A_d/A_s)(w/a)². The direct form of Ỹ = sup − 2Uw/a cancels catastrophically when U² ≫ sup (Ỹ = ζ² ≈ (sup/2U)²): relative error ≈ 4ε·4U²/sup, observed at 8 % of that bound. The harness (12 interfaces, no advection/diffusion) reaches aζ to 2.7e-13 (Newton tolerance 1e-12) for U 1e-3…1e3, from seeds 1e-3 w and 4 w.
- **T4 (i) growth.** dw/dt(0) = A_s a² sup/ℓ to 1e-8 (BE, dt 1e-6). BE vs w_s tanh(w_s t/ℓ): max error 1.8e-2 … 1.1e-3 for dt 0.1 … 0.00625 (ratio 1.97–1.996: first order). N·t90 = atanh(0.9)√8 = 4.164 (numeric 4.1664, BE dt 1e-3), independent of sup (1e-3 … 10), i.e. t90 = 4.16/N = 1.47 ℓ/w_s. **(ii) saturation.** Centre w/w_MLT = 1.00000000; max |w/w_MLT − 1| and |F_c/F_MLT − 1| on 23 < x < 27 (layer 20–30, α_t = 0.1): 4e-7 (U = 0), 7e-10 (U = 0.7). **(iii) overshoot.** The steady stable-layer profile has a first integral [derivation]: w'² = a²|s|/(2α_tℓ²) + w²/(3α_tℓ²), so w reaches 0 at finite distance Δ0 = ℓ√(3α_t) asinh(w_e/q), q = a√(3|s|/2), and falls to w_e/e after Δe = ℓ√(3α_t)[asinh(w_e/q) − asinh(w_e/(eq))]; for w_e ≪ q this is the balance of diffusion with the buoyancy damping a²|s|/ℓ, Δ ≈ ℓ√(2α_t) w_e/(a√|s|), for w_e ≫ q it saturates at ℓ√(3α_t) (exponential tail with ℓ√(3α_t) once w ≳ q). Numerics (1000 cells, h = 0.05) vs this: Δe/ℓ = 0.419/0.418 (α_t 0.1, |s| 0.1), 0.719/0.719 (0.3), 1.306/1.305 (1), 2.255/2.255 (3); |s| = 1: 0.184/0.184, 0.309/0.309, 0.551/0.551, 0.944/0.944; |s| = 10: 0.050/0.056, 0.078/0.082, 0.139/0.134, 0.220/0.220 (the 1/e distance is 1–4 cells there; worst deviation 11 %). Fitted exponents: Δe ∝ α_t^{0.44–0.50} (analytic 1/2), Δe ∝ |s|^{−0.46 … −0.49} at α_t = 0.1, 1. w_e/w_s at the unstable edge falls from 0.78 to 0.45 with growing α_t or |s|. α_t = 0: no penetration. Table: `t4_overshoot.csv`. **(iv) spherical.** The μ = 2 rows at R = 1000ℓ reproduce the plane run to 2.4e-4 (curvature ℓ/R).
- **T5 Ledoux** (ΔlnP = −0.1, ∇ = 0.45 > ∇_ad = 0.4 of the frozen gas, Δlnμ = −0.01, φ = 1; equilibrium Γ1,eq = 1/0.9, δ_eq = 1.5): `blob_sup` gives sup(f = 0) = −0.0500, sup(f = 1) = +0.1667 (Schwarzschild ∇ − ∇_ad = +0.05). f = 0 from w = 0: w stays exactly 0 (200 steps). f = 0 seeded w0 = 0.6: decays to 0 in finite time, analytic t* = (ℓ/w_c)atan(w0/w_c) = 5.429 ℓ/a (w_c = a√|sup|); BE reaches 0 at 6.00, 5.70, 5.60, 5.48, 5.45, 5.44, 5.432 for dt 0.2 … 0.001 (error ≈ 2–3 dt, first order), stays exactly 0, and follows the tan solution to 1.3e-4 (dt 1e-3). f = 1: converges to w_s = 0.408 (final w/w_s = 0.99999975; BE vs tanh 1.5e-3 at dt 0.02). There is no semiconvection model: the layer is convective or not solely through f (sup sign).
- **T6 transit** (L = 10ℓ unstable layer, uniform ṁ = ρu inward, w = 0 inflow, upwind): relative deviation of the downstream-edge amplitude from w_s tanh(w_s d/(uℓ)) (d = L − h, the BE-path distance) for 10 | 30 | 100 cells: N L/u = 0.1: 5.8e-5 | 2.0e-5 | 6.2e-6; 0.3: 5.2e-4 | 1.8e-4 | 5.5e-5; 1: 5.3e-3 | 1.8e-3 | 5.6e-4; 3: 2.5e-2 | 8.3e-3 | 2.5e-3; 10: 8.6e-3 | 1.6e-3 | 3.8e-4; 30: 4.3e-5 | 2.7e-7 | 7.9e-9; 100: 7.4e-9 | 5.5e-14 | 1.8e-14. Max 2.5e-2 | 8.3e-3 | 2.5e-3: first-order convergence (order 1.00). The mirrored problem (outward flow, mirrored layer) gives the same amplitude to 1e-15 (sign convention of dmdt checked). The e budget closes in every step (1.6e-16).
- **T7 supersonic** (sup = 1 everywhere, 100 cells, w = 0 inflow, no gate): max w/w_s = 0.0987 (u = 10 w_s, L = ℓ), 0.288 (L = 3ℓ), 0.0330/0.0987 (u = 30), 0.0099/0.0297 (u = 100), against the transit estimate tanh(L/(uℓ)) = 0.0997, 0.291, 0.0333, 0.0997, 0.0100, 0.0300 (always just below it); F_c/F_MLT equals w/w_s (U = 0). α_t = 0.1 changes max w by ≤ 1e-3 relative. Newton converges in 6 iterations at CFL 5 and at CFL = nz (no kinks).
- **T8 energy.** Random non-uniform grids (μ = 0, 2; cen 0 … 1; mass fluxes of both signs; deliberately violated continuity mzn ≠ mzo; random dVdt; α_p = 2/3, α_t = 0.3): |Σ_k Q_kΔt − ΔΣme − ΔtΣPt dVdt|/scale ≤ 5.1e-17 and |Σ gas + ΣQ|/Σ|Q| ≤ 5.9e-17. In T4 (transient, 0.5 steps) 1.1e-16, in T6 1.6e-16. A steady uniform MLT state has Q = 0 exactly; a growing w has Q > 0.
- **MESA pointwise ODE.** Our Riccati equation mapped to ξ0, ξ1, ξ2 as above: Jermyn's (12) (tanh, with the coth branch for w0 above the plateau; their paper does not state it) evaluated literally agrees with RK4 to ≤ 3.5e-14 (4 parameter sets: U 0 … 0.5, A_s = A_d = 1 and 4/3, α_D = 2 ⇒ A_d = 8/3) at T = 0.1, 1, 5 τ (τ = ℓ/w_plateau), and its plateau −(J + ξ1)/(2ξ2)/√(3/2) equals our root (−A_sUa + √((A_sUa)² + A_sA_d a² sup))/A_d to 1e-13 (and aζ for A_s = A_d). Our BE (n equal substeps, exact Newton): relative error vs the exact ODE, case A (U 0.5, w0 = 0.01 w): T = τ: 1.9e-1 (n = 1), 4.9e-2 (4), 1.1e-2 (16), 2.8e-3 (64); T = 5τ: 4.8e-2, 8.6e-4, 3.2e-6, 4.4e-8; case B (w0 = 3 w, above the plateau) T = τ: 7.5e-1 (n = 1), 1.6e-1, 3.2e-2, 7.3e-3. i.e. our implicit step at Δt ≳ τ is far from MESA's exact-in-step solution (MESA's step is exact for frozen coefficients; ours is BE and relies on the substeps of the outer solver); the long-time limit coincides. Decay (13): re-derived [derived] w = −[ξ1 + |J| tan(θ0 − |J|t/4)]/(2ξ2), tan θ0 = (w0 − c)/q', agrees with RK4 to 1.5e-14 (U = 0, 0.5) and reaches 0 at t* = 1.953 and 1.2189 ℓ/a; BE reaches 0 within 1–2.6 dt. For overdamped decay (J² > 0 with sup < 0, e.g. U = 2, sup = −0.05, both roots negative) (12) goes negative and must be clipped to 0 (not stated in the paper); exact t* = 0.9038, BE hit at 0.906 (dt 1e-3).

### M3a.5 Issues for the integration [interpretation]

1. **α_p default** (above): α_p = 1 is Kuhfuß-isotropic for our w; test 0, 2/3, 1.
2. **Active-set Newton at CFL ≫ 1.** The FB rows are semismooth; the active set (w = 0 where the stable region ahead of the layer is) moves one cell per Newton iteration where the flow crosses cells faster than Δt allows. T6 with 100 cells across L and Δt = domain transit time (CFL 300) needed ≥ 80 iterations (failure at 80; dev run, then switched to CFL ≤ 5 where ≤ 26 iterations suffice; with sup > 0 everywhere, T7, CFL = nz still takes 6). In TITAN the hydro step is limited by the sound speed, not by u only, so CFL ≫ 1 for the (subsonic) w is the normal case in a hot core; if the Newton count is a problem, take the w rows from the previous iterate's active set (a few outer iterations) or limit Δt by the layer transit time.
3. **No exact zero.** Newton lands on w = 0 only to the tolerance (observed −1e-17 … −1e-30 in the scalar local solves; exactly 0.0 in the grid cases of T5/T6). Nothing in F_c, P_t (w²), or the diffusion needs the sign; the integration should not clamp the residual, only post-process (store max(w, 0)) between steps.
4. **Cancellation in Ỹ.** For U² ≫ sup, F_c = A_c ρc_PT w(ℓ/2H_P)(sup − 2Uw/a) loses 4εU²/sup relative; the radiative flux dominates there, so the absolute error is negligible, but diagnostics F_c/F_MLT in that regime should use Ỹ = (A_d/A_s)(w/a)² at steady state.
5. **Newton tolerance and step size.** A relative tolerance 1e-12 on φ cannot be met when Δt ≫ τ_relax (noise ε Δt|R|/w); the harness stops on stagnation of the step (1e-13 w) and φ ≤ 1e-8 relative. The BE one-step solution is not exact for Δt ≳ τ (tables above).
6. **Boundary CVs and gas weights.** Interface 1 and nz+1 are Dirichlet w = 0 here; their Q (the flux into the wall) goes entirely to the adjacent zone. In TITAN's free-surface/inflow boundary variants the outer face flux of the boundary CV is nonzero and must be added to the budget.
7. **Docstring sign in M2's `blob_sup`.** The derivation in the docstring (and D1) gives sup = (∇ − ∇_ad) + (φ/δ)∇_μ with ∇_μ = Δlnμ/ΔlnP; with |ΔlnP_hs| = −ΔlnP the algebra yields sup = (∇ − ∇_ad) − (φ/δ)∇_μ (the Ledoux criterion ∇ > ∇_ad + (φ/δ)∇_μ). The code is right (T5: a stabilising Δlnμ < 0 lowers sup from +0.05 to −0.05); only the printed "+" is a typo. Not edited (M2 text frozen).
8. **ℓ, a, U at interfaces and ρ ℓ at zone centres** are inputs (`TDCBack`); the integration must supply zone ρ ℓ for the diffusion flux and the interface sup from `blob_sup` with the iterate (differentiated).
