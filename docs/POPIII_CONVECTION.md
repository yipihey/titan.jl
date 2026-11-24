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
- [read/derivation] In hydrostatic layers, where Δln P = Δln P_hs and Γ₁,b is the frozen value, sup = ∇ − ∇_ad − (φ/δ)∇_μ (Ledoux: unstable if ∇ > ∇_ad + (φ/δ)∇_μ; sign corrected 22 Nov, see M3a item 7).
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
2. **Δlnμ.** The `Δlnμ` argument of `blob_sup` is the remainder of the actual density contrast, Δlnρ − α_fr ΔlnP + δ_fr ΔlnT (φ_fr = 1). Then A = Δlnρ − ΔlnP/Γ₁,b exactly; α_fr and δ_fr cancel in A. Reasons: (i) the linearisation α ΔlnP − δ ΔlnT has second-order errors of the size of the signal where ΔlnT ~ 1 (ΔlnT = −0.7 and −1.3 at interfaces 40, 41); (ii) a μ from the zones' EOS would count the fast species twice (they follow ρ, T in the CK expansion and are already in α, δ). In a hydrostatic layer it is the same sup = ∇ − ∇_ad − (φ/δ)∇_μ [derivation, `blob_sup` docstring].
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
2. **Active-set Newton at CFL ≫ 1.** The FB rows are semismooth; the active set (w = 0 where the stable region ahead of the layer is) moves one cell per Newton iteration where the flow crosses cells faster than Δt allows. T6 with 100 cells across L and Δt = domain transit time (CFL 300) needed ≥ 80 iterations (one run hit the 80-iteration cap; development run; the committed T6 uses CFL ≤ 5 and then has no Newton failure; with sup > 0 everywhere, T7, CFL = nz still takes 6 iterations). In TITAN the hydro step is limited by the sound speed, not by u only, so CFL ≫ 1 for the (subsonic) w is the normal case in a hot core; if the Newton count is a problem, take the w rows from the previous iterate's active set (a few outer iterations) or limit Δt by the layer transit time.
3. **No exact zero.** Newton lands on w = 0 only to the tolerance (observed −1e-17 … −1e-30 in the scalar local solves; exactly 0.0 in the grid cases of T5/T6). Nothing in F_c, P_t (w²), or the diffusion needs the sign; the integration should not clamp the residual, only post-process (store max(w, 0)) between steps.
4. **Cancellation in Ỹ.** For U² ≫ sup, F_c = A_c ρc_PT w(ℓ/2H_P)(sup − 2Uw/a) loses 4εU²/sup relative; the radiative flux dominates there, so the absolute error is negligible, but diagnostics F_c/F_MLT in that regime should use Ỹ = (A_d/A_s)(w/a)² at steady state.
5. **Newton tolerance and step size.** A relative tolerance 1e-12 on φ cannot be met when Δt ≫ τ_relax (noise ε Δt|R|/w); the harness stops on stagnation of the step (1e-13 w) and φ ≤ 1e-8 relative. The BE one-step solution is not exact for Δt ≳ τ (tables above).
6. **Boundary CVs and gas weights.** Interface 1 and nz+1 are Dirichlet w = 0 here; their Q (the flux into the wall) goes entirely to the adjacent zone. In TITAN's free-surface/inflow boundary variants the outer face flux of the boundary CV is nonzero and must be added to the budget.
7. **Docstring sign in M2's `blob_sup`.** The derivation in the docstring (and D1) gives sup = (∇ − ∇_ad) + (φ/δ)∇_μ with ∇_μ = Δlnμ/ΔlnP; with |ΔlnP_hs| = −ΔlnP the algebra yields sup = (∇ − ∇_ad) − (φ/δ)∇_μ (the Ledoux criterion ∇ > ∇_ad + (φ/δ)∇_μ). The code is right (T5: a stabilising Δlnμ < 0 lowers sup from +0.05 to −0.05); only the printed "+" is a typo. Not edited (M2 text frozen).
8. **ℓ, a, U at interfaces and ρ ℓ at zone centres** are inputs (`TDCBack`); the integration must supply zone ρ ℓ for the diffusion flux and the interface sup from `blob_sup` with the iterate (differentiated).

## M3b — TDC in the implicit Newton system (22 Nov 2026, `TDC=1`, provisional, default off)

Labels: [run] [test] [read] [derivation] [estimate] [interpretation]. Logs in `docs/logs/subconv/m3b/`; run scripts `scripts/subconv/m3brun.sh` (one run), `m3blane.sh` (a queue on one 8-core slot),
analysis `m3bsnap.jl`, `m3bnwt.py`, `m3betot.py`. Branch `popiii-subconv`, TitanCore: see each log (restart tests: `TitanCore.jl` of aace2a0 in a scratch worktree, never committed; IC tests: the branch's).

### M3b.1 What was added [read]

| piece | where |
|---|---|
| `TDC=1` ⇒ `meqn` +1 (needs `COUPLED=1`, `TITAN_NOLIB=1`); with MGROUPS behind the groups | `TitanFortran.jl:115` |
| Switch, iterate `WTN`, old level `WTO`, `SAVEW`, coefficients `TDCP`, face weight `TDC_CEN`, newton hooks `tdc_alloc!`, `tdc_newstep!` (old := iterate), `tdc_reset!`, `tdc_update!` (w ← w + δ·unom), `tdc_dmax!`/`tdc_weight` (|δw| against max(w, 0.01 c_s): a weight, nothing is clamped) | `TitanResidFast.jl:928–963` |
| `_tdc!`: the block of `residual!` (before the GRID rows): the w rows φ_k = fb(w_k, F_k) (`TDC.tdc_F!`) into `R[iw, ·]`, F_c → `fcv`, zone P_t → `W.tPt`, exchange Q → `W.tgx`, D_mix → `W.tD` | `TitanResidFast.jl:973`; call `:1612` |
| code 8 for the w column (linear, scale `unom`, base `WTN`); `convcoef` skipped (no τ weights); MLT block skipped | `:349`, `:375`, `:1772–1773` |
| momentum row k: `+ rmu_k (P_t,k − P_t,k−1)`, θ-centred like p_g | `:1744` |
| gas-energy row of zone j: `+ ½(Q_j + Q_j+1)` (boundary CVs wholly to their zone) | `:1971` |
| species rows 8, 9 of a coupled zone: `+ F_k+1 − F_k`, F_k = −rmu_k ρ_k D_k (x_k − x_k−1)/(r_c,k − r_c,k−1) | `:2008` |
| TR-BDF2: `WTO` saved/restored with the old arrays; the BDF2 stage's old w is w* = w_γ + B (w_γ − w_n) | `Trbdf2.jl:60–62, 227` |
| Newton–Krylov: `WTN` in the snapshot (`extras`), code 8 active on interfaces a … b+1, weight of w, `realized` = applied as given | `NewtonKrylov.jl:63, 148, 169, 250` |
| `TDCPAR`, install, `CPLH` hooks, checkpoint `tdc`, restart, rebuild reset | `popiii_m2c.jl:2002–2022, 844, 1067` |
| ETOTMON: E_t = Σ m_k w_k² and Σ Q Δt into `etot_monitor` (`TDCMON`) | `TitanResidFast.jl:2303`; `ETOTLOG=<file>` per-stage lines `popiii_m2c.jl:2056` |
| CONVDIAG records `tdc_w, tdc_PtP, tdc_Dmix, tdc_Q, tdc_Fres, tdc_a, tdc_U, e_tdc` (+ `fconv`, `vconv`, `fpot`, `vpot`, `sup` as for MLT); layer scalars E_t, P_t/P at the post-shock zone; runrec exports `cv_w, cv_PtP, cv_Dmix, cv_Q`, `cv_Et, cv_PtP_ps, cv_PtP_ps3` | `popiii_m2c.jl:453, 512`; `runrec.jl:95, 99` |
| `tdc_exchange(…; Ptc)`: the P_t of the P_dV term may be the momentum row's θ-centred P_t | `TDC.jl:257` |
| NKCHECK: `NKEQ=<eq,…>` row-by-row listing | `NewtonKrylov.jl:742` |

### M3b.2 Choices [interpretation, with the reasons]

1. **w is a linear unknown with the scale of u** (`unom`), last index `iw = neqn` (10 in the production stack). Dirichlet rows φ = w at the two boundary interfaces (w_base = 0 there), trivial rows (δ = 0) in the ghosts, and for an interface without gravity or with an unphysical thermodynamic state (φ = w ⇒ w = 0, no flux; the step is rejected by the other rows in the second case, as for MLT).
2. **Everything at the iterate, nothing held.** The background (ρ, T, P = p_gas + f E_r, g, H_P, ℓ, κ = `avchi`, sup, δ_b, c_P,b) is evaluated inside the residual from the Newton iterate with Duals; ℓ, g, H_P, ΔlnP_hs exactly as the MLT block (MLT_RCAP, MLT_HSDISC); sup, δ_b, c_P,b from M2's `buoy_blob` (so TDC contains the D1 buoyancy of `MLTBUOY=1`, whatever that switch is). The equilibrium fraction keeps M2's a-based turnover time f = ℓ/(ℓ + t_chem a) (no fixed point in w; a ≥ w for sup ≤ 1, so f is a lower bound) [interpretation; an f with w would add one more column per row and is a candidate for M4].
3. **Dual thermodynamics only where a row or flux can depend on them (exact, not a gate).** The zones' Float64 thermodynamics screen every interface (sup's value is the true one); the Dual thermodynamics are added where sup > 0 or w, wo ≠ 0 at k−1, k, k+1. Elsewhere w = 0 and F > 0 (sup ≤ 0, nothing turbulent near: F = Δt (A_s a² |sup|/ℓ) > 0): fb(0, F) has ∂φ/∂F = 0, so no partial of F, hence none of sup, enters the row, and fcv = ρ c_P T w ℓ Ỹ/2H_P has w = 0 as the factor.
4. **Zone quantities**: ρ_j is the new-level zone density; ℓ_j is the mean of the two interface ℓ (the one-sided value at the first and last zone); the zone-centre area is ((r_j + r_j+1)/2)^μ (θ-centred radii).
5. **w is backward Euler** in every stage (also the θ = ½ stage of TR-BDF2; as the species with `CPL_THETA = 1`), with the interface-CV mass fluxes, masses and donor-cell faces of the momentum row (`TDCCEN = 0`). The BDF2 stage's old state is the linear combination of w (w* = w_γ + B(w_γ − w_n)); it can be negative for a decaying w (then φ is the FB row with w_o < 0; not clamped). **w is not in DTERR's error estimate or in `combos`** (the estimator's filter solves keep WTO = y*); `DTSCHEME` ≠ trbdf2 stops with an error.
6. **P_t**: α_p ρ ē with the θ-centred combination of the new and old level (old: ρ_o, w_o), discretised as the p_g term of the momentum row; the same P_t and dV/dt = rmu_j+1 u_j+1 − rmu_j u_j (θ-centred, the factor of the gas row's p_g term) enter Q (compatible energy). α_p = 1 [M3a.2: isotropic for our w].
7. **Gas coupling**: Q_k = realised exchange of the interface CV (TDC.tdc_exchange); the gas row of zone j gets + ½(Q_j + Q_j+1), the boundary CVs wholly their zone, so Σ_j (gas source) = −Σ_k Q_k exactly.
8. **Species**: only the coupled species 8, 9 get the diffusion flux (D = α_D ℓ w, interface ρ the mean of the two zones, new level); the flux through the domain boundaries is zero. The other network species are not mixed here (they are operator split; their mixing with the converged D_mix is not implemented) [open]. The hot-core `mix!` is untouched.
9. **ETOTMON**: the rows carry + ΣQ; the monitor books E_t = Σ_k m_k w_k² instead (Q Δt = ΔE_t + Δt Σ P_t dV): imb → imb + ΔE_t − Δt ΣQ, thr → thr + |ΔE_t| + Δt Σ|Q|.
10. **α_m (eddy viscosity) is not implemented** (`TDCPAR` 7th entry must be 0). Boundary CVs are Dirichlet (w = 0), no PEXT.
11. **Checkpoint/restart**: the checkpoint's aux carries `tdc = (WTN, WTO, a, b)`; a restart restores them when TDC is on and the zone range matches, and starts from w = 0 for a checkpoint without them (or from a different mesh). A restart with TDC off ignores them (and resets neqn to 9). A RAD rebuild or rezoning (`rebuild!`) resets w to 0 (the mesh changes; REZONE = 0 in production).

### How to run TDC

`TDC=1` in the environment of the zie1 stack (needs `HOT=1 COUPLED=1 TITAN_NOLIB=1 DTERR=… NKRYLOV=1`; `RAD=1`); `TDCPAR=A_s,A_d,A_c,α_p,α_t,α_D,α_m` (default `1,1,1,1,0.01,1,0` = the B2 best estimate); `TDCCEN` face weight (0 = donor); `CONVDIAG=1 RUNREC=<dir>` records the diagnostics, `ETOTMON=1 ETOTLOG=<file>` the total-energy ledger, `NWTRACE=<file>` the Newton iterations. `scripts/subconv/m3brun.sh <name> <steps> RESTART=rp_k650 TDC=1 CONVDIAG=1`.

### M3b.3 Tests [test/run; logs in `docs/logs/subconv/m3b/`]

**(a) Bit-identity, TDC unset** [test, `bitid_ic40.txt`]: IC, zie1 settings, 40 steps, SNAPEVERY=5, HEAD b2c9799 against bcd248e (the commit before M3b; scratch worktree), the branch's TitanCore in both: 10 snapshots, 300 fields, 0 differences. `neqn` = 9 with TDC unset (the log of the first line of each run). The restart path was not re-run for bit-identity (restarts of M1/M2 are exact only with aace2a0's TitanCore).

**(b) NKCHECK** [run; `nkcheck_*.txt`; restarts of 20 TDC steps from rp_k650 (`tdc` checkpoint, w ≈ 2–4×10⁵ cm/s at interfaces 39–41, see the table of pilot (d)) and rp_k290; `NWPROBE=1e5 NWPROBEAT=0`, ε = 1e-5; TitanCore of aace2a0]. Max rel = |Jv_FD − J_b w|/Σ|J_b w| per row, all zones | post-shock zones (650: 34–44; 290: 33–43), TDC (`nk650t`, `nk290t`) against the TDC-unset baseline from rp_k650/rp_k290 (`nk650b`, `nk290b`):

| eq | 650 TDC all | 650 TDC 34–44 | 650 base all | 650 base 34–44 | 290 TDC all | 290 TDC 33–43 | 290 base all | 290 base 33–43 |
|---|---|---|---|---|---|---|---|---|
| r | 5.0e-6 | 4.4e-6 | 5.0e-6 | 4.9e-6 | 5.0e-6 | 4.7e-6 | 5.0e-6 | 4.0e-6 |
| m | 3.8e-6 | 3.8e-6 | 3.9e-6 | 3.1e-6 | 4.5e-6 | 4.5e-6 | 3.6e-6 | 3.6e-6 |
| ρ | 4.9e-6 | 3.9e-6 | 4.7e-6 | 4.5e-6 | 4.9e-6 | 4.4e-6 | 4.5e-6 | 3.3e-6 |
| **u** | 9.8e-6 | 9.8e-6 | 4.2e-6 | 3.0e-6 | 3.5e-6 | 3.5e-6 | 4.2e-6 | 3.9e-6 |
| **T** | 7.9e-6 | 7.9e-6 | 1.7e-2 (zone 1) | 3.4e-6 | 4.6e-6 | 4.5e-6 | 2.6e-2 (zone 1) | 3.8e-6 |
| E_r | 7.9e-6 | 7.9e-6 | 9.1e-4 | 9.1e-4 | 7.6e-5 | 7.6e-5 | 1.1e-5 | 1.1e-5 |
| F_r | 6.9e-2 | 5.0e-2 | 17 | 17 | 0.20 | 0.20 | 0.51 | 0.15 |
| **x_H⁺** | 9.2e-4 | 3.7e-6 | 2.7e-4 | 4.9e-6 | 1.6e-3 | 1.6e-3 | 2.5e-3 | 2.5e-3 |
| **x_H₂** | 3.9e-6 | 3.9e-6 | 3.8e-6 | 2.2e-6 | 4.3e-6 | 4.3e-6 | 3.5e-6 | 3.5e-6 |
| **w** (row 10) | see below | see below | – | – | see below | see below | – | – |

[interpretation] The rows with new terms (T: F_c and Q; u: P_t; species: mixing) agree with the FD product at the FD noise floor (≲1e-5) in the post-shock zones, and wherever they are not (E_r, F_r, x_H⁺ at 1e-3: closure and chemistry rows, zone 1 of T) the TDC-unset baseline shows the same (docs/POPIII_NEWTON_KRYLOV.md). The baseline's T mismatch in zone 1 (1.7e-2) is absent in the 650 TDC restart (a different state). Row probes (`NKROW`, one unit vector per column, `nkcheck_nk650t.txt`, `nkcheck_nk290t.txt`, `nkdbg`): the w rows of interfaces 40, 41 (650) and 34 (290) match the band row in every column to ≤1e-5 absolute differences relative to entries of 1e4 (e.g. `column T zone 39: FD −5.1637e4, band −5.1637e4`); the centre rows (interface 2, T zone 1–2, ρ, u, m at k115 state dbg_k153): rel ≤ 3e-6 in zones 1–8.
**The w rows** (`NKEQ=10`, `nkcheck_nkeq650.txt`, `nkcheck_nkeq290.txt`): where w ≥ 7×10³ cm/s (650: zones 39–42, 290: zones 34–36) rel = 1.7e-6 … 1.5e-5 (the largest, 1.5e-5, in the interface with the smallest w, 7.4e3, next to a w = 0 interface), where w = 0 and F ≫ 0 and the neighbours have w = 0 (650: zones 36–38) 1e-10. **Not at the noise floor: the rows whose w is ≪ the FD step** (ε v unom ≈ 10 cm/s): 650 zones 43–75, 290 zones 37–42, … (max rel 2×10⁴ … 3×10⁶). There the FB row has its kink at w ≈ −F₀/κ ~ 10⁻⁷ cm/s from the base point (F₀ = Δt A_s a²|sup|/ℓ ≈ 50 cm/s, κ = ∂F/∂w ≈ 6×10⁸ for U ~ 10¹⁰), so the central difference straddles it and returns the average of the two one-sided slopes: the signature is rel = 0.5000 exactly in the rows with 0 < w ≪ 10 (650 zones 44–46: 0.5000, w = 7e-4, 2e-5, 7e-8) and arbitrary large values at w = 0 with F > 0. The state columns of these rows are exactly 0 in J_b (φ = fb(0, F > 0) does not depend on F) and the w column is 1 [derivation above, choice 3]; no FD with a step ≪ 10⁻⁷ cm/s exists in double precision, so these rows are **verified by the derivation, not by FD** [limitation]. Iteration 2 (`nkit2_650`, `nkit2_290`, `NWPROBE=1e6`, NKCHECK=2): same picture, the other equations 1e-5 … 3e-4 (u at zone 31: 3.3e-4, x_H⁺ 1.5e-2 at 290 where Newton's iterate is far from the root; the baseline at iteration 2 was not run).

**(c) Energy** [run; `etot.log` per stage, `pilots_summary.txt`; ETOTMON: gas + radiation + kinetic + (TDC) E_t; Σ over the accepted TR and BDF stages of the 200 steps]: |Σ imb|/Σ thr = 6.1e-3 (TDC) against 5.6e-3 (baseline) at k650, 5.8e-3 against 5.7e-3 at k115 (100 steps), max |G| 1.0e-2 against 1.0e-2 (k650); k290: 3.8e-3 against 3.6e-3, k2395 4.5e-3 (TDC). [interpretation] The ETOTMON imbalance of the stack is dominated by the existing ~5e-3 of throughput (radiation/chemistry rows; docs/POPIII_ENERGY.md); TDC adds no measurable drift (E_t is in the ledger, Σ E_t is 1e43–1e44 erg sr⁻¹ against Σ|dE| 6e46–6e47, i.e. the turbulent energy is 1e-4 of the energy exchanged). The exact identity (Σ gas source = −Σ Q, ΔE_t + Δt Σ P_t dV = Δt Σ Q) is a property of the rows, checked in M3a T8 to 5e-17; the monitor's residual cannot resolve it against the stack's own 5e-3.

**(d) Pilots** [run; `pilots_summary.txt`; 200 steps from rp_k650 and rp_k290, 100 from rp_k115 and rp_k2395; α_p = 1, α_t = 0.01, α_D = 1, A_s = A_d = A_c = 1; baseline = TDC unset; B1 = MLTBUOY=1 alone, from `docs/logs/subconv/m2/pilots_summary.txt` (same TitanCore, same stack; not re-run)]:

| restart (M☉), steps | Newton-failed attempts TDC / base / B1 | min snapshot dt [s] TDC / base / B1 | F_c/F_tot in the layer | w in the layer [cm/s] (max) | P_t/P post-shock | w at M > 2 inflow | R*, Ṁ |
|---|---|---|---|---|---|---|---|
| k650 (3.2), 200 | 29 / 33 / 30 | 1.1e5 / 1.2e5 / 1.5e5 | 0.995 → 1 (655–680), 0.2–0.4 after R* jumps 278 → 359 R☉ (695), 1.0 at 785–790, then 0 at 820 | 1.4e5 → 4.0e5 (interface 41/42); 0 from 820 (layer dissolves) | 1.6–13 % (6.9 % at 660) | 2e-6 … 7e-4 | 278 → 359 R☉ at 695 (as baseline/B1), Ṁ 2.10e-2 → 1.58e-2 M☉/yr (TDC, 850) |
| k290 (0.29), 200 | 30 / 32 / 18 (100 steps) | 1.7e4 / 5.8e3 / 2.2e4 | 1 at interfaces 33–36 (B1: 1 at 35–36) | 1.7e4 → 2.05e5 (interface 36) | 1.6–2.8 % | 1 … 17 (max over M > 2 interfaces); 4×10⁴ at the sonic interface 37 (M = 1.0) | 59.9 R☉, Ṁ 3.24e-2 → 2.77e-2 |
| k115 (0.017), 100 | 26 / 12 / 13 | **8 s** / 9.5e4 / 6.9e4 | 1 (interfaces 2, 25–28) | 1.4e5 – 2.0e5 | 0.4–1.8 % | **1.2e3 … 2.7e4** (M > 2; up to 16 % of the layer's w) | 4.6 → 10 R☉, Ṁ 7.9e-2 → 7.7e-2 |
| k2395 (34), 100 | 21 / – / – | 4.7e5 | 1 at interface 7 (core), no layer | 1.7e2 (core) | – (no shock) | 1.7e4 at interface 43 (M = 4.3, F_c = 2×10⁵ erg s⁻¹ cm⁻², ≤ 10⁻⁶ F_tot) | 215 R☉, Ṁ 1.1e-2 |

[interpretation] **k650:** the layer is turbulent from step 655 with F_c/F_tot ≈ 1, as in B1 (0.999), P_t/P up to 13 % (a dynamical effect: α_p = 1 is the Kuhfuß isotropic value; with 0.666 it would be 2/3 of that); the cycle of B1 (layer destroyed by the shock jumping outward, rebuilt) is reproduced; after step 820 the layer is not rebuilt (w → 0 in 5 steps: sup < 0 everywhere) as in B1 at 835. Newton failures are not more frequent than in the baseline (29 vs 33). **w in the supersonic inflow:** at k650 w at interfaces with M > 2 stays ≤ 7×10⁻⁴ cm/s (10⁻⁹ of w_s ≈ 10⁵) and at k290 ≤ 17 cm/s, without any gate (T7 of M3a in the full code); the sonic interface at k290 (M = 1.0) carries 4×10⁴ (20 % of the layer's w), one interface further in (M = 3.5) 3.6 cm/s. At k115 (a weaker shock, 0.017 M☉, cold inflow with U ≈ 0.1–2) w penetrates the supersonic inflow: 1.2×10³–2.7×10⁴ cm/s at M = 2–5 [not ≪ w_s there]. **k2395:** the stellar core's one marginally unstable interface (sup = 1.9e-6) has F_c = 1.01 F_MLT (local limit, below); in the supersonic hot inflow (M = 4.3, U = 4.4, T = 1.2×10⁴ K) the buoyancy is unstable (sup = 0.014) and w = 1.7×10⁴ (1.3 % of w_s = a√sup): the flux is 10⁻⁶ of F_tot; the TDC without a τ gate does not remove it, the radiative damping U is only 4. No spurious convection of consequence at 34 M☉. **k115: a failure of the Newton iteration** (open problem 1): from step 154 the TR stage fails at every Δt down to 26 s (dt collapse to 8 s), 26 failed attempts against 12 in the baseline, and 6.7 Newton iterations per attempt against 4.0. See M3b.5.
Cost [run]: 10 steps from rp_k650, wall time including the start-up: 77 s (TDC) against 37 s (baseline) (`t650_10`, `b650_10`); in the 200 steps of k650 the zone thermodynamics took 1976 s Float64 (1.6×10⁶ evaluations, 1.2 ms each: there is no τ gate, every interface is screened) plus 1423 s with derivatives (`BUOYT`, `pilots_summary.txt`), for a run of about 100 min.

**(e) Local MLT limit in the full code** [run; `local_mlt_limit.txt`, `m3blocal.jl`]: F_c(TDC)/F_MLT of the same sup, ℓ, U at the interface (`cd.fpot`), subsonic interfaces carrying > 1 % of the flux: k2395 core (hydrostatic, steady, sup = 1.9e-6): **1.0107**; k650 layer at step 810 (slowly decaying, sup = 0.066): 0.9856; k290 at step 490: 1.59 (interface 35 with sup = 4×10⁻³ carries the w advected and diffused from the interface 36 with sup = 0.05: not a local state). So the ratio is 1 where advection and diffusion are negligible (1.1 % in the core: the residual of the BE step, A_d/A_s = 1 and Ỹ with U = 10⁻¹⁴).

**(f)** `popiii_regress.jl` only reads stored snapshots and `noh_test.jl` cannot be started here (TitanFortran library not loadable; MLT off): neither exercises the changed code [read, as M2]; not run.

**(g) Newton behaviour** [run; `m3bnwt.py` from `NWTRACE`]: iterations per converged attempt (mean/median/max): k650 TDC 2.98/2/24 against base 2.63/2/17; k115 6.68/5/24 against 3.96/3/23; k2395 1.84/2/3; k290 TDC 8.0/3/25 against base 1.99/2/19 (>10 iterations: 144 of 400 attempts against 1). In k290 the long attempts (144 of 400 with > 10 iterations) are limited by the shock-zone rows (the largest raw corrections of a sampled long attempt, `tdc290_200/nwt.txt` attempt 248: T 0.07–0.11, E_r 0.5, the species ±0.5 (capped) at zone 37, the sonic point; the w correction there is 4×10⁻³–6×10⁻³ of unom and not the largest). I did not see the "one cell per iteration" creep of the active set (M3a issue 1) as the limiter in the samples read, but I did not separate it from the shock-zone iterations statistically **The TDC raises the Newton iterations per attempt at k290 four-fold (8.0 against 1.99) with the same number of failures (30 against 32), and at k115 by 1.7× (6.7 against 4.0); at k650 by 13 % (2.98 against 2.63)** [run]. I have not isolated the cause; the long attempts are not limited by the w row itself.

### M3b.4 The k115 Newton failure [run; `nkcheck_nkdbg.txt`, `dbg_*` runs summarised here; not solved]

State: the replay from rp_k115 with TDC; at step 154 (R* = 10 R☉, the whole 0.017 M☉ star convective, w ~ 10⁵ cm/s in the layer and 10³–3×10⁴ cm/s in the cold supersonic inflow, interfaces 30–33) every TR stage fails after 25 iterations at every Δt from 1.6×10⁴ s down to 26 s (5 retries, then a converged step at Δt = 26 s … 7 s; 26 failed attempts in 100 steps, 12 in the baseline; `tdc115_100/nwt.txt` attempts 88–92). Restart from the checkpoint written after step 153 (`dbg115`, `CKDUMP=150:153`) reproduces it (`dbg_v0`: the attempts of step 154 fail at 1.6×10⁴ … 244 s again). Findings:
- **Not the Jacobian.** NKCHECK at that state (`nkdbg`: dt 1.5e4, iteration 1): all equations ≤ 5e-6 in zones 1–8 (centre), the w rows of zones 2–8 ≤ 1.3e-6, the row probes of interface 2 match column by column; globally E_r 1.4e-2 and F_r 0.34 in the radiation rows of zone 32–33 (the closure rows, as in the baseline).
- **Not P_t and not the e_t diffusion.** `TDCPAR` with α_p = 0 (`dbg_v1`) and α_t = 0 (`dbg_v3`): the same five failures at step 154 and the same iteration counts. With α_D = 0 (`dbg_v2`) the failure moves to step 153 (two failed attempts, then converged at 397 s): the species mixing matters but is not the only contribution.
- **Where the iteration fails.** The raw corrections of the failing attempts (`NWTRACE`): F_r (equation 7) of zones 34–100 (−1 … −200), E_r and T of zones 29–32 (0.1–0.6), x_H⁺ of zone 30 pinned at the ±0.5 trust-region cap, and w and u of interface 2 and ρ of zone 28 growing by ×1.5 per iteration (w: 1.8e-3 → 2.5e-2 of unom in 14 iterations), `dmax` stuck at 2.15 with the trust-region factor 0.01–0.1: a nonlinear divergence of the whole system from a state from which even the Δt → 0 limit does not converge (the old state is the start and the residual at Δt = 26 s is still not reducible). The baseline from the same checkpoint was not possible (the state is TDC's).
- [interpretation, not shown] The candidates are (i) the supersonic cold inflow at M = 2–5 where TDC has w ~ 10³–10⁴ cm/s with U ~ 0.1–4 (no gate) driving the species and radiation rows through F_c, Q, D_mix, (ii) the chemistry's switch inside `buoy_blob` (f = 0 where the equilibrium table is unavailable: `BUOYEQFAIL`, 1.6×10⁵ interfaces in this run) making sup discontinuous in the iterate. (ii) is a property of M2 (the frozen-blob fallback), is not smoothed here, and would show as a jump of sup between iterates; I have not measured it.

### M3b.5 Open problems [interpretation]

1. The k115 failure above (dt collapse to 8 s, 26 failed attempts); k650, k290, k2395 have none.
2. w is not in DTERR's error estimate (`combos`, `old_from!`): the step size does not see w's own truncation error; for the unstable-layer growth (t90 = 4.2/N) the step control comes from T, ρ, x only. Adding w means a conserved combination (m w) in `combos` and a weight in `wnorm`.
3. FD cannot verify the stiff w rows (above); the derivation says state columns = 0 and w column = 1; a complex-step or a dual-number Jacobian of the FB row alone against `TDC.tdc_rows!` (M3a T8 did this on a harness: 0.0 difference) is the available check.
4. Cost: every interface is screened with Float64 zone thermodynamics at every residual evaluation (1.2 ms per zone); caching per iterate/stage (the zone thermodynamics depend only on (ρ_j, T_j, x_j) of one zone, which the Krylov products change by 1e-5) is the obvious saving, unimplemented.
5. The equilibrium fraction f is M2's a-based one (not w-based); the fallback f = 0 where the table is unavailable is a discontinuity (≥ 1.3×10⁵ interfaces per 100 steps).
6. Not implemented: α_m; species mixing of the non-coupled network species and of the hot-core `mix!` by D_mix; PEXT; DTSCHEME ≠ trbdf2; bit-identity of a *restart* with TDC unset against the commit before.
7. Rezoning/RAD rebuild reset w to 0 (never exercised: REZONE = 0).
8. The k2395 supersonic hot inflow (M = 4.3) is unstable in the buoyancy sense (sup = 0.014) with U = 4: w = 1.7×10⁴ (1.3 % of w_s) there, F_c ≤ 10⁻⁶ F_tot: harmless here, the physical suppression (Henyey thin-blob correction, M3a open item) is still missing.

## M3c — The k115 Newton failure is not in the TDC rows (23 Nov 2026; no solver change)

Labels: [run] [test] [read] [interpretation] [open]. Logs `docs/logs/subconv/m3dbg/` (`attempts.txt`: every Newton attempt of every run below; `interfaces_fA.trace.gz`: per-interface w, sup, F_c, D_mix, eqfail for each Float64 residual evaluation of the failing step; `nkscan_vscan_att3.txt`: NKSCAN of the first failing attempt; `debug_knobs.patch`: the temporary DEBUG knobs, removed again). Scripts `scripts/subconv/m3dbgrun.sh` (restart from a checkpoint, scratch worktree with TitanCore.jl of aace2a0), `m3dbgsync.sh`, `m3dbgsummary.sh`. Branch code at d891b47 for all runs.

### M3c.1 Method [run]

Two checkpoints written from the M3b replay (`SAVE=… CKDUMP=…`, TDC, rp_k115): `m3d_k154` (the state before step 154, the M3b failure) and `m3e_k140` (before step 140). All variants restart from the *same* state and compare the Newton attempts of the same step (the controller's dt sequence is identical until the first difference). The full TDC run is reproduced exactly (`fA`, `hA`: step 154 fails at dt 1.56e4, 3.9e3, 9.7e2, 2.4e2, 61 s with 26 iterations each, converges at 15 s; the original `tdc115_100` and `dbg_v0` sequences).

### M3c.2 Hypotheses H1–H5 against the fixed state [run; `attempts.txt`]

| test (state before step 154, 1 step) | result at dt = 1.56e4 … 61 s | 
|---|---|
| reference (`fA`) | 26-iteration failure at all five dt |
| **H1** interfaces switching between table-available and unavailable | not seen in any of 151 residual evaluations of the step: `eqfail` occurs only at interfaces 86–100 (ρ < 1e-20), where w = 0, sup < 0 and no row or flux depends on it, and no interface changes state between iterates [`interfaces_fA.trace.gz`] |
| **H2** production τ weight on F_c (earlier-trajectory run `vtau`, from dbg_k153) | same failures (the layer is optically thick, weight 1 where F_c ≠ 0) |
| **H3** the centre: F_c off on interfaces 2–25 (`vE1`, trajectory) | same failures. At the failing state the Newton corrections are not at the centre: weighted b is dominated by E_r, F_r of zones 31–37 (`nkscan_vscan_att3.txt`: 6@32 −175, 6@31 +116, 7@33 −78, then decaying outward); the M3b reading "w, u at interface 2" was the largest *raw* entries of tiny rows |
| **H4** F_c off at interface 29 (`fC`) or at 30…100 (`fD`) | identical iteration counts (26 ×5, then 13, 6) |
| H4 D_mix off (α_D = 0, `fB`) | identical failures (converges at 15 s in 20 instead of 13 iterations) |
| H4 all couplings off: no F_c, no Q, α_p = α_D = 0 (`fG2`) | identical failures |
| **TDC unset** (`fG1`: restart without TDC, production MLT, neqn 9) | **identical failures, identical iteration counts (26 ×5, then 13, 6)** |
| the same at the onset of the dt collapse, state before step 140 (`hA` TDC, `hG1` TDC unset) | both: step 141 fails at 1.8e5, 4.6e4, 1.1e4 s, converges at 2.9e3 s (11 iterations): identical |

[interpretation] **None of the TDC rows or couplings (F_c, Q, P_t, D_mix, the w rows themselves) is the cause of the failure at the state, at step 154 or at its onset, step 141.** The earlier ablations that appeared to cure it (`vAc0`, `vE2`, `vE4`–`vE6`: F_c off in part of the layer, started from dbg_k153) change the trajectory from step 153 on (the controller then starts step 154 at dt 3.2e3 against 1.56e4 and with a different state) and do not test the state; the fixed-state runs above do. The "α_D = 0 moves it to step 153" of M3b.4 is the same trajectory effect.

### M3c.3 Where it fails [run, NWTRACE/NKSCAN; interpretation]

- The Krylov solve is accurate (rel. residual ≤ 5e-3, 1–2 products; NKCHECK ≤ 5e-6 for the TDC rows in M3b). The Newton direction is the problem: the first correction at dt = 1.56e4 s is 253 in the weighted norm (E_r of zones 31–32 changes by ~5× their value, with sign alternation E_r +, −, − … and F_r alternating in zones 32–36) and the radius-0.5 global step is λ ≈ 0.003–0.01. The line scan along the step: T/T0 = 0.67 at λ = 1/256 and an explosion (T/T0 = 1.5e4) at λ = 1/128 where E_r(32) would become non-positive. The iteration walks along a trust-region-limited path of the radiation block and does not converge in 25 iterations. The initial weighted correction scales as dt^1.7 → dt^0.25 (0.75, 8.6, 31, 76, 160, 253 for dt = 15, 61, 2.4e2, 9.7e2, 3.9e3, 1.56e4 s).
- **The closure.** With the Eddington factor f held during the Newton solve (`NWFREEZE=f`, diagnostic of docs/POPIII_NEWTON_KRYLOV.md; `hF2`; all closure quantities `hF1`) step 141 converges at dt = 4.6e4 s in 18–19 iterations and at 2.3e4 in 18–19 (unfrozen: needs 2.9e3 s). So the nonlinear dependence of f (the formal solution at the iterate; non-local, with the branch decisions `lim`, `sgn` of the limiter flipping along the step) on a photosphere that is a one-zone cliff is what the globalised Newton cannot cross. That is the *known* radiation-closure Newton wall of docs/POPIII_NEWTON_KRYLOV.md §1 (the wall that made the Krylov solve necessary), not a convection effect.
- **The state.** From step 141 (`tdc115_100` snapshots) the TDC run sits at a photosphere (first zone with κ < 1e-2) fixed at zone 30 with E_r(29)/E_r(30) = 60–98 and T(29) = 8.7–10 kK → T(30) = 3.05–3.17 kK, opacity 19–70 → 1.3e-4–2.3e-4 cm² g⁻¹, F_c/F_tot = 0.9–0.97 at the interface below; the baseline has the same size of cliff but it moves outward (jp 29 → 35, R* growth) and cycles (E_r ratio 6 ↔ 125, step 120–215). The TDC photosphere does not move for 75 steps (the dt is 10³–10⁴ times smaller) so the cliff stays at the ratio where f is stiff. At steps 153–155 the radiation flux in zones 30–36 rings (F_r = +1.15e9, −1.26e9, +5.7e9, −1.8e8, −1.8e8, +3.1e8 at step 155 against a smooth 1.0–2.3e9 at steps 145, 150, 205): an odd-even mode of the closure that the large steps cannot damp.
- The baseline fails in the same way on this cycle (steps 132, 141, 143, 157 at dt 1e6, 1e6, 4e5, 1e6) and recovers after one factor-3–4 reduction; the TDC state needs 60× (141) and then 10⁴× (154).

### M3c.4 Conclusion and what was not done [interpretation; open]

- The cause of the k115 collapse is **not** in the TDC rows (H1–H5 refuted at the fixed state; TDC-unset restart reproduces it to the iteration), but in the **radiation-closure Newton wall at the one-zone photosphere**, reached because the TDC star does not do the baseline's photosphere cycle (the cliff stays at zone 30 where f is stiff). TDC influences the *state* through F_c = 0.9–1 F_tot below the photosphere; its rows are not part of the failing block.
- **No solver change was made, nothing to verify against bit-identity or the k290/k650 pilots** (the branch is at d891b47 + docs). The DEBUG knobs (`TDCDBG=log,tau,noFc,noFcR,noQ`, per-interface trace) are in `debug_knobs.patch` only.
- A principled fix belongs to the closure/globalisation work, not to M3: f's dependence on the iterate at τ ~ 1–10 cliffs (smoothing of the limiter/branch decisions of the formal solution, or a globalisation that follows the f-dependence: freezing f in the Newton solve is a lagged input, not acceptable as a method, M3 §1 of docs/POPIII_NEWTON_KRYLOV.md), after which the k115 TDC pilot should be re-run. Not attempted here: it changes production radiation transfer.
- [open] `NITER=100/300` did not change the 26-iteration cap in these restarts (`hN1`, `hN2`): the iteration cap that the trace shows is not set by `NITER`; whether a larger cap would let the trust-region path converge was not tested.
- [open] Why the TDC photosphere does not cycle (F_c carrying the flux up to τ ≈ 10 with a steady surface) was not investigated; B1 (MLTBUOY=1, τ weight, MLT_NOINFALL) ran the same restart without collapse and may differ exactly here.

## M3d — The TDC Newton wall: the species' local initial guess omitted the turbulent mixing (23 Nov 2026; unconditional inside `TDC=1`, the switch `TDCMIXGUESS` was removed)

Labels: [run] [test] [read] [derivation] [interpretation] [open]. Logs `docs/logs/subconv/m3wall/`; scripts `scripts/subconv/m3wall_*` (`m3wall_run.sh` restart from a checkpoint with `NWTRACE`, `m3wall_variants.sh` one process, one variant per probe,
`m3wall_att.py`, `m3wall_attshow.sh`, `m3wall_dtnwt.py` trace readers). Code of the tests: scratch worktrees of the branch (TitanCore.jl of aace2a0, never committed); DEBUG knobs in `debug_knobs.patch` only (removed again).

### M3d.1 The problem and the method [run]

M4 lane B2 (TDC=1, `TDCPAR=1,1,1,1,0.01,1,0`, zie1 stack) advances at Δt = 1.6e4–7.3e4 s with a Newton-failure rejection every ~4 steps (dt ÷4, then regrowth ×2 per step to the same ~7e4). Start state of all tests:
`m2c_B2S1_auto0` (B2 step 600, 0.29 M☉, R* = 59.9 R☉; the photospheric cliff of the layer is at zone 37: T 11 281 → 3 157 K, χ 122 → 3.4e-4 cm⁻¹, x_H⁺ 0.36 → 4.5e-6 → 3.8e-11 over zones 36 → 38,
`profile_step605.txt`, `xHp_zone36_38.txt`). Tools: `NWTRACE` (every Newton iteration: largest weighted entry, λ, trust radius), `NWPROBE` (one trapezoidal-stage attempt from the live state at fixed Δt; stage = 0.586 Δt; 26 iterations = failed),
`TDCDBG=noFc,mltFc,noQ,wfrz` and per-probe `TDCVARS` (DEBUG, scratch only): F_c → 0, F_c → the local-MLT value of the same sup, U, ℓ at the iterate (`kw_mlt_flux`), exchange Q → 0, the w row replaced by w = w_old; α_p, α_t, α_D through `TDCPAR`.

### M3d.2 Findings [run]

1. **The wall is a window in Δt, not a w-rows problem** (`probes.txt`; iterations of the stage Newton at probe Δt):

| probe Δt [s] | 3e4 | 6e4 | 1e5 | 2e5 | 3e5 | 6e5 | 1e6 | 2e6 |
|---|---|---|---|---|---|---|---|---|
| TDC (M3b guess) | F | F | F | F | F | 3 | 4 | – |
| TDC unset, same state | 5 | 5 | 6 | 10 | F | F | F | – |
| TDC, α_D = 0 | 5 | 5 | 7 | 11 | F | F | F | F |
| **TDC with the M3d guess** | F* | 2 | 2 | 2 | 2 | 3 | 4 | 8 |

(F = 26 iterations, not converged; F* = first probe after a restart, see M3d.5.) The control step (B2 itself) fails at 6.4e4–7.3e4 s: the same threshold. The Newton is *not* made hard by the large step (at 6e5–1e6 s it converges in 3–4 iterations).
2. **No TDC row is the cause** (`variants_dt3e4.txt`, probe Δt 3e4, 11 variants in one process): F_c off or → local MLT value, Q off, α_p = 0, α_t = 0 (each alone) all still fail (26); **α_D = 0 converges in 5, and every variant that has α_D = 0 converges in 5**
(incl. α_p = α_D = 0, Q off, w frozen at w_old; F_c → 0). The w unknown is not involved: with the w row replaced by w = w_old (`wfrz`) the failure is the same, the iterations never have a w entry as the largest weighted correction (8 failed
attempts, 200 Newton iterations: x_H⁺/x_H₂ at zone 37, F_r at zones 38 and 74; `largest_entries_failed_attempts.txt`), and |δw| ≤ 5e-3 `unom` per iteration with the weight 1/max(|w|, 0.01 c_s) of `tdc_dmax!`/`tdc_weight`, i.e. the scaling of
w in the convergence measure and in the trust-region measure is not limiting.
3. **First correction** (`first_iteration.txt`): the Newton direction of the first iteration from the start state is |δ₀|_W = 2.7 (TDC unset, TDC with α_D = 0) but **136** with α_D = 1: the largest weighted entry is the x_H⁺ row of zone 37, 143 at weight 0.24, i.e. δ ln x_H⁺ ≈ 6×10²,
and the global step fraction λ = 0.019 (0.49 without D_mix). The trust region and the Armijo search scale the *whole* vector by one λ, so E_r, T, F_r and everything else of the system creep at 2 % per iteration (median λ of the 8 failed attempts 0.125): the step is lost to a single log-variable entry.
4. **Where that entry comes from** (`xHp_row_terms.txt`, the terms of the x_H⁺ row of zone 37 at the first residual evaluations, Δt 3e4): the iterate the Newton starts from has x_H⁺ = 3.2e-7 against the old/converged 4.5e-6, because Newton's initial species are the **local solution of the zone's two
species rows alone** (`cpl_local!`, `TitanResidFast.cpl_local_all!`: neighbours' x held, 2×2 Newton in ln x), and that solve contained advection and chemistry but **not the turbulent mixing term** that the residual rows carry since M3b (`F_{k+1} − F_k`, D = α_D ℓ w). In zone 37 the mixing influx from zone 36
(x_H⁺ 0.355; D = ℓ w = 8.0e11 × 1.46e5 = 1.2e17 cm² s⁻¹ at interface 37) is 1.4e22 g s⁻¹ against storage −7.2e19, chemistry +3.9e19 and advection 9e17: 200× the largest other term of the row. Without the mixing the guess relaxes to 3e-7 (14× below the mixing-fed solution; [interpretation] recombination with t_rec ≪ Δt); in the
row the x_H⁺ source from the neighbour does not depend on the zone's own x, so Newton in ln x from below sees δ ln x ≈ ½((x_eq/x)² − 1) [derivation for P = k x²], hundreds, against the ln 14 = 2.6 that is needed. [interpretation] This is a property of the log variable with a source independent of the zone's own value (the same would hold for the advective inflow of a species from an ionised neighbour), exposed here because the
only large x-independent source into a cold, recombined zone is the D_mix flux across the cliff.
5. **No physics requires Δt ~ 1e5 s** (`timescales_step605.txt`, step 605, interfaces 33–37 carry w = 8.5e3 … 2.0e5 cm/s): turnover ℓ/a = 5.7e5 … 1.1e6 s; KW growth time ℓ/(a√sup) = 4.2e6 s (36), 2.4e6 s (37); w advection across the CV Δr/|u| = 6.6e6 s (36), 2.1e6 s (37), 1.3e6 s (38); ℓ/w ≥ 4.4e6 s;
diffusion Δr²/(α_t ℓ w) ≥ 8e8 s; radiative loss ℓ/(2Ua) = 2.5e6 s at the active interface 37 (≤ 700 s at 38–41, where w ≈ 0 and the damping of the linear w term is merely stiff, not a time scale of the solution). The smallest TDC time scale of the layer is 4e5 s, i.e. steps of 1e6 s are resolved and the 7e4 s steps are numerical, not physical.

### M3d.3 Fix [read; code `cf584c8`]

The local initial guess solves the same rows as the residual: `PopIIIMaterial.cpl_local!` (julia/PopIIIMaterial.jl:2213) adds `+ (c_k + c_k+1) x − c_k+1 x_{k+1} − c_k x_{k−1}` with the coefficients c_j = rmu_j ρ̄_j α_D ℓ_j w_j/Δr_c,j of `TitanResidFast.tdc_mix_coef` (julia/TitanResidFast.jl:970; w the iterate `WTN`,
ℓ, r, rmu of the last residual evaluation, ρ the predicted one; zero at the boundaries). The hooks of w (`tdc_newstep!`, `tdc_reset!`) now run before the species hooks (julia/popiii_m2c.jl:2014) so that the guess sees the old w also after a failed attempt. Nothing in the rows, the Jacobian, the weights or the converged solution changes (Newton's start only); switch `TDCMIXGUESS=0` = M3b guess (inventory, docs/POPIII_SWITCH_INVENTORY.md).
The two options the task names as candidates, w scaling and a missing derivative, are not involved (point 2, M3b.3(b)); the error estimate is discussed in M3d.5.

### M3d.4 Tests [test/run]

- **Probes** (table above, row "M3d guess"): 2 iterations at probe Δt 6e4 … 3e5, 3 at 6e5, 4 at 1e6, 8 at 2e6; first correction |δ₀|_W = 0.017 (against 136).
- **200 steps from the checkpoint** (`fix200`, `dt_history_fix200.txt`, `attempts_*.txt`; same start, same cores): over the first 40 steps ΣΔt = 1.06e7 s against 1.18e6 s (×9), median Δt 1.7e5 against 3.4e4, 3 Newton-failed attempts against 8 (7.5 against 20 per 100 steps); over 200 steps ΣΔt = 7.6e8 s = 24.0 yr (before: 0.037 yr in the 40 steps, ≈ 0.2 yr per 200 steps extrapolated; B2 production 2 yr per 400 steps, 60 failures per 350 steps), 22 failed attempts (11 per 100 steps),
median Δt 3.3e5 s, 25-step medians up to 9.4e6 s (steps 650–674) and 7.8e6 s (750–774) with maxima 3.7e7 s; wall 3.2e3 s for 200 steps (15.8 s/step against ≈ 25 s/step in B2).
- **Bit-identity, TDC unset** (`bitid_ic40.txt`): IC, zie1 settings, 40 steps, SNAPEVERY=5, HEAD 00326c0 against HEAD + fix: 10 snapshots, 300 fields, 0 differences.
- **NKCHECK:** no residual row, derivative, weight or measure was changed (initial guess and hook order only), so the J·v results of M3b.3(b) stand unchanged; not repeated.

### M3d.5 What is left, and what was not done [run; interpretation; open]

- **A second wall remains, and it is not in the TDC rows** (`second_wall_k616.txt`, `second_wall_k800.txt`). With the fix the Newton still fails at isolated Δt in the layer's cliff zone: 22 failed attempts in 200 steps, each retry (Δt ÷ 4) converges in 2 iterations, and the failing attempts show a **period-2 cycle**
at the cliff zone (x_H2 ±0.3 … 0.5 at the cap, T ±3 %, E_r ±0.2, x_H⁺ ±16 %, alternating every iteration, λ 0.37/0.79; `fix200/nwt.txt` attempts 36, 49). It is dt-non-monotonic: from the state before step 617 the stage Newton converges at probe Δt 2e5 (6 it) and 8e5 (5) but fails at 4e5; from the state of step 800 (cliff now at zones 38–39) it fails
at 1e5 and 1.6e6 and converges in 3, 3, 5 iterations at 2e5, 4e5, 8e5. With TDC unset from the same states: 12 it at 2e5 and failures at ≥ 4e5 (k616), 9/11/13 iterations at 1e5/2e5/4e5 and failures at ≥ 8e5 (k800), i.e. the TDC path is not worse than the production MLT path from the same state. Single TDC ingredients change which Δt fail
(k616: α_p = 0, Q off, w frozen, F_c → local MLT each converge at 4e5 in 10, 10, 8, 3 iterations, while `ap0,aD0,noQ,wfrz` and α_D = 0, α_t = 0 alone fail), so there is no coupling to blame and no systematic effect: the signature is that of the closure/chemistry bistability at the one-zone photosphere of M3c (Eddington-factor branch decisions, x_H₂ switching),
which the baseline passes with one ÷3–4 retry as well (M3c.3). It is not part of this fix [interpretation]. Its cost: two of the eight 25-step blocks (675–699, 775–799) still run at median Δt 9e4 and 5e4 s (the controller grows ×1.5 per step into the failing Δt ≈ 1e5 s window, fails, retreats ÷4), the other six at 1.6e5 … 9.4e6 s (`dt_history_fix200.txt`); the failing attempts of 775–799 are again x_H₂ of the cliff zone (99 of 100 first-entries at zone 39).
- **The error estimate.** [read, M3b.3(5)] w is still not in DTERR's estimate. With the fix Δt reaches 2.7e7 s while the layer's w is quasi-steady (ℓ/a = 6e5 … 1e6 s ≪ Δt): backward Euler (L-stable) follows the steady target, which moves only as T, ρ, x of the zones move, and those are in the estimate; between snapshots 5 steps apart w changed by ≤ 10 % at interface 37 (27 % at the weak interface 35 in the first window at Δt ≈ 1e6 s, ≤ 3 % at 620–635, Δt 4e4 … 2e5 s; `wchg_fix200.txt`); the layer moves outward with the shock (interfaces 33–37 at step 605, 38–39 at step 800; R* 60 → 129 R☉ over the 200 steps). I have no failure or accuracy loss that points to w; adding m·w to `combos` is therefore left open (not needed for this wall), but the quasi-steady regime means a w-only truncation error would first matter
at Δt ≳ ℓ/a with a *changing* source, which the T, x estimate bounds from the other side [interpretation, untested].
- **First step after a restart** has no residual workspace yet, so its guess carries no mixing (`tdc_mix_coef` returns 0 until a Float64 residual has been evaluated; v3 probe 1: 26 iterations; fix200 attempt 1: 16 iterations, converged). Harmless in the 200-step run; a guess from the old state's ℓ would need the equation of state at the start, not done.
- **D_mix magnitude** (α_D = 1: D = ℓ w, 1.2e17 cm² s⁻¹ at the cliff, mixing time of zone 37 5.6e6 s ≫ Δt) is a modelling input not tested here; with the fixed guess the α_D = 0 and α_D = 1 stage Newtons both converge at k800 (aD0: 9, α_D = 1: 3 iterations at 2e5).
- **Denormals:** in dissolved parts of the layer w decays exponentially to |w| ~ 1e-320 (`wmax`/`wchg`; also slightly negative values: fb(w, F) with F > 0 gives w ≈ 0⁻ is allowed, nothing is clamped); not a convergence problem, a possible slowness of the Float64 arithmetic [open, not measured].
- Not done: the complementary test of the same fix on the k115/k290/k650 pilots of M3b.3(d) (the k115 failure is the closure wall of M3c); a 200-step *unfixed* run (the before numbers are the 40-step restart `base` and the B2 lane).

### M3d.6 Conclusion [interpretation]

The Newton wall at Δt ≈ 7e4 s of the TDC path was a **starting-point defect, not a property of the TDC equations**: Newton's initial species were computed from rows without the turbulent mixing flux, and the mixing-fed, recombined cliff zone starts 14× too low in x_H⁺, which a log-variable Newton step answers with δ ln x ~ 6×10², and the
shared trust-region fraction collapses to 2 %. Putting the same term in the local guess (a Newton-start change only) lowers the first correction 136 → 0.017, converges the stage Newton in 2–4 iterations over probe Δt = 6e4 … 1e6 s, and moves the 200-step run from Δt ≈ 3e4 s (20 failed attempts per 100 steps) to a median 3.3e5 s with phases of 1e7 s (24 yr in 200 steps, 11 failed attempts per 100 steps).

## M4 — Bracket of the accretion history (22–24 Nov 2026)

### Set-up [run]

**Lanes.** All lanes start from replay checkpoint `rp_k290`: zie1, 5 yr after protostar formation, 0.29 M☉, 100 zones.
- Lane runner: `scripts/subconv/m4run.sh`, which runs segments of 400 steps until M* ≥ 34 M☉.
- Code: `wt_m4`, with the TitanCore of aace2a0 (as for every replay restart, M1/M3b).
  - B1 ran on 4c026bc. The TDC lanes ran on dcb9bb8, which includes M3d and predates the decide() upwind fix of 4636f84. That fix changes only the Jacobian, not the residual.
- Settings: zie1's run.env plus `CONVDIAG=1`.
- Comparison: `scripts/subconv/m4compare.py`, giving `docs/logs/subconv/m4/compare_final_zie1.txt` (with its .png files) and `compare_final_zie3.txt`.

**Why the 0-yr phase is not bracketed.** TDC from k115 meets the radiation-closure wall of M3c, which is not a convection effect.

| case | definition | reached |
|---|---|---|
| B0 | production (replay zie1rp; no convection in the layer, M1) | 33.8 M☉ (step 2400) |
| B1 | local MLT with the M2 buoyancy and thermodynamics (`MLTBUOY=1`) | 34.1 |
| B2 | TDC best estimate, `TDCPAR=1,1,1,1,0.01,1,0` (α = 1.5, α_p = 1, α_t = 0.01, α_D = 1) | 45.9 |
| B3a / B3b | B2 with α = 0.75 / 3.0 | 49.3 / 46.3 |
| B3c | B2 without overshoot diffusion (α_t = 0) | 45.8 |
| B4 | B2 without turbulent pressure (α_p = 0) | 45.8 |
| B5 | efficient-convection proxy: B2 with α = 10 (ℓ ≫ H_P, so ∇ → ∇_ad, and D_mix ∝ ℓw large) | 7.2 (still running) |
| B0z3 / B2z3 | B0 / B2 on zie3's 298-zone mesh, from zie3rp k775 (0.29 M☉) | 1.3 / 0.6 (running) |

B5 approximates "forced adiabatic + full mixing" with large α; it does not literally impose ∇ = ∇_ad [interpretation].

### Results [run, `compare_final_zie1.txt`; spread = (max − min)/B0 over B0–B5 at fixed M*]

**The layer itself** (cv_ scalars):
- In every convective case it carries the whole flux, F_conv/F_tot ≈ 1, up to about 2 M☉ (≈ 50 yr). It is intermittent after that and has mostly disappeared by 3–10 M☉. That is consistent with M1: the layer is unstable only early.
- Turbulent pressure behind the shock, P_t/P: 2–17 % in B2, 3–11 % in B3b (2α), ≤ 5 % in B3a, and 22–59 % in B5.

**Accretion:**
- Ṁ through 10⁴ R☉ is unchanged: spread ≤ 0.36 % in all cases, at every M*.
- Ṁ at the shock [P] varies by 4–10 %.

**Infall and entropy:**
- Peak infall speed u_min: spread ≤ 13 % without B5, and up to 23 % with B5 at 3–5 M☉ (faster infall onto a smaller star).
- Post-shock entropy proxy s_ps [P]: spread ≤ 6 %.

**The star:**
- R*: B1–B4 differ by at most one zone (the mesh has 29 % steps in r). B5 is two zones smaller at 3–5 M☉ (−40 %).
- T_c rises with convection efficiency: +2–5 % at 2–30 M☉ (B2 +5.4 % at 30 M☉; B5 +4.9 % at 5 M☉).
- Deuterium: L_D at 10–30 M☉ differs by factors of 2–4 between cases, and ⟨X_D⟩ inside R* by 2 % at 30 M☉. Burning starts at the same mass; its rate is sensitive to T_c.
- No case reaches pp or CNO onset (L > 1e-3 L_ph) by 46 M☉, so these onset masses are not bracketed.

**Photosphere:**
- B0's L_ph scatters by about 1 dex (rolling σ log L_ph ≈ 1.0–1.4) through 20 M☉, so the pointwise L_ph and T_eff differences are noise-dominated [P].
- The TDC cases are quieter at 10–20 M☉ (σ 0.015–0.6). Convection in the layer damps the photospheric fluctuations of the 100-zone runs [run; mechanism not isolated].

**Resolution** [run, `compare_final_zie3.txt`, 0.29–0.5 M☉ only]:
- On the 298-zone mesh, TDC also convects with F_conv/F = 1 and P_t/P 1.3–1.8 %. Its effect is smaller than at 100 zones: R* −1.5 %, Ṁ_shock +1.8 %, s_ps −0.3 %, T_c +0.04 %.
- The zie3 mesh itself (B0z3 vs B0) moves R* by more than the whole convection bracket does (46 vs 60 R☉ at 0.3 M☉).

### Conclusion [interpretation]

In 1D, convection in the sub-shock layer is real, and with the gate open it is efficient early on (≲ 2 M☉). It barely changes the accretion history:
- the mass inflow is unchanged;
- R* moves by about one zone;
- infall speed changes by ~10 %, and by ~20 % in the efficient limit;
- T_c rises by a few %, which shifts the deuterium-burning rate by factors of 2–4.

The effect at 100 zones is of the same order as, or smaller than, the effect of the resolution itself, and on the 298-zone mesh it is smaller still.

**What 1D cannot bracket:**
- corrugation of the shock by plumes;
- lateral transport;
- turbulent pressure fluctuations in the jump.

These need the multi-D calibration of M5. **Open:** B5 and the zie3 pair to 34 M☉ (running); the T10 resolution test at ~1000 effective zones; the 0-yr phase (blocked by the M3c closure wall).

## M5 — Plan for the multi-D calibration (written plan only; not started, needs Tom's agreement)

**Target.** A 2D wedge, later 3D, of the 111-yr state (rp_k655: R* = 278 R☉, 3.2 M☉). Radial range about 240–700 R☉: the unstable layer at 0.6–0.8 R*, the shock, and about 1.5 R* of inflow.

**Resolution and domain.**
- ≥ 20–30 zones per H_P in the layer. H_P/r ≈ 0.06–0.2 (M1), so Δr/r ≈ 3e-3: about 400 radial zones.
- 90° in θ with comparable aspect ratio, i.e. about 500 θ zones. Periodic or reflecting in θ.

**Physics.**
- The run's EOS: a PopIIIEqEOS table in the equilibrium limit, which M1 shows applies to the layer (t_chem ≪ 1/N).
- Gray flux-limited diffusion with the run's Rosseland/Planck means (G2/StellarOpacity tables).
- Gravity from the enclosed mass.

**Boundary conditions.**
- Inflow at the outer boundary from the 1D profile (ρ, u, T).
- Inner boundary at 240 R☉: a hydrostatic reflecting wall with the 1D luminosity.

**Calibration quantities:**
- the convective flux and the depth of the mixed layer, against α (A_s, A_d);
- the turbulent pressure in the post-shock state, against α_p;
- overshoot towards the shock, against α_t;
- the mixing of the ionisation fraction, against α_D;
- shock corrugation amplitude and SASI-like modes (1D cannot represent these);
- convergence at 2× resolution.

**Host code** (Tom to choose):
- (a) **NestedTurbulence** (Julia/KernelAbstractions, GPU, Float32 per the A6000 rule). It needs gray FLD and a tabulated EOS.
- (b) An external radiation-hydro code with FLD, e.g. Athena++ (implicit radiation) or a MESA-independent setup.

[estimate] Cost: with (a), 2D at 400×500 is about 2e5 cells; to run ~10 turnovers (≈ 10⁷ s) at an FLD-limited Δt is hours on the GPU.
