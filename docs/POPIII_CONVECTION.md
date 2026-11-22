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
