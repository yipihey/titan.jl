# Pop III sub-shock convection: M1 b–d diagnosis (22 Nov 2026)

Labels: [run] = measured in a run of this study; [read] = read from code or from a log; [estimate] = computed from the logged numbers with a stated approximation; [interpretation] = inference, not shown directly.
All numbers are in `docs/logs/subconv/m1diag/` (file names in brackets); scripts are `scripts/subconv/diag_*.{jl,sh,py}`. Units: k_B/baryon/yr for entropy rates, yr for times (t_chem in seconds), r/R* with R* the first interface from the centre with u < −2 km/s, as in `convdiag_snap`.

## 0. Replay versus production [run]

`diag_prodcmp.jl` compares every field of every common snapshot with `isequal` (bitwise). zie1 replay (`zie1rp`, commit aace2a0) versus production `zie1` (`snaps_000000001_2321108.jls`): 482 common snapshots (steps 1…2400, every 5), 14 460 fields, **0 differences** [step0_prodcmp.txt]. The replay is bit-identical to production.

## 1. Restarts with CONVDIAG=1 and bit-identity [run]

Setup: zie1 (zie3) `run.env` verbatim with SNAPEVERY=5, `CONVDIAG=1`, `RESTART=rp_kNNN`, MAXSTEPS=200 (100 for zie3), own CKDIR, 8 threads (`diag_epoch.sh`, `diag_chain.sh`).

**Finding that changes the setup: the code at c523e91 does not continue the replay bit-identically; the code at aace2a0 does.** [restart_code_version_test.txt, titancore_aace2a0_to_c523e91.diff]
- Restart from `rp_k115`, pure c523e91 (the pinned copy), 200 steps: 906 of 1200 fields differ at all 40 snapshots against the replay, already at the first snapshot (step 120: t differs at the 10⁻⁹ level, dt 4.4e5 s vs 1.4e5 s). The restarted run then fails Newton solves and advances with dt ≈ 3e2–1.7e3 s against the replay's 1e5–8e6 s, so the state diverges completely (max relative T difference 17 at step 315). The M1 test `rp650off` (also c523e91, restart k650) differs from the replay in 84 of 120 fields as well. The M1 on/off bit-identity test compared c523e91 with itself, so it did not show this.
- Restart from `rp_k115`, code aace2a0, CONVDIAG unset, 10 steps: 60 fields, 0 differences.
- Restart from `rp_k115`, c523e91 with only `julia/TitanCore.jl` taken from aace2a0 (a scratch worktree), CONVDIAG=1, 10 steps: 60 fields, 0 differences.
- [read] The only non-diagnostic solver difference between aace2a0 and c523e91 in `julia/` is the `predict!` else-branch of `TitanCore.jl` (commit b482bdc, 22 Nov 06:34, log-extrapolated widths and zone masses; `julia_diffstat_aace2a0_to_c523e91.txt`). b482bdc is an ancestor of c523e91 but not of aace2a0. Why a branch that `popiii_m2c.jl` leaves unused with `EXTRAP_FRAC=0.3` changes the trajectory was not investigated. [interpretation] Either the branch is reached after a restart (EXTRAP_FRAC left at Inf), or a second effect of the same file is involved; this needs a look before M2 builds on the logmass-based branch.
- All diagnostic runs below therefore use c523e91 with `TitanCore.jl` from aace2a0 (`wt_diag_b`, uncommitted scratch copy). With it:

| epoch | restart | snapshots compared | fields | differences | reference |
|---|---|---|---|---|---|
| zie1 k115 | step 115 | 40 (120…315) | 1200 | 0 | zie1rp [zie1_k115_bitid_vs_replay.txt] |
| zie1 k295 | step 295 | 40 (300…495) | 1200 | 0 | zie1rp [zie1_k295_bitid_vs_replay.txt] |
| zie1 k655 | step 655 | 40 (660…855) | 1200 | 0 | zie1rp [zie1_k655_bitid_vs_replay.txt] |
| zie1 k2395 | step 2395 | 1 (2400; the replay ended there) | 30 | 0 | zie1rp [zie1_k2395_bitid_vs_replay.txt] |
| zie1 k2395 | step 2395 | 40 (2400…2595) | 1200 | 0 | production zie1 [zie1_k2395_bitid_vs_production.txt] |
| zie3 k190 | step 190 | 20 (195…290) | 600 | 0 | zie3rp [zie3_k190_bitid_vs_replay.txt] |
| zie3 k775 | step 775 | 5 (780…800; the replay ended there) | 150 | 0 | zie3rp [zie3_k775_bitid_vs_replay.txt] |

The solver fields of the CONVDIAG restarts are bit-identical to the replay (and, for k2395, to production) wherever they overlap; the diagnostics therefore describe the production trajectory.

## 2. The layer, epoch by epoch

Layer = the contiguous run of interfaces just inside the shock with Mach < 0.5 and N²_L < 0 or N²_eq < 0 (`diag_layer.jl series`, `*_series.txt`; resolution in `resolution.txt`).

| epoch | zones | R* [R☉] | layer interfaces | layer width [H_P] | zones per H_P at the layer | max(−N²)t_dyn² | frozen potential flux (all gates open) |
|---|---|---|---|---|---|---|---|
| zie1 k115 (0.017 M☉) | 100 | 5.99–59.95 (grows 10×) | none at steps 120–125, then 1–3 | 2.5–4.2 | 0.6–0.9 | 0.01–1.6 | 0 (∇ < ∇_ad,fr) |
| zie1 k295 (0.29 M☉) | 100 | 60–100 | 1–3 | 3.2–4.2 | 0.5–0.9 | 0.2–4.4 | 0 at steps 300–415; 1e5 F_tot at steps ≥ 420 |
| zie1 k655 (3.2 M☉) | 100 | 278–359 | 2, then 1, then none after step 765 | 3.3–3.8 | 0.4–0.75 | 0.8–11 | 0.07–700 F_tot [zie1_k655_series.txt] |
| zie1 k2395 (34 M☉) | 100 | 215 | none | – | – | – | – |
| zie3 k190 (0.03 M☉) | 298 | 7.7–16.7 | 1–4 | 0.5–0.9 | 1.9–2.3 | 0.03–0.7 | 0 |
| zie3 k775 (0.29 M☉) | 298 | 46–51 | 1 | 0.75–0.8 | 1.3 | 0.9–1.4 | 0 |

- [run] zie1 k2395 (34 M☉): all subsonic interfaces outside the core have N² > 0 up to the shock (interface 41); the unstable subsonic interfaces counted in `zie1_k2395_gates.txt` (7 of them) are core interfaces at r/R* < 3e-4, T = 1.8e5 K, table branch, where MLT is active with ∇ ≈ ∇_ad (`zie1_k2395_s2400.txt`). The gas inside the shock is hot (T ≥ 5e4 K) and nearly ideal, with ∇_ad,fr = 0.25–0.28. There is no sub-shock layer at that epoch.
- [run] zie1 k115 at step 120–125 has no layer either: R* = 6–8 R☉, and the only unstable interfaces are at Mach 0.19–0.8 (the Mach-0.5 criterion excludes them).
- [read] The layer mass printed in `layer[]` is the mass between the two outer interfaces, so it is 0 for a one-interface layer.
- [run] The mesh is logarithmic with ratio 1.2915 (zie1) and 1.0890 (zie3) in this region, so the layer is 1–3 cells wide at both resolutions. zie1 has 0.4–0.9 cells per H_P, zie3 1.3–2.3 (`resolution.txt`).

## 3. (a) Is MLT active in the N² < 0 zones, and which gate blocks it? [run]

No. In every layer interface at every epoch and snapshot:
- branch = 1 (coupled network zone), `w_th = 0` (∇_ad,used = 1e30), `w_tau = 1`, `w_inf = 1`, F_conv = 0.
- Counts over all subsonic (Mach < 0.5) interfaces inside the shock with N²_L < 0 or N²_eq < 0 [`*_gates.txt`, last line]: of the unstable interfaces, the number with `w_tau < 1` is 0 and with `w_inf < 1` is 0 at all epochs. The only gate that closes is the thermodynamic one: `w_th = 0` in 71/93 (k115), 82/82 (k295), 38/74 (k655), 0/280 (k2395, core only), 57/84 (zie3 k190), 73/113 (zie3 k775).
- The unstable interfaces with `w_th = 1` (branch 2) are core interfaces (zie1 k115 j = 2 at r/R* = 3e-4; zie3 k190 j = 34 at r/R* = 0.011; both T ≈ 8e4 K, ∇ = ∇_ad,used within 1e-3, F_conv/F_tot = 1) [nonnetwork_unstable.txt]; they are not part of the sub-shock layer. In k655 the 38 blocked interfaces are the layer interfaces of steps 660–765; the other 36 unstable interfaces are two core interfaces (j = 7, 8, r/R* ≈ 1e-4, T = 1.25e5 K, branch 2, `w_th = 1`, F_conv/F_tot = 1) at steps 770–855, after the layer has disappeared (at step 770 N² > 0 in all subsonic interfaces below the shock) [zie1_k655_gates.txt, nonnetwork_unstable.txt].
- [run] Unblocked (all gates open, frozen ∇_ad, c_P, δ in place of the gate, `fpot_th`): F_pot/F_tot = 58 and 124 (step 655), 101 and 11.6 (step 700) in zie1 k655; 0 in k115, zie3 k190 and k775, and in k295 until step ≈ 415, because ∇ < ∇_ad,fr there (see (c)); 1.1–1.3e5 in k295 from step 420.
- The gate is therefore not an optical-depth or infall effect: at 0.017–3.2 M☉ the layer has τ ≥ 1e3 (`tau`, k655 interface 42: 1.3e3; k115 interface 33: 5e8) and the infall factor is 1. This confirms the M0/M1 reading, for four more epochs and for zie3.

## 4. (b) What sets the negative entropy gradient? [run]/[estimate]/[interpretation]

Definitions as in `diag_layer.jl`: zone rates c_i = −e_i μ m_H/(k_B T m_zone) yr in k_B/baryon/yr (m_zone = ρ ΔV per sr), heating positive; `sum` = Σ c_i; the inversion Δs across interface j is Δs = c_P N² Δr/(g_eff δ) in k_B/baryon with c_P,fr (so the Δs_eq column is a lower bound in magnitude); dΔs/dt = rate(zone j) − rate(zone j−1).

#### zie1, epoch k115 (0.017 M☉ at the checkpoint), step 200 (`zie1_k115_s200.txt`)

step 200 t=2.423805e+05 yr  zones 100  R*=35.94 Rsun (shock interface 34)  layer interfaces 32..33  layer[]=0.5995,0.7743,0.0137,0.6772,0,0,2,35.94

Interfaces of the layer (and the shock interface):

| j | r/R* | T[K] | Mach | N2_eq | N2_L | nab | nab_rad | nad_fr | nad_eq | w_tau | w_inf | w_th | kap | Fpot/F | dr/HP |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 32 | 0.5995 | 20999.5 | 0.013 | -1.805e-12 | 2.425e-12 | 0.180 | 0.304 | 0.3380 | 0.1442 | 1.00 | 1.00 | 0 | 5.884e+03 | 0.000e+00 | 1.227e+00 |
| 33 | 0.7743 | 16408.8 | 0.042 | -1.274e-12 | 2.174e-12 | 0.153 | 0.378 | 0.3624 | 0.1232 | 1.00 | 1.00 | 0 | 2.575e+03 | 0.000e+00 | 1.426e+00 |
| 34 | 1.0000 | 11860.9 | 0.284 | -1.325e-12 | 1.830e-12 | 0.153 | 2.246 | 0.3867 | 0.1043 | 1.00 | 1.00 | 0 | 3.391e+02 | 0.000e+00 | 1.763e+00 (shock) |

Energy-row terms of the zones adjacent to the layer, as entropy rates (k_B/baryon/yr, heating > 0; `sum` = the sum of the terms, `res/sum` = row residual / sum):

| k | r/R* | T | sum ds/dt | av | adv | rad | chem | res/sum | t_rad[yr] | t_settle[yr] |
|---|---|---|---|---|---|---|---|---|---|---|
| 31 | 0.5318 | 23405.8 | -7.15e-01 | 4.08e-04 | -7.16e-01 | 1.55e-04 | 2.06e-05 | -2.37e-02 | 6.86e+03 | 5.45e-01 |
| 32 | 0.6869 | 18593.1 | -4.34e+00 | 1.67e-02 | -4.35e+00 | 2.39e-04 | 7.40e-05 | -1.18e-02 | 4.27e+02 | 2.69e-01 |
| 33 | 0.8871 | 14224.5 | -2.15e+01 | 5.94e+00 | -2.74e+01 | 1.68e-03 | 4.28e-04 | -2.19e-03 | -1.63e+03 | 9.53e-02 |

Time scales and potential MLT in the layer interfaces:

| j | t_settle | t_N | t_chem[s] | sup_fr | t_turn_fr | Fpot/F | sup_eq | t_turn_eq | Feq_lb/F |
|---|---|---|---|---|---|---|---|---|---|
| 32 | 3.80e-01 | 2.03e-02 | 8.30e-11 | -0.158 | NaN | 0.00e+00 | 0.035 | 6.55e-02 | 3.95e+07 |
| 33 | 1.59e-01 | 2.15e-02 | 9.17e-10 | -0.210 | NaN | 0.00e+00 | 0.029 | 8.76e-02 | 2.14e+06 |

#### zie1, epoch k295 (0.29 M☉), step 400 (`zie1_k295_s400.txt`)

step 400 t=2.423924e+05 yr  zones 100  R*=100 Rsun (shock interface 38)  layer interfaces 36..37  layer[]=0.5995,0.7743,0.05468,2.178,0,0,2,100

Interfaces of the layer (and the shock interface):

| j | r/R* | T[K] | Mach | N2_eq | N2_L | nab | nab_rad | nad_fr | nad_eq | w_tau | w_inf | w_th | kap | Fpot/F | dr/HP |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 36 | 0.5995 | 21176.4 | 0.006 | -7.017e-13 | -1.976e-13 | 0.262 | 1.282 | 0.2940 | 0.1787 | 1.00 | 1.00 | 0 | 2.010e+03 | 0.000e+00 | 1.462e+00 |
| 37 | 0.7743 | 13492.8 | 0.056 | -7.292e-13 | 4.553e-13 | 0.158 | 2.509 | 0.3682 | 0.1058 | 1.00 | 1.00 | 0 | 8.200e+02 | 0.000e+00 | 2.011e+00 |
| 38 | 1.0000 | 7017.2 | 0.698 | -3.097e-14 | 1.160e-13 | 0.281 | 14.362 | 0.3989 | 0.2618 | 1.00 | 1.00 | 0 | 5.438e-02 | 0.000e+00 | 3.927e+00 (shock) |

Energy-row terms of the zones adjacent to the layer, as entropy rates (k_B/baryon/yr, heating > 0; `sum` = the sum of the terms, `res/sum` = row residual / sum):

| k | r/R* | T | sum ds/dt | av | adv | rad | chem | res/sum | t_rad[yr] | t_settle[yr] |
|---|---|---|---|---|---|---|---|---|---|---|
| 35 | 0.5318 | 25577.1 | -7.29e-02 | 2.04e-05 | -7.32e-02 | 2.73e-04 | 7.11e-06 | -3.17e-01 | 6.20e+03 | 3.21e+00 |
| 36 | 0.6869 | 16775.7 | -2.19e+00 | 3.19e-02 | -2.22e+00 | 4.42e-04 | 3.55e-05 | -2.32e-02 | 1.19e+03 | 1.07e+00 |
| 37 | 0.8871 | 10209.9 | 2.87e+00 | 2.16e+01 | -1.26e+01 | -6.10e+00 | 8.31e-04 | -6.05e+00 | -1.44e+00 | 1.45e-01 |

Time scales and potential MLT in the layer interfaces:

| j | t_settle | t_N | t_chem[s] | sup_fr | t_turn_fr | Fpot/F | sup_eq | t_turn_eq | Feq_lb/F |
|---|---|---|---|---|---|---|---|---|---|
| 36 | 1.87e+00 | 7.13e-02 | 4.91e-10 | -0.032 | NaN | 0.00e+00 | 0.084 | 1.07e-01 | 1.26e+06 |
| 37 | 2.67e-01 | 4.70e-02 | 7.12e-08 | -0.211 | NaN | 0.00e+00 | 0.052 | 1.26e-01 | 1.90e+04 |

#### zie1, epoch k655 (3.2 M☉), step 700 (`zie1_k655_s700.txt`)

step 700 t=2.425643e+05 yr  zones 100  R*=359.4 Rsun (shock interface 43)  layer interfaces 41..42  layer[]=0.5995,0.7743,0.6797,8.214,0,100.9,2,359.4

Interfaces of the layer (and the shock interface):

| j | r/R* | T[K] | Mach | N2_eq | N2_L | nab | nab_rad | nad_fr | nad_eq | w_tau | w_inf | w_th | kap | Fpot/F | dr/HP |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 41 | 0.5995 | 45327.5 | 0.001 | -2.115e-13 | -2.147e-13 | 0.485 | 0.643 | 0.3103 | 0.2962 | 1.00 | 1.00 | 0 | 1.015e+00 | 1.009e+02 | 1.388e+00 |
| 42 | 0.7743 | 17267.4 | 0.013 | -2.141e-13 | -5.054e-13 | 0.316 | 24326.680 | 0.1533 | 0.2046 | 1.00 | 1.00 | 0 | 3.368e+03 | 1.155e+01 | 2.435e+00 |
| 43 | 1.0000 | 4825.6 | 1.140 | 2.073e-14 | 2.109e-14 | 0.029 | 0.029 | 0.3980 | 0.3923 | 1.00 | 1.00 | 0 | 3.879e-04 | 0.000e+00 | 1.857e+01 (shock) |

Energy-row terms of the zones adjacent to the layer, as entropy rates (k_B/baryon/yr, heating > 0; `sum` = the sum of the terms, `res/sum` = row residual / sum):

| k | r/R* | T | sum ds/dt | av | adv | rad | chem | res/sum | t_rad[yr] | t_settle[yr] |
|---|---|---|---|---|---|---|---|---|---|---|
| 40 | 0.5318 | 61830.3 | -6.77e-03 | 8.25e-11 | -9.67e-03 | 2.90e-03 | 1.60e-06 | 1.01e+00 | 7.04e+02 | 4.18e+01 |
| 41 | 0.6869 | 28824.7 | -1.35e-01 | 7.33e-05 | -1.28e-01 | -6.72e-03 | 7.57e-06 | -1.20e-01 | 2.81e+02 | 1.90e+01 |
| 42 | 0.8871 | 5710.2 | -8.39e+00 | 4.12e+01 | -2.67e+00 | -4.69e+01 | 1.26e-03 | 2.95e-02 | -5.15e-01 | 1.51e+00 |

Time scales and potential MLT in the layer interfaces:

| j | t_settle | t_N | t_chem[s] | sup_fr | t_turn_fr | Fpot/F | sup_eq | t_turn_eq | Feq_lb/F |
|---|---|---|---|---|---|---|---|---|---|
| 41 | 3.49e+01 | 6.84e-02 | 1.74e-10 | 0.175 | 2.12e-01 | 1.01e+02 | 0.189 | 2.02e-01 | 1.16e+02 |
| 42 | 3.01e+00 | 4.46e-02 | 1.65e-04 | 0.163 | 2.12e-01 | 1.16e+01 | 0.111 | 2.53e-01 | 6.63e+00 |


Selected rows of `*_summary.txt` (per layer interface: Δs, dΔs/dt and its terms, ∇_rad/∇, κ at the interface and its maximum):

| epoch / step | j | Δs_L | Δs_eq (lower bound) | dΔs/dt | of which av | adv | rad | ∇_rad/∇ | κ_j [cm²/g] | max κ in layer ±3 |
|---|---|---|---|---|---|---|---|---|---|---|
| zie1 k655, 700 | 41 | −1.10 | −1.07 | −0.13 | 7e-5 | −0.12 | −9.6e-3 | 1.3 | 1.0 | 3.4e3 (j 42) |
| zie1 k655, 700 | 42 | −9.14 | −3.77 | −8.26 | +41.2 | −2.54 | −46.9 | 7.7e4 | 3.4e3 | 3.4e3 |
| zie1 k655, 660 | 41 | −7.48 | −2.97 | −2.00 | +29.1 | −1.24 | −29.8 | 5.6e3 | 2.6e3 | 2.6e3 |
| zie1 k295, 400 | 36 | −0.19 | −0.47 | −2.12 | +0.032 | −2.15 | +1.7e-4 | 4.9 | 2.0e3 | 2.0e3 |
| zie1 k295, 400 | 37 | +0.59 | −0.37 | +5.06 | +21.5 | −10.4 | −6.1 | 15.9 | 8.2e2 | 2.0e3 (j 36) |
| zie1 k115, 200 | 32 | +0.35 | −0.13 | −3.62 | +0.016 | −3.64 | +8e-5 | 1.7 | 5.9e3 | 1.45e4 (j 29) |
| zie1 k115, 200 | 33 | +0.52 | −0.13 | −17.1 | +5.9 | −23.0 | +1.4e-3 | 2.5 | 2.6e3 | 1.45e4 (j 29) |
| zie3 k190, 250 | 86 | +0.29 | −0.03 | −60.8 | +10.7 | −71.6 | +6e-3 | 1.4 | 3.7e3 | 1.1e4 (j 83) |
| zie3 k775, 800 | 100 | +0.54 | −0.07 | −44.1 | +41.8 | −85.9 | +1.9e-2 | 2.3 | 6.0e2 | 1.9e3 (j 97) |

(The `res/sum` column of the energy tables shows that the row residual of the converged step is 0.1–1 of the sum in the layer's interior zones at k655 and k295 (|sum| ≲ 1e-2); only rates ≳ 1 k_B/baryon/yr, i.e. the shock-adjacent zones, are well above the solver residual.)

Quantitative reading:
1. [run] In the interior zones of the layer the heating/cooling sources are negligible: |e_rad| is 1e-3–1e-5 of |e_adv| (k655 zone 41: rad −6.7e-3, adv −0.128; k295 zone 36: rad 4e-4, adv −2.22; k115 zone 32: rad 2.4e-4, adv −4.35; zie3 k190 zone 85: rad 1.5e-3, adv −28.8). The radiative cooling time of a layer zone, e·m_zone/|e_rad|, is 2e2–7e3 yr there (|t_rad| in the tables), against t_settle = H_P/|u| of 0.03–40 yr. The interior of the layer is adiabatic on the infall time; radiation does not change its entropy.
2. [run] Where the sources are not negligible is the zone adjacent to the shock (zone j_out of the layer, and the zones beyond). At k655 (step 700) the viscous heating +41.2 and the radiative cooling −46.9 nearly cancel and leave a net source −5.7; adv adds −2.7: sum −8.4 k_B/baryon/yr, with t_settle = 1.5 yr. The product, 8.4 × 1.5 ≈ 13 k_B/baryon, is of the order of the inversion Δs_L = −9.1 at interface 42 [estimate]. At k295 (step 400, zone 37) and zie3 (zone 100) the net source is a heating (av +21.6 vs rad −6.1; av +42.3 vs rad 0.02), at k115 zone 33 av +5.9 with no radiative term.
3. [run] The adv term of each layer zone is the largest term everywhere in the layer and always negative (the zone loses entropy). For an inflow through an Eulerian mesh it equals (ṁ/m)(s_up − s_zone), i.e. it measures the inversion itself, with t_adv ≈ t_settle × (dr/H_P) = 0.1–1 yr; it is [interpretation] a transport term (the profile is carried inward), not a source.
4. [run] Opacity: κ rises outward by 3 decades across the layer at k655 (1.0 → 3.4e3 cm²/g from interface 41 to 42, the H ionization edge at T = 1.7e4 K) and has an interior maximum inside the layer at the other epochs (k115: 1.45e4 at j = 29; k295: 2.0e3 at j = 36). ∇_rad (from the flux actually present and the κ MLT uses) exceeds ∇ in the layer: 1.3 and 7.7e4 (k655, step 700), 4.9 and 15.9 (k295), 1.7 and 2.5 (k115), 1.4 (zie3 k190), 2.3 (zie3 k775). The radiative flux profile L_rad(j) = 4πr²F_rad at k655, step 700 is non-monotonic: 2.8e3, 1.6e3, 2.6e3, 2.8e4 L☉ at j = 40, 41, 42, 43 [read, `zie1_k655_s700_radbudget.txt`], a dip at the opacity step. The zone sums of the energy-row term e_rad (4π × per-sr value) are +3.3e2, −2.7e2, +2.9e2, +2.6e4 L☉ (net cooling positive) in zones 39–42: in the shock-adjacent zone 42 this equals the flux difference L(43) − L(42) = 2.6e4 L☉, but in the two interior zones the e_rad terms (≈ ±3e2 L☉) are about a factor 3–4 smaller than the flux differences (−1.2e3 and +1.0e3 L☉). [interpretation] The mismatch (frame of `frad` versus the lab-frame radiation row) was not resolved; either way the interior radiative exchange is ≲ 3e2 L☉ against ~5e3 L☉-equivalent for the advective term of zone 41, so the conclusion below does not depend on it.

Conclusion (b) [interpretation]:
- The data do not support the "radiative flux blocked by the ionization opacity, with convergence heating inside" mechanism as the cause of the inversion: ∇_rad ≫ ∇ is present (the flux would need a steeper gradient than the layer has), but the radiative energy exchange of the layer interior is 1e-3–1e-5 of the advective term, and t_rad ≫ t_settle.
- The data also do not support "the outer part of the layer is radiatively cooled faster than the inner part" in general: only at k655 is the net source of the shock-adjacent zone negative (rad −46.9 exceeding av +41.2); at k295, k115 and in zie3 it is a net heating.
- What the data show is that the layer interior is adiabatic and carries an entropy profile with s decreasing outward across 1–3 cells: |Δs| = 0.03–3.8 k_B/baryon per interface at equilibrium chemistry (lower bounds, c_P,fr) and up to 9 (k655, frozen), which is transported inward at the infall speed; the entropy of a parcel is then set where it passes the shock and the shock-adjacent zones, where av and rad are of order 10–200 and nearly cancel. A decreasing entropy of the gas leaving the shock zone with time (freshly accreted gas arriving with lower entropy than earlier gas, now deeper) would produce exactly this profile. This was not tested: it needs the entropy of fixed mass coordinates against time (Lagrangian tracking over the 200 snapshots), which I propose as the next diagnostic. A decisive test of the "no local source" statement is an M2 run in which the layer is convective: if the inversion is inherited from the shock zone, convection removes it at once; if the shock zone keeps producing it, it reappears in about t_settle.

## 5. (c) Schwarzschild versus Ledoux, frozen versus equilibrium chemistry [run]

N²_eq uses ∇_ad, δ of the equilibrium table EOS; N²_L the frozen network/CK EOS plus the μ term (`nab_L`). Layer interfaces from the tables above:

| epoch / step | interfaces | N²_eq < 0 | N²_L < 0 | survives frozen chemistry? |
|---|---|---|---|---|
| zie1 k115, 200 | 32, 33 | both (−1.8e-12, −1.3e-12) | neither (+2.4e-12, +2.2e-12) | no |
| zie1 k295, 400 | 36, 37 | both (−7.0e-13, −7.3e-13) | 36 only (−2.0e-13; 37: +4.6e-13) | one of two |
| zie1 k655, 700 | 41, 42 | both (−2.1e-13, −2.1e-13) | both (−2.1e-13, −5.1e-13) | yes |
| zie3 k190, 250 | 86 | −4.3e-12 | +1.9e-11 | no |
| zie3 k775, 800 | 98, 100, 101 (zig-zag in sign) | −6.3e-13, −1.7e-12, −4.5e-12 | none (+3.1e-12, +5.6e-12, +4.5e-12) | no |

- [run] N²_S ≈ N²_L at k655 (interface 42: −3.0e-13 vs −5.1e-13: the μ gradient adds −2.0e-13), so the Ledoux correction is secondary.
- [run] The two limits differ in the layer because ∇_ad,fr = 0.30–0.39 while ∇_ad,eq = 0.09–0.18 (H ionization): ∇ = 0.13–0.26 lies between them at k115, k295, zie3 (equilibrium-unstable, frozen-stable), and above both at k655 (∇ = 0.485, 0.316 vs ∇_ad,fr = 0.310, 0.153).
- [run] t_chem is 8e-11–1.7e-4 s in the layer interfaces against t_N = 0.007–0.07 yr (2e5–2e6 s) and t_turn ≳ 0.05 yr. [estimate] Thermodynamics at equilibrium composition (f = t_turn/(t_turn + t_chem) ≈ 1 in the M2 blob interpolation) is the relevant limit; the frozen N²_L is the wrong one for the stability of the layer, and it only matters at k655, where both agree. [interpretation] The layer is unstable at every epoch from 0.017 M☉ to 3.2 M☉ in the sense relevant for convection; at the three low-mass cases the instability exists only because the ionization changes the adiabatic gradient.

## 6. (d) zie3 (298 zones) at the same mass [run]

#### zie3 (298 zones), epoch k190 (0.03 M☉), step 250 (`zie3_k190_s250.txt`)

step 250 t=2.361924e+05 yr  zones 298  R*=14.07 Rsun (shock interface 87)  layer interfaces 86..86  layer[]=0.9183,0.9183,0,0.3853,0,0,1,14.07

Interfaces of the layer (and the shock interface):

| j | r/R* | T[K] | Mach | N2_eq | N2_L | nab | nab_rad | nad_fr | nad_eq | w_tau | w_inf | w_th | kap | Fpot/F | dr/HP |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 86 | 0.9183 | 17116.7 | 0.097 | -4.313e-12 | 1.923e-11 | 0.135 | 0.194 | 0.3650 | 0.1188 | 1.00 | 1.00 | 0 | 3.718e+03 | 0.000e+00 | 5.162e-01 |
| 87 | 1.0000 | 15079.1 | 0.234 | -1.714e-11 | 4.188e-11 | 0.140 | 0.436 | 0.3720 | 0.1125 | 1.00 | 1.00 | 0 | 1.946e+03 | 0.000e+00 | 5.574e-01 (shock) |

Energy-row terms of the zones adjacent to the layer, as entropy rates (k_B/baryon/yr, heating > 0; `sum` = the sum of the terms, `res/sum` = row residual / sum):

| k | r/R* | T | sum ds/dt | av | adv | rad | chem | res/sum | t_rad[yr] | t_settle[yr] |
|---|---|---|---|---|---|---|---|---|---|---|
| 85 | 0.8807 | 17852.6 | -2.87e+01 | 1.04e-01 | -2.88e+01 | 1.47e-03 | 1.15e-04 | -2.91e-03 | -7.40e+02 | 3.37e-02 |
| 86 | 0.9591 | 16380.8 | -8.96e+01 | 1.08e+01 | -1.00e+02 | 7.60e-03 | 2.17e-04 | 5.54e-05 | -2.42e+02 | 2.13e-02 |

Time scales and potential MLT in the layer interfaces:

| j | t_settle | t_N | t_chem[s] | sup_fr | t_turn_fr | Fpot/F | sup_eq | t_turn_eq | Feq_lb/F |
|---|---|---|---|---|---|---|---|---|---|
| 86 | 2.92e-02 | 7.22e-03 | 1.61e-10 | -0.230 | NaN | 0.00e+00 | 0.016 | 5.24e-02 | 4.60e+06 |

#### zie3 (298 zones), epoch k775 (0.29 M☉), step 800 (`zie3_k775_s800.txt`)

step 800 t=2.361965e+05 yr  zones 298  R*=46.42 Rsun (shock interface 101)  layer interfaces 100..100  layer[]=0.9183,0.9183,0,1.343,0,0,1,46.42

Interfaces of the layer (and the shock interface):

| j | r/R* | T[K] | Mach | N2_eq | N2_L | nab | nab_rad | nad_fr | nad_eq | w_tau | w_inf | w_th | kap | Fpot/F | dr/HP |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 0.9183 | 13708.1 | 0.128 | -1.732e-12 | 5.573e-12 | 0.126 | 0.289 | 0.3829 | 0.1052 | 1.00 | 1.00 | 0 | 6.030e+02 | 0.000e+00 | 7.524e-01 |
| 101 | 1.0000 | 11217.6 | 0.494 | -4.504e-12 | 4.482e-12 | 0.143 | 0.517 | 0.3910 | 0.0925 | 1.00 | 1.00 | 0 | 1.491e+02 | 0.000e+00 | 9.286e-01 (shock) |

Energy-row terms of the zones adjacent to the layer, as entropy rates (k_B/baryon/yr, heating > 0; `sum` = the sum of the terms, `res/sum` = row residual / sum):

| k | r/R* | T | sum ds/dt | av | adv | rad | chem | res/sum | t_rad[yr] | t_settle[yr] |
|---|---|---|---|---|---|---|---|---|---|---|
| 99 | 0.8807 | 14806.0 | -1.66e+01 | 4.81e-01 | -1.71e+01 | 3.62e-03 | 1.07e-05 | 5.19e-04 | -2.28e+02 | 8.63e-02 |
| 100 | 0.9591 | 12610.3 | -6.07e+01 | 4.23e+01 | -1.03e+02 | 2.28e-02 | 3.85e-05 | 7.29e-04 | -1.35e+02 | 3.50e-02 |

Time scales and potential MLT in the layer interfaces:

| j | t_settle | t_N | t_chem[s] | sup_fr | t_turn_fr | Fpot/F | sup_eq | t_turn_eq | Feq_lb/F |
|---|---|---|---|---|---|---|---|---|---|
| 100 | 5.49e-02 | 1.34e-02 | 6.00e-09 | -0.257 | NaN | 0.00e+00 | 0.021 | 1.09e-01 | 1.16e+04 |


- Same gate: branch 1, `w_th = 0`, `w_tau = w_inf = 1`, F_conv = 0, in all 20 snapshots of both zie3 epochs [zie3_k190_gates.txt, zie3_k775_gates.txt].
- The instability is present at the finer resolution with the same magnitude: N²_eq = −4.3e-12 (0.03 M☉) and −1.7e-12 (0.29 M☉), max(−N²)t_dyn² = 0.03–0.7 and 0.9–1.4; equilibrium-unstable and frozen-stable at both, as in zie1 k115/k295. So it is not an artefact of zie1's 0.4–0.9 cells per H_P.
- Resolution: zie3 has 1.3–2.3 zones per H_P in the layer (1.3 at 0.29 M☉), zie1 0.4–0.9 (`resolution.txt`: 0.41–0.86 at k115–k655). The layer is one to four cells wide, 0.5–0.9 H_P in zie3 and 2.5–4.2 H_P in zie1; at 2 zones per H_P the sign of N²_eq alternates from interface to interface (zie3 k775: 97 +, 98 −, 99 +, 100 −, 101 −), so the resolved structure is a zig-zag of the same amplitude, not a smooth inversion. [interpretation] The layer is not resolved at either resolution; its depth is a few cells.
- Time scales in zie3: t_turn (equilibrium potential) is 0.05–0.11 yr and t_settle at the outer interfaces 0.03 (k190, j = 86) and 0.055 (k775, j = 100), so t_turn/t_settle = 1.8–2.0 at the outermost layer interface, 0.65 at k775 j = 98; in zie1 k295/k655 it is 0.006–0.5. [estimate] In zie3 the gas crosses the outer layer interfaces before one turnover of the potential MLT; the infall gate would matter there, while in zie1's layer convection could turn over many times per crossing.

## 7. Time scales in the layer (all from the tables above) [run]/[estimate]

| epoch / step | j | t_settle | t_N | t_chem [s] | t_turn (eq., potential) | t_turn/t_settle | F_pot/F_tot frozen | F_eq,lb/F_tot |
|---|---|---|---|---|---|---|---|---|
| zie1 k115, 200 | 32 / 33 | 0.38 / 0.16 | 0.020 / 0.021 | 8e-11 / 9e-10 | 0.066 / 0.088 | 0.17 / 0.55 | 0 | 4e7 / 2e6 |
| zie1 k295, 400 | 36 / 37 | 1.87 / 0.27 | 0.071 / 0.047 | 5e-10 / 7e-8 | 0.107 / 0.126 | 0.06 / 0.47 | 0 | 1.3e6 / 1.9e4 |
| zie1 k655, 700 | 41 / 42 | 34.9 / 3.0 | 0.068 / 0.045 | 1.7e-10 / 1.7e-4 | 0.20 / 0.25 (frozen 0.21) | 0.006 / 0.08 | 101 / 11.6 | 116 / 6.6 |
| zie3 k190, 250 | 86 | 0.029 | 0.0072 | 1.6e-10 | 0.052 | 1.8 | 0 | 4.6e6 |
| zie3 k775, 800 | 100 | 0.055 | 0.013 | 6e-9 | 0.109 | 2.0 | 0 | 1.2e4 |

(t in yr. t_turn = ℓ/v with the KW velocity v = ℓ √(gδ/8H_P) ζ, ζ = √(sup) in the limit U → 0 used for the equilibrium column — U ≤ 1.6e-3 ≪ √sup = 0.2–0.4 at these interfaces; F_eq,lb uses c_P,fr and is a lower bound since c_P,eq ≥ c_P,fr. [read] `t_turn`/`vpot` of the `cd` fields are 0/NaN in the network zones because ∇_ad,used = 1e30; they are recomputed here with the frozen and equilibrium thermodynamics.)

- t_chem ≪ t_N ≲ t_turn throughout.
- t_turn is 1.5–8 t_N, and t_turn/t_settle < 1 except in zie3. With all gates open the potential MLT flux exceeds the total flux by a factor 7 to 1e7 and the MLT velocities are 2–6e5 cm/s.

## 8. Summary

- [run] (a) MLT is off in every sub-shock unstable interface from 0.017 to 3.2 M☉ (zie1) and at 0.03/0.29 M☉ (zie3): `w_th = 0` (branch 1, ∇_ad,used = 1e30) is the only closed gate; `w_tau = w_inf = 1`. At 34 M☉ there is no sub-shock layer.
- [run] (c) The layer is Schwarzschild-unstable with equilibrium chemistry at all epochs (N²_eq < 0, t_chem ≪ t_N); with frozen chemistry it is stable at 0.017 and 0.29 M☉ (both resolutions) and unstable at 3.2 M☉.
- [interpretation] (b) The interior of the layer is adiabatic on the infall time and radiation cannot transport entropy there (t_rad = 2e2–7e3 yr against t_settle 0.03–40 yr); the inversion is carried from the shock zone, not built by radiative exchange or flux convergence inside the layer. The sign of the net source of the shock-adjacent zone is epoch-dependent. The test that would settle the origin is the Lagrangian entropy history.
- [run] (d) zie3 shows the same instability at the same magnitude; the layer is 0.5–0.9 H_P wide at 1.3–2.3 zones per H_P (zie1: 0.4–0.9 zones per H_P).
- [run] Process finding: the pinned c523e91 does not continue the replay; `TitanCore.jl` from aace2a0 does (section 1).

## Log index

`step0_prodcmp.txt`; `restart_code_version_test.txt`; `titancore_aace2a0_to_c523e91.diff`; `julia_diffstat_aace2a0_to_c523e91.txt`; per epoch E in {zie1_k115, zie1_k295, zie1_k655, zie1_k2395, zie3_k190, zie3_k775}: `E_bitid_vs_replay.txt`, `E_series.txt`, `E_gates.txt`, `E_summary.txt`, `E_s<step>.txt` (interface, energy and time-scale tables), `zie1_k2395_bitid_vs_production.txt`, `resolution.txt`, `nonnetwork_unstable.txt`. Scripts: `diag_prodcmp.jl`, `diag_epoch.sh`, `diag_chain.sh`, `diag_rtest.sh`, `diag_analyze.sh`, `diag_layer.jl`, `diag_resol.jl`, `diag_nonnet.jl`, `diag_cmpstate.jl`, `diag_md.py`. Run directories: `/zpool/nvme/data/tabel_scratch/titan_popiii/subconv/tests/m1diag/`.
