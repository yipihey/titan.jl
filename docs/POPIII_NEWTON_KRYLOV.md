# Newton-Krylov: the Newton iteration made exact against everything `update!` refreshes

Branch `nkrylov` (from popiii-s6 a4837c7 + the `nwdiag` diagnostics). Code: `julia/NewtonKrylov.jl`; hook `TitanDriver.KRYLOV[]` / `TitanDriver.newton_linear!()`
(used by `TitanDriver.step!` and the three Newton solves of `Trbdf2.newton0!`); switch `NKRYLOV` in `julia/popiii_m2c.jl`. All runs below:
ERLAB=1, 100 zones, the z100mm settings, restart copies of the production checkpoints, `/zpool/nvme/data/tabel_scratch/titan_popiii/runs/t_nkNN`,
notes `runs/t_nk_notes.md`.

## 0. Switch

`NKRYLOV=1[,<eta>[,<mmax>[,<eps>[,<central>]]]]` (provisional, **default off = bit-identical**, T5).
`eta` linear tolerance relative to the band correction (1e-3), `mmax` Krylov dimension (30), `eps` perturbation size of the finite-difference product (1e-5, in the
weighted norm), `central` 1 (default) central / 0 one-sided difference (test parameters). With `NWPROBE` the Krylov solve is switched on at the probe state only.
`NKCHECK=<iter>` (diagnostic, with `NWPROBE`): matrix-vector check at Newton iteration `<iter>` of the first attempt (T1), then exit.

Counters: the DTERR trace line gets ` | NK <Krylov its per Newton iteration>, <mv per step>, nonconv, fallback`; the run summary prints
`NKRYLOV ...: N Newton linear solves, M matrix-vector products, nonconvergences, fallbacks, state-drift events (max), seconds`; with `NWTRACE` every Newton iteration
gets a line `NKR att .. it .. krylov its .. converged .. rel residual .. |delta0|_W .. |delta|_W .. drift ..`.

## 1. Why (the nwdiag evidence, `runs/t_nw_notes.md`)

The production run stalls at a Newton wall: linear convergence with a contraction that grows with Δt (0.35 at step Δt 3e3 s, 0.81 at 1.5e4, 0.93 at 4e4); the 25-iteration
cap is reached at Δt ≈ 1.2e4 s, the median failed Δt of the production run (1.21e4). The closure quantities `eddfac!` writes in every `update!` (Eddington factor f = K/J,
then h = `HEDD`, the flux-tie weight `TIEW`) enter the residual as plain Float64, so ∂f/∂x is missing in the assembled block-band Jacobian.

| measurement (nwdiag, step-41 state of t_nw06/07/25/26/29) | result |
|---|---|
| probe iterations, step Δt 3e3 / 6e3 / 1e4 / 1.5e4 / 2e4 / 4e4 s | 9 / 14 / 21 / 31 / 41 / not converged in 60 (contraction 0.35 ... 0.93) |
| full FD (everything `update!` refreshes) against the assembled Jacobian, F_r row, column r | up to 1.0 × diagonal; entries at offsets −2, −3 where the assembled matrix is structurally 0 |
| same with f frozen / h frozen / tie frozen / all closure frozen | 0.27 / 0.79 / 0.93 / gone (col x_H+ 0.011 remains) |
| `NWFREEZE=f,g,h,hx,fx,gx,tie` (closure held over the Newton): iterations at Δt 3e3 ... 1e6 | 13, 13, 13, 13, 12, 11 |
| `NWFREEZE=f` only: iterations at 1.5e4, 4e4, 1e6 | 12, 10, 10 (unfrozen 31, >60, −) |
| `CPLJLAG=0` | 10, 16, 32, >60: not the cause |

Holding the closure fixed over the solve is a lagged input, not acceptable as a method. The Krylov solve below makes the Newton iteration exact instead.

## 2. Method

Unknowns are the code's correction vector δ (`rhs`; logarithmic variables where the code uses them, u and F scaled by `unom`/`frnom`, species in coupled zones). R(x) is the
residual after a full `update!` at x. J_b is the assembled band Jacobian, δ0 = −J_b⁻¹R(x) the correction `matslv` has just written. `krylov!` (called by `newton_linear!` right
after `matslv`) solves, by left-preconditioned GMRES,

    A δ = P δ0,    A v = P J_b⁻¹ (J v),    J v ≈ [R(x ⊕ εv) − R(x ⊖ εv)] / 2ε,

then the unchanged `maxdel!`, `update!` follow (clamps, convergence measure, closing `lboos` iteration, ERESTOL untouched).

* **x ⊕ εv**: εv is written into `rhs` and `TD.r(:update)` is called, i.e. the Newton update through its own code path (EOS, opacity, `eddfac!`, tie weights, species), then a
  fresh `context()` and the Float64 `residual!` at zero correction. The state is then restored from an exact snapshot (every array of the Julia state `TF.JSTATE` except the
  Jacobian and `rhs`, `HEDD`, `HXEDD`, `FXEDD`, `TIEW`, `GXEDD`, the species iterate) and the residual evaluated at the restored state is compared with R(x): **the drift is
  0.0 exactly** (`drift` in the NKR lines; T2).
* **P**: the entries `update!` applies are the unknowns: interface unknowns r, m, u, F of zones ngrs..ngre+1, zone unknowns ρ, T, E_r of ngrs..ngre, species in coupled zones.
  The ghost entries are slaved to them by the boundary conditions inside `update!`; P zeroes them. With an exact band Jacobian A = P and GMRES ends after one product with δ = Pδ0.
  The three Dirichlet rows of the boundary interfaces (u at ngrs and ngre+1, F at ngrs) have R = δ (the correction itself, not a function of the state): J v is 0 there, and so
  is every Krylov vector's entry, since δ0 has 0 there.
* **Branch decisions** (`decide`: limiter sign, upwind direction, viscosity switch): the assembled Jacobian differentiates the branches the base state takes (the Dual's value
  decides). The Float evaluations record them at R(x) (`TitanResidFast.FREEZE = 1`) and replay them in the perturbed evaluations (`FREEZE = 2`). Other regime changes inside the
  residual (clamps, `max`/`min` in the closure, the slope limiter inside `eddfac!`) are not under `decide`; the FD sees them only if a perturbation of 1e-5 crosses one.
* **GMRES**: no restart, Krylov dimension `mmax` = 30, modified Gram-Schmidt, Givens rotations. Norm: the weighted 2-norm ‖W·‖₂ with the weights of the DTERR error norm
  `Trbdf2.wnorm` (ρ, T relative; E_r against E + efl aT⁴; species against x + xfl; u against max(|u|, c_s); m against the zone mass; r against the zone width on a moving mesh;
  F against max(|F|, ffl cE)). `maxdel!`'s own measure (max of |δr|, |δm|, |δρ|, |δT|, |δE| and the weighted species) ignores u and F; the unweighted 2-norm is useless
  (the opaque-interior F corrections are relative changes of 1e2 of a flux 1e-4 cE). Stop when the preconditioned residual is ≤ η ‖Pδ0‖_W; if `mmax` is reached the last iterate
  (the best) is used and counted as a non-convergence.
* **ε**: see §4 (a deviation from "max|εv| = 1e-7").
* Fallback: if anything in the Krylov solve throws or is non-finite, the band correction stays and the event is counted (`fallback`).

## 3. T1: matrix-vector check

State: the stuck state (checkpoint `runs/t_nw_tools/nw0`, z100mm step ≈ 43,500), step 41 of the process, first trapezoidal-stage attempt at step Δt 1.5e4 s (TR stage 8.79e3), Newton
iteration 2; a random v (fixed LCG seed) over the 876 active unknowns; J v by central differences with the default ε; J_b w the band matrix applied to the change `update!`
actually made (including the ghost entries, which follow the boundary conditions). Per row, |Jv − J_b w| / Σ_j |J_b,ij w_j|. 888 rows are compared; the 4 Dirichlet rows
(R = δ) are left out. Files: `runs/t_nk25/seg_t_nk25.out` (closure held: `NWFREEZE=f,g,h,hx,fx,gx,tie`), `runs/t_nk26/seg_t_nk26.out` (not held).
The AD residual agrees with the Float64 residual to 0.000e+00 in both.

| | median | max | rows > 1e-4 | rows > 1e-2 |
|---|---|---|---|---|
| closure frozen (t_nk25) | 1.6e-6 | 2.92 (F_r, zone 16) | 63 | 18 |
| not frozen (t_nk26) | 1.7e-6 | 2.85 (F_r, zone 16) | 92 | 20 |

**The requirement (≲ 1e-4 per row with the closure held) is not met: 63 rows exceed it.** They are:

* **F_r rows of zones 15–24**: 2.92, 2.55, 2.27, 2.10, 0.97, 0.76 × rowscale for zones 16, 17, 19, 18, 15, 20; 0.20 (23), 0.18 (21), 0.035 (22), 0.027 (24). A row probe of the worst
  row (F_r, zone 16; one unit vector per active column, 876 products, `NKCHECK row probe` in `seg_t_nk25.out`) shows that the whole mismatch is the **x_H⁺ column of the zones 15 and 16**:
  finite-difference −2.456e31 and −2.411e31 against the band −2.280e30 and −2.715e30 (9–11 ×). Every other column of the row agrees: T −1.2644e30 / −1.2648e30,
  ρ +1.0399e30 / +1.0400e30, E_r +5.9277e29 / +5.9277e29, F_r +7.2344e29 / +7.2344e29 (FD / band, zone 16). The assembled ∂F_r/∂x_H⁺ (through the opacity of the coupled zones: NLOIT/NLOPAC
  slope refresh) is a factor ~10 too small in zones 15–24.
* the rows T (zones 16–23, up to 1.4e-2 at zone 17, 9.7e-3 at 18), u (zones 16–22, up to 3.0e-2 at zone 18), x_H⁺ and x_H₂ (zones 21–28, up to 4.1e-2 at x_H⁺ zone 27, 1.7e-2 at x_H₂ zone 27), E_r (zone 22: 1.5e-2): 1e-3 … 4e-2.
  Same zones; by the F_r result they are the x_H⁺ (and x_H₂) opacity/rate dependence the band Jacobian underestimates.
* without the freeze additionally F_r rows of zones 27–34 (`t_nk26`): 0.33 (29), 0.30 (28), 0.047 (27), 0.033 (33), 0.019 (32), 0.0062 (31) × rowscale: the closure (f, h, tie) that `update!` refreshes
  and the band Jacobian lacks, in the zones of the photosphere (the nwdiag NWFD columns).

These rows are what the band Jacobian lacks after the closure is held (measured). (Inference, not tested by a switch: the dt-independent linear rate of about 0.61, largest correction E_r in zone 22, is located in the same region (zones 16–28) as these rows.)

Dirichlet rows (u at ngrs and ngre+1, F at ngrs): R is the correction itself (`R[iu, kb(ngrs)] = X(iu, ngrs)`), J v = 0 there, J_b w = w.


## 4. Choice of ε and of the difference

One-sided probes (iterations at step Δt 3e3, 1.5e4, 4e4, 1e5, 1e6 [1e7] s from the stuck state, NITER = 60; plain Newton: 10, 31, >60, ...). "Weighted" = ε measured as the largest
weighted entry of εv (final), "plain" = largest entry of εv.

| variant | run | iterations |
|---|---|---|
| one-sided, ε = 1e-7, plain, unweighted 2-norm in GMRES | t_nk12 | no convergence in 60 at 3e3 ... 1e6 (dmax 3.5e-3 ... 9.5e-2), 19 at 1e7 |
| one-sided, ε = 1e-7, plain | t_nk14 | 36, 10, 7, 40, 38, no convergence |
| one-sided, 1e-6, plain | t_nk19 | 4, 6, 46, 13, no, no |
| one-sided, 3e-6, plain | t_nk20 | 23, 14, no, 13, no |
| one-sided, 1e-5, plain | t_nk21 | 6, 7, 15, 8, 34 |
| one-sided, 1e-4, plain | t_nk24 | 3, 7, 4, 10, 42 |
| central, 1e-5, plain | t_nk22 | 7, 7, 9, 10, 58 |
| **central, 1e-5, weighted (final)** | t_nk27 | **3, 4, 4, 4, 4, 5** (1e7 s) |

Reason for the deviation from "one-sided, max|εv| ≈ 1e-7": the residual carries round-off noise that comes from the formal solution of `eddfac!` (1 − (1 + τ)e^{−τ} in optically thin zones),
so f and h are noisy at the 1e-6 level; divided by a perturbation of 1e-7 that is a Jacobian error of order 10 in the F_r rows of the outer zones, and GMRES stagnates (relative residual 0.4–0.7
after 30 iterations, `NKR` lines of `t_nk14`, `t_nk17`). The sensitivity of J v (T1 state, random v, `NKCHECK eps x..` in `seg_t_nk25/26.out`; worst row relative to the row scale, 2-norm in brackets),
central difference, final ε = 1e-5:

| ε × | 0.01 | 0.1 | 1 | 10 | 100 |
|---|---|---|---|---|---|
| closure frozen: change of J v vs ε = 1e-5 | 1.8e-6 (1.1e-9) | 1.8e-6 (9.0e-10) | – | 5.1e-3 (9.0e-8) | 1.8e-2 (9.1e-6) |
| closure not frozen | 3.5e-2 (8.7e-10) | 2.0e-3 (8.5e-10) | – | 1.9e-2 (8.5e-8) | 2.6e-2 (8.6e-6) |

With the closure held J v is stable to 2e-6 for smaller ε (a smooth function); with it the F_r rows of the photosphere zones are noisy at the 3 % level at ε = 1e-7 (and 0.2 % at 1e-6), and
the truncation error is 0.5 % at 1e-4. ε = 1e-5 sits between. (One-sided at the same state, `t_nk10/11`: worst-row change 25 × rowscale at 1e-9, 0.06 at 1e-8, 0.27 at 1e-6 and 1e-5 for the unfrozen
closure; 2.6e-5 … 2.6e-3 for the frozen one.) The 2-norm changes are ≤ 9e-6 everywhere: the noise is confined to few rows, which the weighted norm and the GMRES minimisation tolerate.
Drift: the residual after the restore equals the residual before **exactly** (0.0, every Newton iteration of every run: `drift` in the NKR lines, summary "state-drift events 0").


## 5. T2: probes

Stuck state = `runs/t_nw_tools/nw0`, step 41 of the process, `NITER=60`, NWPROBE (a trapezoidal-stage attempt from the live state, retry-like; Krylov solve on at the probe only).
Plain Newton (`NKRYLOV` off) on this tree: `t_nk01`/`t_nk04` (ERLAB=1 / not): 10 / 14 / 31 iterations at step Δt 3e3 / 6e3 / 1.5e4 s (the nwdiag reference was 9 / 14 / 31: the 3e3 probe is one
iteration longer on this tree, the other two reproduce). `runs/t_nk27` (`NKRYLOV=1`, default parameters), raw `runs/t_nk27/seg_t_nk27.out`, trace `runs/t_nk27/newtrace.txt`:

| step Δt [s] (TR stage) | Newton iterations | final dmax | Krylov its per Newton iteration (incl. the closing iteration) | wall per attempt [s] | dmax sequence |
|---|---|---|---|---|---|
| 3e3 (1.76e3) | 3 | 2.5e-9 | 5, 5, 4, 5 | 6.9 (first call, includes compilation) | 5.4e-2, 1.2e-3, 6.4e-7 |
| 1.5e4 (8.79e3) | 4 | 1.9e-10 | 5, 5, 6, 2, 4 | 2.3 | 1.1e-1, 6.0e-3, 1.7e-5, 2.1e-9 |
| 4e4 (2.34e4) | 4 | 1.4e-9 | 6, 6, 5, 5, 2 | 2.3 | 1.6e-1, 9.9e-3, 5.0e-5, 1.6e-8 |
| 1e5 (5.86e4) | 4 | 1.4e-9 | 6, 6, 5, 5, 2 | 2.3 | 1.8e-1, 6.4e-3, 2.3e-5, 7.4e-9 |
| 1e6 (5.86e5) | 4 | 1.5e-8 | 6, 6, 5, 5, 2 | 2.4 | 1.9e-1, 1.1e-2, 1.4e-4, 2.9e-7 |
| 1e7 (5.86e6) | 5 | 3.8e-8 | 5, 6, 6, 5, 2, 3 | 2.7 | 3.9e-1, 3.2e-2, 4.2e-4, 1.9e-5, 2.1e-7 |

Gate (converged at every Δt up to 1e6 in ≤ 8 Newton iterations): **passed** (3–4), also at 1e7 (5). Convergence is superlinear (the successive ratios fall: 0.02 → 5e-4 → …).
Summary line: 30 Newton linear solves, 140 matrix-vector products (4.67 per solve), 0 GMRES non-convergences, 0 fallbacks, 0 drift events. Plain Newton took 31 iterations in 2.5 s at 1.5e4 (cheap
iterations): the gain of the Krylov solve is the iteration count (and with it the step size the Newton can take), not the cost of an attempt.

Later state (z100mm step 45,500, `m2c_z19z100mmS12`, copied; step 41 of the process), `runs/t_nk28` (NKRYLOV) and `runs/t_nk29` (plain, same tree):

| step Δt [s] | plain (t_nk29) | NKRYLOV (t_nk28) | Krylov its per Newton iteration | dmax sequence (NKRYLOV) |
|---|---|---|---|---|
| 1.5e4 | 39 | 4 | 6, 5, 5, 5, 2 | 1.6e-1, 3.8e-2, 1.1e-3, 1.5e-6 |
| 1e5 | 40 | 5 | 6, 5, 5, 5, 4, 5 | 1.9e-1, 5.7e-2, 1.2e-2, 1.6e-4, 2.0e-7 |
| 1e6 | 28 | 5 | 5, 5, 5, 4, 5, 5 | 3.9e-1, 1.9e-1, 4.6e-2, 1.9e-3, 5.6e-6 |

(The task's reference for this state, 19 / 16 / 24 iterations, is not reproduced by the plain run on this tree at the same probe point: 39 / 40 / 28; the table compares like with like.) Summary: 17 solves, 82 products (4.82 per solve), no non-convergence.


## 6. T3: 400 steps

400 steps from the stuck checkpoint `nw0` (step 43,500), ERLAB=1, NITER = 25 (unchanged), 4 threads. `runs/t_nk30` (`NKRYLOV=1`), control `runs/t_nk31` (this tree, NKRYLOV off; it reproduces
nwdiag's `t_nw02`: 103 Newton failures, 0 error rejections), `runs/t_nk41` (NKRYLOV=1 with the step capped, `DTERRMAX=1.2e4`, for the equal-time comparison). Statistics from `seg_*.out`
("DTERR statistics"), `perf.csv`, `export/track.csv` (made with `julia/runrec.jl`).

| | control (t_nk31) | NKRYLOV (t_nk30) | NKRYLOV, dt ≤ 1.2e4 (t_nk41) |
|---|---|---|---|
| Newton failures (per 100 steps) | 103 (25.8) | 11 (2.8) | 1 (0.25) |
| error-estimate rejections | 0 | 12 | 0 |
| mean accepted Δt, steps 43,550–43,900 [s] | 5.6e3 | 2.6e6 | ≤ 1.2e4 (all steps at the cap) |
| simulated time over the 400 steps (end − control's end) | 0.0625 yr in steps 43,550–43,900 | +35.87 yr more than the control (timen 7.507832e12 vs 7.506700e12 s) | 0.133 yr (steps 43,550–43,900) |
| time reached (track.csv, step 43,900) [yr] | 237854.884 | 237890.744 | 237854.966 |
| wall per step [s] (perf.csv, step 43,900) | 3.10 (1239 s) | 7.29 (2917 s) | 4.37 (1747 s) |
| wall per simulated yr [s] (400 steps; the control's 400-step time estimated as 400/350 × 0.0625 yr) | ~1.7e4 | 81 | 1.15e4 (400 × 1.2e4 s = 0.152 yr) |
| Krylov counters | – | 5952 solves (14.9 per step), 23487 matrix-vector products (3.95 per solve), 0 non-convergences, 0 fallbacks, 0 drift events, 2078 s in the Krylov solve (71 % of the wall) | see seg_t_nk41.out |

**The step runs away, as expected**: with the Newton wall gone Δt grows to 3e7 s (steps 43,550 …) and is then limited by the error estimate (rtol 0.1): the limiting variable of the accepted
steps is x_H₂ (250 steps) and x_H⁺ (137), in zones 27–31 (the photosphere; the `DTERR jstep` lines, `max 0.5 … 38` of the rejected: x_H₂ zone 31), E_r in 9 and T in 2. The 11 Newton failures are at Δt ≥ 9e5 s
(one restart-artefact failure at 5.3e3 in step 1 as in every restart); the 12 error rejections (all x_H₂ zone 31, err 1.0–1.4) shrink Δt by ×0.8. One `bound` event (mesh motion bound).
No fallback, no non-converged GMRES solve in 5952 solves.

Does the solution stay sane? Measured (track.csv, `runs/t_nk30/export/track.csv`): M* 1.373 → 2.609 M☉ in 35.9 yr (mean accretion rate 0.0344 M☉/yr; the control 0.0376 M☉/yr over its 0.06 yr; the task's 0.036),
R* 107.5 → 138.0, r_ph 115.3 → 148.2, T_c 128377 → 133106 K, all monotone. L_ph: 4320 (step 43,550), 4319 (43,600), **7802 (43,650)**, 6330, 6463, 6432, 6432, 6437: one excursion of +80 % at step 43,650
that relaxes within 11 yr; not diagnosed (the control, which spans 0.06 yr, cannot show it).

Solution at equal simulated time (`t_nk41`, Δt capped near the control's Δt, against `t_nk31`; track.csv interpolated in t to the control's last output, t = 237854.8844 yr; and 237854.84, .86):

| quantity | t = .8400 ctl / NK / rel | t = .8600 | t = .8844 |
|---|---|---|---|
| M* [M☉] | 1.37296 / 1.37296 / 1.0e-6 | 1.37369 / 1.37369 / 1.7e-6 | 1.37457 / 1.37458 / 2.5e-6 |
| R* | 107.526 / 107.539 / 1.2e-4 | 107.549 / 107.571 / 2.0e-4 | 107.589 / 107.617 / 2.6e-4 |
| r_ph | 115.347 / 115.360 / 1.1e-4 | 115.374 / 115.396 / 1.9e-4 | 115.418 / 115.446 / 2.4e-4 |
| L_ph [L☉] | 5654.6 / 5679.1 / 4.3e-3 | 5624.7 / 5667.5 / 7.6e-3 | 5599.3 / 5666.3 / 1.2e-2 |
| T_c [K] | 128384 / 128381 / −1.8e-5 | 128387 / 128383 / −3.4e-5 | 128391 / 128385 / −4.5e-5 |

The solutions agree to a few 1e-4 in the radii and 1e-6 in M*; L_ph differs by 0.4–1.2 % (L_ph drifts 0.3 % per 0.01 yr in the control and wobbles ±0.4 % in the capped run: 5690, 5662, 5674, 5664, 5654 …).


## 7. T4: from the initial conditions

From the initial condition, the standard settings of the reference run `runs/t_ic_l` (its `phaseB.env`, `run.env`, 128 IC shells, 100 zones, 4 threads, `ERLAB=1 EDIAG=2`) plus `NKRYLOV=1`, `runs/t_nk40`
(launcher `runs/t_nk_tools/nkic.sh`), stopped by hand at step 1100 (M* = 0.2038 M☉). Reference: `t_ic_l` (M* = 0.2003 at step 1300). Both from `perf.csv` (every 100 steps; ndiv = cumulative Newton
divergences, col. 13; fail_tr, col. 16):

| | step | t [yr] | wall [s] | ndiv | fail_tr | T_c [K] | M* [M☉] | R* |
|---|---|---|---|---|---|---|---|---|
| t_ic_l (plain) | 1100 | 2.41418443e5 | 2557.8 | 53 | 47 | 1.0882e5 | 0.1661 | 28.56 |
| t_ic_l (plain) | 1300 | 2.41419022e5 | 2926.8 | 60 | 53 | 1.1111e5 | 0.2003 | 33.20 |
| t_nk40 (NKRYLOV) | 1000 | 2.40562020e5 | 5771.5 | 32 | 21 | 1.0905e5 | 0.1914 | 29.93 |
| t_nk40 (NKRYLOV) | 1100 | 2.40562236e5 | 6330.2 | 43 | 32 | 1.1017e5 | 0.2038 | 33.33 |

Result: **passes core formation** (T_c 1.10e5 K, M* 0.2038 at step 1100, R* 33.3). 1100 steps against 1300 for the reference (the Newton wall no longer shrinks Δt),
**but 6330 s against ~2930 s wall** (2.2 ×): every Krylov iteration is one central-difference product (two `update!` + context + residual evaluations); the DTERR lines give 2.6 Krylov iterations per Newton iteration and ~35 matrix-vector products
per step on average (mean of the 55 trace lines; last interval 2.69 and 23.0). Newton failures: 43 against 60 (53 at step 1100): not the 25-fold drop of the stuck state — the remaining failures are 10 in the
first 61 steps (Δt = 2e10–5.4e11 s, the start of the collapse; the one in step 15 with a GMRES fallback), 11 at the core formation (steps 1023–1089, Δt 3.6e3–6.5e5 s) and 22 in between (Δt 1.5e4–2.3e7 s), not the closure wall of the stuck state. 0 GMRES
non-convergences in the trace lines; **68 fallbacks** to the band correction (`non-finite J v`: a perturbed evaluation gave a non-finite residual; one `SingularException`), concentrated in the intervals
of the first steps. The summary line of this run is missing (stopped by hand).

The run reaches the same stage at a different time: t(M* ≈ 0.20) = 2.40562e5 yr against 2.41419e5 yr (−0.36 %). The two runs take different Δt through the first 60 steps (they differ in which attempts fail);
whether the time difference is the time-step error of the collapse phase or something else is **not examined** (it would need the plain run at the same steps).


## 8. T5: bit-identity with the switch off

60 steps from the stuck checkpoint (`nw0`), ERLAB=1, `DTERRTRACE=1`, 4 threads, NKRYLOV not set: the unmodified base tree (`/home/tabel/Projects/titan.jl`, popiii-s6 a4837c7, read only) `runs/t_nk32` and this tree `runs/t_nk33`.
Compared: the decompressed final checkpoint state (`zcat ck/m2c_nk3?S1.tstate.gz | md5sum`: 97e00ce1ce80cb88de8e8a7328f7ee00 in both), the `.aux` files (md5 7350d11ed7eb7123e9f9c41d53364611 in both),
and the 60 `DTERR jstep` lines (md5 a737fa41ef456b4a09a6db67f55ffe81 in both; `diff` empty). Both runs: 14 Newton failures in the 60 steps. **Bit-identical.**


## 9. Open issues

RESULTS_OPEN
