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
and the 60 `DTERR jstep` lines (md5 a737fa41ef456b4a09a6db67f55ffe81 in both; `diff` empty). Both runs: 14 Newton failures in the 60 steps. **Bit-identical.** Repeated on the final commit (after the cleanup commit, `runs/t_nk34`): the same three md5 sums. The final code with the switch on, probes at 1.5e4 and 1e6 s (`runs/t_nk35`): 4 and 4 Newton iterations.


## 9. Open issues

1. **(Found and fixed on branch `jacfix`, §10.) The second inexactness is located, not removed.** With the closure held, the band Jacobian still lacks ∂F_r/∂x_H⁺ in zones 15–24 (a factor ~10 in the entries of zones k−1, k) and is 0.1–4 % off in the
   T, u, E_r, x_H⁺, x_H₂ rows of zones 16–28 (T1). Candidates are the opacity of the coupled zones (NLOIT/NLOPAC slope refresh, `CPLJLAG` reuse); not tested which. Whether this is the dt-independent 0.61 rate
   is an inference from the location (the largest correction E_r, zone 22). The Krylov solve absorbs it (3–5 Newton iterations at every Δt) but it would be cheaper to fix the assembled entry.
2. **Finite-difference noise.** `eddfac!` evaluates 1 − (1 + τ)e^{−τ}-type expressions in optically thin zones; f and h carry ~1e-6 round-off, so J v in the F_r rows of the outer zones needs ε ≳ 1e-5 (central)
   and the one-sided ε ≈ 1e-7 of the first design does not work (§4). A cancellation-free evaluation (expm1) would change `eddfac!` at the round-off level (not bit-identical), so it was not touched.
3. **Cost.** In the 400-step run the Krylov solve is 71 % of the wall (2078 of 2917 s; 15 linear solves and 59 products per step because the steps are 460 × longer and each retry repeats them). Per simulated year it
   is ~200 × cheaper (T3), per step 2.4 × dearer; from the initial condition (T4) the wall to M* = 0.2 is 2.2 × longer, not shorter (step count 1100 against 1300). Each linear solve takes a full snapshot of
   the Julia state (all arrays of `TF.JSTATE` except the Jacobian and rhs) per solve and restores it per product and refactors the band matrix with pivoting; both could be cut (snapshot only what `update!` writes; reuse the matslv factorisation) — not done.
4. **Fallbacks.** In T4, 68 linear solves fell back to the band correction (`non-finite J v`: a perturbed update gave a non-finite residual in the start of the collapse; one SingularException). In T3 none.
5. **T4 timing of core formation** differs by −0.36 % from the reference (§7), not diagnosed. L_ph shows one +80 % excursion at step 43,650 in the free-running T3 (§6), not diagnosed.
6. **Not run**: the `TitanDriver.step!` path (θ-scheme, non-DTERR) uses the same `newton_linear!` but was not exercised with the switch on (it is bit-identical with the switch off by construction: same calls);
   `NITER` smaller than 25 with NKRYLOV; `eta`/`mmax` variations (defaults 1e-3 and 30: no non-convergence in 5952 + 30 + 17 solves of T2/T3, 0 in T4's trace lines; the Krylov iterations per solve are 2.6–5).
7. For the production run (not changed by this work): with NKRYLOV the Newton failures per 100 steps fall from 25.8 to 2.8 and the mean accepted Δt rises from 5.6e3 to 2.6e6 s from the stuck state, so a restart of
   z100mm with `NKRYLOV=1` is the experiment to do; the step then runs to the error-estimate limit (x_H₂/x_H⁺ in the photosphere zones 27–31) and the accreted mass per wall-hour rises ~200 ×.


## 10. Jacobian defect found by the matrix-vector check

Branch `jacfix` (from `nkrylov` 0607bad). Runs `t_jf01`–`t_jf26` and `t_jb01`–`t_jb03` (`jb` = the unfixed tree, launcher `runs/t_nk_tools/jbrun.sh`; `jf` = this tree, `jfrun.sh`), plus `t_nk43`.
`NKCHECK` gained `NKROW=eq,zone[;eq,zone…]` (diagnostic: the rows of the column-by-column probe; zone counted from the inner boundary as in the listing).

### Cause

Four quantities that enter the residual had no derivative, or a wrong one (all in the zones 15–28, the x_H⁺ ≈ 0.99–0.9990 photosphere/ionisation front; the dependence is local, within the ±2-zone band):

1. **x_H⁺ column of the F_r row (the finding of T1).** The F_r row depends on x_H⁺ only through the opacity (avchi·F_r). The opacity of a coupled zone is `PopIIIMaterial.nlo_it`; its derivatives
   in (ρ, T, x_H⁺, x_H₂) of the non-LTE ratios R were differences (`PopIIIMaterial.jl:622` and following, 0607bad): **one-sided, relative step 1e-2 in x_H⁺**. In these zones the neutral fraction
   n_H/n = 1 − x_H⁺ is ~1e-3 (d[:H] = 0.0011 at x_H⁺ = 0.9989): a step of 1 % in x_H⁺ crosses x_H⁺ = 1, `max(d[:H], 0) = 0`, and the difference is a different, flat branch. The measured slope
   d ln χ/d ln x_H⁺ was −6.2 (zone 15), −7.6 (16), −45 (17) against −67, −67, −174 for the central difference of the Float64 residual at h = 1e-6…1e-4 (and −33 at h = 1e-2, the same function seen
   over the step: the dependence has a scale of 1e-3 in ln x_H⁺). Dual numbers have the same answer as the 1e-6 difference. The electron floor of χ (σ_T n_e/ρ) was also differentiated by hand
   (`ned = nebase·x/x_v`, no T or ρ dependence; `PopIIIMaterial.jl:666`).
2. **The gas–radiation exchange cap** `excap` (`TitanResidFast.jl:759`): `min(1, EXCAP/(c val(κ_E) val(ρ) Δt))`, `val` = no derivative ("a constant of the iteration"). Where it is below 1 (opaque zones; EXCAP = 1e8 in the
   production settings) the energy rows' exchange term f·(cκ_E E − 4πκ_P B)ρ dV then lacked ∂f/∂κ_E: the product f·κ_E·E·ρ = EXCAP·E/(cΔt) does not depend on κ_E at all, the held factor made it depend on κ_E with the full c·E·ρ. This was the E_r row: the x_H⁺ column of E_r, zone 18 was +2.8e41 in the band against −4.8e45 (FD).
3. **Photo-rate temperature dependence in the CK expansion** (`ChemNet.jl:745`): `ck_eval` took the photo-ionisation/recombination rates at `ForwardDiff.value(T)`, so every Jacobian through the CK expansion (the coupled zones' local
   Jacobians `cpl_vj`, and the opacity slopes) lacked ∂(photo-rates)/∂T. For the composition in zone 17: d ln n_e/d ln T = 0.031 (Dual) against 0.117 (difference); with the rates frozen at the base T the difference gives 0.0307, so this is the whole defect.
   This is what the 1e-3…4e-2 mismatch of the T, u, E_r, x_H⁺, x_H₂ rows of zones 16–28 was (T diagonal of the T row in zone 18: 30 % off).
4. (not a defect of the old Jacobian, a prerequisite of the fix) `G2Opacity.jl:199` formed the ratio of two occupation probabilities as `occupation(ns)/max(occupation(n), 1e-300)`; in Duals the quotient rule squares the denominator, which underflows for
   probabilities of 1e-100 (the dense zones): NaN partials. The ratio is now taken from the logarithms (`Z0Physics.ln_occupation`); the Float64 value is the same to round-off.

### Fix

* `PopIIIMaterial.nlo_slopes`: one ForwardDiff evaluation (JTag Duals, four partials in ln ρ, ln T, ln x_H⁺, ln x_H₂) through the same `nlo_comp` + `nlo_ratio` that give the value, at the iterate's star-light
  factors (q, t); gives all slopes of R_E, R_P, R_R **and of n_e** exactly (the electron floor is now differentiated through the same expression). Replaces six extra evaluations of R by differences. (For `G2OPAC` ≠ 2 the
  means are Float64-only and the old differences remain; the production setting is `G2OPAC=2`.) The slopes in ln q and ln t of SLIMP (differences, `PopIIIMaterial.jl` SLIMP block) are untouched.
* `TitanResidFast.excap`: no `val`: the cap factor is differentiated with κ_E and ρ.
* `ChemNet.ck_eval`: `phc(T)` with the Dual T; `PopIIIMaterial.photo_cache` caches the rates on the last Float64 T and evaluates a Dual T directly (the photo-rate tables are Catmull–Rom in ln T, differentiable).
* `Z0Physics.ln_occupation`, `G2Opacity._means_fast` (above), and `n2_lte` without `Float64(T)`.

Diff stat 0607bad..jacfix: 6 files, +86 −34 (ChemNet 1, G2Opacity 7, NewtonKrylov 34 (`NKROW`), PopIIIMaterial 67, TitanResidFast 4, z0/atomic 6).

### Before / after

Stuck state `nw0` (t_nk30 checkpoint), `ERLAB=1 NITER=60`. Before = tree nkrylov, after = this tree.

a. `NKCHECK`, closure frozen (`NWFREEZE=f,g,h,hx,fx,gx,tie`, `NWPROBE=1.5e4 NKRYLOV=1`; iteration 2 as in T1):

| | before (t_nk25) | after (t_jf18) | after, `CPLJLAG=0` (t_jf23, iteration 1 state, see below) |
|---|---|---|---|
| max row mismatch | 2.920 | 4.86e-3 | 4.9e-6 |
| median | 1.61e-6 | 1.52e-6 | 1.38e-6 |
| rows > 1e-4 / > 1e-2 | 63 / 18 | 18 / 0 | 0 / 0 |
| F_r zone 16, x_H⁺ zone 15: FD / band | −2.4562e31 / −2.2803e30 | −4.8352e31 / −4.8352e31 (|Δ| 1.8e25, 3.7e-7) | |
| F_r zone 16, x_H⁺ zone 16: FD / band | −2.4112e31 / −2.7153e30 | −3.4525e31 / −3.4525e31 (|Δ| 7.5e25, 2.2e-6) | |

(The two columns are at slightly different states: the Newton step of iteration 1 differs with the Jacobian, so the FD values in the rows differ between the columns. `AD residual vs Float64 residual: 0` in all.) The 18 rows > 1e-4 of t_jf18 (largest x_H₂ zone 26, 4.9e-3; then x_H₂ zones 22, 24, x_H⁺ zones 29, 26, 22) are the **deliberate reuse of the local chemistry Jacobian** (`CPLJLAG`, default 3e-3,0.1:
the Jacobian of a coupled zone is recomputed only when ρ or T moved by more than 3e-3 or x by more than 0.1; "the Newton becomes inexact, the converged state is unchanged", `PopIIIMaterial.jl` CPLJLAG): with `CPLJLAG=0` they are gone (t_jf23: 0 rows > 1e-4, max 4.9e-6 = the FD noise floor of eps 1e-5, median 1.4e-6 as before). Not changed (a design decision of the production settings).
Contribution of each part (iteration-1 state, frozen closure, max / rows > 1e-4 / rows > 1e-2): t_jb01 (before) 2.92 / 57 / 22; + slopes and occupation (t_jf09) 0.277 / 60 / 12; + EXCAP (t_jf10) 3.5e-2 / 54 / 2; + photo-rates (t_jf16, `CPLJLAG` default) **4.9e-6 / 0 / 0**.
(With the default `CPLJLAG` the iteration-1 matrix is presumably the freshly computed one, so the lag shows only from iteration 2; not checked.)

b. `NKCHECK`, closure live (no `NWFREEZE`), iteration 2: before t_nk26 max 2.853, 92 rows > 1e-4, 20 > 1e-2; after t_jf19 max 0.326, 59, 5. The 5 rows > 1e-2 are F_r zones 27–33 (0.326, 0.305, 0.048, 0.032, 0.018: the closure rows 27–34 of the Eddington factors,
expected to remain; the Krylov solver's job), then F_r zones 60, 64 (8.8e-3, 6.9e-3: the finite-difference noise of f, h in the thin outer zones, §Open issues 2).

c. Plain Newton at the stuck state (`NWPROBE=3e3,1.5e4,4e4,1e5,1e6`), iterations (61 = no convergence in NITER = 60):

| Δt [s] | 3e3 | 1.5e4 | 4e4 | 1e5 | 1e6 |
|---|---|---|---|---|---|
| plain, before (t_nk43) | 10 | 31 | 61 (dmax 3.3e-5) | 61 (4.4e-4) | 61 (1.9e-2) |
| plain, after (t_jf20) | 11 | 35 | 61 (9.1e-5) | 61 (2.1e-3) | 61 (5.3e-2) |
| NKRYLOV, before (t_nk27) | 3 | 4 | 4 | 4 | 4 |
| NKRYLOV, after (t_jf21) | 3 | 3 | 4 | 3 | 4 |
| plain, closure frozen, before (t_jb02) | 13 | 13 | 13 | 12 | 11 (final dmax ~5e-6) |
| plain, closure frozen, after (t_jf22) | 9 | 8 | 9 | 9 | 9 (final dmax 3e-8…6e-8) |

**The missing x_H⁺ dependence is not the cause of the dt-independent contraction rate of plain Newton: with it fixed, plain Newton converges no better (11, 35, fail, fail, fail).** With the closure frozen the exact band Jacobian
gives the 8–9 iterations at every Δt (before 11–13), i.e. what stops plain Newton at Δt ≥ 4e4 s is the closure (`eddfac!` f, g, h, tie weights that `update!` rewrites: the F_r rows 27–34 of b), consistent with T1's attribution; the frozen problem is still
not quadratic (8–9 iterations from the probe's first iterate; not investigated further).

d. 400 steps from the stuck state (`nw0`), `ERLAB=1`:

| | control t_nk31 (before) | t_jf25 (after) | t_nk30 NK (before) | t_jf26 NK (after) |
|---|---|---|---|---|
| Newton failures (`DTERR statistics`) | 103 | 103 | 11 | **51** |
| error rejections | 0 | 0 | 12 | 5 |
| mean step over steps 43550–43900 [s] | 5637 | 5592 | 2.58e6 | 3.52e6 |
| simulated time of steps 43550–43900 [yr] | 0.0625 | 0.0620 | 28.57 | 39.08 |
| wall [s] | 1239.4 | 1071.3 | 2916.8 | 3034.6 |
| Krylov: solves / products | | | 5952 / 23487 | 6501 / 25858 |

The plain run is unchanged in failures and step size (103, 5.6e3 s); its wall was 14 % shorter (not examined whether that is the cheaper opacity slopes — one 4-partial Dual evaluation instead of six differences — or run-to-run variation on the shared machine).
**The NK run has more Newton failures after the fix (51 against 11).** The two NK trajectories separate (after 28.6 and 39.1 yr; M* 2.608 and 3.114 M☉; M* gain over steps 43550–43900 per year: 0.0341 and 0.0337 M☉/yr); the failures in t_jf26 cluster at the large steps (dt 5e6–8e7 s at steps 43552, 43650, 43705, 43749, 43768) and in one episode at step 43769 (dt 2.6e4 → 0.1 s in 12 halvings, then recovered).
Why the corrected Jacobian does not reduce them, and whether the larger steps the controller takes are the cause, is **not diagnosed**; one 400-step trajectory each is no statistic.

e. NKRYLOV unset: the code path does not touch `NewtonKrylov.jl` (only `check!`, run with `NKCHECK`, changed there). **The Jacobian fix does change the default path** (intended): the assembled band matrix and thereby the Newton iterates differ; the residual function is the same
(its Float64 value changes only at round-off: the occupation ratio is taken from logarithms). Converged state after 50 steps from `nw0` (t_jf24 against t_jb03 = tree nkrylov, `export/track.csv`, step 43550; both runs have 12 Newton failures in the 50 steps): M* 1.37230817 / 1.37230810 (rel 4.9e-8), R* 107.510463 / 107.510458 (4.4e-8), r_ph 4.6e-8, T_c 3.3e-8,
L_ph 5697.0226 / 5697.0154 (1.3e-6). After 350 further steps (t_jf25 against t_nk31, step 43900): M* 1.0e-5, R* 3.2e-4, L_ph 1.7e-2 (between steps 43600 and 43850: R* ≤ 9e-6, L_ph 2e-6…8.5e-4; the last output, at step 43900, is the one where the two runs differ most; not examined).

### What remains

* Closure rows: F_r zones 27–34 (and the closure's coupling into E_r, T) with the closure live: 5 rows > 1e-2 (§b). Unchanged by design (`update!` rewrites f, g, h, tie after the Newton correction); the Krylov solver handles it.
* `CPLJLAG` (default 3e-3,0.1): 18 rows > 1e-4 (max 4.9e-3) at iteration ≥ 2; deliberate.
* FD noise of f, h in the outer zones (F_r zones 60, 64: 7–9e-3 with the closure live; below 1e-4 with the closure frozen).
* The SLIMP slopes in ln q, ln t are still differences (relative steps 1e-2, 1e-3, `PopIIIMaterial.jl` SLIMP block): the E_r and F_r columns of the opacity agree with FD to 2e-6 in the frozen test, so no defect was seen, but they are not exact Duals. `cpl_vjq`'s t-derivative is a central difference (ln t, 1e-3).
* Plain Newton at Δt ≥ 4e4 s still does not converge (§c); the cause is not the Jacobian entries tested here.
* NK with the fixed Jacobian: more failures in the one 400-step run (§d), not diagnosed.
* Not run: T4 (from the initial condition) and the production restart.


## 11. NKLS: backtracking line search on the preconditioned residual (switch, default off)

What: `NKLS=1[,trigger[,lmin]]` (popiii_m2c.jl, next to NKRYLOV; defaults trigger 1e-2, lmin 1/64). After GMRES returned the step x, the level function
T(lambda) = || W P M^-1 R(x_k + lambda x) ||_2 (same weighted norm and band-LU preconditioner as GMRES) is evaluated with the TRUE residual at the trial state
(`perturbed(...; freeze = 0)`: update! applied, branch decisions re-evaluated, state restored in a `finally`). lambda = 1, 1/2, ... >= lmin; the first lambda with
T(lambda) <= (1 - lambda/4) T(0) and T finite is taken; a non-finite residual or an exception counts as a rejected trial. If none is accepted the finite trial of smallest T
is used if T < T(0), else lmin (counted as lsfail). It engages only if max|W x| > trigger. `blk` receives lambda*x on the active entries. maxdel (including its clamp
for dmax > dtol) and update are unchanged and act on lambda*x.
Why: the trace of t_jf29 attempt 26 shows the full Newton step raising ||M^-1 R||_W by a factor 30 with nothing to stop it.
FREEZE semantics (TitanResidFast.decide): 0 returns the branch condition and touches no record; 1 records, 2 replays; so freeze = 0 neither records nor consumes.
T(0) is computed from S.R0 and agrees with ||W delta0|| to a relative 1.5e-10 at most over all engaged iterations (printed per line as T0rel; summary line).
Trace: the "NKR att" line ends with ` ls lambda .. T0 .. T .. trials .. scnorm0 .. scnorm .. T0rel ..` (or ` ls off`); summary_text() adds the counters.
conv (Newton tolerance on dmax) is 1e-3 in these runs (run.env). dmax is an unweighted measure, so a damped step is not excluded from satisfying dmax <= conv by the
trigger alone: at lambda = 1/64 and max|Wx| = 1e-2 the raw correction is 1.6e-4 if the weights are 1. Not observed to matter, not guarded.

Measured (restart at step 62,300 of z100nk_ep, ERLAB=1 NKRYLOV=1 DTERR=0.03, 80 steps, 4 cores; output /zpool/nvme/data/tabel_scratch/titan_popiii/runs/t_ls0N):
* NKLS unset, 20 steps: identical DTERR lines to t_jf27 (the same code path before the change).
* NKLS=1 (t_ls02): 9 Newton failures (dt 5.0e3, 1.75e4, 3.07e4, 2.98e4, 2.4e4, 3.18e4, 3.28e4, 3.39e4, 6.18e4 s) against 11 without (t_jf27); 457 line searches engaged,
  189 backtracked, 81 failed, 1271 extra residual evaluations; 12.5 min wall incl. start-up.
* NKLS=1 NITER=60 (t_ls03): the same 9 failures; 769 engaged, 469 backtracked, 194 failed; 16.0 min wall.
* In the failing attempt 26 (jstep 62310, 1.027e4 s) |delta0|_W stays at 2.3e-3 to 2.6e-3 while |x|_W is 0.1 to 3.4; lambda goes to lmin or all trials are non-finite
  (T NaN) in 13 of 25 iterations, and T does not fall below 2.2e-3: no descent along the GMRES direction.
