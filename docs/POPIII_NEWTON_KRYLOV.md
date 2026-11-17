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

RESULTS_T1

## 4. Choice of ε and of the difference

RESULTS_EPS

## 5. T2: probes

RESULTS_T2

## 6. T3: 400 steps

RESULTS_T3

## 7. T4: from the initial conditions

RESULTS_T4

## 8. T5: bit-identity with the switch off

RESULTS_T5

## 9. Open issues

RESULTS_OPEN
