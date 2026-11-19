# Handoff: Pop III accretion solver (titan.jl), 19 Nov 2026

For an agent continuing this research program on another machine. Read `docs/POPIII_GOAL.md` first.

## Goal (Tom Abel)
Fast, stable code for very long Pop III accretion simulations: protostar formation, final accretion, HII-region breakout.
Ionisation, H⁻ photodetachment and H₂/LW must be correct. Fully consistent and conserving; clean code without stale switches.

## Working rules Tom has set
- **100 zones** is the default for runs and tests. It is not a hard rule: convergence tests at 50 or 200 zones are fine.
- **No lagged inputs.** Everything is evaluated at the Newton iterate, with derivatives.
- **Jacobian consistency is a gate.** The FD J·v check (`NKCHECK=2`) must pass.
- **No clamps, damping or new switches.** Use general methods, not patches, and delete the old path when you replace it.
- **Physics over cost.**
- **Ask Tom before any method or config change to a production run.**
- **Cheap agents and short reports.** Implement through cheap agents (Sonnet), keep reports short, and never poll runs with sleep loops: run them in the background and wait.

## Code: branch `popiii-combined` (this commit)
It holds, in order of merge:
1. Old-level consistency (`peq!` refresh), which gave SOD 464 → 397 steps with smaller L1 error.
2. The new SOD baseline (397 steps), with L1 checks (`julia/SodReference.jl`).
3. Energy (iecons):
   - Compatible internal-energy rows, used when `ETOTAL=0`.
   - Total-energy monitor `ETOTMON`.
   - Heat booking `CQB`/`YQB` for zones that change coupled status.
   - Conversion of old checkpoints on restart.
   - Node-mass gravity force, which is the exact −∂W/∂r. The exact dW/dt is booked in the monitor and the total row. See `docs/POPIII_ENERGY.md` §23.
4. Lagrangian first interface `ICMLOG=m1`: the IC mesh is logarithmic in enclosed mass.
5. Newton globalisation (nkglob):
   - Updates of ρ, T, E_r on a C¹ positive map; species update in log space.
   - Steihaug–Toint trust region and Armijo on the weighted residual.
   - Newton difficulty feeds into the dt control.
   - Clean retry reset.
   - CPLCLAMP, NKLS and the NK band fallback are deleted. See `docs/POPIII_NEWTON_KRYLOV.md`.
6. Lagged-input removal from branch `lagfix` up to 7fe4791:
   - A3: starlight, columns, Lyα, all at the iterate.
   - A6: HD.
   - The **smooth deep-zone weight** `sl_deepw(τ)`, which replaces the hard flag `deep = τ ≥ 3`.
   - lagfix commits after 7fe4791 (the ALI/EP-factor model) are **not** included: they raised Newton failures from 15 to 49.
7. Stall diagnosis and fixes (stalldiag), with diagnostics `SD*`, all off by default:
   - Closure derivatives in the band Jacobian, `CLOSJ=1` (default in `popiii_stack.jl`).
   - Numerically stable thin-zone formal-solution terms.
   - Textbook Steihaug boundary step.

Regression: `julia/popiii_regress.jl` accepts with SOD 397, SUBC 587 and SUPC 713 steps. Further tests: `noh_test.jl` (production AV, TR-BDF2), `selfgrav_energy_test.jl` and `popiii_tests.jl`.

## The main open problem: stall at M* ≈ 3–3.5 M☉
Every 100- and 128-zone run stalled at M* ≈ 3–3.5 M☉, with dt falling from 3e4 s to about 70 s:
- z100mm and z100nk (100 zones, IC abn02red, GRIDEQ moving mesh);
- z128nk (128 zones).

All three runs are stopped.

What has been established, with runs t_sd01–t_sd89 on the old machine:
1. **Cause 1, found and fixed.** The hard flag `deep = τ ≥ 3` was set from the previous step's state, and zone 30 sat at τ ≈ 3. The flag switched off the three-level heat and the H bound-free opacity, which feeds back on τ. The result was a dt-independent 4-step limit cycle in E_r, T and x_H+ of zones 30–33. With `sl_deepw` the cycle is gone.
2. **Cause 2, partly fixed.** The Eddington/tie closures (f, h, g, tie weight) were refreshed every iterate but not differentiated in the band Jacobian, so P·J_b⁻¹·J was far from I and GMRES stalled. With `CLOSJ`, the stalled state runs 300 steps with 8 Newton failures instead of 28, and dt grows to 4e5 s.
3. **Still open.** From the pre-stall checkpoint (S6, production settings) the run collapses at step 85,638 (30 Newton failures in 190 steps). The largest corrections are in the **F_r rows of the outer envelope, zones 90–101**.
   - NKCHECK row F_r(101) has mismatch 0.144 in its E_r(100) column (FD −984 vs band −855).
   - The outer-boundary Eddington ratio g (`P.geddr`, from `eddfac!` in TitanCore.jl, used with `vef_esurf` in the F_r row at ngre+1 in TitanResidFast.jl) **sits at a clamp, 0.50000101**.
   - **Next task:** find the clamp and the missing E_r(100) term, and replace them with a smooth, differentiable closure (a Dual, see `closure_lift`). Then rerun S6 past t = 2.379230e5 yr.
   - A local agent was doing exactly this on branch `stalldiag` when this note was written. Its uncommitted work is not in this commit, so redo the task here.
4. Other hard flags held from the start-of-step state are listed in `PopIIIMaterial.jl`. They did not change in the stall runs, so they are not urgent:
   - the photosphere index `kph` and `SLSTAR`;
   - `ph`/`jph`;
   - the `wq ≥ 1` branches, which are continuous;
   - `kHon` (MGROUPS only).

## States to reproduce (`handoff/states/`)
| Checkpoint | Step | t (yr) | dt (s) | Notes |
|---|---|---|---|---|
| `m2c_z19z100nkS6` | 85,500 | 237,916.4 | ~3e4 | before the stall |
| `m2c_z19z100nkSTALL` | ~88,900 | | ~70 | stalled, in the 4-cycle on the old code |

- `production.env` holds the production settings of z100nk; machine paths are removed.
- `run.env` holds the mesh and run settings.
- These checkpoints have no controller state, so every restart logs "no solver-controller state", and its first step fails once. That is expected.
- To restart, run `julia -t 4 --project=<env> julia/popiii_m2c.jl lines 128` with the env vars `RESTART=m2c_z19z100nkS6`, `CKDIR=<dir with the files>`, `RUNREC`/`POPIII_OUT=<run dir>` and `EQEOS_CACHE=<file>`. Launchers: `scripts/popiii_run/iceulrun.sh`.

## Open items, after the stall
- **Combined-code validation:**
  - 300-step restart;
  - a start from the IC with `ICMLOG=1e-8` (core formation);
  - Noh, SOD, SUBC and SUPC.
  - Then a new 128-zone run, z128m2: `ICMLOG=1e-8` from the IC.
- **ALI/EP factor:** make the full ALI line-cooling solution smooth and repeatable from call to call. It changed 12–190% between calls in some zones on lagfix. Only then reintroduce the ALI model.
- **Energy:** the remaining total-energy gap at the accretion shock (G ~1e-2, zones 27–34) comes from mass crossing a moving mesh with an upwinded interface density. Tom accepted it for now as a shock-capturing error. The follow-up idea is to let the mesh follow the gas through the shock zones.
- **Queued:**
  - A7/A8: delete CPL_BAD and make the PSW regime switch smooth (G2BLEND).
  - A11.
  - B11: the DTERR estimator through the same operator as the solve.
  - Delete SUBZONE and make SLABFLUX/SLABTIE the default.
  - Remove the dual-energy switches (step 4) once convergence is at least as good.
  - Possibly add a Noh regression.
- **Known:** the F_r(101) Jacobian mismatch, now the current task. With `CLOSJ`, about 2.7e-4 remains in the photosphere F_r rows, from the undifferentiated r dependence of the sweep geometry and from columns outside the ±1 window.

## Explorer
The Tom-facing run explorer is the artifact https://claude.ai/artifact/EiGp93oWTQ9KAkhwKRAMsZ. Its build scripts lived in the old session's scratchpad and are not in the repo.
