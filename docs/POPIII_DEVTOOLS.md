# Dev loop tools (branch devtools, 24 Nov 2026)

Scripts in `scripts/dev/`. They extend `scripts/proto/gate_run.sh` (the runner; one added hook, `GATE_PRINT`) and the P5 harness
(`prod_env/{build_prod.sh,run_test.sh,cmp.sh,metrics.py}`); they do not replace them.

## mkimage.sh `<worktree> [grey|mg|name]`
Builds the TitanStack image (docs/POPIII_PRECOMPILE.md) of a tree into `$DATA/images/<branch>-<sha>[-dirty-<hash>][-mg]/{env,depot}` and prints
`export JLENV=... JULIA_DEPOT_PATH=... TITAN_PKG=1` (`eval "$(mkimage.sh WT)"`). The env is a copy of `prod_env/env` whose Manifest points at the tree.
- Reuse: same name and `Base.isprecompiled` true = no-op (about 0.3 s). A dirty `julia/` (diff or untracked files) is built under `-dirty-<hash>`.
- `mg` is a separate image: `MGROUPS` is a load-time key (conditional JuliaSolver include), so an MGROUPS=4 run on the grey image falls back to the include path.
- A tree that is not a git checkout (a `git archive` scratch tree) carries `.titan_sha` and is named `tree-<sha>`. The gate runs use such trees, so a
  worktree that moves on (popiii-proto did, during this work) cannot invalidate an image.
- Build wall: 440-450 s for two images built in parallel on nodes 3 and 4 (about 350 s of `build.jl` plus about 93 s of dependency precompile in the fresh depot).

## cklib (`$DATA/cklib/`, README there)
Copies of the checkpoints and the run.env/phaseB.env/launch.env of zie1 and zie3, laid out for `gate_run.sh` (`REFENV=cklib CKSRC=cklib/<run>/ck RESTART=<name>`):
zie1 `z19zie1S6_auto1` (23300), `c5_dA` (23420), `z19zie1S16` (63334); zie3 `fp3a1` (5740), `z19zie3S8` (32000), `z19zie3S9` (36000). The IC needs no file.

## gates.sh `<worktree> [--nodes 3,4] [--only 1,2,3,4,5] [--ref COMMIT|DIR] [--include-path] [--nocache]`
The five default merge gates of POPIII_PROTO_GATE.md against the previous gated commit (default: the tree's first parent, as a `git archive` scratch tree).

| gate | jobs (run label) | pass rule |
|---|---|---|
| 1 | zie1 S6_auto1, 5 steps, `NKCHECK=2 NKCHECKN=1 NKROWPROBE=0`, grey (`nkg_5`) and `MGROUPS=4` (`nkm_5`) | 2-norm rel (eps x1) at most 1.15 x the documented 1.5e-10 / 2.306e-7 |
| 2 | zie1 IC, 40 steps, MGROUPS unset, this tree and the reference (`ic_40`) | step lines (cmp.sh filter) and decompressed tstate identical |
| 3 | zie1 S16 and zie3 S8, 50 steps, `ETOTMON=1` (`s16_50`, `s8_50`) | completes, 50 steps, identical to the reference commit's run |
| 4 | the same two windows | ledger imb/thr (ledger.py of proto-gate) at most 1.05 x the reference commit's |
| 5 | zie3 fp3a1, 40 steps, `MGROUPS=4` (`mgz3_40`) | 40 accepted, dt median at least 0.8 x 9e3 s |

Jobs go to the nodes statically (heaviest first), one julia process per node at a time, 8 threads on cpus 8n..8n+7. The reference runs (gates 2-4) are cached in
`$DATA/dev/gcache/<tree name>/<label>` and reused while the reference commit stays the same; `--nocache` and `--include-path` use a fresh directory. Output: one table
(gate, check, value, reference, PASS/FAIL), `<out>/summary.txt`, `<out>/gates.log`. Gate 1/5 references are the documented constants, not re-run.
Check that no run prints "falling back to the include path" (the image's baked keys differ): that would show as a slow run, not as a failure.

## Several runs in one julia process: not byte-identical
`inproc.jl` + `inproc_test.sh` run `PopIIIRun.main(ARGS)` repeatedly in one process (environment reset to the start state per job, stack defaults re-applied).
Test on 3eb90bc with the TitanStack image: job a1 = zie1 IC, 3 steps; b2 = zie3 fp3a1 restart, 20 steps; c3 = a1 again; against three separate processes.
- a1 (first in the process): step lines and tstate byte-identical.
- b2: tstate byte-identical, but one step line differs (diagnostic fields: 2.08 vs 2.06, 0.2 vs 0.3, 40 vs 44, +6.29e-04 vs +3.84e-05), so state carried over from a1 touches the solver.
- c3: NOT identical: the second IC run does not start at step 0; it begins at jstep 5740+ (b2's step counter) and ends after 2 step lines ("run ended (3, :maxsteps); jstep 5743").
Leaking driver/solver globals seen in the output: TitanCore scalars (`jstep`, `timen`, `tmax`: the IC path does not reset them), the cumulative counters printed at the end
("ALI line transfer: 5763 calls" = `NALICALL`, the NKRYLOV summary, the DTERR statistics dict with `steps 5763`), `atexit` summaries registered once per call
(printed three times at exit), and whatever cached solver state produces the differing step-line diagnostics of b2 (module-level Refs/caches of `PopIIIRun`/`TitanCore`/`NewtonKrylov`:
the `const ... = Ref(...)` set of PopIIIRun.jl, `TX.*` caches). Not fixed. `gates.sh --inproc` therefore exits with an explanation unless `GATES_ALLOW_INPROC=1`.
Timing note (shared box, load 25-75): in-process 483 s for the three jobs against 902 s separate; the saving is the repeated start-up (image load ~1 s plus IC and compile).

## First measurements (3eb90bc, ref 6875c7e, nodes 3 and 4, box load 25-75)
- Image builds: grey and mg in parallel 448 s / 447 s (a single image alone: 455 s); second call: 0.3 s no-op.
- `gates.sh --nodes 3,4`, all five gates, reference runs included: 2600 s wall (job walls with the image: ic_40 264 s, s16_50 394 s, s8_50 533 s, nkm_5 660 s, mgz3_40 1406 s, the three reference jobs 263/400/about 530 s).
  A second call reuses the cached reference runs (only the changed jobs run): 1453 s here, because gate 5's MG job was re-run with `PRINTEVERY=1` (perf.csv rows).
- Same gates serially on the include path (estimate): the image saves about 200 s per job (measured: ic_40 467 s include path vs 264 s image), so about 9 jobs x 200 s = 1800 s more,
  and sum of job walls about 4800 s with the image, i.e. about 6600 s serial on the include path against 2600 s (first run) and 1453 s (cached reference).
- Table: gate 1 grey 1.537e-10 and MG 2.306e-07 PASS; gate 2 identical (17 step lines, tstate 10530280 B); gate 3 both restarts 50 steps, identical to 6875c7e; gate 4 ETOTMON 2.127e-03 / 1.150e-03, equal to 6875c7e;
  gate 5 40 accepted, dt median 9.195e+03 s, 1 tr failure: all PASS.
