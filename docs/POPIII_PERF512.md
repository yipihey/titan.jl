# PERF512: where the 512-zone fixed-mesh Pop III run spends its time, and what was made faster

Branch `perf512` (from `popiii-s6` 79d4dc1), 15 Nov 2026. All measurements are short restarts of the z512p
configuration (phase B from the 0.05 M☉ checkpoint `z19f128c_A3_1`, remapped by `envbuild!` onto the 508-zone
logarithmic mesh, Δln r = 0.049, 0.1 R☉ to the cloud edge; z512p's own `phaseB.env` + `run.env`, i.e. COUPLED, DTERR,
NLTE=hyb, LINEFAST, HOT, …), 300 steps (jstep 27710 → 28010), each from a `git archive` tree in its own run directory
`D/t_perf<N>` (D = /zpool/nvme/data/tabel_scratch/titan_popiii/runs), `numactl --physcpubind=96-99|100-103|96-103
--membind=12`. Tools: `TIMING=1` (perf.csv: wall time per 100 steps of the hydro+radiation step, chemistry, H₂ lines,
GC and allocation), `CPLPROF=1` (wall time and call counts of matgen, matslv, maxdel, update, peq; coupled-zone
Jacobian count), the `Profile` stdlib (`CPLPROFILE=27760 CPLPROFILEN=200`: flat and call-tree dumps,
`cplprofile_*.txt` in the run directory), the ALI/EP call times of the log. Scripts: /zpool/nvme/data/tabel_scratch/titan_perf512
(`prun.sh` launcher, `cmpstate.jl` state comparison, `bench.jl` harness).

Bit-identity test: the final checkpoint of each run (all 591 TITAN state arrays, `cmpstate.jl`, and the serialized material /
driver state `.aux`, byte compare) against the unmodified `popiii-s6` tree on the same checkpoint and the same 300 steps
(t_perf2). The 300 steps contain 9 Newton failures, 1 error rejection, 26 ALI calls and 6 checkpoints-worth of
bookkeeping, so the failure/retry paths are exercised.

## 1. The observed "2.4 s/step hydro+radiation" is mostly the first report window

perf.csv windows of the unmodified code (t_perf2, 4 threads, MATTHREADS=4) — seconds per 100 steps:

| window (jstep) | step (hydro+rad) | chemistry | H₂ lines | GC | alloc |
|---|---|---|---|---|---|
| 27710–27800 (90 steps, incl. JIT compile and 3 Newton failures of the remap) | 234.9 | 10.1 | 96.4 | 11.7 | 108 GB |
| 27800–27900 | 123.0 | 4.4 | 107.1 | 9.3 | 90 GB |
| 27900–28000 | 99.8 | 4.3 | 59.4 | 7.7 | 75 GB |

The 235.6 s/100 steps of z512p's log is the first window after the start: ~110–130 s of it are JIT compilation
(first window 235 s vs 100–123 s afterwards; the same first-window excess appears at 128 zones: 141 s vs 26 s) plus the
three diverging attempts of the remapped first step (3 × 25 Newton iterations). Its second window was 123.7 s, i.e.
1.0–1.2 s/step of hydro+radiation, plus 0.6–1.1 s/step of H₂ lines: a total of 1.6–2.3 s/step.

## 2. Breakdown of a step (unmodified code, 512 zones, window 27900–28000)

Per accepted step: one TR-BDF2 step = TR stage Newton + BDF stage Newton (each iterations + one "boost" iteration) +
error estimate (3 extra matgen+matslv: estimate!, refilter!) — on average **14.6 matgen, 11.6 Newton updates** per step
(1456 / 1156 per 100 steps). Newton failures: 9 in 300 steps (each 25 iterations, ~6 % of all iterations), 1 error
rejection — the step controller wastes little.

| routine (CPLPROF) | calls/step | ms/call | s/100 steps | per zone per call |
|---|---|---|---|---|
| matgen (residual + coloured forward-AD Jacobian, material in Duals, MATTHREADS) | 14.6 | 19.6 | 28.6 | 39 µs |
| matslv (BlockSolver.Penta, 9×9 static blocks) | 14.6 | 1.7 | 2.5 | 3.4 µs |
| update (correction, **Float64 EOS/opacity fill**, derived arrays) | 11.6 | 36.8 | 42.5 | 72 µs |
| peq (incl. the coupled species' local Newton `cpl_local_all!`, EOS/opacity, eddfac) | 1 | 160 | 16.0 | — |
| rest (cpl_prepare!, predict!, estimate bookkeeping, save/restore) | | | ~10 | |

Profile of the step (main thread, 200 steps):
- **fill_eos! of update!/peq: 47 % of the step.** For the ~400 COUPLED zones the Float64 EOS is `cpl_zone`, a full CK
  expansion of the reduced network (`ChemNet.ck_eval`: ~40 rate coefficients, each reverse rate an equilibrium constant
  `lnKeq` built from partition functions of up to 61 levels), evaluated **serially** at every Newton update. (The Jacobian's
  material sweep was already threaded; this Float64 refresh was not.)
- matgen 23 %: the residual's own 54-partial Dual arithmetic (serial), the threaded material sweep, `context()` (old-state
  nuclear source `source_old!`, MLT `convcoef`, serial).
- `cpl_local_all!` 10 % (per coupled zone a 2×2 Newton of its species rows with ForwardDiff, serial).
- eddfac (formal solution over tangent rays, O(N²)): 22 ms/call, 2 calls/step, ~3 % of the unmodified step. matslv 5 %. cpl_prepare! 3 %.

Answers to the specific questions:
- **Jacobian**: analytic forward-mode AD (ForwardDiff), stencil-coloured with 6 colours × neqn 9 = 54 partials; no finite
  differences. Material per zone per matgen: one evaluation in 54-wide Duals for the table/network zones; for the coupled
  zones the 7 local outputs (p, e, S_H⁺, S_H₂, Λ, cq, yq) come from a cached local 4-partial Jacobian, recomputed only when
  (ρ, T) moved by > 3e-3 or x by > 0.1 (CPLJLAG): 29 local Jacobians per matgen of ~400 coupled zones (93 % reuse,
  2.2 M reuses per 100 steps), and lifted into the 54-wide Duals by the chain rule.
- Per Newton iteration and coupled zone: one Float64 CK expansion (update), ~0.07 local Jacobians, two cache-hit lifts.
  Per step and coupled zone: up to 12 local 2×2 Newton iterations (2 CK expansions each) in `cpl_local_all!`, 1 in
  cpl_prepare!, 1 in cpl_finish!.
- Block-tridiagonal/pentadiagonal solve: 1.7 ms per call, O(N), 5 %. Not a target.
- Table lookups (PopIIIEqEOS, StellarOpacity): inside the threaded sweeps, cheap.
- H₂ lines: EP call 0.21 s every step, ALI call 7.8 s (2 ALI iterations) every ~11 steps (26 calls / 300 steps).
  EP: escape-probability static integral (42 E₂ per line and zone) and two serial loops; ALI: the LINEFAST formal solution
  over n tangent rays × n zones × bins per thick line — O(N²) — threaded over lines (scales 3.96× from 1 to 4 threads).
- Output/snapshots: 0.3–0.4 s per 100 steps.

## 3. Scaling with the number of zones (unmodified code, 4 threads; the same restart on 128, 254 and 508-zone log meshes)

| quantity | 128 | 254 | 508 | exponent |
|---|---|---|---|---|
| matgen ms/call | 4.6 | 9.4 | 19.6 | N^1.05 |
| update ms/call | 8.8 | 18.1 | 36.8 | N^1.03 |
| matslv ms/call | 0.43 | 0.89 | 1.7 | N^1.0 |
| peq ms/step | 38 | 75 | 160 | N^1.04 |
| hydro+rad step s/100 steps (window 3) | 26.3 | 48.8 | 99.8 | N^0.97 |
| EP call s | 0.045 | 0.11 | 0.21 | N^1.1 |
| ALI call s | 0.60 | 2.06 | 7.8 | **N^1.85** |
| ALI calls per 300 steps | 66 | 50 | 26 | |
| Newton updates per step | 13.0 | 11.5 | 11.6 | |

The hydro step is O(N): time per zone per Newton iteration is constant (matgen ~37–39 µs, update ~70 µs per zone).
The only superlinear part is the ALI formal solution (O(N²) over tangent rays; eddfac's ray work is O(N²) too:
22 ms/call at 512, ~3 % of the unmodified step, ~9 % of the new one). At 512 zones the ALI is called less often (the adaptive interval grows), which keeps lines/step
~O(N^1.5) in these 300 steps.

## 4. Thread scaling (unmodified code, 512 zones, s/100 steps, window 27900–28000)

| threads (julia -t = MATTHREADS) | step | lines | ALI call | GC |
|---|---|---|---|---|
| 1 | 128.0 | 214.6 | 30.5 s | 12.5 |
| 4 | 99.8 | 59.4 | 7.8 s | 7.7 |
| 8 | 92.2 | 32.1 | 4.0 s | 7.2 |

The lines (ALI, EP population solves, escape probabilities) parallelize almost perfectly; the hydro step only 1.28× from
1 to 4 threads, because the large serial parts above (fill_eos!, cpl_local_all!, context, the residual) did not run on
the threads at all.

## 5. Implemented (pure performance, bit-identical)

All MATTHREADS-threaded per-zone loops go through `TitanResidFast.zone_sweep!(f, ks, ok)`: zones from a shared atomic
counter, the serial loop's semantics kept exactly (ZONE[] left at the last zone of the serial order, or restored where the
serial loop did not set it; the exception of the lowest failing zone rethrown after all tasks ended), serial when
MATTHREADS ≤ 1 or already inside a threaded sweep.

| commit | change | effect |
|---|---|---|
| dbe5512 | `fill_eos!`/`fill_opac!` of NetGas threaded (the coupled zones' Float64 CK expansions of every update!/peq); `cpl_local_all!` threaded (results applied together, as before); cpl_prepare!'s coupled-zone values precomputed on threads (the serial loop then reads them from the zone cache; failures left to the serial loop) | update 36.8 → 12.4 ms, peq 160 → 67 ms; step 99.8 → 61.5 s/100 |
| afd300f | `ChemNet.lnKeq` (Float64 T): per-thread memo of the species terms ln(Z_tr Z_int) − E0/kT at the last T (same term, same summation order); EP line-luminosity loops zone-outer (same recursion, cache-friendly); `incident()` Jbg behind a function barrier (boxed captures) | lnKeq 2.5 µs → 0.18 µs; step 61.5 → 55.2, lines −4 s/100 |
| 51a2dd6 | `source_old!` and the residual's coupled-zone rate loop threaded | small (≤ 1 s/100) |
| 7868862 | `convcoef`'s MLT optical-depth opacities precomputed on threads for materials whose Float64 opacity has no side effects (`opac_pure`, NetGas unless NLOIT), with the serial loop's own (stale) ZONE[] | matgen 18.6 → 15.3 ms/call; step 55.2 → 49.3 s/100 (with 51a2dd6) |
| ff8a374 | `CPLPROFILEN` (profile length) and a call-tree dump beside the flat profile | diagnostic only |

Bit-identity: t_perf3 (dbe5512), t_perf5 (afd300f), t_perf6 (51a2dd6), t_perf7 (51a2dd6 + EPTAB=1), t_perf15 (7868862): all 591 state
arrays and the .aux identical to t_perf2 (unmodified popiii-s6) after 300 steps; also at 128 and 254 zones (t_perf10 vs
t_perf8, t_perf11 vs t_perf9), and unmodified -t 1 vs -t 4 (t_perf4 vs t_perf2).

Measured (4 threads, 512 zones, s/100 steps, window 27900–28000):

| | step | chem | lines | total | s/step |
|---|---|---|---|---|---|
| popiii-s6 (t_perf2) | 99.8 | 4.3 | 59.4 | 163.9 | 1.64 |
| perf512 (afd300f, t_perf5) | 55.2 | 3.6 | 55.4 | 114.5 | 1.15 |
| perf512 + EPTAB=1 (t_perf7) | 53.3 | 3.7 | 46.9 | 104.2 | **1.04** |
| perf512 final (7868862, t_perf15) | 49.3 | 3.6 | 54.4 | 107.6 | 1.08 |

Hydro+radiation step **2.0×** faster at 512 zones (99.8 → 49.3 s/100; per call: update 36.8 → 9.3 ms, peq 160 → 57 ms, matgen 19.6 → 15.3 ms) (also at 128 zones: 26.3 → 14.7 s/100; 254: 48.8 → 24.8), whole run 1.4–1.6×.
`EPTAB=1` is an existing switch (Nov 2026, documented "the same to the last bit"; confirmed here: t_perf7 identical) that
z512p's run-directory phaseB.env lacked although the repository's phaseB.env sets it: EP call 0.21 → 0.11 s.
**Recommendation: set EPTAB=1 in the run env** (no code change needed).

Threads with the new code (51a2dd6, EPTAB off; s/100 steps, window 27900–28000):

| | step | lines | ALI call | total s/step |
|---|---|---|---|---|
| popiii-s6, 4 threads (t_perf2) | 99.8 | 59.4 | 7.8 s | 1.64 |
| popiii-s6, 8 threads (t_perf12) | 92.2 | 32.1 | 4.0 s | 1.29 |
| perf512, 4 threads (t_perf6, profiled) | 57.5 | 58.3 | 7.7 s | 1.20 |
| perf512, 8 threads (t_perf13) | 40.8 | 27.9 | 4.0 s | **0.72** |

With the serial parts on threads the step now gains 1.4× from 4 to 8 threads (the unmodified code 1.08×); 8 threads +
EPTAB=1 should be ~0.65 s/step, 2.5× the production configuration. All of these runs are bit-identical to t_perf2.

## 6. Ranked further opportunities (expected gain × risk)

Pure performance (same numbers):
1. **More threads for the lines** (−t 8 / MATTHREADS=8): ALI scales linearly with threads; after the changes above the hydro
   step also profits. Risk none. Gain see §5 (8 threads).
2. **CK expansion cost** (ChemNet `ck_eval`/`ck_expand`/`dxdt_k`): still the largest CPU consumer of the step (now on threads).
   Untyped `Reaction.k`/`Species.Zint` Function fields (dynamic dispatch per reaction), per-call vectors (`ks`, `x`, `phi`,
   `fullvec`, `dd`), NamedTuple merges for tl3: a typed, allocation-free rewrite with the same arithmetic order would cut it
   2–5× and most of the 0.7 GB/step allocation (GC 6.5–7.7 s/100 steps, stop-the-world). Medium effort, low risk if checked
   bit for bit as here.
3. **The residual's serial Dual arithmetic** (~1/2 of matgen after threading the material): thread the per-zone row
   loops of `residual!` (rows of zone k need only k±2), or evaluate with fewer partials per zone. Large effort.
4. **eddfac**: 22 ms per call at 508 zones (bench), 2 calls per step (peq and the BDF predictor) = ~9 % of the new hydro
   step, O(N²). The inward and outward sweeps compute the same `exp(-dτ)` and path lengths for each (ray, segment), and the
   geometry (3 sqrt per segment) is fixed on the Eulerian mesh: caching them halves it bit for bit; the rays are independent
   and can be threaded (ray-major sweep into a (ray × node) buffer, moments summed afterwards in the same order).
5. JIT compilation (~110–130 s per process start, the first-window excess): a sysimage or precompile workload. Only matters for
   short segments.
6. Type instabilities: `PM.HOT`, `TX.CPLW`, `CPLCACHE`, `CPL_ZONE`, `OPAC` are `Ref{Any}`/`Vector{Any}`; the coupled-zone path
   returns through `Any` (heap-allocated 54-wide Duals). Per-zone dynamic dispatch, a few % of matgen. Measured (`bench.jl`, `code_typed`): `material_gas0` and `material_opac`
   of NetGas infer `Tuple{Any,Any,Any}`, `cpl_g` returns `Any`; `eddfac!`, `update!`, `Penta.matslv!` are type-stable.
   Allocation per call (bench, end-of-run state): matgen 54 MB with fresh coupled-zone Jacobians, 17 MB cached; update
   17 MB; peq 73 MB; cpl_prepare! 35 MB — the CK expansions and the Any-returning per-zone tuples.

Changes of the numerics (not adopted; would need a default-off switch and a physics check):
7. **ALI cadence** (ALICAP / ALIGROW / ALISHRINK, ntol of the ALI iteration): the ALI is the single largest item after
   the changes (≈ 0.4–0.7 s/step at 4 threads, O(N²)). Halving its frequency saves ~0.2–0.3 s/step.
8. `LINEXC2=16` (resonance window 16 instead of 38 thermal widths²; documented ≲ 3e-9 change of the cooling): fewer bins
   per segment in the ALI and EP.
9. CPLJLAG (coupled-zone Jacobian reuse) wider; CONV / NITER (Newton tolerance, 11.6 iterations per step); the error
   estimator's 3 extra matgens per step (20 % of matgen) with the first filter only (DTERRFILT=1).
10. Fewer colours (5 instead of 6) in the AD Jacobian (TitanResidFast note: aliases a 1e-8 term).
11. Note (latent): `convcoef`'s optical-depth loop calls `material_opac` without setting ZONE[]; for NetGas with NLOFIX the
    electron-scattering floor then uses the composition of the last zone of the preceding sweep for every k. Kept as is
    (bit-identity); setting ZONE[] = k there changes τ for the MLT weight slightly.

## 7. PERF-RED (branch `perf-red` from popiii-s6 c69c49c): the abn02red from-IC run

abnred8 (IC=abn02red, ICEUL, 300 network zones, full COUPLED, DTERR, NLOPAC, VEFIMP, LINEIMP, NLOIT, LOCALJ, MATTHREADS=8) cost
62–76 s per step with 150–230 GB allocated per step. Allocation profile (`CPLPROFALLOC=1e-4` over steps 1–3, t_perfr1): **97 % of
the bytes come from `ChemNet.photorates`, the direct frequency quadrature (n = 200 per level)**, called by `photorates_lj` (LOCALJ)
whenever T_c of the local diluted blackbody is outside the tables (log T_c < 3: no star yet at the start, every zone), i.e. in every
CK expansion's photo-rates and in every `nlo_comp` of the NLOIT opacity (whose photo-rates are uncached, 3–6 evaluations per
Jacobian refresh). The loop boxed every intermediate because the cross section σ(ν) (an untyped field of `Photo`) and J(ν) were
called through untyped values. Per call (isolated test): 132 kB and 120 µs → 4 kB and 9.7 µs.

Fix (30ad384): the quadrature behind a function barrier (`_photoRA`), the same arithmetic in the same order; bit-identical
(360 old/new comparisons incl. Dual T; 6-step ICEUL start t_perfr2 vs c69c49c t_perfr3: all 591 state arrays and the .aux identical).

| steps 10–20 (8 threads) | step s/step | alloc GB/step | GC s/step |
|---|---|---|---|
| c69c49c (abnred8) | 75.6 | 229 | 15.3 |
| perf-red (t_perfr5) | 9.9 | 10.4 | 0.8 |

Remaining: matgen 0.49 s and 0.49 GB per call (127 per 10 steps); `_nlmf` (NLOIT non-LTE means) ~1.2 GB/step; Base array allocations
~4 GB/step of unresolved origin. Changes of the numerics (report only): extend the photorates_lj tables to T_c < 1e3 K (or skip the
local-J term where W B(T_c) ≪ the CMB), cache NLOIT's photo-rates per zone and T.
