# titan.jl run-time switch audit 2: what the clean code base keeps (branch `cleanup-audit`, 16 Nov 2026)

Second audit of every ENV switch, made for a clean code base that drops the earlier experiments. The first audit is
`docs/POPIII_SWITCH_INVENTORY.md` (14 Nov, 326df84). This one starts from `popiii-s6` at **1782928** (GRIDEQ merged) and also covers
`g2-domain` (00a21f1, not merged). The IonFS branches (`ionfs-implicit`, `ionfs-v2`) are still being worked on, so they are only
listed. Stage 1 of the clean-up (§8) is implemented on this branch. Nothing is merged into `popiii-s6`.

> **Status, 16 Nov evening.** The PI decided the open questions of §9. Stages 2 and 3 and two result-changing adoptions are now
> implemented (§10): the fmove mesh bound in DTERR, and G2BLEND + H2PF as the defaults. The configuration that remains is described
> in docs/POPIII_STACK.md. The tables in §2–§8 describe the code at 1782928 and stage 1.

## 0. Scope, reference stacks, method

- **Switches.** Every `get(ENV, …)`, `ENV[…]` and `haskey(ENV, …)` in `julia/*.jl` and `julia/z0/*.jl`, plus GRIDEQ's `geqv(...)`
  (GEQC, GEQW, GEQX) and the `${X:-…}` variables of `scripts/popiii_run/*.sh`. That gives **465 names**.
  **345 of them are on the production path**: `popiii_m2c.jl` and everything it includes. 14 of those 345 are shell variables of
  the launch scripts. The other 120 names are read only by test, benchmark or old-stage scripts (§3).
- **Reference stacks** (read-only; later entries win, composed as `iceulrun.sh` composes them: `phaseB.env`, then COMMON, then
  `run.env`):
  - **z300mm** (MAIN): taken from the recorded `env` lines of `runs/z300mm/meta_000006450_1524405.txt` (segment S3, 117 keys).
    GEQC and GEQW come from `run.env`, because `recenvkeys` does not see `geqv` reads.
  - **z300g2** (fixed-mesh reference): the same stack without `GRIDEQ`, `GEQRAMP`, `GEQC`, `GEQW` and `FLUSHOUT`.
  - **z100mm** (100-zone from-IC attempt): the z300 stack with a 100-zone mesh (`ENVEUL=0,0.2326,…`), started with `ICEUL=1`.
    Its last try ends with `GEQRAMP=0 GRIDEQ=0`. All three tries failed at step 0.
- **Line numbers** in the tables are at 1782928 (before stage 1).
- **Verdicts.** **KEEP-STACK**: used by the stack and needed. **KEEP-TOOL**: diagnostics, I/O or performance worth keeping.
  **MAKE-DEFAULT**: a general method the stack always uses; it becomes the default and the switch goes. **REMOVE-SUPERSEDED**: an
  older patch that a general method replaced; the note names the method. **REMOVE-DEAD**: an experiment never adopted, or a path
  no longer reachable. **UNCLEAR**: the note says what decides it.

### 0.1 Findings about the stack itself (no code change needed to see them)

1. **`CONV=1e-3` in every run.env has no effect.** `popiii_m2c.jl` never reads CONV; only `order_test`, `popiii_m1`, `s5_run` and
   `z0_jac_check` do. The z300mm Newton converges at the checkpoint's `conv = 1e-5` (read from `m2c_z19z300mmS2.tstate.gz`). If
   1e-3 was intended, it was never applied.
2. **Three run.env keys are read by nothing:** `FLUXTIECORE=1`, `NLOPACJTAU=1,10` and `PSWDEEP=1e3`. Their code was removed earlier.
3. **With `DTERR=0.1`, all of TITAN's Δt rules are diagnostic only** (`Trbdf2.jl:198`, `step_err!`). That covers `timstp!` with
   STOL and CDYN, and every `DT_HOOK`: CPLXSTOL, DTRADMAX, DTENV and the GRIDEQ fmove cap (GEQC[5]). It also covers the
   line-cooling cap and STOLCHEM. So the stack's `CDYN=3e5 STOL=0.3 FAILTFAC=0.5 CPLXSTOL=3 STOLCHEM=1 DTRADMAX=0 ERTAUMIN=3000`
   do not change the run. There is one exception. The HOT nuclear cap is also a DT_HOOK, but `StarMaterial.burn!` applies it to
   `dtime` directly (`StarMaterial.jl:359`), after `step_err!` has set Δt, so it still acts. The GRIDEQ fmove cap, on the other
   hand, is inactive in z300mm. That may matter for the moving mesh (§9).
4. **With `LINEIMP=1`, the IMPCOOL slope (a second transfer call), DTCOOL, TCOOLTH and the t_cool cap are not used.**
5. **With `NLOIT=1`, NLOIMPT, NLODAMP, CPLNLOX, NLOFIX and NLOTMAX have no effect in network zones** (NLOIT docstring). The stack
   still sets `NLOIMPT=1 NLOFIX=1`.
6. **Of the table/network rules, only PSWNH is active:** `PSWITCH=1e30 PSWTAU=0 PSWSTAY=0`. PSWNH uses PSWHYS and PSWTMIN.
   PSWTEX=5000 is read only through the two rules that are off. `PSWRAMP=10` acts only when a restart moves the table edge.
7. **`XADV2` is used although no stack sets it.** Its default is the value of SL3, which is 1.
8. **`CPLBAND=` (empty) couples every network zone.** The split chemistry patches (CHEMSUB, CHEMCKHSS, ECHEMSPLIT=0, BALSPLIT=0)
   then act only where a zone is not coupled.
9. **The 1782928 code reproduces production bit for bit.** A restart of `z19z300mmS2` (step 6400) for 100 steps on a `git archive`
   of 1782928 gave exactly z300mm's own `S3_auto1` at step 6500: 591/591 state arrays and the .aux identical (§8). The z300mm meta
   says `git_dirty_julia yes`, but the dirty files did not change the result.

## 1. Counts

| Verdict | Production-path switches | of which removed in stage 1 |
|---|---|---|
| KEEP-STACK | 83 | – |
| KEEP-TOOL | 82 | – |
| MAKE-DEFAULT | 33 | – |
| REMOVE-SUPERSEDED | 66 | – |
| REMOVE-DEAD | 51 | 44 |
| UNCLEAR | 33 | – |
| total with a verdict | 348 (345 read at 1782928, plus H2PF and G2BLEND from g2-domain, plus CONV) | 44 |

The REMOVE-SUPERSEDED group retires against these general methods: **DTERR** (G3: CDYN, STOL, FAILTFAC, CPLXSTOL, STOLCHEM, DTRADMAX,
ERTAUMIN, ERTHICK, DTENV, CONV); **LINEIMP** (IMPCOOL, IMPALWAYS, DTCOOL, TCOOLTH, LINEFRAC); **NLOIT** (NLOIMPT, NLOFIX, NLODAMP,
CPLNLOX, NLOTMAX); **G2 dense-interior table** (PSWITCH, PSWTAU, PSWSTAY, PSWTEX, CHEMFBDENSE, CHEMTMIN; and with G2BLEND also
PSWNH, PSWHYS, PSWTMIN); **global COUPLED** (CPLBAND, CHEMSUB, CHEMSUBMAX, CHEMCKHSS, ECHEMSPLIT, BALSPLIT); **VEFIMP**
(ERTAUCONV); **ICEUL + GRIDEQ** (G7: the 15 REZ* switches, RADFROM, REFINE, MSTOP, ENVMS0, ENVDT0, ENVSPLIT, ENVMERGE,
ENVSINKTAU, TRESET, ETOTGRAVGRID); **ChemNet** (CHEM=rhfc02 with RHOEQ and EQEOS); **SL3** (SLH); **AVZONE** (AVQL1).

## 2. Every switch on the production path (verdicts)

Columns: default (code at 1782928); value in the z300mm stack (`–` = not set; differences of z300g2 / z100mm in brackets); file:line of the read at 1782928; verdict; note. **S1** = removed in stage 1 (commit on `cleanup-audit`).


### I/O, run control, performance

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| SAVE | unset | `z19z300mmS3` | popiii_m2c.jl:2712, popiii_m2c.jl:2719 … | checkpoint name (end + every CKEVERY) | KEEP-TOOL |  |
| CKEVERY | `"500"` | `100` | popiii_m2c.jl:1933 | checkpoint interval (steps) | KEEP-TOOL |  |
| CKDIR | `joinpath(TF.ROOT, "l…` | `/zpool/nvme/data/tabel_s…` | popiii_m2c.jl:758 | checkpoint directory | KEEP-TOOL |  |
| RESTART | unset | `z19z300mmS2` | popiii_m2c.jl:1565, popiii_m2c.jl:1566 … | restart from a checkpoint (CKENV consistency check) | KEEP-TOOL |  |
| POPIII_OUT | `joinpath(TF.ROOT, "l…` | `/zpool/nvme/data/tabel_s…` | popiii_m2c.jl:2377, popiii_m2c.jl:2381 … | log/snapshot directory | KEEP-TOOL |  |
| RUNREC | `""` | `/zpool/nvme/data/tabel_s…` | popiii_m2c.jl:419 | run recording (snaps_*.jls, perf.csv, meta) | KEEP-TOOL |  |
| TAG | `""` | `_z300mm` | popiii_m2c.jl:413, popiii_m2c.jl:47 | log/snapshot file-name suffix | KEEP-TOOL |  |
| SNAPS | `"0"` | `0` | popiii_m2c.jl:546 | in-memory profile snapshots every 0.25 dex n_c | KEEP-TOOL | stack sets the default 0; drop from env |
| SNAPEVERY | `"0"` | `50` | popiii_m2c.jl:547 | snapshot every n steps (RUNREC) | KEEP-TOOL |  |
| PRINTEVERY | `"500"` | `100` | popiii_m2c.jl:1933 | status-line interval | KEEP-TOOL |  |
| MAXSTEPS | `"20000"` | `1000` | popiii_m2c.jl:2460 | steps per process (segment length; iceulrun.sh SEG) | KEEP-TOOL |  |
| DTSTOP | `"0"` | `1e-9` | popiii_m2c.jl:1932 | stop when Δt collapses below this (protects the checkpoints) | KEEP-TOOL |  |
| TITAN_NOLIB | `"0"` | `1` | TitanFortran.jl:100, TitanFortran.jl:105 … | pure-Julia backend, no libtitan | MAKE-DEFAULT | every Pop III run sets 1; make it the default of popiii_m2c (keep the Fortran path only for TitanDriver bit tests) |
| TITAN_MGR | unset | `600` | TitanFortran.jl:101, TitanFortran.jl:99 | grid array size (needs NOLIB) | KEEP-TOOL |  |
| TITAN_BITWISE | `"0"` | – | TitanDriver.jl:142 | bitwise comparison mode against the Fortran (FIX_CLOCK/FIX_ENERGY off) | KEEP-TOOL | only for the porting harness; keep while the Fortran parity tests exist |
| MSG | `"0"` | `1` | popiii_m2c.jl:387, popiii_m2c.jl:921 | Newton messages to stdout | KEEP-TOOL |  |
| DEBUG | `"0"` | – | popiii_m2c.jl:2716 | print the backtrace of an exception | KEEP-TOOL |  |
| FLUSHOUT | `"0"` | `1` (g2: –) | popiii_m2c.jl:586 | flush stdout after the hydro step | KEEP-TOOL |  |
| NANSCAN | `"0"` | – | popiii_m2c.jl:586 | scan the state for non-finite values after hydro and chemistry | KEEP-TOOL |  |
| TIMING | `"0"` | `1` | popiii_m2c.jl:1927 | per-step wall-time split (perf.csv) | KEEP-TOOL |  |
| DTDIAG | `"0"` | `1` | popiii_m2c.jl:2585, popiii_m2c.jl:2634 | which variable/zone limits Δt, on the status line | KEEP-TOOL |  |
| CHEMDIAG | `"0"` | – | popiii_m2c.jl:1927 | chemistry ΔT/Δu statistics per status line | KEEP-TOOL |  |
| MATTHREADS | `"0"` | `8` | TitanResidFast.jl:660 | threads of the per-zone material evaluation in the residual | KEEP-TOOL |  |
| CHEMTHREADS | `"0"` | `1` | PopIIIMaterial.jl:1256 | threaded per-zone chemistry | KEEP-TOOL |  |
| CHEMSCHED | `"zone"` | – | PopIIIMaterial.jl:1263 | zone scheduling of the threaded chemistry (zone\|block); results identical | KEEP-TOOL |  |
| LINESCHED | `"zone"` | – | H2Lines.jl:604 | zone scheduling of the threaded line transfer; results identical | KEEP-TOOL |  |
| ZONESWEEP | `"1"` | – | TitanResidFast.jl:872 | perf512: zone-sweep residual assembly (default on; =0 the old order, bit-identical) | KEEP-TOOL | default already on; the =0 branch can go once perf512 is accepted (stage 2) |
| LNKMEMO | `"1"` | – | ChemNet.jl:165 | fix-race: memoised ln k of the rates (default on; =0 the direct sum, bit-identical) | KEEP-TOOL | as ZONESWEEP |
| THCHAIN | `"1"` | – | popiii_m2c.jl:812 | threaded chain of the hot-branch EOS (default on) | KEEP-TOOL | as ZONESWEEP |
| CPLPROF | `"0"` | – | popiii_m2c.jl:2338 | profiling of matgen/matslv/... (CPLPROF), Profile samples and allocation reports (perf512) | KEEP-TOOL |  |
| CPLPROFEVERY | `"100"` | – | popiii_m2c.jl:2342 | profiling of matgen/matslv/... (CPLPROF), Profile samples and allocation reports (perf512) | KEEP-TOOL |  |
| CPLPROFILE | unset | – | popiii_m2c.jl:2351 | profiling of matgen/matslv/... (CPLPROF), Profile samples and allocation reports (perf512) | KEEP-TOOL |  |
| CPLPROFILEN | `"200"` | – | popiii_m2c.jl:2352 | profiling of matgen/matslv/... (CPLPROF), Profile samples and allocation reports (perf512) | KEEP-TOOL |  |
| CPLPROFALLOC | `"0"` | – | popiii_m2c.jl:2356 | profiling of matgen/matslv/... (CPLPROF), Profile samples and allocation reports (perf512) | KEEP-TOOL |  |
| LINEGPU | unset | – | h2lines_gpu.jl:45, popiii_m2c.jl:19 | H₂ line transfer on the GPU (h2lines_gpu.jl) | UNCLEAR | never used by a production run; keep only if a GPU line kernel is still wanted (A6000 shared) |
| LINEGPUFLOAT | `"32"` | – | h2lines_gpu.jl:35 | H₂ line transfer on the GPU (h2lines_gpu.jl) | UNCLEAR | never used by a production run; keep only if a GPU line kernel is still wanted (A6000 shared) |
| EQEOS_CACHE | unset | `/zpool/nvme/data/tabel_s…` | PopIIIEqEOS.jl:139 | cache file of the equilibrium (SCVH95 + He + x_e) table | KEEP-TOOL |  |
| POPIII_LIT | `joinpath(homedir(), …` | – | ContinuumOpacity.jl:71, H2Levels.jl:44 … | literature-data directory | KEEP-TOOL |  |
| NUCNET_REACLIB | unset | – | NucNet.jl:31, NucNet.jl:8 | REACLIB file of NucNet | KEEP-TOOL |  |
| ABNFINCSV | `joinpath(TF.ROOT, "l…` | – | popiii_m2c.jl:104 | digitised ABN02 Fig. 2 tables of the abn02fin/abn02red ICs | KEEP-TOOL |  |
| ABNREDCSV | `joinpath(TF.ROOT, "l…` | – | popiii_m2c.jl:105 | digitised ABN02 Fig. 2 tables of the abn02fin/abn02red ICs | KEEP-TOOL |  |
| S5DIR | `"/tmp/s5"` | – | popiii_m2c.jl:1112 | start from an S5 isolated-star state | UNCLEAR | not part of the from-IC line; keep only if the S5 star line stays (decides: PI) |
| STARFROM | unset | – | popiii_m2c.jl:1559, popiii_m2c.jl:1560 | start from an S5 isolated-star state | UNCLEAR | not part of the from-IC line; keep only if the S5 star line stays (decides: PI) |

### set-up, IC, physics model

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| IC | `"p100"` | `abn02red` | popiii_m2c.jl:413, popiii_m2c.jl:47 … | initial condition (abn02red in the stack) | KEEP-STACK | make abn02red the default; p100/abn02z19/abn02fin ICs: stage 3 |
| MCLOUD | `z19 ? "1e4" : fin ? …` | `7137` | popiii_m2c.jl:240 | cloud mass [M☉] | KEEP-STACK |  |
| TCMB | `"0"` | `54.5` | PopIIIMaterial.jl:1686, popiii_m2c.jl:1524 … | CMB temperature (cooling floor, line background, boundary) | KEEP-STACK |  |
| CHEM | `"rhfc02"` | `netmodernli` | popiii_m2c.jl:351 | material: RHFC02 PrimordialGas or ChemNet net* | MAKE-DEFAULT | stack always netmodernli; PrimordialGas branch (RHOEQ, EQEOS) goes with it (stage 3) |
| RHOEQ | `"5e-11"` | – | popiii_m2c.jl:353 | equilibrium switch / EOS of PrimordialGas (unused with CHEM=net*) | REMOVE-SUPERSEDED | superseded by ChemNet (CHEM=netmodernli) |
| EQEOS | `"scvh"` | – | popiii_m2c.jl:353 | equilibrium switch / EOS of PrimordialGas (unused with CHEM=net*) | REMOVE-SUPERSEDED | superseded by ChemNet (CHEM=netmodernli) |
| K5 | `"pss"` | – | popiii_m2c.jl:334 | choice of rate k5 of the RHFC02 network | UNCLEAR | only the RHFC02 network reads it? check ChemRates; with netmodernli remove |
| NLTE | `"1"` | `hyb` | popiii_m2c.jl:551 | H₂ line source function 1/0/ep/ali/hyb | MAKE-DEFAULT | stack hyb; ep/ali/1/0 branches are tests |
| NBIN | `"30"` | – | popiii_m2c.jl:752 | frequency bins of the line transfer | KEEP-STACK |  |
| RAD | `"0"` | `1` | popiii_m2c.jl:43 | radiation E_r, F_r on (RA04 CIE cooling off) | MAKE-DEFAULT | always 1 from the IC with ICEUL; RAD=0 + RADFROM = the old stage chain (stage 3) |
| HOT | `"0"` | `1` | popiii_m2c.jl:22 | stellar EOS, nuclear network and MLT in the hot core | KEEP-STACK | stack 1; make default with the HOTT1/HOTT2 blend kept |
| HOTT1 | `"3e5"` | – | popiii_m2c.jl:842 | blend window PopIII EOS → stellar EOS | KEEP-STACK | defaults used by the stack |
| HOTT2 | `"9e5"` | – | popiii_m2c.jl:842 | blend window PopIII EOS → stellar EOS | KEEP-STACK | defaults used by the stack |
| DH | `"2.5e-5"` | – | PopIIIMaterial.jl:1684, popiii_m2c.jl:828 | D/H for NucNet and HD chemistry | KEEP-STACK |  |
| HDCHEM | `"0"` | `1` | PopIIIMaterial.jl:1683 | D⁺/HD trace network and HD cooling | MAKE-DEFAULT | stack 1 |
| HECHEM | `"0"` | `1` | ChemNet.jl:74, PopIIIEqEOS.jl:42 | He, He⁺, He²⁺ network species; He table EOS | MAKE-DEFAULT | stack 1 |
| CHEM_YHE | `string(1 / 12)` | – | ChemNet.jl:32 | He/H by number (0 = pure-H benchmarks) | KEEP-STACK |  |
| H2PPD | `"0"` | – | ChemRates.jl:108 | H₂⁺ photodissociation by detailed balance | UNCLEAR | off in the stack; physics option, PI decides (default off = no H₂⁺ photodissociation) |
| STARLIGHT | `"0"` | `1` | popiii_m2c.jl:511 | star light below 13.6 eV in network zones (H⁻, LW) | MAKE-DEFAULT | stack 1 |
| SL3 | `"0"` | `1` | popiii_m2c.jl:1211, popiii_m2c.jl:358 … | three-level H atom (Peebles C, Sobolev Lyα) in network zones | MAKE-DEFAULT | stack 1; also sets the XADV2 default |
| SLHM | `"1"` | – | popiii_m2c.jl:513 | H⁻ photodetachment by star light (default 1) | KEEP-STACK |  |
| SLH | `"0"` | – | popiii_m2c.jl:514 | H(n≥2) photoionisation via the network ("energy-unsafe") | REMOVE-SUPERSEDED | superseded by SL3 (three-level atom); stage 2 |
| SLSHIELD | `"static"` | `sobolev` | popiii_m2c.jl:515 | LW shielding static/sobolev/none/off | MAKE-DEFAULT | stack sobolev |
| LYABB | `"0"` | `1` | PopIIIMaterial.jl:1107 | bound–bound energy exchange of the three-level atom (missing term) | MAKE-DEFAULT | stack 1: a correctness fix |
| TL3RECKT | `"0"` | – | PopIIIMaterial.jl:1013 | recombination kinetic energy ⟨hν⟩−χ₂ instead of 0.75 kT | UNCLEAR | old audit: a correctness fix (LTE consistency); off in the stack; PI decides MAKE-DEFAULT vs remove |
| TL3NEFULL | `"0"` | – | ChemNet.jl:846 | n_e of the three-level atom from charge conservation, not n_H x_H⁺ | UNCLEAR | correctness option never adopted; test with the stack, then default or remove |
| TL3DEEP | `"0"` | `1` | PopIIIMaterial.jl:1072 | G2: three-level atom also in the deep network zones | MAKE-DEFAULT | stack 1 |
| LINEINC | `"0"` | `1` | popiii_m2c.jl:566 | H₂ levels see the incident CMB + diluted photosphere | MAKE-DEFAULT | stack 1 |
| LINECOLL | `"0"` | `1` | popiii_m2c.jl:558 | gas term = net collisional exchange at ALI populations | MAKE-DEFAULT | stack 1 |
| LINESRC | `"linear"` | `zone` | H2Lines.jl:238 | zone-constant line source | MAKE-DEFAULT | stack zone |
| LINEQUAD | `"mu"` | `p2` | H2Lines.jl:242 | flux quadrature trapezoid in μ² | MAKE-DEFAULT | stack p2 |
| LINECORE | `"1"` | `8` | popiii_m2c.jl:658 | transparent merged-core zones for impact parameters | KEEP-STACK | stack 8 |
| LTAUCUT | `"10"` | – | popiii_m2c.jl:658 | H₂ line transfer only outside continuum τ (default 10) | KEEP-STACK |  |
| LTRAN | `"1"` | – | popiii_m2c.jl:910 | 1 = formal solution for f, g; 2 = Eddington approximation | UNCLEAR | stack default 1; =2 is a test path |
| EPGRAD | `"thermal"` | – | H2Lines.jl:602 | Sobolev gradient over one thermal width (default thermal) | KEEP-STACK | default already the fix; the other values are dead |
| DM | `"1"` | – | popiii_m2c.jl:163 | ABN02 dark-matter potential on (default 1) | KEEP-STACK |  |
| DMCORE | `ic == "abn02red" ? "…` | – | popiii_m2c.jl:177 | DM potential: cored inner profile (default on for abn02red) / which ABN02 DM line (z19\|final) | KEEP-STACK |  |
| DMLINE | `"z19"` | – | popiii_m2c.jl:114 | DM potential: cored inner profile (default on for abn02red) / which ABN02 DM line (z19\|final) | KEEP-STACK |  |
| XE0 | `"2e-4"` | – | popiii_m2c.jl:218, popiii_m2c.jl:312 | initial x_e / flat x_e (abn02fin/red recombination) | KEEP-STACK | XE0 in CKENV |
| XEFLAT | `"0"` | – | popiii_m2c.jl:218 | initial x_e / flat x_e (abn02fin/red recombination) | KEEP-STACK | XE0 in CKENV |
| LRJOIN | `"auto"` | – | popiii_m2c.jl:272 | join of the hydrostatic interior with the ABN02 n(r) (IC set-up) | KEEP-STACK |  |
| DLNJOIN | `"0.15"` | – | popiii_m2c.jl:276 | join of the hydrostatic interior with the ABN02 n(r) (IC set-up) | KEEP-STACK |  |
| ICEUL | `"0"` | – (z100: 1) | popiii_m2c.jl:1476 | start from the IC directly on the fixed log mesh (iceulrun.sh) | KEEP-STACK | the from-IC path of z100mm; make the default IC path once the Lagrangian stages go (stage 3) |
| ICXMZ | `"1"` | – | grideq_driver.jl:87, popiii_m2c.jl:1540 | ICEUL/GRIDEQ: zero the exterior mass coordinate at the start (default 1) | KEEP-STACK |  |
| ALPHAMLT | `"1.5"` | `1.5` | popiii_m2c.jl:916 | mixing-length α (stack = default 1.5) | KEEP-STACK | drop from env |
| COLDMLT | `"0"` | `1e4` | popiii_m2c.jl:917 | MLT also in hot table zones above T; MLT_NOINFALL | KEEP-STACK | stack 1e4 |
| Z0_TABLE | `""` | – | Z0Physics.jl:17 | `ext` = extended Z=0 opacity table | UNCLEAR | unset everywhere; decide which table is the one |
| OPACSENS | `"1,1"` | – | popiii_m2c.jl:902 | κ_cond, κ_rad sensitivity factors (test) | KEEP-TOOL |  |
| ENVSINKTAU | `"1e6"` | – | popiii_m2c.jl:1139 | relaxation time of the inner sink (r_in > 0) | REMOVE-SUPERSEDED | the stack has a regular centre (ENVEUL r_in = 0): no sink; G7 |

### chemistry, COUPLED (G1/G6)

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| HDINIT | `"eq"` | – | PopIIIMaterial.jl:1685 | initial HD/D⁺ state | KEEP-TOOL |  |
| COUPLED | `"0"` | `1` | TitanFortran.jl:104, popiii_m2c.jl:1265 | x_H⁺, x_H₂ as unknowns 8, 9 of the hydro Newton | MAKE-DEFAULT | stack 1 (all network zones: CPLBAND empty) |
| CPLBAND | `""` | `(empty)` | PopIIIMaterial.jl:1805 | which zones are coupled (rest split) | REMOVE-SUPERSEDED | stack empty = all network zones; the split CHEMCK path stays only for table/hot zones; superseded by global COUPLED |
| CPLCLAMP | `"1"` | `3` | TitanResidFast.jl:241 | clamp of the species correction per Newton iteration (ln x) | KEEP-STACK | stack 3 (G6 guard; keep until a globalised Newton exists) |
| CPLBADDELTA | `"100"` | – | TitanResidFast.jl:297 | decouple a zone for the rest of a failed step if \|δ ln x\| > value | KEEP-STACK | default 100 is active in the stack (G6 guard) |
| CPLXSTOL | `string(log(3.0))` | `3` | popiii_m2c.jl:1266 | Δt cap from the ln x_H⁺ change of coupled zones | REMOVE-SUPERSEDED | stack 3, but a DT_HOOK: with DTERR diagnostic only (Trbdf2.jl:198); superseded by DTERR |
| CPLIMPQ | `"0"` | – | TitanResidFast.jl:239 | three-level heats of the iterate in the coupled rows (experimental) | REMOVE-DEAD | never adopted; stage 2 |
| CPLTHETA | `"1"` | – | TitanResidFast.jl:247 | time weight of the coupled chemistry (default 1 = BE in every stage) | KEEP-STACK | default used; make the constant |
| CPLLOCAL | `"1"` | – | popiii_m2c.jl:2427 | predictor: zone species rows solved alone (default 1) | KEEP-STACK | default used |
| CPLJLAG | `"3e-3,0.1"` | – | PopIIIMaterial.jl:1885 | reuse of the coupled-zone Jacobian blocks (perf) | KEEP-STACK | default used |
| CPLTRACEX | `"0"` | – | TitanResidFast.jl:303 | hold trace species of coupled zones (no correction) | REMOVE-DEAD | S1: removed |
| CPLDBG | `"0"` | – | TitanResidFast.jl:241 | print large species corrections | REMOVE-DEAD | S1: removed |
| CPLDBGZ | `"0"` | – | PopIIIMaterial.jl:1799 | cooling terms of one coupled zone (diagnostic) | KEEP-TOOL |  |
| CPLDBGN | `"20"` | – | PopIIIMaterial.jl:1799 | cooling terms of one coupled zone (diagnostic) | KEEP-TOOL |  |
| CPLJCHECK | unset | – | popiii_m2c.jl:2430, popiii_m2c.jl:2431 | finite-difference Jacobian check of the coupled rows / of the radiation+energy rows, once | KEEP-TOOL |  |
| CPLJCHH | `"1e-6"` | – | TitanResidFast.jl:1742 | finite-difference Jacobian check of the coupled rows / of the radiation+energy rows, once | KEEP-TOOL |  |
| SZJCHECK | unset | – | popiii_m2c.jl:2434, popiii_m2c.jl:2435 | finite-difference Jacobian check of the coupled rows / of the radiation+energy rows, once | KEEP-TOOL |  |
| SZJCHH | `"1e-6"` | – | TitanResidFast.jl:1772 | finite-difference Jacobian check of the coupled rows / of the radiation+energy rows, once | KEEP-TOOL |  |
| CPLSCAN | unset | – | popiii_m2c.jl:2449, popiii_m2c.jl:2450 | thermal-bistability scan of coupled zones, then exit | REMOVE-DEAD | S1: removed |
| CPLSCANNLO | `"0"` | – | popiii_m2c.jl:2320 | thermal-bistability scan of coupled zones, then exit | REMOVE-DEAD | S1: removed |
| CHEMCK | `"0"` | `1` | ChemNet.jl:743 | reduced network (H₂, H⁺ slow; fast species algebraic) | MAKE-DEFAULT | stack 1 (implied by COUPLED) |
| CHEMCKSWEEPS | `"3"` | – | ChemNet.jl:743 | Gauss–Seidel sweeps of the CK fast species (3) | KEEP-STACK |  |
| CHEMCKHSS | `"0"` | `100` | ChemNet.jl:759 | H⁺ quasi-steady state in the split CK step | REMOVE-SUPERSEDED | stack 100, but only the split step reads it; with COUPLED everywhere only table/hot zones are split; superseded by COUPLED (verify) |
| CHEMCKADAPT | `""` | – | ChemNet.jl:750 | adaptive BE sub-steps inside the CK step | REMOVE-DEAD | S1: removed |
| CHEMSUB | `"0"` | `1.3,3` | PopIIIMaterial.jl:1276 | sub-step a split three-level step changing T/n_e by > f | REMOVE-SUPERSEDED | stack 1.3,3; split path only (CHEMCKHSS skips it); superseded by COUPLED |
| CHEMSUBMAX | `"16"` | – | PopIIIMaterial.jl:1277 | depth of CHEMSUB | REMOVE-SUPERSEDED | with CHEMSUB |
| CHEMSUBCYC | `"0"` | – | ChemNet.jl:609 | Anninos sub-cycled chemistry (step_sub) | REMOVE-DEAD | S1: removed |
| CHEMSUBMAX_SC | `"200000"` | – | ChemNet.jl:610 | Anninos sub-cycled chemistry (step_sub) | REMOVE-DEAD | S1: removed |
| CHEMLOGX | `"0"` | – | ChemNet.jl:596 | chemistry Newton in ln x | REMOVE-DEAD | S1: removed |
| CHEMROBUST | `"0"` | – | ChemNet.jl:579 | non-finite Jacobian = failed attempt; its debug print | REMOVE-DEAD | S1: removed |
| CHEMNANDBG | `"0"` | – | ChemNet.jl:704 | non-finite Jacobian = failed attempt; its debug print | REMOVE-DEAD | S1: removed |
| CHEMTRCLIP | `"0"` | – | ChemNet.jl:588 | per-species positivity clip | REMOVE-DEAD | S1: removed |
| CHEMSTAG | `"0"` | – | ChemNet.jl:716 | accept stagnated / abandon slow chemistry attempts early | REMOVE-DEAD | S1: removed |
| CHEMGIVEUP | `"0"` | – | ChemNet.jl:722 | accept stagnated / abandon slow chemistry attempts early | REMOVE-DEAD | S1: removed |
| CHEMREPLAY | `"0"` | – | PopIIIMaterial.jl:1283 | replay deep CHEMSUB zones under solver variants | REMOVE-DEAD | S1: removed |
| CHEMREPLAYTRACE | `"0"` | – | PopIIIMaterial.jl:1344 | replay deep CHEMSUB zones under solver variants | REMOVE-DEAD | S1: removed |
| CHEMSTATS | `"0"` | – | ChemNet.jl:709 | chemistry call/Jacobian counters / Newton trace of the first failed steps | KEEP-TOOL |  |
| CHEMTRACE | `"0"` | – | ChemNet.jl:726 | chemistry call/Jacobian counters / Newton trace of the first failed steps | KEEP-TOOL |  |
| CHEMZPROF | `"0"` | – | PopIIIMaterial.jl:1368 | per-zone chemistry cost profile | KEEP-TOOL |  |
| ECHEMSPLIT | `"0"` | `0` | popiii_m2c.jl:1258 | hydro advects e − Σx E0 | REMOVE-SUPERSEDED | stack 0 (run.env overrides 1); excluded by COUPLED; superseded by G1 |
| XADV2 | `get(ENV, "SL3", "0")` | – | popiii_m2c.jl:1211 | 2nd-order (van Leer) correction of the species advection | MAKE-DEFAULT | NOT set, but its default is SL3 (=1): the stack uses it; make it the constant |
| XADVU | `"0"` | – | popiii_m2c.jl:1249 | e/x mismatch after the advection (1 diag, 2 reset e) | REMOVE-DEAD | S1: removed |
| BALSPLIT | `"0"` | `0` | popiii_m2c.jl:520 | balanced splitting of chemistry heat | REMOVE-SUPERSEDED | stack 0 (= default); superseded by COUPLED |
| STOLCHEM | `"0"` | `1` | popiii_m2c.jl:1940 | re-evaluate stol after chemistry | REMOVE-SUPERSEDED | stack 1, skipped under DTERR (popiii_m2c STOLCHEM && !ERRCTL); superseded by DTERR |
| DXMAX | `"4e-3"` | – | popiii_m2c.jl:866 | HOT: composition cap / energy fraction / burn–hydro source-iteration tolerances | KEEP-STACK | defaults used by the stack (HOT=1); G1 would retire them, not yet |
| EFRAC | `"5"` | – | popiii_m2c.jl:866 | HOT: composition cap / energy fraction / burn–hydro source-iteration tolerances | KEEP-STACK | defaults used by the stack (HOT=1); G1 would retire them, not yet |
| SPLITTOL | `"3e-3"` | – | popiii_m2c.jl:869 | HOT: composition cap / energy fraction / burn–hydro source-iteration tolerances | KEEP-STACK | defaults used by the stack (HOT=1); G1 would retire them, not yet |
| SPLITIT | `"4"` | – | popiii_m2c.jl:869 | HOT: composition cap / energy fraction / burn–hydro source-iteration tolerances | KEEP-STACK | defaults used by the stack (HOT=1); G1 would retire them, not yet |
| ZSPLITTOL | `"3e-3"` | – | popiii_m2c.jl:869 | HOT: composition cap / energy fraction / burn–hydro source-iteration tolerances | KEEP-STACK | defaults used by the stack (HOT=1); G1 would retire them, not yet |
| S5NUM | `"0"` | – | popiii_m2c.jl:872 | S5 quasi-static numerics for STARFROM stars | UNCLEAR | with STARFROM |
| S5TOL | `"0.03"` | – | popiii_m2c.jl:878 | S5 quasi-static numerics for STARFROM stars | UNCLEAR | with STARFROM |

### material (G2), table ↔ network

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| SCVHXE | `"0"` | `1` | PopIIIEqEOS.jl:54 | G2: SCVH95 + x_e table (EQEOS_CACHE) | KEEP-STACK | stack 1; MAKE-DEFAULT once the cache is in POPIII_LIT |
| HM88 | `"0"` | `1` | ChemNet.jl:494 | G2: Hummer–Mihalas occupation probabilities in the network EOS | KEEP-STACK | stack 1; MAKE-DEFAULT |
| G2OPAC | `"0"` | `2` | G2Opacity.jl:230, G2Opacity.jl:36 | G2: network-composition continuum opacity (2 = fast form) | KEEP-STACK | stack 2; MAKE-DEFAULT |
| G2NFAST | `"400"` | – | G2Opacity.jl:134 | G2OPAC=2: frequency points of the fast means | KEEP-STACK |  |
| NLOHBFDEEP | `"0"` | `1` | PopIIIMaterial.jl:417 | G2: H bound–free in κ_P of the deep network zones | KEEP-STACK | stack 1; MAKE-DEFAULT |
| H2PF | `0` | – | ChemNet.jl:76 (g2-domain) | g2-domain: H₂ partition function from the RACPPK level sum (branch, unmerged) | UNCLEAR | not in the stack yet; adopt together with G2BLEND or drop |
| G2BLEND | `0` | – | PopIIIMaterial.jl:435 (g2-domain) | g2-domain: smooth network/table blend by the non-ideality s (branch, unmerged) | UNCLEAR | if adopted it overrides PSWITCH/PSWNH/PSWTAU/PSWTEX/PSWSTAY: those become REMOVE-SUPERSEDED with it |
| PSWITCH | `"1e4"` | `1e30` | popiii_m2c.jl:366 | table where p_gas ≥ PSWITCH | REMOVE-SUPERSEDED | stack 1e30 = rule off; superseded by PSWNH (G2: table only in the dense interior); in CKENV: translate (stage 2) |
| PSWTAU | `"0"` | `0` | popiii_m2c.jl:1025 | table where T ≥ PSWTEX and τ ≥ PSWTAU | REMOVE-SUPERSEDED | stack 0 = off (run.env overrides phaseB.env 3); superseded by G2 |
| PSWSTAY | `"1"` | `0` | popiii_m2c.jl:538 | stay rule of table zones (T ≥ PSWTEX) | REMOVE-SUPERSEDED | stack 0 = off; superseded by G2 |
| PSWTEX | `"0"` | `5000` | popiii_m2c.jl:1025 | table zone stays on table until T < PSWTEX | REMOVE-SUPERSEDED | stack 5000 but only read through PSWTAU/PSWSTAY (both off): no effect; superseded by G2 |
| PSWHYS | `"1"` | `0.3` | popiii_m2c.jl:1023 | hysteresis factor of the table rules | KEEP-STACK | stack 0.3: still used by the PSWNH rule (n_H ≥ PSWHYS·PSWNH to stay) |
| PSWTMIN | `"0"` | `4500` | popiii_m2c.jl:1024 | table also needs T ≥ PSWTMIN (4500) | KEEP-STACK | still used by the PSWNH rule; part of the dense-interior table rule |
| PSWNH | `"0"` | `1e19` | PopIIIMaterial.jl:822 | table where n_H ≥ PSWNH (non-ideal interior) | KEEP-STACK | stack 1e19: the only active table rule (G2); replaced by G2BLEND if that is adopted |
| PSWRAMP | `"0"` | `10` | PopIIIMaterial.jl:815 | gradual table→network hand-over, one zone per n steps | UNCLEAR | stack 10; only acts at a restart that moves the table edge (RAMPK0); decide by a restart test with it unset |
| PSWAGREE | `"0"` | – | popiii_m2c.jl:531 | table only where table and three-level steady state agree | REMOVE-DEAD | S1: removed (with PSWTAUDEEP, xss_tl3, T_ss_at_u) |
| PSWRAD | `"0"` | – | PopIIIMaterial.jl:867 | enter the table only if T_gas ≈ T_rad | REMOVE-DEAD | S1: removed |
| CHEMFBDENSE | `"1"` | – | PopIIIMaterial.jl:902 | failed network step at p ≥ pswitch → table fallback | REMOVE-SUPERSEDED | with PSWITCH=1e30 never reached (p ≥ 1e30); stage 2 with PSWITCH |
| CHEMFBNET | `"0"` | `16` | PopIIIMaterial.jl:1045 | failed CK step: retry on the full network up to n times | KEEP-STACK | stack 16 |
| CHEMTMIN | `"0"` | – | PopIIIMaterial.jl:1364 | chemistry step below T floor counts as failed → table | REMOVE-SUPERSEDED | off; G2/G6 patch never set by the stack; stage 2 |
| NLOPAC | `"0"` | `1` | popiii_m2c.jl:530 | network zones: non-LTE opacity (table × ratio) | MAKE-DEFAULT | stack 1; base of NLOIT |
| NLOIT | `"0"` | `1` | PopIIIMaterial.jl:485 | network-zone opacities from the Newton iterate (G2, lags removed) | MAKE-DEFAULT | stack 1 |
| NLOITW | `"5000,6000"` | – | PopIIIMaterial.jl:486 | NLOIT: base blend window (table below 5000 K, continuum code above 6000 K) / slope reuse | KEEP-STACK | defaults used |
| NLOITJ | `"0,0"` | – | PopIIIMaterial.jl:487 | NLOIT: base blend window (table below 5000 K, continuum code above 6000 K) / slope reuse | KEEP-STACK | defaults used |
| NLOIMPT | `"0"` | `1` | PopIIIMaterial.jl:446 | lagged-NLO patches: T extrapolation of the ratio / χ floor σ_T n_e/ρ | REMOVE-SUPERSEDED | stack 1, but inactive in network zones with NLOIT (PopIIIMaterial.jl NLOIT docstring); superseded by NLOIT; verify with a bit test (stage 2) |
| NLOFIX | `"0"` | `1` | PopIIIMaterial.jl:740 | lagged-NLO patches: T extrapolation of the ratio / χ floor σ_T n_e/ρ | REMOVE-SUPERSEDED | stack 1, but inactive in network zones with NLOIT (PopIIIMaterial.jl NLOIT docstring); superseded by NLOIT; verify with a bit test (stage 2) |
| NLODAMP | `"1"` | – | PopIIIMaterial.jl:321 | under-relaxation / x_H⁺ slope / T cut of the lagged NLO factors | REMOVE-SUPERSEDED | unset (defaults), inactive with NLOIT; superseded by NLOIT |
| CPLNLOX | `"1"` | – | PopIIIMaterial.jl:1798 | under-relaxation / x_H⁺ slope / T cut of the lagged NLO factors | REMOVE-SUPERSEDED | unset (defaults), inactive with NLOIT; superseded by NLOIT |
| NLOTMAX | `"Inf"` | – | PopIIIMaterial.jl:747 | under-relaxation / x_H⁺ slope / T cut of the lagged NLO factors | REMOVE-SUPERSEDED | unset (defaults), inactive with NLOIT; superseded by NLOIT |
| NLOPACJ | `"0"` | `1` | popiii_m2c.jl:533 | κ_E weighted with B_ν(T_c) of the star light | MAKE-DEFAULT | stack 1; also used by NLOIT (nlo_ratio, PopIIIMaterial.jl:576) |
| CHICAP | `"0"` | – | PopIIIMaterial.jl:464 | cap χ (diagnostic switch, "not physics") | KEEP-TOOL |  |
| LOCALJ | `"0"` | `1` | PopIIIMaterial.jl:1028 | network photo-processes see W B(T_c) instead of the CMB-only field | MAKE-DEFAULT | stack 1 |
| LOCALJWMIN | `"0"` | `1e-3` | PopIIIMaterial.jl:1033 | LOCALJ: minimum dilution | KEEP-STACK | stack 1e-3 |
| LOCALJTAB | `"1"` | – | PopIIIMaterial.jl:1030 | LOCALJ photo-rate tables (default on) / table size | KEEP-TOOL | performance |
| LJTABN | `"1000"` | – | ChemNet.jl:352 | LOCALJ photo-rate tables (default on) / table size | KEEP-TOOL | performance |
| EPTAB | `"0"` | `1` | H2Lines.jl:629 | tabulated EP line-cooling kernels | MAKE-DEFAULT | stack 1 (and repo phaseB.env) |
| LINEFAST | `"0"` | `1` | H2Lines.jl:431 | fast resonance window of the line transfer / its width (38, 16 in the repo phaseB.env) | MAKE-DEFAULT | stack LINEFAST=1, LINEXC2 default 38 |
| LINEXC2 | `"38"` | – | H2Lines.jl:432 | fast resonance window of the line transfer / its width (38, 16 in the repo phaseB.env) | MAKE-DEFAULT | stack LINEFAST=1, LINEXC2 default 38 |
| NLOITDIAG | `"0"` | – | PopIIIMaterial.jl:675 | NLOIT: per-zone opacity diagnostics | KEEP-TOOL |  |

### star light J_ν closure (G4)

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| SLBLEND | `"0"` | `1` | popiii_m2c.jl:532 | star-light colour temperature blended in τ / Eddington–Barbier / dilution from max(c(E_r−E_cmb), F_r) | KEEP-STACK | stack 1; patches of the J_ν closure until a per-zone spectrum exists (G4): keep for now |
| SLEB | `"0"` | `1` | popiii_m2c.jl:535 | star-light colour temperature blended in τ / Eddington–Barbier / dilution from max(c(E_r−E_cmb), F_r) | KEEP-STACK | stack 1; patches of the J_ν closure until a per-zone spectrum exists (G4): keep for now |
| SLWFLUX | `"0"` | `1` | PopIIIMaterial.jl:749 | star-light colour temperature blended in τ / Eddington–Barbier / dilution from max(c(E_r−E_cmb), F_r) | KEEP-STACK | stack 1; patches of the J_ν closure until a per-zone spectrum exists (G4): keep for now |
| SLIMP | `"0"` | `2` | PopIIIMaterial.jl:1117 | star light at the Newton iterate (2) | KEEP-STACK | stack 2; MAKE-DEFAULT |
| SLFREEZE | `"0"` | – | PopIIIMaterial.jl:1147 | freeze the star-light parameters (test) | REMOVE-DEAD | never set; stage 2 |
| SUBZONE | `"0"` | `1,2` | popiii_m2c.jl:1972, popiii_m2c.jl:1973 | interface optical depth across an unresolved front (two-layer sub-zone form) | KEEP-STACK | stack 1,2; MAKE-DEFAULT |

### H₂ lines

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| LINEIMP | `"0"` | `1` | popiii_m2c.jl:580 | line cooling of each zone at the iterate (lag removed) | MAKE-DEFAULT | stack 1 |
| LIFMAX | `"10"` | – | popiii_m2c.jl:580 | LINEIMP: bound of the ALI/EP factor (10) | KEEP-STACK |  |
| IMPCOOL | `"1"` | – | popiii_m2c.jl:342, popiii_m2c.jl:344 | line-cooling slope from a second transfer call | REMOVE-SUPERSEDED | unset (default IMPCOOL=1), not used with LINEIMP (popiii_m2c LINEIMP docstring); superseded by LINEIMP. IMPCOOL in CKENV: translate |
| IMPALWAYS | `"0"` | – | popiii_m2c.jl:344 | line-cooling slope from a second transfer call | REMOVE-SUPERSEDED | unset (default IMPCOOL=1), not used with LINEIMP (popiii_m2c LINEIMP docstring); superseded by LINEIMP. IMPCOOL in CKENV: translate |
| DTCOOL | `"1.0"` | – | popiii_m2c.jl:343 | Δt ≤ DTCOOL t_cool | REMOVE-SUPERSEDED | no effect with LINEIMP and DTERR; superseded; in CKENV |
| TCOOLTH | `"0"` | `1` | popiii_m2c.jl:349 | t_cool from 3p/2ρ (only the cap uses it) | REMOVE-SUPERSEDED | stack 1, no effect with LINEIMP/DTERR |
| LINEFRAC | `"0"` | `0` | popiii_m2c.jl:1978 | reuse Λ until LINEFRAC × t_cool | REMOVE-SUPERSEDED | stack 0 (= default) |
| NALI | `"20"` | – | popiii_m2c.jl:551 | hybrid NLTE: ALI refresh interval and its adaptation | KEEP-STACK | stack ALIWMIN 1e-3, ALIGROW 0.1, ALISHRINK 0.5 (≠ defaults: make those the defaults) |
| ALIWMIN | `"0"` | `1e-3` | popiii_m2c.jl:551 | hybrid NLTE: ALI refresh interval and its adaptation | KEEP-STACK | stack ALIWMIN 1e-3, ALIGROW 0.1, ALISHRINK 0.5 (≠ defaults: make those the defaults) |
| ALIGROW | `"0.01"` | `0.1` | popiii_m2c.jl:551 | hybrid NLTE: ALI refresh interval and its adaptation | KEEP-STACK | stack ALIWMIN 1e-3, ALIGROW 0.1, ALISHRINK 0.5 (≠ defaults: make those the defaults) |
| ALISHRINK | `"0.05"` | `0.5` | popiii_m2c.jl:551 | hybrid NLTE: ALI refresh interval and its adaptation | KEEP-STACK | stack ALIWMIN 1e-3, ALIGROW 0.1, ALISHRINK 0.5 (≠ defaults: make those the defaults) |
| ALICAP | `"0"` | – | popiii_m2c.jl:551 | hybrid NLTE: ALI refresh interval and its adaptation | KEEP-STACK | stack ALIWMIN 1e-3, ALIGROW 0.1, ALISHRINK 0.5 (≠ defaults: make those the defaults) |
| LLDIAG | `"0"` | – | popiii_m2c.jl:715 | line-luminosity closure / per-zone consistency / LINEINC diagnostics | KEEP-TOOL |  |
| LINECONS | `"0"` | – | popiii_m2c.jl:723 | line-luminosity closure / per-zone consistency / LINEINC diagnostics | KEEP-TOOL |  |
| LINEINCDIAG | `"0"` | – | popiii_m2c.jl:566 | line-luminosity closure / per-zone consistency / LINEINC diagnostics | KEEP-TOOL |  |
| LYART | `"0"` | – | LyaTransfer.jl:26 | Lyα transfer for the three-level atom (LyaTransfer.jl) instead of Sobolev escape | UNCLEAR | off everywhere; a physics method, not a patch: keep if Lyα transfer is planned, otherwise delete LyaTransfer.jl (~200 lines) |
| LYATAU | `"20"` | – | LyaTransfer.jl:27 | Lyα transfer for the three-level atom (LyaTransfer.jl) instead of Sobolev escape | UNCLEAR | off everywhere; a physics method, not a patch: keep if Lyα transfer is planned, otherwise delete LyaTransfer.jl (~200 lines) |
| LYAGAMMA | `"1"` | – | LyaTransfer.jl:28 | Lyα transfer for the three-level atom (LyaTransfer.jl) instead of Sobolev escape | UNCLEAR | off everywhere; a physics method, not a patch: keep if Lyα transfer is planned, otherwise delete LyaTransfer.jl (~200 lines) |
| LYAH | `"0.2"` | – | LyaTransfer.jl:36 | Lyα transfer for the three-level atom (LyaTransfer.jl) instead of Sobolev escape | UNCLEAR | off everywhere; a physics method, not a patch: keep if Lyα transfer is planned, otherwise delete LyaTransfer.jl (~200 lines) |
| LYAVMAX | `"2e5"` | – | LyaTransfer.jl:36 | Lyα transfer for the three-level atom (LyaTransfer.jl) instead of Sobolev escape | UNCLEAR | off everywhere; a physics method, not a patch: keep if Lyα transfer is planned, otherwise delete LyaTransfer.jl (~200 lines) |
| LYANSUB | `"24"` | – | LyaTransfer.jl:98 | Lyα transfer for the three-level atom (LyaTransfer.jl) instead of Sobolev escape | UNCLEAR | off everywhere; a physics method, not a patch: keep if Lyα transfer is planned, otherwise delete LyaTransfer.jl (~200 lines) |
| LYASUBTOP | `"30"` | – | LyaTransfer.jl:98 | Lyα transfer for the three-level atom (LyaTransfer.jl) instead of Sobolev escape | UNCLEAR | off everywhere; a physics method, not a patch: keep if Lyα transfer is planned, otherwise delete LyaTransfer.jl (~200 lines) |
| LYANCORE | `"6"` | – | LyaTransfer.jl:69 | Lyα transfer for the three-level atom (LyaTransfer.jl) instead of Sobolev escape | UNCLEAR | off everywhere; a physics method, not a patch: keep if Lyα transfer is planned, otherwise delete LyaTransfer.jl (~200 lines) |
| LYAMAXIT | `"40"` | – | LyaTransfer.jl:168 | Lyα transfer for the three-level atom (LyaTransfer.jl) instead of Sobolev escape | UNCLEAR | off everywhere; a physics method, not a patch: keep if Lyα transfer is planned, otherwise delete LyaTransfer.jl (~200 lines) |
| LYATOL | `"1e-3"` | – | LyaTransfer.jl:168 | Lyα transfer for the three-level atom (LyaTransfer.jl) instead of Sobolev escape | UNCLEAR | off everywhere; a physics method, not a patch: keep if Lyα transfer is planned, otherwise delete LyaTransfer.jl (~200 lines) |
| LYANG | `"4"` | – | LyaTransfer.jl:169 | Lyα transfer for the three-level atom (LyaTransfer.jl) instead of Sobolev escape | UNCLEAR | off everywhere; a physics method, not a patch: keep if Lyα transfer is planned, otherwise delete LyaTransfer.jl (~200 lines) |

### radiation (G4)

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| VEFIMP | `"0"` | `1` | popiii_m2c.jl:2065, popiii_m2c.jl:2066 | Eddington factors at every Newton iterate, Auer form, CMB as background | MAKE-DEFAULT | stack 1 |
| VEFBGCHK | `"0"` | – | popiii_m2c.jl:2081, popiii_m2c.jl:2082 | check of the VEF background split against the full formal solution | KEEP-TOOL |  |
| EDDGMU0 | `"0"` | `1` | popiii_m2c.jl:2060 | outer g integrated to μ = 0 | MAKE-DEFAULT | stack 1; a correctness fix |
| EDDSCAT | `"0"` | `1` | popiii_m2c.jl:2053 | formal solution: χ above κ_P treated as scattering | KEEP-STACK | stack 1; G4/G2 patch — keep until the opacity split is consistent |
| EDDREAL | `"0"` | – | popiii_m2c.jl:1945 | enforce f ≥ (F/cE)² above τ | REMOVE-DEAD | S1: removed (superseded by VEFIMP) |
| EDDTAU | `"0.3"` | – | popiii_m2c.jl:1946 | enforce f ≥ (F/cE)² above τ | REMOVE-DEAD | S1: removed (superseded by VEFIMP) |
| EDDRELAX | `"1"` | – | popiii_m2c.jl:1946 | enforce f ≥ (F/cE)² above τ | REMOVE-DEAD | S1: removed (superseded by VEFIMP) |
| EDDCHI | `"0"` | – | popiii_m2c.jl:2097 | χ ≤ κ_P + σ_e in the formal solution | REMOVE-DEAD | S1: removed |
| EDDDBG | unset | – | popiii_m2c.jl:2108, popiii_m2c.jl:2109 | formal-solution J, f vs E, F per zone | REMOVE-DEAD | S1: removed |
| FLUXTIE | unset | `0.1,1` | popiii_m2c.jl:2054, popiii_m2c.jl:2055 | thin zones: F = h c E (h of the iterate with VEFIMP) | KEEP-STACK | stack 0.1,1 |
| FLUXTIEDIAG | unset | – | TitanCore.jl:821, TitanCore.jl:822 … | print h, F/cE, f per zone | REMOVE-DEAD | S1: removed |
| ERTAUCONV | `"1.0"` | `0` | popiii_m2c.jl:911 | causality floor of E_r in thin zones | REMOVE-SUPERSEDED | stack 0 = off (phaseB.env); superseded by VEFIMP; delete the floor code |
| ERTAUMIN | `"1.0"` | `3000` | popiii_m2c.jl:912 | E_r of the thin envelope / thick zones in the Δt rule | REMOVE-SUPERSEDED | Δt-rule only: diagnostic under DTERR; superseded by DTERR |
| ERTHICK | unset | – | popiii_m2c.jl:1962 | E_r of the thin envelope / thick zones in the Δt rule | REMOVE-SUPERSEDED | Δt-rule only: diagnostic under DTERR; superseded by DTERR |
| DTRADMAX | `"1500"` | `0` | popiii_m2c.jl:919 | Δt ≤ DTRADMAX (DT_HOOK) | REMOVE-SUPERSEDED | stack 0 = off; superseded by DTERR |
| ETASWTH | unset | `1e-3` | popiii_m2c.jl:1964 | internal-energy row also in radiation-dominated / flux-tied zones | KEEP-STACK | stack 1e-3 / 0.5 (G5 patches tied to FLUXTIE; keep while FLUXTIE is used) |
| ETASWTIE | unset | `0.5` | popiii_m2c.jl:1965 | internal-energy row also in radiation-dominated / flux-tied zones | KEEP-STACK | stack 1e-3 / 0.5 (G5 patches tied to FLUXTIE; keep while FLUXTIE is used) |
| ERDIAG | unset | – | popiii_m2c.jl:2598 | radiation field / profile print at given steps | KEEP-TOOL |  |
| ERDIAGZ | unset | – | popiii_m2c.jl:2603 | radiation field / profile print at given steps | KEEP-TOOL |  |
| PROFAT | unset | – | popiii_m2c.jl:2608 | radiation field / profile print at given steps | KEEP-TOOL |  |
| EXCAP | unset | `1e8` | popiii_m2c.jl:1968 | scale gas–radiation exchange so c κ_E ρ Δt ≤ EXCAP | KEEP-STACK | stack 1e8 (round-off guard); G5 reformulation would retire it |
| OLDOPAC | `"0"` | – | popiii_m2c.jl:1969 | COUPLED: re-evaluate start-of-step opacities with the Newton material (θ-scheme old level consistent) | UNCLEAR | a consistency fix never adopted; test with the stack: MAKE-DEFAULT or remove |

### energy (G5)

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| ETOTAL | `"0"` | `1` | popiii_m2c.jl:1967 | conservative total-energy row (gas row where e/(e+KE) < ETASW) | MAKE-DEFAULT | stack 1 |
| ETASW | `"0.1"` | `0.6` | popiii_m2c.jl:1967 | dual-energy threshold | KEEP-STACK | stack 0.6 |
| ETOTGRAVGRID | `"0"` | – | popiii_m2c.jl:1966 | gravity work on moving-face mass in the total-energy row (Lagrangian zones) | REMOVE-SUPERSEDED | unset; Lagrangian stages only (stage 3) |
| TRESET | unset | – | popiii_m2c.jl:1815 | one-time set of zone T at a restart (repair) | REMOVE-SUPERSEDED | never set by the stack; G7 one-run (stage 3) |

### hydro, viscosity

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| AVZONE | `"0"` | `2` | popiii_m2c.jl:1454, popiii_m2c.jl:1546 … | viscous length = AVZONE × zone width | MAKE-DEFAULT | stack 2; AVQL1 then unused |
| AVQL1 | `"0.15"` | – | popiii_m2c.jl:374 | viscous length ql1 r̄ | REMOVE-SUPERSEDED | unused when AVZONE > 0 (TitanResidFast: ql = AVZONE·Δr); in CKENV |
| AVCQ2 | `"1.0"` | – | popiii_m2c.jl:374 | quadratic viscosity coefficient (1.0) | KEEP-STACK | default used; in CKENV |
| AVTENSION | `"1"` | `0` | popiii_m2c.jl:379, popiii_m2c.jl:380 | 0 = no viscous tension in shrinking zones | KEEP-STACK | stack 0 ≠ default 1: make 0 the default |
| AVCAP | `"0"` | – | popiii_m2c.jl:383 | cap viscous length at AVCAP × r̄ | REMOVE-DEAD | never set; stage 2 (residual) |
| SHOCKRAD | `"0"` | – | popiii_m2c.jl:378 | sub-grid radiative shock (viscous heat → radiation) | REMOVE-DEAD | never set by a production run; stage 2 (residual) |
| SHOCKDU | unset | – | popiii_m2c.jl:384 | sub-grid radiative shock (viscous heat → radiation) | REMOVE-DEAD | never set by a production run; stage 2 (residual) |
| ADVREL | `"0"` | – | TitanResidFast.jl:134 | unconverged relative mass flux carries momentum and energy consistently | UNCLEAR | phase A only; on the moving mesh (GRIDEQ) it may matter: test |
| UPFREEZE | `"0"` | – | TitanResidFast.jl:601 | freeze upwind direction from Newton iteration n (for GRIDWIN) | REMOVE-DEAD | S1: removed |

### mesh (G7)

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| ENVEUL | unset | `0,0.0793995,0.0793995,1e…` (z100: 0,0.232584376,0.23…) | popiii_m2c.jl:1133 | fixed log mesh: r_in, Δln r, (regular centre: first interface, Δln r below r_cb) | KEEP-STACK | stack ENVEUL=0,0.0794,0.0794,1e30 ENVCORE=0.1,0.0794,1e30; the remap at a restart (lgrid 3 → 2) is stage 3 |
| ENVCORE | `"0.1,0.08,10"` | `0.1,0.0793995,1e30` (z100: 0.1,0.232584376,1e…) | popiii_m2c.jl:1136 | fixed log mesh: r_in, Δln r, (regular centre: first interface, Δln r below r_cb) | KEEP-STACK | stack ENVEUL=0,0.0794,0.0794,1e30 ENVCORE=0.1,0.0794,1e30; the remap at a restart (lgrid 3 → 2) is stage 3 |
| ENVLIN | `"1"` | `1` | popiii_m2c.jl:1137 | zones of equal Δr between r = 0 and the first interface | KEEP-STACK | stack 1 (= default; run.env sets 10 then 1) |
| ENVCONS | `"0"` | `1` | popiii_m2c.jl:1138 | mass-conserving ENVEUL remap / IC fill | MAKE-DEFAULT | stack 1 |
| ENVMS0 | unset | – | popiii_m2c.jl:1570 | reference shell mass of the sink at a restart | REMOVE-SUPERSEDED | regular centre: unused; G7 (stage 3) |
| ENVDT0 | `"100"` | – | popiii_m2c.jl:1457, popiii_m2c.jl:1550 … | first Δt after a remap/split/merge or ICEUL start | REMOVE-SUPERSEDED | stage 3 with the remaps (ICEUL keeps its own 1e-3 t_ff default) |
| ENVSPLIT | unset | – | popiii_m2c.jl:1785 | one-time split/merge of fixed-mesh zones at a restart | REMOVE-SUPERSEDED | superseded by GRIDEQ (G7); stage 3 |
| ENVMERGE | unset | – | popiii_m2c.jl:1786 | one-time split/merge of fixed-mesh zones at a restart | REMOVE-SUPERSEDED | superseded by GRIDEQ (G7); stage 3 |
| DTENV | unset | – | popiii_m2c.jl:1773 | Δt ≤ C × dynamical time over r1…r2 (DT_HOOK) | REMOVE-SUPERSEDED | repo phaseB.env sets it, the stack does not; DT_HOOK = diagnostic under DTERR |
| GRIDEQ | unset | `1` (g2: –; z100: 0) | grideq_driver.jl:35 | TITAN implicit grid equation on the whole mesh (moving mesh) | KEEP-STACK | z300mm 1, z300g2/z100mm off |
| GEQC | `[100.0, 1.5, 30.0, 0…` | `100,1.5,30,0.5,0.5,0.1,1…` (g2: –) | grideq_driver.jl:38 | GRIDEQ: constants (gmax, α, κ, fv, fmove, τ, ctol …) / monitor weights / start-up ramp | KEEP-STACK |  |
| GEQW | `[1.0, 1.0, 1.0, 1.0,…` | `1,1,1,1,1,0,0,3` (g2: –) | grideq_driver.jl:37 | GRIDEQ: constants (gmax, α, κ, fv, fmove, τ, ctol …) / monitor weights / start-up ramp | KEEP-STACK |  |
| GEQRAMP | unset | `1e-3,300,10,20` (g2: –; z100: 0) | grideq_driver.jl:40 | GRIDEQ: constants (gmax, α, κ, fv, fmove, τ, ctol …) / monitor weights / start-up ramp | KEEP-STACK |  |
| GEQX | `[0.0, 1e-3]` | – | grideq_driver.jl:39 | GRIDEQ: lagged feature-tracking term (default off) | KEEP-STACK | default used |
| GEQN | `string(nsh)` | – | popiii_m2c.jl:1576 | GRIDEQ at a remap: number of zones | UNCLEAR | only with a remap restart (stage 3) |
| ICMLOG | unset | – | popiii_m2c.jl (iceulbuild!) | IC start mesh logarithmic in enclosed mass from m1 [Msun] (N = GEQN or nsh zones); replaces the rule's starting mesh; first interface as the default start (GRIDEQ keeps it fixed), interfaces 2..N log in mass from m1 to the cloud mass | KEEP-STACK | start from the IC only |
| GEQREMAP | `"0"` | – | grideq_driver.jl:55 | GRIDEQ remap of the state ("kept but not recommended") | REMOVE-DEAD | never set; stage 2 |
| GEQDBG | `"0"` | – | TitanResidFast.jl:965, TitanResidFast.jl:966 | GRIDEQ debug print of the monitor | KEEP-TOOL |  |
| GRIDMAP | unset | – | popiii_m2c.jl:1154 | planned moving mesh (GridMap.jl) | REMOVE-DEAD | S1: removed (superseded by GRIDEQ) |
| GMJUMP | `"0"` | – | popiii_m2c.jl:1162 | planned moving mesh (GridMap.jl) | REMOVE-DEAD | S1: removed (superseded by GRIDEQ) |
| GMTAU | unset | – | popiii_m2c.jl:1164 | planned moving mesh (GridMap.jl) | REMOVE-DEAD | S1: removed (superseded by GRIDEQ) |
| GMWIN | unset | – | popiii_m2c.jl:1177 | planned moving mesh (GridMap.jl) | REMOVE-DEAD | S1: removed (superseded by GRIDEQ) |
| GMSMOOTH | `"1"` | – | popiii_m2c.jl:1189 | planned moving mesh (GridMap.jl) | REMOVE-DEAD | S1: removed (superseded by GRIDEQ) |
| GRIDWIN | unset | – | popiii_m2c.jl:1825 | Dorfi–Drury window on the fixed mesh | REMOVE-DEAD | S1: removed (superseded by GRIDEQ) |
| REZONE | `"0"` | `0` | popiii_m2c.jl:1019 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZEVERY | `"20"` | `20` | popiii_m2c.jl:1019 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZJSPLIT | `"0.3"` | `Inf` | popiii_m2c.jl:1020 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZJMERGE | `"0.05"` | – | popiii_m2c.jl:1020 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZDMMIN | `"1e-7"` | – | popiii_m2c.jl:1021 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZQMAX | `"2.5"` | – | popiii_m2c.jl:1021 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZNMAX | `"260"` | – | popiii_m2c.jl:1022 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZFRAC | `"0.01"` | – | popiii_m2c.jl:1022 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZKEEPC | `"5"` | – | popiii_m2c.jl:1022 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZDR | `"0"` | – | popiii_m2c.jl:1032 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZPSPLIT | `"1"` | – | popiii_m2c.jl:1065 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZSLOPE | `""` | – | popiii_m2c.jl:1033 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZWIN | unset | – | popiii_m2c.jl:1026 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZNULL | `"0"` | – | popiii_m2c.jl:1026 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| REZDBG | `"0"` | – | popiii_m2c.jl:1052 | conservative mass rezoning of the Lagrangian mesh (Rezone.jl) | REMOVE-SUPERSEDED | REZONE=0 in the stack (REZEVERY/REZJSPLIT set but unread); Lagrangian stages only; stage 3 |
| RADFROM | unset | – | popiii_m2c.jl:1561, popiii_m2c.jl:1563 … | stage A2: radiation rebuild from a hydro-only checkpoint / static pre-refinement | REMOVE-SUPERSEDED | superseded by ICEUL (G7); stage 3 |
| REFINE | `:-1e-3,0.1,0.03` | – | fullrun.sh:32, popiii_m2c.jl:1564 | stage A2: radiation rebuild from a hydro-only checkpoint / static pre-refinement | REMOVE-SUPERSEDED | superseded by ICEUL (G7); stage 3 |
| MSTOP | `"0"` | – | popiii_m2c.jl:1931 | stop at M★ (A3 → B hand-over) | REMOVE-SUPERSEDED | stage chain; stage 3 |
| TSTOP | `"6000"` | `2e8` | popiii_m2c.jl:753 | stop at central T / n | KEEP-STACK | stack 2e8 / 1e27 (never reached); keep as run limits, make those the defaults |
| NSTOP | `"1e18"` | `1e27` | popiii_m2c.jl:753 | stop at central T / n | KEEP-STACK | stack 2e8 / 1e27 (never reached); keep as run limits, make those the defaults |

### time step (G3)

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| DTERR | unset | `0.1` | Trbdf2.jl:358, Trbdf2.jl:359 | TR-BDF2 embedded error estimate is the only Δt control | MAKE-DEFAULT | stack 0.1 |
| DTERRCTL | `"0.9,0.2,5"` | – | Trbdf2.jl:359 | DTERR controller / PI / Newton-failure factor / second filter / norm | KEEP-STACK | defaults used |
| DTERRPI | `"1"` | – | Trbdf2.jl:361 | DTERR controller / PI / Newton-failure factor / second filter / norm | KEEP-STACK | defaults used |
| DTERRFAIL | `"0.25"` | – | Trbdf2.jl:361 | DTERR controller / PI / Newton-failure factor / second filter / norm | KEEP-STACK | defaults used |
| DTERRFILT | `"2"` | – | Trbdf2.jl:362 | DTERR controller / PI / Newton-failure factor / second filter / norm | KEEP-STACK | defaults used |
| DTERRNORM | `"rms"` | – | Trbdf2.jl:362 | DTERR controller / PI / Newton-failure factor / second filter / norm | KEEP-STACK | defaults used |
| DTERRBE | `"0"` | – | Trbdf2.jl:362 | DTERR test options (BE retry, Δt bound) | KEEP-TOOL |  |
| DTERRMAX | `"Inf"` | – | Trbdf2.jl:362 | DTERR test options (BE retry, Δt bound) | KEEP-TOOL |  |
| DTERRDIAG | `"0"` | – | Trbdf2.jl:368 | DTERR stage-derivative diagnostics / trace interval | KEEP-TOOL |  |
| DTERRTRACE | `"20"` | – | Trbdf2.jl:362 | DTERR stage-derivative diagnostics / trace interval | KEEP-TOOL |  |
| INTEGRATOR | `"theta"` | – | popiii_m2c.jl:32 | TR-BDF2 (default) or TITAN θ scheme | MAKE-DEFAULT | default trbdf2 used; the θ branch is TITAN parity only; in CKENV |
| THETA | `"0.55"` | – | popiii_m2c.jl:399 | θ of the TITAN scheme (0.55) | UNCLEAR | TR-BDF2 sets θ per stage; check that thet enters nothing else (predict!, totnrg) before removing; in CKENV |
| CDYN | `"0"` | `3e5` | popiii_m2c.jl:400 | Δt ≤ CDYN × dynamical time | REMOVE-SUPERSEDED | stack 3e5, read only by timstp! = diagnostic under DTERR; superseded by DTERR |
| STOL | `"0.05"` | `0.3` | popiii_m2c.jl:1936 | TITAN max relative change per step | REMOVE-SUPERSEDED | stack 0.3; timstp! only = diagnostic under DTERR |
| FAILTFAC | unset | `0.5` | popiii_m2c.jl:1942 | Δt factor after a failed attempt | REMOVE-SUPERSEDED | stack 0.5; DTERR uses DTERRFAIL; superseded |
| NITER | unset | `25` | popiii_m2c.jl:1935 | Newton iterations per attempt | KEEP-STACK | stack 25 (used by newton!) |
| NTRY | unset | `60` | popiii_m2c.jl:1941 | attempts before giving up | KEEP-STACK | stack 60 (DTERR rejections count against it) |
| DTMAXTEST | unset | – | popiii_m2c.jl:2412, popiii_m2c.jl:2413 | test: Δt ≤ value (Δt → 0 limit) | KEEP-TOOL |  |
| STEPKX | `"0"` | – | popiii_m2c.jl:2594 | which equation sets the step / fastest infall / energy terms of zones z1..z2 / star light of zones z1..z2 | KEEP-TOOL |  |
| UDIAG | `"0"` | – | popiii_m2c.jl:2629 | which equation sets the step / fastest infall / energy terms of zones z1..z2 / star light of zones z1..z2 | KEEP-TOOL |  |
| ZDIAG | unset | – | popiii_m2c.jl:2463, popiii_m2c.jl:2464 | which equation sets the step / fastest infall / energy terms of zones z1..z2 / star light of zones z1..z2 | KEEP-TOOL |  |
| SLZ | unset | – | popiii_m2c.jl:2485 | which equation sets the step / fastest infall / energy terms of zones z1..z2 / star light of zones z1..z2 | KEEP-TOOL |  |
| FORCECHECK | `"0"` | – | popiii_m2c.jl:390 | initial force balance per interface | KEEP-TOOL |  |
| TABAUDIT | `"1"` | – | popiii_m2c.jl:2240, popiii_m2c.jl:2260 … | table vs network LTE audit of the table zones, then exit | KEEP-TOOL | G2 audit tool; drop with the table (G2BLEND) |
| TABAUDNH | `"17"` | – | popiii_m2c.jl:2260 | table vs network LTE audit of the table zones, then exit | KEEP-TOOL | G2 audit tool; drop with the table (G2BLEND) |
| JACCHK | unset | – | popiii_m2c.jl:1992, popiii_m2c.jl:1993 | FD Jacobian check of absolute zone z at step k, iteration iter | KEEP-TOOL |  |
| JACDUMP | unset | – | popiii_m2c.jl:2038, popiii_m2c.jl:2039 | dump the block Jacobian at step k | REMOVE-DEAD | S1: removed |
| NEWTKX | `"0"` | – | popiii_m2c.jl:1980 | one-off Newton debug prints | REMOVE-DEAD | S1: removed (NEWTSEC had zones 112–126 hard-wired) |
| NEWTSEC | unset | – | popiii_m2c.jl:2128, popiii_m2c.jl:2129 | one-off Newton debug prints | REMOVE-DEAD | S1: removed (NEWTSEC had zones 112–126 hard-wired) |
| NEWTDBG | unset | – | popiii_m2c.jl:2182, popiii_m2c.jl:2183 | one-off Newton debug prints | REMOVE-DEAD | S1: removed (NEWTSEC had zones 112–126 hard-wired) |
| NEWTFLIP | `"0"` | – | popiii_m2c.jl:2189 | one-off Newton debug prints | REMOVE-DEAD | S1: removed (NEWTSEC had zones 112–126 hard-wired) |
| TL3DBG | `"0"` | – | PopIIIMaterial.jl:957 | three-level / CHEMSUB debug prints | KEEP-TOOL |  |
| ZONEDBG | `"0"` | – | PopIIIMaterial.jl:959 | chemistry/table step of one zone | REMOVE-DEAD | S1: removed |
| TSTPDBG | `"0"` | – | TitanCore.jl:200 | T/ρ of one zone at every peq!/timstp! | REMOVE-DEAD | S1: removed |
| IFDBG | unset | – | popiii_m2c.jl:2492, popiii_m2c.jl:2493 | zones of a table/network interface after the chemistry | REMOVE-DEAD | S1: removed |
| IFDBGN | `"4"` | – | popiii_m2c.jl:2339 | zones of a table/network interface after the chemistry | REMOVE-DEAD | S1: removed |
| DAMPNEWTON | `"0"` | – | Trbdf2.jl:71 | damped/clipped Newton of the TR-BDF2 stage solves (S6) | REMOVE-DEAD | S1: removed (never adopted; DTERR + VEFIMP took over) |
| DAMPALPHA | `"0.3"` | – | Trbdf2.jl:72 | damped/clipped Newton of the TR-BDF2 stage solves (S6) | REMOVE-DEAD | S1: removed (never adopted; DTERR + VEFIMP took over) |
| DAMPTRIG | `"0.1"` | – | Trbdf2.jl:73 | damped/clipped Newton of the TR-BDF2 stage solves (S6) | REMOVE-DEAD | S1: removed (never adopted; DTERR + VEFIMP took over) |
| DAMPTRACE | `"0"` | – | Trbdf2.jl:75 | damped/clipped Newton of the TR-BDF2 stage solves (S6) | REMOVE-DEAD | S1: removed (never adopted; DTERR + VEFIMP took over) |
| CONV | `"1e-5"` | `1e-3` | order_test.jl:42, popiii_m1.jl:96 … | Newton convergence tolerance — read only by test scripts (order_test, popiii_m1, s5_run, z0_jac_check) | REMOVE-SUPERSEDED | the stack sets CONV=1e-3 but popiii_m2c never reads it: TITAN default conv applies. Delete from run.env (or wire it up deliberately) |

### ionising ray (M3)

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| RAY | `"0"` | – | popiii_m2c.jl:484 | photon-conserving ionising ray from the photosphere (IonRay) | UNCLEAR | off in the stack; IONFS/IONFSIMP (branches ionfs-implicit, ionfs-v2, in progress) build on it: decide after those land |
| RAYTMIN | `"1.5e4"` | – | popiii_m2c.jl:486 | photon-conserving ionising ray from the photosphere (IonRay) | UNCLEAR | off in the stack; IONFS/IONFSIMP (branches ionfs-implicit, ionfs-v2, in progress) build on it: decide after those land |
| RAYXMIN | `"1e-6"` | – | popiii_m2c.jl:507 | photon-conserving ionising ray from the photosphere (IonRay) | UNCLEAR | off in the stack; IONFS/IONFSIMP (branches ionfs-implicit, ionfs-v2, in progress) build on it: decide after those land |
| RAYTEST | unset | – | popiii_m2c.jl:487 | ray stress test (fixed Q, T_eff) | REMOVE-DEAD | remove after the IonFS merge (avoid conflicts now) |

### scripts/popiii_run

| Switch | Default | z300mm | Where | What | Verdict | Note |
|---|---|---|---|---|---|---|
| JLENV | `:?set JLENV to a jul…` | – | fullrun.sh:29, iceulrun.sh:10 | iceulrun.sh/fullrun.sh: julia project, data root, IC shells, threads, numactl prefix, segment length | KEEP-TOOL |  |
| DATA | `:-/zpool/nvme/data/t…` | – | fullrun.sh:29, iceulrun.sh:10 | iceulrun.sh/fullrun.sh: julia project, data root, IC shells, threads, numactl prefix, segment length | KEEP-TOOL |  |
| NZ | `:-200` | – | fullrun.sh:32, iceulrun.sh:13 | iceulrun.sh/fullrun.sh: julia project, data root, IC shells, threads, numactl prefix, segment length | KEEP-TOOL |  |
| TH | `:-4` | – | fullrun.sh:32, iceulrun.sh:13 | iceulrun.sh/fullrun.sh: julia project, data root, IC shells, threads, numactl prefix, segment length | KEEP-TOOL |  |
| PIN | `:-` | – | fullrun.sh:32, iceulrun.sh:13 | iceulrun.sh/fullrun.sh: julia project, data root, IC shells, threads, numactl prefix, segment length | KEEP-TOOL |  |
| SEG | `:-100000` | – | fullrun.sh:32, iceulrun.sh:13 | iceulrun.sh/fullrun.sh: julia project, data root, IC shells, threads, numactl prefix, segment length | KEEP-TOOL |  |
| RUN | unset | – | fullrun.sh:42, fullrun.sh:66 … | run name (fullrun.sh/iceulrun.sh argument) | KEEP-TOOL |  |
| A1X | `:-` | – | fullrun.sh:52 | fullrun.sh: extra settings of stages A1/A2/first B segment; number of B segments | REMOVE-SUPERSEDED | stage chain (fullrun.sh), superseded by iceulrun.sh; stage 3 |
| A2X | `:-` | – | fullrun.sh:56 | fullrun.sh: extra settings of stages A1/A2/first B segment; number of B segments | REMOVE-SUPERSEDED | stage chain (fullrun.sh), superseded by iceulrun.sh; stage 3 |
| B1X | `:-` | – | fullrun.sh:64 | fullrun.sh: extra settings of stages A1/A2/first B segment; number of B segments | REMOVE-SUPERSEDED | stage chain (fullrun.sh), superseded by iceulrun.sh; stage 3 |
| NSEG | `:-0` | – | fullrun.sh:32 | fullrun.sh: extra settings of stages A1/A2/first B segment; number of B segments | REMOVE-SUPERSEDED | stage chain (fullrun.sh), superseded by iceulrun.sh; stage 3 |

## 3. Switches outside the production path (120 names; by script family)

None of these is read by `popiii_m2c.jl` or its includes, so none of them can affect the stack. The verdict applies to the whole
family of scripts (see §5).

| Family (files) | Switches | Verdict |
|---|---|---|
| S4/S5 isolated stars (`s4_a/b/c/d`, `s4_star`, `s4_common`, `s4_test`, `s5_b`, `s5_run`) | ALPHA ATM BURNT CAPDIAG DBGTO DETAIL DT0 DTMAX DTWHY DUMPNOCORE DUMPREJ EVERY FCOUPLING GEOM GROW HSDISC INT ITERTEST ITMAX JACTEST KCAP LAGRHO M MFLOOR MIX MIXTEST N NANREPORT NENV NOCONV NORULE NUC NUCSTEP NSTEPS OUT PREDICT PROFILE PROFILE_AT PROFILE_RELAX PROFILE_X RCAP REGIONCAP RELAX RESDUMP RHO RSUN S4DBG SMAXF STOPZAMS T0 TEND TMAX TRACE VSPLIT ZONEDUMP | UNCLEAR: keep only if the S5 star line (STARFROM) stays |
| Quasi-star / LRD line (`qstar_run`, `QuasiStarSetup`, `qstar_analyze`, `stepdiag`, `newtrace`, `infall_tests`, `test_run_convention`, `lrd_metrics`) | QS_COMPTON QS_CONV QS_DTRIG QS_FREEZE QS_GASROW QS_GRIDTAU QS_GRIDW QS_META QS_MLT_OLD QS_MSG QS_SAVE QS_SET QS_TRBDF2 QS_XSCALE SCAN RUN_CHECK LRD_R LRD_STAR | UNCLEAR: a separate science line; not Pop III |
| Porting/parity harness (`etotal_run`, `phase3/5/6_run`, `phase5_check`, `newton_run`, `energy_ledger`, `order_test`, `z0_jac_check`, `av_test`, `compton_box`) | TITAN_SET DN_SET ETOT_VERBOSE RUNTAG FOLLOW_JULIA PARITY TITAN_NCELL AMP CADV NLIST NREF SCHEMES AVC QL1P QL1S THET NOCOMPTON CONV | REMOVE-DEAD: the port is done; keep one Fortran-parity regression if TITAN_BITWISE stays |
| Early Pop III models (`popiii_m1`, `popiii_m2_onezone`, `popiii_regimes`, `popiii_tests`) | DECADES DUMPAT DUMPRHO FT Q RLIM XIMAX DLNRHO NMAX POPIII_PROJECT POPIII_SERIAL LI | REMOVE-DEAD (m1/onezone/regimes); `popiii_tests` UNCLEAR (is it still a regression suite?) |
| Spectra and analysis (`closed_form`, `line_profile`, `nlte_spectrum`, `snapspec`, `snapfull`, `SpecFull`, `ionfront_*`, `hdchem_check`, `stellareos_test`, `tierc`) | CF_FUV CF_L CF_NLTE CF_STAR CF_TEFF LP_NLTE LP_NZ LP_TAG LP_XI SPECMAX HDLINES OUTDIR SPEC CKPATH EKPATH HELMDIR TIERC_REF | KEEP-TOOL (`snapspec`/`snapfull` feed the run explorer's spectra); the rest UNCLEAR, one by one |

## 4. Dead code paths (reached only through removable switches)

Removed in stage 1 (§8): `Trbdf2.damped_iterate!`, `newton_damped!`, `merit`, `rowscales!`, `scan`, `worst`, `psave!/prestore!`;
`ChemNet.step_sub`, `cdrates`, `attempt_log` (inside `step`), the adaptive loop of `step_ck`; `PopIIIMaterial.chemreplay`, `xss_tl3`,
`T_ss_at_u`; `popiii_m2c.cplscan`, `xadvu!`, `gridmap_step!`, `gridmap_line`, `gridwin_line`, `gmwin_kw`, `gmtau_kw` and all of
`GridMap.jl`; `TitanResidFast.gridwin_rows!` and `upwind`; `TitanCore.tstpdbg`.

Still in the code, to go in later stages:

| Code | Reached only through | Lines (≈) | Stage |
|---|---|---|---|
| `Rezone.jl` + `popiii_m2c.rezone!`, `REZ`, `prerefine!` | REZONE=1, REFINE (Lagrangian stages) | 256 + 90 | 3 |
| `popiii_m2c.radbuild!` | RADFROM (stage A2) | 50 | 3 |
| `envbuild!`, the lgrid 3 → 2 branch | RESTART of a Lagrangian checkpoint with ENVEUL | 105 (shared helpers `envlr` and the fill stay for ICEUL) | 3 |
| `envsplit!`, `envmerge!` | ENVSPLIT, ENVMERGE | 200 | 3 |
| `treset!` | TRESET | 40 | 3 |
| `starfrom!` + `s5_run.jl` family | STARFROM | 25 + ~1000 | 4 (if the S5 line goes) |
| PrimordialGas branch (`PM.PrimordialGas`, RHOEQ, EQEOS, `PrimordialChem` RHFC02 network used only as material) | CHEM=rhfc02 | ~100 in PopIIIMaterial (PrimordialChem stays for its rates and `cool_floor`) | 3 |
| table/network switch rules (the PSWITCH, PSWTAU, PSWSTAY and PSWTEX branches; CHEMFBDENSE; CHEMTMIN) | PSWITCH < 1e30, PSWTAU > 0, PSWSTAY | ~60 | 2 |
| lagged-NLO patches (`nlopc` damping, `nlopc_x`, the NLOFIX floor, the NLOTMAX cut, NLOIMPT extrapolation) in network zones | NLOIT=0 | ~150 | 2 |
| IMPCOOL slope (second transfer call), DTCOOL cap, LINEFRAC reuse | LINEIMP=0 | ~80 | 2 |
| split chemistry patches: CHEMSUB tree, CHEMCKHSS, `echem_pre!/echem_post!` (ECHEMSPLIT), BALSPLIT, the CPLBAND parser | COUPLED=0 or CPLBAND set | ~250 | 2 |
| TITAN θ integrator path (`TD.STEPFN` default, THETA) | INTEGRATOR=theta | small; TitanDriver keeps it for the Fortran parity | 2 |
| `TitanResidFast.FREEZE` / `DEC` branch freezing | only `DampedNewton.jl` (stepdiag, newton_run: not production) | ~30 in the residual | 4 |
| `StellarOpacityLow.jl`, `TitanResidRef.jl` | nothing includes them | 101 + 546 | 4 |
| `TitanResid.jl`, `JuliaSolver.jl`, `PortHarness.jl`, `DampedNewton.jl` | porting harness only | 546 + 214 + 69 + 126 | 4 |
| `LyaTransfer.jl` | LYART=1 (never set) | 222 | 5 (PI) |
| `h2lines_gpu.jl`, `H2LinesKA.jl` | LINEGPU (never set) | ~90 + 493 | 5 (PI) |

## 5. Legacy stage scripts, env files and test/benchmark scripts

- **Stage chain.** `scripts/popiii_run/fullrun.sh` (A1 → A2 RADFROM/REFINE → A3 MSTOP → B remap) and `phaseA.env`,
  `stageA1.env`, `stageA2.env`, `stageA3.env`: superseded by `iceulrun.sh` (ICEUL). Stage 3.
- **`scripts/popiii_run/phaseB.env` (repo copy) is stale.** It is not the stack: it has `CDYN=1e5`, `DTENV=1000,1e5,0.2`,
  `LINEXC2=16` and `ENVCONS=1`, and none of the stack's run.env physics. Every current run keeps its own `phaseB.env`, which differs
  from the repo copy. Stage 2 should replace it with one stack file holding only the settings that do not become defaults (§7).
- **Old run scripts in `scripts/`:** `phase3_runs*.sh`, `phase4_damped.sh`, `phase5_*.sh`, `jfast_runs.sh`, `jfull_runs.sh`,
  `etot_*.sh`, `energy_sweep*.sh`, `legs.sh`, `driver_all.sh`, `run_all.sh`, `run_bench.sh`, `check_exe.sh` (177 lines). They drive
  the porting harness. REMOVE-DEAD in stage 4, together with the harness `.jl` files. `scripts/g2audit/` (G2 audit) is KEEP-TOOL
  until G2BLEND is decided.
- **Julia scripts outside the production path:** 115 files, 13 462 lines.
  - **REMOVE-DEAD, ≈ 5 000 lines (stage 4).** The porting/parity harness (`phase2_run`, `phase3_run`, `phase5_check`,
    `phase6_run`, `verify_phase1`, `TitanResid`, `TitanResidRef`, `JuliaSolver`, `PortHarness`, `resid_check`, `fast_check`,
    `solver_check`, `harness_eos`, `fd_jacobian`, `driver_all`, `legs`, `legcompare`, `Metrics`, `etot*`, `energy_*`, `order_test`,
    `av_test`, `compton_box`, `newton_run`, `DampedNewton`, `eddington_newton*`, `gasrow_check`, `extend_grid`, `to_adaptive`), the
    early models (`popiii_m1`, `popiii_m2_onezone`, `popiii_regimes`), and `StellarOpacityLow`.
  - **UNCLEAR (PI), ≈ 2 600 lines.** The S4/S5 and quasi-star lines.
  - **KEEP-TOOL.** The physics unit tests (`h2lines_test`, `h2levels_test`, `he_chem_test`, `he_eos_test`, `nucnet_test`,
    `stellareos_test`, `stellaropacity_*`, `modelatom_test`, `subzone_test`, `spectrum_test`, `nlte_test`, `hdchem_check`,
    `lopac5_check`, `h2crit_check`, `ionfront_*`) and the post-processing (`runrec`, `snapspec`, `snapfull`, `SpecFull`, `Spectrum`,
    `NLTEH`, `ModelAtom`, `nlte_*`, `export_report_data`, `line_profile`, `closed_form`). Each test should be run once before
    stage 4: any that no longer runs is REMOVE-DEAD.
- **Stage 1 left no test or benchmark script orphaned.** No script in `scripts/` or `julia/` referenced any of the stage-1
  switches or symbols (checked with grep, §6).

## 6. Dependency check of the REMOVE items

**Static.** For each of the 44 stage-1 switches:
1. It is not in the recorded env of z300mm or z300g2 (meta files) and not in the env text of z300mm, z300g2 or z100mm
   (`phaseB.env`, `run.env`, `launch.env`). Checked by script.
2. Its code ran only behind `if <switch>` (or `X[] > 0`, `!== nothing`) with the default branch kept verbatim. The few rewrites:
   - `tab = ptab || (ag === nothing ? X : ag)` with `ag ≡ nothing` became `tab = ptab || X`, which has the same short-circuit.
   - `upwind(uf, hold, …)` with `uf ≡ nothing` became `decide(:up, …)`.
   - `LOGX[] ? attempt_log(…) : attempt(…)` became `attempt(…)`.
   - `fT > 0 || (adv …; return …)` with `fT ≡ 0` became the direct `adv` and return.
3. No removed symbol is referenced anywhere else in `julia/`, apart from `DampedNewton.jl`'s own DTRIG, which is a different
   module. The check was a script over the diff (`symcheck.py`). It found one real reference: `gm_tau`, defined in the GRIDMAP
   block and used by `grideq_driver.jl`. It was restored verbatim. The first test tree (b5cf21b, without it) failed at once with
   `UndefVarError: gm_tau`.
4. No CKENV key was touched (`IC TCMB AVZONE AVQL1 AVCQ2 CHEM NLTE NALI NBIN K5 MCLOUD LRJOIN DLNJOIN DM PSWITCH IMPCOOL DTCOOL
   INTEGRATOR THETA XE0 RAD`), so restarts of the existing checkpoints are unaffected.

**Evidence.** The restart and ICEUL runs in §8.

For the later stages, the REMOVE-SUPERSEDED items that the stack *does* set need their value translated. These are PSWITCH=1e30,
PSWTAU=0, PSWSTAY=0, PSWTEX, NLOIMPT=1, NLOFIX=1, CDYN, STOL, FAILTFAC, CPLXSTOL, STOLCHEM, DTRADMAX=0, ERTAUMIN, ERTAUCONV=0,
TCOOLTH=1, CHEMSUB, CHEMCKHSS, CPLBAND=, ECHEMSPLIT=0, BALSPLIT=0, LINEFRAC=0, REZONE=0, REZEVERY, REZJSPLIT and CONV. For each
of them, §0.1 gives the static reason the setting has no effect on the stack, and each stage's verification below checks it.
Caveat for the .aux part of the bit test in stage 2: NLOIMPT still shapes the lagged NLOPC factors (`nlopc`: qPref, qRref). With
NLOIT those factors are no longer used for network-zone opacities, but they are saved in the checkpoint (`lag.NLOPC`). Removing
NLOIMPT can therefore change the .aux while leaving all 591 state arrays identical. Compare the state, and compare the .aux
field by field.
PSWITCH is the only one in CKENV, and `restart!` already allows it to differ when ENVEUL is set (`allow = ("PSWITCH",)`).
IMPCOOL, DTCOOL, INTEGRATOR, THETA and AVQL1 are in CKENV but unset in the stack. When their code goes, the CKENV tuple loses
them, and `restart!` must skip `aux.env` keys that are no longer in CKENV.

## 7. Clean-up plan (safe, reviewable stages)

Every stage is verified the same way, on a `git archive` of the stage's commit against one of the previous commit:
- **(R)** restart z300mm's `z19z300mmS2` for 100 steps. The state at step 6500 must equal z300mm's own `S3_auto1`: 591 state
  arrays and the .aux.
- **(I)** start from the IC (`ICEUL=1`) and compare the final checkpoints of both trees. Until the cold-IC failure of the G2
  stacks is fixed (§9.9), this uses the zdiagB stack.
- **(G)** stages 2 and 3 also repeat (R) for z300g2 (no GRIDEQ).

The env given to the new tree is the stack with the removed and now-default switches deleted. The comparison is bit for bit, except
where a stage changes a default on purpose. No stage below does that: a superseded patch is removed only where §0.1 shows it has no
effect.

| Stage | Content | Files | Lines removed (≈) | Verification |
|---|---|---|---|---|
| **1 (done, §8)** | 44 dead or experimental switches and one-off debug probes the stacks never set | ChemNet, PopIIIMaterial, TitanCore, TitanResidFast, Trbdf2, popiii_m2c, GridMap.jl (deleted) | **1 141 deleted, 17 added** | (R) + (I) |
| **2** | Make the stack's general methods defaults: COUPLED, CHEMCK, DTERR=0.1, VEFIMP, LINEIMP, NLOIT, NLOPAC, SLIMP=2, SUBZONE=1,2, EDDGMU0, ETOTAL, AVZONE=2, AVTENSION=0, STARLIGHT, SL3, LYABB, LOCALJ, TL3DEEP, HM88, SCVHXE, G2OPAC=2, NLOHBFDEEP, HECHEM, HDCHEM, EPTAB, LINEFAST, the NLTE=hyb line settings, ALI*, the PSWNH dense-interior table. Then delete what they supersede: the table rules (PSWITCH, PSWTAU, PSWSTAY, PSWTEX, CHEMFBDENSE, CHEMTMIN), the lagged-NLO patches, the IMPCOOL/DTCOOL/TCOOLTH/LINEFRAC line lag, the split-chemistry patches (CPLBAND, CHEMSUB, CHEMSUBMAX, CHEMCKHSS, ECHEMSPLIT, BALSPLIT), the G3 caps that DTERR made diagnostic (CDYN, STOL, FAILTFAC, CPLXSTOL, STOLCHEM, DTRADMAX, ERTAUMIN, ERTHICK, DTENV), ERTAUCONV, SLH, CPLIMPQ, AVCAP, SHOCKRAD, GEQREMAP, SLFREEZE, and the default-on bit-identical perf switches' old branches (ZONESWEEP, LNKMEMO, THCHAIN). `restart!` skips aux.env keys that left CKENV. Replace the repo `phaseB.env` with the minimal stack file and drop CONV, FLUXTIECORE, NLOPACJTAU and PSWDEEP. | popiii_m2c, PopIIIMaterial, ChemNet, TitanResidFast, TitanCore, Trbdf2, H2Lines, scripts/popiii_run/phaseB.env | 900–1 200 | (R) for z300mm and z300g2; (I). Translate the env: delete every switch that became a default or was removed. Expect bit-identical, because no removed branch was active (§0.1); any difference is a finding to resolve before merge. |
| **3** | One run from the IC (G7): drop the stage chain and remaps. fullrun.sh, phaseA/stageA*.env, RADFROM/radbuild!, REFINE/prerefine!, REZONE*/Rezone.jl, MSTOP, ENVMS0, ENVDT0, ENVSPLIT/ENVMERGE, TRESET, the lgrid 3 → 2 branch of envbuild!, ENVSINKTAU, ETOTGRAVGRID, the Lagrangian IC path (ICEUL becomes the only start), CHEM=rhfc02/PrimordialGas (RHOEQ, EQEOS, K5 if unused), the IC choices other than abn02red/abn02fin (PI), RAD=0 | popiii_m2c, Rezone.jl (deleted), PopIIIMaterial, scripts/popiii_run | 900 + 256 | (R): the restarted runs never take these paths. (I) with `ICEUL=1` (it becomes the default). The older Lagrangian-lineage checkpoints (t_c128c, c128, …) then no longer restart; they stay on disk. |
| **4** | Drop the old harness and side lines: the §5 REMOVE-DEAD julia scripts, scripts/*.sh, `TitanResid(Ref)`, `JuliaSolver`, `PortHarness`, `DampedNewton` + residual FREEZE/DEC, `StellarOpacityLow`; S4/S5 and the quasi-star lines if the PI agrees | julia/*.jl (≈ 40–60 files), scripts/ | 5 000–8 000 | (R) + (I) (only the production files' includes change); the KEEP-TOOL unit tests run once |
| **5 (after the open branches land)** | RAY/RAYTEST vs IonFS; G2BLEND/H2PF (if adopted, PSWNH, PSWHYS and PSWTMIN go too); LyaTransfer; LINEGPU/H2LinesKA; OLDOPAC, TL3RECKT, TL3NEFULL, H2PPD (adopt or remove) | IonRay, PopIIIMaterial, ChemNet, LyaTransfer.jl, h2lines_gpu.jl | 300–1 000 | per decision; adoption changes results deliberately (a physics change: its own test run, not a bit test) |

## 8. Stage 1, as implemented (commit dee8c76 on `cleanup-audit`)

**Removed: 44 switches** (all REMOVE-DEAD, none set by any of the three stacks). **7 files, 1 141 lines deleted, 17 added.**

| File | Switches removed | Lines |
|---|---|---|
| Trbdf2.jl | DAMPNEWTON DAMPALPHA DAMPTRIG DAMPTRACE (damped stage Newton, its trace and statistics) | −109 |
| ChemNet.jl | CHEMSUBCYC CHEMSUBMAX_SC (sub-cycler `step_sub`, `cdrates`), CHEMLOGX (`attempt_log`), CHEMROBUST CHEMNANDBG, CHEMTRCLIP, CHEMSTAG, CHEMGIVEUP, CHEMCKADAPT; the internal TRACEALL (set only by chemreplay) | −236 / +5 |
| PopIIIMaterial.jl | CHEMREPLAY CHEMREPLAYTRACE (`chemreplay`), PSWAGREE (+ PSWTAUDEEP, `xss_tl3`, `T_ss_at_u`), PSWRAD, ZONEDBG | −158 / +1 |
| popiii_m2c.jl | EDDREAL EDDTAU EDDRELAX, EDDCHI, EDDDBG, NEWTKX, NEWTSEC, NEWTDBG NEWTFLIP, JACDUMP, XADVU (`xadvu!`), CPLSCAN CPLSCANNLO (`cplscan`), IFDBG IFDBGN, FLUXTIEDIAG, PSWAGREE wiring, GRIDMAP GMJUMP GMTAU GMWIN GMSMOOTH (`gridmap_step!`, `gridmap_line`; `gm_tau` kept for GRIDEQ), GRIDWIN (`gridwin_line`, its Δt hook) | −425 / +2 |
| TitanResidFast.jl | GRIDWIN (`gridwin_rows!`, `gridmoves`), UPFREEZE (`upwind`), CPLTRACEX, CPLDBG | −66 / +8 |
| TitanCore.jl | TSTPDBG (`tstpdbg`), FLUXTIEDIAG | −18 / +1 |
| GridMap.jl | (whole file; only GRIDMAP used it) | −129 |

Kept on purpose, although never set: OLDOPAC, TL3NEFULL, TL3RECKT, H2PPD (correctness options; PI decides), RAYTEST (avoids a
conflict with the IonFS branches), REZNULL/REZDBG (they go with Rezone.jl in stage 3), and the diagnostics marked KEEP-TOOL.

**Verification** (test trees are `git archive` of 1782928 and of the stage-1 commit under `runs/t_cl1`; cores 72–79, membind 9,
8 threads; the env is composed exactly as `iceulrun.sh` composes it from z300mm's `phaseB.env` + COMMON + `run.env`; the only
other change is EQEOS_CACHE, pointed at a copy of the same file):

1. **(R) restart, `runs/t_cl1`.** Restart of `z19z300mmS2` (z300mm, step 6400, GRIDEQ moving mesh, 290 zones) for 100 steps on
   each tree. The step-6500 checkpoints of the base tree, of the stage-1 tree and of z300mm itself (`S3_auto1`, copied before it
   was overwritten) agree: **591/591 state arrays identical, the .aux identical field by field and byte for byte**. This holds for
   both the `_auto1` and the final checkpoint. The snapshot records (steps 6450 and 6500) are bit-identical. The run logs are
   identical line for line once timing lines and the checkpoint name are filtered out. The 100 steps exercise the full stack:
   DTERR with 19 Newton-failure rejections and a FSAL restart, COUPLED in 243 zones, NLOIT, LINEIMP, hybrid ALI (363 calls),
   VEFIMP, SUBZONE, SLIMP, GRIDEQ and HOT.
2. **(I) ICEUL start.** With `ICEUL=1` from the IC (abn02red, 128 IC shells, 290-zone mesh):
   - **z300mm stack (`runs/t_cl2`) and z300g2 stack (`runs/t_cl3`).** Both trees fail identically at step 0: NaN in zone 1 in the
     first Newton iteration, 60 attempts, no checkpoint written. Their outputs are identical line for line once timing lines are
     filtered out. The failure is already in 1782928 (§9.9).
   - **zdiagB stack (`runs/t_cl4`).** This is the z300g2 stack without the G2 material, SUBZONE and SLIMP, i.e. the last stack that
     ran from the IC. Each tree ran 30 steps to t = 1.03e12 s, through the IC build (ICEUL envbuild, ENVCONS fill), DTERR (one
     rejection), COUPLED, NLOIT, LINEIMP and VEFIMP. The final checkpoints agree: **591/591 state arrays identical, the .aux
     identical byte for byte**, and the logs are identical.
   - Until the cold-IC failure is fixed, the (I) check of every later stage uses the zdiagB stack. When it is fixed, (I) moves to
     the z300mm stack.

First attempt (recorded for honesty): a stage-1 tree without `gm_tau` (b5cf21b) failed at once in `grideq_step!`. `gm_tau` was
defined in the GRIDMAP block. It was restored verbatim before the verified commit. Its old tree and output are in
`runs/t_cl1/outdated/`.

**Merge outlook.** Merging stage 1 into the open branches gives two trivial conflicts (checked with `git merge-tree`):
- `ionfs-v2` adds NEWTONERR next to the removed damped-Newton block in `Trbdf2.newton!`. Keep both.
- `g2-domain` adds the G2BLEND line after the `tab = …` line that stage 1 simplified. Keep both.

## 9. Open questions for the PI

1. **CONV.** The stack's `CONV=1e-3` was never applied (Newton runs at `conv = 1e-5`). Delete it from the env, or make popiii_m2c
   read it? Reading it changes results.
2. **GRIDEQ's fmove Δt cap (GEQC[5] = 0.5) is inactive under DTERR,** like every DT_HOOK. Is that intended for the moving mesh,
   or should DTERR honour a mesh-motion bound?
3. **Side lines.** Keep or drop the S4/S5 isolated-star line (STARFROM, s5_run, ≈ 1 000 lines) and the quasi-star/LRD line
   (≈ 1 600 lines)?
4. **Never-adopted correctness options.** Adopt (results change) or delete: TL3RECKT, TL3NEFULL, OLDOPAC, H2PPD, ADVREL.
5. **Methods on hold.** LyaTransfer (LYART, 222 lines) and the GPU line kernels (LINEGPU, H2LinesKA, ≈ 600 lines): keep?
6. **G2BLEND/H2PF (g2-domain).** If adopted, they replace the last table rule (PSWNH, PSWHYS, PSWTMIN) and TABAUDIT. Decide
   before stage 2 finalises the material defaults.
7. **RAY vs IonFS.** Does the IonFS work replace the driver's RAY path (RAY, RAYTMIN, RAYXMIN, RAYTEST)? IonRay's module
   functions stay either way.
8. **ICs.** Which ICs stay as supported starts: abn02red only, or also abn02fin, abn02z19 and p100? Older production lines
   (t_c128c, t_nolte2, c128, …) stop being restartable after stage 3; their checkpoints stay on disk.
9. **No current stack starts from the IC.** With `ICEUL=1` on z300mm's 290-zone mesh, both the z300mm stack and the z300g2
   stack fail at step 0: NaN in zone 1 in the first Newton iteration, 60 attempts. This happens identically on 1782928 and on
   stage 1, so it predates the clean-up. It is the same failure as all three z100mm tries. The zdiagB stack (z300g2 without
   SCVHXE, HM88, G2OPAC, NLOHBFDEEP, TL3DEEP, the PSW*/PSWNH table settings, SUBZONE and SLIMP) does start (§8 (I)). So one of the
   G2/SUBZONE/SLIMP additions breaks the cold IC. Bisecting it is the next step for the from-IC goal (G7); it is not part of this
   clean-up.

## 10. Stages 2 and 3 and the result-changing adoptions (PI decisions, 16 Nov)

**PI decisions**, carried out on `cleanup-audit`:
1. CONV is dropped from the env.
2. DTERR respects the moving mesh's motion bound.
3. The S4/S5 line and the quasi-star/LRD line are kept.
4. TL3RECKT, TL3NEFULL, OLDOPAC, H2PPD, ADVREL and SLFREEZE are deleted, together with the other never-adopted patches.
5. LyaTransfer and the GPU line kernels are kept.
6. G2BLEND and H2PF become the default once they are merged and the IC start works.
7. The RAY path is kept.
8. The ICs kept are abn02red (main), abn02fin and abn02z19.

The branch first merged `popiii-s6` 38f886e (which includes 7547bab slimp-ic) and `g2-domain` 00a21f1. One conflict in
`chem_core!` was resolved by keeping stage 1's simplified rule and adding the G2BLEND line. The result is **cedfdc6 = M0**, the
reference for the bit tests.

### 10.1 Commits

| Commit | Content | Lines | Changes results? |
|---|---|---|---|
| ddd0c55, cedfdc6 | merges of popiii-s6 and g2-domain | – | no (M0 reproduces z300mm, below) |
| **52a947f** stage 2a | `julia/popiii_stack.jl`: the z300mm stack as the defaults of popiii_m2c (a run sets only its mesh/output). Hard-wired methods: COUPLED, LINEIMP, NLOIT; DTERR required. Removed as superseded: IMPCOOL/IMPALWAYS/DTCOOL/TCOOLTH/LINEFRAC (LINEIMP), the lagged NLOPAC factors and NLOIMPT/NLOFIX switch/NLODAMP/CPLNLOX/NLOTMAX (NLOIT), ECHEMSPLIT and CPLBAND (COUPLED everywhere), STOLCHEM/CPLXSTOL/DTRADMAX/DTENV (DT hooks, diagnostic only under DTERR), the ERTAUCONV floor (VEFIMP). Removed as never adopted: TL3RECKT, TL3NEFULL, OLDOPAC, H2PPD, ADVREL, SLFREEZE, CPLIMPQ, AVCAP, SHOCKRAD/SHOCKDU, GEQREMAP. `restart!` skips checkpoint keys whose switch is gone. | +90 −429 | no |
| **1c03845** stage 3 | one run from the IC: the start is the IC on the fixed log mesh (formerly ICEUL=1; the switch is gone) or a RESTART of a fixed-mesh checkpoint, and STARFROM is kept. Removed: RADFROM (radbuild!), REFINE, REZONE/REZ* and Rezone.jl, the lgrid 3 → 2 remap (envbuild!), ENVSPLIT, ENVMERGE, TRESET, MSTOP, the RHFC02 PrimordialGas material (CHEM=rhfc02, RHOEQ, EQEOS), and the ICs p100, abn02 and abn02hs. Removed from scripts/popiii_run: fullrun.sh, phaseA.env, stageA1–3.env; iceulrun.sh is the launcher and phaseB.env is empty. | +37 −852 | no |
| **8e49124** (ii) | `Trbdf2.DT_BOUNDS`: hard bounds that DTERR respects; GRIDEQ's fmove cap is registered there | +15 −2 | **yes** |
| **ea2bd2a** (ii) | G2BLEND always on (the table only where the non-ideality weight is complete) and H2PF a stack default. Removed: PSWITCH (NetGas pswitch = Inf), PSWNH, PSWHYS, PSWTMIN, PSWTEX, PSWSTAY, PSWTAU, CHEMFBDENSE. PSWRAMP is kept. | +14 −65 | **yes** |
| **283c43e** | `GRIDEQ=off` for the fixed mesh (the stack default is 1); `docs/POPIII_STACK.md` | +56 −1 | no |

CONV needs no code change: popiii_m2c never read it. It is gone from the minimal stack file.

### 10.2 Bit-identity of stages 2a and 3 (env translated to the new defaults)

The test trees are `git archive` exports under `runs/t_cl5`. They ran on cores 72–79 (membind 9) with 8 threads.
- **Old env** (on M0): z300mm's `phaseB.env` + COMMON + `run.env`, as iceulrun.sh composes them.
- **New env**: an empty `phaseB.env` and `run.env` = `SNAPEVERY=50`, which equals the stack's SNAPEVERY and only overrides the
  test launcher's COMMON. For z100mm the new `run.env` adds its mesh: `ENVEUL=0,0.232584376,0.232584376,1e30
  ENVCORE=0.1,0.232584376,1e30`.
- `envcmp.py` lists every key whose effective value differs between the two. The keys that new code still reads are BALSPLIT=0
  (its default), CONV (read only by test scripts) and REZ* (unused with REZONE=0). Every other difference is a removed switch.

| Test | M0 (old env) | stage 2a (new env) | stage 3 (new env) |
|---|---|---|---|
| (R) restart z19z300mmS2 → step 6500 (100 steps, GRIDEQ, 290 zones) | = z300mm's own S3_auto1: 591/591 arrays, .aux bytes identical | = M0 and = z300mm: 591/591 arrays, .aux identical field by field | = z300mm S3_auto1: 591/591 arrays, .aux identical field by field |
| (I) start from the IC with the z100mm stack (100-zone mesh, GRIDEQ from the IC, 30 steps to t = 7.0e12 s) | runs (the slimp-ic fix makes the from-IC start work, §9.9 resolved) | = M0: 591/591 arrays, .aux identical field by field | = M0: 591/591 arrays, .aux identical field by field |

The .aux files of the new tree are not identical byte for byte. The only difference is the iteration order of the `env` Dict
(the same 8 keys and values; the defaults are injected into ENV in a different order). Every other field serializes to the same
bytes.

### 10.3 Result-changing adoptions: test runs

Each restart below starts from z300mm's `z19z300mmS2` (step 6400) and runs 200 steps with the new stack env. The reference is
z300mm itself over the same 200 steps (`seg_z19z300mmS3.out`, `perf.csv`), which is bit-identical to M0. The from-IC runs use the
z100mm stack for 100 steps. Their reference is z100mm's own production run (on 7547bab, old env). The per-step Δt is taken from
the DTDIAG "line cap … dtime before" line, and the rejections from the DTERR reject lines.

| Run | Code | Steps | Δt median / 10 % / 90 % [s] | steps with Δt < 1 s | Newton failures (error rejections) | Σ Δt [s] | at the last step |
|---|---|---|---|---|---|---|---|
| z300mm (reference) | 1782928 | 6400–6600 | 4.10e3 / 1.31e2 / 9.12e4 | 10 | 22 (0) | 5.24e6 | t 2.3893506e5 yr, M★ 0.1771 M☉, R★ 28.78, T_c 1.024e5 K, wall 788 s / 100 steps |
| fmove in DTERR (`t_cl7/out_rfm`) | 8e49124 | 6400–6600 | 4.10e3 / 1.12e2 / 9.36e4 | 10 | 21 (0) | 5.28e6 | t 2.3893507e5 yr, M★ 0.1773, R★ 28.87, T_c 1.024e5, wall 791 s; **the bound set Δt in 2 of 200 steps** |
| + G2BLEND + H2PF (`t_cl8/out_rg2`) | ea2bd2a | 6400–6600 | 3.51e3 / 1.11e2 / 8.75e3 | 10 | 23 (2) | 7.95e5 | t 2.3893493e5 yr, M★ 0.1689, R★ 27.61, T_c 1.020e5, wall 892 s; network zones 243 → 262 |
| z100mm (reference, from the IC) | 7547bab | 0–100 | – | – | 24 divergences | – | t 2.33435e5 yr, log n_c 8.846, T_c 668.3 K, wall 352 s |
| fmove, from the IC (`t_cl7/out_ifm`) | 8e49124 | 0–100 | – | – | 24 | – | **identical to z100mm** (t, n_c, T_c, divergences): the bound did not bind |
| + G2BLEND + H2PF, from the IC (`t_cl8/out_ig2`) | ea2bd2a | 0–100 | – | – | 22 | – | t 2.34377e5 yr, log n_c 9.003, T_c 685.8 K, wall 350 s |

**fmove in DTERR.** The change is stable. In the current state of z300mm (0.17 M☉, the photosphere front at zone ~78) the mesh
moves slowly: the bound set Δt in 2 of 200 steps, and in none of the first 100 steps from the IC. It did not change the stability
picture. Δt still collapses to ~1e-4 s every ~120 steps, in 10 of 200 steps, as in z300mm. The trace shows E_r of zone 78 limiting
each collapse, followed by recovery within ~20 steps. **This collapse is a radiation/front effect, not mesh motion.** It is the
main remaining obstacle to "fast and stable".

**G2BLEND + H2PF.** The change is stable. From the z300mm checkpoint, the smooth blend moves 19 zones of the dense interior from
the table to the network (243 → 262 network zones; all of them coupled). The restart then runs with a typical Δt like z300mm's
(median 3.5e3 s against 4.1e3 s) and the same number of sub-second steps. But the 1e5 s steps after the collapse recoveries, which
z300mm reached by step 6580, have not come back within 200 steps (90 % quantile 8.8e3 s against 9.1e4 s), so the physical time
covered is 7× smaller. The cost per step is 13 % higher (892 against 788 s per 100 steps). At about the same physical time, the
stellar quantities agree with z300mm to the level of the material change: T_c 1.020e5 against 1.022e5 K, and R★ 27.6 against
27.1–28.8 R☉. M★ (the mass inside the first interface falling faster than 2 km/s) reads 0.1689 against 0.1677 M☉ at step 6500,
because the material moved the infall boundary by about one zone. From the IC the collapse is slightly faster (log n_c 9.00
against 8.85 at step 100, T_c 686 against 668 K), with 22 instead of 24 Newton divergences. This is the H2PF/HM88 network EOS in
the cold gas.

**Open question.** The longer test needed before G2BLEND can be called "faster" is whether its Δt recovers to the 1e5 s steps
over a full segment (1000 steps). The from-IC run would have to be continued until the core forms (≈ 4000 steps for z100mm).
Neither fits in this session's core budget. The adoption is in its own commit (ea2bd2a), so it can be reverted alone.
