# Pop III prototype, Stage 1 (branch popiii-proto, 24 Nov 2026)

## What is merged
subconv-mg 1f9d063 (M1-M6 convection diagnostics, MLTBUOY, TDC), tdc-henyey 9026276 (M7), mgprod 063452c (multigroup, a1a2ff4 true-residual convergence + exact J.v, hydrogen logit 54814fc, eddfac/outer-boundary g fixes), prodpkg cc7b8dc (TitanStack image, PopIIIRun.jl + popiii_run_body.jl), gasheat be26fad (docs 31, julia/photoheat_check.jl).
The cc7b8dc merge moved popiii_m2c.jl into PopIIIRun.jl/popiii_run_body.jl; our 15 hunks (CONVDIAG, MLTBUOY, TDC, ETOTLOG, LRESPTEST, MG radsolve, xn arguments) were ported by patch; CONVDIAG and TDC are typed globals set by the body.

## Stage-1 configuration
Defaults of julia/popiii_stack.jl = zie1/zie3 grey production: their run.env (later duplicate keys win: IC=abn02red, SUBZONE=0, PSWITCH=1e30, CDYN=3e5, ECHEMSPLIT=0, CONV=0.01, NLTE=ali) + NKRYLOV=1, TITAN_NOLIB=1, TITAN_MGR=600; the stack's own TCMB=54.5, CLOSJ=1, GEQRAMP, PRINTEVERY=100 also applied in production (meta files). Changed vs the old stack: CONV 1e-3 -> 0.01, DTERR 0.1 -> 0.03, GRIDEQ 1 -> off, SUBZONE 1,2 -> 0, MATTHREADS 8 -> 16; added NKRYLOV=1, SLABFLUX=1, SLABTIE=1. zie1 and zie3 agree on all physics/solver settings (they differ only in ENVEUL/ENVCORE). Per-run files: scripts/proto/run_zie1.env, run_zie3.env (mesh, IC, MCLOUD). TIMING=1 DTDIAG=1 are launcher settings (gate_run.sh sets them). Radiation grey (MGROUPS unset), convection production MLT; TDC/MLTBUOY/CONVDIAG present, off.
28 run.env keys are read by no code (ALIGROW ALISHRINK ALIWMIN CPLBAND CPLCLAMP CPLXSTOL DTRADMAX ECHEMSPLIT EPTAB ERTAUCONV FLUXTIECORE LINECOLL LINEFRAC LINEIMP LOCALJWMIN NLOFIX NLOIMPT NLOIT NLOPACJTAU NLTE PSWDEEP REZEVERY REZJSPLIT REZONE STOLCHEM TCOOLTH ...): not in the Stage-1 defaults.

## Deletions (this stage)
- ETOTAL/ETASW/ETASWTH/ETASWTIE/ETOTGRAVGRID/ERLAB/ECOMPAT error guard in popiii_run_body.jl and in noh_test.jl, phase6_run.jl, selfgrav_energy_test.jl (dual-energy path already removed 21 Nov, e3d5fab-era; nothing else reads the keys).
- Also, each verified by reading the code path: SLH (pH is always nothing; it already errored under SL3=1, which Stage 1 sets); TIEFRAME env read (TX.TIEFRAME is set unconditionally by the ECOMPAT block).
Inventory updated (docs/POPIII_SWITCH_INVENTORY.md, last section).

## Gate (logs docs/logs/proto/stage1/; runs in /zpool/.../proto/stage1)
Reference = cc7b8dc code (the pkgstart 1ed14cf runs are the same code: `git diff 1ed14cf cc7b8dc` = switch_to_prod.sh only), include path, zie1 run.env + empty phaseB.env, 8 threads, SNAPEVERY=5, compare_snaps.jl (isequal).
1. zie1 from the IC, 40 steps, merged popiii-proto vs cc7b8dc: NOT identical (201 of 300 fields; first differing step 1: u 5e-10, m 9e-14, rho 2e-13, T 3e-13 relative; by step 2 T up to 0.13 relative in a few zones, chaotic growth). Same result before the cc7b8dc merge, and for zie3 (200 of 300).
2. 50-step grey restart from z19zie1S16 (prod_ck copy): NOT identical (219 of 300; step 63335: dt differs, u 1e-7, rho 3e-10, T 2e-10).
Cause (bisected on 1-step IC runs, same config): a1a2ff4 is bit-identical to the reference; 54814fc (hydrogen ionisation unknown = logit ln(x_H+/x_HI), merged in 9a4c8d3) is the first differing commit and the later mgprod commits do not change it. popiii-proto equals mgprod 063452c bitwise (zie1 IC, 3 snapshots, 0 differences; MGROUPS=4 restart from rp_k650, 40 steps, 8 snapshots, BIT-IDENTICAL). The logit changes the grey path (unknown, DTERR/Krylov/maxdel weights), so bit identity with cc7b8dc is not attainable on mgprod-based trees without reverting it. Our convection code (off) adds no difference.
Dropped by instruction (Tom, prune): zie3 restart, image-path (TITAN_PKG) run, NKCHECK run; the TitanStack image for this tree did build (image_build.log, 373 s) but was not exercised. zie3 IC 40-step proto log kept.

## Bit-identity vs production
The Stage-1 differences vs cc7b8dc come from mgprod 54814fc (hydrogen-ionisation logit unknown), which the multigroup session verified is an equivalent change of variable in grey (LG check, /zpool/nvme/data/tabel_scratch/titan_popiii/lg_runs): at fixed dt (SDFIXDT, 2 steps) 76856d0 vs cc7b8dc agree to round-off for CONV 1e-2...1e-8 (zie1 S16: T 1e-10-7e-10, rho ~5e-10, x_H+ <= 3e-9; IC: T <= 2e-12). The free-dt differences are the step controller (min(x_H+, x_HI) DTERR weights change the proposed dt by ~5 %; at the IC cc7b8dc fails Newton at 2.4e7 s and accepts 6e6 s, the logit tree accepts 2.4e7 s). Hence the grey bit-identity reference for gate 2 is popiii-proto itself. Tight-CONV runs need DTERR = 1e3*CONV and NKRYLOV=1,CONV,60.

## Open lag
Grey Lyman-Werner heat GASQ is lagged by one step (exact in MG); measured 7e-10 relative, unchanged.

## Deferred switch-pruning candidates (not done)
BALSPLIT (CHEMB: TitanResidFast, PopIIIMaterial, PopIIIRun, in the checkpoint lag tuple), XADV2 (alias of SL3), CPLLOCAL, THCHAIN, INTEGRATOR, LTRAN, SDLAG/SDFRZ*, STEPKX, ZEROSTEP, GEQ* (inert with GRIDEQ=off), the 28 unread run.env keys above, K5, EXCAP. Not touched: NKCHECK, DTERRDIAG, EIGDIAG/EIGSTEP, SWDIAG, HELDDUMP, CONVDIAG, TITAN_PKG/envguard, MGROUPS, TDC, MLTBUOY.

## Stage 2 (24 Nov 2026): TDC is the only convection path

### Changes
- With HOT=1 (convection on) w is always the last Newton unknown (neqn+1); TDCPAR default 1,1,1,1,0.01,1,0, Henyey losses (M7), network linear-response thermodynamics (M6). `TX.TDC_ON[]` is kept as the internal "w installed" flag (set with CONVECTION in popiii_run_body.jl:783); the env switches are gone.
- Species mixing (StarMaterial.mix!, StarMaterial.jl:~551): D = alpha_D l w from the converged TDC state (was v l/3).
- DTERR: w is a combination (like r, x_me) in `combos`/`old_from!`/`refilter*` and a weighted entry in `wnorm` (Trbdf2.jl): delta_w against rtol*max(|w|, TDC_WFLOOR c_s), the convention of u and of the Newton measure; error name "w".
- config_tdc.env = config_grey.env + CONVDIAG=1 (switch lines gone); COLDMLT removed from popiii_stack.jl and config_grey.env. envguard.jl is a source scanner (no key list): nothing to edit. pc.jl: precompile lines of `_mlt_pre!`/`convcoef` removed.

### Deletions (before-state line numbers of 700e60e)
TitanResidFast.jl: `convcoef` 405-466 and `Ctx.cv`; MLT block of `residual!` 1854-1999 (incl. MLT_HS blend, MLT_AD lagged path, MLT_BUOY branch, tau gate); `_mlt_pre!` 1211-1242; `TAU_MLT`, `MLT_AD`, `MLT_HS`, `MLT_NOINFALL`, `MLT_BUOY`, `MLT_RCAP`, `MLT_HSDISC` 896-955 (RCAP form and the discrete DlnP_hs, `lph`, are unconditional in `_tdc!`); WS fields mnad/mcp/mdl/mkp; SW_MLT* registrations (1465-1466); conv_diag_pass w_tau/w_inf/tau/cvw. PopIIIRun.jl:633-637 (TDC global :38, MLTBUOY, COLDMLT, RCAP/HSDISC settings); popiii_run_body.jl:319-320 (TDC env), TitanFortran.jl:116 (meqn now keyed on HOT=1); PopIIIMaterial.jl: `COLDMLT`, `cold_thermo` (233-247), `material_thermo(::NetGas)` with the 1e30 placeholder (309-318); CONVDIAG fields hs_w, nad_used, nab, w_tau, w_inf, tau, w_th (PopIIIRun.convdiag_snap, runrec.jl cv_* exports; fpot_th := fpot). stepdiag.jl QS_MLT_OLD line. Inventory: docs/POPIII_SWITCH_INVENTORY.md, last section.

### Gate (logs docs/logs/proto/stage2/; runs in /zpool/.../proto/stage2; cores 64-87; Stage-1 = 700e60e tree, same settings, CONVDIAG ETOTMON on)
| gate | Stage 1 (700e60e) | Stage 2 |
|---|---|---|
| T0 tdc_tests t4-t10, mesa | - | all PASS, exit 0 (9,5,3,1,3,16,2,15) |
| 1 NKCHECK grey, front S6_auto1 (AD-FD 2-norm) | 1.5e-10 (ref) | 3.0e-10; non-w rows <= 2e-8 except F_r 3e-3 (as before: noisy F_r rows); w rows: front zones 3e-16, worst w row 1.0 (zone 98, w~0 kink/noise) |
| 1 NKCHECK grey, S16 (iteration 1) | - | 2.9e-10; w worst 1.9 (zone 50) |
| 1 NKCHECK MGROUPS=4, c5_dA | 2.306e-7 (ref) | 4.1e-7 (T/E_g4 6.6e-7 zone 21; w worst 1.0; F_g2 row z101 noise) |
| 3 zie1 S16, 50 steps | 50 acc, 1 Newton rej, dt med 6.55e6 s, t_sim 1.97e8 s, F_c/F_tot 0, 629 s total wall | 50 acc, 1 Newton rej, **dt med 23 s, t_sim 894 s**, F_c/F_tot 1.0, P_t/P 0.042, w_max 4.2e7 cm/s, 2276 s wall |
| 3 zie3 S8, 50 steps | 50 acc, 0 rej, dt med 8049 s, t_sim 4.49e5 s, F_c/F_tot 0, 772 s wall | 50 acc, 0 rej, dt med 8260 s, t_sim 4.81e5 s, F_c/F_tot 0.97, P_t/P 0.012, 1249 s wall |
| 4 ETOTMON |sum imb|/sum thr, S16 / S8 | 2.1e-3 / 1.1e-3 | 6.5e-4 (max G 7e-3) / 1.2e-3 |
| 5 zie3 fp3a1 MG, 40 steps | 40 acc, 0 rej, dt med 9248 s, 1743 s wall | 40 acc, 0 rej, dt med 9312 s, ledger 3.0e-3 (same), 2716 s wall |
Wall includes ~400-700 s start-up (no TITAN_PKG image used); per-step wall not separable here.

### Findings
- Gates 5 and zie3 S8 pass (accepted, dt, ledger). Gate 1 passes on non-w rows (grey 2x the Stage-1 2-norm, MG 1.8x: both dominated by w rows where w ~ 0 or at an FB kink; w rows are exact at the front zones).
- **zie1 S16 collapses to dt ~ 10-60 s** (Stage 1: 6.5e6 s). Not caused by the w DTERR term (TDCWERR=0 run: dt median 33 s) nor by the mixing patches (bridge/tmax off: 21 s). The checkpoint has no w (starts from w = 0) and Stage 1's MLT had no convective flux at all in this state (F_c/F_tot 0) whereas TDC convects the core (F_c/F_tot ~ 1, w to 4e7 cm/s, limiting variables w / u / F_r / E_r of central zones): a restart transient or a real timescale. 50 steps cannot tell; needs a longer S16 run from a TDC state (w carried in the checkpoint).

### Open items
- zie1 S16 dt collapse above (longer run from a TDC checkpoint; carry w in checkpoint is implemented, `aux.tdc`).
- mix! `bridge=4`/`tmax=0.02` patches kept: switching them off changed nothing on S16 but S16 is in the transient, no test on a state with a moving core boundary (gate 3 zie3 has no hot core effect). Delete once a TDC state is tested.
- Non-coupled network species (Y of the stellar network outside the coupled x_H+, x_H2): not mixed by TDC rows, only by mix!.
- Legacy MLT drivers break: s5_run.jl, qstar_run.jl (TX.MLT_* assignments), QuasiStarSetup/ProblemSetup/extend_grid/newtrace (CONVECTION without w), popiii_regress SUBC, stellareos_test uses only material_thermo (kept). scripts/subconv diag_*.jl, abnfig.py read removed cv_* fields.
- No TITAN_PKG image built for this commit; NKCHECK w rows at the FB kink are not FD-testable (tol.json 'w rows info').
