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

## Stage 2
Switch the convection to TDC (+ Henyey) as default and gate it on the zie1 front; re-run the gate against the prodpkg commit the MG session names.
