# The supported Pop III configuration

titan.jl supports one configuration for Pop III protostar formation, accretion and H II break-out. It starts from the ABN02 cloud
and runs as one simulation from the initial condition on one mesh. This configuration is the set of defaults of
`julia/popiii_m2c.jl`, listed in `julia/popiii_stack.jl`. The 16 Nov 2026 clean-up made it so
(docs/POPIII_SWITCH_AUDIT2.md, docs/POPIII_GOAL.md). A run sets only what is its own, such as the mesh, the cloud or the output.

## Launch from the initial condition

```bash
export JLENV=<julia project with the repository's packages>
R=/zpool/nvme/data/tabel_scratch/titan_popiii/runs/<run>; mkdir -p $R
cp <an existing eqeos_he_xe.jls> $R/      # the SCVH95 + He + x_e table cache (built and written there if missing)
cat > $R/launch.env <<'X'
NZ=128
TH=8
PIN="numactl --physcpubind=<8 cores> --membind=<their node>"
SEG=1000
X
echo "" > $R/run.env                       # nothing: the 290-zone stack mesh; or e.g. the 100-zone mesh below
scripts/popiii_run/iceulrun.sh <run>       # segments of SEG steps; the first starts from the IC, each later one restarts the last
scripts/popiii_run/iceulrun.sh <run> <checkpoint>   # continue from a checkpoint
```

- The 100-zone mesh (z100mm) goes in `run.env`: `ENVEUL=0,0.232584376,0.232584376,1e30 ENVCORE=0.1,0.232584376,1e30`.
- A fixed mesh (no moving mesh) goes in `run.env`: `GRIDEQ=off`. The value `GRIDEQ=0` is something else: it selects the
  prescribed-target mode.
- The other ICs are `IC=abn02fin` and `IC=abn02z19`. Set them in `run.env`, together with a matching `MCLOUD`.

`iceulrun.sh` writes the files of each segment into the run directory: `procedure.log`, `seg_*.out`, the m2c log, `snaps_*.jls`,
`meta_*.txt` with the effective settings, `perf.csv`, and the checkpoints in `ck/`. A restart checks the settings of the
checkpoint (CKENV). It refuses Lagrangian checkpoints and checkpoints without radiation, i.e. those from the old stage chain.

## The stack

These are the defaults (`julia/popiii_stack.jl`). The methods marked † are no longer switches: they are the only code path.

| Part | Method | Settings |
|---|---|---|
| backend | pure-Julia TITAN | `TITAN_NOLIB=1 TITAN_MGR=600` |
| IC | ABN02 red output (9 Myr), 7137 M☉, T_CMB(z = 19) | `IC=abn02red MCLOUD=7137 TCMB=54.5`, radiation from the IC (`RAD=1`) |
| mesh | log mesh from a regular centre (one central zone to 0.1 R☉, Δln r 0.0794, 290 zones), mass-conserving fill; TITAN's implicit grid equation (GRIDEQ mode 1, iterate-only monitor) | `ENVEUL ENVCORE ENVLIN=1 ENVCONS=1 GRIDEQ=1 GEQRAMP GEQC GEQW` |
| material | ChemNet network (H, He, Li, D) with HM88 occupation probabilities and the H2PF partition function, the three-level H atom (TL3DEEP), network continuum opacity at the iterate (NLOIT†, G2OPAC=2, NLOHBFDEEP); SCVH95+x_e table only where non-ideal (G2BLEND†, smooth blend); stellar EOS + nuclear network + MLT in the hot core | `CHEM=netmodernli HM88 H2PF SL3 TL3DEEP NLOPAC NLOPACJ G2OPAC=2 NLOHBFDEEP SCVHXE HOT=1 COLDMLT=1e4` |
| chemistry | coupled implicit: x_H⁺, x_H₂ as Newton unknowns in every network zone (COUPLED†), reduced CK network, HD, He | `CHEMCK CHEMCKHSS CHEMSUB CHEMFBNET=16 CPLCLAMP=3 HDCHEM HECHEM` |
| H₂ lines | ALI transfer every step (no escape-probability or hybrid mode; warm start from the previous step's populations), line cooling at the iterate (LINEIMP†), grey continuum absorption of line photons returned to the gas, incident CMB + photosphere | `LINEINC LINESRC=zone LINEQUAD=p2 LINECORE=8 LINEFAST` |
| star light | H⁻, LW, three-level photo-processes from the photosphere, at the iterate | `STARLIGHT SLBLEND SLEB SLWFLUX SLSHIELD=sobolev LYABB LOCALJ SLIMP=2` |
| radiation | grey E_r, F_r; Eddington factors at every Newton iterate, Auer sphericality form, CMB as background (VEFIMP); sub-zone front closure; flux tie in the thin outer zones | `VEFIMP EDDGMU0 EDDSCAT SUBZONE=1,2 FLUXTIE=0.1,1` |
| energy | total-energy row (dual-energy switch), exchange cap | `ETOTAL ETASW=0.6 ETASWTH ETASWTIE EXCAP=1e8` |
| hydro | zone-width viscosity, no tension | `AVZONE=2 AVTENSION=0` |
| step control | TR-BDF2 with its embedded error estimate (DTERR, required), respecting the mesh-motion bound (DT_BOUNDS) | `DTERR=0.1 NITER=25 NTRY=60` (TITAN's own rule, `CDYN`, `STOL`, is diagnostic only) |
| run control | | `TSTOP=2e8 NSTOP=1e27 DTSTOP=1e-9 CKEVERY=100 SNAPEVERY=50 PRINTEVERY=100 MSG=1 FLUSHOUT=1` |

Diagnostics and performance switches (DTDIAG, TIMING, CPLPROF, JACCHK, CPLJCHECK, SZJCHECK, NANSCAN, …) are listed in
docs/POPIII_SWITCH_AUDIT2.md §2. The options kept off by default are the RAY ionising ray, IonFS (branch), LyaTransfer (LYART)
and the GPU line kernels (LINEGPU).
