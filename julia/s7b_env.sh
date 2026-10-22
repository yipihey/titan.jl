#!/bin/zsh
# S7b envelope process (popiii_m2c.jl, ABN02 z = 19 sink run from the protostar-formation snapshot) with the S7b hooks. Usage: s7b_env.sh <tag> [KEY=VAL ...]
# The snapshot (core2) is not committed: SNAP=<path> (default lib/popiii/snaps_abn02z19_lines_n200_core2.jls).
cd "$(dirname "$0")"
tag=$1; shift
export TITAN_NOLIB=1 IC=abn02z19 TCMB=54.5 AVZONE=2 INTEGRATOR=trbdf2 CHEM=netatomicli NLTE=ep XFLOOR=1e-12 CHEMFALLBACK=1 SNAPS=0
export RESTART=${SNAP:-../lib/popiii/snaps_abn02z19_lines_n200_core2.jls} RESIDX=0 SINK=1 SPLIT=2 RSINKAU=5 RSINKP=1 RSMREF=1 RSMFLOOR=0.01 CDYN=0.03 TMAXYR=1e5 NSTOP=1e30 TSTOP=1e30 MAXSTEPS=600000
export S7B=1 S7B_TAG=$tag TAG=_s7b_$tag
for kv in "$@"; do export "$kv"; done
exec julia -t ${JT:-4} --project=env popiii_m2c.jl lines 200
