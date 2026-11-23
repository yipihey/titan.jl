#!/bin/bash
# M3b test runs: the zie1 settings, MAXSTEPS steps, SNAPEVERY=5, from the IC or from a replay checkpoint (RESTART=rp_kNNN: its files are COPIED from zie1rp/ck into the run's own ck dir).
# Usage: m3brun.sh <name> <MAXSTEPS> [extra env words ...]   (e.g. RESTART=rp_k650 TDC=1 CONVDIAG=1 NKRYLOV=1 NWPROBE=1.0e5 NWPROBEAT=0 NKCHECK=1)
# Env: CODE = the worktree whose julia/ is run (default: the one this script lives in; a scratch worktree with TitanCore.jl of aace2a0 for the restart tests, see docs M3b).
# PINC = cores (default 64-71), TH = threads (8). Output in $T/m3b/<name>.
set -u
DATA=/zpool/nvme/data/tabel_scratch/titan_popiii; T=$DATA/subconv/tests/m3b; N=$1; MAX=$2; shift 2
WT=${CODE:-$(cd "$(dirname "$0")/../.." && pwd)}; R=$T/$N; mkdir -p $R/ck
[ -f $R/eqeos_he_xe.jls ] || cp $DATA/runs/zie1/eqeos_he_xe.jls $R/
for w in "$@"; do case $w in RESTART=rp_k*) for e in aux tstate.gz; do [ -f $R/ck/m2c_${w#RESTART=}.$e ] || cp $DATA/subconv/runs/zie1rp/ck/m2c_${w#RESTART=}.$e $R/ck/; done;; esac; done
RUNENV="$(cat $DATA/runs/zie1/run.env | sed 's/SNAPEVERY=50/SNAPEVERY=5/g')"
cd $WT/julia
env RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_$N $RUNENV "$@" MAXSTEPS=$MAX SNAPEVERY=5 \
  numactl --physcpubind=${PINC:-64-71} --localalloc julia -t ${TH:-8} --project=$DATA/subconv/jlenv popiii_m2c.jl lines 128 > $R/run.out 2>&1 < /dev/null
echo "done $N: $(tail -1 $R/run.out | cut -c1-200)"
