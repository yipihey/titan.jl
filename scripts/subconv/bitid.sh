#!/bin/bash
# M1 bit-identity test: the zie1 settings from the IC (or a restart), MAXSTEPS steps with SNAPEVERY=5, CONVDIAG unset or =1.
# Usage: bitid.sh <name> <MAXSTEPS> [extra env words ...]     (extra e.g. CONVDIAG=1  or  RESTART=rp_k650 SAVE=x)
# Output in $T/<name>; code from the worktree this script lives in. Pinned to 8 cores (80-87), 8 threads.
set -u
DATA=/zpool/nvme/data/tabel_scratch/titan_popiii; T=$DATA/subconv/tests; N=$1; MAX=$2; shift 2
WT=$(cd "$(dirname "$0")/../.." && pwd); R=$T/$N; mkdir -p $R/ck
[ -f $R/eqeos_he_xe.jls ] || cp $DATA/runs/zie1/eqeos_he_xe.jls $R/
RUNENV="$(cat $DATA/runs/zie1/run.env | sed 's/SNAPEVERY=50/SNAPEVERY=5/g')"
cd $WT/julia
env RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_$N $RUNENV "$@" MAXSTEPS=$MAX SNAPEVERY=5 \
  numactl --physcpubind=80-87 --localalloc julia -t 8 --project=$DATA/subconv/jlenv popiii_m2c.jl lines 128 > $R/run.out 2>&1 < /dev/null
echo "done $N: $(tail -1 $R/run.out | cut -c1-200)"
