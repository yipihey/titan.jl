#!/bin/bash
# M3c debug run: restart from the k115 TDC checkpoint dbg_k153 (step 153), MAXSTEPS steps. Code = scratch worktree $CODE (TitanCore.jl of aace2a0; julia/*.jl else synced from the branch by m3dbgsync.sh).
# Usage: m3dbgrun.sh <name> <MAXSTEPS> [env words ...]   Output $T/m3dbg/<name>. PINC cores (default 64-71), TH threads (8).
set -u
DATA=/zpool/nvme/data/tabel_scratch/titan_popiii; T=$DATA/subconv/tests/m3dbg; N=$1; MAX=$2; shift 2
WT=${CODE:-$T/wt}; R=$T/$N; mkdir -p $R/ck
[ -f $R/eqeos_he_xe.jls ] || cp $DATA/runs/zie1/eqeos_he_xe.jls $R/
for e in aux tstate.gz; do [ -f $R/ck/m2c_${CK:-dbg_k153}.$e ] || cp $T/ck/m2c_${CK:-dbg_k153}.$e $R/ck/; done
RUNENV="$(cat $DATA/runs/zie1/run.env | sed 's/SNAPEVERY=50/SNAPEVERY=5/g')"
cd $WT/julia
env RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_$N $RUNENV TDC=1 RESTART=${CK:-dbg_k153} NWTRACE=$R/nwt.txt "$@" MAXSTEPS=$MAX SNAPEVERY=5 \
  numactl --physcpubind=${PINC:-64-71} --localalloc julia -t ${TH:-8} --project=$DATA/subconv/jlenv popiii_m2c.jl lines 128 > $R/run.out 2>&1 < /dev/null
echo "done $N: $(tail -1 $R/run.out | cut -c1-200)"
