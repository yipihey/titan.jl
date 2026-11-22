#!/bin/bash
# M1 diagnosis: restart from a replay checkpoint with CONVDIAG=1, run settings of the production run, code from the PINNED worktree wt_diag.
# Usage: diag_epoch.sh <prodrun> <replayrun> <ckstep> <MAXSTEPS> [nz]     e.g. diag_epoch.sh zie1 zie1rp 115 200 128
set -u
DATA=/zpool/nvme/data/tabel_scratch/titan_popiii; S=$DATA/subconv; P=$DATA/runs/$1; Q=$S/runs/$2; K=$3; MAX=$4; NZ=${5:-128}
WT=${WT:-$S/wt_diag_b}; CORES=${CORES:-80-87}; R=$S/tests/m1diag/$1_k$K; mkdir -p $R/ck
cp $Q/ck/m2c_rp_k$K.aux $Q/ck/m2c_rp_k$K.tstate.gz $R/ck/ 2>/dev/null; cp -n $Q/eqeos_he_xe.jls $R/
cd $WT/julia
env RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_d$1_k$K $(cat $P/run.env | sed 's/SNAPEVERY=50/SNAPEVERY=5/g') \
  CONVDIAG=1 RESTART=rp_k$K SAVE=x MAXSTEPS=$MAX SNAPEVERY=5 \
  numactl --physcpubind=$CORES --localalloc julia -t 8 --project=$S/jlenv popiii_m2c.jl lines $NZ > $R/run.out 2>&1 < /dev/null
echo "done $1 k$K: $(tail -1 $R/run.out | cut -c1-200)" > $R/done
