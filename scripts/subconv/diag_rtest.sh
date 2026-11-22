#!/bin/bash
# restart-reproducibility test: diag_rtest.sh <name> <worktree> <prodrun> <replayrun> <ckstep> <MAXSTEPS> <cores> [extra env words]   (CONVDIAG only if given in the extra words)
set -u
DATA=/zpool/nvme/data/tabel_scratch/titan_popiii; S=$DATA/subconv; N=$1; WT=$2; P=$DATA/runs/$3; Q=$S/runs/$4; K=$5; MAX=$6; CORES=$7; shift 7
R=$S/tests/m1diag/rt_$N; mkdir -p $R/ck; cp $Q/ck/m2c_rp_k$K.aux $Q/ck/m2c_rp_k$K.tstate.gz $R/ck/; cp -n $Q/eqeos_he_xe.jls $R/
cd $WT/julia
env RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_rt$N $(cat $P/run.env | sed 's/SNAPEVERY=50/SNAPEVERY=5/g') "$@" \
  RESTART=rp_k$K SAVE=x MAXSTEPS=$MAX SNAPEVERY=5 numactl --physcpubind=$CORES --localalloc julia -t 8 --project=$S/jlenv popiii_m2c.jl lines 128 > $R/run.out 2>&1 < /dev/null
echo "done $N" > $R/done
