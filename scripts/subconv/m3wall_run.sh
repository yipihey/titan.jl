#!/bin/bash
# m3wall run: <name> <MAXSTEPS> [env words]. Restart from ck/m2c_B2S1_auto0 (B2 step ~600), TDC=1 TDCPAR default of B2. Code = $CODE (default scratch wt).
set -u
DATA=/zpool/nvme/data/tabel_scratch/titan_popiii; T=$DATA/subconv/tests/m3wall; N=$1; MAX=$2; shift 2
WT=${CODE:-$T/wt}; R=$T/$N; mkdir -p $R/ck
[ -f $R/eqeos_he_xe.jls ] || cp $DATA/runs/zie1/eqeos_he_xe.jls $R/
cp -n $T/ck/m2c_B2S1_auto0.* $R/ck/
RUNENV="$(cat $DATA/runs/zie1/run.env | sed 's/SNAPEVERY=50/SNAPEVERY=5/g')"
cd $WT/julia
env RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_$N $RUNENV CONVDIAG=1 TDC=1 TDCPAR=1,1,1,1,0.01,1,0 RESTART=B2S1_auto0 SAVE=$N NWTRACE=$R/nwt.txt "$@" MAXSTEPS=$MAX \
  numactl --physcpubind=${PINC:-72-79} --localalloc julia -t ${TH:-8} --project=$DATA/subconv/jlenv popiii_m2c.jl lines 128 > $R/run.out 2>&1 < /dev/null
echo "done $N: rejects $(grep -c 'DTERR reject.*Newton' $R/run.out) $(tail -1 $R/run.out | cut -c1-160)"
