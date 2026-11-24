#!/bin/bash
# M6 copy of mgint_run.sh with CORES=<cpulist> (default 112-119); M3e (MGROUPS x convection): one run of the zie1 stack (docs/POPIII_CONVECTION.md M3e). Cores ${CORES:-112-119}, 8 threads, one job at a time.
# Usage: mgint_run.sh <name> <worktree> <MAXSTEPS> [env words ...]    e.g. RESTART=rp_k650 CKSRC=<dir with m2c_<name>.{aux,tstate.gz}> MGROUPS=4 TDC=1 CONVDIAG=1
# Output in $BITT (default .../subconv/tests/mgint)/<name>; code from <worktree>. IC run unless RESTART=<checkpoint name> (the checkpoint files are COPIED from CKSRC, never written there).
set -u
DATA=/zpool/nvme/data/tabel_scratch/titan_popiii; T=${BITT:-$DATA/subconv/tests/mgint}; N=$1; WT=$2; MAX=$3; shift 3
R=$T/$N; mkdir -p $R/ck
[ -f $R/eqeos_he_xe.jls ] || cp $DATA/runs/zie1/eqeos_he_xe.jls $R/
RUNENV="$(cat $DATA/runs/zie1/run.env | sed 's/SNAPEVERY=50/SNAPEVERY=5/g')"
EW=(); CKSRC=
for w in "$@"; do case "$w" in CKSRC=*) CKSRC=${w#CKSRC=};; *) EW+=("$w");; esac; done
RN=; for w in "${EW[@]}"; do case "$w" in RESTART=*) RN=${w#RESTART=};; esac; done
if [ -n "$RN" ]; then for f in aux tstate.gz; do [ -f $R/ck/m2c_$RN.$f ] || cp $CKSRC/m2c_$RN.$f $R/ck/; done; fi
cd $WT/julia
echo "===== START $(date) $N wt=$(git -C $WT rev-parse --short HEAD) $*" > $R/run.out
env RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_$N $RUNENV "${EW[@]}" MAXSTEPS=$MAX SNAPEVERY=5 \
  numactl --physcpubind=${CORES:-112-119} --localalloc julia -t 8 --project=$DATA/subconv/jlenv popiii_m2c.jl lines 128 >> $R/run.out 2>&1 < /dev/null
echo "done $N: $(tail -1 $R/run.out | cut -c1-200)"
