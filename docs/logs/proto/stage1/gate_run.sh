#!/bin/bash
# Stage-1 gate runner (docs/POPIII_PROTO.md): <steps> steps of the zie1 or zie3 configuration from the IC (or a restart), SNAPEVERY=5.
# Usage: gate_run.sh <label> <tree> <zie1|zie3> <steps> <ref|proto> [extra VAR=val ...]
#   ref   : the production launch of iceulrun.sh: phaseB.env + run.env copies in $REFENV/<run> (proto/refenv, copied from the production run directories on 24 Nov before they switched to prodpkg; launch.env NZ/TH/SEG likewise), in <tree> (a pkgstart checkout)
#   proto : the Stage-1 defaults (julia/popiii_stack.jl) + scripts/proto/run_<run>.env only
# Output: $OUT/<label>/ (run.out, snaps_*, ck/). Env: CKSRC (checkpoint dir for RESTART=name, default prod_ck/<run>), EQEOS, REFENV, PINC (cores, default 112-119), TH (default 8), JLENV, OUT.
set -u
L=$1; TREE=$2; RUN=$3; N=$4; MODE=$5; shift 5
P=/zpool/nvme/data/tabel_scratch/titan_popiii; OUT=${OUT:-$P/proto/stage1}; R=$OUT/$L; mkdir -p $R/ck
JLENV=${JLENV:-$P/proto/jlenv}; RR=${REFENV:-$P/proto/refenv}/$RUN
[ -f $R/eqeos_he_xe.jls ] || cp ${EQEOS:-$P/prod_ck/$RUN/eqeos_he_xe.jls} $R/
for w in "$@"; do case $w in RESTART=*) cp -n ${CKSRC:-$P/prod_ck/$RUN}/m2c_${w#RESTART=}.* $R/ck/;; esac; done   # restart checkpoint files are COPIED
. $RR/launch.env     # NZ, TH, PIN, SEG of the production run; PIN is replaced below
case $MODE in
  ref)   B="$(grep -v '^#' $RR/phaseB.env)"; RUNENV="$(cat $RR/run.env)";;
  proto) B=""; RUNENV="$(grep -v '^#' $(dirname "$0")/run_$RUN.env | tr '\n' ' ')";;
esac
cd $TREE/julia
T0=$(date +%s.%N)
env $B RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_$RUN $RUNENV "$@" MAXSTEPS=$N SNAPEVERY=5 \
  numactl --physcpubind=${PINC:-112-119} --localalloc julia -t ${TH:-8} --project=$JLENV popiii_m2c.jl lines $NZ > $R/run.out 2>&1 < /dev/null
echo "===== EXIT $? wall $(echo "$(date +%s.%N) - $T0" | bc) =====" >> $R/run.out
