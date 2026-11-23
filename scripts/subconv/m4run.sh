#!/bin/bash
# M4 bracket run: one case from a replay checkpoint, in segments of SEG steps (each restarts the previous segment's final checkpoint),
# until M* >= MSTOP or NSEGMAX segments. Resumable: rerunning continues after the last completed segment (procedure.log).
# Code: CODE (scratch worktree of the branch commit with julia/TitanCore.jl of aace2a0, as for all replay-checkpoint restarts; see docs M3b).
# Usage: m4run.sh <case> <start rp_kNNN> <cores> [env words ...]
set -u
DATA=/zpool/nvme/data/tabel_scratch/titan_popiii; T=$DATA/subconv/m4; C=$1; START=$2; PINC=$3; shift 3
SEG=${SEG:-400}; MSTOP=${MSTOP:-34}; NSEGMAX=${NSEGMAX:-40}; WT=${CODE:?set CODE}
R=$T/$C; mkdir -p $R/ck; log() { echo "$(date '+%F %T') $*" >> $R/procedure.log; }
[ -f $R/eqeos_he_xe.jls ] || cp $DATA/runs/zie1/eqeos_he_xe.jls $R/
for e in aux tstate.gz; do [ -f $R/ck/m2c_$START.$e ] || cp $DATA/subconv/runs/zie1rp/ck/m2c_$START.$e $R/ck/; done
RUNENV="$(cat $DATA/runs/zie1/run.env)"
n=$(grep -c 'segment done' $R/procedure.log 2>/dev/null); ck=$START; [ "$n" -gt 0 ] && ck=$(grep 'segment done' $R/procedure.log | tail -1 | awk '{print $5}')
log "lane start: case $C from $ck, code $(git -C $WT rev-parse --short HEAD) (+TitanCore aace2a0), cores $PINC, env: $*"
cd $WT/julia
while [ $n -lt $NSEGMAX ]; do
  n=$((n+1)); save=${C}S$n; out=$R/seg_$save.out
  env RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_$C $RUNENV CONVDIAG=1 SNAPEVERY=5 "$@" \
      RESTART=$ck SAVE=$save MAXSTEPS=$SEG setsid numactl --physcpubind=$PINC --localalloc julia -t 8 --project=$DATA/subconv/jlenv popiii_m2c.jl lines 128 > $out 2>&1 < /dev/null
  end=$(grep -E 'run ended|stopped' $out | tail -1 | cut -c1-200)
  ms=$(grep -E '^ *hot: t' $out | tail -1 | sed -E 's/.*m\* ([0-9.eE+-]+) Msun.*/\1/')
  if ! echo "$end" | grep -q ':maxsteps'; then log "segment FAILED $save: $end (M* $ms)"; exit 1; fi
  log "segment done $save : $end (M* $ms)"; ck=$save
  python3 -c "import sys; sys.exit(0 if float('$ms')>=$MSTOP else 1)" 2>/dev/null && { log "reached M* $ms >= $MSTOP"; break; }
done
