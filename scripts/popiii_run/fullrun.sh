#!/bin/bash
# From-scratch Pop III run with recording (docs/POPIII_RUNS.md). Phase A: cloud collapse on the Lagrangian mesh with rezoning,
# from the IC (or from checkpoint $2) to n_c = NSTOPA. Phase B: restart remapped onto the fixed accretion mesh (ENVEUL) with
# the production flags, in segments of SEG steps. Every attempt saves under a fresh name (step numbers repeat across
# attempts); procedure.log records each start, its settings and how it ended.
# Usage: fullrun.sh <run> [restart checkpoint for phase A | "B:<checkpoint>" to resume phase B] [threads]
HERE=$(cd "$(dirname "$0")" && pwd); JLENV=${JLENV:?set JLENV to a julia project with the repository packages}; DATA=${DATA:-/zpool/nvme/data/tabel_scratch/titan_popiii}
RUN=$1; FROM=${2:-}; TH=${3:-4}; SEG=100000; NSTOPA=2.55e23
R=$DATA/runs/$RUN; CK=$R/ck; mkdir -p $CK
COMMON="RUNREC=$R POPIII_OUT=$R CKDIR=$CK SNAPEVERY=20 PRINTEVERY=100 SNAPS=0 TIMING=1 DTDIAG=1 MSG=1 TAG=_$RUN"
cd $HERE/../../julia
log() { echo "$(date '+%F %T') $*" >> $R/procedure.log; }
run() { local out=$1; shift; env "$@" setsid julia -t $TH --project=$JLENV popiii_m2c.jl lines 200 > $out 2>&1 < /dev/null; }
fresh() { local b=$1 n=1; while [ -e $CK/m2c_${b}${n}.aux ] || [ -e $CK/m2c_${b}${n}_auto0.aux ]; do n=$((n+1)); done; echo ${b}${n}; }
if [[ "$FROM" != B:* ]]; then
  save=$(fresh z19${RUN}A); out=$R/phaseA_$save.out
  log "phase A start ${FROM:+from $FROM }as $save: $(cat $HERE/phaseA.env) $COMMON"
  run $out $(cat $HERE/phaseA.env) $COMMON NSTOP=$NSTOPA SAVE=$save MAXSTEPS=10000000 ${FROM:+RESTART=$FROM}
  if ! grep -q 'stopped: T_c' $out; then log "phase A FAILED ($save): $(grep -E 'stopped|run ended' $out | tail -2 | tr '\n' ' ')"; exit 1; fi
  log "phase A done ($save): $(grep 'stopped: T_c' $out)"; ck=$save
else
  ck=${FROM#B:}
fi
B="$(cat $HERE/phaseB.env)"
while true; do
  save=$(fresh z19${RUN}B); out=$R/phaseB_$save.out
  ms0=$(grep -h -oE 'M_shell0 [0-9.e+-]+' $R/m2c_*_${RUN}.log 2>/dev/null | tail -1 | awk '{print $2}')
  log "phase B start from $ck (ENVMS0 ${ms0:-from remap}) as $save: $B $COMMON"
  run $out $B $COMMON ${ms0:+ENVMS0=$ms0} RESTART=$ck SAVE=$save MAXSTEPS=$SEG
  end=$(grep -E 'run ended|stopped' $out | tail -1)
  if ! echo "$end" | grep -q ':maxsteps'; then log "phase B FAILED ($save): $(grep -E 'stopped|run ended' $out | tail -2 | tr '\n' ' ')"; exit 1; fi
  log "phase B segment done ($save): $end"; ck=$save
done
