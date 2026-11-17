#!/bin/bash
# The Pop III run launcher (docs/POPIII_STACK.md): one run from the initial condition on the fixed log mesh (GRIDEQ moving mesh), the
# physics stack of julia/popiii_stack.jl, in segments of SEG steps (each segment restarts the previous one's final checkpoint with
# the same settings; the first one starts from the IC).
#
# Usage: iceulrun.sh <run> [<checkpoint to continue>] [threads]
# <run directory>/phaseB.env (optional; default: phaseB.env next to this script, which sets nothing): settings that differ from the
#   stack; <run directory>/run.env: the run's own settings (mesh ENVEUL/ENVCORE, MCLOUD, ...), appended after it (later entries win);
#   <run directory>/launch.env: NZ (IC shells for the tabulated profile set-up, default 128), TH, PIN, SEG
HERE=$(cd "$(dirname "$0")" && pwd); JLENV=${JLENV:?set JLENV to a julia project with the repository packages}; DATA=${DATA:-/zpool/nvme/data/tabel_scratch/titan_popiii}
RUN=$1; ck=${2:-}; R=$DATA/runs/$RUN; CK=$R/ck; mkdir -p $CK
[ -f $R/launch.env ] && . $R/launch.env
TH=${3:-${TH:-8}}; SEG=${SEG:-1000}; NZ=${NZ:-128}; PIN=${PIN:-}
sf() { if [ -f $R/$1 ]; then cat $R/$1; else cat $HERE/$1; fi; }
COMMON="RUNREC=$R POPIII_OUT=$R CKDIR=$CK EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_$RUN $([ -f $R/run.env ] && cat $R/run.env)"
cd $HERE/../../julia
log() { echo "$(date '+%F %T') $*" >> $R/procedure.log; }
run() { local out=$1; shift; env "$@" setsid $PIN julia -t $TH --project=$JLENV popiii_m2c.jl lines $NZ > $out 2>&1 < /dev/null; }
fresh() { local b=$1 n=1; while [ -e $CK/m2c_${b}${n}.aux ] || [ -e $CK/m2c_${b}${n}_auto0.aux ]; do n=$((n+1)); done; echo ${b}${n}; }
B="$(sf phaseB.env | grep -v '^#')"
log "launch (iceulrun.sh): $NZ IC shells, $TH threads${PIN:+, $PIN}, commit $(git -C $HERE rev-parse --short HEAD 2>/dev/null)${ck:+, from $ck}"
while true; do
  save=$(fresh z19${RUN}S); out=$R/seg_$save.out
  if [ -z "$ck" ]; then start="(from the IC)"; st=; else start="RESTART=$ck"; st=$start; fi
  log "segment start ($start) as $save: $B $COMMON"
  run $out $B $COMMON $st SAVE=$save MAXSTEPS=$SEG
  end=$(grep -E 'run ended|stopped' $out | tail -1)
  if ! echo "$end" | grep -q ':maxsteps'; then log "segment FAILED ($save): $(grep -E 'stopped|run ended' $out | tail -2 | tr '\n' ' ')"; exit 1; fi
  log "segment done ($save): $end"; ck=$save
done
