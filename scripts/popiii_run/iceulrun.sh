#!/bin/bash
# One-stage Pop III run from the initial condition on the fixed log-Eulerian mesh (ICEUL=1, popiii_m2c.jl): no Lagrangian stages,
# no remap, RAD=1 and the phase-B physics stack from t = 0, in segments of SEG steps (each segment restarts the previous one's
# final checkpoint with the same settings; the first one builds the IC with ICEUL=1).
#
# Usage: iceulrun.sh <run> [<checkpoint to continue>] [threads]
# <run directory>/phaseB.env: the physics stack (default: phaseB.env next to this script); <run directory>/run.env: the run's own
#   settings (mesh ENVEUL/ENVCORE, MCLOUD, ...), appended after it; <run directory>/launch.env: NZ (IC shells for the tabulated
#   profile set-up, default 200), TH, PIN, SEG
# iceul.env next to this script: the run.env of z100ic3 (14 Nov), the first run from the IC past the protostar's swelling to the main
#   sequence (100 zones; t_c128c's run.env without ENVMS0, + MCLOUD ENVCONS, the mesh, and the from-IC levers; later entries override
#   earlier ones, e.g. CPLBAND= and ECHEMSPLIT=0 for COUPLED in every network zone). Use with t_c128c's phaseB.env as <run>/phaseB.env.
HERE=$(cd "$(dirname "$0")" && pwd); JLENV=${JLENV:?set JLENV to a julia project with the repository packages}; DATA=${DATA:-/zpool/nvme/data/tabel_scratch/titan_popiii}
RUN=$1; ck=${2:-}; R=$DATA/runs/$RUN; CK=$R/ck; mkdir -p $CK
[ -f $R/launch.env ] && . $R/launch.env
TH=${3:-${TH:-4}}; SEG=${SEG:-100000}; NZ=${NZ:-200}; PIN=${PIN:-}
sf() { if [ -f $R/$1 ]; then cat $R/$1; else cat $HERE/$1; fi; }
COMMON="RUNREC=$R POPIII_OUT=$R CKDIR=$CK SNAPEVERY=20 PRINTEVERY=100 SNAPS=0 TIMING=1 DTDIAG=1 MSG=1 TAG=_$RUN $([ -f $R/run.env ] && cat $R/run.env)"
cd $HERE/../../julia
log() { echo "$(date '+%F %T') $*" >> $R/procedure.log; }
run() { local out=$1; shift; env "$@" setsid $PIN julia -t $TH --project=$JLENV popiii_m2c.jl lines $NZ > $out 2>&1 < /dev/null; }
fresh() { local b=$1 n=1; while [ -e $CK/m2c_${b}${n}.aux ] || [ -e $CK/m2c_${b}${n}_auto0.aux ]; do n=$((n+1)); done; echo ${b}${n}; }
B="$(sf phaseB.env)"
log "launch (iceulrun.sh): $NZ IC shells, $TH threads${PIN:+, $PIN}, commit $(git -C $HERE rev-parse --short HEAD 2>/dev/null)${ck:+, from $ck}"
while true; do
  save=$(fresh z19${RUN}S); out=$R/seg_$save.out
  if [ -z "$ck" ]; then start="ICEUL=1"; else start="RESTART=$ck"; fi
  log "segment start ($start) as $save: $B $COMMON"
  run $out $B $COMMON $start SAVE=$save MAXSTEPS=$SEG
  end=$(grep -E 'run ended|stopped' $out | tail -1)
  if ! echo "$end" | grep -q ':maxsteps'; then log "segment FAILED ($save): $(grep -E 'stopped|run ended' $out | tail -2 | tr '\n' ' ')"; exit 1; fi
  log "segment done ($save): $end"; ck=$save
done
