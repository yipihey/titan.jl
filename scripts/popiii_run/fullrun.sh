#!/bin/bash
# From-scratch Pop III run with recording. Phase A: cloud collapse on the Lagrangian mesh with rezoning, from the IC to
# n_c = NSTOPA (the 0.47 M_sun protostar; the old lineage switched at log n_c 23.41). Phase B: restart remapped onto the fixed
# accretion mesh (ENVEUL) with the production flags, in segments of SEG steps, each from the newest checkpoint.
# Usage: fullrun.sh <run> [threads]
HERE=$(cd "$(dirname "$0")" && pwd); JLENV=${JLENV:?set JLENV to a julia project with the repository packages}; DATA=${DATA:-/zpool/nvme/data/tabel_scratch/titan_popiii}
RUN=$1; TH=${2:-4}; SEG=100000; NSTOPA=2.55e23
R=$DATA/runs/$RUN; CK=$R/ck; mkdir -p $CK
COMMON="RUNREC=$R POPIII_OUT=$R CKDIR=$CK SNAPEVERY=20 PRINTEVERY=100 SNAPS=0 TIMING=1 DTDIAG=1 TAG=_$RUN"
cd $HERE/../../julia
log() { echo "$(date '+%F %T') $*" >> $R/procedure.log; }
run() { local out=$1; shift; env "$@" setsid julia -t $TH --project=$JLENV popiii_m2c.jl lines 200 > $out 2>&1 < /dev/null; }
if [ ! -f $CK/m2c_z19${RUN}A.aux ]; then
  log "phase A start: IC -> n_c $NSTOPA; env: $(cat $HERE/phaseA.env) $COMMON"
  run $R/phaseA.out $(cat $HERE/phaseA.env) $COMMON NSTOP=$NSTOPA SAVE=z19${RUN}A MAXSTEPS=10000000
  grep -q 'stopped: T_c .*NSTOP' $R/phaseA.out || { log "phase A did not reach NSTOP: $(grep -E 'stopped|ERROR|run ended' $R/phaseA.out | tail -2)"; exit 1; }
  log "phase A done: $(grep 'stopped: T_c' $R/phaseA.out)"
fi
B="$(cat $HERE/phaseB.env)"
ck=z19${RUN}A; seg=1
while true; do
  save=z19${RUN}B$seg
  ms0=$(grep -h -oE 'M_shell0 [0-9.e+-]+' $R/m2c_*_${RUN}.log 2>/dev/null | tail -1 | awk '{print $2}')
  log "phase B segment $seg from $ck (ENVMS0 ${ms0:-from remap}) as $save; env: $B $COMMON"
  run $R/phaseB$seg.out $B $COMMON ${ms0:+ENVMS0=$ms0} RESTART=$ck SAVE=$save MAXSTEPS=$SEG
  end=$(grep -E 'run ended|stopped' $R/phaseB$seg.out | tail -1)
  log "segment $seg ended: $end"
  echo "$end" | grep -q ':maxsteps' || exit 1
  ck=$save; seg=$((seg+1))
done
