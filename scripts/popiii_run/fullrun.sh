#!/bin/bash
# From-scratch Pop III run with recording (docs/POPIII_RUNS.md). Every stage is its own julia process and saves under a fresh
# name (step numbers repeat across attempts); procedure.log records each start, its settings and how it ended.
#
#   A1  cloud collapse: hydrodynamics, chemistry, H2 line transfer, RA04 CIE continuum cooling (RAD=0), from the IC to
#       n_c = 3.5e18 cm^-3                                                              phaseA.env + stageA1.env
#   A2  radiation rebuild from the A1 checkpoint (RADFROM: E_r, F_r and the StellarOpacity continuum replace the CIE
#       escape-factor term) with static pre-refinement of the future core (REFINE); second collapse and core formation,
#       to T_c = 5e4 K                                                                  phaseA.env + stageA2.env
#   A3  core growth on the Lagrangian mesh, dynamical time-step limit relaxed, moving refinement window ahead of the
#       star, to M* = 0.5 Msun                                                          phaseA.env + stageA3.env
#   B   remap onto the fixed accretion mesh (ENVEUL), production flags, segments of SEG steps        phaseB.env
#
# RAD=1 must never be set from the IC without the rebuild: it switches the CIE cooling off and no radiation field exists.
#
# Usage: fullrun.sh <run> [<stage>:<checkpoint>] [threads]
#   A1:<ck> continue A1 from its checkpoint      R:<ck>  start A2 (rebuild) from a hydro-only checkpoint
#   A2:<ck> continue A2                          A3:<ck> start or continue A3        B:<ck> start or continue B
# A2X="<settings>" (environment): extra settings of the rebuild stage (they override stageA2.env)
# A1X="<settings>" (environment): the same for stage A1 (e.g. A1X="NSTOP=1e14" A2X="DTRADMAX=0": continuum transport from log n_c 14)
# <run directory>/run.env: settings of this run for every stage and every restart (e.g. MCLOUD=2e4); they override the stage files
HERE=$(cd "$(dirname "$0")" && pwd); JLENV=${JLENV:?set JLENV to a julia project with the repository packages}; DATA=${DATA:-/zpool/nvme/data/tabel_scratch/titan_popiii}
RUN=$1; FROM=${2:-}; TH=${3:-4}; SEG=100000; REFINE=${REFINE:-1e-3,0.1,0.03}
R=$DATA/runs/$RUN; CK=$R/ck; mkdir -p $CK
COMMON="RUNREC=$R POPIII_OUT=$R CKDIR=$CK SNAPEVERY=20 PRINTEVERY=100 SNAPS=0 TIMING=1 DTDIAG=1 MSG=1 TAG=_$RUN $([ -f $R/run.env ] && cat $R/run.env)"
cd $HERE/../../julia
log() { echo "$(date '+%F %T') $*" >> $R/procedure.log; }
run() { local out=$1; shift; env "$@" setsid julia -t $TH --project=$JLENV popiii_m2c.jl lines 200 > $out 2>&1 < /dev/null; }
fresh() { local b=$1 n=1; while [ -e $CK/m2c_${b}${n}.aux ] || [ -e $CK/m2c_${b}${n}_auto0.aux ]; do n=$((n+1)); done; echo ${b}${n}; }
# stage <name> <pattern that must end the log> <env...>: run one phase-A stage; sets ck to its final checkpoint
stage() {
  local name=$1 pat=$2; shift 2
  save=$(fresh z19${RUN}_${name}_); out=$R/stage${name}_$save.out
  log "stage $name start as $save: $A $* $COMMON"
  run $out $A "$@" $COMMON SAVE=$save MAXSTEPS=10000000
  if ! grep -q "$pat" $out; then log "stage $name FAILED ($save): $(grep -E '^stopped|run ended|^ERROR' $out | tail -2 | tr '\n' ' ')"; exit 1; fi
  log "stage $name done ($save): $(grep "$pat" $out | tail -1)"; ck=$save
}
A="$(cat $HERE/phaseA.env)"; st=A1; ck=
case "$FROM" in "") ;; A1:*|R:*|A2:*|A3:*|B:*) st=${FROM%%:*}; ck=${FROM#*:} ;; *) echo "unknown start $FROM" >&2; exit 2 ;; esac
if [ $st = A1 ]; then
  stage A1 'stopped: T_c' $(cat $HERE/stageA1.env) ${A1X:-} ${ck:+RESTART=$ck}; st=R
fi
if [ $st = R ]; then
  stage A2 'stopped: T_c' $(cat $HERE/stageA2.env) RADFROM=$ck REFINE=$REFINE ${A2X:-}
  grep -q 'RAD: rebuilt' $out || { log "stage A2 ($save) has no radiation rebuild"; exit 1; }; st=A3
elif [ $st = A2 ]; then
  stage A2 'stopped: T_c' $(cat $HERE/stageA2.env) RESTART=$ck; st=A3
fi
if [ $st = A3 ]; then
  stage A3 'stopped: M_star' $(cat $HERE/stageA3.env) RESTART=$ck
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
