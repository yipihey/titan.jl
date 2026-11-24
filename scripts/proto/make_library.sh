#!/bin/bash
# make_library.sh : write the PROTOTYPE's own checkpoint library (replaces the zie1 replay checkpoints of states.txt). NOT run by the gate; run once the prototype
# branch (popiii-proto) exists:  CODE=<prototype worktree> scripts/proto/make_library.sh [--config F] [--cores "a-b c-d"]
#
# For each case (z100 = zie1's 100-zone mesh, z298 = zie3's 298-zone mesh) the prototype code is replayed FROM THE IC with the Stage-1 grey configuration and
# full checkpoints are dumped (CKDUMP) at the standard epochs:
#   z100: steps 115 (0 yr: protostar just formed), 295 (5 yr), 655 (111 yr), 2395 (1.76 kyr)   = the zie1 replay epochs (rp_k115/295/655/2395)
#   z298: the same EPOCHS mass-matched: the steps at which the M* of zie1's epochs is reached (0.017, 0.29, 3.2, 34 Msun); the step numbers of zie3 are not
#         those of zie1 (finer mesh -> more steps), so they are looked up in zie3's production track (runs/zie3/perf.csv: first row with mstar >= the mass).
# Resumable (as scripts/popiii_run/iceulrun.sh / subconv/replay.sh): a rerun restarts the newest _auto checkpoint, skips dumps that exist, and stops when the last
# epoch is dumped. A checkpoint of a step the run never reaches (e.g. a different dt history of the prototype) is written for the NEAREST step: the dump range is
# +-CKWIN steps (default 3) around the epoch, the file closest to the epoch is renamed ..._k<epoch>.* by `pick` at the end.
# Output: $LIB/<case>/ck/m2c_lib_<case>_k<step>.{aux,tstate.gz}, and the line block for states.txt printed to stdout (append it, delete the old lines).
# The gate itself needs: name, ckdir, ckname, epoch, case, nkz, nwdt -- nkz and nwdt are set from the table below (check them against the first NKCHECK run).
set -u
HERE=$(cd "$(dirname "$0")" && pwd); . $HERE/common.sh
CFG=$HERE/config_grey.env; LIB=${LIB:-$DATA/proto/library}; CKWIN=${CKWIN:-3}
while [ $# -gt 0 ]; do case $1 in --config) CFG=$2; shift 2;; --cores) GATE_CORES=$2; shift 2;; *) die "unknown arg $1";; esac; done
read -ra SLOTS <<< "$GATE_CORES"
EPOCHS_MSUN=(0.017 0.29 3.2 34); EPOCH_NAMES=(0yr 5yr 111yr 1.76kyr); STEPS_Z100=(115 295 655 2395)
NKZ_Z100=(20:30 33:43 34:44 1:10)                  # zones probed by NKCHECK at those epochs (docs M3b/M4); z298 mesh: scale ~3x, to be set from the first runs
NKZ_Z298=(60:90 100:130 100:135 1:30)
NWDT=(1.5e4 1e5 1e5 1e5)

z298_steps() {   # the first step of zie3's production track with M* >= each mass (perf.csv has mstar every PRINTEVERY steps; refine with CKWIN)
  python3 - "$DATA/runs/zie3/perf.csv" "${EPOCHS_MSUN[@]}" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
for m in map(float, sys.argv[2:]):
    r = next((x for x in rows if float(x["mstar"]) >= m), None)
    print(r["step"] if r else "NA")
PY
}

lib_case() {   # case cores steps...
  local cs=$1 cores=$2; shift 2; local steps=("$@") R=$LIB/$cs
  mkdir -p $R/ck; local log=$R/procedure.log; lg() { echo "$(date '+%F %T') $*" >> $log; }
  local last=${steps[-1]} dump="" s
  for s in "${steps[@]}"; do dump="$dump${dump:+,}$((s-CKWIN)):$((s+CKWIN))"; done
  while true; do
    ls $R/ck/m2c_lib_${cs}_k$last.aux > /dev/null 2>&1 && { lg "library complete"; break; }
    local auto=$(ls -t $R/ck/m2c_lib_${cs}_auto?.aux 2>/dev/null | head -1) rs=""
    [ -n "$auto" ] && rs="RESTART=$(basename $auto .aux | sed 's/^m2c_//')"
    lg "replay segment ${rs:-from the IC}, CKDUMP=$dump"
    run_pop $R $cores $WT $((last+CKWIN+1)) $CFG $cs $rs SAVE=lib_$cs CKDUMP=$dump CKEVERY=50 CKDIR=$R/ck
    local end=$(grep -E 'run ended|stopped' $R/run.out | tail -1 | cut -c1-200); lg "ended: $end"
    echo "$end" | grep -qE ':maxsteps|:tmax' || { lg "FAILED"; return 1; }
  done
  # pick: for each epoch the dumped file nearest to it is copied to the exact name _k<epoch>
  for s in "${steps[@]}"; do
    [ -f $R/ck/m2c_lib_${cs}_k$s.aux ] && continue
    local best=$(ls $R/ck | sed -nE "s/^m2c_lib_${cs}_k([0-9]+)\.aux$/\1/p" | awk -v t=$s 'BEGIN{b=-1} {d=$1-t; if(d<0)d=-d; if(b<0||d<bd){b=$1;bd=d}} END{print b}')
    [ "$best" -ge 0 ] 2>/dev/null && [ $((best-s)) -le $CKWIN ] && [ $((s-best)) -le $CKWIN ] && for e in aux tstate.gz; do cp $R/ck/m2c_lib_${cs}_k$best.$e $R/ck/m2c_lib_${cs}_k$s.$e; done && lg "epoch $s <- dump k$best"
  done
}

mapfile -t Z298 < <(z298_steps)
lib_case z100 ${SLOTS[0]} "${STEPS_Z100[@]}" &
lib_case z298 ${SLOTS[1]:-${SLOTS[0]}} "${Z298[@]}" &
wait
echo "# --- prototype library, $(date +%F), code $(git -C $WT rev-parse --short HEAD), config $(basename $CFG): replace the zie1rp lines of states.txt by"
for i in 0 1 2 3; do echo "z100_k${STEPS_Z100[$i]}  $LIB/z100/ck  lib_z100_k${STEPS_Z100[$i]}  ${EPOCH_NAMES[$i]}_(${EPOCHS_MSUN[$i]}Msun)  z100  ${NKZ_Z100[$i]}  ${NWDT[$i]}"; done
for i in 0 1 2 3; do echo "z298_k${Z298[$i]}  $LIB/z298/ck  lib_z298_k${Z298[$i]}  ${EPOCH_NAMES[$i]}_(${EPOCHS_MSUN[$i]}Msun)  z298  ${NKZ_Z298[$i]}  ${NWDT[$i]}"; done
