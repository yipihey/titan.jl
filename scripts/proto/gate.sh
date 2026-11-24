#!/bin/bash
# The prototype gate. docs/POPIII_PROTO_GATE.md.
#   gate.sh <T0|T1|T2|T3|all> [--ref metrics.json] [--out DIR] [--config FILE] [--bitref TREE|COMMIT] [--tol tol.json]
#                             [--states a,b] [--cores "120-127 112-119"] [--ref-export [case=]DIR] [--code WORKTREE]
# Writes <out>/<tier>/metrics.json + summary.txt (+ every run's files), then, with --ref, compare_metrics.py <tier>/metrics.json against --ref (a metrics.json of the same tier,
# or a directory holding <tier>/metrics.json) and exits non-zero on any FAIL. Without --ref only the absolute rules of tol.json (unit-test failures, NKCHECK bounds, ledger, bit-identity) are applied.
#   --config     run configuration (default scripts/proto/config_grey.env = Stage 1; config_tdc.env = Stage 2). Never hard-coded elsewhere.
#   --bitref     (T1) reference tree for the bit-identity check: a worktree directory or a commit; first use: prodpkg (the grey production default).
#   --out        default $DATA/proto/gate/<tier>_<UTC date>_<code commit>
#   env: GATE_CORES (slots of 8 cores, max 4 used), GATE_THREADS, GATE_STATES, GATE_PILOT_STEPS, GATE_ETOT_STEPS, GATE_SEG, GATE_MSTOP, GATE_T3_CASES, GATE_TIMEOUT, JLENV, CODE.
set -u
HERE=$(cd "$(dirname "$0")" && pwd); . $HERE/common.sh
declare -A REFEXPORT
TIER=${1:-}; [ -n "$TIER" ] || { sed -n 2,13p $0; exit 2; }; shift
REF=; OUT=; CFG=$HERE/config_grey.env; BITREF=; TOL=$HERE/tol.json
while [ $# -gt 0 ]; do case $1 in
  --ref) REF=$2; shift 2;; --out) OUT=$2; shift 2;; --config) CFG=$(cd "$(dirname "$2")" && pwd)/$(basename "$2"); shift 2;;
  --bitref) BITREF=$2; shift 2;; --tol) TOL=$2; shift 2;; --states) GATE_STATES=$2; export GATE_STATES; shift 2;;
  --cores) GATE_CORES=$2; shift 2;; --code) WT=$(cd "$2" && pwd); shift 2;;
  --ref-export) case $2 in *=*) REFEXPORT[${2%%=*}]=${2#*=};; *) REFEXPORT[z100]=$2; REFEXPORT[z298]=$2;; esac; shift 2;;
  *) die "unknown argument $1";; esac; done
export CFG BITREF
case ${TIER^^} in T0|T1|T2|T3|ALL) ;; *) die "tier must be T0, T1, T2, T3 or all";; esac
TIERS=${TIER^^}; [ $TIERS = ALL ] && TIERS="T0 T1 T2"      # 'all' = the three tiers that fit a working day; T3 is the nightly one
[ -n "$OUT" ] || OUT=$DATA/proto/gate/$(date -u +%Y%m%d)_$(git -C $WT rev-parse --short HEAD)
mkdir -p $OUT; OUT=$(cd $OUT && pwd); export OUT WT
read -ra SLOTS <<< "$GATE_CORES"
echo "gate: tiers $TIERS code $WT ($(git -C $WT rev-parse --short HEAD)) config $CFG out $OUT cores [$GATE_CORES] $(stamp)" | tee -a $OUT/gate.log
for f in tier_t0 tier_t1 tier_t2 tier_t3; do . $HERE/$f.sh; done
rc=0
for t in $TIERS; do
  echo "gate: $t start $(stamp)" | tee -a $OUT/gate.log
  case $t in T0) tier_t0;; T1) tier_t1;; T2) tier_t2;; T3) tier_t3;; esac
  d=$OUT/${t,,}; echo "gate: $t done $(stamp): $d/metrics.json" | tee -a $OUT/gate.log
  r=-; if [ -n "$REF" ]; then if [ -d "$REF" ]; then r=$REF/${t,,}/metrics.json; else r=$REF; fi; fi
  python3 $HERE/compare_metrics.py $d/metrics.json $r $TOL > $d/compare.txt 2>&1; c=$?
  grep -vE "^(PASS|INFO)" $d/compare.txt | tail -40; tail -1 $d/compare.txt
  [ $c -ne 0 ] && rc=1
done
echo "gate: finished $(stamp), exit $rc" | tee -a $OUT/gate.log
exit $rc
