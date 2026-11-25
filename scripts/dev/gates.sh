#!/bin/bash
# gates.sh <worktree> [--nodes 3,4] [--only 1,2,3,4,5] [--ref COMMIT|DIR] [--inproc] [--include-path] [--nocache] [--out DIR] [--mode proto|ref]
# The five default merge gates of docs/POPIII_PROTO_GATE.md for <worktree>, against the reference commit (default: the first parent of the
# tree's HEAD = the previous gated commit on popiii-proto). Extends scripts/proto/gate_run.sh (the runner) and the P5 harness (prod_env/cmp.sh,
# metrics.py); needs scripts/dev/mkimage.sh (TitanStack image per commit) and the checkpoint library $DATA/cklib (README there).
#   gate 1  NKCHECK=2 grey and MGROUPS=4, zie1 S6_auto1 (5 steps)        refs 1.5e-10 / 2.306e-7 (doc constants)
#   gate 2  bit-identity IC 40 steps, MGROUPS off, vs reference commit    (step lines + decompressed tstate)
#   gate 3  late restarts zie1 S16 and zie3 S8, 50 steps                  (completes; identical to the reference commit's run)
#   gate 4  ETOTMON ledger over those two 50-step windows                 (no worse than the reference commit)
#   gate 5  MG window zie3 fp3a1, 40 steps, MGROUPS=4                     (all accepted, dt median ~9e3 s)
# Runs on one julia process per NUMA node of --nodes (default 3,4; node n = cpus 8n..8n+7, 8 threads). Jobs are assigned to nodes
# statically; results of the REFERENCE commit are cached in $DATA/dev/gcache/<tree>/ and reused (--nocache: everything in a fresh dir).
#   --inproc        a node's jobs of one tree and image variant run in ONE julia process (scripts/dev/inproc.jl); default: one process per job
#   --include-path  no TITAN_PKG image: the include path (the old speed; for timing comparisons)
#   --mode          proto (default): popiii_stack defaults + scripts/proto/run_<run>.env; ref: the production run.env copies in $DATA/cklib
# Prints one summary table (gate, value, reference, PASS/FAIL), exit status 1 if any gate FAILs. Logs: <out>/gates.log.
set -u
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd $HERE/../.. && pwd)
DATA=${DATA:-/zpool/nvme/data/tabel_scratch/titan_popiii}; CK=$DATA/cklib
WT=; NODES=3,4; ONLY=1,2,3,4,5; REF=; INPROC=0; INCL=0; NOCACHE=0; OUTARG=; MODE=proto
while [ $# -gt 0 ]; do case $1 in
  --nodes) NODES=$2; shift 2;; --only) ONLY=$2; shift 2;; --ref) REF=$2; shift 2;; --inproc) INPROC=1; shift;; --include-path) INCL=1; shift;;
  --nocache) NOCACHE=1; shift;; --out) OUTARG=$2; shift 2;; --mode) MODE=$2; shift 2;; -h|--help) sed -n 2,19p $0; exit 0;;
  -*) echo "gates: unknown option $1" >&2; exit 2;; *) WT=$1; shift;; esac; done
[ -n "$WT" ] || { sed -n 2,4p $0; exit 2; }
[ $INPROC = 1 ] && [ "${GATES_ALLOW_INPROC:-0}" != 1 ] && { echo "gates: --inproc is disabled: several PopIIIRun.main calls in one process are NOT byte-identical to separate processes (docs/POPIII_DEVTOOLS.md, driver globals leak). GATES_ALLOW_INPROC=1 overrides, for experiments only." >&2; exit 2; }
WT=$(cd $WT && pwd); IFS=, read -ra NODE <<< "$NODES"; NN=${#NODE[@]}
has() { case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }
TNAME=$($HERE/mkimage.sh $WT name) || exit 2
# ---- reference tree: a git-archive scratch copy (real files) of the commit, unless a directory is given
NEEDREF=0; for g in 2 3 4; do has $g && NEEDREF=1; done
RT=; RNAME=-
if [ $NEEDREF = 1 ]; then
  if [ -d "$REF/julia" ]; then RT=$(cd $REF && pwd)
  else
    TSHA=$(cat $WT/.titan_sha 2>/dev/null || git -C $WT rev-parse HEAD); RSHA=$(git -C $REPO rev-parse --short ${REF:-$TSHA^1}) || exit 2
    RT=$DATA/dev/trees/$RSHA
    if [ ! -f $RT/.titan_sha ]; then mkdir -p $RT; git -C $REPO archive $RSHA | tar -x -C $RT; echo $RSHA > $RT/.titan_sha; fi
  fi
  RNAME=$($HERE/mkimage.sh $RT name) || exit 2
fi
STAMP=$(date -u +%Y%m%d_%H%M%S)
if [ $NOCACHE = 1 ] || [ $INCL = 1 ] || [ $INPROC = 1 ]; then BASEOUT=${OUTARG:-$DATA/dev/runs/${STAMP}_$TNAME}; TOUT=$BASEOUT/tested; ROUT=$BASEOUT/ref
else BASEOUT=${OUTARG:-$DATA/dev/runs/${STAMP}_$TNAME}; TOUT=$DATA/dev/gcache/$TNAME; ROUT=$DATA/dev/gcache/$RNAME; fi
mkdir -p $BASEOUT $TOUT $ROUT; LOG=$BASEOUT/gates.log
echo "gates: tree $WT ($TNAME) ref ${RNAME} nodes $NODES gates $ONLY mode $MODE inproc=$INPROC include-path=$INCL out $BASEOUT $(date '+%F %T')" | tee -a $LOG
GT0=$(date +%s)
# ---- jobs:  id|which(T/R)|run|steps|variant|extra words     (label = <name>_<steps>; the run dir is <out>/<label>)
JOBS=()
add() { JOBS+=("$1"); }
NK="RESTART=z19zie1S6_auto1 NKCHECK=2 NKCHECKN=1 NKROWPROBE=0"
has 1 && { add "nkm|T|zie1|5|mg|$NK MGROUPS=4"; add "nkg|T|zie1|5|grey|$NK"; }
has 5 && add "mgz3|T|zie3|40|mg|RESTART=fp3a1 MGROUPS=4 PRINTEVERY=1"
has 2 && { add "ic|T|zie1|40|grey|SAVE=sv"; add "ic|R|zie1|40|grey|SAVE=sv"; }
if has 3 || has 4; then
  add "s16|T|zie1|50|grey|RESTART=z19zie1S16 SAVE=sv ETOTMON=1 ETOTLOG=@R@/etot.log"; add "s16|R|zie1|50|grey|RESTART=z19zie1S16 SAVE=sv ETOTMON=1 ETOTLOG=@R@/etot.log"
  add "s8|T|zie3|50|grey|RESTART=z19zie3S8 SAVE=sv ETOTMON=1 ETOTLOG=@R@/etot.log";   add "s8|R|zie3|50|grey|RESTART=z19zie3S8 SAVE=sv ETOTMON=1 ETOTLOG=@R@/etot.log"
fi
# ---- run directory / cache check of a job
jinfo() {  # job -> sets J_ID J_W J_RUN J_N J_VAR J_EXTRA J_TREE J_OUT J_LAB J_DIR
  IFS='|' read -r J_ID J_W J_RUN J_N J_VAR J_EXTRA <<< "$1"
  if [ $J_W = T ]; then J_TREE=$WT; J_OUT=$TOUT; else J_TREE=$RT; J_OUT=$ROUT; fi
  J_LAB=${J_ID}_$J_N; J_DIR=$J_OUT/$J_LAB; J_EXTRA=${J_EXTRA//@R@/$J_DIR}
}
trash() { [ -e "$1" ] && { mkdir -p $DATA/dev/outdated; mv "$1" $DATA/dev/outdated/$(date +%s)_$$_$(basename "$1"); }; return 0; }   # never rm
done_ok() { grep -aq "^===== EXIT 0 " $1/run.out 2>/dev/null; }
env_for() {   # job: exports JLENV / JULIA_DEPOT_PATH / TITAN_PKG of the image (or the include-path env)
  if [ $INCL = 1 ]; then export JLENV=$DATA/proto/jlenv; unset TITAN_PKG
  else eval "$(MKIMAGE_NODE=$2 $HERE/mkimage.sh $J_TREE $J_VAR 2>>$LOG)" || return 1; fi
  export JLENV; return 0
}
runjob() {  # job node
  jinfo "$1"; local node=$2 cpus="$((8*$2))-$((8*$2+7))"
  if [ $NOCACHE = 0 ] && [ $INCL = 0 ] && [ $INPROC = 0 ] && done_ok $J_DIR; then echo "gates: $J_LAB ($J_W) cached" >> $LOG; return 0; fi
  trash $J_DIR
  env_for "$1" $node || { echo "gates: image build failed for $J_LAB" >> $LOG; return 1; }
  local t0=$(date +%s)
  ( export OUT=$J_OUT CKSRC=$CK/$J_RUN/ck REFENV=$CK EQEOS=$CK/eqeos_he_xe.jls PINC=$cpus TH=8
    $REPO/scripts/proto/gate_run.sh $J_LAB $J_TREE $J_RUN $J_N $MODE $J_EXTRA )
  echo "gates: $J_LAB ($J_W) node $node $(( $(date +%s) - t0 )) s $(tail -1 $J_DIR/run.out | cut -c1-60)" >> $LOG
}
runqueue() {  # node jobs...   (one process per job)
  local node=$1; shift; local j; for j in "$@"; do runjob "$j" $node; done
}
runqueue_inproc() {  # node jobs... : per (tree, variant) one julia process
  local node=$1; shift; local j k grp
  local keys; keys=$(for j in "$@"; do IFS='|' read -r a w c d v e <<< "$j"; echo "$w|$v"; done | awk '!s[$0]++')
  for k in $keys; do
    local jf=$BASEOUT/inproc_${k//|/_}_n$node.jobs; : > $jf; local first=
    for j in "$@"; do
      IFS='|' read -r a w c d v e <<< "$j"; [ "$w|$v" = "$k" ] || continue
      jinfo "$j"; [ -z "$first" ] && first=$j
      trash $J_DIR; mkdir -p $J_DIR
      GATE_PRINT=$jf OUT=$J_OUT CKSRC=$CK/$J_RUN/ck REFENV=$CK EQEOS=$CK/eqeos_he_xe.jls $REPO/scripts/proto/gate_run.sh $J_LAB $J_TREE $J_RUN $J_N $MODE $J_EXTRA
    done
    jinfo "$first"; env_for "$first" $node || { echo "gates: image build failed ($k)" >> $LOG; continue; }
    local t0=$(date +%s)
    ( cd $J_TREE/julia && INPROC_TREE=$J_TREE numactl --physcpubind=$((8*node))-$((8*node+7)) --localalloc julia -t 8 --project=$JLENV $REPO/scripts/dev/inproc.jl $jf ) >> $LOG 2>&1
    echo "gates: in-process queue $k node $node $(( $(date +%s) - t0 )) s" >> $LOG
  done
}
# ---- static assignment of the jobs to the nodes (longest first, to the least loaded node); weights ~ expected wall
declare -A W=([ic]=8 [s16]=8 [s8]=8 [mgz3]=8 [nkm]=5 [nkg]=5); LOAD=(); Q=(); for ((i=0;i<NN;i++)); do LOAD[$i]=0; Q[$i]=; done
SORTED=$(for j in "${JOBS[@]}"; do IFS='|' read -r a w c d v e <<< "$j"; echo "${W[$a]} $j"; done | sort -s -k1,1nr | cut -d' ' -f2-)
while IFS= read -r j; do
  [ -n "$j" ] || continue; IFS='|' read -r a w c d v e <<< "$j"
  [ $NOCACHE = 0 ] && [ $INCL = 0 ] && [ $INPROC = 0 ] && [ $w = R ] && { jinfo "$j"; done_ok $J_DIR && continue; }
  b=0; for ((i=1;i<NN;i++)); do [ ${LOAD[$i]} -lt ${LOAD[$b]} ] && b=$i; done
  LOAD[$b]=$(( ${LOAD[$b]} + ${W[$a]} )); Q[$b]+="$j"$'\n'
done <<< "$SORTED"
for ((i=0;i<NN;i++)); do
  mapfile -t JQ <<< "$(printf '%s' "${Q[$i]}" | sed '/^$/d')"; [ -n "${JQ[0]:-}" ] || continue
  echo "gates: node ${NODE[$i]} queue: $(for j in "${JQ[@]}"; do IFS='|' read -r a w c d v e <<< "$j"; echo -n "$a/$w/$v "; done)" | tee -a $LOG
  if [ $INPROC = 1 ]; then runqueue_inproc ${NODE[$i]} "${JQ[@]}" & else runqueue ${NODE[$i]} "${JQ[@]}" & fi
done
wait
echo "gates: all jobs finished, wall $(( $(date +%s) - GT0 )) s" | tee -a $LOG
python3 $HERE/gates_summary.py $TOUT $([ $NEEDREF = 1 ] && echo $ROUT || echo -) $ONLY | tee $BASEOUT/summary.txt
rc=${PIPESTATUS[0]}
echo "gates: total wall $(( $(date +%s) - GT0 )) s, exit $rc" | tee -a $LOG
exit $rc
