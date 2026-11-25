#!/bin/bash
# inproc_test.sh <worktree> [node-sep] [node-inproc]   feasibility of several PopIIIRun.main calls in one julia process (docs/POPIII_DEVTOOLS.md).
# Jobs: (1) zie1 IC 3 steps, (2) zie3 restart fp3a1 20 steps, (3) zie1 IC 3 steps again; each run (a) as separate processes and (b) in one process
# (scripts/dev/inproc.jl). Compares per job the step lines (prod_env/cmp.sh filter) and the decompressed tstate byte for byte; also (1) vs (3).
set -u
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd $HERE/../.. && pwd); WT=$(cd ${1:?worktree} && pwd); NS=${2:-3}; NI=${3:-4}
DATA=${DATA:-/zpool/nvme/data/tabel_scratch/titan_popiii}; CK=$DATA/cklib; BASE=$DATA/dev/inproc_test_$(date -u +%H%M%S); mkdir -p $BASE/sep $BASE/inp
eval "$(MKIMAGE_NODE=$NS $HERE/mkimage.sh $WT grey)" || exit 2
JOBS=("a1 zie1 3 SAVE=sv" "b2 zie3 20 RESTART=fp3a1 SAVE=sv" "c3 zie1 3 SAVE=sv")
common() { export CKSRC=$CK/$2/ck REFENV=$CK EQEOS=$CK/eqeos_he_xe.jls TH=8; }
inp() {
  local jf=$BASE/inp.jobs; : > $jf
  for j in "${JOBS[@]}"; do read -r L RUN N X1 X2 <<< "$j"; ( common x $RUN; export OUT=$BASE/inp GATE_PRINT=$jf; $REPO/scripts/proto/gate_run.sh $L $WT $RUN $N proto $X1 $X2 ); done
  T0=$(date +%s)
  ( cd $WT/julia && INPROC_TREE=$WT numactl --physcpubind=$((8*NI))-$((8*NI+7)) --localalloc julia -t 8 --project=$JLENV $HERE/inproc.jl $jf ) > $BASE/inp.log 2>&1
  echo "inproc total wall $(( $(date +%s) - T0 )) s" >> $BASE/inp.log
}
sep() {
  T0=$(date +%s)
  for j in "${JOBS[@]}"; do read -r L RUN N X1 X2 <<< "$j"; ( common x $RUN; export OUT=$BASE/sep PINC=$((8*NS))-$((8*NS+7)); $REPO/scripts/proto/gate_run.sh $L $WT $RUN $N proto $X1 $X2 ); done
  echo "separate total wall $(( $(date +%s) - T0 )) s" > $BASE/sep.log
}
sep & inp & wait
cat $BASE/sep.log; tail -n 3 $BASE/inp.log
python3 - $BASE <<'EOF'
import sys, re, gzip, os
B = sys.argv[1]; F = re.compile(r"^ +[0-9]+ +[0-9.]+ +[0-9.e+]+ |^DTERR|^ *divergence")
def sl(m, l): return [x for x in open(f"{B}/{m}/{l}/run.out", errors="replace").read().splitlines() if F.match(x)]
def ts(m, l): return gzip.open(f"{B}/{m}/{l}/ck/m2c_sv.tstate.gz").read()
for l in ("a1", "b2", "c3"):
    try:
        a, b = sl("sep", l), sl("inp", l); ta, tb = ts("sep", l), ts("inp", l)
        print(f"{l}: step lines {len(a)} vs {len(b)}: {'IDENTICAL' if a == b else 'DIFFER'}; tstate {len(ta)} vs {len(tb)} B: {'IDENTICAL' if ta == tb else 'DIFFER'}")
        if a != b:
            for i, (x, y) in enumerate(zip(a, b)):
                if x != y: print("   first difference, step line", i, "\n   sep:", x[:160], "\n   inp:", y[:160]); break
    except Exception as e: print(l, "ERROR", e)
try: print("in-process a1 vs c3 (same job twice):", "IDENTICAL" if sl("inp", "a1") == sl("inp", "c3") and ts("inp", "a1") == ts("inp", "c3") else "DIFFER")
except Exception as e: print("a1 vs c3 ERROR", e)
EOF
echo "results in $BASE"
