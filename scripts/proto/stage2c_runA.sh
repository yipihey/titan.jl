#!/bin/bash
# Task A (Stage 2c): TDC (HEAD code) from the zie1 S16 copy, segmented and resumable, until dt is steady or 2000 steps / 24 h.
# usage: stage2c_runA.sh  (env: SEGSTEPS=250 MAXTOT=2000 WALLMAX=86400 PINC=64-71). Re-running resumes from the checkpoint tdcA in A/ck.
P=/zpool/nvme/data/tabel_scratch/titan_popiii; OUT=$P/proto/stage2c; R=$OUT/A; mkdir -p $R
TREE=/home/tabel/Projects/titan.jl/.claude/worktrees/proto; S=$TREE/scripts/proto
IMG=$P/images/popiii-proto-6241f5a                 # scripts/dev/mkimage.sh $TREE at 6241f5a (HEAD julia code; later commits touch docs/scripts only, so the image stays valid)
export JLENV=$IMG/env JULIA_DEPOT_PATH=$IMG/depot:/home/tabel/.julia TITAN_PKG=1
export OUT REFENV=$P/cklib CKSRC=$P/cklib/zie1/ck EQEOS=$P/cklib/eqeos_he_xe.jls PINC=${PINC:-64-71} SNAPEVERY=10
SEGSTEPS=${SEGSTEPS:-250}; MAXTOT=${MAXTOT:-2000}; WALLMAX=${WALLMAX:-86400}; T0=$(date +%s)
START=63334                                         # zie1 S16 checkpoint step
while :; do
  seg=$(ls $R/run_s*.out 2>/dev/null | wc -l)
  last=$(python3 $S/stage2c_steady.py $R --table 2>/dev/null | awk 'NF==4{k=$1} END{print k+0}'); [ "$last" = 0 ] && last=$START
  done_steps=$((last - START)); el=$(( $(date +%s) - T0 ))
  [ -f $R/STOP ] && { echo stop-file; break; }
  [ $done_steps -ge $MAXTOT ] && { echo "2000 steps"; break; }
  [ $el -gt $WALLMAX ] && { echo "wall cap"; break; }
  [ $seg -gt 0 ] && python3 $S/stage2c_steady.py $R >/dev/null && { echo steady; break; }
  if [ -f $R/ck/m2c_tdcA.aux ]; then RS=tdcA; else RS=z19zie1S16; fi
  n=$((MAXTOT - done_steps)); [ $n -gt $SEGSTEPS ] && n=$SEGSTEPS
  RUNOUT=run_s$(printf %02d $seg).out PRINTEVERY=10 $S/gate_run.sh A $TREE zie1 $n proto RESTART=$RS SAVE=tdcA CKEVERY=50 CONVDIAG=1
done
echo "runA finished: last step $last, wall $el s" > $R/runA.done
