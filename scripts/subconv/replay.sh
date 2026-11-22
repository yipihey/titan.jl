#!/bin/bash
# Replay of a production run's first segment from the IC, with full checkpoints at listed steps (CKDUMP), into its own run
# directory. Code: a worktree at the production commit (WT). Settings: the production run.env verbatim (+ CKDUMP, SAVE, MAXSTEPS).
# Resumable: if a checkpoint <SAVE>_auto{0,1} exists, the newest is restarted (restarts are bit-identical, popiii_m2c S6-A).
# Usage: replay.sh <prodrun> <replayrun> <WT> <MAXSTEPS> <CKDUMP> [PIN]
set -u
DATA=/zpool/nvme/data/tabel_scratch/titan_popiii; S=$DATA/subconv
P=$DATA/runs/$1; R=$S/runs/$2; WT=$3; MAX=$4; DUMP=$5; PIN=${6:-"numactl --physcpubind=72-79 --localalloc"}
mkdir -p $R/ck; CK=$R/ck; log() { echo "$(date '+%F %T') $*" >> $R/procedure.log; }
PIN0=$PIN; . $P/launch.env; PIN=$PIN0
COMMON="RUNREC=$R POPIII_OUT=$R CKDIR=$CK EQEOS_CACHE=$R/eqeos_he_xe.jls TIMING=1 DTDIAG=1 TAG=_$2 $(cat $P/run.env)"
[ -f $R/eqeos_he_xe.jls ] || cp $P/eqeos_he_xe.jls $R/
last=$(ls -t $CK/m2c_rp_auto?.aux 2>/dev/null | head -1); st=; [ -n "$last" ] && st="RESTART=$(basename $last .aux | sed 's/^m2c_//')"
done_=$(ls $CK 2>/dev/null | grep -c _k)
log "replay of $1 (commit $(git -C $WT rev-parse --short HEAD)), MAXSTEPS=$MAX CKDUMP=$DUMP ${st:-from the IC}; $done_ dump files present; $PIN"
cd $WT/julia
env $COMMON $st SAVE=rp CKDUMP=$DUMP CKEVERY=50 MAXSTEPS=$MAX setsid $PIN julia -t ${TH:-8} --project=$S/jlenv popiii_m2c.jl lines ${NZ:-128} > $R/seg_rp_$(date +%s).out 2>&1 < /dev/null
log "replay ended: $(grep -E 'run ended|stopped' $R/seg_rp_*.out | tail -1)"
