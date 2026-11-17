#!/bin/bash
# exprun.sh <tree> <rundir> <name> <maxsteps> [extra env...]: zdiagB's S6 settings from checkpoint S5 (copied into <rundir>/ck),
# RUNREC snapshots every 5 steps, checkpoints every 100 steps; cores 96-103.
TREE=$1; R=$2; NAME=$3; MS=$4; shift 4
JLENV=/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/jlenv
S6=$(cat /tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/g2/s6env.txt)
cd $TREE/julia
echo "$(date '+%F %T') start $NAME tree $TREE maxsteps $MS extra: $*" >> $R/procedure.log
env $S6 RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck TAG=_$NAME MATTHREADS=8 SNAPEVERY=5 CKEVERY=100 RESTART=${RST:-z19zdiagBS5} SAVE=$NAME MAXSTEPS=$MS TAGRUN=g2exp "$@" \
    numactl --physcpubind=104-111 --membind=13 julia -t 8 --project=$JLENV popiii_m2c.jl lines 128 > $R/$NAME.out 2>&1 < /dev/null
echo "$(date '+%F %T') end $NAME exit $? $(grep 'run ended' $R/$NAME.out | tail -1)" >> $R/procedure.log
