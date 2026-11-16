#!/bin/bash
# bitrun.sh <tree> <rundir> <name> [extra env...]: restart zdiagB's S5 checkpoint for MAXSTEPS steps with zdiagB's S6 settings.
TREE=$1; R=$2; NAME=$3; shift 3
JLENV=/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/jlenv
S6=$(cat /tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/g2/s6env.txt)
cd $TREE/julia
echo "$(date '+%F %T') start $NAME tree $TREE extra: $*" >> $R/procedure.log
env $S6 RUNREC=$R POPIII_OUT=$R CKDIR=$R/ck TAG=_$NAME MATTHREADS=8 RESTART=z19zdiagBS5 SAVE=$NAME MAXSTEPS=${MAXSTEPS:-8} TAGRUN=g2bit "$@" \
    numactl --physcpubind=104-111 --membind=13 julia -t 8 --project=$JLENV popiii_m2c.jl lines 128 > $R/$NAME.out 2>&1 < /dev/null
echo "$(date '+%F %T') end $NAME exit $?" >> $R/procedure.log
