#!/bin/bash
# run2.sh <label> <tree> <cores> <run zie1|zie3> <ckname> <steps> [extra VAR=val...]  -> stage2/<label>
P=/zpool/nvme/data/tabel_scratch/titan_popiii; S=$P/proto/stage2
L=$1; TREE=$2; CORES=$3; RUN=$4; CK=$5; N=$6; shift 6
export OUT=$S REFENV=$P/cklib CKSRC=$P/cklib/$RUN/ck EQEOS=$P/cklib/eqeos_he_xe.jls PINC=$CORES
G=/home/tabel/Projects/titan.jl/.claude/worktrees/proto/scripts/proto/gate_run.sh
$G $L $TREE $RUN $N proto RESTART=$CK ETOTMON=1 ETOTLOG=$S/$L/etot.log CONVDIAG=1 "$@"
echo done > $S/$L.done
