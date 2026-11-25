#!/bin/bash
# run.sh <label> <cores> <steps> [extra]
P=/zpool/nvme/data/tabel_scratch/titan_popiii; S=$P/proto/s16diag
L=$1; CORES=$2; N=$3; shift 3
export OUT=$S REFENV=$P/cklib CKSRC=$P/cklib/zie1/ck EQEOS=$P/cklib/eqeos_he_xe.jls PINC=$CORES SNAPEVERY=${SNAPEVERY:-1}
TREE=/home/tabel/Projects/titan.jl/.claude/worktrees/proto
$S/gate_run_snap.sh $L $TREE zie1 $N proto RESTART=z19zie1S16 CONVDIAG=1 "$@"
echo done > $S/$L.done
