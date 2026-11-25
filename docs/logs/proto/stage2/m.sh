#!/bin/bash
# m.sh <rundir> <prefix>: metrics of a stage-2 run
P=/home/tabel/Projects/titan.jl/.claude/worktrees/proto/scripts/proto; J=/zpool/nvme/data/tabel_scratch/titan_popiii/proto/jlenv
python3 $P/parse_run.py $1 $2; python3 $P/ledger.py $1/etot.log $2.etot
julia --startup-file=no --project=$J $P/pilot_snap.jl $1 $2 2>&1 | tail -20
grep -c "^DTERR jstep" $1/run.out | sed "s/^/$2.dterr_lines /"
