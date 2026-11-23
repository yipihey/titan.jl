#!/bin/bash
# copy the branch's julia/*.jl (except TitanCore.jl) into the scratch worktree of the M3c tests
S=$(cd "$(dirname "$0")/../.." && pwd); W=/zpool/nvme/data/tabel_scratch/titan_popiii/subconv/tests/m3dbg/wt
for f in $S/julia/*.jl; do b=$(basename $f); [ $b = TitanCore.jl ] || cmp -s $f $W/julia/$b || { cp $f $W/julia/$b; echo synced $b; }; done
