#!/bin/bash
P=/zpool/nvme/data/tabel_scratch/titan_popiii; I=/zpool/nvme/data/tabel_scratch/titan_popiii/images/tree-700e60e
export JLENV=$I/env JULIA_DEPOT_PATH=$I/depot:/home/tabel/.julia TITAN_PKG=1
export OUT=$P/proto/stage2c REFENV=$P/cklib CKSRC=$P/cklib/zie1/ck EQEOS=$P/cklib/eqeos_he_xe.jls PINC=64-71 SNAPEVERY=1 DIAGB=/home/tabel/Projects/titan.jl/.claude/worktrees/proto/docs/logs/proto/stage2c/diagB.jl
/home/tabel/Projects/titan.jl/.claude/worktrees/proto/scripts/proto/gate_run.sh b $P/proto/stage2c/wt700 zie1 1 ref RESTART=z19zie1S16 CONVDIAG=1
