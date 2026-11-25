#!/bin/bash
# Stage 2c task A: 40 steps with SNAPEVERY=1 from the relaxed TDC checkpoint (breathing amplitude / period at fine cadence)
P=/zpool/nvme/data/tabel_scratch/titan_popiii; IMG=$P/images/popiii-proto-6241f5a
export JLENV=$IMG/env JULIA_DEPOT_PATH=$IMG/depot:/home/tabel/.julia TITAN_PKG=1
export OUT=$P/proto/stage2c REFENV=$P/cklib CKSRC=$P/cklib/zie1/S16_tdc EQEOS=$P/cklib/eqeos_he_xe.jls PINC=72-79 SNAPEVERY=1
RUNOUT=run_s00.out PRINTEVERY=1 /home/tabel/Projects/titan.jl/.claude/worktrees/proto/scripts/proto/gate_run.sh Afine /home/tabel/Projects/titan.jl/.claude/worktrees/proto zie1 40 proto RESTART=z19zie1S16_tdc CONVDIAG=1
