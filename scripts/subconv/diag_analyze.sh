#!/bin/bash
# diag_analyze.sh <name e.g. zie1_k115> <replaysnaps.jls> <steps...>: bit-identity vs the replay + series + tables at the given steps -> docs/logs/subconv/m1diag/
S=/zpool/nvme/data/tabel_scratch/titan_popiii/subconv; H=$(cd "$(dirname "$0")" && pwd); L=$H/../../docs/logs/subconv/m1diag; N=$1; RP=$2; shift 2
SN=$(ls $S/tests/m1diag/$N/snaps_*.jls); J="julia --project=$S/jlenv"
cd $S/wt_diag/julia
$J $H/diag_prodcmp.jl $RP $SN > $L/${N}_bitid_vs_replay.txt 2>&1
$J $H/diag_layer.jl $SN series > $L/${N}_series.txt 2>&1
$J $H/diag_layer.jl $SN gates > $L/${N}_gates.txt 2>&1
$J $H/diag_layer.jl $SN summary > $L/${N}_summary.txt 2>&1
for st in "$@"; do $J $H/diag_layer.jl $SN $st 3 > $L/${N}_s$st.txt 2>&1; done
