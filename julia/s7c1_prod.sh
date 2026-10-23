#!/bin/zsh
# S7c-1 two-way coupled production (ABN02 z = 19): envelope (popiii_m2c.jl sink run) + protostar (s6b_run.jl) exchanging at barriers, 3 threads each (2 processes).
# Usage: s7c1_prod.sh <tag> <RT_LYMAN 0|1> <RT_LYA 0|1> [LYASURV 0|1] [KEY=VAL ... for both processes]
# Star: PROTGAMMA = 0.25 (s7b_star.sh), Perf records lib/popiii/perf/s7c1_<tag>_*; envelope: S7B_PERFPREFIX = s7c1.
cd "$(dirname "$0")"
tag=$1; export RT_LYMAN=$2 RT_LYA=$3 S7B_LYASURV=${4:-0}; if [ $# -ge 4 ]; then shift 4; else shift $#; fi
export S7B_PERFPREFIX=s7c1 S7B_MSYNC=3 S7B_TCMAX=${S7B_TCMAX:-60000} JT=3 SNAP=${SNAP:-/tmp/s7c1_ck/core2.jls}
./s7b_env.sh $tag "$@" > ../lib/popiii/m2c_abn02z19_lines_n200_s7b_$tag.log 2>&1 &
sleep 20
SCRIPT=s7c1_perf.jl PERF=s7c1_${tag}_star ./s7b_star.sh $tag CKEVERY=100 "$@" > /tmp/s7c1_${tag}_star.out 2>&1 &
wait
