#!/bin/zsh
# S7c-2 one-way production runs (ABN02 z = 19, F1 set-up of S7b: core2 snapshot, S6b track as the source, feedback from M_pt = 3 Msun).
# Usage: s7c2_prod.sh <tag> <RT_LYMAN 0|1> <RT_LYA 0|1>     (tags used: c2L = Lyman wings only, c2A = trapped Lyα only, c2LA = both)
cd "$(dirname "$0")"
tag=$1
export RT_LYMAN=$2 RT_LYA=$3 S7B_PERFPREFIX=s7c2
export SNAP=${SNAP:-/private/tmp/w/snaps/snaps_abn02z19_lines_n200_core2.jls}
export S7B_TRACK=../lib/popiii/s6b_p3_z19f.csv S7B_MSYNC=3 S7B_NPROF=12
./s7b_env.sh $tag > ../lib/popiii/m2c_abn02z19_lines_n200_s7b_$tag.log 2>&1
