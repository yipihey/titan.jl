#!/bin/zsh
# S8a two-way coupled production (ABN02 z = 19, S7c-1 scheme): envelope (popiii_m2c.jl sink run) + protostar (s6b_run.jl), 3 threads each (2 processes).
# Usage: s8a_prod.sh <tag> <lya: mc|mt08|none> [KEY=VAL ... for both processes]
#   mc    RT_LYA=mc   Monte Carlo Lyα (LyaMC.jl), RT_LYMAN=1
#   mt08  RT_LYA=1    the approximate trapped force of S7c-2 (MT08 eq. 22), RT_LYMAN=1
#   none  RT_LYA=0    no Lyα force (RT_LYMAN=1)
# Estimator fix (Mdot -> 0 when the sink starves): S7B_STARVE (default 3), the star follows a falling Mdot: S7B_SLEWDN=0.2, barrier spacing <= S7B_DTMAX yr.
cd "$(dirname "$0")"
tag=$1; lya=$2; shift 2
export RT_LYMAN=1
case $lya in mc) export RT_LYA=mc;; mt08) export RT_LYA=1;; *) export RT_LYA=0;; esac
export S7B_PERFPREFIX=s8a S7B_MSYNC=3 S7B_TCMAX=${S7B_TCMAX:-40000} JT=3 SNAP=${SNAP:-/tmp/s7c1_ck/core2.jls} S7B_SLEWDN=${S7B_SLEWDN:-0.2} S7B_DTMAX=${S7B_DTMAX:-60} S7B_NPROF=${S7B_NPROF:-12}
./s7b_env.sh $tag "$@" > ../lib/popiii/m2c_abn02z19_lines_n200_s7b_$tag.log 2>&1 &
sleep 20
SCRIPT=s7c1_perf.jl PERF=s8a_${tag}_star ./s7b_star.sh $tag CKEVERY=100 "$@" > /tmp/s8a_${tag}_star.out 2>&1 &
wait
