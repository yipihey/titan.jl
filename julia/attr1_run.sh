#!/bin/zsh
# ATTR-1: one coupled (envelope + protostar) run of the FULL-1 chain with single configuration switches.
# Usage: attr1_run.sh <name> <snap: fb|core2> <thickcont: 1|0> <starvecut: 3|0> <mdotstop> [lya: 0|mc] ["K=V K=V" extra envelope env] ["K=V" extra env of both processes]
#   fb    = formation snapshot of E2E-A (/tmp/attr1/fb_c2.jls, FULL-1);  core2 = S6b core2 snapshot (/tmp/attr1/core2.jls, S7b..S8a)
#   thickcont 1 = THICKCONT=3 TEQMAX=1e30 in the envelope (FULL-1);  0 = defaults (THICKCONT=0, TEQMAX=1.2e4, as s7b_env.sh)
#   starvecut 3 = S7B_STARVECUT=3 (FULL-1) / 0 off;  mdotstop = S7B_MDOTSTOP [Msun/yr] (FULL-1 1e-3, S8a 0)
cd "$(dirname "$0")"
name=$1; snap=$2; thick=$3; cut=$4; mstop=$5; lya=${6:-0}; xenv=$7; both=$8
tag=attr1_$name
W=../lib/popiii/e2e_work/$tag; mkdir -p $W
src=/tmp/attr1/fb_c2.jls; [[ $snap == core2 ]] && src=/tmp/attr1/core2.jls
printf 'stage=formation\nsnapshot=%s\nreused=attr1\n' $src > $W/formation.done
extra="$xenv $both"; [[ $thick == 0 ]] && extra="THICKCONT=0 TEQMAX=1.2e4 $xenv $both"
export E2E_FB=1 E2E_TAG=$tag E2E_PERF=$tag E2E_THREADS=${E2E_THREADS:-4} E2E_LYA=$lya E2E_STARVECUT=$cut E2E_MDOTSTOP_FB=$mstop E2E_HI21EVERY=1000000 E2E_STARFREEZE=1
export E2E_EXTRA_COUPLED_ENV="$extra" E2E_EXTRA_COUPLED_STAR="$both"
exec julia --project=env e2e_chain.jl coupled
