#!/bin/zsh
# S7c-1 validation of the S6a-fix / S6b results with the changed defaults (EXTRAP 0, ESTAU 1, LSTARQ 1e-3): constant-Mdot protostar, settings of s6b_pv6 (Mdot 1.1e-3) and s6b_pt5 (Mdot_fid 4.41e-3).
# Usage: s7c1_val.sh <tag> <Mdot Msun/yr> <Mend Msun> [KEY=VAL ...]   -> lib/popiii/s6b_s7c1_<tag>.{log,csv}, lib/popiii/perf/s7c1_val_<tag>_*
cd "$(dirname "$0")"
tag=$1; md=$2; me=$3; shift 3
export TITAN_NOLIB=1 OUT=s7c1_$tag MDOT=$md MEND=$me PRECTHICK=1 BCSMOOTH=0.03 BCSLOPE=8 BCFALLBACK=1 ESLIM=1.05 TOL=0.05 NIN=80 QC=1.083 QSPLIT=0.3 NENV=70 PSI=0 MAXSTEPS=400000 EVERY=200 CSVEVERY=500 PERF=s7c1_val_$tag
for kv in "$@"; do export "$kv"; done
exec julia -t ${JT:-3} --project=env s7c1_perf.jl
