#!/bin/zsh
# S7c-1: single-process restart of the S6b/S7b z19 star from a checkpoint with the tabulated Mdot(M*) (reproduces the precursor stall at ~32 Msun).
# Usage: s7c1_repro.sh <tag> <ckdir> <ckname> [KEY=VAL ...]
cd "$(dirname "$0")"
tag=$1; ckd=$2; ckn=$3; shift 3
export TITAN_NOLIB=1 OUT=s7c1_$tag PRECTHICK=1 BCSMOOTH=0.03 BCSLOPE=8 BCFALLBACK=1 ESLIM=1.05 TOL=0.05 NIN=80 QC=1.083 QSPLIT=0.3 NENV=70 PSI=0
export MDOTFILE=${MDOT:-../lib/popiii/s6c__z19a2_mdot_smooth.csv} MDOTRAMP=4.41e-3 M0GUESS=0.22 MDOTRAMPP=2.5 MAXSTEPS=${MAXSTEPS:-400} EVERY=10 CSVEVERY=200 CKEVERY=100000
export CKDIR=$ckd RESTART=$ckn
for kv in "$@"; do export "$kv"; done
exec julia -t 3 --project=env s6b_run.jl
