#!/bin/zsh
# S7b protostar process (s6b_run.jl, ABN02 z = 19 star fed by the collapse table until the coupling starts, then by the envelope process). Usage: s7b_star.sh <tag> [KEY=VAL ...]
# MDOT = Mdot(M*) table of the collapse for phase :pre (default lib/popiii/s6c__z19a2_mdot_smooth.csv); settings of the S6b Part 3 / E2E-A protostar stage.
cd "$(dirname "$0")"
tag=$1; shift
export TITAN_NOLIB=1 OUT=s7b_$tag PRECTHICK=1 BCSMOOTH=0.03 BCSLOPE=8 BCFALLBACK=1 ESLIM=1.05 TOL=0.05 NIN=80 QC=1.083 QSPLIT=0.3 NENV=70 PSI=0
export MDOTFILE=${MDOT:-../lib/popiii/s6c__z19a2_mdot_smooth.csv} MDOTRAMP=4.41e-3 M0GUESS=0.22 MDOTRAMPP=2.5 MAXSTEPS=600000 EVERY=50 CSVEVERY=200 CKEVERY=300 CKDIR=/tmp/s6_s7b_$tag SAVE=s7b_${tag}_p3
export S7B_STAR=1 S7B_TAG=$tag
for kv in "$@"; do export "$kv"; done
exec julia -t 4 --project=env s6b_run.jl
