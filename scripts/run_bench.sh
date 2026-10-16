#!/bin/bash
# usage: run_bench.sh NN NAME [exe] [outdir_root]
# runs TITAN benchmark NN (titan.inNN) in reference_runs/NAME with a 600 s cap
set -u
NN=$1; NAME=$2
EXE=${3:-$(cd "$(dirname "$0")/.." && pwd)/fortran/titan.x}
ROOT=${4:-$(cd "$(dirname "$0")/.." && pwd)/reference_runs}
SRC=$(cd "$(dirname "$0")/.." && pwd)/fortran
D=$ROOT/$NAME; mkdir -p $D; cd $D
cp $SRC/titan.in$NN titan.in
[ -n "${PREDUMP:-}" ] && cp "$PREDUMP" .
t0=$(date +%s.%N)
timeout 600 $EXE < /dev/null > stdout.txt 2> stderr.txt
rc=$?
t1=$(date +%s.%N)
echo "wall=$(echo "$t1 - $t0" | bc) exit=$rc" > time.txt
cat time.txt
