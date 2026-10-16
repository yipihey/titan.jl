#!/bin/bash
# Reproduce all TITAN 1.2 reference runs and CSV extractions.
# Build first:  cd $(cd "$(dirname "$0")/.." && pwd)/fortran && make codes XFLAGS="-fno-automatic -finit-local-zero"
set -u
T=$(cd "$(dirname "$0")/.." && pwd)
R=$T/scripts/run_bench.sh
RR=$T/reference_runs
$R 03 SOD; $R 04 BLST; $R 05 SDOV
$R 06 HEAT
mkdir -p $RR/COOL && cp $RR/HEAT/HEAT.dmp $RR/COOL/COOL.dmp   # COOL restarts from HEAT dump 4
$R 07 COOL
$R 08 SUBC; $R 09 SUPC
$R 10 RLST                                                   # leg 1: steps 1..2001
mkdir -p $RR/RLST_leg2 && cd $RR/RLST_leg2 && cp ../RLST/RLST.dmp . && cp $T/fortran/titan.in titan.in
t0=$(date +%s.%N); timeout 600 $T/fortran/titan.x </dev/null >stdout.txt 2>stderr.txt; rc=$?
t1=$(date +%s.%N); echo "wall=$(echo "$t1 - $t0" | bc) exit=$rc" > time.txt   # leg 2: restart dump 19 / step 1801
# extraction: final model + every archived model
S=$T/scripts/read_dump.py
cd $RR
ex(){ d=$1; f=$2; g=$3; n=$4
  python3 -I $S $d/$f $g $n $d/${f%.dmp}_final.csv
  last=$(python3 -I -c "import sys;sys.path.insert(0,'$T/scripts');import read_dump as r,numpy as np;print(np.frombuffer(r.read_record('$d/$f',2)[:8],'<i4')[0])")
  mkdir -p $d/dumps
  for j in $(seq 1 $last); do python3 -I $S $d/$f $g $n $d/dumps/$(printf "%s_d%03d.csv" ${f%.dmp} $j) $j >/dev/null; done; }
ex SOD SOD.dmp 4 100; ex BLST BLST.dmp 4 200; ex SDOV SDOV.dmp 5 100
ex HEAT HEAT.dmp 5 105; ex COOL COOL.dmp 5 105
ex SUBC SUBC.dmp 5 100; ex SUPC SUPC.dmp 5 100; ex RLST_leg2 RLST.dmp 5 100
