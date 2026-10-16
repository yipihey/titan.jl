#!/bin/bash
# Run all benchmarks with a given executable and byte-compare .dmp/.out to reference_runs.
EXE=$1; OUT=$2; T="$(cd "$(dirname "$0")/.." && pwd)"; R=$T/scripts/run_bench.sh
for p in "03 SOD" "04 BLST" "05 SDOV" "06 HEAT" "08 SUBC" "09 SUPC" "10 RLST"; do set -- $p
  nice $R $1 $2 $EXE $OUT >/dev/null
  if [ $2 = HEAT ]; then mkdir -p $OUT/COOL; cp $OUT/HEAT/HEAT.dmp $OUT/COOL/COOL.dmp; nice $R 07 COOL $EXE $OUT >/dev/null; fi
done
mkdir -p $OUT/RLST_leg2; cd $OUT/RLST_leg2; cp ../RLST/RLST.dmp .; cp $T/fortran/titan.in .; nice timeout 600 $EXE </dev/null >stdout.txt 2>&1
for n in SOD BLST SDOV HEAT COOL SUBC SUPC RLST RLST_leg2; do f=${n%_leg2}
  a=$(cmp -s $OUT/$n/$f.dmp $T/reference_runs/$n/$f.dmp && echo same || echo DIFF)
  b=$(cmp -s $OUT/$n/$f.out $T/reference_runs/$n/$f.out && echo same || echo DIFF)
  echo "$n dmp:$a out:$b"; done
