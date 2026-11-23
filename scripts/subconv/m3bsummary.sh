#!/bin/bash
# M3b: summary of the pilot runs (run.out, perf, snapshots): for each <name> in tests/m3b: steps, Newton-failure rejections, error rejections, min accepted dt (snapshot steps),
# Newton iteration statistics (nwt.txt), energy ledger (etot.log), and the snapshot table of m3bsnap.jl. Usage: m3bsummary.sh <name> ...  (runs julia for the table)
T=/zpool/nvme/data/tabel_scratch/titan_popiii/subconv/tests/m3b; S=$(cd "$(dirname "$0")" && pwd); J="julia --project=/zpool/nvme/data/tabel_scratch/titan_popiii/subconv/jlenv"
for n in "$@"; do
  d=$T/$n; echo "=== $n"
  grep -m1 'run ended' $d/run.out | cut -c1-160
  echo "Newton-failure rejections: $(grep -c 'DTERR reject.*Newton failure' $d/run.out); error rejections: $(grep -c 'DTERR reject.*err ' $d/run.out)"
  grep -m1 'DTERR statistics' $d/run.out | cut -c1-300
  grep -m1 'MLTBUOY zone thermodynamics' $d/run.out | cut -c1-200
  [ -f $d/nwt.txt ] && python3 $S/m3bnwt.py $d
  [ -f $d/etot.log ] && python3 $S/m3betot.py $d
  $J $S/m3bsnap.jl $d 5 2>&1 | head -60
done
