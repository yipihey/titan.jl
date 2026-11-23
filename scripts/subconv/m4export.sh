#!/bin/bash
# m4export.sh <case>...   export M4 cases (B0, B1, B2, ...) to $M4/exports/<case>/{frames.f32,frames.json,track.csv}
# B0 = production replay zie1rp (all steps exported; m4compare.py restricts to step >= 290). Lane dirs are only read.
set -u
S=/zpool/nvme/data/tabel_scratch/titan_popiii/subconv
WT=$(cd "$(dirname "$0")/../.." && pwd)
for c in "$@"; do
  if [ "$c" = B0 ]; then d=$S/runs/zie1rp; else d=$S/m4/$c; fi
  mkdir -p $S/m4/exports/$c
  echo "== $c ($d)"
  (cd $WT/julia && numactl --physcpubind=96-99 --localalloc julia -t 2 --project=$S/jlenv $WT/scripts/subconv/m4extra.jl $d $S/m4/exports/$c) 2>&1 | grep -v "^┌\|^└" | tail -5
done
