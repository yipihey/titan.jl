#!/bin/bash
# M3b pilot lane: runs the listed jobs one after the other on one 8-core slot (code: CODE, a scratch worktree with TitanCore.jl of aace2a0).
# Usage: m3blane.sh <cores> <jobfile>   jobfile line: <name> <steps> <env words...>   (see m3brun.sh). ETOTLOG, NWTRACE are set per job into the run directory.
set -u
S=$(cd "$(dirname "$0")" && pwd); T=/zpool/nvme/data/tabel_scratch/titan_popiii/subconv/tests/m3b
while read -r name steps rest; do
  [ -z "${name:-}" ] && continue
  case $name in \#*) continue;; esac
  mkdir -p $T/$name
  PINC=$1 ETOTLOG=$T/$name/etot.log NWTRACE=$T/$name/nwt.txt $S/m3brun.sh $name $steps $rest
done < $2
