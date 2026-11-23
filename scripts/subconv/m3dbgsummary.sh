#!/bin/bash
# M3c: one line per run directory of $T/m3dbg: the Newton attempts (jstep, theta, dt, ok/iterations) from nwt.txt. Usage: m3dbgsummary.sh <name>...
T=/zpool/nvme/data/tabel_scratch/titan_popiii/subconv/tests/m3dbg
for v in "$@"; do echo "== $v"; grep -h "NWT head\|NWT end" $T/$v/nwt.txt | sed 's/NWT head att \([0-9]*\) jstep \([0-9]*\) thet \([0-9.]*\) dtime \([0-9.e+-]*\).*/\2 th\3 dt\4/;s/NWT end att [0-9]* ok \(.*\) iter \(.*\)/\1\/\2/' | paste -d' ' - - ; done
