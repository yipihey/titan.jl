#!/bin/bash
# variants at fixed dt, one process: probe i uses variant i of TDCVARS
cd /zpool/nvme/data/tabel_scratch/titan_popiii/subconv/tests/m3wall
export PINC=${PINC:-72-79}
N=$1; DT=$2; shift 2
V="$1"; shift
NV=$(echo "$V" | awk -F'|' '{print NF}')
P=$(python3 -c "print(','.join(['$DT']*$NV))")
./run.sh $N 2 TDCVARS="$V" NWPROBE=$P NWPROBEAT=0 "$@"
