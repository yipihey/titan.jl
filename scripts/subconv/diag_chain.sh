#!/bin/bash
# sequential epochs: diag_chain.sh <prodrun> <replayrun> <MAXSTEPS> <nz> k1 k2 ...
H=$(dirname "$0"); a=$1; b=$2; m=$3; z=$4; shift 4
for k in "$@"; do bash $H/diag_epoch.sh $a $b $k $m $z; done
