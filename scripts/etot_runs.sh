#!/bin/bash
# Phase 4 item 2: internal-energy (std) vs conservative total-energy (etot[, dual switch]) rows
cd "$(dirname "$0")/../julia"
export TITAN_BITWISE=1 TITAN_SET="jstepe=3000"
J="nice /usr/local/bin/julia --startup-file=no --project=env etotal_run.jl"
for a in "std" "etot" "etot 0.1" "etot 0.3"; do $J SUBC 8 $a; $J SUPC 9 $a; $J SOD 3 $a; done
echo ETOT DONE
