#!/bin/bash
# Phase 4 item 4: steady-state timing fortran vs jfast, tier (b)/(c) for jfast
cd "$(dirname "$0")/../julia"
export TITAN_BITWISE=1
J="nice /usr/local/bin/julia --startup-file=no --project=env phase3_run.jl"
for m in fortran jfast; do $J SUBC 8 $m; $J SUPC 9 $m; TITAN_SET="jstepe=3000" $J SOD 3 $m; done
$J SUBC 8 jfast tierb; $J SUPC 9 jfast tierb; TITAN_SET="jstepe=3000" $J SOD 3 jfast tierb
nice /usr/local/bin/julia --startup-file=no tierc.jl jfast
echo JFAST DONE
