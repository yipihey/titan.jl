#!/bin/bash
cd "$(dirname "$0")/../julia"
export TITAN_BITWISE=1          # compare against the bitwise-TITAN reference trajectories
J="nice /usr/local/bin/julia --startup-file=no --project=env phase3_run.jl"
$J SUPC 9 jfull; TITAN_SET="jstepe=3000" $J SOD 3 jfull
$J SUBC 8 jfull tierb; $J SUPC 9 jfull tierb; TITAN_SET="jstepe=3000" $J SOD 3 jfull tierb
nice /usr/local/bin/julia --startup-file=no tierc.jl jfull
echo JFULL DONE
