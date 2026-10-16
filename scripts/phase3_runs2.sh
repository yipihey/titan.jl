#!/bin/bash
cd "$(dirname "$0")/../julia"
J="nice /usr/local/bin/julia --startup-file=no --project=env phase3_run.jl"
$J SUBC 8 jall; $J SUPC 9 jall; TITAN_SET="jstepe=3000" $J SOD 3 jall
$J SUBC 8 jall tierb; $J SUPC 9 jall tierb
TITAN_FIX_CLOCK=1 $J SUBC 8 fortran > /dev/null 2>&1; mv ../lib/p2_runs/SUBC_fortran ../lib/p2_runs/SUBC_fixclock; mv ../lib/p3_runs/run_SUBC_fortran.log ../lib/p3_runs/run_SUBC_fixclock.log
$J SUBC 8 fortran
nice /usr/local/bin/julia --startup-file=no tierc.jl jall
echo P3B DONE
