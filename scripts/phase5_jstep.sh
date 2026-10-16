#!/bin/bash
# Phase 5b: whole step in Julia (jstep) vs fixed all-Fortran reference (std); tier (c)
cd "$(dirname "$0")/../julia"
J="nice julia --startup-file=no --project=env"
$J phase3_run.jl SUBC 8 jstep; $J phase3_run.jl SUPC 9 jstep; TITAN_SET="jstepe=3000" $J phase3_run.jl SOD 3 jstep
TIERC_REF=std $J tierc.jl jstep
echo JCORE DONE
