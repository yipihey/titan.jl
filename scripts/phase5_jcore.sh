#!/bin/bash
# Phase 5: full runs with Julia core routines (jcore) vs fixed all-Fortran reference (std); tier (c)
cd "$(dirname "$0")/../julia"
J="nice julia --startup-file=no --project=env"
$J phase3_run.jl SUBC 8 jcore; $J phase3_run.jl SUPC 9 jcore; TITAN_SET="jstepe=3000" $J phase3_run.jl SOD 3 jcore
TIERC_REF=std $J tierc.jl jcore
echo JCORE DONE
