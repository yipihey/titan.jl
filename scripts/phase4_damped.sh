#!/bin/bash
# Phase 4 item 5: damped Newton vs TITAN retry, both on the full Julia stack, fixes ON; tier (c) vs std
cd "$(dirname "$0")/../julia"
J="nice julia --startup-file=no --project=env"
$J phase3_run.jl SUBC 8 std; $J phase3_run.jl SUPC 9 std; TITAN_SET="jstepe=3000" $J phase3_run.jl SOD 3 std
for p in "SOD 3" "SUBC 8" "SUPC 9"; do TITAN_SET="jstepe=3000" $J newton_run.jl $p titan,damped 0; done
TIERC_REF=std $J tierc.jl damped_nr
TIERC_REF=std $J tierc.jl jfast_nr
echo DAMPED DONE
