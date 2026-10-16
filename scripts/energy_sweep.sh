#!/bin/bash
# Phase-3 energy convergence sweep (sequential, nice'd). Results: lib/p3_runs/energy_*.log
cd "$(dirname "$0")/../julia"
J="nice /usr/local/bin/julia --startup-file=no energy_ledger.jl"
for P in "SUBC 8" "SUPC 9" "SOD 3"; do
  $J $P lgrid3 lgrid=3
  $J $P stol05 stol=0.05
  $J $P stol025 stol=0.025
  $J $P conv7 conv=1e-7
  TITAN_NCELL=200 $J $P n200
done
$J SUBC 8 thet075 thet=0.75
$J SOD 3 thet1 thet=1.0
TITAN_NCELL=200 $J SUBC 8 n200stol05 stol=0.05
echo SWEEP DONE
