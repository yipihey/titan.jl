#!/bin/bash
cd "$(dirname "$0")/../julia"
J="nice /usr/local/bin/julia --startup-file=no --project=env phase3_run.jl"
for P in "SUBC 8" "SUPC 9" "SOD 3"; do
  $J $P jsolve; $J $P jall; $J $P jall tierb
done
echo P3 DONE
