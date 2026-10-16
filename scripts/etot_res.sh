#!/bin/bash
cd "$(dirname "$0")/../julia"
export TITAN_BITWISE=1 TITAN_SET="jstepe=6000"
J="nice /usr/local/bin/julia --startup-file=no --project=env etotal_run.jl"
for n in 200 400; do for f in std etot; do TITAN_NCELL=$n RUNTAG=_n$n $J SUBC 8 $f; done; done
TITAN_NCELL=200 RUNTAG=_n200 $J SOD 3 std; TITAN_NCELL=200 RUNTAG=_n200 $J SOD 3 etot
echo RES DONE
