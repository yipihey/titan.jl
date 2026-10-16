#!/bin/bash
cd "$(dirname "$0")/../julia"
J="nice /usr/local/bin/julia --startup-file=no legs.jl"
for M in orig fixed; do
  $J SOD titan.in03 $M; $J BLST titan.in04 $M; $J SDOV titan.in05 $M; $J HEAT titan.in06 $M
  $J COOL titan.in07 $M HEAT; $J SUBC titan.in08 $M; $J SUPC titan.in09 $M
  $J RLST titan.in10 $M; $J RLST_leg2 titan.in $M RLST
done
echo LEGS DONE
