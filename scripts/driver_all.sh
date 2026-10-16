#!/bin/bash
cd "$(dirname "$0")/../julia"
J="nice /usr/local/bin/julia --startup-file=no driver_all.jl"
$J BLST titan.in04; $J SDOV titan.in05; $J HEAT titan.in06
$J COOL titan.in07 HEAT        # COOL restarts from the reference HEAT dump (copied as COOL.dmp below)
$J RLST titan.in10
$J RLST_leg2 titan.in RLST     # leg 2 restarts from the reference RLST leg-1 dump
echo DRIVER_ALL DONE
