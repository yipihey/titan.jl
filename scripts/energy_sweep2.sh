#!/bin/bash
# re-run step-limited sweep cases with jstepe raised so all reach tmax
cd "$(dirname "$0")/../julia"
J="nice /usr/local/bin/julia --startup-file=no energy_ledger.jl"
$J SUBC 8 stol05 stol=0.05 jstepe=20000
$J SUBC 8 stol025 stol=0.025 jstepe=20000
$J SUPC 9 stol05 stol=0.05 jstepe=20000
$J SUPC 9 stol025 stol=0.025 jstepe=20000
$J SOD 3 stol05 stol=0.05 jstepe=20000
$J SOD 3 stol025 stol=0.025 jstepe=20000
$J SOD 3 lgrid3 lgrid=3 jstepe=20000
TITAN_NCELL=200 $J SOD 3 n200 jstepe=20000
$J SOD 3 thet1 thet=1.0 jstepe=20000
$J SOD 3 conv7 conv=1e-7 jstepe=20000
TITAN_NCELL=200 $J SUBC 8 n200stol05 stol=0.05 jstepe=20000
TITAN_NCELL=300 $J SUBC 8 n300 jstepe=20000
TITAN_NCELL=300 $J SUPC 9 n300 jstepe=20000
$J SUBC 8 lgrid3 lgrid=3 jstepe=20000
echo SWEEP2 DONE
