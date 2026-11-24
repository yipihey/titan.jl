#!/bin/bash
# M6: worst |J v(AD) - J v(FD)|/row scale per equation from the "NKCHECK AD row eq <n> zone <z>" lines of a run.out (NKEQ listing), zones Z1..Z2.
# Usage: lresp_nksummary.sh <run.out> <Z1> <Z2> <MG|OFF>      names: eq 4 u, 5 T, 8 z_H, 9 x_H2; MG: 10 E_g1, 11 F_g1, 12 E_g2, … 17 F_g4, 18 w; OFF: 10 w
f=$1; z1=$2; z2=$3; mode=$4
grep -a "NKCHECK AD row" "$f" | awk -v z1=$z1 -v z2=$z2 -v mode=$mode '
BEGIN { nm[4]="u"; nm[5]="T"; nm[8]="z_H"; nm[9]="x_H2"
  if (mode=="MG") { nm[10]="E_g1"; nm[11]="F_g1"; nm[12]="E_g2"; nm[13]="F_g2"; nm[14]="E_g3"; nm[15]="F_g3"; nm[16]="E_g4"; nm[17]="F_g4"; nm[18]="w" } else { nm[10]="w" } }
{ eq=$5; z=$7+0; rel=$0; sub(/.* rel /,"",rel); split(rel,a," "); r=a[1]
  if (z>=z1 && z<=z2) { if (r!=r || r=="NaN") { nan[eq]++ } else if (!(eq in w) || r+0>w[eq]) { w[eq]=r+0; wz[eq]=z } } }
END { n=asort_idx=0; for (e=4;e<=18;e++) if ((e in w) || (e in nan)) printf("%s %.1e(z%d)%s | ", nm[e], (e in w)?w[e]:-1, wz[e], (e in nan)?"+"nan[e]"NaN":"") ; print "" }'
