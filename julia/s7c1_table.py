#!/usr/bin/env python3
"""S7c-1: table of the two-way coupled run from the envelope record (s7b_<tag>_env.csv): python3 s7c1_table.py <tag> [Mpt ...]"""
import csv, sys
r = list(csv.DictReader(open(f"../lib/popiii/s7b_{sys.argv[1]}_env.csv"))) if False else list(csv.DictReader(open(f"/Users/tabel/Projects/titan.jl/.claude/worktrees/agent-a5ea8b87f9981e927/lib/popiii/s7b_{sys.argv[1]}_env.csv")))
Ms = [float(x) for x in sys.argv[2:]] or [10, 20, 30, 50, 70, 100, 130]
print("Mpt   t[yr]  Mdot[Msun/yr]  M*   R*[Rsun]  L*[Lsun]  Teff*[K]  Q*[1/s]   R_II/R*  F/(Mdot v)  edge[AU]")
for M in Ms:
    x = next((x for x in r if float(x['Mpt']) >= M), None)
    if x is None: continue
    f = lambda k: float(x[k])
    print(f"{f('Mpt'):6.1f} {f('tc_yr'):7.0f} {f('Mdot_win'):9.2e} {f('Mstar'):6.1f} {f('Rstar_Rsun'):7.1f} {f('Lstar'):9.2e} {f('Teffstar'):8.0f} {f('Qstar'):9.2e} {f('RII_over_R'):6.2f} {f('Fsum_over_Mdotv'):8.3f} {f('r_edge_AU'):8.0f}")
x = r[-1]; print("last:", x['j'], "t_c", x['tc_yr'], "Mpt", x['Mpt'])
