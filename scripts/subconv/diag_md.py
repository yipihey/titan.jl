#!/usr/bin/env python3
"""Markdown tables of the layer interfaces/zones from a diag_layer.jl log: diag_md.py <log> -> stdout (layer rows only)."""
import sys, re
L = open(sys.argv[1]).read().split('\n')
hd = L[0]
sec = 0; rows = {0: [], 1: [], 2: []}
for l in L[1:]:
    if l.startswith('# energy'): sec = 1; continue
    if l.startswith('# time scales'): sec = 2; continue
    if l.startswith('#') or not l.strip() or l.strip().split()[0] in ('j', 'k'): continue
    rows[sec].append(l)
print(hd[:200]); print()
def g(x): return x
def tab(cols, rs, pick):
    print('| ' + ' | '.join(cols) + ' |'); print('|' + '---|' * len(cols))
    for l in rs:
        f = l.split(); mark = l.rstrip().endswith('*') or '<-shock' in l
        if not mark: continue
        f = [x for x in f if x not in ('*', '<-shock')]
        print('| ' + ' | '.join(f[i] for i in pick) + (' (shock)' if '<-shock' in l else '') + ' |')
# interface table indices: j r/R* T rho Mach br N2_eq N2_S N2_L nab nab_rad nab_L nad_fr nad_eq w_tau w_inf w_th kap Fpot/F Fconv/F dr/HP dlnT/dr L
tab(['j', 'r/R*', 'T[K]', 'Mach', 'N2_eq', 'N2_L', 'nab', 'nab_rad', 'nad_fr', 'nad_eq', 'w_tau', 'w_inf', 'w_th', 'kap', 'Fpot/F', 'dr/HP'], rows[0], [0, 1, 2, 4, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 20])
print()
# energy rows: k r/R* T rho mu net av kin adv dfe rad conv chem src res/net |rad|/|tot| t_rad t_settle dr/HP
tab(['k', 'r/R*', 'T', 'sum ds/dt', 'av', 'adv', 'rad', 'chem', 'res/sum', 't_rad[yr]', 't_settle[yr]'], rows[1], [0, 1, 2, 6, 7, 9, 11, 13, 15, 17, 18])
print()
# time: j t_settle t_N t_chem U sup_fr v_fr t_turn_fr Fpot sup_eq v_eq t_turn_eq Feq_lb tau
print('| j | t_settle | t_N | t_chem[s] | sup_fr | t_turn_fr | Fpot/F | sup_eq | t_turn_eq | Feq_lb/F |'); print('|' + '---|' * 10)
lay = [l.split()[0] for l in rows[0] if l.rstrip().endswith('*')]
for l in rows[2]:
    f = l.split()
    if f[0] in lay: print('| ' + ' | '.join(f[i] for i in [0, 1, 2, 3, 5, 7, 8, 9, 11, 12]) + ' |')
