#!/usr/bin/env python3
"""S7c-1: ratio of two protostar tracks (s6b csv) at equal M: python3 s7c1_compare.py new.csv old.csv [M1 M2 ...]. No dependencies."""
import sys, csv, bisect

def load(p):
    rows = list(csv.DictReader(open(p)))
    f = lambda k: [float(r[k]) for r in rows]
    return dict(t=f("t_yr"), M=f("M"), R=f("R_Rsun"), L=f("Lstar_Lsun"), Tc=f("Tc"), fD=f("fD"))

def at(d, key, M):
    m, y = d["M"], d[key]
    i = bisect.bisect_left(m, M)
    if i <= 0: return y[0]
    if i >= len(m): return y[-1]
    w = (M - m[i-1]) / (m[i] - m[i-1]) if m[i] > m[i-1] else 0.0
    return (1 - w) * y[i-1] + w * y[i]

def first(d, key, thr, below=False):
    for M, v in zip(d["M"], d[key]):
        if (v < thr) if below else (v > thr): return M
    return float("nan")

a, b = load(sys.argv[1]), load(sys.argv[2])
Ms = [float(x) for x in sys.argv[3:]] or [0.5, 1, 3, 5, 8, 10, 15, 20, 30, 40, 50]
print("M      R_new/R_old  L_new/L_old  Tc_new/Tc_old   (R_new, R_old)")
for M in Ms:
    if M > min(a["M"][-1], b["M"][-1]): continue
    print(f"{M:6.1f} {at(a,'R',M)/at(b,'R',M):10.3f} {at(a,'L',M)/at(b,'L',M):11.3f} {at(a,'Tc',M)/at(b,'Tc',M):12.3f}   ({at(a,'R',M):.1f}, {at(b,'R',M):.1f})")
for nm, d in (("new", a), ("old", b)):
    R = d["R"]; i = R.index(max(R))
    print(f"{nm}: R peak {R[i]:.1f} Rsun at M = {d['M'][i]:.2f}; fD<0.99/0.5/0.1 at M = {first(d,'fD',0.99,True):.2f}/{first(d,'fD',0.5,True):.2f}/{first(d,'fD',0.1,True):.2f}; "
          f"L*>1e2/1e3/1e4 at M = {first(d,'L',1e2):.2f}/{first(d,'L',1e3):.2f}/{first(d,'L',1e4):.2f}; last M = {d['M'][-1]:.2f}, R = {R[-1]:.2f}")
