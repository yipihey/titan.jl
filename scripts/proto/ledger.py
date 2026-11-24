#!/usr/bin/env python3
"""ledger.py <etot.log> <prefix> : ETOTMON ledger of a run (ETOTLOG lines 'jstep kind dt imb thr G dE fin fout'; accepted stages = last 'tr' and
the 'bdf' after it, per jstep; as scripts/subconv/m3betot.py). Keys <prefix>.imb_over_thr (= |sum imb|/sum thr), .max_G, .median_G, .stages, .steps."""
import sys, statistics as st
f, pre = sys.argv[1], sys.argv[2]; last = {}
try: lines = open(f).readlines()
except OSError: lines = []
for ln in lines:
    x = ln.split()
    if len(x) < 9: continue
    j = int(x[0]); kind = x[1]; v = [float(t) for t in x[2:]]
    if kind == "tr": last[j] = {"tr": v}
    elif kind == "bdf" and j in last: last[j]["bdf"] = v
imb = thr = 0.0; G = []
for j, s in last.items():
    for k in ("tr", "bdf"):
        if k in s: imb += s[k][1]; thr += s[k][2]; G.append(abs(s[k][3]))
if G:
    print(pre + ".imb_over_thr", abs(imb) / thr); print(pre + ".max_G", max(G)); print(pre + ".median_G", st.median(G))
print(pre + ".stages", len(G)); print(pre + ".steps", len(last))
