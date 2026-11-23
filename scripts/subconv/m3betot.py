#!/usr/bin/env python3
"""M3b: total-energy drift of a run with ETOTMON=1 ETOTLOG=<run>/etot.log (lines: jstep kind dt imb thr G dE fin fout; erg/sr, one per converged stage of every attempt).
Accepted stages = the last 'tr' and the 'bdf' that follows it, per jstep. Reports sum(imb), sum(thr), |sum imb|/sum thr, median/max |G| (G = imb/thr per stage), and the
cumulative imbalance as a fraction of the total gas+radiation+kinetic(+turbulent) energy change budget sum|dE|.   Usage: m3betot.py <rundir> ..."""
import sys, statistics as st
for d in sys.argv[1:]:
    last = {}
    for ln in open(d + "/etot.log"):
        f = ln.split()
        if len(f) < 9: continue
        j = int(f[0]); kind = f[1]; v = [float(x) for x in f[2:]]
        if kind == "tr": last[j] = {"tr": v}
        elif kind == "bdf" and j in last: last[j]["bdf"] = v
    imb = thr = dE = 0.0; G = []; n = 0
    for j, s in sorted(last.items()):
        for k in ("tr", "bdf"):
            if k in s:
                dt, i, t, g, de = s[k][:5]; imb += i; thr += t; dE += abs(de); G.append(abs(g)); n += 1
    print("%-14s accepted stages %4d (%d steps): sum imb %+.4e  sum thr %.4e  |sum imb|/sum thr %.3e  sum|dE| %.4e  |sum imb|/sum|dE| %.3e  median|G| %.2e  max|G| %.2e" % (
        d.rstrip('/').split('/')[-1], n, len(last), imb, thr, abs(imb) / thr, dE, abs(imb) / dE, st.median(G), max(G)))
