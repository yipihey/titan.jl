#!/usr/bin/env python3
"""stage2c_steady.py <rundir> [--table]: dt trajectory of a segmented run from the DTERR lines (one per 10 steps) of <rundir>/run_s*.out.
Exit 0 = steady (last three 100-step windows have median dt within a factor 1.5 and 1e3 s < dt, or dt >= 5e6 s), 1 = not yet. Prints the last step."""
import re, sys, glob, statistics as st
rx = re.compile(r"DTERR jstep (\d+) t (\S+) yr: dt (\S+) err (\S+) .*?rej (\d+) fail (\d+)")
rows = {}
for f in sorted(glob.glob(sys.argv[1] + "/run_s*.out")):
    for l in open(f, errors="replace"):
        m = rx.search(l)
        if m: rows[int(m[1])] = (float(m[3]), int(m[5]), int(m[6]))
ks = sorted(rows)
if "--table" in sys.argv:
    for k in ks: print(k, *rows[k])
if not ks: print("nolines"); sys.exit(1)
dts = [rows[k][0] for k in ks]; print("last step", ks[-1], "n", len(ks), "dt", dts[-1])
w = [st.median(dts[i:i+10]) for i in range(len(dts) - 30, len(dts), 10)] if len(dts) >= 30 else None
if w: print("last 3 window medians", w)
steady = w is not None and (max(w) / min(w) < 1.5 and min(w) > 1e3 or min(w) >= 5e6)
sys.exit(0 if steady else 1)
