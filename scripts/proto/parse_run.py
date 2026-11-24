#!/usr/bin/env python3
"""parse_run.py <rundir> <prefix> : step-control / Newton / wall metrics of a pilot run directory (run.out, nwt.txt from NWTRACE, perf.csv, wall.txt).
Keys: <prefix>.ended_maxsteps (1 if 'run ended: (N, :maxsteps)'), .steps, .newton_fail_rej (DTERR reject ... Newton failure), .err_rej (DTERR reject ... err),
.nwt_attempts/.nwt_fail_frac/.nwt_it_mean/.nwt_it_median/.nwt_it_max/.nwt_gt10 (from nwt.txt), .wall_total_s, .wall_per_step (slope of perf.csv
wall_s vs step if >= 2 rows, else wall_total/steps), .trueres_refused (Newton convergences refused by the true-residual test, if printed)."""
import os, re, statistics as st, sys
d, pre = sys.argv[1], sys.argv[2]
out = []
txt = open(os.path.join(d, "run.out"), errors="replace").read() if os.path.exists(os.path.join(d, "run.out")) else ""
m = re.search(r"run ended: \((\d+), :(\w+)\)", txt)
steps = int(m.group(1)) if m else 0
out += [("ended_maxsteps", int(bool(m and m.group(2) == "maxsteps"))), ("steps", steps)]
out += [("newton_fail_rej", len(re.findall(r"DTERR reject[^\n]*Newton failure", txt))), ("err_rej", len(re.findall(r"DTERR reject[^\n]* err ", txt)))]
m = re.search(r"true residual[^\n]*refused[^\n]*?(\d+)", txt, re.I)
nw = os.path.join(d, "nwt.txt")
if os.path.exists(nw):
    ok, bad = [], 0
    for ln in open(nw, errors="replace"):
        if ln.startswith("NWT end"):
            q = re.search(r"ok (\w+) iter (\d+)", ln)
            if q:
                if q.group(1) == "true": ok.append(int(q.group(2)))
                else: bad += 1
    n = len(ok) + bad
    if n:
        out += [("nwt_attempts", n), ("nwt_fail_frac", bad / n)]
        if ok: out += [("nwt_it_mean", st.mean(ok)), ("nwt_it_median", st.median(ok)), ("nwt_it_max", max(ok)), ("nwt_gt10", sum(1 for i in ok if i > 10))]
w = float(open(os.path.join(d, "wall.txt")).read()) if os.path.exists(os.path.join(d, "wall.txt")) else None
if w is not None: out.append(("wall_total_s", w))
pf = os.path.join(d, "perf.csv"); rows = []
if os.path.exists(pf):
    L = open(pf).read().splitlines()
    if len(L) > 1:
        h = L[0].split(",")
        for ln in L[1:]:
            x = ln.split(",")
            if len(x) == len(h): rows.append((int(float(x[h.index("step")])), float(x[h.index("wall_s")])))
if len(rows) >= 2 and rows[-1][0] > rows[0][0]: out.append(("wall_per_step", (rows[-1][1] - rows[0][1]) / (rows[-1][0] - rows[0][0])))
elif w is not None and steps: out.append(("wall_per_step", w / steps))
for k, v in out: print("%s.%s %s" % (pre, k, repr(v) if isinstance(v, float) else v))
