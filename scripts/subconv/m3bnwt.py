#!/usr/bin/env python3
"""M3b: Newton statistics of runs with NWTRACE (lines 'NWT head att N jstep J ...', 'NWT att N it I ...', 'NWT end att N ok true|false iter K').
Usage: m3bnwt.py <rundir> [<rundir> ...]   -> per run: attempts, converged/failed, mean/median/max iterations of converged attempts, share with > 10 iterations,
and the steps' count of attempts. Also the accepted dt statistics from the run's own `snaps` are in m3bsnap.jl."""
import sys, re, statistics as st
for d in sys.argv[1:]:
    ok_it = []; fail_it = []; n = 0; jsteps = set()
    for line in open(d + "/nwt.txt"):
        if line.startswith("NWT head"):
            m = re.search(r"jstep (\d+)", line); jsteps.add(int(m.group(1)))
        elif line.startswith("NWT end"):
            m = re.search(r"ok (\w+) iter (\d+)", line); n += 1
            (ok_it if m.group(1) == "true" else fail_it).append(int(m.group(2)))
    print("%-14s attempts %4d (steps %3d): converged %4d, mean it %.2f, median %.1f, max %d, >10 it: %d | failed %3d (%.1f%%)" % (
        d.rstrip('/').split('/')[-1], n, len(jsteps), len(ok_it), st.mean(ok_it) if ok_it else 0, st.median(ok_it) if ok_it else 0, max(ok_it) if ok_it else 0,
        sum(1 for i in ok_it if i > 10), len(fail_it), 100.0 * len(fail_it) / max(n, 1)))
