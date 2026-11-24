#!/usr/bin/env python3
"""nkparse.py <run.out> <prefix> [grey|mg] : metrics from the NKCHECK listing of a run (NewtonKrylov.check!).
Keys <prefix>.<eq>.all / .zones = per-equation max relative AD/FD-vs-band J.v mismatch (all zones / the NKZ range),
<prefix>.max, .median, .rows_gt1e-4, .rows_gt1e-2 (global), <prefix>.ad_resid (AD vs Float64 residual), <prefix>.nan_rows (NaN entries), <prefix>.found (1/0).
mg mode: equation codes 10..17 are the group rows E_g1,F_g1,..., 18 is w (docs/POPIII_CONVECTION.md M3e)."""
import re, sys
f, pre = sys.argv[1], sys.argv[2]; mode = sys.argv[3] if len(sys.argv) > 3 else "grey"
MG = {"10": "E_g1", "11": "F_g1", "12": "E_g2", "13": "F_g2", "14": "E_g3", "15": "F_g3", "16": "E_g4", "17": "F_g4", "18": "w"}
OFF = {"10": "w"}
out = []; found = 0; nan = 0
num = lambda s: float(s)
for ln in open(f, errors="replace"):
    if not ln.startswith("NKCHECK"): continue
    m = re.match(r"NKCHECK eq (\S+)\s*: max rel (\S+) \(zone (\d+)\) \| zones (\d+)\.\.(\d+): max rel (\S+)", ln)
    if m:
        nm = m.group(1).replace("+", "p"); nm = (MG if mode == "mg" else OFF).get(nm, nm if not nm.isdigit() else "eq" + nm)
        for tag, v in (("all", m.group(2)), ("zones", m.group(6))):
            x = num(v)
            if x != x: nan += 1; continue
            out.append(("%s.%s.%s" % (pre, nm, tag), x))
        found = 1; continue
    m = re.match(r"NKCHECK ne .*max ([\d.eE+-]+), median ([\d.eE+-]+), rows > 1e-4: (\d+), > 1e-2: (\d+)", ln)
    if m:
        out += [(pre + ".max", num(m.group(1))), (pre + ".median", num(m.group(2))), (pre + ".rows_gt1e-4", int(m.group(3))), (pre + ".rows_gt1e-2", int(m.group(4)))]; continue
    m = re.match(r"NKCHECK AD residual vs Float64 residual: max \|diff\| / max \|R\| = (\S+)", ln)
    if m: out.append((pre + ".ad_resid", num(m.group(1))))
out += [(pre + ".found", found), (pre + ".nan_rows", nan)]
for k, v in out: print(k, repr(v) if isinstance(v, float) else v)
