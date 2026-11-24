#!/usr/bin/env python3
"""finalize.py <tier-dir> <tier> : collect every '*.kv' file ('key value' lines) of a tier directory into metrics.json (flat key -> number; the last value of a key wins,
non-numeric values are dropped with a warning) and write summary.txt (human-readable, grouped by the key's first two dot-fields).
Also merges <tier-dir>/notes.txt (SKIP / info lines written by the tiers) into the summary."""
import glob, json, math, os, sys
d, tier = sys.argv[1], sys.argv[2]
M = {}
for f in sorted(glob.glob(os.path.join(d, "**", "*.kv"), recursive=True)):
    for ln in open(f):
        p = ln.split()
        if len(p) != 2: continue
        try: v = float(p[1])
        except ValueError: print("finalize: non-numeric", ln.strip(), file=sys.stderr); continue
        if math.isnan(v) or math.isinf(v): M[p[0]] = v
        else: M[p[0]] = int(v) if v == int(v) and abs(v) < 1e15 and "." not in p[1] and "e" not in p[1].lower() else v
clean = {k: (None if isinstance(v, float) and (math.isnan(v) or math.isinf(v)) else v) for k, v in M.items()}
json.dump(clean, open(os.path.join(d, "metrics.json"), "w"), indent=1, sort_keys=True)
L = ["tier %s: %d metrics" % (tier, len(M)), ""]
nf = os.path.join(d, "notes.txt")
if os.path.exists(nf): L += open(nf).read().rstrip().splitlines() + [""]
grp = {}
for k in sorted(M): grp.setdefault(".".join(k.split(".")[:2]), []).append(k)
for g in sorted(grp):
    L.append("[%s]" % g)
    for k in grp[g]: L.append("  %-60s %s" % (k, "%.6g" % M[k] if isinstance(M[k], float) else M[k]))
open(os.path.join(d, "summary.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L[:3 if len(L) > 3 else None]))
