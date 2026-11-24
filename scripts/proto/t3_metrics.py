#!/usr/bin/env python3
"""t3_metrics.py <export dir of the new case> <prefix> <reference export dir | -> <smin of the reference> : accretion-history observables at fixed M* (1, 3, 10, 30 Msun)
from a runrec export (track.csv, frames.json/f32; scripts/subconv/m4extra.jl) with the series() of m4compare.py. Keys <prefix>.<qty>.M<x> and, with a reference,
<prefix>.<qty>.M<x>.reldiff = (new-ref)/|ref| at the same M*, plus <prefix>.M_final / .frames."""
import importlib.util, os, sys
import numpy as np
spec = importlib.util.spec_from_file_location("m4c", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "subconv", "m4compare.py"))
m4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m4)
new, pre, ref, smin = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
QTY = ["Rstar", "Mdot1e4", "Mdot_sh", "umin", "u_pk", "s_ps", "T_ps", "Tc", "Teff", "Lph", "cv_Fmax", "cv_PtPps"]
TARGETS = [1.0, 3.0, 10.0, 30.0]
def ser(path, smin_):
    d, c = os.path.split(path.rstrip("/"))
    cc = m4.load(d, c, smin_)
    if cc is None: return None
    cc["sel"] = np.where(cc["tr"]["step"] >= smin_)[0]
    return m4.series(cc) if len(cc["sel"]) else None
S = ser(new, 0)
if S is None: print(pre + ".frames 0"); sys.exit(0)
print(pre + ".frames", len(S["M"])); print(pre + ".M_final", float(np.nanmax(S["M"])))
R = ser(ref, smin) if ref != "-" and os.path.isdir(ref) else None
for q in QTY:
    v = m4.at_x(S["M"], S[q], TARGETS)
    r = m4.at_x(R["M"], R[q], TARGETS) if R is not None else None
    for i, t in enumerate(TARGETS):
        if np.isfinite(v[i]): print("%s.%s.M%g %.8g" % (pre, q, t, v[i]))
        if r is not None and np.isfinite(v[i]) and np.isfinite(r[i]) and r[i] != 0: print("%s.%s.M%g.reldiff %.6g" % (pre, q, t, (v[i] - r[i]) / abs(r[i])))
