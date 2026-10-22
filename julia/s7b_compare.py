# S7b: compare coupled runs (halving test of the exchange interval). Usage: python3 -I s7b_compare.py <tagA> <tagB> [t_c list in yr]
# Reads lib/popiii/s6b_s7b_<tag>.csv (protostar series, absolute star time) + the 'coupling starts' line of the star stdout/log for the clock zero, and s7b_<tag>_env.csv (barrier records).
import sys, csv, re, os
LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib", "popiii")
YR = 3.15576e7
def star(tag):
    rows = list(csv.DictReader(open(os.path.join(LIB, f"s6b_s7b_{tag}.csv"))))
    t0 = None
    for fn in (f"/tmp/s7b_{tag}_star.out", os.path.join(LIB, f"s6b_s7b_{tag}.log")):
        if os.path.exists(fn):
            for l in open(fn):
                m = re.search(r"coupling starts at t = ([0-9.eE+-]+) s", l)
                if m and "star" in l: t0 = float(m.group(1)) / YR
            if t0 is not None: break
    return rows, t0
def interp(xs, ys, x):
    if x <= xs[0]: return ys[0]
    if x >= xs[-1]: return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x: lo = mid
        else: hi = mid
    w = (x - xs[lo]) / (xs[hi] - xs[lo]); return (1 - w) * ys[lo] + w * ys[hi]
def series(tag, key):
    rows, t0 = star(tag)
    t = [float(r["t_yr"]) - t0 for r in rows]; return t, [float(r[key]) for r in rows]
tags = sys.argv[1:3]; tcs = [float(x) for x in sys.argv[3:]] or [10, 20, 30, 40, 60, 80]
keys = [("M", "M"), ("R_Rsun", "R*"), ("Lstar_Lsun", "L*"), ("Teff", "Teff,ph"), ("Tc", "Tc")]
print("t_c [yr]   " + "   ".join(f"{lab:>9s}({a}/{b}) rel.diff" for _, lab in keys for a, b in [tags]))
for tc in tcs:
    out = []
    for key, lab in keys:
        v = []
        for tag in tags:
            t, y = series(tag, key); v.append(interp(t, y, tc))
        out.append(f"{lab}: {v[0]:.5g} / {v[1]:.5g}  {(v[1]/v[0]-1) if v[0] else float('nan'):+.2e}")
    print(f"{tc:7.1f}   " + "   ".join(out))
# envelope side: accreted mass and rate at the barriers
for tag in tags:
    rows = list(csv.DictReader(open(os.path.join(LIB, f"s7b_{tag}_env.csv"))))
    print(tag, "barriers:", len(rows), " last t_c = %.2f yr, M_acc = %.3f, Mdot = %.4e, M_pt = %.3f" % (float(rows[-1]["tc_yr"]), float(rows[-1]["Macc"]), float(rows[-1]["Mdot_win"]), float(rows[-1]["Mpt"])))
