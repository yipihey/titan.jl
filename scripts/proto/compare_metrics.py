#!/usr/bin/env python3
"""compare_metrics.py new.json ref.json tol.json [--only PREFIX] : PASS/WARN/FAIL per key; exit 1 if any FAIL.
ref.json may be '-' (only the absolute rules, 'upper'/'exact0', apply). tol.json: {"rules": [{"pat": glob, "mode": ..., "tol": x, "warn": y, "atol": a, "note": ".."}, ...]}
First matching rule per key wins (fnmatch on the full key). Modes:
  absupper : |new| <= tol (a signed relative difference to the production reference).
  upper : new <= tol (absolute bound; WARN above warn). Needs no reference.      lower : new >= tol.
  zero  : new must be 0 (counts of failed tests, exit codes, differing fields). Needs no reference.
  rel   : |new-ref| <= warn*|ref| PASS, <= tol*|ref| WARN, else FAIL (warn defaults tol/2; atol is an absolute floor on the allowed difference).
  rise  : one-sided: new <= ref*(1+warn) + atol PASS, <= ref*(1+tol)+atol WARN, else FAIL (cost / count metrics: only an increase is bad).
  info  : never fails (printed as INFO).
A key in new with no rule is INFO; a ref key missing in new is WARN 'missing' (FAIL for keys matching a rule marked "required": true)."""
import fnmatch, json, math, sys
a = [x for x in sys.argv[1:] if not x.startswith("--")]
only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else ""
if only: a.remove(only)
new = json.load(open(a[0])); ref = {} if a[1] == "-" else json.load(open(a[1])); rules = json.load(open(a[2]))["rules"]
def rule(k):
    for r in rules:
        if fnmatch.fnmatchcase(k, r["pat"]): return r
cnt = {"PASS": 0, "WARN": 0, "FAIL": 0, "INFO": 0}
def show(st, k, msg): cnt[st] += 1; print("%-5s %-58s %s" % (st, k, msg))
for k in sorted(set(new) | set(ref)):
    if only and not k.startswith(only): continue
    r = rule(k); v = new.get(k); w = ref.get(k)
    if k not in new:
        if r and r.get("required") and r["mode"] not in ("upper", "zero", "lower"): show("FAIL", k, "missing in new (ref %s)" % w)
        elif r and r["mode"] != "info": show("WARN", k, "missing in new (ref %s)" % w)
        continue
    if r is None or r["mode"] == "info": show("INFO", k, "%s%s" % (v, "" if w is None else "  (ref %s)" % w)); continue
    m = r["mode"]; tol = r.get("tol", 0.0); wn = r.get("warn"); atol = r.get("atol", 0.0)
    if v is None or (isinstance(v, float) and math.isnan(v)): show("FAIL", k, "NaN"); continue
    if m == "absupper": show("FAIL" if abs(v) > tol else "WARN" if abs(v) > (tol / 2 if wn is None else wn) else "PASS", k, "%.4g (|.| <= %.3g)" % (v, tol))
    elif m == "upper":
        wn = tol if wn is None else wn
        show("FAIL" if v > tol else "WARN" if v > wn else "PASS", k, "%.4g (<= %.3g)" % (v, tol))
    elif m == "lower": show("FAIL" if v < tol else "PASS", k, "%.4g (>= %.3g)" % (v, tol))
    elif m == "zero": show("PASS" if v == 0 else "FAIL", k, "%.4g (must be 0)" % v)
    elif w is None: show("INFO", k, "%.6g (no reference)" % v)
    elif m == "rel":
        wn = tol / 2 if wn is None else wn; dd = abs(v - w); sc = abs(w)
        st = "PASS" if dd <= wn * sc + atol else "WARN" if dd <= tol * sc + atol else "FAIL"
        show(st, k, "%.6g vs ref %.6g (rel %.2e, tol %.2g)" % (v, w, dd / sc if sc else float("inf"), tol))
    elif m == "rise":
        wn = tol / 2 if wn is None else wn
        st = "PASS" if v <= w * (1 + wn) + atol else "WARN" if v <= w * (1 + tol) + atol else "FAIL"
        show(st, k, "%.6g vs ref %.6g (%+.1f%%, limit +%.0f%%)" % (v, w, 100 * (v / w - 1) if w else float("inf"), 100 * tol))
    else: show("FAIL", k, "unknown mode %s" % m)
print("== %d PASS, %d WARN, %d FAIL, %d INFO" % (cnt["PASS"], cnt["WARN"], cnt["FAIL"], cnt["INFO"]))
sys.exit(1 if cnt["FAIL"] else 0)
