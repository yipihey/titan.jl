#!/usr/bin/env python3
"""gates_summary.py <tested-out> <ref-out|-> <gates e.g. 1,2,3,4,5>  : the summary table of scripts/dev/gates.sh.
Evaluates the five default merge gates of docs/POPIII_PROTO_GATE.md from the run directories (job label <name>_<steps>):
  1 NKCHECK=2 (AD vs FD J.v, 2-norm) grey and MGROUPS=4, zie1 S6_auto1 : worse than the documented value by > 15 % FAILS (grey 1.5e-10, MG 2.306e-7)
  2 bit-identity MGROUPS unset vs the reference commit: step lines (cmp.sh's filter) + decompressed tstate of the 40-step IC run
  3 late restarts zie1 S16, zie3 S8, 50 steps: completes, 50 steps, step lines + tstate identical to the reference commit's run
  4 ETOTMON ledger |sum imb|/sum thr over the same 50-step windows: must not exceed the reference commit's (+5 % noise margin)
  5 MG window zie3 fp3a1 40 steps: all 40 accepted, median dt >= 0.8 * 9e3 s (063452c: all accepted, ~9e3 s)
Exit status 1 if any row FAILs. The step-line filter and the tstate comparison are those of prod_env/cmp.sh; the ledger is ledger.py of scripts/proto (proto-gate)."""
import sys, os, re, csv, gzip, statistics as S

T, R, GATES = sys.argv[1], sys.argv[2], set(sys.argv[3].split(","))
NK_GREY, NK_MG, NK_SLACK = 1.5e-10, 2.306e-7, 1.15
MG_DT_REF, MG_DT_SLACK = 9e3, 0.8
STEPRE = re.compile(r"^ +[0-9]+ +[0-9.]+ +[0-9.e+]+ |^DTERR|^ *divergence")

def txt(d, lab, f="run.out"):
    p = os.path.join(d, lab, f)
    return open(p, errors="replace").read() if os.path.exists(p) else None

def exit_ok(d, lab):
    t = txt(d, lab); return t is not None and re.search(r"^===== EXIT 0 ", t, re.M) is not None

def wall(d, lab):
    t = txt(d, lab)
    m = re.findall(r"^===== EXIT \d+ wall ([\d.]+)", t or "", re.M)
    return float(m[-1]) if m else None

def steplines(d, lab):
    t = txt(d, lab); return None if t is None else [l for l in t.splitlines() if STEPRE.match(l)]

def tstate(d, lab, name="sv"):
    p = os.path.join(d, lab, "ck", "m2c_%s.tstate.gz" % name)
    if not os.path.exists(p): return None
    with gzip.open(p, "rb") as f: return f.read()

def perf(d, lab):
    p = os.path.join(d, lab, "perf.csv")
    if not os.path.exists(p): return None
    rows = list(csv.DictReader(open(p)))
    if len(rows) < 2: return None
    dts = [float(r["dt_s"]) for r in rows[1:]]
    return dict(n=len(rows) - 1, med=S.median(dts), ftr=int(rows[-1]["fail_tr"]) - int(rows[0]["fail_tr"]),
                fbd=int(rows[-1]["fail_bdf"]) - int(rows[0]["fail_bdf"]))

def ended(d, lab):
    """steps of a finished run: 'run ended: (N, :maxsteps)'"""
    m = re.findall(r"^run ended: \((\d+), :maxsteps\)", txt(d, lab) or "", re.M)
    return int(m[-1]) if m else None

def dterr(d, lab):
    """accepted steps, dt median and Newton-failure delta from the DTERR trace lines (the run needs PRINTEVERY=1)"""
    L = [l for l in (txt(d, lab) or "").splitlines() if l.startswith("DTERR jstep")]
    if not L: return None
    dts = [float(re.search(r": dt (\S+)", l).group(1)) for l in L]
    fl = [int(re.search(r" fail (\d+)", l).group(1)) for l in L]; rj = [int(re.search(r" rej (\d+)", l).group(1)) for l in L]
    return dict(n=len(L), med=S.median(dts), fail=fl[-1] - fl[0], rej=rj[-1] - rj[0])

def ledger(d, lab):
    p = os.path.join(d, lab, "etot.log")
    if not os.path.exists(p): return None
    last = {}
    for ln in open(p):
        x = ln.split()
        if len(x) < 9: continue
        j = int(x[0]); v = [float(t) for t in x[2:]]
        if x[1] == "tr": last[j] = {"tr": v}
        elif x[1] == "bdf" and j in last: last[j]["bdf"] = v
    imb = thr = 0.0; n = 0
    for s in last.values():
        for k in ("tr", "bdf"):
            if k in s: imb += s[k][1]; thr += s[k][2]; n += 1
    return abs(imb) / thr if thr > 0 and n else None

def nk2norm(d, lab):
    t = txt(d, lab) or ""
    m = re.search(r"NKCHECK AD J v vs FD \(eps x1\).*2-norm rel (\S+)", t)
    return float(m.group(1)) if m else None

def same(lab):
    if R == "-": return None, "no reference"
    a, b = steplines(T, lab), steplines(R, lab)
    if a is None or b is None or not a: return None, "run missing"
    sa, sb = tstate(T, lab), tstate(R, lab)
    if sa is None or sb is None: return None, "tstate missing"
    ok = a == b and sa == sb
    return ok, "%d step lines, tstate %d B %s" % (len(a), len(sa), "identical" if ok else ("STEPS DIFFER" if a != b else "TSTATE DIFFERS"))

rows = []
def add(g, name, value, ref, ok): rows.append((g, name, value, ref, ok))

if "1" in GATES:
    for nm, lab, ref in (("NKCHECK grey zie1 S6_auto1", "nkg_5", NK_GREY), ("NKCHECK MGROUPS=4 zie1 S6_auto1", "nkm_5", NK_MG)):
        v = nk2norm(T, lab)
        add(1, nm, "%.3e" % v if v is not None else "no result", "%.3e (+15%%)" % ref, None if v is None else v <= NK_SLACK * ref)
if "2" in GATES:
    ok, msg = same("ic_40"); add(2, "bit-identity IC 40 steps, MGROUPS off", msg, "ref commit", ok)
if "3" in GATES:
    for nm, lab in (("restart zie1 S16, 50 steps", "s16_50"), ("restart zie3 S8, 50 steps", "s8_50")):
        n = ended(T, lab); ok, msg = same(lab)
        done = exit_ok(T, lab) and n is not None and n >= 50
        add(3, nm, "%s; %s" % ("%d steps" % n if n is not None else "no result", msg), "50 steps, = ref", None if n is None else (done and ok is not False))
if "4" in GATES:
    for nm, lab in (("ETOTMON zie1 S16 window", "s16_50"), ("ETOTMON zie3 S8 window", "s8_50")):
        v, r = ledger(T, lab), (ledger(R, lab) if R != "-" else None)
        add(4, nm, "%.3e" % v if v is not None else "no result", "%.3e" % r if r is not None else "none",
            None if v is None else (True if r is None else v <= 1.05 * r + 1e-12))
if "5" in GATES:
    p = perf(T, "mgz3_40"); n5 = ended(T, "mgz3_40")
    if p is None: add(5, "MG window zie3 fp3a1, 40 steps", "no result", "40 acc, dt~9e3", None)
    else: add(5, "MG window zie3 fp3a1, 40 steps", "%s accepted, dt med %.3e s, fails tr/bdf %d/%d" % (n5, p["med"], p["ftr"], p["fbd"]), "40 acc, dt~9e3 s",
              exit_ok(T, "mgz3_40") and n5 == 40 and p["med"] >= MG_DT_SLACK * MG_DT_REF)

fail = 0
print("%-4s %-34s %-62s %-18s %s" % ("gate", "check", "value", "reference", ""))
for g, nm, v, r, ok in rows:
    st = "FAIL(no result)" if ok is None else ("PASS" if ok else "FAIL")
    if ok is not True: fail = 1
    print("%-4d %-34s %-62s %-18s %s" % (g, nm, v, r, st))
sys.exit(fail)
