#!/usr/bin/env python3
"""ATTR-1 summary: per run (lib/popiii/s7b_<tag>_env.csv) the final M_pt, the time of the last sink event, the star mass and the end time,
plus the time at which M_pt first reaches a few masses.  Usage: attr1_summary.py <csv-or-tag> ..."""
import csv, os, sys
L = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib", "popiii")
MS = [30, 52, 78, 102, 108, 114, 121, 128, 135, 143]
def summ(arg):
    f = arg if os.path.isfile(arg) else os.path.join(L, "s7b_%s_env.csv" % arg)
    r = list(csv.DictReader(open(f)))
    mpt = [float(x["Mpt"]) for x in r]; tc = [float(x["tc_yr"]) for x in r]
    last = max(i for i in range(len(r)) if i == 0 or mpt[i] > mpt[i - 1] + 1e-6)
    out = "%-34s M_pt = %7.2f  t_last_event = %7.0f yr  M* = %7.2f  t_end = %7.0f yr  Mdot_end = %.1e  rows %d" % (os.path.basename(f)[4:-8], mpt[-1], tc[last], float(r[-1]["Mstar"]), tc[-1], float(r[-1]["Mdot_win"]), len(r))
    tt = []
    for m in MS:
        t = next((tc[i] for i in range(len(r)) if mpt[i] >= m - 0.5), None)
        tt.append("%.0f" % t if t is not None else "-")
    return out + "\n    t(M_pt >= " + ",".join(map(str, MS)) + ") = " + " ".join(tt)
for a in [x for x in sys.argv[1:] if not x.startswith("--")]:
    print(summ(a))

# ---- table of the ATTR-1 runs: python3 attr1_summary.py --table > ../lib/popiii/attr1_summary.csv
RUNS = [  # tag, lya, snapshot, THICKCONT/TEQMAX, STARVECUT+MDOTSTOP, other, note
    ("abn02z19_fb", "mc", "E2E-A c2", "3/1e30", "3 + 1e-3", "", "FULL-1"),
    ("attr1_fm", "mc", "E2E-A c2", "3/1e30", "3 + 1e-3", "", "FULL-1 repeated (bit-identical)"),
    ("attr1_fm2", "mc", "E2E-A c2", "3/1e30", "3 + 1e-3", "S7B_MSYNC 2.9", "perturbation of the exchange spacing"),
    ("attr1_fm3", "mc", "E2E-A c2", "3/1e30", "3 + 1e-3", "3x packets", "S8A_NPK_* x3"),
    ("attr1_fma", "mc", "core2", "3/1e30", "3 + 1e-3", "", "revert (a)"),
    ("attr1_fma5", "mc", "core2", "3/1e30", "3 + 1e-3", "S8A_NEVERY 5", "revert (a), MC force held 5 barriers instead of 25"),
    ("attr1_fmb", "mc", "E2E-A c2", "off", "3 + 1e-3", "", "revert (b)"),
    ("attr1_sm", "mc", "core2", "off", "off", "", "S8a-like with MC and EDDCAP from the start (stopped at 13 kyr)"),
    ("d4", "mc", "core2", "off", "off", "EDDCAP after restart at 100", "S8a d4"),
    ("attr1_f0", "none", "E2E-A c2", "3/1e30", "3 + 1e-3", "star frozen at its failure (145.5)", "FULL-1 without Lya"),
    ("attr1_f0d", "none", "E2E-A c2", "3/1e30", "off", "star frozen (145.5)", "revert (d) (stopped at 14.3 kyr)"),
    ("attr1_f0a", "none", "core2", "3/1e30", "3 + 1e-3", "star alive", "revert (a)"),
    ("attr1_f0b", "none", "E2E-A c2", "off", "3 + 1e-3", "star alive", "revert (b)"),
    ("attr1_s0b", "none", "core2", "3/1e30", "off", "star alive", "S8a-like + (b) (stopped at TCMAX 22 kyr, still accreting)"),
    ("attr1_s0", "none", "core2", "off", "off", "star frozen (144.8)", "S8a d5 reproduced (stopped at 17 kyr)"),
    ("d5", "none", "core2", "off", "off", "star frozen (144.8)", "S8a d5"),
]
if "--table" in sys.argv:
    print("run,lya,snapshot,thickcont_teqmax,starvecut_mdotstop,other,M_pt,t_last_event_yr,M_star,t_end_yr,note")
    for tag, lya, snap, th, cut, oth, note in RUNS:
        f = os.path.join(L, "s7b_%s_env.csv" % ("abn02z19_fb" if tag == "abn02z19_fb" else tag.replace("attr1_", "attr1_")))
        r = list(csv.DictReader(open(f))); mpt = [float(x["Mpt"]) for x in r]; tc = [float(x["tc_yr"]) for x in r]
        last = max(i for i in range(len(r)) if i == 0 or mpt[i] > mpt[i - 1] + 1e-6)
        print("%s,%s,%s,%s,%s,%s,%.2f,%.0f,%.2f,%.0f,%s" % (tag, lya, snap, th, cut, oth, mpt[-1], tc[last], float(r[-1]["Mstar"]), tc[-1], note))
    sys.exit(0)
