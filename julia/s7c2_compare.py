"""S7c-2: time at which M_pt first reaches given masses in the one-way runs (S7b F0/F1 and S7c-2 c2L/c2A/c2LA), final mass, halt time.
   python3 s7c2_compare.py [tags ...]   (reads lib/popiii/s7b_<tag>_env.csv)"""
import csv, sys, os
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib", "popiii")
tags = sys.argv[1:] or ["F0", "F1", "c2L", "c2A", "c2LA"]
masses = [10, 20, 30, 50, 70, 90, 100, 120, 140, 160]
print("tag   " + " ".join("%8s" % ("M=%d" % m) for m in masses) + "   M_final  t_last[yr]  t_halt[yr]  r_edge_last[AU]")
for t in tags:
    f = os.path.join(D, "s7b_%s_env.csv" % t)
    if not os.path.exists(f):
        continue
    r = list(csv.DictReader(open(f)))
    out = []
    for m in masses:
        x = next((float(q["tc_yr"]) for q in r if float(q["Mpt"]) >= m), None)
        out.append("%8s" % ("-" if x is None else "%.0f" % x))
    Mf = float(r[-1]["Mpt"])
    th = next((float(q["tc_yr"]) for q in r if float(q["Mpt"]) >= Mf - 1e-9), float("nan"))
    print("%-5s %s   %7.2f  %10.0f  %10.0f  %10.0f" % (t, " ".join(out), Mf, float(r[-1]["tc_yr"]), th, float(r[-1]["r_edge_AU"])))
print()
for t in tags:
    f = os.path.join(D, "s7b_%s_env.csv" % t)
    if not os.path.exists(f):
        continue
    r = list(csv.DictReader(open(f)))
    print(t, "Mdot_win [Msun/yr] at M_pt = 10, 30, 50, 70, 90:", " ".join("%.2e" % next((float(q["Mdot_win"]) for q in r if float(q["Mpt"]) >= m), float("nan")) for m in (10, 30, 50, 70, 90)))
