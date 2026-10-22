# S7b: table of the coupled run (star csv + envelope barrier csv). Usage: python3 -I s7b_summary.py <tag> [M list]
import sys, csv, os, re
LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib", "popiii")
tag = sys.argv[1]; Ms = [float(x) for x in sys.argv[2:]] or [4, 10, 20, 30, 40, 60, 80, 100, 120, 150, 180, 200, 230, 250]
env = list(csv.DictReader(open(os.path.join(LIB, f"s7b_{tag}_env.csv"))))
star = list(csv.DictReader(open(os.path.join(LIB, f"s6b_s7b_{tag}.csv")))) if os.path.exists(os.path.join(LIB, f"s6b_s7b_{tag}.csv")) else []
f = lambda r, k: float(r[k])
print(f"{'Mstar':>7s} {'t_c[yr]':>9s} {'Mdot':>9s} {'R*':>7s} {'L*':>9s} {'Linf':>9s} {'Tph':>6s} {'Teff*':>7s} {'Q*':>9s} {'Qesc,in':>9s} {'RII/R*':>7s} {'Q/Qcrit':>9s} {'kLWmax':>9s} {'a/g max':>9s} {'F/Mdot v':>9s} {'xe max':>8s} {'lit':>4s}")
last = -1
for r in env:
    M = f(r, "Mstar")
    for m in Ms:
        if M >= m and last < m:
            last = m
            print(f"{M:7.2f} {f(r,'tc_yr'):9.1f} {f(r,'Mdot_win'):9.2e} {f(r,'Rstar_Rsun'):7.1f} {f(r,'Lstar'):9.2e} {f(r,'Linf'):9.2e} {f(r,'Tph'):6.0f} {f(r,'Teffstar'):7.0f} {f(r,'Qstar'):9.2e} {f(r,'Qesc_in'):9.2e} {f(r,'RII_over_R'):7.2f} {f(r,'Q_over_Qcrit'):9.2e} {f(r,'kLW_max'):9.2e} {f(r,'arad_over_g_max'):9.2e} {f(r,'Fsum_over_Mdotv'):9.2e} {f(r,'xe_max'):8.1e} {int(f(r,'n_lit')):4d}")
            break
print("barriers", len(env), " last: M* = %.2f, M_pt = %.2f, Mdot = %.3e, t_c = %.1f yr" % (f(env[-1], "Mstar"), f(env[-1], "Mpt"), f(env[-1], "Mdot_win"), f(env[-1], "tc_yr")))
print("max Q/Qcrit over the run: %.3e ; max Qesc_in: %.3e ; max a/g: %.3e ; max kLW: %.3e ; max Fsum/(Mdot v_ff): %.3e" % (max(f(r, "Q_over_Qcrit") for r in env), max(f(r, "Qesc_in") for r in env), max(f(r, "arad_over_g_max") for r in env), max(f(r, "kLW_max") for r in env), max(f(r, "Fsum_over_Mdotv") for r in env)))
