# S7b: breakout accretion rate of the H II region in spherical free-fall accretion (OI02 eq. 5-7) along the stored/offline star track:
# Mdot_break(M) = the rate below which R_II (Q = alpha Mdot'^2/(8 pi G m_p^2 M) ln(R_II/R*) + Mdot'/m_p) reaches r_s, for the ionising rate Q*(M) of the stellar surface.
# Usage: python3 -I s7b_break.py <s7b_offline csv> [alpha]
import sys, csv, math
G = 6.67430e-8; mp = 1.67262192e-24; Msun = 1.98847e33; YR = 3.15576e7; AU = 1.495978707e13; Rsun = 6.957e10
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0e-13       # case B at ~1.5e4 K
rows = list(csv.DictReader(open(sys.argv[1])))
print(f"{'M':>7s} {'R*[Rsun]':>9s} {'r_s[AU]':>8s} {'Q*[1/s]':>10s} {'Mdot_env':>10s} {'Q/Qcrit':>9s} {'Mdot_break':>11s} {'env/break':>9s}")
for r in rows:
    M = float(r["M"]); Md = float(r["Mdot"]); Q = float(r["Q_star"]); R = float(r["R_Rsun"]) * Rsun; rs = float(r["rs_AU"]) * AU
    L = math.log(rs / R)
    if L <= 0: continue
    # solve Q = alpha x^2/(8 pi G mp^2 M) L + x/mp for x = Mdot [g/s] (positive root)
    A = alpha * L / (8 * math.pi * G * mp**2 * M * Msun); B = 1 / mp
    x = (-B + math.sqrt(B * B + 4 * A * Q)) / (2 * A)
    xb = x * YR / Msun
    print(f"{M:7.2f} {float(r['R_Rsun']):9.1f} {float(r['rs_AU']):8.1f} {Q:10.2e} {Md:10.2e} {float(r['Q_over_Qcrit']):9.2e} {xb:11.2e} {Md / xb:9.1f}")
