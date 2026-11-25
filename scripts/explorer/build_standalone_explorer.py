#!/usr/bin/env python3
"""
build_standalone_explorer.py

Refresh the embedded self-contained JSON payload inside popiii_explorer_v82.html
from the latest Titan.jl validation/accretion/episodic/supermassive RHD CSVs & profiles.
Computes exact accretion shock radius R_shock (peak deceleration -du/cs inside r < 0.6 r_max)
for display on the HOY12 Fig. 1 Kippenhahn diagram and radial profile panels.
"""
import argparse, glob, json, math, os, re
import numpy as np

DEFAULT_RUNS = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs"
DEFAULT_HTML = "/home/tabel/.gemini/antigravity/brain/91f66bf7-1694-46da-a2b7-de91d0662618/popiii_explorer_v82.html"

def load_csv(p):
    with open(p) as f:
        hdr = f.readline().strip().split(",")
    a = np.loadtxt(p, delimiter=",", skiprows=1)
    if a.size == 0:
        return None
    if a.ndim == 1:
        a = a[None, :]
    return {h: a[:, i] for i, h in enumerate(hdr)}

def r3(x):
    return None if not np.isfinite(x) else float(f"{x:.3g}")

def rt(x):
    return None if not np.isfinite(x) else round(float(x), 5)

def euv_rate(Rph_rsun, Teff):
    if Teff <= 500 or Rph_rsun <= 0:
        return 0.0
    x0 = 157807.0 / Teff
    if x0 > 600:
        return 0.0
    s = sum(math.exp(-n * x0) * (x0 * x0 / n + 2 * x0 / (n * n) + 2 / (n * n * n)) for n in range(1, 5))
    return (4 * math.pi * (Rph_rsun * 6.957e10) ** 2 * 2 * math.pi * (1.380649e-16 * Teff) ** 3
            / (6.62607015e-27 ** 3 * 2.99792458e10 ** 2) * s)

def find_shock_radius(rc, u, cs):
    """Locate primary accretion shock front: strongest inward deceleration -du/cs inside r < 0.6 r_max."""
    n = len(rc)
    if n < 5:
        return None
    rmax = 0.60 * rc[-1]
    du = np.diff(u)
    csm = np.maximum(0.5 * (cs[:-1] + cs[1:]), 1e-3)
    comp = -du / csm  # positive where infalling gas (u < 0) decelerates sharply
    mask = (rc[:-1] < rmax) & (u[1:] < -2.0) & (comp > 0.4)
    if not np.any(mask):
        return None
    cand = np.where(mask)[0]
    kbest = cand[int(np.argmax(comp[cand]))]
    return float(0.5 * (rc[kbest] + rc[kbest + 1]))

SPECS = [
    ("hoy12_md1e1", ["hoy12_md1e1_puls", "hoy12_md1e1_zoom"], "HOY12 0.1 M☉/yr (2→48 M☉ + pulsation zoom)", "#2a78d6", 0.1),
    ("hoy12_md3e1", ["hoy12_md3e1_puls", "hoy12_md3e1_zoom"], "HOY12 0.3 M☉/yr (2→108 M☉ + pulsation zoom)", "#eb6834", 0.3),
    ("sms_md1e0_gr", ["sms_md1e0_gr"], "Supermassive 1.0 M☉/yr + 1D spherical GR", "#7d5fd3", 1.0),
    ("episodic_short", ["episodic_short"], "Episodic 40yr quiet / 8yr burst @ 0.2 M☉/yr", "#1baf7a", 1e-4),
    ("episodic_10yr", ["episodic_10yr"], "Episodic 10yr quiet / 5yr burst @ 0.2 M☉/yr", "#eda100", 1e-4),
]

def build(runs_dir, html_path):
    pdir = os.path.join(runs_dir, "profiles")
    out = {}
    shells = [1, 3, 10, 30, 100, 300, 1000]
    for key, tags, desc, col, md_def in SPECS:
        rows = {}
        for t in tags:
            p = os.path.join(runs_dir, f"{t}.csv")
            if not os.path.exists(p) or os.path.getsize(p) < 200:
                continue
            c = load_csv(p)
            if c is None:
                continue
            for i in range(len(c["step"])):
                st = int(c["step"][i])
                rows[st] = {k: float(v[i]) for k, v in c.items()}
        if not rows:
            continue
        sts = sorted(rows.keys())
        idx_g = np.unique(np.linspace(0, len(sts) - 1, min(80, len(sts)), dtype=int))
        idx_z = np.arange(max(0, len(sts) - 260), len(sts))
        idx = np.unique(np.concatenate([idx_g, idx_z]))
        trk = []
        for ii in idx:
            r = rows[sts[ii]]
            seuv = r["S_EUV_s"] if "S_EUV_s" in r else euv_rate(r["R_ph"], r["Teff"])
            xi = (r["max_2Gm_rc2"] if "max_2Gm_rc2" in r
                  else 2 * 6.6743e-8 * (r["M_star"] * 1.98847e33) / (max(r["R_star"], 0.1) * 6.957e10 * 2.99792458e10 ** 2))
            trk.append([
                sts[ii], rt(r["t_yr"]), r3(r["M_star"]), r3(r["R_star"]), r3(r["R_ph"]),
                r3(r["Teff"]), r3(r["L_ph"]), r3(r.get("max_dlnrho_dlnr", 0.0)), r3(r.get("max_Gamma_Edd", 0.0)),
                r3(xi), r3(seuv), r3(r.get("Mdot_cur", md_def))
            ])
        pfiles = {}
        for t in tags:
            for fn in sorted(glob.glob(os.path.join(pdir, f"{t}_s*.csv"))):
                m = re.search(r"_s(\d+)_", os.path.basename(fn))
                if m:
                    pfiles[int(m.group(1))] = fn
        psts = sorted(pfiles.keys())
        pidx = np.unique(np.linspace(0, len(psts) - 1, min(22, len(psts)), dtype=int)) if psts else []
        profs = []
        for pi in pidx:
            st = psts[pi]
            pf = load_csv(pfiles[st])
            if pf is None:
                continue
            rc = np.maximum(pf["r_Rsun"], 1e-6)
            m = np.maximum(pf["m_Msun"], 1e-12)
            lr = np.log10(rc)
            n = len(rc)
            Mst = rows[st]["M_star"] if st in rows else float(m[-1])
            Rst = rows[st]["R_star"] if st in rows else float(rc[-1])
            Rph = rows[st]["R_ph"] if st in rows else Rst
            Rsh = find_shock_radius(rc, pf["u_kms"], pf["cs_kms"])
            sh = []
            for ms in shells:
                if Mst >= ms * 1.02:
                    jj = np.searchsorted(m, ms)
                    if 0 < jj < n:
                        f = (ms - m[jj - 1]) / max(m[jj] - m[jj - 1], 1e-30)
                        sh.append(r3(lr[jj - 1] + f * (lr[jj] - lr[jj - 1])))
                    else:
                        sh.append(None)
                else:
                    sh.append(None)
            jj70 = np.searchsorted(m, 0.7 * Mst)
            r70 = r3(lr[min(jj70, n - 1)]) if Mst > 0.5 else None
            Lrad, Lconv = pf["L_rad_Lsun"], pf["L_conv_Lsun"]
            fc = Lconv / np.maximum(np.abs(Lrad) + np.abs(Lconv), 1e-30)
            cb = []
            in_c = False
            r0c = 0.0
            for k in range(n):
                if rc[k] > Rst * 1.05:
                    break
                if fc[k] > 0.02 and not in_c:
                    in_c = True
                    r0c = lr[k]
                elif fc[k] <= 0.02 and in_c:
                    in_c = False
                    cb.append([r3(r0c), r3(lr[k - 1])])
            if in_c:
                cb.append([r3(r0c), r3(min(lr[k], np.log10(max(Rst, 0.1))))])
            ki = np.unique(np.linspace(0, n - 1, min(40, n), dtype=int))
            rho = np.maximum(pf["rho"], 1e-30)
            dlnr = pf["dlnrho_dlnr"] if "dlnrho_dlnr" in pf else np.zeros(n)
            xi = (pf["xi_GR"] if "xi_GR" in pf
                  else 2 * 6.6743e-8 * (m * 1.98847e33) / (rc * 6.957e10 * 2.99792458e10 ** 2))
            profs.append({
                "st": st, "t": rt(rows[st]["t_yr"]) if st in rows else 0,
                "M": r3(Mst), "R": r3(Rst), "Rph": r3(Rph), "Rsh": r3(Rsh) if Rsh else None,
                "Te": r3(rows[st]["Teff"]) if st in rows else 0,
                "L": r3(rows[st]["L_ph"]) if st in rows else 0,
                "sh": sh, "r70": r70, "cb": cb,
                "lr": [r3(lr[k]) for k in ki],
                "lT": [r3(np.log10(max(pf["T"][k], 10))) for k in ki],
                "lD": [r3(np.log10(rho[k])) for k in ki],
                "u":  [r3(pf["u_kms"][k]) for k in ki],
                "G":  [r3(np.log10(max(pf["Gamma_Edd"][k], 1e-12))) for k in ki],
                "fc": [r3(fc[k]) for k in ki],
                "dr": [r3(np.clip(dlnr[k], -8, 12)) for k in ki],
                "xi": [r3(np.log10(max(xi[k], 1e-12))) for k in ki],
            })
        out[key] = {"desc": desc, "col": col, "trk": trk, "profs": profs}
        print(f"Embedded {key:15s}: M*={rows[sts[-1]]['M_star']:7.2f} Msun, R*={rows[sts[-1]]['R_star']:7.1f} Rsun, t={rows[sts[-1]]['t_yr']:6.2f} yr, profs={len(profs)}")

    payload = json.dumps(out, separators=(",", ":"), ensure_ascii=False)
    with open(html_path) as f:
        html = f.read()
    s0 = html.index("const DATA = ")
    s1 = html.index(";\n", s0)
    html = html[:s0] + "const DATA = " + payload + html[s1:]
    with open(html_path, "w") as f:
        f.write(html)
    print(f"Wrote standalone artifact {html_path} ({len(html)/1024:.1f} KB)")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-dir", default=DEFAULT_RUNS)
    ap.add_argument("--html", default=DEFAULT_HTML)
    args = ap.parse_args()
    build(args.runs_dir, args.html)
