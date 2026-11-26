#!/usr/bin/env python3
"""
export_landmarks_data.py

Packages all Titan.jl landmark benchmark CSVs (lm01..lm08), exact ABN02 Fig. 2
digitized reference tables, Pop III protostar/pulsation/episodic/GR RHD runs,
and live 1D RHD H II/He III/LW breakout data into compact JSON files under
report/app/data/ for the interactive Titan.jl Laboratory & Quarto book widgets.
"""
import csv, glob, json, math, os, re, shutil
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
LM_DIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
RUNS_DIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs"
ABN_RED = os.path.join(ROOT, "lib", "popiii", "abn02red", "abn02_red_fig2.csv")
ABN_FIN = os.path.join(ROOT, "lib", "popiii", "abn02fin", "abn02_final_fig2.csv")
OUT_DATA = os.path.join(ROOT, "report", "app", "data")
OUT_FIGS = os.path.join(ROOT, "report", "figures")
BRAIN_DIR = "/home/tabel/.gemini/antigravity/brain/91f66bf7-1694-46da-a2b7-de91d0662618"

os.makedirs(OUT_DATA, exist_ok=True)
os.makedirs(OUT_FIGS, exist_ok=True)

def r4(x):
    if x is None or not math.isfinite(x):
        return None
    if x == 0.0:
        return 0.0
    return float(f"{x:.4g}")

def load_csv(path, stride=1, max_pts=None):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    hdr = [h.strip() for h in rows[0]]
    data = [r for r in rows[1:] if r and r[0].strip()]
    if max_pts and len(data) > max_pts:
        idx = np.unique(np.linspace(0, len(data) - 1, max_pts).round().astype(int))
        data = [data[i] for i in idx]
    elif stride > 1:
        data = data[::stride]
    out = {}
    for ci, h in enumerate(hdr):
        vals = []
        is_num = True
        for r in data:
            v = r[ci].strip()
            try:
                vals.append(r4(float(v)))
            except ValueError:
                is_num = False
                vals.append(v)
        out[h] = vals
    return out

def load_abn(path):
    tab = {}
    with open(path) as f:
        for l in f:
            if l.startswith("#") or l.startswith("curve") or not l.strip():
                continue
            parts = l.strip().split(",")
            k, x, y = parts[0], float(parts[1]), float(parts[2])
            tab.setdefault(k, {"x": [], "y": []})
            tab[k]["x"].append(r4(x))
            tab[k]["y"].append(r4(y))
    return tab

def dump(name, obj):
    p = os.path.join(OUT_DATA, name)
    with open(p, "w") as f:
        json.dump(obj, f, separators=(",", ":"))
    sz = os.path.getsize(p) / 1024.0
    print(f"  wrote {p} ({sz:.1f} KiB)")

def export_shocks_collapse():
    subc = load_csv(os.path.join(LM_DIR, "lm01_subc.csv"), max_pts=160)
    supc = load_csv(os.path.join(LM_DIR, "lm01_supc.csv"), max_pts=160)
    jruns = {}
    for key, fn in [
        ("J0", "lm02_collapse_J0.csv"),
        ("J10", "lm02_collapse_J10.csv"),
        ("J30", "lm02_collapse_J30.csv"),
        ("J100", "lm02_collapse_J100.csv"),
        ("J1e3", "lm02_collapse_J1e3_T1e4.csv"),
        ("J1e4", "lm02_collapse_J1e4_T1e5.csv"),
    ]:
        jruns[key] = load_csv(os.path.join(LM_DIR, fn), max_pts=180)
    h2lines = load_csv(os.path.join(LM_DIR, "lm02_h2_lines_J0.csv"), max_pts=160)
    dump("lm_shocks_collapse.json", {
        "ensman": {"subc": subc, "supc": supc},
        "omukai": {"runs": jruns, "h2lines": h2lines},
    })

def export_feedback_popiii():
    halos = {
        "3e5": load_csv(os.path.join(LM_DIR, "lm03_halo_3e5.csv")),
        "1e6": load_csv(os.path.join(LM_DIR, "lm03_halo_1e6.csv")),
        "5e6": load_csv(os.path.join(LM_DIR, "lm03_halo_5e6.csv")),
    }
    profs = {
        "0Myr":    load_csv(os.path.join(LM_DIR, "lm03_prof_1e6_0Myr.csv")),
        "0p25Myr": load_csv(os.path.join(LM_DIR, "lm03_prof_1e6_0p25Myr.csv")),
        "0p8Myr":  load_csv(os.path.join(LM_DIR, "lm03_prof_1e6_0p8Myr.csv")),
        "2p5Myr":  load_csv(os.path.join(LM_DIR, "lm03_prof_1e6_2p5Myr.csv")),
    }
    schaerer = load_csv(os.path.join(LM_DIR, "lm04_schaerer_zams.csv"))
    raiter   = load_csv(os.path.join(LM_DIR, "lm04_raiter_density.csv"))
    dump("lm_feedback_popiii.json", {
        "whalen": {"halos": halos, "profs_1e6": profs},
        "schaerer": schaerer,
        "raiter": raiter,
    })

def export_exotics_pisn():
    mcoll = load_csv(os.path.join(LM_DIR, "lm05_sms_mcoll_curve.csv"))
    gr_all = load_csv(os.path.join(LM_DIR, "lm05_sms_gr_collapse_tracks.csv"))
    gr_cases = {}
    for i, c in enumerate(gr_all["case"]):
        gr_cases.setdefault(c, {"t_tdyn": [], "rho_rat": [], "umin_c": [], "R_Rsun": []})
        for k in ("t_tdyn", "rho_rat", "umin_c", "R_Rsun"):
            gr_cases[c][k].append(gr_all[k][i])
    qs_all = load_csv(os.path.join(LM_DIR, "lm06_quasistar_seq.csv"))
    qs_05 = {k: [] for k in qs_all}
    for i, f in enumerate(qs_all["f_bh"]):
        if abs(f - 0.05) < 1e-3:
            for k in qs_all:
                qs_05[k].append(qs_all[k][i])
    qs_rhd = load_csv(os.path.join(LM_DIR, "lm06_quasistar_rhd_prof.csv"), max_pts=140)
    ds_seq = load_csv(os.path.join(LM_DIR, "lm06_darkstar_seq.csv"))
    seds   = load_csv(os.path.join(LM_DIR, "lm06_jwst_seds.csv"), max_pts=120)

    # Pair-instability 2D grid (reshape neatly into xs, ys, z)
    pm = load_csv(os.path.join(LM_DIR, "lm07_pair_instability_map.csv"))
    xs = sorted(set(pm["log10_rho"]))
    ys = sorted(set(pm["log10_T"]))
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    z_p = [None] * (len(xs) * len(ys))
    z_n = [None] * (len(xs) * len(ys))
    for i in range(len(pm["log10_rho"])):
        ix = xi[pm["log10_rho"][i]]
        iy = yi[pm["log10_T"][i]]
        z_p[iy * len(xs) + ix] = pm["Gamma1_pairs"][i]
        z_n[iy * len(xs) + ix] = pm["Gamma1_nopairs"][i]

    pisn_lc  = load_csv(os.path.join(LM_DIR, "lm07_pisn_breakout_lc.csv"), max_pts=150)
    ppisn_lc = load_csv(os.path.join(LM_DIR, "lm07_ppisn_collision_lc.csv"), max_pts=150)
    dump("lm_exotics_pisn.json", {
        "sms": {"mcoll": mcoll, "gr": gr_cases},
        "quasistar": {"seq": qs_05, "rhd": qs_rhd},
        "darkstar": ds_seq,
        "seds": seds,
        "pisn": {
            "pair_map": {"xs": xs, "ys": ys, "gamma1_pairs": z_p, "gamma1_nopairs": z_n},
            "breakout_lc": pisn_lc,
            "ppisn_lc": ppisn_lc,
        },
    })

def export_abn02():
    dump("lm_abn02.json", {
        "red": load_abn(ABN_RED),
        "fin": load_abn(ABN_FIN),
        "at_1e12": load_csv(os.path.join(LM_DIR, "lm08_abn02_at_1e12.csv")),
        "core_birth": load_csv(os.path.join(LM_DIR, "lm08_abn02_core_birth.csv")),
        "collapse_track": load_csv(os.path.join(LM_DIR, "lm08_abn02_collapse_track.csv"), max_pts=180),
        "protostar": load_csv(os.path.join(LM_DIR, "lm08_abn02_protostar_track.csv"), max_pts=160),
        "photoevap_polar": load_csv(os.path.join(LM_DIR, "lm08_abn02_photoevap_polar.csv")),
        "photoevap_1d": load_csv(os.path.join(LM_DIR, "lm08_abn02_photoevap_1d.csv")),
    })

def export_popiii_runs():
    catalog = [
        ("hoy12_md1e1_full", "hoy12_md1e1_puls", "HOY12 0.1 M☉/yr (2→48 M☉) + κ-mechanism pulsations (P=0.82 yr, +50% inversion)"),
        ("hoy12_md3e1_full", "hoy12_md3e1_puls", "HOY12 0.3 M☉/yr (2→109 M☉) + pulsations (P=1.44 yr, +26% inversion)"),
        ("episodic_cycle",   "episodic_cycle",   "Episodic accretion burst cycle (6 yr quiet @ 1e-3 / 4 yr burst @ 0.25 M☉/yr)"),
        ("episodic_short",   "episodic_short",   "Episodic quiescent contraction (12 yr @ 1e-4 M☉/yr; R* 1313→645 R☉)"),
        ("hoy12_md3e2",      "hoy12_md3e2",      "HOY12 0.03 M☉/yr (1→11.3 M☉) full 1D Eulerian RHD + Stage-2 TDC"),
    ]
    pdir = os.path.join(RUNS_DIR, "profiles")
    runs_out = {}
    for key, csv_stem, desc in catalog:
        cpath = os.path.join(RUNS_DIR, f"{csv_stem}.csv")
        if not os.path.exists(cpath):
            continue
        ts = load_csv(cpath, max_pts=160)
        # Collect up to 14 evenly spaced radial profiles for interactive click-to-inspect
        pfiles = sorted(glob.glob(os.path.join(pdir, f"{csv_stem}_s*.csv"))) if os.path.isdir(pdir) else []
        snaps = []
        if pfiles:
            idx = np.unique(np.linspace(0, len(pfiles) - 1, min(14, len(pfiles))).round().astype(int))
            for i in idx:
                pf = pfiles[i]
                m = re.search(r"_s(\d+)_", os.path.basename(pf))
                st = int(m.group(1)) if m else i
                pr = load_csv(pf, max_pts=110)
                keep = {k: pr[k] for k in ("r_Rsun", "m_Msun", "rho", "T", "u_kms", "Gamma_Edd", "L_tot_Lsun") if k in pr}
                snaps.append({"step": st, "prof": keep})
        keep_ts = {k: ts[k] for k in ("step", "t_yr", "M_star", "R_star", "R_ph", "Teff", "L_ph", "L_Edd", "Tc", "max_dlnrho_dlnr", "max_Gamma_Edd", "Mdot_cur") if k in ts}
        runs_out[key] = {"desc": desc, "ts": keep_ts, "snaps": snaps}

    # Add Live 1D RHD H II / He III / LW Breakout run
    hii_path = os.path.join(RUNS_DIR, "hii_breakout_rhd.csv")
    hii = load_csv(hii_path, max_pts=160) if os.path.exists(hii_path) else None
    dump("popiii_runs.json", {"runs": runs_out, "hii_breakout": hii})

def copy_figures():
    figs = [
        "figL1_ensman_and_omukai_collapse.png",
        "figL2_whalen04_and_schaerer02.png",
        "figL3_sms_quasistar_darkstar.png",
        "figL4_pisn_shock_breakout_and_ppisn.png",
        "figL5_abn02_end_to_end.png",
        "fig7_dynamical_hii_breakout_and_lines.png",
    ]
    for fn in figs:
        src = os.path.join(BRAIN_DIR, fn)
        dst = os.path.join(OUT_FIGS, fn)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            print(f"  copied figure {fn} ({os.path.getsize(dst)/1024:.0f} KiB)")

def main():
    print("Exporting Titan.jl Landmark & Pop III datasets for GitHub Pages...")
    export_shocks_collapse()
    export_feedback_popiii()
    export_exotics_pisn()
    export_abn02()
    export_popiii_runs()
    copy_figures()
    print("Done.")

if __name__ == "__main__":
    main()
