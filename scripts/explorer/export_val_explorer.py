#!/usr/bin/env python3
"""
export_val_explorer.py

Export Titan.jl validation / accretion / episodic / supermassive RHD runs
(from val_const_mdot.jl CSVs + radial profile CSVs) into the exact Pop III Run Explorer
binary format (runs.json, <run>/frames.json, <run>/frames.b64.txt) compatible with
scripts/explorer/panel_index_v82.html.

Usage:
  python3 scripts/explorer/export_val_explorer.py \
      --runs-dir /zpool/nvme/data/tabel_scratch/titan_popiii/val/runs \
      --out-dir  /zpool/nvme/data/tabel_scratch/titan_popiii/val/explorer
"""
import argparse, base64, datetime, glob, json, math, os, re, shutil, time
import numpy as np

MSUN = 1.98847e33
RSUN = 6.957e10
YR   = 3.15576e7
MP   = 1.6726219e-24

RUN_CATALOG = [
    ("hoy12_md1e1_full", ["hoy12_md1e1_puls", "hoy12_md1e1_zoom"],
     "HOY12 MD1e1 (0.1 M☉/yr, 2→48 M☉) + high-cadence pulsation zoom (P=0.82 yr, Q=0.062 d, +50% density inversion)"),
    ("hoy12_md3e1_full", ["hoy12_md3e1_puls", "hoy12_md3e1_zoom"],
     "HOY12 MD3e1 (0.3 M☉/yr, 2→108.5 M☉) + high-cadence pulsation zoom (P=1.44 yr, Q=0.073 d, +26% density inversion)"),
    ("episodic_cycle", ["episodic_cycle"],
     "Episodic full burst cycle (6 yr quiet @ 1e-3 / 4 yr burst @ 0.25 M☉/yr from 47.4 M☉ supergiant)"),
    ("episodic_short", ["episodic_short"],
     "Episodic quiescent contraction (12 yr @ 1e-4 M☉/yr; R* 1313→645 R☉, outward shock detachment)"),
    ("sms_md1e0_gr", ["sms_md1e0_gr"],
     "Supermassive star formation with 1D spherical TOV/1PN GR corrections (1.0 M☉/yr, 10→10⁴ M☉)"),
    ("hoy12_md3e2", ["hoy12_md3e2"],
     "HOY12 MD3e2 (0.03 M☉/yr, 1→11.3 M☉) full 1D Eulerian RHD + Stage-2 TDC"),
]

FIELDS = [
    ["logr", "log10 r [Rsun]"], ["logm", "log10 m [Msun]"], ["logrho", "log10 rho [g/cm^3]"],
    ["logT", "log10 T [K]"], ["logTrad", "log10 Trad [K]"], ["u", "u [km/s]"],
    ["L", "L_tot [Lsun]"], ["logtau", "log10 tau_R"], ["Gamma", "Gamma_Edd"],
    ["Mdot", "|4pi r^2 rho u| [Msun/yr]"], ["logr_o", "log10 r_outer [Rsun]"], ["logm_o", "log10 m_outer [Msun]"],
    ["logXH", "log10 X_H"], ["logXD", "log10 (D/D0 * 3.75e-5)"],
    ["f_conv", "L_conv / L_tot"], ["f_pturb", "P_turb / P_tot"], ["dlnrho", "dln rho / dln r"],
    ["eps_D", "eps_D [erg/g/s]"], ["eps_H", "eps_H [erg/g/s]"], ["xi_GR", "2Gm / (r c^2)"],
]

SCALARS = [
    "step", "t_yr", "dt_s", "logn_c", "T_c", "M_star", "R_star", "r_ph",
    "T_eff_lum", "L_ph", "u_min", "L_nuc", "nz",
    "max_dlnrho", "max_Gam", "max_PtP", "max_xi_gr", "S_EUV", "Mdot_cur"
]

def load_csv(path):
    with open(path) as f:
        hdr = f.readline().strip().split(",")
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    if arr.size == 0:
        return None
    if arr.ndim == 1:
        arr = arr[None, :]
    return {h: arr[:, i] for i, h in enumerate(hdr)}

def safe_log10(x, floor=-30.0):
    arr = np.asarray(x, dtype=np.float64)
    out = np.full_like(arr, floor)
    m = np.isfinite(arr) & (arr > 0)
    out[m] = np.log10(arr[m])
    return out.astype(np.float32)

def collect_profiles(runs_dir, tags):
    pdir = os.path.join(runs_dir, "profiles")
    if not os.path.isdir(pdir):
        return {}
    by_step = {}
    for tag in tags:
        for fn in sorted(glob.glob(os.path.join(pdir, f"{tag}_s*.csv"))):
            m = re.search(r"_s(\d+)_", os.path.basename(fn))
            if not m:
                continue
            st = int(m.group(1))
            by_step[st] = fn
    return by_step

def build_profile_arrays(path):
    p = load_csv(path)
    if p is None:
        return None
    n = len(p["r_Rsun"])
    rc = np.maximum(p["r_Rsun"], 1e-6)
    dr = np.maximum(p["dr_Rsun"], 1e-8)
    ro = rc + 0.5 * dr
    m  = np.maximum(p["m_Msun"], 1e-12)
    rho = np.maximum(p["rho"], 1e-30)
    T   = np.maximum(p["T"], 10.0)
    Pr  = np.maximum(p["P_rad"], 1e-30)
    Trad = np.power(3.0 * Pr / 7.5657e-15, 0.25)
    u   = p["u_kms"]
    Lrad = p["L_rad_Lsun"]
    Lconv = p["L_conv_Lsun"]
    Ltot = Lrad + Lconv
    Ptot = np.maximum(p["P_gas"] + p["P_rad"] + p["P_turb"], 1e-30)
    f_conv = np.clip(Lconv / np.maximum(np.abs(Ltot), 1.0), -1.0, 1.5)
    f_pturb = np.clip(p["P_turb"] / Ptot, 0.0, 1.0)
    mdot = 4.0 * np.pi * (rc * RSUN) ** 2 * rho * np.abs(u * 1e5) / (MSUN / YR)
    eps_D = p["eps_D_cgs"] if "eps_D_cgs" in p else np.zeros(n)
    eps_H = p["eps_H_cgs"] if "eps_H_cgs" in p else np.zeros(n)
    xi_GR = p["xi_GR"] if "xi_GR" in p else (2.0 * 6.6743e-8 * (m * MSUN) / ((ro * RSUN) * 2.99792458e10 ** 2))
    out = np.zeros((len(FIELDS), n), dtype=np.float32)
    vals = [
        safe_log10(rc), safe_log10(m), safe_log10(rho), safe_log10(T), safe_log10(Trad),
        u.astype(np.float32), Ltot.astype(np.float32), safe_log10(np.maximum(p["tau_R"], 1e-12)),
        p["Gamma_Edd"].astype(np.float32), mdot.astype(np.float32),
        safe_log10(ro), safe_log10(m),
        safe_log10(np.maximum(p["X_H"], 1e-20)),
        safe_log10(np.maximum(p["D_over_D0"] * 3.75e-5, 1e-20)),
        f_conv.astype(np.float32), f_pturb.astype(np.float32), p["dlnrho_dlnr"].astype(np.float32),
        eps_D.astype(np.float32), eps_H.astype(np.float32), xi_GR.astype(np.float32),
    ]
    for j, v in enumerate(vals):
        out[j, :] = np.nan_to_num(v, nan=0.0, posinf=0.0, neginf=0.0)
    return out, n, float(np.min(u))

def synth_profile_from_row(row, n=32):
    # Lightweight 1D fallback for time steps between full profile dumps so high-cadence pulsation
    # time series & FFT have EVERY single step while profile scrubber snaps smoothly to nearest profile!
    r_s = max(row["R_star"], 1.0)
    r_p = max(row["R_ph"], r_s)
    m_s = max(row["M_star"], 0.1)
    rc = np.geomspace(max(0.02 * r_s, 0.5), r_p * 1.5, n)
    ro = rc * 1.08
    m  = m_s * np.minimum(1.0, (rc / r_s) ** 1.5)
    rho = np.full(n, max(row["rho_c"], 1e-12)) * np.power(rc / rc[0], -2.0)
    T = np.full(n, max(row["Tc"], 1e4)) * np.power(rc / rc[0], -0.6)
    out = np.zeros((len(FIELDS), n), dtype=np.float32)
    out[0] = safe_log10(rc); out[1] = safe_log10(m); out[2] = safe_log10(rho)
    out[3] = safe_log10(T);  out[4] = safe_log10(T)
    out[6] = row["L_ph"];    out[8] = row.get("max_Gamma_Edd", 0.0)
    out[10] = safe_log10(ro); out[11] = safe_log10(m)
    return out, n, float(row.get("u_ph_kms", 0.0))

def export_run(key, tags, desc, runs_dir, out_dir, max_frames=500):
    dfs = []
    for t in tags:
        fn = os.path.join(runs_dir, f"{t}.csv")
        if os.path.exists(fn) and os.path.getsize(fn) > 200:
            d = load_csv(fn)
            if d is not None:
                dfs.append(d)
    if not dfs:
        return None
    # Stitch cleanly by strictly increasing step
    all_keys = sorted(set().union(*[d.keys() for d in dfs]))
    rows = []
    last_step = -1
    for d in dfs:
        n = len(d["step"])
        for i in range(n):
            st = int(round(d["step"][i]))
            if st <= last_step:
                continue
            last_step = st
            rows.append({k: float(d[k][i]) if k in d else float("nan") for k in all_keys})
    if len(rows) < 2:
        return None

    profs = collect_profiles(runs_dir, tags)
    prof_steps = sorted(profs.keys())
    # Keep ALL profile steps + uniformly subsampled time-series rows up to max_frames
    keep_idx = set()
    step_to_idx = {int(round(r["step"])): idx for idx, r in enumerate(rows)}
    for ps in prof_steps:
        if ps in step_to_idx:
            keep_idx.add(step_to_idx[ps])
    # Add high-cadence points (especially last 250 steps for pulsation FFT!)
    if len(rows) <= max_frames:
        keep_idx.update(range(len(rows)))
    else:
        tail_n = min(220, max_frames // 2)
        keep_idx.update(range(len(rows) - tail_n, len(rows)))
        head_idx = np.linspace(0, len(rows) - tail_n - 1, max_frames - tail_n, dtype=int)
        keep_idx.update(int(i) for i in head_idx)
    indices = sorted(keep_idx)

    # Cache loaded profiles so nearby time-series steps share exact nearest full profile
    prof_cache = {}
    def nearest_prof(st):
        if not prof_steps:
            return None
        pos = np.searchsorted(prof_steps, st)
        cand = []
        if pos < len(prof_steps): cand.append(prof_steps[pos])
        if pos > 0: cand.append(prof_steps[pos - 1])
        best = min(cand, key=lambda s: abs(s - st))
        if best not in prof_cache:
            prof_cache[best] = build_profile_arrays(profs[best])
        return prof_cache[best]

    rdir = os.path.join(out_dir, key)
    os.makedirs(rdir, exist_ok=True)
    frames_meta = []
    chunks = []
    off = 0
    for idx in indices:
        r = rows[idx]
        st = int(round(r["step"]))
        pr = nearest_prof(st)
        if pr is None:
            pr = synth_profile_from_row(r)
        arr, nz, umin = pr
        raw = arr.tobytes()
        chunks.append(raw)
        lognc = math.log10(max(r.get("rho_c", 1e-10), 1e-30) / MP)
        sv = [
            st, r.get("t_yr", 0.0), r.get("dt_s", 0.0), lognc, r.get("Tc", float("nan")),
            r.get("M_star", float("nan")), r.get("R_star", float("nan")), r.get("R_ph", float("nan")),
            r.get("Teff", float("nan")), r.get("L_ph", float("nan")), umin, r.get("L_nuc", 0.0), nz,
            r.get("max_dlnrho_dlnr", float("nan")), r.get("max_Gamma_Edd", float("nan")),
            r.get("max_Pt_over_P", float("nan")), r.get("max_2Gm_rc2", float("nan")),
            r.get("S_EUV_s", float("nan")), r.get("Mdot_cur", float("nan")),
        ]
        sv = [None if (isinstance(v, float) and not math.isfinite(v)) else v for v in sv]
        frames_meta.append({"off": off, "n": nz, "s": sv})
        off += len(raw)

    bin_all = b"".join(chunks)
    with open(os.path.join(rdir, "frames.json"), "w") as f:
        json.dump({"fields": FIELDS, "scalars": SCALARS, "frames": frames_meta}, f)
    with open(os.path.join(rdir, "frames.b64.txt"), "w") as f:
        f.write(base64.b64encode(bin_all).decode("ascii"))

    last = rows[-1]
    return {
        "key": key, "desc": desc, "process": False, "frames": True, "status": "stopped",
        "step": int(round(last["step"])), "t_yr": last.get("t_yr"), "dt_s": last.get("dt_s"),
        "M_star": last.get("M_star"), "R_star": last.get("R_star"),
        "zones": frames_meta[-1]["n"], "threads": 2,
        "rate_Msun_yr": last.get("Mdot_cur") if math.isfinite(last.get("Mdot_cur", float("nan"))) else None,
        "perf_mtime": int(time.time()),
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-dir", default="/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs")
    ap.add_argument("--out-dir",  default="/zpool/nvme/data/tabel_scratch/titan_popiii/val/explorer")
    ap.add_argument("--html",     default=os.path.join(os.path.dirname(__file__), "panel_index_v82.html"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    shutil.copy2(args.html, os.path.join(args.out_dir, "index.html"))
    manifest_runs = []
    for key, tags, desc in RUN_CATALOG:
        entry = export_run(key, tags, desc, args.runs_dir, args.out_dir)
        if entry:
            manifest_runs.append(entry)
            print(f"Exported {key:20s}: M*={entry['M_star']:.2f} Msun, R*={entry['R_star']:.1f} Rsun, step={entry['step']}")
    now = int(time.time())
    doc = {"updated": now, "updated_iso": datetime.datetime.fromtimestamp(now).astimezone().isoformat(timespec="minutes"), "runs": manifest_runs}
    with open(os.path.join(args.out_dir, "runs.json"), "w") as f:
        json.dump(doc, f, indent=1, allow_nan=False)
    print(f"Wrote {os.path.join(args.out_dir, 'runs.json')} ({len(manifest_runs)} runs) + index.html")

if __name__ == "__main__":
    main()
