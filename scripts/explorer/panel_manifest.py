#!/usr/bin/env python3
# Writes panel/runs.json for the PopIII explorer: one entry per run that has a panel subdirectory or a live julia process
# (exact environ line TAG=_<run>), with the latest state from perf.csv and rates over the last ~hour of wall time.
import os, json, time, glob, datetime
S = "/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad"
D = "/zpool/nvme/data/tabel_scratch/titan_popiii/runs"
P = os.path.join(S, "panel")

DESC = {
    "zeul25": "TEST (20 Nov): as zeul1 (fixed log Eulerian mesh, code 06f2bd6) but 2500 zones: dln r 0.0102 from 0.01 Rsun to the IC edge; 16 threads cores 0-15",
    "zie3": "RESOLUTION TEST (21 Nov, Tom): zie1 at 3x resolution, 298 zones on the fixed Eulerian mesh (dln r 0.0853 vs 0.256), otherwise zie1's settings and code (iedef aace2a0), from the IC, cores 0-7: does the protostar's mass accretion history change? (the 100-zone fixed mesh may pin the I-front)",
    "zie1": "PRODUCTION (21 Nov, Tom: internal energy only): zeulnc with the dual-energy path deleted (branch iedef aace2a0): compatible internal-energy rows in every zone, no ETOTAL/ETASW* switches, total energy only monitored (ETOTMON=1); no CMB, from the IC, 100 zones, cores 32-39. zeulnc stalled at the photosphere front from an ETASWTIE row flip (13/93 Newton failures in a 40-step replay; 0 with internal energy only)",
    "zeulnc": "STOPPED (21 Nov, Tom: CMB off, it affects nothing here): zeul1 (fixed log Eulerian mesh, 100 zones, code 06f2bd6) with TCMB unset, from the IC; zeul1 lost the CMB inside the envelope (E_r < a T_CMB^4, spurious +2e4 L_sun) and its outer-boundary F_r row drove most Newton failures; cores 32-39",
    "zeul5": "TEST (20 Nov): as zeul1 (fixed log Eulerian mesh, code 06f2bd6) but 500 zones (496): dln r 0.0512 from 0.01 Rsun to the IC edge; to resolve the shock and photosphere; cores 40-47",
    "zeul1": "TEST (20 Nov): fixed logarithmic Eulerian mesh, no moving mesh: 100 zones, regular centre, first interface 0.01 Rsun, dln r 0.256 to the IC edge; NLTE=ali; code 06f2bd6; 8 threads",
    "zbase3": "BASELINE (20 Nov): 100 zones from the IC, NLTE=ali, mesh with iterate relaxation time smin(kappa t_cross, 0.2 t_dyn), iterate tau term, normalised monitor mixture a=0.4 (floor 6/decade), centre 1e-4 Msun; code dynrelax 06f2bd6; 8 threads",
    "z200a": "TEST (20 Nov): 200 zones from the IC, NLTE=ali, iterate-based mesh relaxation time (smin(kappa t_cross, 0.2 t_dyn)), tau monitor term wtau=3 window to 1e8, gmax 1000, resolution floor 6/decade, centre 1e-4 Msun; code dynrelax 27e6efd; 8 threads",
    "zbase2": "BASELINE attempt (20 Nov): zbase1 settings + NLTE=ali + resolution floor + kappa_f=0; failed at step 208 (dt collapse, mesh lagged the collapse); code 02549c7",
    "zbase1": "BASELINE (19 Nov): from the IC (IC=abn02red), complete z100nk settings with CONV=0.01 (Newton tolerance; DTERR 0.03) and ICMLOG=1e-4 (log-in-mass IC mesh, centre zone 1e-4 M☉), 100 zones GRIDEQ moving mesh, NKRYLOV + CLOSJ; code popiii-combined 5f833e8; 8 threads",
    "c128": "Split-chemistry production line, 141.5 to 272.6 M☉ (stopped 14 Nov; continued by t_c128c)",
    "t_c128c": "Former production (stopped 15 Nov at 555.44 M☉): τ≥3 table layout; superseded by t_nolte2",
    "t_es1": "Side test (option c): production layout with warm-front zones 97–110 (1.8e5–8.9e6 R☉) split ×3, from t_nolte2 at 555.5 M☉",
    "z512ic": "512 zones, purely logarithmic (Δlog₁₀ r 0.0213) from 0.1 R☉ to the cloud edge, fixed Eulerian, from the IC; performance study",
    "z512p": "Same 512-zone log mesh (508 after rescale to 143 pc) from the 0.05 M☉ protostar (f128c stage A3); performance study",
    "abnfin4": "IC=abn02fin smoke test with z100ic3 settings unchanged (core on the network, PSWTMIN=4500): failed at step 95 (only the step-50 snapshot)",
    "abnfin6": "IC=abn02fin smoke test, 96 uniform zones, core on the table (PSWTMIN=0 PSWITCH=1): 1500 steps, second core at 6.7 yr, M* 0.43 M☉",
    "z100mm": "From the IC with 100 zones on the GRIDEQ moving mesh (IC=abn02red, 7137 M☉ to 1e9 R☉) + G2 network; interface opacity and flux tie from the slab solution (SLABFLUX/SLABTIE, SUBZONE off); star light at the Newton iterate (SLIMP=2); full accurate stack, 8 threads",
    "z300mm": "MAIN from-IC run (16 Nov): GRIDEQ moving mesh (iterate-only monitor, 290 zones) + G2 network through the photosphere + SUBZONE + SLIMP, full accurate stack; continues t_mm24; IC=abn02red, 8 threads",
    "z300g2": "Fixed-mesh reference (16 Nov): zdiagB → E5 config: G2 network through the photosphere (SCVHXE, HM88, G2OPAC=2, table only n_H ≥ 1e19) + SUBZONE=1,2 + SLIMP=2; IC=abn02red, 291 zones, full accurate stack, 8 threads",
    "t_sz31": "zdiagB continued from S5 with SUBZONE=1,2 (front-interface optical depth) + SLIMP=2 (star light W and T_c at the Newton iterate); dt ≥ 1e4 s through the photosphere front",
    "zdiagB": "MAIN from-IC run (15 Nov): IC=abn02red (ABN02 red 9-Myr output, z=19 DM, recombined x_e), 7137 M☉ to 1e9 R☉; one central zone to 0.1 R☉ + 290 log zones (Δlog₁₀ r 0.0345); full accurate stack incl. VEFIMP, 8 threads",
    "zdiagA": "Diagnostic: z300red with VEFIMP=0 (10 linear centre zones kept); stopped after the collapse (log n_c 17)",
    "z300red": "MAIN from-IC run (15 Nov): IC=abn02red (ABN02 red 9-Myr output from the vector data, z=19 DM, recombined x_e), 7137 M☉ to 1e9 R☉; 300 fixed zones (10 linear to 0.1 R☉ + 290 log); most accurate stack, 8 threads",
    "z300acc": "From IC=abn02fin (recombined x_e), cloud cut at 1e9 R☉ (7372 M☉); 300 fixed zones: 10 linear 0–0.1 R☉ + 290 log (Δlog₁₀ r 0.0345) to 1e9 R☉; most accurate stack (full COUPLED, DTERR, NLOPAC, VEFIMP, LINEIMP, NLOIT, LOCALJ), 8 threads",
    "z100fin": "From the ABN02 FINAL output (IC=abn02fin, Fig. 2 vector data, n_c 3e13), 2e4 M☉ to 63.9 pc; z100s50 mesh (50 zones 0.1–10 R☉ + 49 out); core on the table (PSWTMIN=0 PSWITCH=1)",
    "z100s50b": "z100s50 mesh (50 zones 0.1–10 R☉ + 49 to the edge, fixed Eulerian) started from a 0.05 M☉ protostar (f128c stage-A3 checkpoint)",
    "z100s50": "From the IC, fixed 100-zone Eulerian mesh: 50 zones 0.1–10 R☉ (Δlog₁₀ r 0.040) + 49 zones to 166.7 pc (Δlog₁₀ r 0.181); z100ic3 settings",
    "z100e": "Fixed 100-zone log-Eulerian grid 0.1 R☉–143 pc (Δlog₁₀ r 0.110, no motion), production stack, from 0.05 M☉",
    "z100a": "Same fixed 100-zone grid from the end of stage A2 (just after the second collapse, 2e-4 M☉)",
    "t_tlh": "Short check of TL3RECKT=1 H2PPD=1 from production at 399 M☉ (control = t_c128c)",
    "t_gs_ic9": "Moving mesh GRIDSEG 32/64/32 (128 zones) from the IC, no remap",
    "t_gs_ic13": "Moving mesh GRIDSEG 24/52/24 (100 zones) from the IC, no remap",
    "t_te1": "TR-BDF2 error control (DTERR 0.05) + chemistry coupled in every zone, fixed 100-zone grid from 0.05 M☉",
    "t_tau3b": "Production with the τ ≥ 3 table layout but LOCALJ=0 (LOCALJ breaks the far cloud)",
    "t_nolte": "Candidate: table off everywhere outside n_H ≥ 1e16 (PSWTAU=0 + TL3DEEP), from production at 480 M☉",
    "t_nolte2": "PRODUCTION since 15 Nov: table only at n_H ≥ 1e16 (PSWTAU=0 + TL3DEEP + CPLTRACEX), from the old production at ~507 M☉",
    "z100ic3": "NATIVE: one run from the IC, fixed 100-zone log grid, no remap; TR-BDF2 error control (DTERR 0.1), chemistry coupled everywhere, NLOPAC=0",
    "t_cplwide": "Production with the coupled band extended over the warm front (CPLBAND r=8.5–1e10 R☉) + CPLTRACEX",
    "t_damp": "Production + DAMPNEWTON=1 (zone-wise damped Newton) against the warm-front T limit cycle",
    "f128c": "From the ABN02 IC with the production stack; phase B on the old 64-zone mesh failed at 0.40 M☉ (stopped)",
    "f128d": "From the IC, phase B at 0.05 M☉ with production settings + coupled stack (380 zones, ENVCONS)",
    "f128u": "Control for f128d: same, coupled solver off",
    "f128e": "From the IC, phase B at 0.05 M☉ with the repo-default settings (DTENV cap, CDYN=1e5) + coupled stack",
    "z100": "100-zone line from the initial conditions (moving mesh GRIDMAP=99, radius+mass arc length)",
    "z100h": "100-zone line (GRIDMAP=99) from base64 at 0.045 M☉",
    "t_ps0": "z100h with PSWTMIN=0 (H₂ dissociation layer on the equilibrium table)",
    "t_b224c": "b224 from 5.86 M☉ with the 100-zone line's chemistry speed-ups, 4 threads",
    "t_b224c8": "t_b224c on 8 threads, per-zone chemistry profile on",
    "t_xu1": "Test from t_b224c8 at 9.2 M☉: control with XADVU=1 (mismatch diagnostic only)",
    "t_xu2": "Test from 9.2 M☉: XADVU=2 (energy reset to hydro T after species advection)",
    "t_ec": "Test from 9.2 M☉: ECHEMSPLIT=1 (thermal/chemical energy split)",
    "t_gt2": "Moving-mesh test (GRIDMAP=99, GMTAU=0.4, PSWTMIN=0) from 0.045 M☉",
    "b224": "b224 baseline (224 zones, base64 from 0.045 M☉)",
    "base64": "base64 (64 zones, from the initial conditions, stopped at 1.78 M☉)",
    "nohd64": "nohd64 (base64 without HD)", "nohd64c": "nohd64 continuation",
    "f128": "f128 (full6 on 128 zones, from 14.3 M☉)",
    "full4": "full4 (radiation from log n_c 14)", "full3": "full3 (radiation from log n_c 18.5)",
    "full6": "full6 (2×10⁴ M☉ cloud, shock-only viscosity)", "c146": "c146 (146 zones, 106 to 141.5 M☉)",
    "fx3": "fx3 (286 zones, 102.7 to 105.9 M☉)", "full2": "full2 (valid below log n_c 14)", "full1": "full1",
    "full5": "full5", "ref": "reference history",
}
HIST_ORDER = ["zie1", "zie3", "zeulnc", "zeul25", "zeul5", "zbase3", "zeul1", "z200a", "zbase2", "zbase1", "b224", "base64", "nohd64", "f128", "full4", "full3", "full6", "c146", "fx3", "full2", "ref"]

def live_tags():
    tags = set()
    for e in os.listdir("/proc"):
        if not e.isdigit(): continue
        try:
            with open(f"/proc/{e}/comm") as f:
                if f.read().strip() != "julia": continue
            with open(f"/proc/{e}/environ", "rb") as f:
                for l in f.read().split(b"\0"):
                    if l.startswith(b"TAG=_"): tags.add(l[5:].decode(errors="replace"))
        except Exception: pass
    return tags

def num(x):
    try: return float(x)
    except Exception: return None

def perf_stats(fn):
    with open(fn) as f: lines = f.read().strip().split("\n")
    hdr = lines[0].split(","); rows = [l.split(",") for l in lines[1:] if l and not l.startswith("step")]
    rows = [r for r in rows if len(r) == len(hdr)]
    if not rows: return None
    ix = {h: i for i, h in enumerate(hdr)}; g = lambda r, k: num(r[ix[k]]) if k in ix else None
    last = rows[-1]
    # walk back about an hour of wall time (wall_s is cumulative within one process, pid column)
    wall = 0.0; j = len(rows) - 1
    while j > 0 and wall < 3600:
        a, b = rows[j - 1], rows[j]
        dw = g(b, "wall_s") - g(a, "wall_s") if a[ix["pid"]] == b[ix["pid"]] else g(b, "wall_s")
        wall += max(dw or 0, 0); j -= 1
    m0, m1, t0, t1 = g(rows[j], "mstar"), g(last, "mstar"), g(rows[j], "t_yr"), g(last, "t_yr")
    rate_h = (m1 - m0) / (wall / 3600) if wall > 60 and m0 is not None else None
    rate_yr = (m1 - m0) / (t1 - t0) if t1 and t0 and t1 > t0 else None
    # ndiv over the last 1000 steps (cumulative within a process)
    s_last = g(last, "step"); nd = 0; k = len(rows) - 1
    while k > 0 and s_last - g(rows[k - 1], "step") <= 1000:
        a, b = rows[k - 1], rows[k]
        d = g(b, "ndiv") - g(a, "ndiv") if a[ix["pid"]] == b[ix["pid"]] else g(b, "ndiv")
        nd += max(d or 0, 0); k -= 1
    rec = rows[-10:]
    def avg(key):
        v = [g(r, key) / max(g(r, "nsteps") or 100, 1) * 100 for r in rec if g(r, key) is not None]
        return round(sum(v) / len(v), 2) if v else None
    return dict(step=int(s_last), t_yr=g(last, "t_yr"), dt_s=g(last, "dt_s"), M_star=m1, R_star=g(last, "rstar"),
                zones=int(g(last, "nz") or 0), threads=int(g(last, "nthreads") or 0),
                rate_Msun_h=None if rate_h is None else round(rate_h, 4),
                rate_Msun_yr=None if rate_yr is None else float(f"{rate_yr:.3e}"), rate_window_s=round(wall),
                ndiv_1000=int(nd), chem_s_100=avg("t_chem"), hydro_s_100=avg("t_step"), lines_s_100=avg("t_lines"),
                dtlim=last[ix["dtlim"]] if "dtlim" in ix else None, perf_mtime=int(os.path.getmtime(fn)))

def main():
    tags = live_tags(); now = time.time()
    keys = set(d for d in os.listdir(P) if os.path.isdir(os.path.join(P, d)))
    keys |= set(t for t in tags if os.path.isdir(os.path.join(D, t)))
    keys |= set(t for t in os.listdir(D) if os.path.exists(os.path.join(D, t, "perf.csv")) and now - os.path.getmtime(os.path.join(D, t, "perf.csv")) < 86400 and not t.startswith("t_pf"))      # every run of the last 24 h (tests included)
    SKIP = set(open(os.path.join(os.path.dirname(P), "panel_skip.txt")).read().split()) if os.path.exists(os.path.join(os.path.dirname(P), "panel_skip.txt")) else set()
    keys -= SKIP      # pruned from the explorer (artifact file limit)
    out = []
    for k in sorted(keys):
        e = dict(key=k, desc=DESC.get(k, k), process=k in tags, frames=os.path.exists(os.path.join(P, k, "frames.json")))
        fn = os.path.join(D, k, "perf.csv")
        if os.path.exists(fn):
            try: e.update(perf_stats(fn) or {})
            except Exception as ex: e["err"] = str(ex)
        age = now - e.get("perf_mtime", 0)
        e["status"] = ("running" if age < 1800 else "stalled") if e["process"] else "stopped"
        out.append(e)
    def order(e):
        if e["status"] != "stopped": return (0, 0 if e["key"] == "t_nolte2" else (1 if e["key"] in ("zie1", "zie3", "zeulnc", "zeul25", "zeul5", "zbase3", "zeul1", "zbase2", "zbase1", "z300mm", "z100mm") else 2), -(e.get("M_star") or 0))
        return (1, HIST_ORDER.index(e["key"]) if e["key"] in HIST_ORDER else 99, e["key"])
    out.sort(key=order)
    doc = dict(updated=int(now), updated_iso=datetime.datetime.fromtimestamp(now).astimezone().isoformat(timespec="minutes"), runs=out)
    import math
    fin = lambda v: (None if isinstance(v, float) and not math.isfinite(v) else
                     {k: fin(x) for k, x in v.items()} if isinstance(v, dict) else [fin(x) for x in v] if isinstance(v, list) else v)
    doc = fin(doc)      # browsers' JSON.parse rejects NaN/Infinity (a NaN rate emptied the explorer's run list, 14 Nov)
    with open(os.path.join(P, "runs.json.tmp"), "w") as f: json.dump(doc, f, ensure_ascii=False, indent=1, allow_nan=False)
    os.replace(os.path.join(P, "runs.json.tmp"), os.path.join(P, "runs.json"))
    for e in out:
        if e["status"] != "stopped": print(e["status"], e["key"], e.get("M_star"), e.get("rate_Msun_h"), e.get("chem_s_100"))

main()
