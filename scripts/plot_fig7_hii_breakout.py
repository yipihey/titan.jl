#!/usr/bin/env python3
"""
plot_fig7_hii_breakout.py

Generate publication-grade 6-panel (2x3) Figure 7:
  Live Dynamical 1D RHD H II / He III / LW Region Breakout, Champagne Shock
  Accretion Reversal, and Time-Resolved l-Resolved Non-LTE Atomic Line Signatures.
"""
import glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

RUNS = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs"
PDIR = os.path.join(RUNS, "profiles")
OUT  = "/home/tabel/.gemini/antigravity/brain/91f66bf7-1694-46da-a2b7-de91d0662618/fig7_dynamical_hii_breakout_and_lines.png"

def load_csv(path):
    with open(path) as f:
        hdr = f.readline().strip().split(",")
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    if arr.ndim == 1:
        arr = arr[None, :]
    return {h: arr[:, i] for i, h in enumerate(hdr)}

def load_profiles(prefix="hii_breakout_rhd_s"):
    profs = []
    for fn in sorted(glob.glob(os.path.join(PDIR, f"{prefix}*.csv"))):
        m = re.search(r"_s(\d+)_t([\d\.]+)\.csv$", os.path.basename(fn))
        if m:
            profs.append((int(m.group(1)), float(m.group(2)), load_csv(fn)))
    return sorted(profs, key=lambda x: x[1])

def main():
    d = load_csv(os.path.join(RUNS, "hii_breakout_rhd.csv"))
    trapped_p = os.path.join(RUNS, "hii_trapped_1e3.csv")
    if not os.path.exists(trapped_p):
        trapped_p = os.path.join(RUNS, "hii_trapped_rhd.csv")
    dt = load_csv(trapped_p) if os.path.exists(trapped_p) else None
    dt3 = load_csv(os.path.join(RUNS, "hii_trapped_rhd.csv")) if os.path.exists(os.path.join(RUNS, "hii_trapped_rhd.csv")) else None
    profs = load_profiles()

    t = d["t_yr"]
    tmax = max(t[-1], 600.0)

    plt.rcParams.update({
        "font.size": 9.5, "axes.titlesize": 10.5, "axes.labelsize": 9.5,
        "legend.fontsize": 7.8, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
    })
    fig, axs = plt.subplots(2, 3, figsize=(15.5, 9.2), dpi=180)
    fig.patch.set_facecolor("white")

    # ── Panel (a): Front Trajectories (Breakout vs Trapped) ──
    ax = axs[0, 0]
    ax.axvspan(150, 600, color="#fff3bf", alpha=0.55, label=r"KH contraction ($5\to 95.5$ kK)")
    ax.plot(t, np.maximum(d["R_LW_pc"], 2e-3), color="#5f3dc4", lw=2.0, label=r"$R_{\rm LW}$ ($\mathrm{H}_2$ dissoc., breakout)")
    m_hii = d["R_HII_pc"] > 0
    if np.any(m_hii):
        ax.plot(t[m_hii], d["R_HII_pc"][m_hii], color="#e03131", lw=2.3, label=r"$R_{\rm HII}$ (breakout, $\dot{M}_0=1.5\times10^{-4}$)")
    m_he3 = d["R_HeIII_pc"] > 0
    if np.any(m_he3):
        ax.plot(t[m_he3], d["R_HeIII_pc"][m_he3], color="#1971c2", lw=2.1, label=r"$R_{\rm HeIII}$ (nested He III)")
    m_sh = (d["R_shock_pc"] > 0) & (d["R_shock_pc"] < 0.35)
    if np.any(m_sh):
        ax.plot(t[m_sh], d["R_shock_pc"][m_sh], color="#d9480f", lw=1.9, ls="-.", label=r"$R_{\rm shock}$ (core champagne)")
    m_sh2 = d["R_shock_pc"] >= 0.35
    if np.any(m_sh2):
        ax.scatter(t[m_sh2], d["R_shock_pc"][m_sh2], color="#d9480f", s=24, marker="^", zorder=5, label=r"Outer D-type H II shock")
    if dt is not None:
        m_th = dt["R_HII_pc"] > 0
        if np.any(m_th):
            ax.plot(dt["t_yr"][m_th], dt["R_HII_pc"][m_th], color="#868e96", lw=1.7, ls="--",
                    label=r"$R_{\rm HII}$ trapped ($\dot{M}_0=10^{-3}$ OI02)")
    ax.set_yscale("log")
    ax.set_xlim(0, tmax)
    ax.set_ylim(1.1e-3, 1.9)
    ax.set_xlabel("Time $t$ [yr]")
    ax.set_ylabel("Radius [pc]")
    ax.set_title("(a) H II / He III / LW Breakout vs Trapped Regime", fontweight="bold")
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, ncol=1)

    # ── Panel (b): Dynamical Accretion Reversal & Champagne Outflow Velocity ──
    ax = axs[0, 1]
    ax.axvspan(150, 600, color="#fff3bf", alpha=0.55)
    ax.axhline(0.0, color="#343a40", lw=1.0, ls=":")
    ax.plot(t, d["Mdot_in_Msunyr"], color="#1971c2", lw=2.0, label=r"$\dot{M}(r_{\rm in})$ (inner sink)")
    ax.plot(t, d["Mdot_0p02pc_Msunyr"], color="#e03131", lw=2.2, label=r"$\dot{M}(0.02\,$pc) (envelope)")
    if dt is not None:
        ax.plot(dt["t_yr"], dt["Mdot_0p02pc_Msunyr"], color="#868e96", lw=1.5, ls="--", label=r"Trapped $\dot{M}$ ($10^{-3}$)")
    ax.set_yscale("symlog", linthresh=5e-5)
    ax.set_yticks([-1e-3, -1e-4, 0.0, 1e-4, 5e-4])
    ax.set_yticklabels([r"$-10^{-3}$", r"$-10^{-4}$", "0", r"$+10^{-4}$", r"$+5\times10^{-4}$"])
    ax.set_xlim(0, tmax)
    ax.set_xlabel("Time $t$ [yr]")
    ax.set_ylabel(r"Mass Flux $\dot{M}$ [$M_\odot\,\mathrm{yr}^{-1}$] (+ outflow / $-$ infall)")
    ax2 = ax.twinx()
    ax2.plot(t, d["u_max_kms"], color="#2b8a3e", lw=2.0, ls="-.", label=r"$u_{\max}$ [km/s]")
    ax2.set_ylabel(r"Peak Outflow Speed $u_{\max}$ [km s$^{-1}$]", color="#2b8a3e")
    ax2.tick_params(axis="y", colors="#2b8a3e")
    ax2.set_ylim(-5, max(35.0, float(np.max(d["u_max_kms"])) * 1.18))
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", frameon=True, framealpha=0.92)
    ax.set_title("(b) Accretion Reversal & Champagne Outflow", fontweight="bold")

    # ── Panel (c): Radial Champagne Shock Evolution u(r) & n_H(r) ──
    ax = axs[0, 2]
    ax2 = ax.twinx()
    if profs:
        # Select 4-5 informative epochs spanning initial infall -> breakout -> champagne shock
        targets = [0.0, 450.0, 950.0, 2200.0, profs[-1][1]]
        chosen = []
        seen = set()
        for tt in targets:
            bp = min(profs, key=lambda p: abs(p[1] - tt))
            if bp[0] not in seen:
                seen.add(bp[0])
                chosen.append(bp)
        cols = ["#495057", "#1971c2", "#2b8a3e", "#e67700", "#c92a2a"]
        for ci, (st, ty, pf) in enumerate(chosen):
            col = cols[ci % len(cols)]
            r_pc = pf["r_pc"]
            ax.plot(r_pc, pf["u_kms"], color=col, lw=2.0, label=f"$t={ty:.0f}$ yr")
            ax2.plot(r_pc, pf["nH"], color=col, lw=1.2, ls="--", alpha=0.72)
    ax.axhline(0.0, color="#343a40", lw=0.9, ls=":")
    ax.set_xscale("log")
    ax2.set_yscale("log")
    ax.set_xlim(2.2e-3, 1.25)
    ax.set_xlabel("Radius $r$ [pc]")
    ax.set_ylabel(r"Radial Velocity $u(r)$ [km s$^{-1}$] (solid)")
    ax2.set_ylabel(r"Density $n_{\rm H}(r)$ [cm$^{-3}$] (dashed)", color="#555")
    ax2.set_ylim(1e1, 3e7)
    ax.legend(loc="upper right", frameon=True, framealpha=0.92, title="Epoch")
    ax.set_title("(c) Swept-Up Shock Shell & Core Evacuation", fontweight="bold")

    # ── Panel (d): Nested Thermal & Chemical Structure (Final Snapshot) ──
    ax = axs[1, 0]
    if profs:
        st_f, ty_f, pf = profs[-1]
        r_pc = pf["r_pc"]
        ax.plot(r_pc, np.maximum(pf["xHp"], 1e-8), color="#e03131", lw=2.3, label=r"$x_{\mathrm{H}^+}$ (H II)")
        ax.plot(r_pc, np.maximum(pf["xHepp"] / 0.0833, 1e-8), color="#1971c2", lw=2.1, label=r"$x_{\mathrm{He}^{++}} / y_{\rm He}$ (He III)")
        ax.plot(r_pc, np.maximum(pf["xHep"] / 0.0833, 1e-8), color="#0ca678", lw=1.6, ls=":", label=r"$x_{\mathrm{He}^+} / y_{\rm He}$ (He II)")
        # Also show H2 front during breakout (t ≈ 400 yr) so LW dissociation front is clearly visible
        p_lw = min(profs, key=lambda p: abs(p[1] - 400.0))
        ax.plot(p_lw[2]["r_pc"], np.maximum(p_lw[2]["xH2"], 1e-10), color="#7048e8", lw=1.8, ls="--",
                label=rf"$x_{{\mathrm{{H}}_2}}$ ($t={p_lw[1]:.0f}$ yr LW)")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(2.2e-3, 1.25); ax.set_ylim(1e-7, 2.0)
        ax.set_xlabel("Radius $r$ [pc]")
        ax.set_ylabel("Species Fractional Abundance")
        ax2 = ax.twinx()
        ax2.plot(r_pc, pf["T"], color="#d9480f", lw=1.9, alpha=0.85, label=r"$T_{\rm gas}$ [K]")
        ax2.set_yscale("log"); ax2.set_ylim(1e2, 8e4)
        ax2.set_ylabel("Gas Temperature $T$ [K]", color="#d9480f")
        ax2.tick_params(axis="y", colors="#d9480f")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="lower left", frameon=True, framealpha=0.92, ncol=2)
        ax.set_title(f"(d) Nested He III → H II → LW PDR ($t={ty_f:.0f}$ yr)", fontweight="bold")

    # ── Panel (e): Time-Resolved l-Resolved Non-LTE Atomic Line Lightcurves ──
    ax = axs[1, 1]
    ax.axvspan(150, 600, color="#fff3bf", alpha=0.55)
    lines = [
        ("L_Lya_Lsun",     r"H I Ly$\alpha$ (1216 Å)",    "#e03131", 2.2, "-"),
        ("L_2g_Lsun",      r"H I $2\gamma$ cont.",        "#f08c00", 1.8, "--"),
        ("L_Ha_Lsun",      r"H I H$\alpha$ (6563 Å)",     "#c2255c", 2.0, "-"),
        ("L_Hb_Lsun",      r"H I H$\beta$ (4861 Å)",      "#1971c2", 1.9, "-"),
        ("L_HeI10830_Lsun",r"He I $\lambda 10830$",       "#0ca678", 1.8, "-."),
        ("L_HeII1640_Lsun",r"He II $\lambda 1640$",       "#5f3dc4", 2.2, "-"),
        ("L_HeI5876_Lsun", r"He I $\lambda 5876$",        "#2b8a3e", 1.5, ":"),
        ("L_HeII4686_Lsun",r"He II $\lambda 4686$",       "#364fc7", 1.7, "--"),
        ("L_Paa_Lsun",     r"H I Pa$\alpha$ (1.87 $\mu$m)","#862e9c", 1.5, ":"),
    ]
    for colk, lab, c, lw, ls in lines:
        y = d[colk]
        m = y > 1.0
        if np.any(m):
            ax.plot(t[m], y[m], color=c, lw=lw, ls=ls, label=lab)
    ax.set_yscale("log")
    ax.set_xlim(200, tmax)
    ax.set_ylim(10.0, 6e5)
    ax.set_xlabel("Time $t$ [yr]")
    ax.set_ylabel(r"Line Luminosity $L_{\rm line}$ [$L_\odot$]")
    ax.set_title(r"(e) Non-LTE Atomic Line Lightcurves ($l$-Resolved)", fontweight="bold")
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, ncol=2, fontsize=7.2)

    # ── Panel (f): Diagnostic Pop III Line Ratios vs Time ──
    ax = axs[1, 2]
    ax.axvspan(150, 600, color="#fff3bf", alpha=0.55)
    m = d["L_Hb_Lsun"] > 50.0
    if np.any(m):
        hb = d["L_Hb_Lsun"][m]
        tm = t[m]
        ax.plot(tm, d["L_Ha_Lsun"][m] / hb, color="#c2255c", lw=2.2, label=r"H$\alpha$ / H$\beta$")
        ax.plot(tm, d["L_HeI10830_Lsun"][m] / hb, color="#0ca678", lw=2.0, label=r"He I $\lambda 10830$ / H$\beta$")
        ax.plot(tm, d["L_HeII1640_Lsun"][m] / hb, color="#5f3dc4", lw=2.2, label=r"He II $\lambda 1640$ / H$\beta$")
        ax.plot(tm, d["L_HeI5876_Lsun"][m] / hb, color="#2b8a3e", lw=1.7, ls="-.", label=r"He I $\lambda 5876$ / H$\beta$")
        ax.plot(tm, d["L_HeII4686_Lsun"][m] / hb, color="#364fc7", lw=1.8, ls="--", label=r"He II $\lambda 4686$ / H$\beta$")
    if dt3 is not None:
        mt = dt3["L_Hb_Lsun"] > 50.0
        if np.any(mt):
            r1640_trap = float(dt3["L_HeII1640_Lsun"][mt][-1] / dt3["L_Hb_Lsun"][mt][-1])
            r10830_trap = float(dt3["L_HeI10830_Lsun"][mt][-1] / dt3["L_Hb_Lsun"][mt][-1])
            ax.axhline(r1640_trap, color="#5f3dc4", lw=1.1, ls=":", alpha=0.75)
            ax.axhline(r10830_trap, color="#0ca678", lw=1.1, ls=":", alpha=0.75)
            ax.text(tmax * 0.97, r10830_trap * 1.10, f"trapped He I 10830/Hβ = {r10830_trap:.2f}",
                    ha="right", va="bottom", fontsize=7.2, color="#0ca678")
            ax.text(tmax * 0.97, r1640_trap * 0.88, f"trapped He II 1640/Hβ = {r1640_trap:.2f}",
                    ha="right", va="top", fontsize=7.2, color="#5f3dc4")
    ax.set_yscale("log")
    ax.set_xlim(300, tmax)
    ax.set_ylim(2e-2, 12.0)
    ax.set_xlabel("Time $t$ [yr]")
    ax.set_ylabel(r"Line Luminosity Ratio ($L_{\rm line} / L_{\mathrm{H}\beta}$)")
    ax.set_title("(f) JWST/Pop III Diagnostic Ratios (Breakout vs Trapped)", fontweight="bold")
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, ncol=2, fontsize=7.4)

    fig.suptitle(
        r"Titan.jl 1D RHD + Option A Multigroup IonRay + $l$-Resolved ModelAtom: $100\,M_\odot$ Pop III H II / He III / LW Breakout & Accretion Reversal",
        fontsize=12.2, fontweight="bold", y=0.985
    )
    plt.tight_layout(rect=[0, 0, 1, 0.962])
    plt.savefig(OUT, dpi=180)
    print(f"Wrote {OUT}")

if __name__ == "__main__":
    main()
