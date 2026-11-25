#!/usr/bin/env python3
"""
plot_figL1_ensman_omukai.py

Publication-grade 4-panel (2x2) Figure L1:
  (a) Ensman (1994) Subcritical (6 km/s) Radiative Shock (T_gas, T_rad, rho/rho0)
  (b) Ensman (1994) Supercritical (20 km/s) Radiative Shock + Stahler-Shu-Taam precursor
  (c) Omukai & Nishi (1998) + Omukai (2001) Primordial Collapse T(n_H) across 22 decades & J_crit bifurcation
  (d) Ripamonti & Abel (2004) Emergent ExoMol H2 Ro-Vibrational Line & CIE Continuum Luminosities
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
OUT = "/home/tabel/.gemini/antigravity/brain/91f66bf7-1694-46da-a2b7-de91d0662618/figL1_ensman_and_omukai_collapse.png"

def load_csv(p):
    with open(p) as f:
        hdr = f.readline().strip().split(",")
    a = np.loadtxt(p, delimiter=",", skiprows=1)
    if a.ndim == 1:
        a = a[None, :]
    return {h: a[:, i] for i, h in enumerate(hdr)}

def main():
    sub = load_csv(os.path.join(DIR, "lm01_subc.csv"))
    sup = load_csv(os.path.join(DIR, "lm01_supc.csv"))
    j0  = load_csv(os.path.join(DIR, "lm02_collapse_J0.csv"))
    h2l = load_csv(os.path.join(DIR, "lm02_h2_lines_J0.csv"))

    plt.rcParams.update({
        "font.size": 9.5, "axes.titlesize": 10.5, "axes.labelsize": 9.5,
        "legend.fontsize": 7.8, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
    })
    fig, axs = plt.subplots(2, 2, figsize=(13.6, 9.2), dpi=180)
    fig.patch.set_facecolor("white")

    # (a) Ensman 1994 SUBC
    ax = axs[0, 0]
    x = (sub["r_cm"] - 8.0e11) / 1e10
    ax.plot(x, sub["T_gas"], color="#e03131", lw=2.2, label=r"$T_{\rm gas}$ (Zel'dovich spike)")
    ax.plot(x, sub["T_rad"], color="#1971c2", lw=2.0, ls="--", label=r"$T_{\rm rad} = (E_r / a)^{1/4}$")
    ax.set_xlabel(r"Position $r - r_0$ [$10^{10}$ cm]")
    ax.set_ylabel("Temperature [K]")
    ax.set_title(r"(a) Ensman (1994) Subcritical Shock ($v_p = 6\,$km s$^{-1}$)", fontweight="bold")
    ax2 = ax.twinx()
    ax2.plot(x, sub["M_comp"], color="#495057", lw=1.4, ls=":", label=r"$\rho / \rho_0$")
    ax2.set_ylabel(r"Compression $\rho / \rho_0$", color="#495057")
    ax2.set_ylim(0, max(6.0, float(np.max(sub["M_comp"])) * 1.15))
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper right", frameon=True, framealpha=0.92)

    # (b) Ensman 1994 SUPC
    ax = axs[0, 1]
    x = (sup["r_cm"] - 8.0e11) / 1e10
    ax.plot(x, sup["T_gas"], color="#e03131", lw=2.2, label=r"$T_{\rm gas}$ ($T_{\rm pre} \approx T_{\rm post}$)")
    ax.plot(x, sup["T_rad"], color="#1971c2", lw=2.0, ls="--", label=r"$T_{\rm rad}$ (radiative precursor)")
    ax.set_xlabel(r"Position $r - r_0$ [$10^{10}$ cm]")
    ax.set_ylabel("Temperature [K]")
    ax.set_title(r"(b) Ensman (1994) Supercritical Shock ($v_p = 20\,$km s$^{-1}$)", fontweight="bold")
    ax2 = ax.twinx()
    ax2.plot(x, sup["M_comp"], color="#495057", lw=1.4, ls=":", label=r"$\rho / \rho_0$")
    ax2.set_ylabel(r"Compression $\rho / \rho_0$", color="#495057")
    ax2.set_ylim(0, max(7.0, float(np.max(sup["M_comp"])) * 1.15))
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper right", frameon=True, framealpha=0.92)

    # (c) Omukai & Nishi (1998) + Omukai (2001) J_crit
    ax = axs[1, 0]
    curves = [
        ("J0",        r"$J_{21}=0$ (Pop III $\mathrm{H}_2$ track)",         "#1971c2", 2.4, "-"),
        ("J10",       r"$J_{21}=10$ ($T_{\rm rad}=10^4$ K)",                "#0ca678", 1.7, "-."),
        ("J30",       r"$J_{21}=30$ ($T_{\rm rad}=10^4$ K)",                "#f08c00", 1.8, "--"),
        ("J100",      r"$J_{21}=100$ ($T_{\rm rad}=10^4$ K)",               "#e8590c", 1.9, "-"),
        ("J1e3_T1e4", r"$J_{21}=10^3$ ($T_{\rm rad}=10^4$ K, DCBH)",        "#e03131", 2.2, "-"),
        ("J1e4_T1e5", r"$J_{21}=10^4$ ($T_{\rm rad}=10^5$ K, hard)",        "#7048e8", 1.8, ":"),
    ]
    for tag, lab, col, lw, ls in curves:
        p = os.path.join(DIR, f"lm02_collapse_{tag}.csv")
        if os.path.exists(p):
            c = load_csv(p)
            ax.plot(c["nH"], c["T"], color=col, lw=lw, ls=ls, label=lab)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1.0, 1e22); ax.set_ylim(80.0, 1.5e5)
    ax.set_xlabel(r"Central Hydrogen Number Density $n_{\rm H}$ [cm$^{-3}$]")
    ax.set_ylabel(r"Gas Temperature $T$ [K]")
    ax.set_title(r"(c) Omukai & Nishi (1998) Collapse & Omukai (2001) $J_{\rm crit}$", fontweight="bold")
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=7.5)
    # Physical regime annotations
    ax.annotate("Loitering\n(~200 K)", xy=(1e4, 210), xytext=(8e1, 110),
                arrowprops=dict(arrowstyle="->", color="#1971c2", lw=1.0), fontsize=7.6, color="#1971c2")
    ax.annotate("3-body H₂\nheating", xy=(1e10, 1000), xytext=(3e7, 2500),
                arrowprops=dict(arrowstyle="->", color="#1971c2", lw=1.0), fontsize=7.6, color="#1971c2")
    ax.annotate("Atomic H\nDCBH (8 kK)", xy=(1e5, 8000), xytext=(1e2, 2.2e4),
                arrowprops=dict(arrowstyle="->", color="#e03131", lw=1.0), fontsize=7.6, color="#e03131")
    ax.annotate("Hydrostatic\nCore", xy=(1e20, 2e4), xytext=(3e16, 4.5e4),
                arrowprops=dict(arrowstyle="->", color="#1971c2", lw=1.0), fontsize=7.6, color="#1971c2")

    # (d) Ripamonti & Abel (2004) Emergent H2 Lines & CIE Continuum
    ax = axs[1, 1]
    n = h2l["nH"]
    lcols = [
        ("L_H2tot_Lsun", r"Total ExoMol $\mathrm{H}_2$ lines", "#1971c2", 2.3, "-"),
        ("L_cont_Lsun",  r"Continuum (CIE + $\mathrm{H}^-$)", "#e03131", 2.1, "--"),
        ("L_S0_Lsun",    r"$\mathrm{H}_2$ 0–0 S(0) 28.2 $\mu$m", "#0ca678", 1.7, "-"),
        ("L_S1_Lsun",    r"$\mathrm{H}_2$ 0–0 S(1) 17.0 $\mu$m", "#2b8a3e", 1.7, "-."),
        ("L_S2_Lsun",    r"$\mathrm{H}_2$ 0–0 S(2) 12.3 $\mu$m", "#f08c00", 1.6, ":"),
        ("L_S3_Lsun",    r"$\mathrm{H}_2$ 0–0 S(3) 9.7 $\mu$m",  "#d9480f", 1.6, "-."),
        ("L_v1S1_Lsun",  r"$\mathrm{H}_2$ 1–0 S(1) 2.12 $\mu$m", "#7048e8", 1.8, "-"),
    ]
    for ck, lab, col, lw, ls in lcols:
        y = h2l[ck]
        m = y > 1e-6
        if np.any(m):
            ax.plot(n[m], y[m], color=col, lw=lw, ls=ls, label=lab)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1e2, 1e20); ax.set_ylim(1e-4, 5e3)
    ax.set_xlabel(r"Central Hydrogen Number Density $n_{\rm H}$ [cm$^{-3}$]")
    ax.set_ylabel(r"Jeans-Core Emergent Luminosity [$L_\odot$]")
    ax.set_title(r"(d) Ripamonti & Abel (2004) $\mathrm{H}_2$ Ro-Vib Lines & CIE", fontweight="bold")
    ax.legend(loc="lower left", frameon=True, framealpha=0.92, ncol=2, fontsize=7.3)

    fig.suptitle("Titan.jl Landmark Suite (Wave 1): Ensman (1994) Radiative Shocks & Omukai (1998/2001) / Ripamonti & Abel (2004) Collapse",
                 fontsize=11.8, fontweight="bold", y=0.985)
    plt.tight_layout(rect=[0, 0, 1, 0.962])
    plt.savefig(OUT, dpi=180)
    print(f"Wrote {OUT}")

if __name__ == "__main__":
    main()
