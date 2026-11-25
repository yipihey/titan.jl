#!/usr/bin/env python3
"""
plot_figL2_whalen_schaerer.py

Publication-grade 6-panel (2x3) Figure L2:
  Top row — Whalen, Abel & Norman (2004) + Kitayama+04 + Alvarez+06 Minihalo Photoevaporation:
    (a) Core evacuation & D-type shock expansion (n_H(r) & u(r) at 0, 0.25, 0.8, 2.5 Myr in 1e6 Msun halo)
    (b) Ionizing escape fraction f_esc(t) & bound gas fraction f_bound(t) across 3e5, 1e6, 5e6 Msun halos
    (c) Full-lifetime (2.5 Myr) JWST non-LTE atomic line lightcurves (Lyα, 2γ, Hα, Hβ, He II 1640, He I 10830)
  Bottom row — Schaerer (2002/2003) + Raiter+10 Pop III ZAMS & Nebular Grid:
    (d) Pop III ZAMS HR & Ionizing photon rates (Q_H, Q_HeII) vs Schaerer (2002) Table 3 (5 -> 500 Msun)
    (e) Rest-frame Equivalent Widths EW(Lyα), EW(He II 1640), EW(Hβ) vs M_*
    (f) Raiter+10 high-density nebular diagnostics (n_e = 1e2 -> 1e8 cm⁻³: 2γ suppression, He I 10830 boost, Hα/Hβ)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
OUT = "/home/tabel/.gemini/antigravity/brain/91f66bf7-1694-46da-a2b7-de91d0662618/figL2_whalen04_and_schaerer02.png"

def load_csv(p):
    with open(p) as f:
        hdr = f.readline().strip().split(",")
    a = np.loadtxt(p, delimiter=",", skiprows=1)
    if a.ndim == 1:
        a = a[None, :]
    return {h: a[:, i] for i, h in enumerate(hdr)}

def main():
    h3e5 = load_csv(os.path.join(DIR, "lm03_halo_3e5.csv"))
    h1e6 = load_csv(os.path.join(DIR, "lm03_halo_1e6.csv"))
    h5e6 = load_csv(os.path.join(DIR, "lm03_halo_5e6.csv"))
    zams = load_csv(os.path.join(DIR, "lm04_schaerer_zams.csv"))
    dens = load_csv(os.path.join(DIR, "lm04_raiter_density.csv"))

    plt.rcParams.update({
        "font.size": 9.2, "axes.titlesize": 10.2, "axes.labelsize": 9.2,
        "legend.fontsize": 7.5, "xtick.labelsize": 8.4, "ytick.labelsize": 8.4,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
    })
    fig, axs = plt.subplots(2, 3, figsize=(15.5, 9.2), dpi=180)
    fig.patch.set_facecolor("white")

    # (a) Whalen+04 1e6 Msun radial shock & evacuation
    ax = axs[0, 0]; ax2 = ax.twinx()
    profs = [("0Myr", "0 Myr", "#495057"), ("0p25Myr", "0.25 Myr", "#1971c2"),
             ("0p8Myr", "0.8 Myr", "#0ca678"), ("2p5Myr", "2.5 Myr", "#e03131")]
    for suf, lab, col in profs:
        p = os.path.join(DIR, f"lm03_prof_1e6_{suf}.csv")
        if os.path.exists(p):
            d = load_csv(p)
            ax.plot(d["r_pc"], d["nH"], color=col, lw=2.0, label=lab)
            ax2.plot(d["r_pc"], d["u_kms"], color=col, lw=1.3, ls="--", alpha=0.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.35, 150); ax.set_ylim(2e-2, 5e4)
    ax2.set_ylim(-5, 42)
    ax.set_xlabel("Radius $r$ [pc]")
    ax.set_ylabel(r"Density $n_{\rm H}$ [cm$^{-3}$] (solid)")
    ax2.set_ylabel(r"Velocity $u$ [km s$^{-1}$] (dashed)", color="#555")
    ax.legend(loc="upper right", frameon=True, framealpha=0.92, title=r"$10^6\,M_\odot$ Halo")
    ax.set_title("(a) Whalen+04 Minihalo D-Type Shock & Evacuation", fontweight="bold")

    # (b) Alvarez+06 f_esc(t) & Kitayama+04 f_bound(t)
    ax = axs[0, 1]
    for d, lab, col in [(h3e5, r"$3\times10^5\,M_\odot$", "#0ca678"),
                        (h1e6, r"$1\times10^6\,M_\odot$", "#1971c2"),
                        (h5e6, r"$5\times10^6\,M_\odot$", "#e03131")]:
        ax.plot(d["t_Myr"], d["f_esc_H"], color=col, lw=2.1, label=lab + r" ($f_{\rm esc,H}$)")
        ax.plot(d["t_Myr"], d["f_bound"], color=col, lw=1.6, ls="--", alpha=0.85, label=lab + r" ($f_{\rm bound}$)")
    ax.set_xlim(0, 2.5); ax.set_ylim(-0.03, 1.08)
    ax.set_xlabel("Time $t$ [Myr]")
    ax.set_ylabel(r"$f_{\rm esc}$ (solid) & Bound Gas $f_{\rm bound}$ (dashed)")
    ax.legend(loc="center right", frameon=True, framealpha=0.92, ncol=1, fontsize=7.0)
    ax.set_title(r"(b) Alvarez+06 $f_{\rm esc}(t)$ & Kitayama+04 Unbinding", fontweight="bold")

    # (c) Full-lifetime 2.5 Myr JWST line lightcurves (1e6 Msun halo)
    ax = axs[0, 2]
    for ck, lab, col, ls in [("L_Lya", r"Ly$\alpha$", "#e03131", "-"),
                             ("L_2g",  r"$2\gamma$ cont.", "#f08c00", "--"),
                             ("L_Ha",  r"H$\alpha$", "#c2255c", "-"),
                             ("L_Hb",  r"H$\beta$", "#1971c2", "-"),
                             ("L_1640",r"He II $\lambda 1640$", "#5f3dc4", "-"),
                             ("L_10830",r"He I $\lambda 10830$", "#0ca678", "-."),
                             ("L_4686",r"He II $\lambda 4686$", "#364fc7", ":")]:
        y = h1e6[ck]
        m = y > 1.0
        if np.any(m):
            ax.plot(h1e6["t_Myr"][m], y[m], color=col, lw=1.9, ls=ls, label=lab)
    ax.set_yscale("log"); ax.set_xlim(0, 2.5); ax.set_ylim(1e1, 5e5)
    ax.set_xlabel("Time $t$ [Myr]")
    ax.set_ylabel(r"Halo Line Luminosity [$L_\odot$]")
    ax.legend(loc="upper right", frameon=True, framealpha=0.92, ncol=2, fontsize=7.2)
    ax.set_title(r"(c) 2.5 Myr Nebular Line Fading ($10^6\,M_\odot$ Halo)", fontweight="bold")

    # (d) Schaerer (2002) Table 3 comparison: Q_H, Q_HeII & L vs M_*
    ax = axs[1, 0]; ax2 = ax.twinx()
    M = zams["M_Msun"]
    ax.plot(M, zams["logQH"], color="#e03131", lw=2.2, marker="o", ms=4.5, label=r"Titan.jl $\log Q_{\rm H}$")
    ax.plot(M, zams["logQH_S02"], color="#e03131", lw=1.4, ls="--", marker="x", ms=5.0, label=r"Schaerer+02 $\log Q_{\rm H}$")
    ax.plot(M, zams["logQHe2"], color="#5f3dc4", lw=2.2, marker="s", ms=4.5, label=r"Titan.jl $\log Q_{\rm HeII}$")
    ax.plot(M, zams["logQHe2_S02"], color="#5f3dc4", lw=1.4, ls="--", marker="x", ms=5.0, label=r"Schaerer+02 $\log Q_{\rm HeII}$")
    ax2.plot(M, zams["Teff"] / 1e3, color="#0ca678", lw=1.8, ls="-.", label=r"$T_{\rm eff}$ [kK] (Titan)")
    ax2.plot(M, zams["Teff_S02"] / 1e3, color="#0ca678", lw=1.2, ls=":", marker="+", ms=5, label=r"$T_{\rm eff}$ (S02)")
    ax.set_xscale("log"); ax.set_xlim(4.5, 600); ax.set_ylim(38, 52)
    ax2.set_ylim(20, 125)
    ax.set_xlabel(r"Pop III ZAMS Mass $M_*$ [$M_\odot$]")
    ax.set_ylabel(r"$\log_{10} Q_i$ [photons s$^{-1}$]")
    ax2.set_ylabel(r"$T_{\rm eff}$ [kK]", color="#0ca678")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="lower right", frameon=True, framealpha=0.92, ncol=2, fontsize=6.7)
    ax.set_title("(d) Schaerer (2002) Pop III ZAMS Ionizing Grid", fontweight="bold")

    # (e) Rest-frame EWs vs Mass
    ax = axs[1, 1]
    ax.plot(M, np.maximum(zams["EW_Lya_A"], 1.0), color="#e03131", lw=2.2, marker="o", ms=4.5, label=r"EW(Ly$\alpha$)")
    ax.plot(M, np.maximum(zams["EW_Hb_A"], 1.0),  color="#1971c2", lw=2.1, marker="s", ms=4.5, label=r"EW(H$\beta$)")
    ax.plot(M, np.maximum(zams["EW_1640_A"], 0.1),color="#5f3dc4", lw=2.2, marker="^", ms=5.0, label=r"EW(He II $\lambda 1640$)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(4.5, 600); ax.set_ylim(0.5, 5e3)
    ax.set_xlabel(r"Pop III ZAMS Mass $M_*$ [$M_\odot$]")
    ax.set_ylabel(r"Rest-Frame Equivalent Width EW$_0$ [Å]")
    ax.legend(loc="lower right", frameon=True, framealpha=0.92)
    ax.set_title("(e) Pop III Nebular Equivalent Widths (incl. 2γ Cont.)", fontweight="bold")

    # (f) Raiter+10 high-density diagnostics
    ax = axs[1, 2]; ax2 = ax.twinx()
    ne = dens["ne"]
    ax.plot(ne, dens["HeI10830_Hb"], color="#0ca678", lw=2.2, label=r"He I $\lambda 10830$ / H$\beta$")
    ax.plot(ne, dens["Ha_Hb"],       color="#c2255c", lw=2.1, label=r"H$\alpha$ / H$\beta$")
    ax.plot(ne, dens["HeII1640_Hb"], color="#5f3dc4", lw=2.1, label=r"He II $\lambda 1640$ / H$\beta$")
    ax.plot(ne, dens["HeI5876_Hb"],  color="#2b8a3e", lw=1.7, ls="-.", label=r"He I $\lambda 5876$ / H$\beta$")
    ax2.plot(ne, dens["L_2g_Lsun"] / (dens["L_Lya_Lsun"] + dens["L_2g_Lsun"]), color="#f08c00", lw=1.9, ls="--",
             label=r"$2\gamma / (\mathrm{Ly}\alpha + 2\gamma)$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1e2, 1e8); ax.set_ylim(8e-2, 15.0)
    ax2.set_ylim(0, 0.42)
    ax.set_xlabel(r"Nebular Electron Density $n_e$ [cm$^{-3}$]")
    ax.set_ylabel(r"Line Ratio $L_{\rm line} / L_{\mathrm{H}\beta}$")
    ax2.set_ylabel(r"$2\gamma$ Continuum Fraction", color="#f08c00")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", frameon=True, framealpha=0.92, ncol=2, fontsize=6.8)
    ax.set_title("(f) Raiter+10 High-Density $2\gamma$ & He I Diagnostics", fontweight="bold")

    fig.suptitle("Titan.jl Landmark Suite (Wave 2): Whalen+04 / Alvarez+06 Minihalo Photoevaporation & Schaerer+02 / Raiter+10 Pop III Grid",
                 fontsize=11.6, fontweight="bold", y=0.985)
    plt.tight_layout(rect=[0, 0, 1, 0.962])
    plt.savefig(OUT, dpi=180)
    print(f"Wrote {OUT}")

if __name__ == "__main__":
    main()
