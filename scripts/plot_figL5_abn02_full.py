#!/usr/bin/env python3
"""
plot_figL5_abn02_full.py

Publication-grade 6-panel (2x3) Figure L5:
  Abel, Bryan & Norman (2002, Science 295, 93) End-to-End First Star Evolution:
    (a) Stage 1 Validation: 1D RHD Collapse from abn02_red -> Direct Overlay vs Published abn02_final (n_c = 1e12) + Core Birth (1e21)
    (b) Stage 1 Thermal & Chemical Structure T(M_enc) & f_H2(M_enc) vs abn02_final Fig. 2
    (c) Stage 1 -> 2 Bridge: Exact Extracted ABN02 Cosmological Accretion Rate Mdot_ABN02(M_*) vs OI02 Breakout Thresholds
    (d) Stage 2 Protostellar Growth: R_*(M_*), T_eff(M_*), Q_H(M_*) from Core Birth (0.02 Msun) -> Swelling -> ZAMS -> Breakout
    (e) Stage 3 Full-Lifetime ABN02 Minihalo Photoevaporation: R_HII(t), f_esc(t), f_bound(t) (Case B Polar Disk vs Case A Strict 1D)
    (f) Stage 3 Time-Resolved JWST Atomic Line Signatures (Lyα, Hβ, He II 1640, He I 10830) over 2.2-3.2 Myr Lifetime
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
ABN_RED = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/wt/lib/popiii/abn02red/abn02_red_fig2.csv"
ABN_FIN = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/wt/lib/popiii/abn02fin/abn02_final_fig2.csv"
OUT = "/home/tabel/.gemini/antigravity/brain/91f66bf7-1694-46da-a2b7-de91d0662618/figL5_abn02_end_to_end.png"

def load_csv(p):
    with open(p) as f:
        hdr = f.readline().strip().split(",")
    a = np.loadtxt(p, delimiter=",", skiprows=1)
    if a.ndim == 1:
        a = a[None, :]
    return {h: a[:, i] for i, h in enumerate(hdr)}

def load_abn(p):
    tab = {}
    with open(p) as f:
        for l in f:
            if l.startswith("#") or l.startswith("curve") or not l.strip():
                continue
            k, x, y = l.strip().split(",")[:3]
            tab.setdefault(k, ([], []))
            tab[k][0].append(float(x)); tab[k][1].append(float(y))
    return {k: (np.array(v[0]), np.array(v[1])) for k, v in tab.items()}

def main():
    red = load_abn(ABN_RED)
    fin = load_abn(ABN_FIN)
    s12 = load_csv(os.path.join(DIR, "lm08_abn02_at_1e12.csv"))
    sc  = load_csv(os.path.join(DIR, "lm08_abn02_core_birth.csv"))
    pr  = load_csv(os.path.join(DIR, "lm08_abn02_protostar_track.csv"))
    pp  = load_csv(os.path.join(DIR, "lm08_abn02_photoevap_polar.csv"))
    p1  = load_csv(os.path.join(DIR, "lm08_abn02_photoevap_1d.csv"))

    plt.rcParams.update({
        "font.size": 9.2, "axes.titlesize": 10.1, "axes.labelsize": 9.2,
        "legend.fontsize": 7.3, "xtick.labelsize": 8.3, "ytick.labelsize": 8.3,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
    })
    fig, axs = plt.subplots(2, 3, figsize=(15.6, 9.3), dpi=180)
    fig.patch.set_facecolor("white")

    # (a) Radial Density n(r): abn02_red IC -> Titan at 1e12 vs abn02_final -> Hydrostatic Core Birth (1e21)
    ax = axs[0, 0]
    ax.plot(10**red["n_r"][0], 10**red["n_r"][1], color="#868e96", lw=1.8, ls=":", label=r"ABN02 Fig. 2 red IC ($9$ Myr)")
    ax.plot(10**fin["n_r"][0], 10**fin["n_r"][1], color="#111111", lw=2.2, ls="--", label=r"ABN02 Fig. 2 final (Enzo 3D)")
    ax.plot(s12["r_pc"], s12["nH"], color="#1971c2", lw=2.1, label=r"Titan.jl 1D RHD @ $n_c=10^{12}$")
    ax.plot(sc["r_pc"],  sc["nH"],  color="#e03131", lw=2.2, label=r"Titan.jl Core Birth ($n_c\sim10^{21}$)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1e-8, 30.0); ax.set_ylim(1e1, 3e21)
    ax.set_xlabel("Radius $r$ [pc]")
    ax.set_ylabel(r"Hydrogen Density $n_{\rm H}$ [cm$^{-3}$]")
    ax.legend(loc="lower left", frameon=True, framealpha=0.92)
    ax.set_title("(a) Stage 1: ABN02 Collapse to Hydrostatic Core", fontweight="bold")

    # (b) T(M_enc) & f_H2(M_enc) Direct Validation vs abn02_final Fig. 2
    ax = axs[0, 1]; ax2 = ax.twinx()
    ax.plot(10**fin["T_M"][0], 10**fin["T_M"][1], color="#111111", lw=2.0, ls="--", label=r"$T(M)$ ABN02 final")
    ax.plot(s12["m_Msun"], s12["T"], color="#e03131", lw=2.1, label=r"$T(M)$ Titan @ $10^{12}$")
    ax.plot(sc["m_Msun"],  sc["T"],  color="#d9480f", lw=1.7, ls="-.", label=r"$T(M)$ Core Birth")
    ax2.plot(10**fin["fH2_M"][0], 10**fin["fH2_M"][1], color="#495057", lw=1.5, ls=":", label=r"$f_{\mathrm{H}_2}$ ABN02 final")
    ax2.plot(s12["m_Msun"], np.maximum(s12["fH2"], 1e-5), color="#1971c2", lw=1.8, label=r"$f_{\mathrm{H}_2}$ Titan @ $10^{12}$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax2.set_yscale("log")
    ax.set_xlim(1e-3, 5e3); ax.set_ylim(100, 4e4); ax2.set_ylim(1e-4, 1.5)
    ax.set_xlabel(r"Enclosed Gas Mass $M_{\rm enc}$ [$M_\odot$]")
    ax.set_ylabel("Temperature $T$ [K]", color="#e03131")
    ax2.set_ylabel(r"$\mathrm{H}_2$ Mass Fraction $f_{\mathrm{H}_2}$", color="#1971c2")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper right", frameon=True, framealpha=0.92, fontsize=6.9)
    ax.set_title(r"(b) Stage 1: $T(M)$ & $f_{\mathrm{H}_2}(M)$ vs ABN02 Fig. 2", fontweight="bold")

    # (c) Exact ABN02 Accretion Rate Mdot(M_*) vs Omukai & Inutsuka (2002) Breakout Thresholds
    ax = axs[0, 2]
    M = pr["M_Msun"]
    ax.plot(M, pr["Mdot_1d"],   color="#111111", lw=2.3, label=r"Strict 1D spherical $\dot{M}_{\rm 1D}(M_*)$")
    ax.plot(M, pr["Mdot_disk"], color="#0ca678", lw=2.2, ls="--", label=r"3D Disk-regulated $\dot{M}_{\rm disk}$ ($v_{\rm rot}/v_{\rm Kep}\approx 0.4$)")
    ax.plot(M, pr["Mdot_crit_disk"], color="#1971c2", lw=2.0, label=r"OI02 breakout $\dot{M}_{\rm crit}$ (Disk track)")
    ax.plot(M, pr["Mdot_crit_1d"],   color="#e03131", lw=1.9, ls="-.", label=r"OI02 breakout $\dot{M}_{\rm crit}$ (1D track)")
    ax.axhline(4.2e-3, color="#7048e8", lw=1.1, ls=":", label=r"OP03 supergiant threshold ($4.2\times10^{-3}$)")
    i_pol = np.where((pr["Mdot_crit_disk"] > 0.20 * pr["Mdot_disk"]) & (M > 15))[0]
    i_1d  = np.where((pr["Mdot_crit_1d"]   > pr["Mdot_1d"])          & (M > 15))[0]
    mp = M[i_pol[0]] if len(i_pol) else 55.0
    m1 = M[i_1d[0]]  if len(i_1d)  else 220.0
    if len(i_pol):
        ax.scatter([mp], [pr["Mdot_disk"][i_pol[0]]], color="#0ca678", s=60, zorder=6)
        ax.annotate(f"3D Disk/Polar\nbreakout\n{mp:.0f} M☉", xy=(mp, pr["Mdot_disk"][i_pol[0]]),
                    xytext=(mp*0.18, 8e-5), arrowprops=dict(arrowstyle="->", color="#0ca678"), fontsize=7.4, color="#0ca678", fontweight="bold")
    if len(i_1d):
        ax.scatter([m1], [pr["Mdot_1d"][i_1d[0]]], color="#e03131", s=60, zorder=6)
        ax.annotate(f"Strict 1D\nbreakout\n{m1:.0f} M☉", xy=(m1, pr["Mdot_1d"][i_1d[0]]),
                    xytext=(m1*0.30, 1.8e-2), arrowprops=dict(arrowstyle="->", color="#e03131"), fontsize=7.4, color="#e03131", fontweight="bold")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.025, 300); ax.set_ylim(2e-5, 2e-1)
    ax.set_xlabel(r"Protostellar Mass $M_*$ [$M_\odot$]")
    ax.set_ylabel(r"Accretion / Critical Rate [$M_\odot\,\mathrm{yr}^{-1}$]")
    ax.legend(loc="lower left", frameon=True, framealpha=0.92, fontsize=6.6)
    ax.set_title("(c) Stage 2: ABN02 Accretion & Final Mass", fontweight="bold")

    # (d) Protostellar Radius R_*(M_*) & T_eff(M_*) — 3D Disk Track vs Strict 1D Track
    ax = axs[1, 0]; ax2 = ax.twinx()
    ax.plot(M, pr["R_disk_Rsun"], color="#0ca678", lw=2.3, label=r"$R_*$ (3D Disk → ZAMS)")
    ax.plot(M, pr["R_1d_Rsun"],   color="#e03131", lw=2.0, ls="--", label=r"$R_*$ (Strict 1D bloated)")
    ax.plot(M, pr["Rzams_Rsun"],  color="#495057", lw=1.3, ls=":", label=r"Pop III ZAMS $R_{\rm ZAMS}$")
    ax2.plot(M, pr["Teff_disk_K"] / 1e3, color="#1971c2", lw=2.0, ls="-.", label=r"$T_{\rm eff}$ (Disk, kK)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.025, 300); ax.set_ylim(0.5, 350)
    ax2.set_ylim(0, 118)
    ax.set_xlabel(r"Protostellar Mass $M_*$ [$M_\odot$]")
    ax.set_ylabel(r"Stellar Radius $R_*$ [$R_\odot$]", color="#e03131")
    ax2.set_ylabel(r"Effective Temperature $T_{\rm eff}$ [kK]", color="#1971c2")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", frameon=True, framealpha=0.92, fontsize=7.0)
    ax.set_title("(d) Stage 2: Swelling → KH Contraction → ZAMS", fontweight="bold")

    # (e) Stage 3: ABN02 Minihalo Photoevaporation (f_esc, f_bound, R_HII)
    ax = axs[1, 1]; ax2 = ax.twinx()
    ax.plot(pp["t_Myr"], pp["f_esc_H"], color="#0ca678", lw=2.2, label=rf"Disk/Polar ($M_{{\rm f}}\approx {mp:.0f}\,M_\odot$): $f_{{\rm esc}}$")
    ax.plot(pp["t_Myr"], pp["f_bound"], color="#0ca678", lw=1.6, ls="--", label=r"Disk/Polar: $f_{\rm bound}$")
    ax.plot(p1["t_Myr"], p1["f_esc_H"], color="#e03131", lw=2.2, label=rf"Strict 1D ($M_{{\rm f}}\approx {m1:.0f}\,M_\odot$): $f_{{\rm esc}}$")
    ax.plot(p1["t_Myr"], p1["f_bound"], color="#e03131", lw=1.6, ls="--", label=r"Strict 1D: $f_{\rm bound}$")
    ax2.plot(pp["t_Myr"], pp["nH_core"], color="#1971c2", lw=1.6, ls=":", label=r"Core $n_{\rm H}(t)$ [cm$^{-3}$] (Polar)")
    ax.set_xlim(0, max(pp["t_Myr"][-1], p1["t_Myr"][-1])); ax.set_ylim(-0.03, 1.08)
    ax2.set_yscale("log"); ax2.set_ylim(5e-2, 5e4)
    ax.set_xlabel("Post-Breakout Lifetime $t$ [Myr]")
    ax.set_ylabel(r"$f_{\rm esc}$ (solid) & Bound Gas $f_{\rm bound}$ (dashed)")
    ax2.set_ylabel(r"Core Density $n_{\rm H,core}$ [cm$^{-3}$]", color="#1971c2")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="center right", frameon=True, framealpha=0.92, fontsize=6.7)
    ax.set_title("(e) Stage 3: Full-Lifetime ABN02 Minihalo Photoevaporation", fontweight="bold")

    # (f) Stage 3: Time-Resolved JWST Atomic Line Signatures (Polar vs Strict 1D)
    ax = axs[1, 2]
    for ck, lab, col in [("L_Lya", r"Ly$\alpha$", "#e03131"), ("L_Hb", r"H$\beta$", "#1971c2"),
                         ("L_1640", r"He II $\lambda 1640$", "#5f3dc4"), ("L_10830", r"He I $\lambda 10830$", "#0ca678")]:
        m = pp[ck] > 1.0
        if np.any(m):
            ax.plot(pp["t_Myr"][m], pp[ck][m], color=col, lw=2.1, label=lab + " (Polar)")
        m1 = p1[ck] > 1.0
        if np.any(m1):
            ax.plot(p1["t_Myr"][m1], p1[ck][m1], color=col, lw=1.4, ls="--", alpha=0.8, label=lab + " (1D)")
    ax.set_yscale("log"); ax.set_xlim(0, max(pp["t_Myr"][-1], p1["t_Myr"][-1])); ax.set_ylim(5e1, 2e6)
    ax.set_xlabel("Post-Breakout Lifetime $t$ [Myr]")
    ax.set_ylabel(r"Emergent Nebular Luminosity [$L_\odot$]")
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, ncol=2, fontsize=6.5)
    ax.set_title("(f) Stage 3: ABN02 JWST Line Lightcurves", fontweight="bold")

    fig.suptitle("Titan.jl End-to-End ABN02 (Abel, Bryan & Norman 2002): Cosmological Collapse → Protostar → ZAMS → Minihalo Photoevaporation",
                 fontsize=11.6, fontweight="bold", y=0.985)
    plt.tight_layout(rect=[0, 0, 1, 0.962])
    plt.savefig(OUT, dpi=180)
    print(f"Wrote {OUT}")

if __name__ == "__main__":
    main()
