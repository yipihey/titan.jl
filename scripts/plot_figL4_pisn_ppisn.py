#!/usr/bin/env python3
"""
plot_figL4_pisn_ppisn.py

Publication-grade 4-panel (2x2) Figure L4:
  (a) Heger & Woosley (2002) Exact e± Pair-Instability Valley Γ_1(ρ, T) < 4/3 from FreeEnergyEOS
  (b) Kasen, Woosley & Heger (2011) 200 Msun PISN Shock Breakout UV Flash & Early Lightcurve
  (c) PISN Shock Breakout Radial Structure (T, rho, u) at Peak Flash
  (d) Woosley, Blinnikov & Heger (2007) PPISN Two-Shell Collision Superluminous Lightcurve & Forward/Reverse Shocks
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
OUT = "/home/tabel/.gemini/antigravity/brain/91f66bf7-1694-46da-a2b7-de91d0662618/figL4_pisn_shock_breakout_and_ppisn.png"

def load_csv(p):
    with open(p) as f:
        hdr = f.readline().strip().split(",")
    a = np.loadtxt(p, delimiter=",", skiprows=1)
    if a.ndim == 1:
        a = a[None, :]
    return {h: a[:, i] for i, h in enumerate(hdr)}

def main():
    pm   = load_csv(os.path.join(DIR, "lm07_pair_instability_map.csv"))
    plc  = load_csv(os.path.join(DIR, "lm07_pisn_breakout_lc.csv"))
    ppr  = load_csv(os.path.join(DIR, "lm07_pisn_breakout_prof.csv"))
    pp_l = load_csv(os.path.join(DIR, "lm07_ppisn_collision_lc.csv"))
    pp_p = load_csv(os.path.join(DIR, "lm07_ppisn_collision_prof.csv"))

    plt.rcParams.update({
        "font.size": 9.4, "axes.titlesize": 10.3, "axes.labelsize": 9.4,
        "legend.fontsize": 7.6, "xtick.labelsize": 8.4, "ytick.labelsize": 8.4,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
    })
    fig, axs = plt.subplots(2, 2, figsize=(13.8, 9.2), dpi=180)
    fig.patch.set_facecolor("white")

    # (a) Exact FreeEnergyEOS e± Pair Instability Valley Γ_1 < 4/3
    ax = axs[0, 0]
    lr = np.unique(pm["log10_rho"]); lT = np.unique(pm["log10_T"])
    G1 = pm["Gamma1_pairs"].reshape(len(lr), len(lT)).T
    cf = ax.contourf(lr, lT, G1, levels=np.linspace(1.12, 1.56, 23), cmap="RdYlBu", extend="both")
    cs = ax.contour(lr, lT, G1, levels=[4.0 / 3.0], colors=["black"], linewidths=2.2)
    ax.clabel(cs, fmt=r"$\Gamma_1 = 4/3$", fontsize=8.5, inline=True)
    # Overlay schematic Heger & Woosley (2002) core evolution adiabats (ρ_c ∝ T_c³) for 50, 100, 140 Msun He cores
    T_tr = np.linspace(8.45, 9.50, 40)
    for Mhe, shift, lab, col in [(50, 5.2, r"PPISN ($M_{\rm He}=50\,M_\odot$)", "#2b8a3e"),
                                 (100, 4.7, r"PISN ($M_{\rm He}=100\,M_\odot$)", "#c92a2a"),
                                 (140, 4.4, r"Direct BH ($M_{\rm He}=140\,M_\odot$)", "#495057")]:
        rho_tr = shift + 3.0 * (T_tr - 9.0)
        ax.plot(rho_tr, T_tr, color=col, lw=1.8, ls="--", label=lab)
    cb = plt.colorbar(cf, ax=ax, pad=0.02); cb.set_label(r"Adiabatic Index $\Gamma_1$ (`FreeEnergyEOS`)")
    ax.set_xlim(2.0, 6.8); ax.set_ylim(8.4, 9.6)
    ax.set_xlabel(r"$\log_{10} \rho$ [g cm$^{-3}$]")
    ax.set_ylabel(r"$\log_{10} T$ [K]")
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=7.2)
    ax.set_title(r"(a) Heger & Woosley (2002) $e^\pm$ Pair-Instability Valley ($\Gamma_1 < 4/3$)", fontweight="bold")

    # (b) Kasen+11 PISN Shock Breakout UV Flash & Early Lightcurve
    ax = axs[0, 1]; ax2 = ax.twinx()
    ax.plot(plc["t_day"], plc["L_ergs"], color="#e03131", lw=2.3, label=r"Bolometric $L(t)$ [erg s$^{-1}$]")
    ax2.plot(plc["t_day"], plc["Tph_K"] / 1e4, color="#1971c2", lw=1.9, ls="--", label=r"$T_{\rm ph}(t)$ [$10^4$ K]")
    ax.set_yscale("log"); ax.set_ylim(1e41, max(1e46, float(np.max(plc["L_ergs"])) * 3.0))
    ax2.set_yscale("log"); ax2.set_ylim(0.3, max(100.0, float(np.max(plc["Tph_K"])) / 1e4 * 2.0))
    ax.set_xlabel("Time since explosion $t$ [days]")
    ax.set_ylabel(r"Luminosity $L_{\rm bol}$ [erg s$^{-1}$]")
    ax2.set_ylabel(r"Photospheric $T_{\rm ph}$ [$10^4$ K]", color="#1971c2")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper right", frameon=True, framealpha=0.92)
    ax.set_title(r"(b) Kasen+11 $200\,M_\odot$ PISN Shock Breakout Flash ($5\times10^{52}$ erg)", fontweight="bold")

    # (c) PISN Radial Breakout Shock Structure
    ax = axs[1, 0]; ax2 = ax.twinx()
    ax.plot(ppr["r_Rsun"], ppr["u_kms"] / 1e3, color="#e03131", lw=2.2, label=r"Velocity $u(r)$ [$10^3$ km s$^{-1}$]")
    ax2.plot(ppr["r_Rsun"], ppr["T"], color="#1971c2", lw=2.0, ls="--", label=r"Temperature $T(r)$ [K]")
    ax2.plot(ppr["r_Rsun"], ppr["rho"] * 1e14, color="#495057", lw=1.4, ls=":", label=r"Density $\rho \times 10^{14}$")
    ax.set_xscale("log"); ax2.set_yscale("log")
    ax.set_xlabel(r"Radius $r$ [$R_\odot$]")
    ax.set_ylabel(r"Outflow Velocity $u$ [$10^3$ km s$^{-1}$]")
    ax2.set_ylabel(r"$T$ [K] & Scaled $\rho$", color="#1971c2")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="lower left", frameon=True, framealpha=0.92, fontsize=7.4)
    ax.set_title("(c) PISN Radiative Shock Breakout Front & Precursor", fontweight="bold")

    # (d) Woosley, Blinnikov & Heger (2007) PPISN Colliding Shells Lightcurve + Inset
    ax = axs[1, 1]
    ax.plot(pp_l["t_day"], pp_l["L_ergs"], color="#5f3dc4", lw=2.4, label=r"PPISN Collision $L_{\rm bol}(t)$ (10+10 $M_\odot$)")
    ax.axhline(1e44, color="#e03131", lw=1.1, ls="--", label=r"Superluminous SN threshold ($10^{44}$ erg s$^{-1}$)")
    ax.set_yscale("log"); ax.set_ylim(1e41, max(2e45, float(np.max(pp_l["L_ergs"])) * 3.0))
    ax.set_xlabel("Time $t$ [days]")
    ax.set_ylabel(r"Bolometric Luminosity $L_{\rm bol}$ [erg s$^{-1}$]")
    ax.legend(loc="upper left", frameon=True, framealpha=0.92, fontsize=7.4)
    ax.set_title("(d) Woosley+07 PPISN Two-Shell Collision (SLSN)", fontweight="bold")
    ins = ax.inset_axes([0.48, 0.12, 0.44, 0.42])
    ins.plot(pp_p["r_1e15cm"], pp_p["u_kms"] / 1e3, color="#e03131", lw=1.6, label=r"$u$ [$10^3$ km/s]")
    ins2 = ins.twinx()
    ins2.plot(pp_p["r_1e15cm"], np.log10(np.maximum(pp_p["rho"], 1e-18)), color="#1971c2", lw=1.3, ls="--")
    ins.set_xlabel(r"$r$ [$10^{15}$ cm]", fontsize=7); ins.set_ylabel(r"$u$ [$10^3$ km/s]", fontsize=7)
    ins2.set_ylabel(r"$\log_{10}\rho$", fontsize=6.5, color="#1971c2")
    ins.set_title("Shocked Shell @ Peak", fontsize=7.5, fontweight="bold")
    ins.tick_params(labelsize=6.5); ins2.tick_params(labelsize=6.5); ins.grid(True, alpha=0.25)

    fig.suptitle("Titan.jl Landmark Suite (Wave 4): Heger+02 Pair Instability, Kasen+11 PISN Breakout & Woosley+07 PPISN Colliding Shells",
                 fontsize=11.5, fontweight="bold", y=0.985)
    plt.tight_layout(rect=[0, 0, 1, 0.962])
    plt.savefig(OUT, dpi=180)
    print(f"Wrote {OUT}")

if __name__ == "__main__":
    main()
