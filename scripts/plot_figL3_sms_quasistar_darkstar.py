#!/usr/bin/env python3
"""
plot_figL3_sms_quasistar_darkstar.py

Publication-grade 4-panel (2x2) Figure L3:
  (a) Woods+17 / Haemmerlé+18 Supermassive Star GR Collapse Mass Curve M_coll(Mdot)
  (b) Live 1D RHD+GR Dynamical Collapse vs Stability of a 2.5e5 Msun SMS Core (GR_CORR=true)
  (c) Begelman (2008/2010) Quasi-Star & Spolyar/Freese (2008/2016) Dark Star HR Locus + Live 1D RHD QSTAR Inversion
  (d) JWST Rest-Frame UV-Optical-NIR SED Discriminant: 1e5 Msun Quasi-Star (LRD) vs Dark Star vs Accreting SMS vs Pop III ZAMS
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = "/zpool/nvme/data/tabel_scratch/titan_popiii/val/runs/landmarks"
OUT = "/home/tabel/.gemini/antigravity/brain/91f66bf7-1694-46da-a2b7-de91d0662618/figL3_sms_quasistar_darkstar.png"

def load_csv(p):
    with open(p) as f:
        hdr = f.readline().strip().split(",")
    lines = [l.strip().split(",") for l in open(p).readlines()[1:] if l.strip()]
    out = {}
    for i, h in enumerate(hdr):
        try:
            out[h] = np.array([float(r[i]) for r in lines])
        except ValueError:
            out[h] = np.array([r[i] for r in lines])
    return out

def main():
    mc  = load_csv(os.path.join(DIR, "lm05_sms_mcoll_curve.csv"))
    tr  = load_csv(os.path.join(DIR, "lm05_sms_gr_collapse_tracks.csv"))
    qs  = load_csv(os.path.join(DIR, "lm06_quasistar_seq.csv"))
    qp  = load_csv(os.path.join(DIR, "lm06_quasistar_rhd_prof.csv"))
    ds  = load_csv(os.path.join(DIR, "lm06_darkstar_seq.csv"))
    sed = load_csv(os.path.join(DIR, "lm06_jwst_seds.csv"))

    plt.rcParams.update({
        "font.size": 9.4, "axes.titlesize": 10.3, "axes.labelsize": 9.4,
        "legend.fontsize": 7.6, "xtick.labelsize": 8.4, "ytick.labelsize": 8.4,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
    })
    fig, axs = plt.subplots(2, 2, figsize=(13.8, 9.2), dpi=180)
    fig.patch.set_facecolor("white")

    # (a) Woods+17 M_coll(Mdot)
    ax = axs[0, 0]
    ax.plot(mc["Mdot_Msunyr"], mc["Mcoll_Msun"] / 1e5, color="#e03131", lw=2.4, marker="o", ms=5.5,
            label=r"Titan.jl (`FreeEnergyEOS` + `NucNet` + 1PN GR)")
    ax.plot(mc["Mdot_Msunyr"], mc["Mcoll_W17_Msun"] / 1e5, color="#1971c2", lw=2.0, ls="--", marker="s", ms=5.0,
            label=r"Woods et al. (2017) KEPLER reference")
    ax.plot(mc["Mdot_Msunyr"], mc["Mcore_Msun"] / 1e5, color="#7048e8", lw=1.6, ls=":", marker="^", ms=4.5,
            label=r"Convective H-burning core $M_{\rm core}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"Accretion Rate $\dot{M}$ [$M_\odot\,\mathrm{yr}^{-1}$]")
    ax.set_ylabel(r"GR Collapse Mass $M_{\rm coll}$ [$10^5\,M_\odot$]")
    ax.set_ylim(0, max(4.2, float(np.max(mc["Mcoll_Msun"])) / 1e5 * 1.15))
    ax.legend(loc="upper left", frameon=True, framealpha=0.92)
    ax.set_title(r"(a) Woods+17 / Haemmerlé+18 SMS GR Collapse Curve $M_{\rm coll}(\dot{M})$", fontweight="bold")

    # (b) Live 1D RHD+GR Collapse Tracks
    ax = axs[0, 1]; ax2 = ax.twinx()
    cases = [
        ("GR_supercrit_0p7Rc", r"GR ($0.7\,R_{\rm crit}$): COLLAPSE", "#e03131", 2.3, "-"),
        ("GR_subcrit_1p3Rc",   r"GR ($1.3\,R_{\rm crit}$): STABLE",   "#1971c2", 2.0, "--"),
        ("Newt_0p7Rc",         r"Newtonian ($0.7\,R_{\rm crit}$): STABLE", "#0ca678", 1.8, ":"),
    ]
    for nm, lab, col, lw, ls in cases:
        m = tr["case"] == nm
        if np.any(m):
            ax.plot(tr["t_tdyn"][m], tr["rho_rat"][m], color=col, lw=lw, ls=ls, label=lab)
            if nm == "GR_supercrit_0p7Rc":
                ax2.plot(tr["t_tdyn"][m], -tr["umin_c"][m], color="#d9480f", lw=1.7, ls="-.", label=r"Peak infall $|u_{\min}|/c$")
    ymax_r = max(1.75, float(np.max(tr["rho_rat"])) * 1.15)
    ax.set_ylim(0.85, ymax_r)
    ax2.set_ylim(0, max(0.0038, float(np.max(-tr["umin_c"])) * 1.25))
    ax.set_xlabel(r"Time $t / t_{\rm dyn}$")
    ax.set_ylabel(r"Central Density Compression $\rho_c / \rho_{c,0}$")
    ax2.set_ylabel(r"Infall Speed $|u_{\min}|/c$ (GR supercrit)", color="#d9480f")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", frameon=True, framealpha=0.92)
    ax.set_title(r"(b) Live 1D RHD+GR Collapse of $2.5\times10^5\,M_\odot$ SMS Core", fontweight="bold")

    # (c) Begelman Quasi-Star + Spolyar/Freese Dark Star HR & Live RHD QSTAR inset
    ax = axs[1, 0]
    m05 = np.isclose(qs["f_bh"], 0.05)
    ax.plot(qs["Teff_K"][m05], qs["L_Lsun"][m05], color="#e03131", lw=2.3, marker="o", ms=5.5,
            label=r"Begelman+08 Quasi-Stars ($10^4\to 10^6\,M_\odot$, LRD)")
    ax.plot(ds["Teff_K"], ds["L_Lsun"], color="#5f3dc4", lw=2.3, marker="s", ms=5.0,
            label=r"Spolyar/Freese Dark Stars ($10^2\to 10^6\,M_\odot$, 100 GeV)")
    # Add Pop III ZAMS reference line
    ax.scatter([95500.0], [1.3e6], color="#1971c2", s=60, marker="*", zorder=6, label=r"$100\,M_\odot$ Pop III ZAMS")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlim(1.3e5, 3000); ax.set_ylim(1e5, 1e11)
    ax.set_xlabel(r"Effective Temperature $T_{\rm eff}$ [K] (hot $\leftarrow \rightarrow$ cool)")
    ax.set_ylabel(r"Bolometric Luminosity $L$ [$L_\odot$]")
    ax.legend(loc="lower left", frameon=True, framealpha=0.92, fontsize=7.4)
    ax.set_title("(c) Quasi-Star (LRD) & Dark Star HR Locus (+ Live QSTAR RHD)", fontweight="bold")
    # Small inset showing live 1D RHD QSTAR super-Eddington Gamma_Edd profile (placed cleanly in empty upper-left)
    ins = ax.inset_axes([0.12, 0.46, 0.44, 0.36])
    ins.plot(qp["r_AU"], qp["Gamma_Edd"], color="#e03131", lw=1.6, label=r"$\Gamma_{\rm Edd}$")
    ins.axhline(1.0, color="#495057", lw=0.8, ls=":")
    ins.set_xlabel("r [AU]", fontsize=7); ins.set_ylabel(r"$\Gamma_{\rm Edd}$", fontsize=7)
    ins.set_title("Live 1D RHD QSTAR", fontsize=7.5, fontweight="bold")
    ins.tick_params(labelsize=6.5); ins.grid(True, alpha=0.25)

    # (d) JWST Rest-Frame SED Discriminant (1e5 Msun objects)
    ax = axs[1, 1]
    lam = sed["lam_um"]
    ax.plot(lam, sed["lLl_QS_1e5"],   color="#e03131", lw=2.3, label=r"$10^5\,M_\odot$ Quasi-Star ($T_{\rm eff}\approx 6.2$ kK, LRD)")
    ax.plot(lam, sed["lLl_SMS_1e5"],  color="#f08c00", lw=2.0, ls="--", label=r"$10^5\,M_\odot$ Accreting SMS ($T_{\rm eff}\approx 5.2$ kK)")
    ax.plot(lam, sed["lLl_DS_1e5"],   color="#5f3dc4", lw=2.2, ls="-.", label=r"$10^5\,M_\odot$ Supermassive Dark Star")
    ax.plot(lam, sed["lLl_ZAMS_1e5"], color="#1971c2", lw=2.0, ls=":",  label=r"$10^5\,M_\odot$ Pop III ZAMS Cluster ($95.5$ kK)")
    # Mark JWST NIRCam/NIRSpec rest-frame coverage at z = 10 (0.09 -> 0.45 μm)
    ax.axvspan(0.09, 0.45, color="#d0ebff", alpha=0.45, label=r"JWST NIRSpec rest-frame ($z=10$)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.08, 4.5); ax.set_ylim(1e6, 1e10)
    ax.set_xlabel(r"Rest-Frame Wavelength $\lambda_{\rm rest}$ [$\mu$m]")
    ax.set_ylabel(r"Spectral Luminosity $\lambda L_\lambda$ [$L_\odot$]")
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=7.2)
    ax.set_title(r"(d) JWST SED Discriminant ($10^5\,M_\odot$ Exotic Primordial Objects)", fontweight="bold")

    fig.suptitle("Titan.jl Landmark Suite (Wave 3): Woods+17 SMS GR Collapse, Begelman+08 Quasi-Stars & Spolyar/Freese Dark Stars",
                 fontsize=11.6, fontweight="bold", y=0.985)
    plt.tight_layout(rect=[0, 0, 1, 0.962])
    plt.savefig(OUT, dpi=180)
    print(f"Wrote {OUT}")

if __name__ == "__main__":
    main()
