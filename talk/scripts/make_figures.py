#!/usr/bin/env python3
"""All figures of the talk. From talk/: python3 scripts/make_figures.py
Reads report/app/data (the books' data), lib/popiii/*.csv and *.log, and talk/data (scripts/export_data.jl)."""
import json, os, re
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, Wedge, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..", "..")
BOOK = os.path.join(ROOT, "report", "app", "data"); POP = os.path.join(ROOT, "lib", "popiii")
DATA = os.path.join(HERE, "..", "data"); OUT = os.path.join(HERE, "..", "figures")
# the books' palette (report/_brand.yml)
BLUE, SKY, ORANGE, GREEN, INK = "#0072b2", "#56b4e9", "#d55e00", "#009e73", "#1f2328"
MUTED, GRID, PALE = "#5c6370", "#e3e5e8", "#f3f4f6"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 17, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "axes.linewidth": 1.0,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white",
                     "savefig.dpi": 160, "legend.frameon": False, "lines.linewidth": 2.6, "text.color": INK, "mathtext.default": "regular"})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight", pad_inches=0.08); plt.close(fig); print(name)


def book(name):
    return json.load(open(os.path.join(BOOK, name)))


def blank(w, h):
    fig, ax = plt.subplots(figsize=(w, h)); ax.set_axis_off(); ax.set_xlim(0, w); ax.set_ylim(0, h); ax.set_aspect("equal"); return fig, ax


def arrow(ax, a, b, c=BLUE, lw=2.4, ms=18):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=ms, color=c, lw=lw, shrinkA=0, shrinkB=0))


# ------------------------------------------------------------------ the common picture
def infall_picture(ax, cx, cy, title, rows, star, starc):
    for th in np.linspace(0, 2 * np.pi, 12, endpoint=False):
        arrow(ax, (cx + 3.0 * np.cos(th), cy + 3.0 * np.sin(th)), (cx + 1.85 * np.cos(th), cy + 1.85 * np.sin(th)), SKY)
    ax.add_patch(Circle((cx, cy), 1.62, fc="none", ec=ORANGE, lw=4, ls=(0, (1.2, 0.8))))
    ax.add_patch(Circle((cx, cy), 1.25, fc=starc, ec="none"))
    ax.text(cx, cy, star, ha="center", va="center", fontsize=17, color="white", weight="bold")
    ax.text(cx, cy + 3.45, title, ha="center", fontsize=23, weight="bold")
    for i, (k, v) in enumerate(rows):
        y = cy - 3.45 - 0.62 * i
        ax.text(cx - 0.12, y, k, ha="right", color=MUTED, fontsize=17); ax.text(cx + 0.12, y, v, ha="left", fontsize=17, weight="bold")


def fig_two_objects():
    fig, ax = blank(15, 10.2)
    infall_picture(ax, 3.7, 6.2, "A first star, forming", [("mass so far", "1.5 Suns"), ("size", "100 × the Sun"), ("fed at", "0.05 Suns per year"), ("gas arrives at", "65 km/s")], "young\nstar", ORANGE)
    infall_picture(ax, 11.3, 6.2, "A quasi-star", [("mass", "1,000,000 Suns"), ("size", "120,000 × the Sun"), ("fed at", "3 Suns per year"), ("gas arrives at", "1,670 km/s")], "black hole\ninside a star", INK)
    ax.text(7.5, 8.4, "falling gas", ha="center", color=BLUE, fontsize=17)
    ax.text(7.5, 6.2, "shock:\ngas stops,\nheats,\nand shines", ha="center", va="center", color=ORANGE, fontsize=17)
    save(fig, "two_objects")


# ------------------------------------------------------------------ the code in one picture
def fig_shells():
    fig, ax = blank(15, 7.2)
    cx, cy = 3.6, 3.5
    rs = np.geomspace(0.35, 3.2, 8)
    cols = plt.cm.Blues(np.linspace(0.75, 0.12, len(rs)))
    for r, c in list(zip(rs, cols))[::-1]:
        ax.add_patch(Wedge((cx, cy), r, -55, 55, fc=c, ec="white", lw=2))
    ax.add_patch(Wedge((cx, cy), rs[5], -55, 55, width=rs[5] - rs[4], fc=ORANGE, ec="white", lw=2))
    ax.text(cx - 0.2, cy, "centre", ha="right", va="center", color=MUTED)
    ax.text(cx + 1.7, 0.2, "the star is cut into\nnested spherical shells", ha="center", va="center", color=MUTED, fontsize=16)
    y0 = 6.25; ax.text(8.1, y0 + 0.5, "Each shell carries seven numbers", fontsize=22, weight="bold")
    for i, (a, b) in enumerate([("where it is", "radius"), ("how much gas is outside it", "mass"), ("how dense", "density"), ("how fast it moves", "velocity"),
                                ("how hot the gas is", "temperature"), ("how much light it holds", "radiation energy"), ("how much light flows through", "radiation flux")]):
        ax.text(8.1, y0 - 0.68 * i - 0.2, a, fontsize=18); ax.text(14.9, y0 - 0.68 * i - 0.2, b, fontsize=15, color=MUTED, ha="right")
        ax.plot([8.1, 14.9], [y0 - 0.68 * i - 0.42] * 2, color=GRID, lw=1)
    arrow(ax, (cx + rs[5] * 0.82, cy + rs[5] * 0.42), (7.9, 5.6), ORANGE, 2)
    save(fig, "shells")


def fig_mesh():
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.2))
    t = np.linspace(0, 1, 7)
    for ax, ttl, txt in zip(axs, ("Moves with the gas", "Fixed in space", "Adapts by itself"),
                            ("Each shell keeps its mass.\nBest inside a star.\nStretches badly in fast infall.",
                             "Shells stay put; gas flows through.\nBest for falling gas\nand a sharp shock.",
                             "Shells crowd into sharp features.\nNear the shock, zones 8,000× thinner\nthan on the fixed mesh.")):
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_title(ttl, fontsize=21, weight="bold", pad=10); ax.set_xticks([]); ax.set_yticks([])
        ax.spines["left"].set_visible(True); ax.spines["bottom"].set_visible(True)
        ax.set_xlabel(txt, fontsize=15, labelpad=8, color=INK)
    for i in range(16):                                              # Lagrangian: lines converge inward with time, spacing opening
        r0 = 0.05 + 0.058 * i
        axs[0].plot(t, r0 * (1 - 0.72 * t**1.5) ** (1 if i < 8 else 1.2 - 0.03 * i) , color=BLUE, lw=1.8)
    for r in np.geomspace(0.04, 1.0, 16):
        axs[1].plot([0, 1], [r, r], color=BLUE, lw=1.8)
    tt = np.linspace(0, 1, 60); sh = 0.72 - 0.32 * tt
    for i, u in enumerate(np.linspace(-1, 1, 19)):                    # adaptive: lines bunch around a moving shock
        axs[2].plot(tt, np.clip(sh + 0.5 * np.sign(u) * np.abs(u)**2.6, 0, 1), color=BLUE, lw=1.6)
    axs[2].plot(tt, sh, color=ORANGE, lw=3.5); axs[2].text(0.68, 0.40, "shock", color=ORANGE, fontsize=17, va="top", ha="left")
    axs[1].plot(tt, 0.32 + 0.3 * tt, color=ORANGE, lw=3.5); axs[1].text(0.60, 0.42, "shock", color=ORANGE, fontsize=17, va="top")
    for ax in axs:
        ax.text(0.0, 1.02, "radius", transform=ax.transAxes, color=MUTED, fontsize=13)
        ax.text(0.99, 0.015, "time →", transform=ax.transAxes, color=MUTED, fontsize=13, ha="right", va="bottom", bbox=dict(fc="white", ec="none", pad=1.5), zorder=6)
    fig.tight_layout(w_pad=2.5); save(fig, "mesh")


def fig_bench():
    b = book("bench_SUPC.json"); F, J = b["fortran"], b["julia"]
    rc = lambda d: 0.5 * (np.array(d["r"][:-1]) + np.array(d["r"][1:]))
    fig, ax = plt.subplots(figsize=(10, 6.2))
    ax.plot(rc(J) / 1e11, J["T"], color=BLUE, label="Gas temperature, Julia"); ax.plot(rc(J) / 1e11, J["Trad"], color=ORANGE, label="Light temperature, Julia")
    ax.plot(rc(F) / 1e11, F["T"], "o", ms=5, mfc="white", mec=INK, mew=1.1, markevery=2, label="Original 1995 Fortran"); ax.plot(rc(F) / 1e11, F["Trad"], "o", ms=5, mfc="white", mec=INK, mew=1.1, markevery=2)
    ax.set_xlabel("position in the test problem [10¹¹ cm]"); ax.set_ylabel("temperature [K]"); ax.grid(axis="y", color=GRID); ax.legend(loc="upper right", fontsize=15)
    ax.ticklabel_format(axis="y", style="plain"); fig.tight_layout(); save(fig, "bench")


# ------------------------------------------------------------------ microphysics
def fig_light():
    rng = np.random.default_rng(3); fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    star = dict(fc=ORANGE, ec="none")
    for ax, t, sub in zip(axs, ("Light escapes", "Some colours are trapped", "All light is trapped"),
                          ("thin gas: each shell cools directly", "thick H₂ lines: follow rays,\ncolour by colour", "deep inside: light diffuses outward,\none average field")):
        ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect("equal"); ax.set_axis_off(); ax.set_title(t, fontsize=21, weight="bold", pad=6)
        ax.text(0, -1.17, sub, ha="center", va="top", fontsize=16); ax.add_patch(Circle((0, 0), 1, fc=PALE, ec=MUTED, lw=1.4))
    for th in np.linspace(0.3, 6.2, 9):
        x, y = rng.uniform(-0.45, 0.45, 2); arrow(axs[0], (x, y), (np.cos(th) * 1.02, np.sin(th) * 1.02), SKY, 2.2, 14); axs[0].plot([x], [y], "o", color=BLUE, ms=6)
    for th in np.linspace(0.2, 6.2, 5):
        arrow(axs[1], (0, 0), (np.cos(th) * 1.02, np.sin(th) * 1.02), SKY, 2.2, 14)
    for k in range(3):
        p = rng.uniform(-0.2, 0.2, 2); q = [p]
        for _ in range(7): q.append(np.clip(q[-1] + rng.normal(0, 0.28, 2), -0.62, 0.62))
        q = np.array(q); axs[1].plot(q[:, 0], q[:, 1], color=ORANGE, lw=1.8); axs[1].plot(q[:, 0], q[:, 1], "o", color=ORANGE, ms=4)
    for k in range(5):
        th = rng.uniform(0, 6.28); q = np.cumsum(rng.normal(0, 0.05, (280, 2)), axis=0); q = q[np.hypot(*q.T).cumsum() >= 0]; h = np.hypot(*q.T); e = np.argmax(h > 0.98) or len(q) - 1; q = q[:e]
        axs[2].plot(q[:, 0], q[:, 1], color=ORANGE, lw=1.1)
    fig.tight_layout(w_pad=1); save(fig, "light")


def fig_ladder():
    fig, ax = plt.subplots(figsize=(15, 6.2)); ax.set_xscale("log"); ax.set_xlim(50, 1.2e9); ax.set_ylim(0, 6.6); ax.set_yticks([])
    for s in ("left",): ax.spines[s].set_visible(False)
    ax.set_xlabel("gas temperature [K]   (the exact boundaries move with density)")
    rows = [(5.45, 50, 2.5e3, BLUE, "H₂ molecules"), (5.45, 2.5e3, 1e4, SKY, "H atoms"), (5.45, 1e4, 3e8, ORANGE, "ionized: free electrons and protons"),
            (4.25, 50, 3e4, SKY, "He atoms"), (4.25, 3e4, 1.5e5, ORANGE, "He⁺"), (4.25, 1.5e5, 3e8, "#a84a00", "He⁺⁺ (bare nuclei)")]
    for y, a, b, c, t in rows:
        ax.add_patch(Rectangle((a, y), b - a, 0.85, fc=c, ec="white", lw=2.5))
        ax.text(np.sqrt(a * b), y + 0.42, t, ha="center", va="center", color="white" if c != SKY else INK, fontsize=15 if b / a > 6 else 12.5, weight="bold")
    ax.text(58, 6.42, "hydrogen", color=MUTED, fontsize=14, va="center"); ax.text(58, 5.2, "helium", color=MUTED, fontsize=14, va="center")
    ax.annotate("", (58, 1), (58, 1))
    for T, t, y in [(1e6, "deuterium burns", 2.95), (1e7, "hydrogen burns\n(proton chains)", 2.95), (1e8, "hydrogen burns\nusing carbon", 2.95)]:
        ax.plot([T, T], [0, y - 0.35], color=GREEN, lw=2); ax.plot([T], [y - 0.35], "o", color=GREEN, ms=10)
        ax.text(T, y - 0.15, t, ha="center", va="bottom", fontsize=14)
    for T, t in [(142, "coldest point of\nthe collapse, 142 K"), (2000, "molecules break up:\ncollapse restarts")]:
        ax.plot([T, T], [0, 1.35], color=INK, lw=1.6); ax.plot([T], [1.35], "o", color=INK, ms=8); ax.text(T, 1.55, t, ha="center", va="bottom", fontsize=14)
    ax.text(1.27e5, 1.55, "star's centre\ntoday, 127,000 K", ha="center", va="bottom", fontsize=14); ax.plot([1.27e5, 1.27e5], [0, 1.35], color=INK, lw=1.6); ax.plot([1.27e5], [1.35], "o", color=INK, ms=8)
    fig.tight_layout(); save(fig, "ladder")


def fig_opacity():
    d = pd.read_csv(os.path.join(DATA, "opacity.csv")); fig, ax = plt.subplots(figsize=(9.8, 6.4))
    for (ρ, c, lab) in [(1e-12, SKY, "10⁻¹² g/cm³, thin gas"), (1e-8, GREEN, "10⁻⁸"), (1e-4, BLUE, "10⁻⁴"), (1.0, ORANGE, "1 g/cm³, a star's interior")]:
        m = np.isclose(d.rho, ρ, rtol=1e-3, atol=0); ax.plot(10**d.logT[m], d.kR[m], color=c, label=lab)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e3, 1e8); ax.set_ylim(1e-6, 3e6); ax.grid(axis="y", color=GRID)
    ax.set_xlabel("temperature [K]"); ax.set_ylabel("how opaque the gas is [cm²/g]"); ax.legend(loc="upper right", fontsize=14)
    ax.text(3e7, 0.8, "light scatters off\nfree electrons", ha="center", va="bottom", fontsize=14, color=MUTED)
    fig.tight_layout(); save(fig, "opacity")


# ------------------------------------------------------------------ quasi-star
def fig_branches():
    b = pd.DataFrame(book("branches.json")); b = b[b.fBH == 0.1]; R0 = book("analytic.json")["R0"]
    f = pd.read_csv(os.path.join(ROOT, "steady", "output", "fold.csv")); f = f[(f.fBH == 0.1) & (f.nroots == 2)]
    b = b.assign(RR=b.R / R0); big = b[b.RR >= 0.5].sort_values("Mdot"); small = b[b.RR < 0.5].sort_values("Mdot")
    M1 = np.r_[big.Mdot.clip(lower=1e-3), f.Mdot]; R1 = np.r_[big.RR, 10**f.lR2 / R0]; M2 = np.r_[small.Mdot, f.Mdot]; R2 = np.r_[small.RR, 10**f.lR1 / R0]
    i = np.argsort(M2); fig, ax = plt.subplots(figsize=(9.6, 6.2))
    ax.plot(M1, R1, "-", color=BLUE, label="Large star"); ax.plot(M2[i], R2[i], "-", color=ORANGE, label="Compact star")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(axis="y", color=GRID); ax.set_xlabel("feeding rate [Suns per year]"); ax.set_ylabel("radius, relative to the unfed star"); ax.legend(loc="lower left", fontsize=15)
    fm = f.Mdot.max(); print(" last two-root rate", fm); ax.axvline(9.5, color=MUTED, lw=1.5, ls="--"); ax.text(8.5, 0.06, "beyond 9–10:\nno steady star", ha="right", va="top", color=MUTED, fontsize=15)
    ax.plot([3], [0.9876], "o", ms=12, color=INK, zorder=5); ax.annotate("the run on the next slide", (3, 0.9876), (0.05, 0.62), fontsize=15, arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
    fig.tight_layout(); save(fig, "branches")


def fig_settle():
    a = book("qstar_A220.json"); h = a["hist"]; R = a["R"]; t = np.array(h["t"]); fig, axs = plt.subplots(2, 1, figsize=(9.2, 6.4), sharex=True)
    mask = t > 0
    axs[0].plot(t[mask] / 1e8, np.array(h["rph"])[mask], color=BLUE); axs[0].axhline(0.9876, color=ORANGE, ls="--", lw=2); axs[0].set_ylabel("radius of the\nsurface  [R₀]")
    axs[0].text(5.95, 0.9876, "steady model", color=ORANGE, ha="right", va="top", fontsize=15); axs[0].set_ylim(0.975, 1.005)
    axs[1].plot(t[mask] / 1e8, np.array(h["Teff"])[mask], color=BLUE); axs[1].axhline(7359, color=ORANGE, ls="--", lw=2); axs[1].set_ylabel("surface\ntemperature [K]"); axs[1].set_ylim(7150, 7450)
    axs[1].set_xlabel("time [10⁸ s, about 3 years per unit]"); axs[1].set_xlim(0, 6)
    for ax in axs: ax.grid(axis="y", color=GRID)
    fig.tight_layout(); save(fig, "settle")


def fig_layer():
    a, e = book("qstar_A220.json"), book("qstar_R2.json"); R = a["R"]
    pick = lambda q, t: min(q["snaps"], key=lambda s: abs(s["t"] - t))
    A = pick(a, 1.0e8); E = pick(e, 1.0e8); print(" layer snapshots", A["t"], E["t"])
    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    def about_shock(s):
        ri = np.array(s["r"]) / R; u = np.array(s["u"]); j = np.argmax(np.diff(u[:len(ri)]) < 0) if False else np.argmin(np.diff(u)) ; n = len(s["T"])
        rc = 0.5 * (ri[:n] + ri[1:n + 1]); return (rc - ri[j + 1]) * 1e3, np.array(s["T"])
    x, y = about_shock(E); ax.plot(x, y, "-o", ms=7, color=ORANGE, label="Fixed mesh: the layer is lost")
    x, y = about_shock(A); ax.plot(x, y, "-", color=BLUE, label="Adaptive mesh: the layer is resolved")
    ax.set_xlim(-4, 2); ax.set_yscale("log"); ax.set_ylim(3e3, 1e9); ax.grid(axis="y", color=GRID)
    ax.set_xlabel("distance from the shock [thousandths of the star's radius]"); ax.set_ylabel("gas temperature [K]"); ax.legend(loc="upper left", fontsize=15)
    fig.tight_layout(); save(fig, "layer")


def fig_spectrum():
    d = book("spectra.json"); runs = {r["name"]: r for r in d["runs"]}; N, Lt = runs["W30_late_ext_nlte"], runs["W30_late_ext"]
    fig, ax = plt.subplots(figsize=(10.6, 6.2))
    ls = N["linespec"]; l = np.array(ls["lam"]); y = np.array(ls["nuLnu"])
    ll = np.array(Lt["lam"]); yl = np.array(Lt["nuLnu"])
    m = (ll > 900) & (ll < 3e4); ax.plot(ll[m], yl[m], color=MUTED, lw=2, label="Simple assumption: atoms in equilibrium")
    m = (l > 900) & (l < 3e4); ax.plot(l[m], y[m], color=BLUE, lw=1.6, label="Full calculation")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1200, 3e4); ax.set_ylim(1e41, 1e45); ax.grid(axis="y", color=GRID)
    ax.set_xlabel("wavelength [Å]"); ax.set_ylabel("power per log wavelength [erg/s]"); ax.legend(loc="lower center", fontsize=15)
    ym = y[(l > 6500) & (l < 6620)].max(); ax.annotate("Hα", (6563, ym), (6563 * 1.32, ym * 1.25), fontsize=17, weight="bold", arrowprops=dict(arrowstyle="-", color=INK, lw=1))
    ax.text(3400, 1.6e43 / 2, "step from\nhydrogen", ha="center", va="top", fontsize=14, color=MUTED)
    fig.tight_layout(); save(fig, "spectrum")


def fig_ionized():
    d = book("spectra.json"); runs = {r["name"]: r for r in d["runs"]}; R = d["R"]; fig, ax = plt.subplots(figsize=(7.6, 6.2))
    for n, c, lab in [("W30_late_ext_nlte", BLUE, "Full calculation"), ("W30_late_ext", MUTED, "Equilibrium")]:
        p = runs[n]["prof"]; r = np.array(p["r"]); x = np.array(p["xH"]); s = np.argsort(r)
        ax.plot(r[s], x[s], color=c, label=lab)
    ax.set_xscale("log"); ax.set_xlim(0.9, 120); ax.set_ylim(0, 1.04); ax.grid(axis="y", color=GRID); ax.set_xlabel("distance from the star [star's radius]"); ax.set_ylabel("fraction of hydrogen ionized"); ax.legend(loc="center left", bbox_to_anchor=(0.2, 0.3), fontsize=15)
    fig.tight_layout(); save(fig, "ionized")


# ------------------------------------------------------------------ first star
def fig_collapse():
    rows = []
    for l in open(os.path.join(POP, "m2c_abn02z19_lines_n200_final2.log")):
        p = l.split()
        if len(p) > 5 and p[0].isdigit():
            try: rows.append([float(x) for x in p[1:5]])
            except ValueError: pass
    a = np.array(rows); n, T, on = a[:, 0], a[:, 2], a[:, 3]
    h = pd.read_csv(os.path.join(DATA, "collapse_hot.csv")); r = h[h.run == "rad"].sort_values("log_n"); n3 = h[h.run == "n3"]; c = h[h.run == "c2"]
    fig, ax = plt.subplots(figsize=(14.5, 7))
    for x0, x1, t, col in [(3, 10.3, "H₂ lines escape freely", "#ffffff"), (10.3, 16.2, "H₂ lines trapped", PALE), (16.2, 21.5, "all light trapped", "#e6e8eb"), (21.5, 24, "a star", "#d9dce0")]:
        ax.axvspan(x0, x1, color=col, lw=0, zorder=0); ax.text((x0 + x1) / 2, 2.6e5, t, ha="center", va="top", color=MUTED, fontsize=15)
    m = n < 19; ax.plot(n[m], T[m], color=BLUE, lw=3, label="This work", zorder=3)
    m = r.log_n >= 19; ax.plot(r.log_n[m], r.Tc[m], color=BLUE, lw=3, zorder=3); ax.plot(n3.log_n, n3.Tc, color=BLUE, lw=3, zorder=3)
    ax.plot(c.log_n.iloc[:1], c.Tc.iloc[:1], "o", color=BLUE, ms=11, zorder=4)
    m = n >= 6.2; ax.plot(n[m], on[m], color=ORANGE, ls="--", lw=2.4, label="Omukai & Nishi 1998", zorder=2)
    ax.set_yscale("log"); ax.set_xlim(3, 24); ax.set_ylim(100, 3e5); ax.grid(axis="y", color="#d5d8dc", lw=0.8)
    ax.set_xlabel("central density  [log of atoms per cm³]"); ax.set_ylabel("central temperature [K]"); ax.legend(loc="center left", bbox_to_anchor=(0.01, 0.62), fontsize=16)
    ax.set_yticks([100, 1e3, 1e4, 1e5]); ax.set_yticklabels(["100", "1,000", "10,000", "100,000"])
    ax.annotate("20 million years\nnear 150 K", (4.2, 142), (6.2, 112), fontsize=15, va="bottom", arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
    ax.annotate("then runaway: the rest\ntakes about a million years", (11, 1100), (12.3, 300), fontsize=15, va="bottom", arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
    ax.annotate("molecules break up:\nsecond collapse", (20.7, 1e4), (16.7, 1.8e4), fontsize=15, va="bottom", arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
    ax.text(23.8, 1.5e4, "young\nstar", fontsize=15, ha="right", va="bottom")
    print(" collapse: first hot point", c.log_n.iloc[0], c.Tc.iloc[0]); fig.tight_layout(); save(fig, "collapse")


def fig_hr():
    fig, ax = plt.subplots(figsize=(9.4, 6.4))
    for f, c, m, lab, sch in [("s5_p100e.csv", ORANGE, 100, "100 Suns", (4.976, 6.11)), ("s5_p15e.csv", BLUE, 15, "15 Suns", (4.759, 4.324))]:
        d = pd.read_csv(os.path.join(POP, f)); d = d[d.t_yr > 1e3] if m == 100 else d[d.t_yr > 5e4]
        ax.plot(d.logTeff, d.logL, color=c, label=f"{lab}: this work"); ax.plot(*sch, "D", ms=12, mfc="white", mec=INK, mew=2, zorder=5)
    ax.plot([], [], "D", ms=10, mfc="white", mec=INK, mew=2, label="Start of hydrogen burning,\nSchaerer 2002")
    ax.set_xlim(5.05, 4.35); ax.set_xlabel("log surface temperature [K]   (hotter ←)"); ax.set_ylabel("log luminosity [Suns]"); ax.grid(axis="y", color=GRID); ax.legend(loc="center right", fontsize=14)
    fig.tight_layout(); save(fig, "hr")


def fig_anatomy():
    p = pd.read_csv(os.path.join(DATA, "profile_conv.csv")); fig, axs = plt.subplots(2, 1, figsize=(14, 6.9), sharex=True)
    j = (p.u < -10).values.argmax(); Rsh = p.r.iloc[j]
    axs[0].plot(p.r, p['T'], color=ORANGE); axs[0].set_yscale("log"); axs[0].set_ylabel("temperature [K]"); axs[0].set_yticks([1e3, 1e4, 1e5]); axs[0].set_yticklabels(["1,000", "10,000", "100,000"])
    axs[1].plot(p.r, -p.u, color=BLUE); axs[1].set_ylabel("falling speed [km/s]"); axs[1].set_xscale("log"); axs[1].set_xlim(0.1, p.r.max()); axs[1].set_xlabel("distance from the centre [radius of the Sun]")
    axs[1].set_ylim(0, 72)
    for ax in axs:
        ax.grid(axis="y", color=GRID); ax.axvspan(0.1, Rsh, color=PALE, lw=0, zorder=0); ax.axvline(Rsh, color=INK, lw=1.8)
    axs[0].text(0.15, 2200, "the star", fontsize=20, weight="bold"); axs[0].text(Rsh * 1.25, 3.6e4, f"shock, {Rsh:.0f} Suns' radii", fontsize=16)
    axs[1].text(3e5, 30, "the whole falling cloud stays on the same mesh,\nout to about 150 light-years", fontsize=16, ha="center", va="bottom")
    fig.tight_layout(h_pad=0.6); save(fig, "anatomy")


def fig_growth():
    d = pd.read_csv(os.path.join(DATA, "accretion.csv")); fig, axs = plt.subplots(1, 2, figsize=(14.5, 5.9))
    for run, c, lab in [("c1", MUTED, "No convection"), ("c2", BLUE, "Convection, restarted at day 1,000"), ("c3", ORANGE, "Convection from the start, finer mesh\n(stopped after 0.2 years)")]:
        x = d[d.run == run]
        axs[0].plot(x.t_yr, x.M, color=c, label=lab, lw=2.6 if run != "c1" else 2); axs[1].plot(x.M, x.R.rolling(5, center=True, min_periods=1).median(), color=c, label=lab, lw=2.6 if run != "c1" else 2)
    x = d[(d.run == "c1") & (d.M > 0.7)]; b = x.M.values[0]
    axs[0].set_xlabel("time since the start of this mesh [years]"); axs[0].set_ylabel("mass of the star [Suns]"); axs[0].grid(axis="y", color=GRID); axs[0].set_ylim(0.45, 1.6)
    axs[0].legend(loc="upper left", fontsize=15)
    axs[1].set_xlabel("mass of the star [Suns]"); axs[1].set_ylabel("radius of the shock [Suns' radii]"); axs[1].grid(axis="y", color=GRID); axs[1].set_ylim(25, 120); axs[1].set_xlim(0.48, 1.55)
    fig.tight_layout(w_pad=3); save(fig, "growth")


def fig_shell():
    a = pd.read_csv(os.path.join(DATA, "profile_noconv.csv")); b = pd.read_csv(os.path.join(DATA, "profile_conv.csv"))
    fig, axs = plt.subplots(1, 2, figsize=(14.5, 5.3), sharex=True)
    axs[0].plot(a.r, a['T'] / 1e3, color=MUTED, label="No convection"); axs[0].plot(b.r, b['T'] / 1e3, color=BLUE, label="Convection")
    axs[1].plot(a.r, a.rho, color=MUTED); axs[1].plot(b.r, b.rho, color=BLUE)
    axs[0].set_ylabel("temperature [thousand K]"); axs[1].set_ylabel("density [g/cm³]"); axs[1].set_yscale("log"); axs[0].legend(fontsize=15, loc="upper right")
    for ax in axs: ax.set_xlim(5, 60); ax.grid(axis="y", color=GRID); ax.set_xlabel("distance from the centre [radius of the Sun]")
    axs[1].set_ylim(5e-7, 5e-4); axs[0].set_ylim(0, 240)
    fig.tight_layout(w_pad=3); save(fig, "shell")


for f in (fig_two_objects, fig_shells, fig_mesh, fig_bench, fig_light, fig_ladder, fig_opacity, fig_branches, fig_settle, fig_layer, fig_spectrum, fig_ionized, fig_collapse, fig_hr, fig_anatomy, fig_growth, fig_shell):
    f()
