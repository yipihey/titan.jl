"""Movie of a recorded Pop III run (RUNREC) from its export (julia/runrec.jl).

    python scripts/popiii_movie.py <rundir> [--x r|m] [--fps 24] [--every 1] [--out movie.mp4]

Nine profile panels against log r (or log m) and a strip of the run's track, one frame per exported snapshot. Axis ranges are
fixed over the whole movie (from percentiles of all frames) so motion is real motion.
"""
import argparse, json, os, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, SEC, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1f9e6e"
plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": "#a8a79f", "axes.labelcolor": SEC, "xtick.color": SEC, "ytick.color": SEC,
                     "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": .6, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "axes.titlesize": 10,
                     "axes.titleweight": "bold", "axes.titlelocation": "left", "legend.frameon": False})


def load(rundir):
    ex = os.path.join(rundir, "export")
    meta = json.load(open(os.path.join(ex, "frames.json")))
    raw = np.fromfile(os.path.join(ex, "frames.f32"), dtype=np.float32)
    names = [f[0] for f in meta["fields"]]
    frames = []
    for fr in meta["frames"]:
        n, o = fr["n"], fr["off"] // 4
        a = raw[o:o + n * len(names)].reshape(len(names), n)
        sc = {k: (float("nan") if v is None else v) for k, v in zip(meta["scalars"], fr["s"])}
        frames.append(({k: a[i] for i, k in enumerate(names)}, sc))
    return frames


# (title, [(field, label, colour)], log-ish?, fixed range or None, reference line)
PANELS = [("Temperature (log K)", [("logT", "gas", BLUE), ("logTrad", "radiation", ORANGE)], None, None),
          ("Density (log g cm⁻³)", [("logrho", None, BLUE)], None, None),
          ("Velocity (km s⁻¹)", [("u", None, BLUE)], None, 0.0),
          ("Luminosity (L☉)", [("L", None, BLUE)], "symlog", None),
          ("Composition", [("xHp", "x(H⁺)", BLUE), ("xH2", "x(H₂)", ORANGE)], (-0.03, 1.03), None),
          ("Optical depth (log τ)", [("logtau", None, BLUE)], None, np.log10(2 / 3)),
          ("Eddington ratio Γ", [("Gamma", None, BLUE)], "log", 1.0),
          ("F / cE", [("FcE", None, BLUE)], "log", 1.0),
          ("Inflow rate Ṁ (M☉ yr⁻¹)", [("Mdot", None, BLUE)], "log", None)]


def steps(f, xk, y):
    """A zone's value across the zone, inner to outer interface (the first zone from its centre)."""
    xo = f.get(xk + "_o")
    if xo is None: return f[xk], y
    xi = np.concatenate([[f[xk][0]], xo[:-1]])
    return np.column_stack([xi, xo]).ravel(), np.repeat(y, 2)


def ranges(frames, xk):
    def pr(vals, lo=1, hi=99):
        v = np.concatenate([x[np.isfinite(x)] for x in vals])
        return (np.percentile(v, lo), np.percentile(v, hi)) if v.size else (0.0, 1.0)
    xr = pr([f[xk] for f, _ in frames] + [f[xk + "_o"] for f, _ in frames if xk + "_o" in f], 0, 100)
    yr = []
    for title, curves, mode, ref in PANELS:
        if isinstance(mode, tuple): yr.append(mode); continue
        vals = [f[c[0]] for f, _ in frames for c in curves]
        if mode == "log":
            v = np.concatenate([x[np.isfinite(x) & (x > 0)] for x in vals])
            if v.size == 0: yr.append((1e-3, 1.0)); continue
            lo, hi = np.percentile(np.log10(v), [1, 99.5])
            yr.append((10 ** (lo - 0.3), 10 ** (hi + 0.3)))
        else:
            lo, hi = pr(vals, 0.5, 99.5); pad = 0.05 * (hi - lo + 1e-12); yr.append((lo - pad, hi + pad))
    return xr, yr


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("rundir"); ap.add_argument("--x", default="r"); ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--every", type=int, default=1); ap.add_argument("--out", default=None); ap.add_argument("--dpi", type=int, default=110)
    a = ap.parse_args()
    frames = load(a.rundir)[::a.every]
    xk = "logr" if a.x == "r" else "logm"; xl = "log r (R☉)" if a.x == "r" else "log m (M☉)"
    xr, yr = ranges(frames, xk)
    T = np.array([s["t_yr"] for _, s in frames]); M = np.array([s["M_star"] for _, s in frames]); nc = np.array([s["logn_c"] for _, s in frames])
    Tc = np.array([s["T_c"] for _, s in frames])
    out = a.out or os.path.join(a.rundir, "export", f"movie_{a.x}.mp4")
    fig = plt.figure(figsize=(15, 9.4), constrained_layout=True)
    gs = fig.add_gridspec(4, 3, height_ratios=[1, 1, 1, 0.55])
    axs = [fig.add_subplot(gs[i // 3, i % 3]) for i in range(9)]
    axt = fig.add_subplot(gs[3, :])
    W, H = fig.canvas.get_width_height(); W, H = int(W * a.dpi / fig.dpi), int(H * a.dpi / fig.dpi)
    fig.set_dpi(a.dpi); W, H = fig.canvas.get_width_height(); W -= W % 2; H -= H % 2
    ff = subprocess.Popen(["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgba", "-s", f"{W}x{H}", "-r", str(a.fps),
                           "-i", "-", "-vf", "format=yuv420p", "-c:v", "libx264", "-crf", "20", out], stdin=subprocess.PIPE)
    i_t = np.arange(len(frames))
    for i, (f, s) in enumerate(frames):
        for ax, (title, curves, mode, ref), (ylo, yhi) in zip(axs, PANELS, yr):
            ax.cla(); ax.set_title(title)
            for k, lab, c in curves:
                y = f[k].astype(float)
                if mode == "log": y = np.where(y > 0, y, np.nan)
                ax.plot(*steps(f, xk, y), color=c, lw=1.6, label=lab)
            if mode == "log": ax.set_yscale("log")
            if mode == "symlog": ax.set_yscale("symlog", linthresh=1.0)
            ax.set_xlim(*xr); ax.set_ylim(ylo, yhi)
            if ref is not None: ax.axhline(ref, color=SEC, lw=1, ls=":")
            if any(c[1] for c in curves): ax.legend(loc="best", fontsize=8.5)
        for ax in axs[6:]: ax.set_xlabel(xl)
        axt.cla(); axt.set_title("Track: central density (left) and protostar mass (right)")
        axt.plot(i_t, nc, color=BLUE, lw=1.4); axt.set_ylabel("log n_c (cm⁻³)", color=BLUE); axt.axvline(i, color=INK, lw=1)
        axt.set_xlim(0, len(frames) - 1); axt.set_xlabel("snapshot")
        ax2 = axt.twinx(); ax2.plot(i_t, M, color=ORANGE, lw=1.4); ax2.set_ylabel("M* (M☉)", color=ORANGE); ax2.grid(False)
        Ms = s["M_star"]
        fig.suptitle(f"step {int(s['step'])}   t = {s['t_yr']:.6e} yr   Δt = {s['dt_s']:.2e} s   log n_c = {s['logn_c']:.2f}   T_c = {s['T_c']:.3g} K"
                     + (f"   M* = {Ms:.3f} M☉   R* = {s['R_star']:.3g} R☉" if Ms == Ms else "")
                     + (f"   r_ph = {s['r_ph']:.3g} R☉" if s['r_ph'] == s['r_ph'] else ""), x=0.01, ha="left", fontsize=11)
        fig.canvas.draw(); buf = np.asarray(fig.canvas.buffer_rgba())[:H, :W]
        ff.stdin.write(buf.tobytes()); ax2.remove()
    ff.stdin.close(); ff.wait(); print(out, len(frames), "frames")


if __name__ == "__main__":
    main()
