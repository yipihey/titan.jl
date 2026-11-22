#!/usr/bin/env python3
"""abnfig.py -- publication-quality PDF figure of a titan.jl Pop III protostar run, styled after
Fig. 2 of Abel, Bryan & Norman 2002 (Science 295, 93).

    python3 -I scripts/abnfig/abnfig.py RUN [--panel P] [--out PATH] [--times auto|t1,t2,...]
                                            [--xaxis abn|r|m] [--n N]

Page 1: radial profiles at N automatically chosen output times (log r column / enclosed-mass column).
Page 2: time histories of the protostar scalars and a summary table at the plotted times.
Page 3 (only if the export has cv_* fields, i.e. a CONVDIAG=1 run): sub-shock convection diagnostics of one frame.
Reads the explorer export P/RUN/frames.json + frames.b64.txt (+ perf.csv, P/runs.json); writes P/RUN/abnfig.pdf.
See README.md in this directory.
"""
import argparse
import base64
import datetime
import json
import math
import os
import subprocess
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib import ticker

DEFAULT_PANEL = "/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/panel"

# physical constants (cgs)
MH = 1.6735575e-24
KB = 1.380649e-16
RSUN_AU = 6.957e10 / 1.495978707e13          # 1 R_sun in AU
LOG_RSUN_AU = math.log10(RSUN_AU)
GAMMA_AD = 5.0 / 3.0
ERG2LSUN = 1.0 / 3.828e33      # L_nuc, L_lines scalars are in erg/s
PROTO_LOGNC = 20.0       # a frame "has a protostar" when M_star is finite > 0 AND log n_c >= this (see README)

# ----------------------------------------------------------------------------------------------- style
PAGE_W, PAGE_H = 8.5, 11.0       # US Letter portrait [in]
FS = 8.5                         # base font size [pt]
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": FS, "axes.labelsize": FS + 0.5, "axes.linewidth": 0.7,
    "xtick.labelsize": FS - 0.5, "ytick.labelsize": FS - 0.5,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 4.0, "ytick.major.size": 4.0, "xtick.minor.size": 2.0, "ytick.minor.size": 2.0,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.45, "ytick.minor.width": 0.45,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "legend.fontsize": FS, "legend.frameon": False, "pdf.fonttype": 42, "axes.unicode_minus": True,
    "lines.solid_capstyle": "round", "lines.dash_capstyle": "butt",
})
INK, MUTED, FAINT = "#1a1a1a", "#555555", "#9a9a9a"

# dash patterns for the time sequence (all distinct in greyscale), in on/off points
DASHES = [(None), (6.0, 1.6), (1.1, 1.3), (7.0, 1.6, 1.4, 1.6), (3.0, 1.4), (9.0, 1.8, 1.4, 1.8, 1.4, 1.8),
          (2.0, 1.0, 2.0, 3.2), (12.0, 2.5)]
PRE_DASHES = [(1.0, 1.6), (4.0, 1.4, 1.0, 1.4, 1.0, 1.4)]
PRE_COLS = ["#7a7a7a", "#b0b0b0"]


def time_colours(n):
    """Sequential, colour-blind-safe, greyscale-monotone (dark -> light) colours for n times, early -> late."""
    cm = plt.get_cmap("plasma")
    if n <= 1:
        return [cm(0.0)]
    return [cm(0.0 + 0.80 * i / (n - 1)) for i in range(n)]


# ----------------------------------------------------------------------------------------------- data
def fnum(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return math.nan
    return x if math.isfinite(x) else math.nan


class Run:
    def __init__(self, panel, key):
        self.key, self.dir = key, os.path.join(panel, key)
        with open(os.path.join(self.dir, "frames.json")) as fh:
            d = json.load(fh)
        self.fields = [f[0] for f in d["fields"]]
        self.fdesc = {f[0]: f[1] for f in d["fields"]}
        self.scal = list(d["scalars"])
        self.meta = d.get("frames", [])
        with open(os.path.join(self.dir, "frames.b64.txt")) as fh:
            self.buf = base64.b64decode(fh.read())
        S = np.full((len(self.meta), len(self.scal)), np.nan)
        for i, f in enumerate(self.meta):
            for j, v in enumerate(f.get("s", [])[:len(self.scal)]):
                S[i, j] = fnum(v)
        self.S = S
        self.nf = len(self.meta)
        self._cache = {}
        self.info = {}
        try:
            with open(os.path.join(panel, "runs.json")) as fh:
                for r in json.load(fh).get("runs", []):
                    if r.get("key") == key:
                        self.info = r
        except (OSError, ValueError):
            pass
        self.perf = read_perf(os.path.join(self.dir, "perf.csv"))

    def s(self, name):
        return self.S[:, self.scal.index(name)] if name in self.scal else np.full(self.nf, np.nan)

    def frame(self, i):
        if i not in self._cache:
            f = self.meta[i]
            n, nf = int(f["n"]), len(self.fields)
            a = np.frombuffer(self.buf, dtype="<f4", count=n * nf, offset=int(f["off"])).astype(np.float64)
            a = a.reshape(nf, n)
            self._cache[i] = {k: a[j] for j, k in enumerate(self.fields)}
        return self._cache[i]


def read_perf(path):
    try:
        with open(path) as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
    except OSError:
        return {}
    if len(lines) < 2:
        return {}
    hdr = lines[0].split(",")
    cols = {h: [] for h in hdr}
    for ln in lines[1:]:
        p = ln.split(",")
        if len(p) != len(hdr) or p[0] == hdr[0]:
            continue
        for h, v in zip(hdr, p):
            cols[h].append(fnum(v))
    out = {h: np.array(v) for h, v in cols.items()}
    if "t_yr" in out and len(out["t_yr"]):
        o = np.argsort(out["t_yr"], kind="stable")
        out = {h: v[o] for h, v in out.items()}
    return out


# ----------------------------------------------------------------------------------------------- derived profiles
def derived(fr):
    """Physical profiles (linear values, float64) from a frame's fields."""
    g = lambda k: fr.get(k, None)
    n = len(fr["logr"])
    nanv = np.full(n, np.nan)
    logXH = g("logXH")
    XH = 10.0 ** logXH if logXH is not None else nanv.copy()
    XH = np.where(np.isfinite(XH) & (XH > 0), XH, 0.76)
    logXHe = g("logXHe4")
    XHe = 10.0 ** logXHe if logXHe is not None else nanv.copy()
    XHe = np.where(np.isfinite(XHe) & (XHe > 0), XHe, np.clip(1.0 - XH, 0.0, 1.0))
    rho = 10.0 ** fr["logrho"]
    T = 10.0 ** fr["logT"]
    xH2 = np.clip(np.nan_to_num(g("xH2") if g("xH2") is not None else nanv, nan=0.0), 0, 0.5)
    xHp = np.clip(np.nan_to_num(g("xHp") if g("xHp") is not None else nanv, nan=0.0), 0, 1.0)
    yHe = XHe / 4.0 / XH                                  # He nuclei per H nucleus
    ntot_per_H = 1.0 - xH2 + 2.0 * xHp + yHe - xHp        # H + H2 + H+ + e- + He  (He taken neutral)
    mu = (1.0 / XH) / ntot_per_H
    cs = np.sqrt(GAMMA_AD * KB * T / (mu * MH)) / 1e5     # km/s
    d = {
        "r": 10.0 ** fr["logr"], "m": 10.0 ** fr["logm"], "nH": XH * rho / MH, "T": T,
        "Trad": 10.0 ** fr["logTrad"] if "logTrad" in fr else nanv,
        "u": fr.get("u", nanv), "cs": cs, "L": fr.get("L", nanv),
        "xH2": g("xH2") if g("xH2") is not None else nanv, "xHp": g("xHp") if g("xHp") is not None else nanv,
        "xHm": 10.0 ** fr["logxHm"] if "logxHm" in fr else nanv,
        "tau": 10.0 ** fr["logtau"] if "logtau" in fr else nanv,
        "Mdot": fr.get("Mdot", nanv), "Gamma": fr.get("Gamma", nanv),
    }
    return d


# ----------------------------------------------------------------------------------------------- time selection
def protostar_mask(run):
    M = run.s("M_star")
    nc = run.s("logn_c")
    return np.isfinite(M) & (M > 0) & (~np.isfinite(nc) | (nc >= PROTO_LOGNC))


def prof_vec(run, i, grid):
    fr = run.frame(i)
    x = fr["logr"]
    o = np.argsort(x)
    out = []
    for k, w in (("logrho", 1.0), ("logT", 1.0)):
        y = fr[k][o]
        ok = np.isfinite(x[o]) & np.isfinite(y)
        out.append(w * np.interp(grid, x[o][ok], y[ok], left=np.nan, right=np.nan) if ok.sum() > 1
                   else np.full(len(grid), np.nan))
    y = fr["u"][o]
    ok = np.isfinite(x[o]) & np.isfinite(y)
    u = np.interp(grid, x[o][ok], y[ok], left=np.nan, right=np.nan) if ok.sum() > 1 else np.full(len(grid), np.nan)
    out.append(np.sign(u) * np.log10(1 + np.abs(u) / 3.0))       # velocity: soft-log, ~dex units
    return np.concatenate(out)


def farthest_point(cands, D, start, n):
    chosen = list(dict.fromkeys(start))
    while len(chosen) < n and len(chosen) < len(cands):
        best, bd = None, -1.0
        for c in cands:
            if c in chosen:
                continue
            dm = min(D[c][s] for s in chosen)
            if dm > bd:
                best, bd = c, dm
        if best is None or bd <= 1e-9:
            break
        chosen.append(best)
    return chosen


def select_times(run, nsel, tsel=None):
    """Return sorted frame indices to plot (pre-protostar first), and t0 (or None)."""
    nf = run.nf
    t = run.s("t_yr")
    ps = protostar_mask(run)
    t0 = float(t[np.argmax(ps)]) if ps.any() else None
    if nf == 0:
        return [], t0
    if tsel:                         # explicit times since t0 [yr] (or since the first frame without protostar)
        ref = t0 if t0 is not None else t[0]
        idx = []
        for v in tsel:
            idx.append(int(np.nanargmin(np.abs(t - (ref + v)))))
        return sorted(dict.fromkeys(idx), key=lambda i: t[i]), t0
    if nf <= nsel:
        return list(range(nf)), t0
    allidx = list(range(nf))
    # subsample candidates for speed
    step = max(1, nf // 400)
    pre = [i for i in allidx if not ps[i]]
    pro = [i for i in allidx if ps[i]]
    lo = np.nanmin([np.nanmin(run.frame(i)["logr"]) for i in (allidx[0], allidx[-1])])
    hi = np.nanmax([np.nanmax(run.frame(i)["logr"]) for i in (allidx[0], allidx[-1])])
    grid = np.linspace(lo, hi, 60)

    def dist_matrix(cands, extra):
        V = {c: prof_vec(run, c, grid) for c in cands}
        D = {}
        for a in cands:
            D[a] = {}
            for b in cands:
                dv = V[a] - V[b]
                ok = np.isfinite(dv)
                dp = math.sqrt(np.mean(dv[ok] ** 2)) if ok.any() else 0.0
                D[a][b] = math.sqrt(dp ** 2 + extra(a, b) ** 2)
        return D

    chosen = []
    npre = 0
    if pre:
        npre = min(len(pre), 2 if pro else nsel)
        nc = run.s("logn_c")
        cp = pre[::max(1, len(pre) // 200)]
        if pre[-1] not in cp:
            cp.append(pre[-1])
        D = dist_matrix(cp, lambda a, b: 0.5 * abs(fnum(nc[a]) - fnum(nc[b])) if np.isfinite(nc[a] - nc[b]) else 0.0)
        chosen += farthest_point(cp, D, [pre[0]] + ([pre[-1]] if npre > 1 and not pro else []), npre)
    if pro:
        npro = nsel - len(chosen)
        cp = pro[::step]
        if pro[-1] not in cp:
            cp.append(pro[-1])
        tau = t - t0
        M = run.s("M_star")
        tpos = tau[pro][tau[pro] > 0]
        ts = max(np.nanmin(tpos) if len(tpos) else 1.0, 1e-4 * max(tau[pro[-1]], 1e-30))
        lt = lambda i: math.log10(max(tau[i], 0) + ts)
        lm = lambda i: math.log10(M[i])

        def extra(a, b):
            return math.sqrt((1.5 * (lm(a) - lm(b))) ** 2 + (0.6 * (lt(a) - lt(b))) ** 2)
        D = dist_matrix(cp, extra)
        chosen += farthest_point(cp, D, [pro[0], pro[-1]], npro)
    return sorted(dict.fromkeys(chosen), key=lambda i: (t[i], i)), t0


# ----------------------------------------------------------------------------------------------- formatting
def sig(v, n=3):
    if not np.isfinite(v):
        return "—"
    if v == 0:
        return "0"
    e = math.floor(math.log10(abs(v)))
    d = max(0, n - 1 - e)
    s = f"{v:.{d}f}"
    return s


def fmt_dt(yr):
    if not np.isfinite(yr):
        return "—"
    a = abs(yr)
    if a == 0:
        return "0 yr"
    if a < 1e3:
        return f"{sig(yr, 3)} yr"
    if a < 1e6:
        return f"{sig(yr / 1e3, 3)} kyr"
    return f"{sig(yr / 1e6, 3)} Myr"


def fmt_sci(v, n=3):
    """10^x mathtext for table entries."""
    if not np.isfinite(v):
        return "—"
    if v == 0:
        return "0"
    e = math.floor(math.log10(abs(v)))
    if -2 <= e <= 3:
        return sig(v, n)
    m = v / 10.0 ** e
    return rf"${m:.{n - 1}f}\!\times\!10^{{{e}}}$"


def curve_label(run, i, t0):
    t = run.s("t_yr")[i]
    M = run.s("M_star")[i]
    nc = run.s("logn_c")[i]
    if t0 is not None and protostar_mask(run)[i]:
        return rf"$t-t_0$ = {fmt_dt(t - t0)},  $M_*$ = {sig(M, 3)} M$_\odot$"
    ncs = rf"$\log_{{10}} n_{{\rm c}}$ = {nc:.2f}" if np.isfinite(nc) else ""
    if t0 is None:                  # no protostar anywhere: time since the first frame + central density
        return rf"$t-t_{{\rm first}}$ = {fmt_dt(t - run.s('t_yr')[0])},  " + ncs
    return rf"pre-protostar: {ncs}  ($t-t_0$ = {fmt_dt(t - t0)})"


def git_commit():
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        out = subprocess.run(["git", "-C", here, "log", "-1", "--format=%h %cd", "--date=short"],
                             capture_output=True, text=True, timeout=5)
        br = subprocess.run(["git", "-C", here, "rev-parse", "--abbrev-ref", "HEAD"],
                            capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            return f"{out.stdout.strip()} ({br.stdout.strip()})"
    except (OSError, subprocess.SubprocessError):
        pass
    return "n/a"


# ----------------------------------------------------------------------------------------------- axes helpers
def decade_ticks(ax, which, lo, hi, maxmaj=6):
    """Major ticks at 10^(k j) with k chosen so that <= maxmaj labels; minor ticks at every decade or at 2..9."""
    axis = ax.xaxis if which == "x" else ax.yaxis
    lo10, hi10 = math.log10(lo), math.log10(hi)
    span = hi10 - lo10
    if which == "y" and span < 2.6:          # few decades: label 1, 2, 5 x 10^k
        maj = [m * 10.0 ** e for e in range(int(math.floor(lo10)) - 1, int(math.ceil(hi10)) + 1) for m in ((1, 2, 5) if span < 1.6 else (1, 3))]
        maj = [v for v in maj if lo <= v <= hi]
        axis.set_major_locator(ticker.FixedLocator(maj))
        axis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=100))

        def g(v, _pos):
            fr = (math.log10(v) - lo10) / span
            if fr < 0.06 or fr > 0.94:
                return ""
            e = math.floor(math.log10(v) + 1e-9)
            m = round(v / 10.0 ** e)
            if -2 <= e <= 4:
                return f"{v:.{max(0, -e)}f}"
            return rf"$10^{{{e}}}$" if m == 1 else rf"${m}\times10^{{{e}}}$"
        axis.set_major_formatter(ticker.FuncFormatter(g))
        axis.set_minor_formatter(ticker.NullFormatter())
        return
    k = max(1, math.ceil(span / maxmaj))
    for kk in (1, 2, 3, 4, 5, 10, 20):
        if kk >= k:
            k = kk
            break
    maj = [10.0 ** e for e in range(int(math.floor(lo10 / k) * k), int(math.ceil(hi10)) + 1, k)]
    axis.set_major_locator(ticker.FixedLocator(maj))
    if span <= 7:
        axis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
    else:
        axis.set_minor_locator(ticker.FixedLocator([10.0 ** e for e in range(int(math.floor(lo10)) - 1,
                                                                             int(math.ceil(hi10)) + 2)]))

    def f(v, _pos):
        if v <= 0:
            return ""
        e = math.log10(v)
        if abs(e - round(e)) > 1e-6:
            return ""
        # suppress labels too close to the axis ends (stacked panels share edges)
        fr = (e - lo10) / span if span > 0 else 0.5
        if (which == "y" and (fr < 0.06 or fr > 0.94)) or (which == "x" and (fr < 0.015 or fr > 0.985)):
            return ""
        e = int(round(e))
        return rf"$10^{{{e}}}$"
    axis.set_major_formatter(ticker.FuncFormatter(f))
    axis.set_minor_formatter(ticker.NullFormatter())


def lin_ticks_trim(ax):
    lo, hi = ax.get_ylim()

    def f(v, _pos):
        fr = (v - lo) / (hi - lo)
        if fr < 0.06 or fr > 0.94:
            return ""
        return f"{v:g}".replace("-", "−")
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(f))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())


def narrow_log_ticks(ax):
    """Log axis spanning < ~1.2 dex: 3-5 'nice' linear tick values, labelled plainly or as m x 10^k."""
    lo, hi = ax.get_ylim()
    vals = [v for v in ticker.MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi) if lo <= v <= hi and v > 0]
    if len(vals) < 2:
        vals = [v for v in ticker.MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
                if lo <= v <= hi and v > 0]
    ax.yaxis.set_major_locator(ticker.FixedLocator(vals))
    e = math.floor(math.log10(max(vals)) + 1e-9) if vals else 0
    step = (vals[1] - vals[0]) if len(vals) > 1 else 10.0 ** e
    nd = max(0, -math.floor(math.log10(step / 10.0 ** e) + 1e-9))

    def f(v, _p):
        if v <= 0:
            return ""
        fr = (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
        if fr < 0.06 or fr > 0.94:              # stacked panels share edges: keep labels off them
            return ""
        if -2 <= e <= 4:
            dd = max(0, -math.floor(math.log10(step) + 1e-9))
            return f"{v:.{dd}f}"
        m = v / 10.0 ** e
        return rf"${m:.{nd}f}\times10^{{{e}}}$"
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(f))
    ax.yaxis.set_minor_locator(ticker.NullLocator())


def wrap_math(s, width):
    """textwrap that never breaks inside $...$ mathtext."""
    out, inm = [], False
    for ch in s:
        if ch == "$":
            inm = not inm
        out.append("\x00" if (ch == " " and inm) else ch)
    return "\n".join(ln.replace("\x00", " ") for ln in textwrap.wrap("".join(out), width))


CORNERS = {"ul": (0.0, 0.16, 0.74, 1.0), "ur": (0.84, 1.0, 0.74, 1.0),
           "ll": (0.0, 0.16, 0.0, 0.26), "lr": (0.84, 1.0, 0.0, 0.26)}


def line_points(ax, extra=()):
    """All plotted line vertices + segment midpoints, densified, in axes coordinates."""
    pts = []
    tr = ax.transData + ax.transAxes.inverted()
    for ln in list(ax.get_lines()) + list(extra):
        xy = np.asarray(ln.get_xydata(), float)
        if xy.ndim != 2 or len(xy) == 0:
            continue
        ok = np.all(np.isfinite(xy), axis=1)
        if ax.get_xscale() == "log":
            ok &= xy[:, 0] > 0
        if ax.get_yscale() == "log":
            ok &= xy[:, 1] > 0
        a = tr.transform(np.where(ok[:, None], xy, 1.0))
        a[~ok] = np.nan
        if len(a) > 1:                     # densify each segment (8 sub-points), keep breaks at NaN
            t = np.linspace(0, 1, 8)[:, None, None]
            seg = a[:-1][None] * (1 - t) + a[1:][None] * t
            a = np.concatenate([a, seg.reshape(-1, 2)])
        pts.append(a[np.all(np.isfinite(a), axis=1)])
    return np.concatenate(pts) if pts else np.zeros((0, 2))


def emptiest(pts, prefs, taken=(), boxes=CORNERS):
    best, bn = None, None
    for c in prefs:
        if c in taken:
            continue
        x0, x1, y0, y1 = boxes[c]
        n = int(np.sum((pts[:, 0] >= x0) & (pts[:, 0] <= x1) & (pts[:, 1] >= y0) & (pts[:, 1] <= y1)))
        if bn is None or n < bn:
            best, bn = c, n
    return best


def plot_xy(ax, x, y, **kw):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if kw.pop("pos", False):
        ok &= (y > 0)
    yy = np.where(ok, y, np.nan)
    o = np.argsort(x)
    return ax.plot(x[o], yy[o], **kw)


def interp_at(x, y, x0, logy=True):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0)
    if logy:
        ok &= y > 0
    if ok.sum() < 2 or not np.isfinite(x0) or x0 <= 0:
        return np.nan
    o = np.argsort(x[ok])
    lx = np.log10(x[ok][o])
    l0 = math.log10(x0)
    if l0 < lx[0] or l0 > lx[-1]:
        return np.nan
    yy = np.log10(y[ok][o]) if logy else y[ok][o]
    v = np.interp(l0, lx, yy)
    return 10.0 ** v if logy else v


def letter(ax, s, loc="ul"):
    x, ha = (0.022, "left") if loc[1] == "l" else (0.978, "right")
    y, va = (0.93, "top") if loc[0] == "u" else (0.07, "bottom")
    ax.text(x, y, s, transform=ax.transAxes, ha=ha, va=va, fontsize=FS + 2.5, fontweight="bold", color=INK, zorder=10)


# ----------------------------------------------------------------------------------------------- page 1
PANELS = {
    # key: (ylabel, quantity, log?, span cap [dex] or None, fixed lower limit)
    "nH": (r"$n_{\rm H}$  [cm$^{-3}$]", True),
    "m": (r"$M_{\rm enc}$  [M$_\odot$]", True),
    "r": (r"$r$  [R$_\odot$]", True),
    "L": (r"$L = 4\pi r^2 F$  [L$_\odot$]", True),
    "tau": (r"$\tau_{\rm R}$ (from outside)", True),
    "Mdot": (r"$\dot M = -4\pi r^2\rho\, v_r$  [M$_\odot$ yr$^{-1}$]", True),
    "T": (r"$T_{\rm gas}$, $T_{\rm rad}$  [K]", True),
    "u": (r"$v_r$, $-c_s$  [km s$^{-1}$]", False),
    "xH2": (r"$x({\rm H_2})$", True),
    "xHp": (r"$x({\rm H^+})$", True),
    "xHm": (r"$x({\rm H^-})$", True),
}
SPAN = {"nH": None, "m": 14, "r": None, "L": 9, "tau": 12, "Mdot": 7, "T": None, "xH2": 9, "xHp": 12, "xHm": 12}


def page_profiles(pdf, run, sel, t0, xmode, styles, header):
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    header(fig, "Radial structure at the selected output times")
    left_q = ["nH", "m", "L", "tau", "Mdot"]
    right_q = ["T", "u", "xH2", "xHp", "xHm"]
    colx = {"abn": ("r", "m"), "r": ("r", "r"), "m": ("m", "m")}[xmode]
    nrow = len(left_q)
    # geometry [in]
    L0, R0, gap = 0.80, 0.80, 0.10
    top, bot = 1.55, 2.05
    pw = (PAGE_W - L0 - R0 - gap) / 2
    ph = (PAGE_H - top - bot) / nrow
    profs = {i: derived(run.frame(i)) for i in sel}
    Rs = run.s("R_star")
    rph = run.s("r_ph")
    Ms = run.s("M_star")

    # x ranges per x-quantity from the plotted frames
    def xr(q):
        v = np.concatenate([profs[i][q] for i in sel]) if sel else np.array([1.0, 10.0])
        v = v[np.isfinite(v) & (v > 0)]
        if not len(v):
            return 1.0, 10.0
        lo, hi = np.min(v), np.max(v)
        if q == "m":
            lo = max(lo, hi * 1e-14)
            # do not waste space on the near-empty innermost zones of pre-protostar frames
            lo_pro = [np.nanmin(profs[i]["m"]) for i in sel if protostar_mask(run)[i]]
            if lo_pro:
                lo = max(lo, 10.0 ** (math.log10(min(lo_pro)) - 1.0))
        return 10.0 ** (math.floor(math.log10(lo) * 2) / 2), 10.0 ** (math.ceil(math.log10(hi) * 2) / 2)
    XR = {"r": xr("r"), "m": xr("m")}

    letters = iter("ABCDEFGHIJKLMNOP")
    tlet = {}
    for col, qs in enumerate((left_q, right_q)):
        xq = colx[col]
        axes = []
        for row, q in enumerate(qs):
            if q == xq:                         # all-r: M_enc(r) stays; all-m: show r(M) instead of M(M)
                q = "r" if xq == "m" else "m"
            x0 = L0 + col * (pw + gap)
            y0 = PAGE_H - top - (row + 1) * ph
            ax = fig.add_axes([x0 / PAGE_W, y0 / PAGE_H, pw / PAGE_W, ph / PAGE_H])
            axes.append(ax)
            ax.set_xscale("log")
            label, islog = PANELS[q]
            if islog:
                ax.set_yscale("log")
            ally = []
            for i in sel:
                p = profs[i]
                st = styles[i]
                x = p[xq]
                y = p[q]
                plot_xy(ax, x, y, pos=islog, color=st["c"], lw=st["lw"], dashes=st["d"] or (None, None),
                        ls="-", zorder=3, solid_capstyle="round")
                ally.append(y)
                if q == "T":
                    plot_xy(ax, x, p["Trad"], pos=True, color=st["c"], lw=0.55, alpha=0.75,
                            dashes=st["d"] or (None, None), zorder=2)
                if q == "u":
                    plot_xy(ax, x, -p["cs"], color=st["c"], lw=0.55, alpha=0.75, dashes=st["d"] or (None, None),
                            zorder=2)
                # protostar radius (filled circle) and photosphere (open square) on the curve
                if st["pro"]:
                    for xv, mk, fc in ((Rs[i] if xq == "r" else Ms[i], "o", st["c"]),
                                       (rph[i] if xq == "r" else interp_at(p["r"], p["m"], rph[i]), "s", "white")):
                        yv = interp_at(x, y, xv, logy=islog)
                        if np.isfinite(yv):
                            ax.plot([xv], [yv], marker=mk, ms=3.6 if mk == "o" else 3.2, mfc=fc, mec=st["c"],
                                    mew=0.8, ls="none", zorder=5)
            # y limits
            yy = np.concatenate([np.asarray(a, float) for a in ally]) if ally else np.array([1.0])
            xx = np.concatenate([np.asarray(profs[i][xq], float) for i in sel]) if sel else np.array([1.0])
            inx = (xx >= XR[xq][0]) & (xx <= XR[xq][1])
            yy = yy[np.isfinite(yy) & inx]
            if islog:
                yy = yy[yy > 0]
                if q == "T":
                    tr = np.concatenate([profs[i]["Trad"] for i in sel])
                    tr = tr[np.isfinite(tr) & (tr > 0) & inx]
                    yy = np.concatenate([yy, tr])
                if len(yy) == 0:
                    yy = np.array([1e-1, 1e1])
                hi = np.max(yy)
                lo = np.min(yy)
                if q in ("xH2",):
                    hi = max(hi, 0.5)
                if q in ("xHp",):
                    hi = max(hi, 1.0)
                if SPAN.get(q):
                    lo = max(lo, hi / 10.0 ** SPAN[q])
                l10lo = math.floor(math.log10(lo) - 0.15)
                l10hi = math.log10(hi) + 0.25
                l10hi = math.ceil(l10hi) if (math.ceil(l10hi) - l10hi) < 0.6 else l10hi
                ax.set_ylim(10.0 ** l10lo, 10.0 ** l10hi)
                decade_ticks(ax, "y", 10.0 ** l10lo, 10.0 ** l10hi, maxmaj=5)
            else:
                if len(yy) == 0:
                    yy = np.array([-1.0, 1.0])
                lo, hi = np.min(yy), np.max(yy)
                hi = max(hi, 0.0)
                pad = 0.08 * (hi - lo if hi > lo else 1.0)
                ax.set_ylim(lo - pad, hi + pad)
                ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
                lin_ticks_trim(ax)
                ax.axhline(0.0, color=FAINT, lw=0.5, zorder=1)
            if q == "tau" and ax.get_ylim()[0] < 2.0 / 3.0 < ax.get_ylim()[1]:
                ax.axhline(2.0 / 3.0, color=FAINT, lw=0.6, ls=(0, (3, 2)), zorder=1)
                ax.text(0.985, 2.0 / 3.0, r"$\tau = 2/3$", transform=ax.get_yaxis_transform(), ha="right",
                        va="bottom", fontsize=FS - 1, color=MUTED)
            ax.set_xlim(*XR[xq])
            decade_ticks(ax, "x", *XR[xq], maxmaj=7)
            ax.set_ylabel(label)
            if col == 1:
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
                ax.tick_params(axis="y", which="both", left=True, right=True, labelleft=False, labelright=True)
            if row < nrow - 1:
                ax.tick_params(axis="x", labelbottom=False)
            pts = line_points(ax)
            c = emptiest(pts, ["ul", "ur", "ll", "lr"])
            tlet[q] = next(letters)
            letter(ax, tlet[q], c)
        xlabel = (r"radius $r$  [R$_\odot$]" if xq == "r" else r"enclosed mass $M_{\rm enc}$  [M$_\odot$]")
        axes[-1].set_xlabel(xlabel)
        if xq == "r":                      # AU along the top of the first panel of an r column
            sec = axes[0].secondary_xaxis("top", functions=(lambda v: v * RSUN_AU, lambda v: v / RSUN_AU))
            lo, hi = XR["r"][0] * RSUN_AU, XR["r"][1] * RSUN_AU
            decade_ticks(sec, "x", lo, hi, maxmaj=7)
            sec.tick_params(axis="x", which="both", direction="in")
            sec.set_xlabel(r"$r$  [AU]", labelpad=3)
            axes[0].tick_params(axis="x", which="both", top=False)
        elif xmode == "abn":
            pass
    # legend + key below the panels
    handles = [Line2D([], [], color=styles[i]["c"], lw=styles[i]["lw"] + 0.2,
                      dashes=styles[i]["d"] or (None, None)) for i in sel]
    labels = [curve_label(run, i, t0) for i in sel]
    ncol = 2
    # column-major fill (time increases down the first column)
    leg_y = (bot - 0.50) / PAGE_H
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, leg_y), ncol=ncol,
               handlelength=4.2, columnspacing=2.5, handletextpad=0.7, labelspacing=0.45, fontsize=FS)
    key = [Line2D([], [], marker="o", ms=3.6, mfc=INK, mec=INK, ls="none"),
           Line2D([], [], marker="s", ms=3.2, mfc="white", mec=INK, mew=0.8, ls="none")]
    key.append(Line2D([], [], color=INK, lw=0.55))
    keylab = [r"protostar radius $R_*$ (accretion shock)", r"photosphere $r_{\rm ph}$ ($\tau_{\rm R}=2/3$)",
              r"thin lines: $T_{\rm rad}$ (%s), $-c_s$ (%s)" % (tlet.get("T", "?"), tlet.get("u", "?"))]
    fig.legend(key, keylab, loc="upper center", bbox_to_anchor=(0.5, leg_y + 0.016), ncol=3, fontsize=FS - 0.5,
               handletextpad=0.3, columnspacing=2.0)
    cap = ("Mass-weighted zone values (cell centres) of the 1D radiation-hydrodynamics run. $n_{\\rm H}=X_{\\rm H}\\rho/m_{\\rm H}$; "
           "$x(\\cdot)$ are abundances per H nucleus; $\\tau_{\\rm R}$ is the grey Rosseland optical depth integrated inward from the "
           "outer boundary; $L\\leq0$ and $\\dot M\\leq 0$ (outflow) are not shown on the log axes. "
           "$c_s=(\\gamma kT/\\mu m_{\\rm H})^{1/2}$ with $\\gamma=5/3$ and $\\mu$ from $x({\\rm H_2})$, $x({\\rm H^+})$ and neutral He "
           "(an estimate; the run's own EOS is not exported).")
    fig.text(0.5, 0.30 / PAGE_H, wrap_math(cap, 185), ha="center", va="bottom", fontsize=FS - 1.5,
             color=MUTED, linespacing=1.3)
    pdf.savefig(fig)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------- page 2
def page_history(pdf, run, sel, t0, styles, header):
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    header(fig, "Time evolution and summary")
    t = run.s("t_yr")
    ps = protostar_mask(run)
    pf = run.perf
    have_pro = t0 is not None
    ref = t0 if have_pro else (t[0] if run.nf else 0.0)
    tau = t - ref
    use = ps if have_pro else np.ones(run.nf, bool)

    def series(fname, pname=None, scale=1.0):
        """(x, y) from perf.csv when it has the column (denser), else from frames."""
        if pname and pf and pname in pf and "t_yr" in pf and len(pf["t_yr"]) > len(t):
            x = pf["t_yr"] - ref
            y = pf[pname] * scale
            if have_pro and "mstar" in pf:
                ok = np.isfinite(pf["mstar"]) & (pf["mstar"] > 0) & (pf.get("log_nc", np.full(len(x), 99)) >= PROTO_LOGNC)
                ok &= x >= 0
                y = np.where(ok, y, np.nan)
            elif pname in ("mstar", "rstar"):
                y = np.full(len(x), np.nan)
            return x, y
        msk = ps if fname in ("M_star", "R_star") else use      # the code's M_star/R_star mean nothing pre-protostar
        return tau, np.where(msk, run.s(fname) * scale, np.nan)

    xs = tau[use]
    xs = xs[np.isfinite(xs) & (xs > 0)]
    if len(xs) == 0:
        xs = np.array([1.0, 10.0])
    xhi = np.max(xs)
    xlo = max(np.min(xs), xhi * 1e-6)
    xlo = 10.0 ** math.floor(math.log10(xlo))
    xhi = 10.0 ** (math.log10(xhi) + 0.08)
    CAT = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]       # Okabe-Ito
    # panels: (letter, ylabel, list of (fname, pname, label, colour, dashes), log)
    P = [
        (r"$M_*$  [M$_\odot$]", [("M_star", "mstar", None, INK, None)], True),
        (r"$T_{\rm eff}$  [K]", [("T_eff_lum", None, r"$(L_{\rm ph}/4\pi r_{\rm ph}^2\sigma)^{1/4}$", CAT[0], None),
                                 ("T_eff_sl", None, r"$T_{\rm eff}$ (run)", CAT[1], (5, 1.6))], True),
        (r"$R_*$, $r_{\rm ph}$  [R$_\odot$]", [("R_star", "rstar", r"$R_*$", INK, None),
                                               ("r_ph", None, r"$r_{\rm ph}$", CAT[0], (5, 1.6))], True),
        (r"$\dot M_{10^4}$  [M$_\odot$ yr$^{-1}$]", [("Mdot_1e4", None, None, INK, None)], True),
        (r"$L$  [L$_\odot$]", [("L_ph", None, r"$L_{\rm ph}$", INK, None),
                               ("L_nuc", None, r"$L_{\rm nuc}$", CAT[1], (5, 1.6)),
                               ("L_lines", None, r"$L_{\rm lines}$", CAT[2], (1.2, 1.3))], True),
        (r"$-u_{\min}$  [km s$^{-1}$]", [("u_min", "umin_kms", None, INK, None)], "neg"),
        (r"$n_{\rm c}$  [cm$^{-3}$]", [("logn_c", "log_nc", None, INK, None)], "exp"),
        (r"$T_{\rm c}$  [K]", [("T_c", "Tc", None, INK, None)], True),
    ]
    nrow, ncol = 4, 2
    L0, R0, gap = 0.80, 0.80, 0.10
    top = 1.40
    ph = 1.30
    pw = (PAGE_W - L0 - R0 - gap) / 2
    letters = iter("ABCDEFGH")
    axs = []
    for k, (ylab, ser, mode) in enumerate(P):
        row, col = k // ncol, k % ncol
        x0 = L0 + col * (pw + gap)
        y0 = PAGE_H - top - (row + 1) * ph
        ax = fig.add_axes([x0 / PAGE_W, y0 / PAGE_H, pw / PAGE_W, ph / PAGE_H])
        axs.append(ax)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ally = []
        for fname, pname, lab, c, d in ser:
            x, y = series(fname, pname, ERG2LSUN if fname in ("L_nuc", "L_lines") else 1.0)
            if mode == "neg":
                y = -y
            elif mode == "exp":
                y = 10.0 ** y
            ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
            if not ok.any():
                continue
            o = np.argsort(x[ok])
            ax.plot(x[ok][o], y[ok][o], color=c, lw=0.9, dashes=d or (None, None), label=lab, zorder=3)
            ally.append(y[ok])
        # mark the plotted profile times
        for i in sel:
            if not (np.isfinite(tau[i]) and tau[i] > 0 and use[i]):
                continue
            fname, pname = ser[0][0], ser[0][1]
            v = run.s(fname)[i]
            v = -v if mode == "neg" else (10.0 ** v if mode == "exp" else v)
            if np.isfinite(v) and v > 0:
                ax.plot([tau[i]], [v], marker="o", ms=4.2, mfc=styles[i]["c"], mec="white", mew=0.6, ls="none",
                        zorder=6)
        if not ally:
            ax.text(0.5, 0.5, "no protostar" if fname in ("M_star", "R_star", "T_eff_lum") and not have_pro
                    else ("single frame: no history" if run.nf < 2 else "no positive data"), transform=ax.transAxes, ha="center", va="center", color=MUTED,
                    fontsize=FS, style="italic")
        yy = np.concatenate(ally) if ally else np.array([1.0, 10.0])
        lo, hi = np.min(yy), np.max(yy)
        l10lo, l10hi = math.log10(lo), math.log10(hi)
        span = l10hi - l10lo
        if span < 0.3:                      # nearly constant: show a 0.3-dex window around it
            mid = 0.5 * (l10lo + l10hi)
            l10lo, l10hi = mid - 0.15, mid + 0.15
        else:
            l10lo, l10hi = l10lo - 0.10 * span, l10hi + 0.14 * span
        ax.set_ylim(10.0 ** l10lo, 10.0 ** l10hi)
        if l10hi - l10lo >= 1.2:
            decade_ticks(ax, "y", 10.0 ** l10lo, 10.0 ** l10hi, maxmaj=5)
        else:
            narrow_log_ticks(ax)
        ax.set_xlim(xlo, xhi)
        decade_ticks(ax, "x", xlo, xhi, maxmaj=8)
        ax.set_ylabel(ylab)
        if col == 1:
            ax.yaxis.set_label_position("right")
            ax.tick_params(axis="y", which="both", labelleft=False, labelright=True)
        if row < nrow - 1:
            ax.tick_params(axis="x", labelbottom=False)
        else:
            ax.set_xlabel(r"$t-t_0$  [yr]" if have_pro else r"$t-t_{\rm first}$  [yr]")
        if col == 1:
            ax.yaxis.tick_right()
            ax.tick_params(axis="y", which="both", left=True, right=True, labelleft=False, labelright=True)
        pts = line_points(ax)
        if any(not str(ln.get_label()).startswith("_") for ln in ax.get_lines()):      # legend first (matplotlib 'best' avoids the curves), then the letter avoids both
            leg = ax.legend(loc="best", fontsize=FS - 1, handlelength=2.6, borderaxespad=0.4, labelspacing=0.25)
            bb = leg.get_window_extent().transformed(ax.transAxes.inverted())
            gx, gy = np.meshgrid(np.linspace(bb.x0, bb.x1, 12), np.linspace(bb.y0, bb.y1, 8))
            pts = np.concatenate([pts, np.repeat(np.c_[gx.ravel(), gy.ravel()], 50, axis=0)])
        letter(ax, next(letters), emptiest(pts, ["ul", "ll", "ur", "lr"]))
    if have_pro:
        note = (rf"$t_0$ = {t0:,.2f} yr (run clock) = first frame with a protostar "
                rf"($M_*>0$, $n_{{\rm c}}\geq10^{{{PROTO_LOGNC:.0f}}}$ cm$^{{-3}}$); $M_*(t_0)$ = {sig(run.s('M_star')[np.argmax(ps)], 3)} M$_\odot$. "
                "Coloured dots: the output times of page 1. Curves from perf.csv where denser, otherwise from the frames.")
    else:
        note = "No protostar formed in the exported frames; time is measured from the first frame."
    fig.text(L0 / PAGE_W, (PAGE_H - top - nrow * ph - 0.50) / PAGE_H, wrap_math(note, 190),
             fontsize=FS - 1.5, color=MUTED, va="top")

    # ---- summary table
    ty = PAGE_H - top - nrow * ph - 0.95          # top of the table [in]
    cols = [("", 0.55), (r"$t-t_0$" if have_pro else r"$t-t_{\rm first}$", 0.70), (r"$M_*$", 0.55), (r"$R_*$", 0.55), (r"$r_{\rm ph}$", 0.55),
            (r"$T_{\rm eff}$", 0.60), (r"$L_{\rm ph}$", 0.80), (r"$\dot M_{10^4}$", 0.85), (r"$\log n_{\rm c}$", 0.55),
            (r"$T_{\rm c}$", 0.80), (r"$-u_{\min}$", 0.55), (r"$r(v_{\min})$", 0.62)]
    units = ["", "", r"M$_\odot$", r"R$_\odot$", r"R$_\odot$", "K", r"L$_\odot$", r"M$_\odot$ yr$^{-1}$", r"cm$^{-3}$", "K",
             r"km s$^{-1}$", r"R$_\odot$"]
    tw = sum(w for _, w in cols)
    x0 = (PAGE_W - tw) / 2
    rowh = 0.205
    xs = [x0]
    for _, w in cols:
        xs.append(xs[-1] + w)
    def T(x, y, s, **kw):
        fig.text(x / PAGE_W, y / PAGE_H, s, **kw)
    T(x0, ty + 0.12, "Summary at the plotted times", fontsize=FS + 0.5, fontweight="bold", va="bottom")
    for j, ((h, w), u) in enumerate(zip(cols, units)):
        T(xs[j] + w / 2, ty - 0.02, h, ha="center", va="top", fontsize=FS)
        T(xs[j] + w / 2, ty - 0.20, u, ha="center", va="top", fontsize=FS - 1.5, color=MUTED)
    def hline(y, lw):
        fig.add_artist(Line2D([x0 / PAGE_W, xs[-1] / PAGE_W], [y / PAGE_H, y / PAGE_H], color=INK, lw=lw))
    hline(ty + 0.04, 0.9)
    hline(ty - 0.40, 0.5)
    y = ty - 0.45
    for i in sel:
        st = styles[i]
        fr = run.frame(i)
        u = fr["u"]
        jmin = int(np.nanargmin(u)) if np.isfinite(u).any() else None
        rvmin = 10.0 ** fr["logr"][jmin] if jmin is not None else np.nan
        tt = t[i] - t0 if have_pro else t[i] - ref
        teff = run.s("T_eff_lum")[i]
        teff = teff if np.isfinite(teff) and teff > 0 else run.s("T_eff_sl")[i]
        vals = [fmt_dt(tt) if (have_pro and ps[i]) else ("pre" if have_pro else fmt_dt(tt)),
                sig(run.s("M_star")[i], 3) if ps[i] else "—",
                fmt_sci(run.s("R_star")[i]) if ps[i] else "—",
                fmt_sci(run.s("r_ph")[i]),
                fmt_sci(teff if teff > 0 else np.nan),
                fmt_sci(run.s("L_ph")[i]), fmt_sci(run.s("Mdot_1e4")[i]),
                f"{run.s('logn_c')[i]:.2f}" if np.isfinite(run.s("logn_c")[i]) else "—",
                fmt_sci(run.s("T_c")[i]), sig(-run.s("u_min")[i], 3), fmt_sci(rvmin)]
        yc = y - rowh / 2
        fig.add_artist(Line2D([(xs[0] + 0.05) / PAGE_W, (xs[1] - 0.08) / PAGE_W], [yc / PAGE_H, yc / PAGE_H],
                              color=st["c"], lw=st["lw"] + 0.2, dashes=st["d"] or (None, None)))
        fig.add_artist(Line2D([(xs[0] + 0.27) / PAGE_W], [yc / PAGE_H], marker="o", ms=3.8, mfc=st["c"], mec="white",
                              mew=0.5, ls="none"))
        for j, v in enumerate(vals):
            T(xs[j + 1] + cols[j + 1][1] / 2, yc, v, ha="center", va="center", fontsize=FS - 0.5)
        y -= rowh
    hline(y - 0.03, 0.9)
    tn = (r"$R_*$: the code's accretion-shock radius (first interface from the centre with $v_r<-2$ km s$^{-1}$); $M_*=M_{\rm enc}(R_*)$. "
          r"$T_{\rm eff}=(L_{\rm ph}/4\pi r_{\rm ph}^2\sigma)^{1/4}$ (run's $T_{\rm eff}$ where $L_{\rm ph}\leq 0$). "
          r"$\dot M_{10^4}$: mass inflow rate through $r=10^4$ R$_\odot$. $-u_{\min}$, $r(v_{\min})$: peak infall speed and where it occurs. "
          r"$n_{\rm c}=\rho_{\rm c}/m_{\rm H}$, $T_{\rm c}$: centre (innermost zone).")
    fig.text(x0 / PAGE_W, (y - 0.10) / PAGE_H, wrap_math(tn, 190), fontsize=FS - 1.5, color=MUTED,
             va="top")
    pdf.savefig(fig)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------- main
# ----------------------------------------------------------------------------------------------- page 3 (CONVDIAG runs only)
def page_convection(pdf, run, sel, t0, header, frame=None):
    """Sub-shock convection diagnostics of one frame (needs the cv_* fields of a CONVDIAG=1 run; see docs/POPIII_CONVECTION.md M1).

    Left column vs r/R_star, right column vs M_enc.  Quantities of buoyancy (nabla, nabla_ad, N^2, Ledoux) are drawn only in
    zones with Mach < 0.5 (the supersonic infall has no meaningful stratification), the gates and fluxes everywhere.
    The accretion shock (R_star: first interface with u < -2 km/s) is the vertical dotted line.  The frame is the one with the
    strongest unstable layer (cv_N2td2) among the selected times, or --cvframe."""
    cand = [i for i in (sel if frame is None else [frame]) if np.isfinite(run.frame(i).get("cv_mach", np.array([np.nan]))).any()]
    if not cand:
        return False
    if frame is None:
        td = run.s("cv_N2td2") if "cv_N2td2" in run.scal else np.full(run.nf, np.nan)
        cand.sort(key=lambda i: (np.nan_to_num(td[i], nan=-1.0), i))
    i = cand[-1]
    fr = run.frame(i)
    Rst = run.s("R_star")[i]
    Mst = run.s("M_star")[i]
    r = 10.0 ** fr["logr"]
    m = 10.0 ** fr["logm"]
    xr = r / Rst if np.isfinite(Rst) and Rst > 0 else r
    xlab_r = r"$r/R_\star$" if np.isfinite(Rst) and Rst > 0 else r"$r$  [R$_\odot$]"
    sub = fr["cv_mach"] < 0.5
    g = lambda k: fr.get(k, np.full(len(r), np.nan))
    selx = np.isfinite(xr) & (xr > (0.1 if np.isfinite(Rst) else 0) ) & (xr < (3.0 if np.isfinite(Rst) else np.inf))
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    t_yr = run.s("t_yr")[i]
    header(fig, f"Sub-shock convection diagnostics (step {int(run.s('step')[i])}, t - t0 = {fmt_dt(t_yr - t0)})")
    axs = [[fig.add_axes([0.11 + 0.45 * c, 0.80 - 0.172 * rw, 0.40, 0.145]) for c in range(2)] for rw in range(5)]
    sym = dict(linthresh=1e-12)

    def put(row, fn, ylab, ylog=False, symlog=False, ylim=None):
        for c, (x, xl) in enumerate(((xr, xlab_r), (m, r"$M_{\rm enc}$  [M$_\odot$]"))):
            ax = axs[row][c]
            fn(ax, x)
            ax.set_xscale("log")
            if symlog:
                ax.set_yscale("symlog", **sym)
            elif ylog:
                ax.set_yscale("log")
            if ylim:
                ax.set_ylim(*ylim)
            if c == 0:
                ax.set_ylabel(ylab)
            if row == 4:
                ax.set_xlabel(xl)
            xs = x[selx]
            if xs.size:
                ax.set_xlim(np.nanmin(xs[xs > 0]), np.nanmax(xs))
            if c == 0 and np.isfinite(Rst) and Rst > 0:
                ax.axvline(1.0, color=FAINT, ls=":", lw=0.9)
            if c == 1 and np.isfinite(Mst) and Mst > 0:
                ax.axvline(Mst, color=FAINT, ls=":", lw=0.9)

    def line(ax, x, y, mask=None, **kw):
        y = np.asarray(y, float).copy()
        if mask is not None:
            y = np.where(mask, y, np.nan)
        plot_xy(ax, x, y, **kw)

    def f1(ax, x):
        line(ax, x, g("cv_nab"), sub, color=INK, lw=1.2, label=r"$\nabla$ (MLT, after HS blend)")
        line(ax, x, g("cv_nab_raw"), sub, color=INK, lw=0.7, ls=(0, (1, 1.5)), label=r"$\nabla$ raw")
        line(ax, x, g("cv_nad_eq"), sub, color="#c0392b", lw=1.0, label=r"$\nabla_{\rm ad}$ equil.")
        line(ax, x, g("cv_nad_fr"), sub, color="#2471a3", lw=1.0, ls=(0, (6, 1.6)), label=r"$\nabla_{\rm ad}$ frozen")
        line(ax, x, g("cv_nab_L"), sub, color="#1e8449", lw=1.0, ls=(0, (3, 1.4)), label=r"$\nabla_{\rm L}$ (Ledoux)")
        line(ax, x, g("cv_nab_rad"), sub, color="#7d3c98", lw=0.8, ls=(0, (7, 1.6, 1.4, 1.6)), label=r"$\nabla_{\rm rad}$")
        ax.set_ylim(-1.0, 4.0)
        if ax is axs[0][0]:
            ax.legend(fontsize=FS - 2, ncol=2, loc="upper right")
    put(0, f1, r"$\nabla$ (Mach $<0.5$ only)")

    def f2(ax, x):
        line(ax, x, g("cv_N2_eq"), sub, color="#c0392b", lw=1.1, label=r"$N^2$ equil. (Schwarzschild)")
        line(ax, x, g("cv_N2_L"), sub, color="#1e8449", lw=1.1, ls=(0, (3, 1.4)), label=r"$N^2_{\rm L}$ frozen Ledoux")
        line(ax, x, g("cv_N2_S"), sub, color="#2471a3", lw=0.9, ls=(0, (6, 1.6)), label=r"$N^2$ frozen Schwarzschild")
        ax.axhline(0, color=FAINT, lw=0.5)
        if ax is axs[1][0]:
            ax.legend(fontsize=FS - 2, loc="upper right")
    put(1, f2, r"$N^2$  [s$^{-2}$] (Mach $<0.5$ only)", symlog=True)

    def f3(ax, x):
        line(ax, x, g("cv_Fconv_Ftot"), None, color=INK, lw=1.2, label=r"$F_{\rm conv}/F_{\rm tot}$ (as used)")
        line(ax, x, g("cv_Fpot_Ftot"), None, color="#c0392b", lw=1.0, ls=(0, (6, 1.6)), label=r"gate-free MLT flux$/F_{\rm tot}$")
        if ax is axs[2][0]:
            ax.legend(fontsize=FS - 2, loc="upper right")
    put(2, f3, r"$F_{\rm conv}/F_{\rm tot}$", symlog=True)

    def f4(ax, x):
        line(ax, x, g("cv_w_tau"), None, color=INK, lw=1.1, label=r"$w_\tau$")
        line(ax, x, g("cv_w_inf"), None, color="#2471a3", lw=1.0, ls=(0, (6, 1.6)), label=r"$w_{\rm inf}$")
        line(ax, x, g("cv_w_th"), None, color="#c0392b", lw=1.0, ls=(0, (3, 1.4)), label=r"$w_{\rm th}$ ($\nabla_{\rm ad}<10^{29}$)")
        line(ax, x, g("cv_hs_w"), None, color="#1e8449", lw=0.9, ls=(0, (1.1, 1.3)), label=r"HS blend weight")
        ax.set_ylim(-0.05, 1.1)
        if ax is axs[3][0]:
            ax.legend(fontsize=FS - 2, ncol=2, loc="center right")
    put(3, f4, "MLT gates")

    def f5(ax, x):
        line(ax, x, g("cv_logt_settle"), None, color=INK, lw=1.1, label=r"$t_{\rm settle}=H_P/|u|$")
        line(ax, x, g("cv_logt_N"), sub, color="#c0392b", lw=1.0, ls=(0, (6, 1.6)), label=r"$1/\sqrt{|N^2_{\rm L}|}$")
        line(ax, x, g("cv_logt_turn"), None, color="#1e8449", lw=1.0, ls=(0, (3, 1.4)), label=r"$\ell/v$")
        line(ax, x, g("cv_logt_chem"), None, color="#7d3c98", lw=1.0, ls=(0, (7, 1.6, 1.4, 1.6)), label=r"$t_{\rm chem}$")
        if ax is axs[4][0]:
            ax.legend(fontsize=FS - 2, ncol=2, loc="upper right")
    put(4, f5, r"$\log_{10}$ timescale [s]")
    ax2 = axs[4][0].twinx()
    ax2.plot(xr[selx][np.argsort(xr[selx])], g("cv_mach")[selx][np.argsort(xr[selx])], color=FAINT, lw=0.8)
    ax2.set_ylabel("Mach", color=FAINT)
    ax2.axhline(0.5, color=FAINT, ls=":", lw=0.5)
    pdf.savefig(fig)
    plt.close(fig)
    return True


def make_header(run, sel):
    info = run.info
    nz = [int(f["n"]) for f in run.meta] or [0]
    commit = git_commit()
    rr = [run.frame(i)["logr"] for i in (sel[:1] + sel[-1:])] if sel else []
    rng = ""
    if rr:
        lo = np.nanmin(np.concatenate(rr))
        hi = np.nanmax(np.concatenate(rr))
        rng = rf"; $r$ = $10^{{{lo:.1f}}}$–$10^{{{hi:.1f}}}$ R$_\odot$"
    meshs = f"{nz[0]} zones" if min(nz) == max(nz) else f"{min(nz)}–{max(nz)} zones"
    t = run.s("t_yr")
    gen = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    status = info.get("status", "")
    desc = info.get("desc", "") or ""
    line2 = (f"mesh: {meshs}{rng};  frames: {run.nf} (steps {int(np.nanmin(run.s('step'))) if run.nf else 0}–"
             f"{int(np.nanmax(run.s('step'))) if run.nf else 0})" + (f";  status: {status}" if status else "")
             + f";  generated {gen};  abnfig @ titan.jl {commit}")

    def header(fig, subtitle):
        fig.text(0.80 / PAGE_W, 1 - 0.42 / PAGE_H, run.key, fontsize=FS + 7, fontweight="bold", va="baseline")
        fig.text(1 - 0.80 / PAGE_W, 1 - 0.42 / PAGE_H, subtitle, fontsize=FS + 1.5, va="baseline", ha="right",
                 style="italic", color=MUTED)
        dw = "\n".join(textwrap.wrap(desc, 135)[:3])
        fig.text(0.80 / PAGE_W, 1 - 0.56 / PAGE_H, dw, fontsize=FS, va="top", linespacing=1.25)
        fig.text(0.80 / PAGE_W, 1 - 1.02 / PAGE_H, line2, fontsize=FS - 1.5, va="top", color=MUTED)
        fig.add_artist(Line2D([0.80 / PAGE_W, 1 - 0.80 / PAGE_W], [1 - 0.48 / PAGE_H] * 2, color=INK, lw=0.8))
    return header


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run")
    ap.add_argument("--panel", default=DEFAULT_PANEL, help="explorer data directory P (default: %(default)s)")
    ap.add_argument("--out", default=None, help="output PDF (default: P/RUN/abnfig.pdf)")
    ap.add_argument("--times", default="auto", help="'auto' or comma-separated t - t0 [yr]")
    ap.add_argument("--xaxis", default="abn", choices=["abn", "r", "m"],
                    help="abn: left column vs r, right column vs M_enc (ABN02 layout); r / m: all panels vs r / M_enc")
    ap.add_argument("--cvframe", type=int, default=None, help="frame index for the convection page (only runs with cv_* fields; default: strongest unstable layer)")
    ap.add_argument("--n", type=int, default=7, help="number of output times for --times auto (default 7, max 8)")
    a = ap.parse_args(argv)
    run = Run(a.panel, a.run)
    if run.nf == 0:
        sys.exit(f"abnfig: {a.run}: no frames")
    tsel = None if a.times == "auto" else [float(v) for v in a.times.split(",") if v.strip()]
    sel, t0 = select_times(run, max(1, min(a.n, 8)), tsel)
    ps = protostar_mask(run)
    pro = [i for i in sel if ps[i]]
    pre = [i for i in sel if not ps[i]]
    cols = time_colours(len(pro))
    styles = {}
    for k, i in enumerate(pro):
        styles[i] = {"c": cols[k], "d": DASHES[k % len(DASHES)], "lw": 1.15, "pro": True}
    for k, i in enumerate(pre):
        styles[i] = {"c": PRE_COLS[k % 2] if pro else time_colours(len(pre))[k],
                     "d": PRE_DASHES[k % 2] if pro else DASHES[k % len(DASHES)], "lw": 1.15, "pro": False}
    out = a.out or os.path.join(run.dir, "abnfig.pdf")
    header = make_header(run, sel)
    with PdfPages(out) as pdf:
        page_profiles(pdf, run, sel, t0, a.xaxis, styles, header)
        page_history(pdf, run, sel, t0, styles, header)
        if any(f.startswith("cv_") for f in run.fields):          # CONVDIAG=1 runs only: runs without cv_* fields are unchanged
            page_convection(pdf, run, sel, t0, header, a.cvframe)
        d = pdf.infodict()
        d["Title"] = f"titan.jl run {run.key}: ABN02-style profiles"
        d["Author"] = "abnfig.py (titan.jl)"
    print(f"abnfig: wrote {out}  ({len(sel)} times: frames {sel})")


if __name__ == "__main__":
    main()
