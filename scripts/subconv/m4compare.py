#!/home/tabel/.local/share/anaconda3/bin/python3
"""m4compare.py [--tag TAG] [--exports DIR] [--out DIR] [--window 20] [--smin 290] CASE...

Compare M4 accretion-history cases (first case given as B0 is the reference; B0 = production replay zie1rp, steps >= smin).
Reads <exports>/<case>/{track.csv,frames.json,frames.f32} written by m4export.sh (runrec.jl export).
Writes <out>/compare_<tag>.txt (tables at fixed M* and fixed t-t0, relative differences vs B0, spread across cases),
<out>/compare_<tag>.png (vs M*) and <out>/compare_<tag>_t.png (vs t-t0).

PROXIES (flagged [P] in the tables): s_ps = ln(T^1.5/rho) of the post-shock zone (cgs, arbitrary constant; ideal-gas monatomic
entropy proxy, ignores H2 dissociation/ionisation/radiation); Mdot_sh = 4 pi r^2 rho |u| from the frame zone just outside
the shock (cell-centred, not interface); u_pk = max|u| infall in [R*, 3R*] zones; sigma(log L), sigma(log Teff) = rolling std over
`window` consecutive frames (depends on snapshot cadence); <X_D>_in = mass-weighted over zones inside R* (from logm_o differences).
"""
import sys, os, json, argparse, datetime
import numpy as np

Rs = 6.957e10; Ms = 1.989e33; yr = 3.156e7
T0 = 242378.256758
MT = [0.3, 0.5, 1, 2, 3, 5, 10, 20, 30]
TT = [5, 10, 20, 50, 100, 200, 500, 1000, 3000]   # t - t0 [yr]

def load(root, case, smin):
    d = os.path.join(root, case)
    tp = os.path.join(d, 'track.csv')
    if not os.path.exists(tp): return None
    tr = np.genfromtxt(tp, delimiter=',', names=True)
    if tr.ndim == 0: tr = tr.reshape(1)
    J = json.load(open(os.path.join(d, 'frames.json')))
    fn = [k for k, _ in J['fields']]
    D = np.fromfile(os.path.join(d, 'frames.f32'), dtype='<f4')
    sel = np.where(tr['step'] >= (smin if case == 'B0' or True else 0))[0]
    # restarted lanes begin at step 290 anyway; smin applies to all
    return dict(tr=tr, J=J, fn=fn, D=D, sel=sel)

def frame_fields(c, i):
    f = c['J']['frames'][i]; n = f['n']; o = f['off'] // 4
    return {k: c['D'][o + j * n:o + (j + 1) * n].astype(float) for j, k in enumerate(c['fn'])}

def nan_div(a, b):
    with np.errstate(all='ignore'): return a / b

def series(c):
    tr = c['tr']; sel = c['sel']; N = len(sel)
    S = {}
    g = lambda k: tr[k][sel] if k in tr.dtype.names else np.full(N, np.nan)
    S['M'] = g('M_star'); S['dt'] = g('t_yr') - T0; S['step'] = g('step')
    S['Rstar'] = g('R_star'); S['rph'] = g('r_ph'); S['Teff'] = g('T_eff_lum'); S['Lph'] = g('L_ph')
    S['Mdot1e4'] = g('Mdot_1e4'); S['umin'] = g('u_min'); S['L_D'] = g('L_D'); S['L_pp'] = g('L_pp'); S['L_CNO'] = g('L_CNO')
    S['Tc'] = g('T_c'); S['lognc'] = g('logn_c')
    for k, kk in [('cv_Fmax', 'cv_FconvFtot_max'), ('cv_rin', 'cv_rin'), ('cv_rout', 'cv_rout'), ('cv_mlayer', 'cv_mlayer'),
                  ('cv_Et', 'cv_Et'), ('cv_PtPps', 'cv_PtP_ps')]:
        S[k] = g(kk)
    S['cv_ext'] = S['cv_rout'] - S['cv_rin']
    # profile-derived
    keys = ['Mdot_sh', 'u_pk', 's_ps', 'XD_in', 'T_ps', 'rho_ps']
    for k in keys: S[k] = np.full(N, np.nan)
    for a, i in enumerate(sel):
        F = frame_fields(c, i)
        r = 10 ** F['logr'] * Rs; Rst = S['Rstar'][a] * Rs
        if not np.isfinite(Rst): continue
        js = np.argmax(r >= Rst * (1 - 1e-9))        # shock zone (first zone with infall), R* = its radius
        jo = js + 1 if js + 1 < len(r) else js
        rho = 10 ** F['logrho']; u = F['u'] * 1e5
        S['Mdot_sh'][a] = 4 * np.pi * r[jo] ** 2 * rho[jo] * abs(u[jo]) * yr / Ms
        m3 = (r >= Rst) & (r <= 3 * Rst) & (u < 0)
        if m3.any(): S['u_pk'][a] = -u[m3].min() / 1e5
        jp = max(js - 1, 0)
        T = 10 ** F['logT'][jp]; S['T_ps'][a] = T; S['rho_ps'][a] = rho[jp]
        S['s_ps'][a] = np.log(T ** 1.5 / rho[jp])
        if 'logXD' in F and np.isfinite(F['logXD']).any():
            mo = 10 ** F['logm_o']; dm = np.diff(np.concatenate([[0.], mo]))
            XD = 10 ** F['logXD']; XD = np.where(np.isfinite(XD), XD, 0.0)
            ins = np.arange(len(r)) < js
            if ins.any() and dm[ins].sum() > 0: S['XD_in'][a] = (XD[ins] * dm[ins]).sum() / dm[ins].sum()
    return S

def rolling_std(x, w):
    out = np.full(len(x), np.nan)
    for i in range(len(x)):
        seg = x[max(0, i - w + 1):i + 1]
        if len(seg) >= max(5, w // 2) and np.isfinite(seg).all(): out[i] = np.std(seg)
    return out

def first_onset(M, L, Lph, thr=1e-3):
    on = np.isfinite(L) & np.isfinite(Lph) & (Lph > 0) & (L > thr * Lph)
    for i in range(len(on)):
        if on[i:i + 3].all() and len(on[i:i + 3]) == 3 or (i == len(on) - 1 and on[i]): return M[i]
    return np.nan

def at_x(x, y, targets):
    """first-crossing linear interpolation of y at x = target; NaN if x never reaches target or starts above it."""
    out = np.full(len(targets), np.nan)
    ok = np.isfinite(x)
    xs = np.maximum.accumulate(np.where(ok, x, -np.inf))
    for k, t in enumerate(targets):
        idx = np.where(xs >= t)[0]
        if len(idx) == 0 or idx[0] == 0: continue
        i = idx[0]; x0, x1 = x[i - 1], x[i]
        if not (np.isfinite(x0) and np.isfinite(x1)) or x1 == x0: continue
        w = (t - x0) / (x1 - x0)
        y0, y1 = y[i - 1], y[i]
        out[k] = y0 + w * (y1 - y0) if np.isfinite(y0) and np.isfinite(y1) else np.nan
    return out

# (key, label, proxy?, transform, log-plot?)
QUANT = [
    ('Rstar', 'R* [Rsun]', 0, None, 1), ('rph', 'r_ph [Rsun]', 0, None, 1), ('Teff', 'T_eff [K]', 0, None, 1), ('Lph', 'L_ph [Lsun]', 0, None, 1),
    ('Mdot1e4', 'Mdot(1e4 Rsun) [Msun/yr]', 0, None, 1), ('Mdot_sh', 'Mdot_shock 4pi r2 rho|u| [Msun/yr]', 1, None, 1),
    ('umin', 'u_min [km/s]', 0, None, 0), ('u_pk', 'max|u| in [R*,3R*] [km/s]', 1, None, 1),
    ('s_ps', 's_ps = ln(T^1.5/rho) post-shock zone', 1, None, 0), ('T_ps', 'T post-shock zone [K]', 0, None, 1), ('rho_ps', 'rho post-shock [g/cc]', 0, None, 1),
    ('XD_in', '<X_D> mass-weighted inside R*', 1, None, 1), ('L_D', 'L_D [Lsun]', 0, None, 1),
    ('sigL', 'sigma(log10 L_ph), rolling', 1, None, 1), ('sigT', 'sigma(log10 Teff), rolling', 1, None, 1),
    ('Tc', 'T_c [K]', 0, None, 1),
    ('cv_Fmax', 'max Fconv/F in layer', 0, None, 0), ('cv_ext', 'layer extent r_out-r_in [R*]', 0, None, 0), ('cv_mlayer', 'layer mass [Msun]', 0, None, 0),
    ('cv_Et', 'integrated E_t (runrec cv_Et)', 0, None, 1), ('cv_PtPps', 'P_t/P post-shock', 0, None, 1),
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cases', nargs='+')
    ap.add_argument('--tag', default=None)
    ap.add_argument('--exports', default='/zpool/nvme/data/tabel_scratch/titan_popiii/subconv/m4/exports')
    ap.add_argument('--out', default=None)
    ap.add_argument('--window', type=int, default=20)
    ap.add_argument('--smin', type=int, default=290)
    ap.add_argument('--ref', default='B0', help='reference case for the relative differences')
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    out = a.out or os.path.join(here, '..', '..', 'docs', 'logs', 'subconv', 'm4')
    os.makedirs(out, exist_ok=True)
    tag = a.tag or 'partial_' + datetime.date.today().strftime('%Y%m%d')
    cases = list(a.cases); data = {}; skipped = []
    for cn in cases:
        c = load(a.exports, cn, a.smin)
        if c is None or len(c['sel']) == 0 or (c['sel'].size and not (c['tr']['step'][c['sel']] >= a.smin).any()):
            skipped.append(cn); continue
        c['sel'] = np.where(c['tr']['step'] >= a.smin)[0]
        if len(c['sel']) == 0: skipped.append(cn); continue
        S = series(c)
        S['sigL'] = rolling_std(np.log10(np.where(S['Lph'] > 0, S['Lph'], np.nan)), a.window)
        S['sigT'] = rolling_std(np.log10(np.where(S['Teff'] > 0, S['Teff'], np.nan)), a.window)
        S['onset_pp'] = first_onset(S['M'], S['L_pp'], S['Lph']); S['onset_CNO'] = first_onset(S['M'], S['L_CNO'], S['Lph'])
        data[cn] = S
    names = [c for c in cases if c in data]
    L = []; P = L.append
    P(f'M4 accretion-history comparison  tag={tag}  window={a.window} frames  steps >= {a.smin}  t0 = {T0} yr')
    P(f'cases with data: {names}; without snapshots yet / skipped: {skipped}')
    for cn in names:
        S = data[cn]; P(f'  {cn}: {len(S["M"])} frames, steps {int(S["step"][0])}..{int(S["step"][-1])}, M* {S["M"][0]:.3f}..{np.nanmax(S["M"]):.3f}, t-t0 {S["dt"][0]:.2f}..{S["dt"][-1]:.2f} yr, '
                    f'median dstep {np.median(np.diff(S["step"])) if len(S["step"])>1 else float("nan"):.0f}')
    P('[P] = proxy / approximate quantity. s_ps: ln(T^1.5/rho) (cgs) of the zone just inside the shock (not a true entropy: no ionisation/H2/radiation).')
    P('Mdot_shock: cell-centred 4 pi r^2 rho |u| at the zone just outside R* (R* = first zone with u<-2 km/s). max|u|: max infall speed in [R*,3R*].')
    P('Onset masses: first M* where L_pp or L_CNO > 1e-3 L_ph for 3 consecutive frames (NaN: not reached). Fluctuations depend on snapshot cadence per case.')
    P('cv_* columns exist only for CONVDIAG=1 lanes (NaN for B0). Rel diff = (case-B0)/|B0| ; spread = (max-min)/|B0| over all cases (incl. B0).')
    P('')
    ref = a.ref if a.ref in data else None
    def table(xname, xkey, targets, fmt_x):
        P('=' * 110); P(f'TABLE vs {xname}'); P('=' * 110)
        for key, lab, prox, _, _ in QUANT:
            vals = {cn: at_x(data[cn][xkey], data[cn][key], targets) for cn in names}
            if all(np.isnan(v).all() for v in vals.values()): P(f'-- {lab}{" [P]" if prox else ""}: all NaN'); continue
            P(f'-- {lab}{" [P]" if prox else ""}')
            P('   case   ' + ' '.join(f'{fmt_x(t):>10s}' for t in targets))
            for cn in names: P(f'   {cn:7s}' + ' '.join(f'{v:10.4g}' for v in vals[cn]))
            if ref:
                for cn in names:
                    if cn == ref: continue
                    with np.errstate(all='ignore'): rd = (vals[cn] - vals[ref]) / np.abs(vals[ref])
                    P(f'   d{cn:6s}' + ' '.join(f'{v:10.3g}' for v in rd))
                V = np.array([vals[cn] for cn in names])
                with np.errstate(all='ignore'): sp = (np.nanmax(V, 0) - np.nanmin(V, 0)) / np.abs(vals[ref])
                P('   spread ' + ' '.join(f'{v:10.3g}' for v in sp))
        P('')
    import warnings; warnings.filterwarnings('ignore')
    table('M* [Msun]', 'M', MT, lambda t: f'{t:g}')
    table('t - t0 [yr]', 'dt', TT, lambda t: f'{t:g}')
    P('=' * 110); P('ONSET MASSES [Msun] (first M* with L > 1e-3 L_ph)'); P('=' * 110)
    P('   case    pp         CNO')
    for cn in names: P(f'   {cn:7s} {data[cn]["onset_pp"]:<10.4g} {data[cn]["onset_CNO"]:<10.4g}')
    P('')
    open(os.path.join(out, f'compare_{tag}.txt'), 'w').write('\n'.join(L) + '\n')
    # ---- plots
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    for xkey, xl, suf, logx in [('M', 'M* [Msun]', '', True), ('dt', 't - t0 [yr]', '_t', True)]:
        nq = len(QUANT); nc = 4; nr = (nq + nc - 1) // nc
        fig, ax = plt.subplots(nr, nc, figsize=(4.2 * nc, 3.0 * nr), squeeze=False)
        for q, (key, lab, prox, _, lg) in enumerate(QUANT):
            A = ax[q // nc][q % nc]
            for cn in names:
                S = data[cn]; y = S[key]
                if np.isfinite(y).any(): A.plot(S[xkey], y, lw=1, label=cn)
            if lg and any(np.isfinite(data[c][key]).any() and (np.nanmin(np.where(np.isfinite(data[c][key]),data[c][key],np.nan)) > 0) for c in names): A.set_yscale('log')
            if logx: A.set_xscale('log')
            A.set_title(lab + (' [proxy]' if prox else ''), fontsize=8); A.tick_params(labelsize=7)
            if q // nc == nr - 1: A.set_xlabel(xl, fontsize=8)
        for q in range(nq, nr * nc): ax[q // nc][q % nc].axis('off')
        ax[0][0].legend(fontsize=7)
        fig.suptitle(f'M4 comparison {tag} (steps>={a.smin}); proxies marked', fontsize=10)
        fig.tight_layout(); fig.savefig(os.path.join(out, f'compare_{tag}{suf}.png'), dpi=110); plt.close(fig)
    # D profile panel
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.2))
    for A, Mt in zip(ax, [1, 3, 10]):
        for cn in names:
            S = data[cn]; c = None
            c = load(a.exports, cn, a.smin); sel = c['sel']
            j = np.where(c['tr']['M_star'][sel] >= Mt)[0]
            if len(j) == 0: continue
            F = frame_fields(c, sel[j[0]])
            if np.isfinite(F['logXD']).any(): A.plot(F['logm'] - np.log10(max(c['tr']['M_star'][sel][j[0]], 1e-9)), F['logXD'], label=cn)
        A.set_title(f'log X_D profile at M*>={Mt}', fontsize=9); A.set_xlabel('log m/M*'); A.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out, f'compare_{tag}_D.png'), dpi=100); plt.close(fig)
    print('wrote', os.path.join(out, f'compare_{tag}.txt'))

if __name__ == '__main__': main()
