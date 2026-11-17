# compare restarts: per-step photospheric fluctuation, dt, cost, physics over a common physical time window
#   cmp.py <tmax_rel or 0> run1 run2 ...      (hist_<run>.csv from ana.jl; perf.csv and seg.out from the run dir)
import sys, numpy as np, re
R = '/zpool/nvme/data/tabel_scratch/titan_popiii/runs/'
tmax = float(sys.argv[1]); runs = sys.argv[2:]
def alt(x):  # fraction of sign alternations of successive increments (1/2 for noise, -> 1 for period-2 ringing)
    dx = np.diff(x); dx = dx[dx != 0]
    return np.mean(np.sign(dx[1:]) != np.sign(dx[:-1])) if len(dx) > 2 else np.nan
def p2(x):  # rms period-2 (zig-zag) component: x_k - (x_{k-1}+x_{k+1})/2
    return np.sqrt(np.mean((x[1:-1] - 0.5 * (x[:-2] + x[2:])) ** 2)) if len(x) > 3 else np.nan
def ac1(x):
    dx = np.diff(x); dx = dx - dx.mean()
    return np.sum(dx[1:] * dx[:-1]) / np.sum(dx * dx) if len(dx) > 2 else np.nan
for run in runs:
    d = np.genfromtxt('hist_%s.csv' % run, delimiter=',', names=True)
    t0 = d['t'][0]; tr = d['t'] - t0
    if tmax > 0: d = d[tr <= tmax]; tr = d['t'] - t0
    lL = np.log10(np.maximum(d['Lem'], 1e-30)); lP = np.log10(np.maximum(d['Lph'], 1e-30)); lT = np.log10(d['Trad_ph']); ldt = np.log10(d['dt'])
    tw = np.diff(d['t'])
    meanL = np.sum(d['Lem'][1:] * tw) / np.sum(tw)
    p = np.genfromtxt(R + run + '/perf.csv', delimiter=',', names=True)
    p = p[p['nsteps'] > 0]
    sps = np.sum(p['t_step']) / np.sum(p['nsteps'])
    nfail = (p['fail_tr'][-1] + p['fail_bdf'][-1]) - (p['fail_tr'][0] + p['fail_bdf'][0])
    nst = p['step'][-1] - p['step'][0]
    seg = open(R + run + '/seg.out', errors='replace').read()
    nrej = len(re.findall(r'DTERR reject jstep (\d+): dt [0-9.e+-]+ err', seg)); nnf = len(re.findall(r'DTERR reject jstep \d+: dt [0-9.e+-]+ Newton failure', seg))
    print('%-7s steps %4d (%d-%d) t %.3e s | <dt> %.3e med dt %.3e | L_em t-mean %.3g min %.3g max %.3g | dlogL_em/step rms %.3f alt %.2f ac1 %+.2f zz %.4f | dlogL_ph rms %.3f alt %.2f zz %.4f | dlogTrad_ph rms %.4f alt %.2f zz %.5f | dlogdt alt %.2f ac1 %+.2f | M* %.5f->%.5f R* %.2f->%.2f rph %.2f->%.2f kph %d->%d | s/step %.2f, s/(1e5 s) %.1f, err-rej %d Newton-fail %d per100 %.1f/%.1f' % (
        run, len(d), d['step'][0], d['step'][-1], tr[-1], tr[-1] / max(len(d) - 1, 1), np.median(d['dt']), meanL, d['Lem'].min(), d['Lem'].max(),
        np.sqrt(np.mean(np.diff(lL) ** 2)), alt(lL), ac1(lL), p2(lL), np.sqrt(np.mean(np.diff(lP) ** 2)), alt(lP), p2(lP), np.sqrt(np.mean(np.diff(lT) ** 2)), alt(lT), p2(lT), alt(ldt), ac1(ldt),
        d['Mst'][0], d['Mst'][-1], d['Rst'][0], d['Rst'][-1], d['rph'][0], d['rph'][-1], d['kph'][0], d['kph'][-1],
        sps, sps * (len(d) - 1) / tr[-1] * 1e5, nrej, nnf, 100 * nrej / max(nst, 1), 100 * nnf / max(nst, 1)))
