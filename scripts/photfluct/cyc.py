# sawtooth statistics of the emergent luminosity: collapses = L_em drop by > fac within one record
import sys, numpy as np
f = sys.argv[1]; fac = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0; smin = float(sys.argv[3]) if len(sys.argv) > 3 else 0
d = np.genfromtxt(f, delimiter=',', names=True); d = d[d['step'] >= smin]
L = d['Lem']; t = d['t']; st = d['step']
idx = [i for i in range(1, len(L)) if L[i - 1] > 0 and L[i] > 0 and L[i - 1] / L[i] > fac and L[i-1] > 50]
print(f, 'records', len(d), 'steps %d-%d' % (st[0], st[-1]), 't span %.3e s' % (t[-1] - t[0]))
for j, i in enumerate(idx):
    s = ' collapse at step %6d t=%.4e  L %.3g -> %.3g  kph %d->%d js %d  rph %.2f M %.4f' % (st[i], t[i], L[i-1], L[i], d['kph'][i-1], d['kph'][i], d['js'][i], d['rph'][i], d['Mst'][i])
    if j > 0:
        k = idx[j - 1]; s += '  | period %d steps %.3e s, <dt> %.3e s, drph %.2f Rsun' % (st[i] - st[k], t[i] - t[k], (t[i] - t[k]) / (st[i] - st[k]), d['rph'][i] - d['rph'][k])
        seg = slice(k, i); Ls = L[seg]
        tw = np.diff(t[k:i + 1]); s += '  L_em mean(t) %.3g max %.3g min %.3g' % (np.sum(L[k:i] * tw) / np.sum(tw), Ls.max(), Ls[Ls > 0].min())
    print(s)
print('median Lacc %.3g, median L_ph %.3g, time-mean L_em %.3g, max/min L_em %.3g/%.3g' % (np.median(d['Lacc']), np.median(d['Lph']),
      np.sum(L[:-1] * np.diff(t)) / (t[-1] - t[0]), L.max(), L[L > 0].min()))
