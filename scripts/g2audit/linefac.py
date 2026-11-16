# Line (+ OPAL-only) content of the table opacity in the network domain: table / G2-continuum on HM88 populations (out_p2w), with s from out_nd3.
import numpy as np
def load(d, fn):
    L = open(d + '/' + fn + '.txt').read().split('\n'); F = L[0].split(); LR = np.array(L[1].split(), float); LT = np.array(L[2].split(), float)
    A = np.fromfile(d + '/' + fn + '.bin').reshape(len(LT), len(LR), len(F)); return {f: A[:, :, k].T for k, f in enumerate(F)}, LR, LT
g, LR, LT = load('out_p2w', 'grid'); h, LR2, LT2 = load('out_nd3', 'gridnd')
s = ((h['Gam'] / 0.3) ** 4 + (h['Deg'] / 0.1) ** 4 + (h['Pneu'] / 0.01) ** 4) ** 0.25
i2 = {round(x, 1): i for i, x in enumerate(LR2)}; j2 = {round(x, 2): j for j, x in enumerate(LT2)}
fR = []; fP = []; pts = []
for i, lr in enumerate(LR):
    for j, lt in enumerate(LT):
        if round(lr, 1) not in i2 or round(lt, 2) not in j2: continue
        ss = s[i2[round(lr, 1)], j2[round(lt, 2)]]
        if not (ss < 1): continue
        a = np.log10(g['kR_t'][i, j] / g['kR_p'][i, j]); b = np.log10(g['kP_t'][i, j] / g['kP_p'][i, j])
        if np.isfinite(a) and np.isfinite(b): fR.append(a); fP.append(b); pts.append((lr, lt))
fR = np.array(fR); fP = np.array(fP); pts = np.array(pts)
for lab, m in (('T < 5000 K', pts[:, 1] < np.log10(5000)), ('5000-6300 K', (pts[:, 1] >= np.log10(5000)) & (pts[:, 1] < 3.8)),
               ('6300 K - 5e4 K', (pts[:, 1] >= 3.8) & (pts[:, 1] < np.log10(5e4))), ('5e4 - 3e5 K', (pts[:, 1] >= np.log10(5e4)) & (pts[:, 1] < np.log10(3e5)))):
    print(f'{lab:16s} n={m.sum():5d}  log10(table/G2cont) kappa_R: p5 {np.percentile(fR[m],5):+.2f} p50 {np.percentile(fR[m],50):+.2f} p95 {np.percentile(fR[m],95):+.2f} max {fR[m].max():+.2f} | kappa_P: p5 {np.percentile(fP[m],5):+.2f} p50 {np.percentile(fP[m],50):+.2f} p95 {np.percentile(fP[m],95):+.2f} max {fP[m].max():+.2f}')
