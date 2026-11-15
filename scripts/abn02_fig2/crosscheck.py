# M(r) from integrating panel A's n(r) (rho = n mu m_H) against panel B's M_enc(r), for three readings of mu; and the
# headline numbers of the final-output profile.  usage: python3 -I crosscheck.py <abn02_final_fig2.csv>
import sys
import numpy as np

rows = [l.strip().split(',') for l in open(sys.argv[1]) if not l.startswith('#') and not l.startswith('curve')]
def cur(k):
    a = np.array([[float(r[1]), float(r[2])] for r in rows if r[0] == k]); a = a[np.argsort(a[:, 0], kind='stable')]
    _, i = np.unique(a[:, 0], return_index=True); return a[i, 0], a[i, 1]
lr, ln = cur('n_r'); lrB, lMB = cur('M_r'); lMH, lfH2 = cur('fH2_M'); lMI, lfHI = cur('fHI_M')
lMT, lT = cur('T_M'); lMv, v = cur('vr_M'); lMc, mcs = cur('mcs_M'); lrd, lrhod = cur('dm_r')
pc, mH, Msun, X, Y = 3.0857e18, 1.6735575e-24, 1.98847e33, 0.76, 0.24
x = np.linspace(-6.5, lr[-1], 40001); r = pc * 10 ** x
n = 10 ** np.interp(x, lr, ln)                                        # (constant inside the first point, r < 1e-6 pc)
MB = lambda xx: np.interp(xx, lrB, lMB)                               # panel B, log M at log r
def mu_of(lm):                                                        # mean particle mass [m_H] for the panel-C composition (neutral)
    fh2 = 10 ** np.interp(lm, lMH, lfH2); fhi = 10 ** np.interp(lm, lMI, lfHI)
    return 1.0 / (fhi + fh2 / 2 + Y / 4)
def Menc(mu):
    dm = 4 * np.pi * r ** 2 * n * mu * mH; M = np.concatenate([[0], np.cumsum(0.5 * (dm[1:] + dm[:-1]) * np.diff(r))])
    return M + 4 / 3 * np.pi * r[0] ** 3 * n[0] * mu[0] * mH if np.ndim(mu) else M + 4 / 3 * np.pi * r[0] ** 3 * n[0] * mu * mH
mu_var = mu_of(MB(np.clip(x, lrB[0], None)))
cases = {'mu=1.22 (atomic, as abn02*)': 1.22, 'mu(f_H2(M)) from panel C, neutral': mu_var, 'mu=1 (n = rho/m_H)': 1.0,
         'mu=2.27 (fully molecular)': 2.27}
print('mu(panel C) at the centre %.3f, at 1 Msun %.3f, at 1e3 Msun %.3f' % (mu_of(-2), mu_of(0), mu_of(3)))
xs = np.array([-5, -4, -3, -2, -1, 0, 1, 2, 2.28])
print('log r[pc]   ' + ' '.join('%7.2f' % q for q in xs))
print('panel B     ' + ' '.join('%7.3f' % MB(q) if q >= lrB[0] else '    -  ' for q in xs))
for k, mu in cases.items():
    lM = np.log10(Menc(mu) / Msun); d = np.interp(xs, x, lM) - np.array([MB(q) for q in xs])
    sel = (x >= lrB[0]) & (x <= lrB[-1]); dd = np.interp(x[sel], x, lM) - MB(x[sel])
    print('%-34s' % k, ' '.join('%+7.3f' % q if xx >= lrB[0] else '    -  ' for q, xx in zip(d, xs)),
          '| d log M over B: mean %+.3f rms %.3f max|%.3f|' % (dd.mean(), np.sqrt((dd ** 2).mean()), abs(dd).max()))
# headline numbers
print('n_c = 10^%.3f cm^-3 at r = 10^%.2f pc (innermost point), max log n %.3f' % (ln[0], lr[0], ln.max()))
for m in (0.0, 3.0, lMB[-1]):
    print('M = 10^%.3f Msun at r = 10^%.3f pc = %.3g pc = %.3g au' % (m, np.interp(m, lMB, lrB), 10 ** np.interp(m, lMB, lrB), 10 ** np.interp(m, lMB, lrB) * 206265))
k = np.argmax(ln < ln[0] - 1); print('n falls 1 dex below n_c at r = 10^%.2f pc, M_enc(B) there 10^%.2f Msun' % (lr[k], MB(lr[k])))
print('panel B mass range 10^%.3f .. 10^%.3f Msun (r 10^%.3f .. 10^%.3f pc); C/D/E mass range 10^%.3f .. 10^%.3f' % (lMB[0], lMB[-1], lrB[0], lrB[-1], lMT[0], lMT[-1]))
print('T: %.0f .. %.0f K; T at M_min %.0f K, at 1 Msun %.0f K' % (10 ** lT.min(), 10 ** lT.max(), 10 ** lT[0], 10 ** np.interp(0, lMT, lT)))
print('f_H2: %.3g (centre) .. %.3g; f_HI at centre %.3g; f_H2+f_HI over the range: %.3f..%.3f' % (10 ** lfH2[0], 10 ** lfH2.min(), 10 ** lfHI[0],
      min(10 ** np.interp(q, lMH, lfH2) + 10 ** np.interp(q, lMI, lfHI) for q in np.linspace(-1.99, 4.15, 200)),
      max(10 ** np.interp(q, lMH, lfH2) + 10 ** np.interp(q, lMI, lfHI) for q in np.linspace(-1.99, 4.15, 200))))
i = np.argmin(v); print('v_r min %.3f km/s at M = 10^%.2f Msun' % (v[i], lMv[i]))
mm = np.linspace(lMc[0], lMc[-1], 400); cs = -np.interp(mm, lMc, mcs); vv = np.interp(mm, lMv, v); unc = np.interp(mm, lMc, mcs) <= -2.999
mach = -vv / cs; j = np.argmax(np.where(unc, -1, mach))
print('Mach |v|/c_s max %.2f at M = 10^%.2f (c_s %.2f km/s, frame-clipped c_s excluded: %s); Mach at 0.01, 0.1, 1 Msun: %s' % (mach[j], mm[j], cs[j],
      'none' if not unc.any() else '10^%.2f..10^%.2f' % (mm[unc].min(), mm[unc].max()), ', '.join('%.2f' % (-np.interp(q, lMv, v) / -np.interp(q, lMc, mcs)) for q in (-2, -1, 0))))
print('DM final: log rho_DM [GeV cm^-3] %.2f at 10^%.2f pc .. %.2f at 10^%.2f pc' % (lrhod[0], lrd[0], lrhod[-1], lrd[-1]))
