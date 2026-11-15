# ABN02 Fig. 2 SECOND output (red dotted, "nine Myrs later"): M(r) from panel A's n(r) against panel B, the central plateau,
# headline numbers, and a comparison with the hand read-off tables of IC=abn02 (julia/popiii_m2c.jl, "9 Myr later").
# usage: python3 -I crosscheck_red.py <abn02_red_fig2.csv>
import sys
import numpy as np

rows = [l.strip().split(',') for l in open(sys.argv[1]) if not l.startswith('#') and not l.startswith('curve')]
def cur(k):
    a = np.array([[float(r[1]), float(r[2])] for r in rows if r[0] == k]); return a[:, 0], a[:, 1]
lr, ln = cur('n_r'); lrB, lMB = cur('M_r'); lMH, lfH2 = cur('fH2_M'); lMI, lfHI = cur('fHI_M'); lMT, lT = cur('T_M'); lMv, v = cur('vr_M')
pc, mH, Msun, Rsun = 3.0857e18, 1.6735575e-24, 1.98847e33, 6.957e10
MB = lambda x: np.interp(x, lrB, lMB)
def mu_of(lm):
    f2 = 10 ** np.interp(lm, lMH, lfH2); f1 = 10 ** np.interp(lm, lMI, lfHI); return 1.0 / (f1 + f2 / 2 + (1 - f1 - f2) / 4)
# plateau
kp = np.where(ln == ln[0])[0].max(); rp = lr[kp]
print('n(r) plateau: n = 10^%.4f to log r = %.3f pc (vertex %d; next vertex %.3f at %.3f); M_B(r_p) = 10^%.3f = %.4f Msun' %
      (ln[0], rp, kp, ln[kp + 1], lr[kp + 1], MB(rp), 10 ** MB(rp)))
for lab, d in (('T', (lMT, lT)), ('f_H2', (lMH, lfH2))):
    x, y = d; k = np.argmax(np.abs(y - y[0]) > 0.01); print('%s flat (within 0.01 dex of its central value 10^%.4f) to M = 10^%.2f Msun' % (lab, y[0], x[k]))
mup = mu_of(MB(rp)); rpc = pc * 10 ** rp
print('mu at the plateau (panel C) %.4f; uniform sphere with n_A mu m_H: M = %.4f Msun vs panel B %.4f (%.3f dex); panel B mean density inside r_p -> n = 10^%.3f at that mu' %
      (mup, 4 / 3 * np.pi * rpc ** 3 * 10 ** ln[0] * mup * mH / Msun, 10 ** MB(rp), MB(rp) - np.log10(4 / 3 * np.pi * rpc ** 3 * 10 ** ln[0] * mup * mH / Msun),
       np.log10(10 ** MB(rp) * Msun / (4 / 3 * np.pi * rpc ** 3) / (mup * mH))))
# M(r) from A
x = np.linspace(-6.5, lr[-1], 40001); r = pc * 10 ** x; n = 10 ** np.interp(x, lr, ln)
def Menc(mu):
    dm = 4 * np.pi * r ** 2 * n * mu * mH; M = np.concatenate([[0], np.cumsum(0.5 * (dm[1:] + dm[:-1]) * np.diff(r))])
    return M + 4 / 3 * np.pi * r[0] ** 3 * n[0] * (mu if np.ndim(mu) == 0 else mu[0]) * mH
muv = mu_of(MB(np.clip(x, lrB[0], None)))
xs = np.array([-2.5, -2.0, -1.5, -1.0, 0.0, 1.0, 2.0, 2.28])
print('log r[pc]        ' + ' '.join('%7.2f' % q for q in xs)); print('panel B log M    ' + ' '.join('%7.3f' % MB(q) for q in xs))
for lab, mu in (('mu = 1.22', 1.22), ('mu(panel C)', muv), ('mu = 1', 1.0)):
    lM = np.log10(Menc(mu) / Msun); d = np.interp(xs, x, lM) - MB(xs)
    sel = (x >= lrB[0]); dd = np.interp(x[sel], x, lM) - MB(x[sel])
    print('%-16s' % lab, ' '.join('%+7.3f' % q for q in d), '| over B: mean %+.3f rms %.3f' % (dd.mean(), np.sqrt((dd ** 2).mean())))
print('panel B range: M 10^%.3f .. 10^%.3f Msun, r 10^%.3f .. 10^%.3f pc; C/D/E: 10^%.3f .. 10^%.3f Msun' % (lMB[0], lMB[-1], lrB[0], lrB[-1], lMT[0], lMT[-1]))
r9 = 1e9 * Rsun / pc; print('1e9 R_sun = %.3f pc: M_B = 10^%.4f = %.1f Msun' % (r9, MB(np.log10(r9)), 10 ** MB(np.log10(r9))))
for m in (0.0, 1.0, 2.0, 3.0): print('M = 10^%.0f Msun at r = %.4g pc' % (m, 10 ** np.interp(m, lMB, lrB)))
print('T %.0f..%.0f K (centre %.0f); f_H2 %.3g (centre) .. %.3g; v_r min %.3f km/s at M = 10^%.2f; v_r at 0.01, 1, 100 Msun: %s' %
      (10 ** lT.min(), 10 ** lT.max(), 10 ** lT[0], 10 ** lfH2[0], 10 ** lfH2.min(), v.min(), lMv[np.argmin(v)], [round(float(np.interp(q, lMv, v)), 3) for q in (-2, 0, 2)]))
# IC=abn02 hand read-off (julia/popiii_m2c.jl, the non-z19 abn02 tables)
lr_h = [-3.0, -1.57, -1.28, -1.08, -0.84, -0.59, -0.35, 0.06, 0.87]; ln_h = [5.83, 5.81, 5.40, 5.25, 4.73, 3.99, 3.47, 2.87, 0.97]
lM_h = [-2.0, 1.5, 2.0, 2.5, 3.0, 3.4, 3.7, 4.0]; T_h = [252.0, 252.0, 229.0, 205.0, 164.0, 355.0, 1120.0, 1000.0]
lf_h = [-2.90, -2.90, -2.98, -3.11, -3.18, -3.37, -3.65, -4.0]; v_h = [-0.03, -0.32, -0.49, -0.19, -0.81, -0.70, -1.73, -0.95]
print('IC=abn02 hand tables vs this extraction (hand - vector):')
print('  log n at log r', lr_h, ':', [round(a - float(np.interp(b, lr, ln)), 2) for a, b in zip(ln_h, lr_h)])
print('  log T at log M', lM_h, ':', [round(np.log10(a) - float(np.interp(b, lMT, lT)), 3) for a, b in zip(T_h, lM_h)])
print('  log f_H2      :', [round(a - float(np.interp(b, lMH, lfH2)), 3) for a, b in zip(lf_h, lM_h)])
print('  v_r [km/s]    :', [round(a - float(np.interp(b, lMv, v)), 3) for a, b in zip(v_h, lM_h)])
print('  vector T at those M:', [round(10 ** float(np.interp(b, lMT, lT)), 0) for b in lM_h])
print('  vector v at those M:', [round(float(np.interp(b, lMv, v)), 3) for b in lM_h])
