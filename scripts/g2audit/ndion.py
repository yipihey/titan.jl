import sys, numpy as np
d = sys.argv[1]
L = open(d + '/gridnd.txt').read().split('\n'); F = L[0].split(); LR = np.array(L[1].split(), float); LT = np.array(L[2].split(), float)
A = np.fromfile(d + '/gridnd.bin').reshape(len(LT), len(LR), len(F)); g = {f: A[:, :, k].T for k, f in enumerate(F)}
dp = np.log(g['p_n'] / g['p_t']); du = (g['u_n'] - g['u_t']) * 10 ** LR[:, None] / g['p_t']
TT = np.broadcast_to(LT[None, :], dp.shape)
for lab, m0 in (('ionised (x_e,table > 0.9)', g['xe_t'] > 0.9), ('T > 2e4 K', TT > np.log10(2e4)), ('molecular/atomic T < 6000 K', TT < np.log10(6000))):
    print(lab)
    for nm, P in (('Gamma_e', g['Gam']), ('P_neu', g['Pneu']), ('D', g['Deg'])):
        for lo, hi in ((1e-3, 1e-2), (1e-2, 0.03), (0.03, 0.1), (0.1, 0.3), (0.3, 1), (1, 3), (3, 10)):
            m = m0 & np.isfinite(dp) & (P >= lo) & (P < hi)
            if m.sum() < 5: continue
            print(f'  {nm:7s} [{lo:g},{hi:g}) n={m.sum():5d} |dlnp| p50 {np.percentile(np.abs(dp[m]),50):.3f} p90 {np.percentile(np.abs(dp[m]),90):.3f} | |du| kT p50 {np.percentile(np.abs(du[m]),50):.3f} p90 {np.percentile(np.abs(du[m]),90):.3f}')
