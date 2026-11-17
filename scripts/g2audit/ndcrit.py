import sys, numpy as np
d = sys.argv[1]; G0, D0, P0 = (float(x) for x in sys.argv[2:5])
L = open(d + '/gridnd.txt').read().split('\n'); F = L[0].split(); LR = np.array(L[1].split(), float); LT = np.array(L[2].split(), float)
A = np.fromfile(d + '/gridnd.bin').reshape(len(LT), len(LR), len(F)); g = {f: A[:, :, k].T for k, f in enumerate(F)}
dp = np.log(g['p_n'] / g['p_t']); du = (g['u_n'] - g['u_t']) * 10 ** LR[:, None] / g['p_t']
dc = np.log(g['cP_n'] / g['cP_t']); dg = np.log(g['G1_n'] / g['G1_t'])
s = ((g['Gam'] / G0) ** 4 + (g['Deg'] / D0) ** 4 + (g['Pneu'] / P0) ** 4) ** 0.25
ok = np.isfinite(dp) & np.isfinite(s)
for lab, m in (('s < 1/3 (network)', s < 1/3), ('1/3 <= s < 1 (blend)', (s >= 1/3) & (s < 1)), ('1 <= s < 3', (s >= 1) & (s < 3)), ('s >= 3', s >= 3)):
    m = m & ok
    q = lambda a: (np.nanpercentile(np.abs(a[m]), 50), np.nanpercentile(np.abs(a[m]), 95), np.nanmax(np.abs(a[m])))
    print(f'{lab:20s} n={m.sum():5d}  |dln p| p50/p95/max {"/".join(f"{v:.3f}" for v in q(dp))}   |du| kT {"/".join(f"{v:.3f}" for v in q(du))}   |dln cP| {"/".join(f"{v:.3f}" for v in q(dc))}   |dln G1| {"/".join(f"{v:.3f}" for v in q(dg))}')
np.save(d + '/s.npy', s)
