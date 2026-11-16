import sys, numpy as np
d = sys.argv[1]
L = open(d + '/gridnd.txt').read().split('\n'); F = L[0].split(); LR = np.array(L[1].split(), float); LT = np.array(L[2].split(), float)
A = np.fromfile(d + '/gridnd.bin').reshape(len(LT), len(LR), len(F)); g = {f: A[:, :, k].T for k, f in enumerate(F)}
ep = np.abs(np.log(g['p_n'] / g['p_t'])); eu = np.abs(g['u_n'] - g['u_t']) * 10 ** LR[:, None] / g['p_t']
ec = np.abs(np.log(g['cP_n'] / g['cP_t']))
E = np.maximum(ep, eu)
Gam, Deg, PHM = g['Gam'], g['Deg'], g['PHM']
ok = np.isfinite(E)
print('points', ok.sum())
for nm, P in (('Gamma_e', Gam), ('D=n_e lam^3/2', Deg), ('P_HM=-ln w1', PHM)):
    print(nm)
    for lo, hi in ((1e-6, 1e-3), (1e-3, 1e-2), (1e-2, 3e-2), (3e-2, 0.1), (0.1, 0.3), (0.3, 1), (1, 3), (3, 100)):
        m = ok & (P >= lo) & (P < hi)
        if m.sum() < 5: continue
        print(f'  [{lo:.0e},{hi:.0e}) n={m.sum():5d}  max(|dlnp|,|du| kT) p50 {np.percentile(E[m],50):.3f} p90 {np.percentile(E[m],90):.3f} max {E[m].max():.3f} | dlncP p90 {np.percentile(ec[m],90):.3f}')
# combined criterion s = max(Gam/G0, Deg/D0, PHM/P0)
for G0, D0, P0 in ((0.1, 0.1, 0.1), (0.2, 0.1, 0.1), (0.3, 0.3, 0.3), (0.1, 0.05, 0.05), (0.05, 0.05, 0.05)):
    s = np.maximum(np.maximum(Gam / G0, Deg / D0), PHM / P0)
    m = ok & (s < 1)
    bad = m & (E > 0.05)
    print(f'G0 {G0} D0 {D0} P0 {P0}: s<1 points {m.sum()}, E p99 {np.percentile(E[m],99):.3f} max {E[m].max():.3f}; points with E>0.05: {bad.sum()}',
          ' at (lr,lT):', [(round(LR[i],1), round(LT[j],2), round(E[i,j],3)) for i, j in zip(*np.where(bad))][:6])
