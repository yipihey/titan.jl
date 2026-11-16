import sys, numpy as np
d = sys.argv[1]; thr = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1
L = open(d + '/gridnd.txt').read().split('\n'); F = L[0].split(); LR = np.array(L[1].split(), float); LT = np.array(L[2].split(), float)
A = np.fromfile(d + '/gridnd.bin').reshape(len(LT), len(LR), len(F)); g = {f: A[:, :, k].T for k, f in enumerate(F)}
dp = np.log(g['p_n'] / g['p_t']); du = (g['u_n'] - g['u_t']) * 10 ** LR[:, None] / g['p_t']
E = np.maximum(np.abs(dp), np.abs(du)); ok = np.isfinite(E)
s = np.maximum(np.maximum(g['Gam'] / thr, g['Deg'] / thr), g['PHM'] / thr); m = ok & (s < 1)
idx = np.argsort(np.where(m, E, -1).ravel())[::-1][:int(sys.argv[3]) if len(sys.argv) > 3 else 20]
for f in idx:
    i, j = np.unravel_index(f, E.shape)
    print(f'lr {LR[i]:5.1f} lT {LT[j]:4.2f}: dlnp {dp[i,j]:+.3f} du {du[i,j]:+.3f} kT  xe n/t {g["xe_n"][i,j]:.3e}/{g["xe_t"][i,j]:.3e} wscvh {g["w_scvh"][i,j]:.2f} G {g["Gam"][i,j]:.2e} D {g["Deg"][i,j]:.2e} P {g["PHM"][i,j]:.2e}')
