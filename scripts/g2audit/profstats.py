import sys, numpy as np
sys.path.insert(0, '/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/g2')
import frames as FR


def load(d):
    L = open(d + '/grid.txt').read().split('\n'); F = L[0].split(); LR = np.array(L[1].split(), float); LT = np.array(L[2].split(), float)
    A = np.fromfile(d + '/grid.bin').reshape(len(LT), len(LR), len(F))
    return {f: A[:, :, k].T for k, f in enumerate(F)}, LR, LT


g, LR, LT = load(sys.argv[1])
fr = FR.J['frames'][-1]
r = FR.get(fr, 'logrho'); T = FR.get(fr, 'logT'); tau = FR.get(fr, 'logtau'); m = FR.get(fr, 'logm')
qs = [('p', 'p_n', 'p_t'), ('cP', 'cP_n', 'cP_t'), ('G1', 'G1_n', 'G1_t'), ('xe(XE)', 'xe_n', 'xe_tfix'), ('kR', 'kR_n', 'kR_t'), ('kP', 'kP_n', 'kP_t'), ('kR_p', 'kR_p', 'kR_t'), ('kP_p', 'kP_p', 'kP_t')]
print('zone logrho logT logtau | ' + ' '.join(f'{q[0]:>7s}' for q in qs) + '  (net/table ratio at nearest grid node)')
sel = [k for k in range(fr['n']) if 3.0 <= T[k] <= 5.48 and -14 <= r[k] <= 0]
for k in sel[::1]:
    i = int(round((r[k] - LR[0]) / 0.1)); j = int(round((T[k] - LT[0]) / 0.02))
    vals = [g[a][i, j] / g[b][i, j] for _, a, b in qs]
    if k % 4 == 0 or 54 <= k <= 64:
        print(f'{k:4d} {r[k]:6.2f} {T[k]:5.2f} {tau[k]:6.2f} | ' + ' '.join(f'{v:7.3g}' for v in vals))
# summary along the profile: fraction of zones with |log ratio| > 5 %, 30 %, x3
print('fractions of the', len(sel), 'profile zones with 1e3<T<3e5 K beyond 5% / 30% / x3:')
for nm, a, b in qs:
    v = []
    for k in sel:
        i = int(round((r[k] - LR[0]) / 0.1)); j = int(round((T[k] - LT[0]) / 0.02)); v.append(abs(np.log(g[a][i, j] / g[b][i, j])))
    v = np.array(v); print(f'  {nm:7s} {np.mean(v > np.log(1.05)):.2f} {np.mean(v > np.log(1.3)):.2f} {np.mean(v > np.log(3)):.2f}')
