import sys, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/g2')
import frames as FR

out = sys.argv[1]; png = sys.argv[2]
L = open(out + '/grid.txt').read().split('\n')
F = L[0].split(); LR = np.array(L[1].split(), float); LT = np.array(L[2].split(), float)
A = np.fromfile(out + '/grid.bin').reshape(len(LT), len(LR), len(F))   # Julia column-major [f, i, j] -> C [j, i, f]
g = {f: A[:, :, k].T for k, f in enumerate(F)}                        # [iρ, iT]
g0 = g
if len(sys.argv) > 3:                                                  # the SCVHXE=0 grid (same fields)
    L0 = open(sys.argv[3] + '/grid.txt').read().split('\n'); F0 = L0[0].split()
    A0 = np.fromfile(sys.argv[3] + '/grid.bin').reshape(len(LT), len(LR), len(F0)); g0 = {f: A0[:, :, k].T for k, f in enumerate(F0)}
if len(sys.argv) > 4: TITLE2 = sys.argv[4]
else: TITLE2 = 'network + G2 additions'
UOFF = None

fr = FR.J['frames'][-1]
prof = (FR.get(fr, 'logrho'), FR.get(fr, 'logT'))
cl = (np.array([-6.418, -6.657]), np.array([3.983, 3.525]))


def lr(a, b):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.log10(a / b)


panels = [
    ('p (network LTE / table)', lr(g['p_n'], g['p_t'])),
    ('Δu ρ/p (network − table) [kT per particle], colour ±1', (g['u_n'] - g['u_t']) * 10 ** LR[:, None] / g['p_t']),
    ('Γ1', lr(g['G1_n'], g['G1_t'])), ('∇_ad', lr(g['nab_n'], g['nab_t'])), ('c_P', lr(g['cP_n'], g['cP_t'])),
    ('x_e (masked where both < 1e-9)', np.where(np.maximum(g['xe_n'], g['xe_t']) > 1e-9, lr(g['xe_n'], g['xe_t']), np.nan)),
    ('κ_R (network own continuum / table)', lr(g['kR_n'], g['kR_t'])), ('κ_P = κ_E (network / table)', lr(g['kP_n'], g['kP_t'])),
    ('κ_R (Z0 lines / table)', lr(g['kR_z'], g['kR_t'])),
]
if np.isfinite(g['kR_p']).any():
    panels += [('κ_R (' + TITLE2 + ') / table', lr(g['kR_p'], g['kR_t'])), ('κ_P (' + TITLE2 + ') / table', lr(g['kP_p'], g['kP_t'])),
               ('x_e network LTE / table with SCVHXE reading', np.where(np.maximum(g['xe_n'], g['xe_tfix']) > 1e-9, lr(g['xe_n'], g['xe_tfix']), np.nan))]
if 'xe_w' in g and np.isfinite(g['xe_w']).any():
    panels += [('x_e Saha + HM88 occupation prob. (Z0) / table SCVHXE', np.where(np.maximum(g['xe_w'], g['xe_tfix']) > 1e-9, lr(g['xe_w'], g['xe_tfix']), np.nan)),
               ('x_e table: SCVH as read (SCVHXE=0) / SCVHXE=1', np.where(g['xe_tfix'] > 1e-9, lr(g0['xe_t'], g['xe_tfix']), np.nan))]
n = len(panels); nc = 3; nr = (n + nc - 1) // nc
fig, axs = plt.subplots(nr, nc, figsize=(5.2 * nc, 4.0 * nr), sharex=True, sharey=True, constrained_layout=True)
lev = [np.log10(1.05), np.log10(1.3), np.log10(3)]
levu = [0.05, 0.3, 3]
for ax, (t, Z) in zip(axs.flat, panels):
    Zc = np.clip(Z, -1.5, 1.5)
    im = ax.pcolormesh(LR, LT, Zc.T, cmap='RdBu_r', vmin=-1, vmax=1, shading='nearest')
    ax.contour(LR, LT, np.abs(Z).T, levels=levu if t.startswith('Δu') else lev, colors=['0.55', '0.25', 'k'], linewidths=[0.6, 0.9, 1.3])
    ax.contour(LR, LT, g['w_scvh'].T, levels=[0.01, 0.99], colors='green', linestyles=[':', '--'], linewidths=0.9)
    ax.axhline(np.log10(3e5), color='0.4', lw=0.6, ls='-.')
    ax.plot(prof[0], prof[1], color='orange', lw=1.2)
    ax.plot(cl[0], cl[1], 'o', ms=4, mfc='none', mec='m')
    ax.set_title(t, fontsize=9)
    ax.set_xlim(LR[0], LR[-1]); ax.set_ylim(LT[0], LT[-1])
for ax in axs[-1]: ax.set_xlabel('log ρ [g cm⁻³]')
for ax in axs[:, 0]: ax.set_ylabel('log T [K]')
for ax in axs.flat[n:]: ax.axis('off')
cb = fig.colorbar(im, ax=axs, shrink=0.6, label='log10 ratio (Δu panel: Δu ρ/p, contours 0.05/0.3/3 kT)')
fig.suptitle('G2: network material in LTE vs equilibrium table.  contours |Δ| = 5 % (grey), 30 % (dark), ×3 (black); green: SCVH weight 0 (dotted) / 1 (dashed); '
             'orange: zdiagB step 2650 profile, magenta: cliff zones 61/62; dash-dot: HOT branch 3e5 K', fontsize=8)
fig.savefig(png, dpi=110)
print('wrote', png)
