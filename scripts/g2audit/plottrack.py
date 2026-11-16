import sys, re, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/tmp/claude-1002/-home-tabel-Projects/3ebf2758-521a-5bd2-af15-1951bd5c13b4/scratchpad/g2')
from cmptrack import load
out = sys.argv[1]; runs = sys.argv[2:]
fig, ax = plt.subplots(4, 1, figsize=(8, 9), sharex=True, constrained_layout=True)
cols = ['#2a6fdb', '#d9822b', '#3a9e5f', '#8e44ad', '#7f7f7f']
for c, spec in zip(cols, runs):
    p, lab = spec.split('=', 1)
    r = load(p); st = [x['step'] for x in r]
    ax[0].semilogy(st, [x['dt'] for x in r], color=c, lw=1.2, label=lab)
    ax[1].semilogy(st, [x['L'] for x in r], color=c, lw=1.0)
    ax[2].plot(st, [x['M'] for x in r], color=c, lw=1.2)
    ax[3].plot(st, [1 if x['on'] else 0 for x in r], color=c, lw=0, marker='o', ms=2)
ax[0].axhline(1e4, color='0.6', lw=0.6, ls='--'); ax[0].axvline(3032, color='0.6', lw=0.6, ls=':')
ax[0].set_ylabel('Δt [s]'); ax[1].set_ylabel('L_ph [L☉]'); ax[2].set_ylabel('M* [M☉]'); ax[3].set_ylabel('cliff on\ntable|network')
ax[3].set_yticks([0, 1]); ax[3].set_yticklabels(['no', 'yes']); ax[3].set_xlabel('step')
ax[0].legend(fontsize=8, loc='lower right'); ax[0].set_title('zdiagB restarts from S5 (step 3000); dotted: zdiagB stopped at 3032', fontsize=9)
fig.savefig(out, dpi=100); print('wrote', out)
