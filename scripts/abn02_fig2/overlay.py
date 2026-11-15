# Overlay of the extracted final-output curves (data -> page via the inverse axis maps of calib.json) on a raster render of
# ABN02 page 4.  usage: python3 -I overlay.py <page png> <dpi> <calib.json> <abn02_final_fig2.csv> <out.png>
import sys, json
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

png, dpi, cal, csv, out = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
C = json.load(open(cal)); im = Image.open(png); s = dpi / 72.0
inv = lambda key: (lambda v: (np.asarray(v) - C[key]['offset']) / C[key]['slope_per_pt'])
ax_of = {'n_r': ('x: radius [pc] (A bottom)', 'y: n [cm^-3] (A left)'), 'dm_r': ('x: radius [pc] (A bottom)', 'y: n [cm^-3] (A left)'),
         'M_r': ('x: radius [pc] (B bottom-edge ticks, labels of A)', 'y: M_enc [Msun] (B left)'),
         'fH2_M': ('x: M_enc [Msun] (E bottom, shared C/D/E)', 'y: mass fractions (C right)'),
         'fHI_M': ('x: M_enc [Msun] (E bottom, shared C/D/E)', 'y: mass fractions (C right)'),
         'T_M': ('x: M_enc [Msun] (E bottom, shared C/D/E)', 'y: T [K] (D right)'),
         'vr_M': ('x: M_enc [Msun] (E bottom, shared C/D/E)', 'y: v [km/s] (E right, linear)'),
         'mcs_M': ('x: M_enc [Msun] (E bottom, shared C/D/E)', 'y: v [km/s] (E right, linear)'),
         'mcs_M_symbols': ('x: M_enc [Msun] (E bottom, shared C/D/E)', 'y: v [km/s] (E right, linear)')}
rows = [l.strip().split(',') for l in open(csv) if not l.startswith('#') and not l.startswith('curve')]
x0, x1, y0, y1 = 130, 500, 45, 400                              # page window [pt]
im = im.crop((int(x0 * s), int(y0 * s), int(x1 * s), int(y1 * s)))
fig = plt.figure(figsize=(im.size[0] / 100, im.size[1] / 100), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(im, extent=(x0 * s, x1 * s, y1 * s, y0 * s)); ax.set_axis_off()
col = {'n_r': '#00c8ff', 'dm_r': '#ff8c00', 'M_r': '#00c8ff', 'fH2_M': '#00c8ff', 'fHI_M': '#ff8c00', 'T_M': '#00c8ff', 'vr_M': '#00c8ff',
       'mcs_M': '#00e000', 'mcs_M_symbols': '#ff00ff'}
for k, (kx, ky) in ax_of.items():
    r = np.array([[float(a[1]), float(a[2])] for a in rows if a[0] == k])
    if len(r) == 0: continue                                     # (the red CSV has no DM / -c_s curves)
    px = inv(kx)(r[:, 0]) * s; py = inv(ky)(r[:, 1]) * s
    ax.plot(px, py, 'o', ms=2.2, mfc='none', mec=col[k], mew=0.7, label=k)
ax.set_xlim(130 * s, 500 * s); ax.set_ylim(400 * s, 45 * s)
ax.legend(loc='lower left', fontsize=7, framealpha=0.85)
which = sys.argv[6] if len(sys.argv) > 6 else 'final output'           # optional 6th argument: which output (title only)
ax.set_title('ABN02 Fig. 2 %s: extracted points (open circles) mapped back through the tick calibration' % which, fontsize=8)
fig.savefig(out, dpi=100)
print('wrote', out)
