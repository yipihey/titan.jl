# ABN02 Fig. 2 (Abel, Bryan & Norman 2002, Science 295, 93; arXiv:astro-ph/0112088), page 4 of the arXiv PDF, as vector
# graphics (pdftocairo -svg -f 4 -l 4) -> the FINAL-output curves (black dotted with circles, z = 18.181164).
# usage: python3 -I extract.py <p4.svg> <pdftotext -bbox page-4 html> <outdir>
import sys, os, re, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from svgparse import parse

svg, bbox, outdir = sys.argv[1:4]
P = parse(svg); byline = {p['line']: p for p in P}
W = [(float(a), float(b), float(c), float(d), t.replace('−', '-')) for a, b, c, d, t in
     re.findall(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>', open(bbox).read())]


def labels(region, log=True):
    """tick labels inside region (x0,x1,y0,y1) of the page: '10' + superscript exponent (log) or plain integers (linear)"""
    x0, x1, y0, y1 = region; out = []
    ws = [w for w in W if x0 <= w[0] and w[2] <= x1 and y0 <= w[1] and w[3] <= y1]
    for w in ws:
        if log and w[4] == '10':
            # the superscript: starts where the '10' ends, its top 2-6 pt above the top of the '10'
            e = [v for v in ws if abs(v[0] - w[2]) < 0.6 and 2 < w[1] - v[1] < 6 and re.fullmatch(r'-?\d+', v[4])]
            if len(e) == 1:
                v = e[0]; out.append((int(v[4]), 0.5 * (w[0] + v[2]), 0.5 * (w[1] + w[3])))
        elif not log and re.fullmatch(r'-?\d+', w[4]):
            out.append((int(w[4]), 0.5 * (w[0] + w[2]), 0.5 * (w[1] + w[3])))
    return sorted(out)


def ticks(line, axis, edge, tol=0.3):
    """major-tick positions of tick path `line` along `axis` ('x': vertical ticks starting on the horizontal edge y = edge)"""
    pos = []
    for s in byline[line]['subpaths']:
        (xa, ya), (xb, yb) = s[0], s[-1]
        if axis == 'x' and abs(xa - xb) < 1e-6 and abs(ya - edge) < tol: pos.append(xa)
        if axis == 'y' and abs(ya - yb) < 1e-6 and abs(xa - edge) < tol: pos.append(ya)
    return np.array(sorted(set(round(p, 4) for p in pos)))


report = {}


def calib(name, tickline, axis, edge, labreg, log=True):
    """linear map page coordinate -> data (log10 or linear): least squares through the labelled major ticks"""
    tk = ticks(tickline, axis, edge); lab = labels(labreg, log)
    rows = []
    for v, xc, yc in lab:
        c = xc if axis == 'x' else yc
        j = int(np.argmin(abs(tk - c))); rows.append((v, tk[j], c - tk[j]))
    v = np.array([r[0] for r in rows], float); t = np.array([r[1] for r in rows])
    (s, o), *_ = np.linalg.lstsq(np.vstack([t, np.ones_like(t)]).T, v, rcond=None)
    res = v - (s * t + o)
    allv = s * tk + o; resall = allv - np.round(allv)           # every major tick must fall on an integer (decade / unit)
    report[name] = dict(labels_value_tickpos_labeloffset=[(int(r[0]), round(float(r[1]), 3), round(float(r[2]), 2)) for r in rows],
                        slope_per_pt=float(s), offset=float(o), label_fit_resid=[round(float(x), 5) for x in res],
                        n_major_ticks=int(len(tk)), major_tick_maxdev=float(abs(resall).max()), major_tick_rms=float(np.sqrt((resall ** 2).mean())))
    return lambda c: s * np.asarray(c) + o


# axis frames (page pt, y down): A x 152.07-330.92, y 190.31-384.69; B y 73.66-190.31 (same x); right column x 330.9-488.7,
# C y 73.66-170.92, D 170.92-268.04, E 268.04-384.69.  Tick paths found by listing the solid black stroked paths.
X_r = calib('x: radius [pc] (A bottom)', 852, 'x', 384.69, (100, 335, 383, 402))
Y_A = calib('y: n [cm^-3] (A left)', 905, 'y', 152.07, (100, 152, 197, 392))
X_rB = calib('x: radius [pc] (B bottom-edge ticks, labels of A)', 1870, 'x', 190.31, (100, 335, 383, 402))
Y_B = calib('y: M_enc [Msun] (B left)', 1882, 'y', 152.07, (100, 152, 60, 197))
X_au = calib('x: r [a.u.] (B top)', 2089, 'x', 73.66, (140, 330, 50, 68))
X_M = calib('x: M_enc [Msun] (E bottom, shared C/D/E)', 1433, 'x', 384.69, (335, 490, 383, 402))
Y_C = calib('y: mass fractions (C right)', 1162, 'y', 488.72, (490, 520, 60, 180))
Y_D = calib('y: T [K] (D right)', 1733, 'y', 488.72, (490, 520, 180, 278))
Y_E = calib('y: v [km/s] (E right, linear)', 1485, 'y', 488.55, (490, 515, 278, 392), log=False)


def minor(name, f, tickline, axis, edge):
    """minor ticks of a log axis: their mapped values must sit at log10(k), k = 1..9, within each decade"""
    tk = ticks(tickline, axis, edge); v = f(tk); fr = v - np.floor(v + 1e-9)
    lk = np.log10(np.arange(1, 11)); dev = np.array([np.abs(lk - x).min() for x in fr])
    report[name]['minor_ticks_n'] = int(len(tk)); report[name]['minor_tick_maxdev'] = float(dev.max())
    report[name]['minor_tick_rms'] = float(np.sqrt((dev ** 2).mean()))


minor('x: radius [pc] (A bottom)', X_r, 851, 'x', 384.69)
minor('y: n [cm^-3] (A left)', Y_A, 904, 'y', 152.07)
minor('y: M_enc [Msun] (B left)', Y_B, 1881, 'y', 152.07)
minor('x: M_enc [Msun] (E bottom, shared C/D/E)', X_M, 1432, 'x', 384.69)
minor('y: mass fractions (C right)', Y_C, 1161, 'y', 488.72)
minor('y: T [K] (D right)', Y_D, 1732, 'y', 488.72)


def pts(line):
    return np.array([q for sp in byline[line]['subpaths'] for q in sp])


def style(line):
    st = byline[line]['style']; return dict(stroke=st.get('stroke'), dash=st.get('stroke-dasharray'), width=st.get('stroke-width'))


def markers(line0, color):
    """centres of the circle markers (closed 4-Bezier subpaths) that directly follow path `line0` in document order"""
    c = []
    for p in P:
        if p['line'] <= line0: continue
        if p['ncurve'] == 4 and p['style'].get('stroke') == color:
            q = np.array(p['subpaths'][0]); c.append(((q[:, 0].max() + q[:, 0].min()) / 2, (q[:, 1].max() + q[:, 1].min()) / 2))
        else:
            break
    return np.array(c)


curves = {}


def add(key, line, fx, fy, xname, yname, mcolor='rgb(0%,0%,0%)'):
    q = pts(line); x = fx(q[:, 0]); y = fy(q[:, 1]); mk = markers(line, mcolor)
    def segdist(m):                                   # distance of a marker centre from the polyline [pt]
        a, b = q[:-1], q[1:]; ab = b - a; t = np.clip(((m - a) * ab).sum(1) / np.maximum((ab ** 2).sum(1), 1e-30), 0, 1)
        return float(np.hypot(*(a + t[:, None] * ab - m).T).min())
    d = [segdist(m) for m in mk]
    curves[key] = dict(svg_line=line, style=style(line), x=x, y=y, xname=xname, yname=yname, page=q, nmark=len(mk),
                       mark_vertex_dist=(min(d), max(d)) if d else None)


add('A_n', 1037, X_r, Y_A, 'log10_r_pc', 'log10_n_cm3')
add('A_DM', 1113, X_r, Y_A, 'log10_r_pc', 'log10_rhoDM_GeV_cm3')
add('B_M', 2005, X_rB, Y_B, 'log10_r_pc', 'log10_Menc_Msun')
add('C_a', 1308, X_M, Y_C, 'log10_Menc_Msun', 'log10_massfrac')
add('C_b', 1369, X_M, Y_C, 'log10_Menc_Msun', 'log10_massfrac')
add('D_T', 1806, X_M, Y_D, 'log10_Menc_Msun', 'log10_T_K')
add('E_v', 1579, X_M, Y_E, 'log10_Menc_Msun', 'vr_km_s')
add('E_mcs', 1641, X_M, Y_E, 'log10_Menc_Msun', 'minus_cs_km_s', mcolor='rgb(58.006287%,0%,82.8125%)')
pur = [p for p in P if '58.006287%' in p['style'].get('fill', '')]
pc_ = np.array([np.array(p['subpaths'][0]).mean(0) for p in pur])
curves['E_mcs_symbols'] = dict(x=X_M(pc_[:, 0]), y=Y_E(pc_[:, 1]), page=pc_, xname='log10_Menc_Msun', yname='minus_cs_km_s')
X_au_of_r = None
for k in ('C_a', 'C_b'):
    c = curves[k]; print(k, 'log f at min M', round(float(c['y'][np.argmin(c['x'])]), 3), ' at max M', round(float(c['y'][np.argmax(c['x'])]), 3))
# (the curves dict is not saved; the CSV below carries everything)
json.dump(report, open(os.path.join(outdir, 'calib.json'), 'w'), indent=1)
for k, v in report.items(): print(k, json.dumps(v))
for k, c in curves.items():
    if 'style' in c:
        print(k, c['svg_line'], c['style'], len(c['x']), 'markers', c['nmark'], 'marker-vertex dist', c['mark_vertex_dist'],
              'x %.3f..%.3f y %.3f..%.3f' % (min(c['x']), max(c['x']), min(c['y']), max(c['y'])))
# consistency of the two radius axes: a.u. axis (top of B) vs pc axis: log10(1 pc / 1 au) = 5.31442
s_pc = report['x: radius [pc] (A bottom)']; s_au = report['x: r [a.u.] (B top)']
xs = np.array([200.0, 300.0])
print('log10(au/pc) from the two axes:', (X_au(xs) - X_r(xs)).round(4), '(expected 5.3144)')
# ---- CSV (long format): one row per path vertex
names = {'A_n': 'n_r', 'A_DM': 'dm_r', 'B_M': 'M_r', 'C_a': 'fH2_M', 'C_b': 'fHI_M', 'D_T': 'T_M', 'E_v': 'vr_M', 'E_mcs': 'mcs_M',
         'E_mcs_symbols': 'mcs_M_symbols'}
with open(os.path.join(outdir, 'abn02_final_fig2.csv'), 'w') as f:
    f.write('# ABN02 (Science 295, 93; astro-ph/0112088) Fig. 2, final output (z = 18.181164, black dotted with circles): vertices of the\n'
            '# vector paths of the arXiv PDF page 4 (pdftocairo -svg), mapped to data by the axis tick marks (README.md, scripts/abn02_fig2).\n'
            '# curve: n_r = log10 n[cm^-3] vs log10 r[pc] (A); dm_r = log10 rho_DM[GeV cm^-3] vs log10 r[pc] (A, final-time DM line);\n'
            '# M_r = log10 M_enc[Msun] vs log10 r[pc] (B); fH2_M, fHI_M = log10 mass fraction vs log10 M_enc[Msun] (C; the second curve is\n'
            '# atomic H, not H+: fH2 + fHI = 0.756 +- 0.005 = X at every mass); T_M = log10 T[K] vs log10 M (D); vr_M = v_r[km/s] vs log10 M (E);\n'
            '# mcs_M = -c_s[km/s] vs log10 M (E, purple dotted path, clipped by the frame at -3); mcs_M_symbols = the filled purple symbols.\n'
            'curve,x,y,page_x_pt,page_y_pt\n')
    for k, nm in names.items():
        c = curves[k]
        for x, y, (px, py) in zip(c['x'], c['y'], c['page']):
            f.write('%s,%.5f,%.5f,%.4f,%.4f\n' % (nm, x, y, px, py))
print('wrote', os.path.join(outdir, 'abn02_final_fig2.csv'))
# ---- optional: the SECOND output ("nine Myrs later", RED dotted with red circles: dasharray 6.66,13.32, stroke rgb(100%,0,0))
#      -> <outdir_red>/abn02_red_fig2.csv (same format; the figure has no DM line and no -c_s at this time)
if len(sys.argv) > 4:
    outred = sys.argv[4]; red = 'rgb(100%,0%,0%)'
    reddot = sorted(p['line'] for p in P if p['style'].get('stroke') == red and p['style'].get('stroke-dasharray') == '6.65974,13.3195')
    print('red dotted paths (SVG line numbers):', reddot)
    assert reddot == [981, 1210, 1254, 1529, 1756, 1949], reddot
    add('R_n', 981, X_r, Y_A, 'log10_r_pc', 'log10_n_cm3', mcolor=red)
    add('R_M', 1949, X_rB, Y_B, 'log10_r_pc', 'log10_Menc_Msun', mcolor=red)
    add('R_a', 1210, X_M, Y_C, 'log10_Menc_Msun', 'log10_massfrac', mcolor=red)
    add('R_b', 1254, X_M, Y_C, 'log10_Menc_Msun', 'log10_massfrac', mcolor=red)
    add('R_T', 1756, X_M, Y_D, 'log10_Menc_Msun', 'log10_T_K', mcolor=red)
    add('R_v', 1529, X_M, Y_E, 'log10_Menc_Msun', 'vr_km_s', mcolor=red)
    rnames = {'R_n': 'n_r', 'R_M': 'M_r', 'R_a': 'fH2_M', 'R_b': 'fHI_M', 'R_T': 'T_M', 'R_v': 'vr_M'}
    for k in rnames:
        c = curves[k]
        print(k, c['svg_line'], c['style'], len(c['x']), 'markers', c['nmark'], 'marker-polyline dist', c['mark_vertex_dist'],
              'x %.3f..%.3f y %.3f..%.3f' % (min(c['x']), max(c['x']), min(c['y']), max(c['y'])))
    with open(os.path.join(outred, 'abn02_red_fig2.csv'), 'w') as f:
        f.write('# ABN02 (Science 295, 93; astro-ph/0112088) Fig. 2, SECOND output ("nine Myrs later", RED dotted with circles): vertices of the\n'
                '# vector paths of the arXiv PDF page 4 (pdftocairo -svg), mapped to data by the axis tick marks (lib/popiii/abn02fin/calib.json;\n'
                '# lib/popiii/abn02red/README.md, scripts/abn02_fig2). curve: n_r = log10 n[cm^-3] vs log10 r[pc] (A); M_r = log10 M_enc[Msun]\n'
                '# vs log10 r[pc] (B); fH2_M, fHI_M = log10 mass fraction vs log10 M_enc[Msun] (C: H2 and neutral atomic H); T_M = log10 T[K]\n'
                '# vs log10 M (D); vr_M = v_r[km/s] vs log10 M (E). The figure has no DM line and no -c_s at this time.\n'
                'curve,x,y,page_x_pt,page_y_pt\n')
        for k, nm in rnames.items():
            c = curves[k]
            for x, y, (px, py) in zip(c['x'], c['y'], c['page']):
                f.write('%s,%.5f,%.5f,%.4f,%.4f\n' % (nm, x, y, px, py))
    print('wrote', os.path.join(outred, 'abn02_red_fig2.csv'))
