// qstar.js - quasi-star runs of the Julia TITAN, replayed in time against the steady envelope they started from.
import { loadData, h, caption, fmt } from '../core.js';
import { LinePlot } from '../viz/plot.js';
import { select, slider, readout } from '../viz/ui.js';

const BOOK = new URL('../../../chapters/07-titan-quasistar.html', import.meta.url).href;
const RUNS = {
  R2: 'infall to 3 M⊙/yr, Eulerian grid (R2)',
  S6: 'as R2, with the inner mass sink (S6)',
  A220: 'adaptive grid, 220 zones, to 6e8 s (A220)',
  A440: 'adaptive grid, 440 zones, to 1.2e8 s (A440)',
  QSTAR: 'hydrostatic test, no infall',
  E1: 'infall before the fixes (fails at 2.3e4 s)',
  A1: 'infall, adaptive, before the fixes (fails at 8 s)',
};
const VARS = {
  T: ['T [K] (gas, and radiation T_r = (E/a)^¼)', true], rho: ['ρ [g cm⁻³]', true], u: ['u [cm s⁻¹]', false],
  Lrad: ['L_rad / L', false], Lconv: ['L_conv / L', false], dr: ['zone width Δr / R', true],
};
const HIST = {
  radii: ['r / R', false, [['rsh', 'shock'], ['rph', 'photosphere (τ_R = 2/3)']]],
  mdot: ['Ṁ [M⊙/yr]', true, [['Mdot_sh', 'just outside the shock'], ['Mdot_out', 'at the outer boundary']]],
  teff: ['T_eff [K]', false, [['Teff', 'T_eff at the photosphere']]],
  tpk: ['T [K]', true, [['Tpk', 'peak gas T behind the shock']]],
  lum: ['L / L', false, [['Lph', 'at the photosphere']]],
};
const ARAD = 7.5657e-15;
const COLORS = ['#0072B2', '#D55E00', '#009E73'];
const centres = (r) => Float64Array.from({ length: r.length - 1 }, (_, i) => 0.5 * (r[i] + r[i + 1]));

export async function mountQstar(host, S, onChange = () => {}, opts = {}) {
  if (opts.run && RUNS[opts.run]) S.run = opts.run;
  if (!RUNS[S.run]) S.run = 'R2';
  if (!HIST[S.qhist]) S.qhist = 'radii';
  const bar = h('div', 'qd-controls', host);
  const info = readout(host, { label: 'snapshot' });
  const box = h('div', 'qd-panel', host);
  const d = h('div', null, box); d.setAttribute('aria-label', 'quasi-star TITAN run profile versus radius');
  const plot = new LinePlot(d, { width: 760, height: 360 });
  caption(host, 'Solid: TITAN (Julia) at the selected time; dashed: the steady envelope it started from. Radius in units of the steady photospheric radius R. Stored run output, not recomputed here.', BOOK + '#sec-titan-infall');
  const hbar = h('div', 'qd-controls', host);
  const hbox = h('div', 'qd-panel', host);
  const hd = h('div', null, hbox); hd.setAttribute('aria-label', 'time history of the infall run');
  const hplot = new LinePlot(hd, { width: 760, height: 260 });
  const hcap = caption(host, 'Time history of the run (every stored step). The vertical line marks the snapshot shown above. Infall runs only.', BOOK + '#sec-titan-infall');
  let D = null;
  const draw = () => {
    if (!D) return;
    const k = Math.round(S.ti * (D.snaps.length - 1));
    const s = D.snaps[k];
    const v = VARS[S.qvar] ? S.qvar : 'T';
    const [ylabel, ylog] = VARS[v];
    const atInterfaces = v === 'u' || v === 'Lconv' || v === 'Lrad';
    const x = (atInterfaces ? s.r : centres(s.r)).map((r) => r / D.R);
    const y = v === 'dr' ? Float64Array.from({ length: s.r.length - 1 }, (_, i) => (s.r[i + 1] - s.r[i]) / D.R) : s[v];
    const series = [{ x, y, label: v === 'dr' ? `zone widths, t = ${fmt(s.t, 3)} s` : `TITAN gas, t = ${fmt(s.t, 3)} s` }];
    if (v === 'T' && s.Er) series.push({ x, y: s.Er.map((E) => Math.pow(E / ARAD, 0.25)), label: 'TITAN radiation T_r', color: '#009E73' });
    const st = D.steady;
    if (v === 'T' || v === 'rho' || v === 'Lconv') {
      const sx = st.r.map((r) => r / D.R), sy = v === 'T' ? st.T : v === 'rho' ? st.rho : st.Lconv;
      series.push({ x: sx, y: sy, label: 'steady model', dash: '5 3', color: '#D55E00' });
    }
    plot.setAxes({ xlabel: 'r / R', ylabel, ylog, title: RUNS[S.run], xlim: [x[0], x[x.length - 1]] });
    plot.setSeries(series);
    plot.draw();
    info.set(`${k + 1} of ${D.snaps.length}: t = ${fmt(s.t, 4)} s (${fmt(s.t / 6.7e7, 3)} t_dyn); last stored snapshot t = ${fmt(D.snaps[D.snaps.length - 1].t, 4)} s`);
    const H = D.hist;
    hbox.style.display = hbar.style.display = hcap.style.display = H ? '' : 'none';
    if (H) {
      const [hy, hlog, cols] = HIST[S.qhist];
      hplot.setAxes({ xlabel: 't [s]', ylabel: hy, ylog: hlog, title: '' });
      hplot.setSeries(cols.map(([c, lab], i) => ({ x: H.t, y: H[c], label: lab, color: COLORS[i] })));
      hplot.setMarkers([{ x: s.t, label: `t = ${fmt(s.t, 3)} s`, color: '#666' }]);
      hplot.draw();
    }
    onChange(S);
  };
  const load = async () => { D = await loadData(`qstar_${S.run}.json`); plot.resetView(); hplot.resetView(); draw(); };
  select(bar, { label: 'run', options: Object.keys(RUNS).map((k) => ({ value: k, label: RUNS[k] })), value: S.run, onChange: (v) => { S.run = v; load(); } });
  select(bar, { label: 'variable', options: Object.keys(VARS).map((k) => ({ value: k, label: VARS[k][0] })), value: S.qvar, onChange: (v) => { S.qvar = v; plot.resetView(); draw(); } });
  slider(bar, { label: 'time', min: 0, max: 1, step: 0.001, value: S.ti, format: (v) => `${Math.round(v * 100)} %`, onInput: (v) => { S.ti = v; draw(); } });
  select(hbar, { label: 'history', options: Object.keys(HIST).map((k) => ({ value: k, label: HIST[k][2].map((c) => c[1]).join(', ') + ` (${HIST[k][0]})` })), value: S.qhist, onChange: (v) => { S.qhist = v; hplot.resetView(); draw(); } });
  await load();
  return { draw };
}
