// qstar.js - quasi-star runs of the Julia TITAN, replayed in time against the steady envelope they started from.
import { loadData, h, caption, fmt } from '../core.js';
import { LinePlot } from '../viz/plot.js';
import { select, slider, readout } from '../viz/ui.js';

const BOOK = new URL('../../../chapters/07-titan-quasistar.html', import.meta.url).href;
const RUNS = { QSTAR: 'hydrostatic test, no infall', E1: 'infall, Eulerian grid, ramped (fails)', A1: 'infall, adaptive grid (fails)' };
const VARS = {
  T: ['T [K]', true], rho: ['ρ [g cm⁻³]', true], u: ['u [cm s⁻¹]', false],
  Lconv: ['L_conv / L', false], Lrad: ['L_rad / L', false],
};
const centres = (r) => Float64Array.from({ length: r.length - 1 }, (_, i) => 0.5 * (r[i] + r[i + 1]));

export async function mountQstar(host, S, onChange = () => {}, opts = {}) {
  if (opts.run && RUNS[opts.run]) S.run = opts.run;
  const bar = h('div', 'qd-controls', host);
  const info = readout(host, { label: 'snapshot' });
  const box = h('div', 'qd-panel', host);
  const d = h('div', null, box); d.setAttribute('aria-label', 'quasi-star TITAN run profile versus radius');
  const plot = new LinePlot(d, { width: 760, height: 360 });
  caption(host, 'Solid: TITAN (Julia) at the selected time; dashed: the steady envelope it started from. Radius in units of the steady photospheric radius R. Stored run output, not recomputed here.', BOOK + '#sec-hydrostatic');
  let D = null;
  const draw = () => {
    if (!D) return;
    const k = Math.round(S.ti * (D.snaps.length - 1));
    const s = D.snaps[k];
    const v = VARS[S.qvar] ? S.qvar : 'T';
    const [ylabel, ylog] = VARS[v];
    const atInterfaces = v === 'u' || v === 'Lconv' || v === 'Lrad';
    const x = (atInterfaces ? s.r : centres(s.r)).map((r) => r / D.R);
    const series = [{ x, y: s[v], label: `TITAN, t = ${fmt(s.t, 3)} s` }];
    const st = D.steady;
    if (v === 'T' || v === 'rho' || v === 'Lconv') {
      const sx = st.r.map((r) => r / D.R), sy = v === 'T' ? st.T : v === 'rho' ? st.rho : st.Lconv;
      series.push({ x: sx, y: sy, label: 'steady model', dash: '5 3', color: '#D55E00' });
    }
    plot.setAxes({ xlabel: 'r / R', ylabel, ylog, title: RUNS[S.run], xlim: [x[0], x[x.length - 1]] });
    plot.setSeries(series);
    plot.draw();
    info.set(`${k + 1} of ${D.snaps.length}: t = ${fmt(s.t, 4)} s (${fmt(s.t / 6.7e7, 3)} t_dyn); last stored snapshot t = ${fmt(D.snaps[D.snaps.length - 1].t, 4)} s`);
    onChange(S);
  };
  const load = async () => { D = await loadData(`qstar_${S.run}.json`); plot.resetView(); draw(); };
  select(bar, { label: 'run', options: Object.keys(RUNS).map((k) => ({ value: k, label: RUNS[k] })), value: S.run, onChange: (v) => { S.run = v; load(); } });
  select(bar, { label: 'variable', options: Object.keys(VARS).map((k) => ({ value: k, label: VARS[k][0] })), value: S.qvar, onChange: (v) => { S.qvar = v; plot.resetView(); draw(); } });
  slider(bar, { label: 'time', min: 0, max: 1, step: 0.001, value: S.ti, format: (v) => `${Math.round(v * 100)} %`, onInput: (v) => { S.ti = v; draw(); } });
  await load();
  return { draw };
}
