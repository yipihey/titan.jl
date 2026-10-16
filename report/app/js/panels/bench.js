// bench.js - TITAN benchmarks: the pure-Julia run replayed in time, with the 1995 Fortran code's final state.
import { loadData, h, caption, fmt } from '../core.js';
import { LinePlot } from '../viz/plot.js';
import { select, slider, readout } from '../viz/ui.js';

const BOOK = new URL('../../../chapters/06-titan-port.html', import.meta.url).href;
const VARS = { rho: ['ρ [g cm⁻³]', true], T: ['T_gas [K]', true], Trad: ['T_rad [K]', true], u: ['u [cm s⁻¹]', false] };
const centres = (r) => Float64Array.from({ length: r.length - 1 }, (_, i) => 0.5 * (r[i] + r[i + 1]));

export async function mountBench(host, S, onChange = () => {}) {
  const bar = h('div', 'qd-controls', host);
  const info = readout(host, { label: 'time' });
  const box = h('div', 'qd-panel', host);
  const d = h('div', null, box); d.setAttribute('aria-label', 'benchmark profile versus radius');
  const plot = new LinePlot(d, { width: 760, height: 360 });
  caption(host, 'Solid: the Julia port at the selected time. Dashed: the original Fortran (TITAN 1.2) final dump. Neither is recomputed here.', BOOK + '#sec-port-method');
  let D = null;
  const draw = () => {
    if (!D) return;
    const k = Math.round(S.ti * (D.series.length - 1));
    const s = D.series[k];
    const [ylabel, ylog] = VARS[S.bvar] || VARS.T;
    const zx = (q) => (S.bvar === 'u' ? q.r : centres(q.r));
    const zy = (q) => q[S.bvar];
    plot.setAxes({ xlabel: 'r [cm]', ylabel, ylog: ylog && S.bvar !== 'u', title: `${D.name}: Julia at t = ${fmt(s.t, 4)} s` });
    plot.setSeries([
      { x: zx(s), y: zy(s), label: `Julia, t = ${fmt(s.t, 4)} s` },
      { x: zx(D.fortran), y: zy(D.fortran), label: `Fortran final, t = ${fmt(D.fortran.t, 4)} s`, dash: '5 3', color: '#D55E00' },
    ]);
    plot.draw();
    info.set(`snapshot ${k + 1} of ${D.series.length}; Julia run: ${D.nsteps_julia} steps to t = ${fmt(D.julia.t, 5)} s`);
    onChange(S);
  };
  const load = async () => { D = await loadData(`bench_${S.bench}.json`); draw(); };
  select(bar, { label: 'problem', options: [{ value: 'SOD', label: 'SOD shock tube' }, { value: 'SUBC', label: 'subcritical radiating shock' }, { value: 'SUPC', label: 'supercritical radiating shock' }], value: S.bench, onChange: (v) => { S.bench = v; load(); } });
  select(bar, { label: 'variable', options: Object.keys(VARS).map((k) => ({ value: k, label: VARS[k][0] })), value: S.bvar, onChange: (v) => { S.bvar = v; plot.resetView(); draw(); } });
  slider(bar, { label: 'time', min: 0, max: 1, step: 0.001, value: S.ti, format: (v) => `${Math.round(v * 100)} %`, onInput: (v) => { S.ti = v; draw(); } });
  await load();
  return { draw };
}
