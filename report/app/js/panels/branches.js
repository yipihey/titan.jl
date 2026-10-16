// branches.js - all steady solutions R(Ṁ): the published table (Julia) and, on request, a live recomputation.
import { loadCore, loadData, h, caption, MODEL, yieldToPaint } from '../core.js';
import { LinePlot } from '../viz/plot.js';
import { select, button, readout } from '../viz/ui.js';

const BOOK = new URL('../../../chapters/04-infall.html', import.meta.url).href;
const GRID = [0, 0.1, 0.3, 0.6, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 9.5, 10, 12, 15, 20, 25];

export async function mountBranches(host, S, onChange = () => {}) {
  const qs = await loadCore();
  const pub = await loadData('branches.json');
  const bar = h('div', 'qd-controls', host);
  const status = readout(host, { label: 'status' });
  const grid = h('div', 'qd-grid2', host);
  const mk = (title, aria) => { const p = h('div', 'qd-panel', grid); h('div', 'qd-title', p, title); const d = h('div', null, p); d.setAttribute('aria-label', aria); return new LinePlot(d, { width: 480, height: 320 }); };
  const pR = mk('Boundary radius of every steady solution', 'boundary radius versus accretion rate');
  const pT = mk('Boundary and observed temperature', 'boundary and observed temperature versus accretion rate');
  caption(host, 'Filled: computed here in the browser; open circles: the Julia table of the book. Ṁ = 0 is drawn at 0.03 on the log axis.', BOOK + '#sec-branches');
  let live = [];
  const mdx = (m) => (m > 0 ? m : 0.03);
  const draw = () => {
    const P = pub.filter((r) => Math.abs(r.fBH - S.fbh) < 1e-9);
    const L = live;
    pR.setAxes({ xlog: true, ylog: true, xlabel: 'Ṁ [M☉/yr]', ylabel: 'R [cm]', xlim: [0.02, 40], ylim: [5e13, 3e16] });
    pR.setSeries([
      { x: L.map((r) => mdx(r.mdot)), y: L.map((r) => r.R), label: 'browser', line: false, points: true },
      { x: P.map((r) => mdx(r.Mdot)), y: P.map((r) => r.R), label: 'Julia table', line: false, points: true, hollow: true, radius: 4 },
    ]);
    pT.setAxes({ xlog: true, ylog: true, xlabel: 'Ṁ [M☉/yr]', ylabel: 'T [K]', xlim: [0.02, 40], ylim: [4000, 8e4] });
    pT.setSeries([
      { x: L.map((r) => mdx(r.mdot)), y: L.map((r) => r.Tb), label: 'T_b (boundary)', line: false, points: true },
      { x: L.map((r) => mdx(r.mdot)), y: L.map((r) => r.Tobs), label: 'T_obs (photosphere)', line: false, points: true },
      { x: P.map((r) => mdx(r.Mdot)), y: P.map((r) => r.Tb), label: 'Julia T_b', line: false, points: true, hollow: true, radius: 4 },
      { x: P.map((r) => mdx(r.Mdot)), y: P.map((r) => r.T_obs), label: 'Julia T_obs', line: false, points: true, hollow: true, radius: 4 },
    ]);
    pR.draw(); pT.draw();
    onChange(S);
  };
  const compute = async () => {
    live = [];
    const t0 = performance.now();
    for (const m of GRID) {
      status.set(`solving Ṁ = ${m} …`); await yieldToPaint();
      const o = qs.solve_all(MODEL.M, S.fbh, m, S.alpha, S.eps, MODEL.logRmin, MODEL.logRmax, MODEL.dlogR);
      for (let i = 0; i + 10 <= o.length; i += 10) live.push({ mdot: m, R: o[i], Tb: o[i + 1], Tobs: o[i + 7] });
      draw();
    }
    status.set(`${live.length} solutions on ${GRID.length} rates in ${((performance.now() - t0) / 1000).toFixed(1)} s`);
  };
  select(bar, { label: 'M_BH/M', options: [{ value: 0.1, label: '0.1' }, { value: 0.01, label: '0.01' }], value: String(S.fbh), onChange: (v) => { S.fbh = +v; live = []; draw(); } });
  button(bar, { label: 'compute all branches in the browser', onClick: () => compute() });
  status.set('published solutions shown; press compute to recompute them here (takes a few seconds)');
  draw();
  return { draw };
}
