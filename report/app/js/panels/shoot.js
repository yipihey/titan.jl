// shoot.js - the shooting function across trial radii: steady stars are its crossings of 1.
import { loadCore, h, caption, fmt, MODEL, yieldToPaint } from '../core.js';
import { LinePlot } from '../viz/plot.js';
import { select, readout } from '../viz/ui.js';

const BOOK = new URL('../../../chapters/04-infall.html', import.meta.url).href;
const MDOTS = [0, 0.1, 0.3, 1, 3, 5, 7, 8, 9, 9.5, 10, 20, 25];

export async function mountShoot(host, S, onChange = () => {}) {
  const qs = await loadCore();
  const bar = h('div', 'qd-controls', host);
  const status = readout(host, { label: 'roots' });
  const box = h('div', 'qd-panel', host);
  const d = h('div', null, box); d.setAttribute('aria-label', 'enclosed mass at the inner radius over black-hole mass versus trial radius');
  const plot = new LinePlot(d, { width: 760, height: 340 });
  caption(host, 'For each trial radius R the envelope is integrated inward from the infall boundary; the curve is m(ε_in R)/M_BH (1 where the mass runs out first). Solutions are where it crosses 1 (vertical lines).', BOOK + '#sec-branches');
  const draw = async () => {
    status.set('scanning…'); await yieldToPaint();
    const t0 = performance.now();
    const sc = qs.shoot_scan(MODEL.M, S.fbh, S.mdot, S.alpha, S.eps, MODEL.logRmin, MODEL.logRmax, MODEL.dlogR);
    const n = sc.length / 2;
    const lR = Float64Array.from({ length: n }, (_, i) => MODEL.logRmin + i * MODEL.dlogR);
    const ratio = Float64Array.from({ length: n }, (_, i) => sc[2 * i + 1]);
    const roots = qs.solve_all(MODEL.M, S.fbh, S.mdot, S.alpha, S.eps, MODEL.logRmin, MODEL.logRmax, MODEL.dlogR);
    const rs = [];
    for (let i = 0; i + 10 <= roots.length; i += 10) rs.push(roots.slice(i, i + 10));
    plot.setAxes({ ylog: true, xlabel: 'log₁₀ R [cm]', ylabel: 'm(ε_in R) / M_BH', title: `Ṁ = ${S.mdot} M☉/yr, M_BH/M = ${S.fbh}` });
    plot.setSeries([{ x: lR, y: ratio, label: 'shooting function', points: true, radius: 2 },
      { x: lR, y: Float64Array.from(lR, () => 1), label: 'solution: = 1', dash: '4 3', color: '#999' }]);
    plot.setMarkers(rs.map((r, k) => ({ x: Math.log10(r[0]), label: `${k + 1}: T_b ${fmt(r[1], 3)} K` })));
    plot.draw();
    status.set(rs.length ? rs.map((r, k) => `${k + 1}: R = ${fmt(r[0])} cm, T_b = ${fmt(r[1], 4)} K, T_obs = ${fmt(r[7], 4)} K`).join('  ·  ') + `  (${(performance.now() - t0).toFixed(0)} ms)`
      : `no steady solution in 10¹⁴–4×10¹⁶ cm (${(performance.now() - t0).toFixed(0)} ms)`);
    onChange(S);
  };
  select(bar, { label: 'Ṁ [M☉/yr]', options: MDOTS.map((v) => ({ value: v, label: String(v) })), value: String(S.mdot), onChange: (v) => { S.mdot = +v; draw(); } });
  select(bar, { label: 'M_BH/M', options: [{ value: 0.1, label: '0.1' }, { value: 0.01, label: '0.01' }], value: String(S.fbh), onChange: (v) => { S.fbh = +v; draw(); } });
  await draw();
  return { draw };
}
