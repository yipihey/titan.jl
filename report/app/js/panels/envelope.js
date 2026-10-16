// envelope.js - a steady quasi-star envelope computed in the browser: all roots of the shooting problem, then the
// profile of the chosen root (temperature, gradients, local Eddington factor, convective Mach number).
import { loadCore, h, caption, fmt, MODEL, yieldToPaint } from '../core.js';
import { LinePlot } from '../viz/plot.js';
import { slider, select, readout } from '../viz/ui.js';

const BOOK = new URL('../../../chapters/03-envelope.html', import.meta.url).href;
const NF = 15; // fields per profile row: r, m, P, T, ρ, ∇, ∇ad, ∇rad, v, cs, κ, Γ, β, H_P, Γ1

export async function mountEnvelope(host, S, onChange = () => {}) {
  const qs = await loadCore();
  const bar = h('div', 'qd-controls', host);
  const status = readout(host, { label: 'solution' });
  const grid = h('div', 'qd-grid2', host);
  const mk = (title, aria) => { const p = h('div', 'qd-panel', grid); h('div', 'qd-title', p, title); const d = h('div', null, p); d.setAttribute('aria-label', aria); return new LinePlot(d, { width: 480, height: 300 }); };
  const pT = mk('Temperature and density', 'temperature and density versus radius');
  const pG = mk('Gradients', 'actual, adiabatic and radiative temperature gradients versus temperature');
  const pE = mk('Local Eddington factor Γ', 'local Eddington factor versus temperature');
  const pM = mk('Convective Mach number and β', 'convective Mach number and gas pressure fraction versus temperature');
  caption(host, 'Computed in your browser: scan of log R on [14, 16.6] in steps of 0.025, bisection of every sign change of the shooting condition, then Dormand–Prince integration inward.', BOOK + '#sec-structure');

  let roots = [];
  const rootSel = select(bar, { label: 'root', options: [{ value: -1, label: 'largest R' }], value: String(S.root), onChange: (v) => { S.root = +v; drawProfile(); } });
  const solve = async () => {
    status.set('solving…'); await yieldToPaint();
    const t0 = performance.now();
    const out = qs.solve_all(MODEL.M, S.fbh, S.mdot, S.alpha, S.eps, MODEL.logRmin, MODEL.logRmax, MODEL.dlogR);
    roots = [];
    for (let i = 0; i + 10 <= out.length; i += 10) roots.push(Array.from(out.slice(i, i + 10)));
    const sel = rootSel.el.querySelector('select'); sel.innerHTML = '';
    const opt = (v, l) => { const o = document.createElement('option'); o.value = v; o.textContent = l; sel.appendChild(o); };
    opt(-1, 'largest R');
    roots.forEach((r, k) => opt(k, `${k + 1}: R = ${fmt(r[0])} cm`));
    sel.value = String(S.root >= roots.length ? -1 : S.root);
    status.set(`${roots.length} steady solution${roots.length === 1 ? '' : 's'} (${(performance.now() - t0).toFixed(0)} ms)`);
    drawProfile();
  };
  const drawProfile = () => {
    if (!roots.length) { for (const p of [pT, pG, pE, pM]) { p.setSeries([]); p.draw(); } onChange(S); return; }
    const k = S.root < 0 || S.root >= roots.length ? roots.length - 1 : S.root;
    const r = roots[k];
    const pr = qs.envelope_profile(MODEL.M, S.fbh, S.mdot, S.alpha, S.eps, r[0]);
    const n = pr[0], f = pr[1];
    const col = (j) => Float64Array.from({ length: n }, (_, i) => pr[4 + i * f + j]);
    const rr = col(0), T = col(3), rho = col(4), nab = col(5), nad = col(6), nrad = col(7), v = col(8), cs = col(9), G = col(11), beta = col(12);
    const rR = Float64Array.from(rr, (x) => x / r[0]);
    pT.setAxes({ xlog: true, ylog: true, xlabel: 'r / R', ylabel: 'T [K],  ρ × 10¹⁷ [g cm⁻³]' });
    pT.setSeries([{ x: rR, y: T, label: 'T' }, { x: rR, y: Float64Array.from(rho, (x) => x * 1e17), label: 'ρ × 10¹⁷' }]);
    pG.setAxes({ xlog: true, ylog: true, xlabel: 'T [K]', ylabel: 'd ln T / d ln P' });
    pG.setSeries([{ x: T, y: nrad, label: '∇_rad' }, { x: T, y: nab, label: '∇' }, { x: T, y: nad, label: '∇_ad', dash: '4 3' }]);
    pE.setAxes({ xlog: true, ylog: true, xlabel: 'T [K]', ylabel: 'Γ = κL / (4πcGm)' });
    pE.setSeries([{ x: T, y: G, label: 'Γ' }, { x: T, y: Float64Array.from(T, () => 1), label: 'Γ = 1', dash: '4 3', color: '#999' }]);
    pM.setAxes({ xlog: true, ylog: true, xlabel: 'T [K]', ylabel: 'v_c / c_s,  β' });
    pM.setSeries([{ x: T, y: Float64Array.from(v, (x, i) => x / cs[i]), label: 'Mach v_c/c_s' }, { x: T, y: beta, label: 'β = P_gas/P' }]);
    for (const p of [pT, pG, pE, pM]) p.draw();
    const st = pr[2] === 1 ? 'reached r_in' : pr[2] === 2 ? 'mass exhausted' : 'failed';
    status.set(`${roots.length} solution(s); shown: R = ${fmt(r[0])} cm, T_b = ${fmt(r[1], 4)} K, τ_in = ${fmt(r[5])}, r_ph/R = ${fmt(r[6] / r[0])}, T_obs = ${fmt(r[7], 4)} K, Γ_κ,b = ${fmt(r[8])}; integration ${st}`);
    onChange(S);
  };
  const re = () => solve();
  select(bar, { label: 'M_BH/M', options: [{ value: 0.1, label: '0.1' }, { value: 0.03, label: '0.03' }, { value: 0.01, label: '0.01' }], value: String(S.fbh), onChange: (v) => { S.fbh = +v; re(); } });
  slider(bar, { label: 'α_MLT', min: 0.3, max: 2.0, step: 0.05, value: S.alpha, onInput: (v) => { S.alpha = v; clearTimeout(slider._t); slider._t = setTimeout(re, 250); } });
  select(bar, { label: 'ε_in', options: [1e-4, 1e-3, 1e-2].map((v) => ({ value: v, label: String(v) })), value: String(S.eps), onChange: (v) => { S.eps = +v; re(); } });
  select(bar, { label: 'Ṁ [M☉/yr]', options: [0, 0.03, 0.3, 1, 3, 5, 9].map((v) => ({ value: v, label: String(v) })), value: String(S.mdot), onChange: (v) => { S.mdot = +v; re(); } });
  await solve();
  return { draw: solve };
}
