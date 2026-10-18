// analytic.js - the simple models of the infall's hydrogen (second book): the closed form, evaluated here, against
// the stored semi-analytic model (julia/analytic_infall.jl) and the full non-LTE solution of run W30 (nlte_run.jl).
import { loadData, h, caption, fmt } from '../core.js';
import { LinePlot } from '../viz/plot.js';
import { select, slider, checkbox, readout } from '../viz/ui.js';

const BOOK = new URL('../../../spectra/03-simple-model.html', import.meta.url).href;
const VIEWS = {
  x: 'ionized fraction of hydrogen x_H',
  dL: 'Hα power per ln r, dL(Hα)/d ln r [erg s⁻¹]',
  SB: 'Hα source function S_L / B_ν(Hα, T_eff) (full solution)',
  fuv: 'x_H with the star’s far UV (< 1300 Å) scaled (semi-analytic model)',
};
const C = { full: '#0072B2', simple: '#009E73', closed: '#D55E00', live: '#CC79A7', bb: '#999999' };
const SIGT = 6.6524e-25;

// closed form: x²/(1 − x) = f² Q₀ q / (n_H α_B(T)), f = W + 3τ* R*²/(4r²) (dilution × Thomson trapping), iterated in τ*
function closedForm(A, mdot, fuvq, trap) {
  const r = A.full.r.map((x) => x * A.R0), T = A.full.T, n = r.length;
  const vff = (x) => Math.sqrt(2 * A.G * A.M / x);
  const nH = r.map((x) => A.X * mdot * A.Msun / A.yr / (4 * Math.PI * x * x * vff(x)) / A.mH);
  const W = (x) => 0.5 * (1 - Math.sqrt(Math.max(1 - (A.Rst / x) ** 2, 0)));
  const aB = (t) => 2.59e-13 * (t / 1e4) ** -0.8;
  const Q0 = (A.p2 * A.G2 + A.p3 * A.G3) * fuvq;
  let tau = new Float64Array(n); const x = new Float64Array(n);
  for (let it = 0; it < 80; it++) {
    for (let k = 0; k < n; k++) {
      const f = W(r[k]) + (trap ? 3 * tau[k] * A.Rst ** 2 / (4 * r[k] ** 2) : 0);
      const a = f * f * Q0 / (nH[k] * aB(T[k]));
      x[k] = (-a + Math.sqrt(a * a + 4 * a)) / 2;
    }
    const tn = new Float64Array(n); let acc = 0;
    for (let k = n - 2; k >= 0; k--) { acc += 0.5 * (x[k] * nH[k] / r[k] ** 2 + x[k + 1] * nH[k + 1] / r[k + 1] ** 2) * SIGT * (r[k + 1] - r[k]); tn[k] = r[k] ** 2 * acc; }
    let d = 0; for (let k = 0; k < n; k++) { d = Math.max(d, Math.abs(tn[k] - tau[k])); tau[k] = 0.5 * (tau[k] + tn[k]); }
    if (!trap || d < 1e-4) break;
  }
  let NHI = 0, tauT = 0;
  for (let k = 0; k < n - 1; k++) {
    const dr = r[k + 1] - r[k];
    NHI += 0.5 * ((1 - x[k]) * nH[k] + (1 - x[k + 1]) * nH[k + 1]) * dr;
    tauT += 0.5 * (x[k] * nH[k] + x[k + 1] * nH[k + 1]) * SIGT * dr;
  }
  return { x, NHI, tauT };
}

export async function mountAnalytic(host, S, onChange = () => {}, opts = {}) {
  const A = await loadData('analytic.json');
  if (opts.view && VIEWS[opts.view]) S.aview = opts.view;
  if (!VIEWS[S.aview]) S.aview = 'x';
  const bar = h('div', 'qd-controls', host);
  const bar2 = h('div', 'qd-controls', host);
  const info = readout(host, { label: 'summary' });
  const box = h('div', 'qd-panel', host);
  const d = h('div', null, box); d.setAttribute('aria-label', 'simple models of the infall against the full non-LTE solution');
  const plot = new LinePlot(d, { width: 760, height: 380 });
  caption(host, 'Full: non-LTE hydrogen with full transfer in TITAN run W30 (stored). Semi-analytic: the same 12-level statistical equilibrium in the star’s diluted, Thomson-trapped light with case B, no transfer (stored). Closed form: x²/(1−x) = f²Q₀/(n_H α_B), evaluated here; Ṁ and the far UV can be changed (the temperatures stay those of W30).', BOOK + '#sec-sm-closed');
  const draw = () => {
    const series = [];
    let ylog = false, xlog = true, ylabel = VIEWS[S.aview];
    if (S.aview === 'x') {
      const cf = closedForm(A, S.amdot, S.afuv, S.atrap);
      series.push({ x: A.full.r, y: A.full.x, label: 'full non-LTE (W30)', color: C.full, width: 2 });
      series.push({ x: A.simple.r, y: A.simple.x, label: 'semi-analytic (star_t0, case B, trapping)', color: C.simple });
      series.push({ x: A.simple_bb.r, y: A.simple_bb.x, label: 'semi-analytic with a 7310 K blackbody star', color: C.bb, dash: '2 3' });
      series.push({ x: A.full.r, y: Array.from(cf.x), label: `closed form (Ṁ = ${fmt(S.amdot, 2)} M⊙/yr, far UV × ${fmt(S.afuv, 2)}${S.atrap ? '' : ', no trapping'})`, color: C.live, dash: '6 3' });
      info.set(`neutral column ${fmt(cf.NHI, 3)} cm⁻², Thomson depth ${fmt(cf.tauT, 3)} (full: 8.0e23 and 5.06 on the same points)`);
      ylog = S.alog;
    } else if (S.aview === 'dL') {
      series.push({ x: A.full.r, y: A.full.dL, label: 'full non-LTE (W30)', color: C.full, width: 2 });
      series.push({ x: A.simple.r, y: A.simple.dL, label: 'semi-analytic', color: C.simple });
      info.set(`L(Hα), local net power: full 1.10e43, semi-analytic ${fmt(A.simple.LHa, 3)} erg s⁻¹`);
      ylog = true;
    } else if (S.aview === 'SB') {
      series.push({ x: A.full.r, y: A.full.SB, label: 'full non-LTE (W30)', color: C.full, width: 2 });
      const x0 = 1.2, y0 = A.full.SB[A.full.r.findIndex((r) => r >= x0)];
      series.push({ x: A.full.r, y: A.full.r.map((r) => y0 * (r / x0) ** -1.2), label: 'r^-1.2 (guide)', color: C.bb, dash: '4 3', width: 1 });
      info.set('S_L ∝ r^-1.2 to about 20 R; L(Hα) ∝ ∫ r v (S_L − J) dr then grows only as r_out^0.3');
      ylog = true;
    } else {
      series.push({ x: A.simple.r, y: A.simple.x, label: 'far UV × 1', color: C.simple, width: 2 });
      const cols = ['#56B4E9', '#E69F00', '#999999'];
      A.fuv.forEach((v, i) => series.push({ x: v.s.r, y: v.s.x, label: `far UV × ${v.f} (L(Hα) ${fmt(v.s.LHa, 2)} erg/s)`, color: cols[i] }));
      info.set(`L(Hα) of the semi-analytic model: ${fmt(A.simple.LHa, 3)} erg s⁻¹ at far UV × 1`);
      ylog = S.alog;
    }
    plot.setAxes({ xlabel: 'r / R', ylabel, xlog, ylog, title: '' });
    plot.setSeries(series);
    plot.draw();
    onChange(S);
  };
  select(bar, { label: 'show', options: Object.keys(VIEWS).map((k) => ({ value: k, label: VIEWS[k] })), value: S.aview, onChange: (v) => { S.aview = v; plot.resetView(); draw(); } });
  checkbox(bar, { label: 'log y', value: !!S.alog, onChange: (v) => { S.alog = v ? 1 : 0; plot.resetView(); draw(); } });
  slider(bar2, { label: 'Ṁ [M⊙/yr]', min: 0.1, max: 30, log: true, value: S.amdot, format: (v) => fmt(v, 2), onInput: (v) => { S.amdot = v; draw(); } });
  slider(bar2, { label: 'far UV ×', min: 1e-3, max: 10, log: true, value: S.afuv, format: (v) => fmt(v, 2), onInput: (v) => { S.afuv = v; draw(); } });
  checkbox(bar2, { label: 'Thomson trapping', value: !!S.atrap, onChange: (v) => { S.atrap = v ? 1 : 0; draw(); } });
  draw();
  return { draw };
}
