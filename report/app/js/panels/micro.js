// micro.js - zero-metallicity microphysics at fixed density: ionization, mean opacities, ∇_ad and Γ₁ along T.
import { loadCore, h, caption } from '../core.js';
import { LinePlot } from '../viz/plot.js';
import { slider } from '../viz/ui.js';

const BOOK = new URL('../../../chapters/02-microphysics.html', import.meta.url).href;

export async function mountMicro(host, S, onChange = () => {}) {
  const qs = await loadCore();
  const bar = h('div', 'qd-controls', host);
  const grid = h('div', 'qd-grid3', host);
  const mk = (title, aria) => { const p = h('div', 'qd-panel', grid); h('div', 'qd-title', p, title); const d = h('div', null, p); d.setAttribute('aria-label', aria); return new LinePlot(d, { width: 420, height: 300 }); };
  const pIon = mk('Ionization', 'ionization fractions versus temperature');
  const pKap = mk('Mean opacities', 'Rosseland and Planck mean opacity versus temperature');
  const pTh = mk('Adiabatic gradient and Γ₁', 'adiabatic gradient and first adiabatic exponent versus temperature');
  caption(host, 'Saha equilibrium of H/He (X = 0.76) with blackbody radiation; continuum mean opacities from the tabulated means (no lines, no molecules).', BOOK + '#sec-saha');

  const n = 260, lt0 = 3.3, lt1 = 6.0;
  const T = Float64Array.from({ length: n }, (_, i) => 10 ** (lt0 + (lt1 - lt0) * i / (n - 1)));
  const draw = () => {
    const rho = 10 ** S.lrho;
    const xH = new Float64Array(n), y1 = new Float64Array(n), y2 = new Float64Array(n);
    const kR = new Float64Array(n), kP = new Float64Array(n), nad = new Float64Array(n), g1 = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const e = qs.eos_state(rho, T[i]);
      xH[i] = e[11]; y1[i] = e[12]; y2[i] = e[13]; nad[i] = e[4]; g1[i] = e[5];
      const k = qs.kappa(rho, T[i]);
      kR[i] = k[2] ? k[0] : NaN; kP[i] = k[2] ? k[1] : NaN;      // only inside the table
    }
    pIon.setAxes({ xlog: true, xlabel: 'T [K]', ylabel: 'fraction', ylim: [0, 1.02] });
    pIon.setSeries([{ x: T, y: xH, label: 'H⁺/H' }, { x: T, y: y1, label: 'He⁺/He' }, { x: T, y: y2, label: 'He²⁺/He' }]);
    pKap.setAxes({ xlog: true, ylog: true, xlabel: 'T [K]', ylabel: 'κ [cm² g⁻¹]' });
    pKap.setSeries([{ x: T, y: kR, label: 'κ_R (with scattering)' }, { x: T, y: kP, label: 'κ_P (absorption)' },
      { x: T, y: Float64Array.from(T, () => 0.350), label: 'Thomson 0.35', dash: '4 3', color: '#999' }]);
    pTh.setAxes({ xlog: true, xlabel: 'T [K]', ylabel: '∇_ad, Γ₁' });
    pTh.setSeries([{ x: T, y: nad, label: '∇_ad' }, { x: T, y: g1, label: 'Γ₁' }]);
    for (const p of [pIon, pKap, pTh]) p.draw();
    onChange(S);
  };
  slider(bar, { label: 'log₁₀ ρ [g cm⁻³]', min: -18, max: -6, step: 0.1, value: S.lrho, onInput: (v) => { S.lrho = v; draw(); } });
  h('span', 'qd-note', bar, 'Gaps in κ: outside the table (log R = log ρ − 3 log T₆ ∈ [−12, 1]).');
  draw();
  return { draw };
}
