// spec.js - emergent continuum spectra of TITAN snapshots (julia/Spectrum.jl, spectra_run.jl), stored, not recomputed.
import { loadData, h, caption, fmt } from '../core.js';
import { LinePlot } from '../viz/plot.js';
import { select, checkbox } from '../viz/ui.js';

const BOOK = new URL('../../../chapters/08-spectra.html', import.meta.url).href;
const VIEWS = {
  sed: 'νL_ν [erg s⁻¹]',
  ratio: 'L_ν / L_ν(steady star at t = 0)',
  T: 'T [K] used for the spectrum (after the Unsöld–Lucy correction)',
  xH: 'H ionized fraction used for the spectrum (Saha, or non-LTE)',
  full: 'νL_ν with the hydrogen lines (non-LTE spectra only) [erg s⁻¹]',
  Halpha: 'Hα from the infall (estimate), F / F_continuum',
  Hbeta: 'Hβ from the infall (estimate), F / F_continuum',
};
const LVAR = [['lte_b1', 'LTE source function, b₂ = 1', '#0072B2', null], ['lte_b10', 'LTE, b₂ = 10', '#0072B2', '5 3'],
  ['lte_b01', 'LTE, b₂ = 0.1', '#0072B2', '1 3'], ['scat_b1', 'pure scattering, b₂ = 1', '#D55E00', null],
  ['scat_b10', 'pure scattering, b₂ = 10', '#D55E00', '5 3'], ['scat_b01', 'pure scattering, b₂ = 0.1', '#D55E00', '1 3']];
const XAX = { lam: 'λ [Å]', E: 'photon energy [keV]' };
const COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7'];
const hP = 6.62607015e-27, kB = 1.380649e-16, cc = 2.99792458e10, SIGMA = 5.670374419e-5;

function blackbody(nu, T, L) {             // νL_ν of a blackbody sphere of luminosity L and temperature T
  const R2 = L / (4 * Math.PI * SIGMA * T ** 4);
  return Array.from(nu, (v) => 4 * Math.PI ** 2 * R2 * v * (2 * hP * v ** 3 / cc ** 2) / Math.expm1(hP * v / (kB * T)));
}

export async function mountSpec(host, S, onChange = () => {}, opts = {}) {
  const D = await loadData('spectra.json');
  const byName = Object.fromEntries(D.runs.map((r) => [r.name, r]));
  if (opts.spec && byName[opts.spec]) S.spec = opts.spec;
  if (!byName[S.spec]) S.spec = 'A220_late_ext';
  if (S.spec2 !== 'none' && !byName[S.spec2]) S.spec2 = 'star_t0';
  if (!VIEWS[S.sview]) S.sview = 'sed';
  const bar = h('div', 'qd-controls', host);
  const bar2 = h('div', 'qd-controls', host);
  const box = h('div', 'qd-panel', host);
  const d = h('div', null, box); d.setAttribute('aria-label', 'emergent spectrum of a TITAN snapshot');
  const plot = new LinePlot(d, { width: 760, height: 380 });
  const tbox = h('div', 'qd-panel', host);
  const table = h('table', null, tbox); table.style.cssText = 'font-size:0.85em;border-collapse:collapse;width:100%';
  caption(host, 'Static-frame LTE continuum (no lines) of stored TITAN snapshots: spherical rays, Feautrier/ALI with coherent Thomson scattering, Unsöld–Lucy correction of the hydrostatic atmosphere below the shock. Infall in LTE at the TITAN temperature, or (limit) fully ionized; beyond the grid a free-fall extension with an assumed T ∝ r^-1/2. The non-LTE spectrum has non-LTE hydrogen in the infall; its line luminosity is listed separately. Stored output, not recomputed here.', BOOK + '#sec-spec-method');
  const opt = (r) => ({ value: r.name, label: r.label });
  const draw = () => {
    const A = byName[S.spec], B = S.spec2 === 'none' ? null : byName[S.spec2];
    const runs = B ? [A, B] : [A];
    const xk = S.sxax === 'E' ? 'E' : 'lam';
    const series = [];
    if (S.sview === 'sed' || S.sview === 'ratio') {
      const ref = byName.star_t0;
      runs.forEach((r, i) => {
        const x = r[xk];
        const y = S.sview === 'sed' ? r.nuLnu : r.nuLnu.map((v, j) => (ref.nuLnu[j] > 0 ? v / ref.nuLnu[j] : NaN));
        series.push({ x, y, label: r.label, color: COLORS[i] });
      });
      if (S.sview === 'sed' && S.sbb) {
        runs.forEach((r, i) => {
          const nu = r.lam.map((l) => cc / (l * 1e-8));
          series.push({ x: r[xk], y: blackbody(nu, r.Teff, r.Ltitan), label: `blackbody T = ${fmt(r.Teff, 4)} K, L = TITAN L`, color: COLORS[i], dash: '5 3', width: 1 });
        });
      }
      const ymax = Math.max(...series.flatMap((s) => Array.from(s.y).filter(Number.isFinite)));
      plot.setAxes({ xlabel: XAX[xk], ylabel: VIEWS[S.sview], xlog: true, ylog: true,
        ylim: S.sview === 'sed' ? [ymax * 1e-12, ymax * 3] : [1e-4, 1e4], title: '' });
    } else if (S.sview === 'full') {
      runs.forEach((r, i) => {
        const ls = r.linespec;
        if (ls && ls.lam.length) series.push({ x: ls.lam, y: ls.nuLnu, label: r.label + ' (continuum + H lines)', color: COLORS[i] });
        else series.push({ x: r.lam, y: r.nuLnu, label: r.label + ' (continuum only)', color: COLORS[i], dash: '4 3' });
      });
      const ym = Math.max(...series.flatMap((q) => Array.from(q.y).filter((v) => Number.isFinite(v) && v > 0)));
      plot.setAxes({ xlabel: 'λ [Å] (vacuum)', ylabel: VIEWS.full, xlog: true, ylog: true, xlim: [900, 30000], ylim: [ym * 1e-4, ym * 3], title: '' });
    } else if (S.sview === 'Halpha' || S.sview === 'Hbeta') {
      // a non-LTE solution has one profile (its own populations and source function); the LTE spectra show the six
      // brackets (b₂ × thermal or scattering source function); a comparison spectrum adds its profiles, thinner
      const add = (r, thin) => {
        const ln = (r.lines || []).find((l) => l.name === S.sview);
        if (!ln) return false;
        if (ln.nlte && ln.nlte.length) series.push({ x: ln.v, y: ln.nlte, label: `non-LTE: ${r.label}`, color: '#009E73', width: thin ? 1 : 2 });
        else LVAR.forEach(([k, lab, col, dash]) => series.push({ x: ln.v, y: ln[k], label: (thin ? '(compare) ' : '') + lab, color: col, dash, width: thin ? 0.8 : 1.5 }));
        return true;
      };
      const ok = add(A, false); if (B) add(B, true);
      plot.setAxes({ xlabel: 'velocity [km/s] (positive = redshift)', ylabel: VIEWS[S.sview] + (ok ? '' : ' — not computed for this spectrum'), xlog: false, ylog: false, title: A.label });
    } else {
      runs.forEach((r, i) => series.push({ x: r.prof.r, y: S.sview === 'T' ? r.prof.T : r.prof.xH, label: r.label, color: COLORS[i] }));
      if (S.sview === 'T') runs.forEach((r, i) => series.push({ x: r.prof.r, y: r.prof.T_titan, label: 'TITAN T (before the correction)', color: COLORS[i], dash: '4 3', width: 1 }));
      plot.setAxes({ xlabel: 'r / R', ylabel: VIEWS[S.sview], xlog: true, ylog: S.sview === 'T', title: '' });
    }
    plot.setSeries(series);
    plot.setMarkers(S.sview === 'T' || S.sview === 'xH' ? runs.map((r, i) => ({ x: r.rph, label: 'photosphere', color: COLORS[i] })) : []);
    plot.draw();
    const head = ['spectrum', 'continuum L / L_TITAN', 'lines (non-LTE)', 'optical 3646–7000 Å', 'FUV 912–2000 Å', 'X 0.5–2 keV', 'X 2–10 keV', 'X 10–40 keV', 'Balmer jump [dex]'];
    const cell = (v) => (v > 1e30 ? fmt(v, 2) : '—');
    table.innerHTML = '<caption style="caption-side:top;text-align:left">band luminosities [erg s⁻¹]</caption><tr>' +
      head.map((x) => `<th style="text-align:left;padding:2px 6px">${x}</th>`).join('') + '</tr>' +
      runs.map((r) => '<tr>' + [r.label, fmt(r.bands.total / r.Ltitan, 4), r.Lline > 0 ? fmt(r.Lline, 2) : '—', cell(r.bands.opt), cell(r.bands.FUV), cell(r.bands.X_05_2),
        cell(r.bands.X_2_10), cell(r.bands.X_10_40), fmt(r.Dbalmer, 3)].map((x) => `<td style="padding:2px 6px">${x}</td>`).join('') + '</tr>').join('');
    onChange(S);
  };
  select(bar, { label: 'spectrum', options: D.runs.map(opt), value: S.spec, onChange: (v) => { S.spec = v; draw(); } });
  select(bar, { label: 'compare with', options: [{ value: 'none', label: 'none' }, ...D.runs.map(opt)], value: S.spec2, onChange: (v) => { S.spec2 = v; draw(); } });
  select(bar2, { label: 'show', options: Object.keys(VIEWS).map((k) => ({ value: k, label: VIEWS[k] })), value: S.sview, onChange: (v) => { S.sview = v; plot.resetView(); draw(); } });
  select(bar2, { label: 'x axis', options: Object.keys(XAX).map((k) => ({ value: k, label: XAX[k] })), value: S.sxax, onChange: (v) => { S.sxax = v; plot.resetView(); draw(); } });
  checkbox(bar2, { label: 'blackbody', value: !!S.sbb, onChange: (v) => { S.sbb = v ? 1 : 0; draw(); } });
  draw();
  return { draw };
}
