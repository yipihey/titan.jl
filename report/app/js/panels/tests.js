// tests.js - the decisive tests of the second book: Monte Carlo line profiles with electron scattering, the Balmer-limit
// region with the full series, axisymmetric infall + outflow profiles, He I λ10830. Stored results (tests.json).
import { loadData, h, caption, fmt } from '../core.js';
import { LinePlot, PALETTE } from '../viz/plot.js';
import { select, checkbox } from '../viz/ui.js';

const BOOK = new URL('../../../spectra/06-tests.html', import.meta.url).href;
const VIEWS = {
  mc: 'Hα / Hβ escaping profile, Monte Carlo (normalized to unit area, per km/s)',
  mcincl: 'Hα by inclination for infall + polar outflow (Monte Carlo)',
  balmer: 'the Balmer-limit region: L_ν / continuum L_ν at 4050 Å',
  axisym: 'Hα / Hβ from infall + polar outflow (ray integration, F/F_c)',
  he: 'He I λ10830, F/F_c',
  cool: 'the cooler star C30 beside W30',
};
const SEL = {
  mc: ['Halpha_es0_res1', 'Halpha_es1_res1', 'Halpha_es1_res1_th90', 'Halpha_es1_res1_th60'],
};

export async function mountTests(host, S, onChange = () => {}, opts = {}) {
  const D = await loadData('tests.json');
  if (opts.view && VIEWS[opts.view]) S.tview = opts.view;
  if (!VIEWS[S.tview]) S.tview = 'mc';
  const bar = h('div', 'qd-controls', host);
  const bar2 = h('div', 'qd-controls', host);
  const box = h('div', 'qd-panel', host);
  const d = h('div', null, box); d.setAttribute('aria-label', 'decisive tests of the fed quasi-star spectrum');
  const plot = new LinePlot(d, { width: 760, height: 380 });
  const note = h('div', 'qd-panel', host); note.style.fontSize = '0.85em';
  caption(host, 'Stored results of julia/mc_escatter.jl, balmer_series.jl, axisym_lines.jl, he10830.jl and nlte_spectrum.jl for the non-LTE solutions of TITAN runs W30 and C30 (the cooler star). Not recomputed here.', BOOK);
  const ctl = [];
  const draw = () => {
    ctl.forEach((c) => c.el.remove()); ctl.length = 0;
    const series = [];
    let xlabel = 'velocity [km/s] (positive = redshift)', ylog = false, xlim = null, txt = '';
    if (S.tview === 'mc') {
      const tags = S.tmcall ? D.mc.map((m) => m.tag) : SEL.mc;
      D.mc.filter((m) => tags.includes(m.tag)).forEach((m, i) => series.push({ x: m.V, y: m.prof, label: m.label, color: PALETTE[i % PALETTE.length] }));
      ctl.push(checkbox(bar2, { label: 'all runs', value: !!S.tmcall, onChange: (v) => { S.tmcall = v ? 1 : 0; draw(); } }));
      ylog = true; xlim = [-4000, 4000];
      txt = D.mc.map((m) => `${m.label}: FWHM ${fmt(m.fwhm, 3)} km/s, centroid ${fmt(m.centroid, 3)} km/s, exponential wing scale blue ${fmt(m.Wb, 3)} / red ${fmt(m.Wr, 3)} km/s`).join('<br>');
    } else if (S.tview === 'mcincl') {
      D.mc.filter((m) => m.profi && m.profi.length).forEach((m, j) => ['0–30°', '30–60°', '60–90°'].forEach((lab, b) =>
        series.push({ x: m.V, y: m.profi[b], label: `${m.label}, inclination ${lab}`, color: PALETTE[(3 * j + b) % PALETTE.length], dash: j ? '5 3' : null })));
      ylog = true; xlim = [-4000, 4000];
      txt = 'Inclination = angle between the line of sight and the outflow axis.';
    } else if (S.tview === 'balmer') {
      const b = D.balmer.find((x) => x.tag === (S.tbal || 'n40_b12')) || D.balmer[0];
      if (b) {
        // L_ν normalized to the continuum at 4050 Å (the ratio to the continuum would mix the jump into the smoothed curves)
        const k0 = b.lam.reduce((m, l, i) => (Math.abs(l - 4050) < Math.abs(b.lam[m] - 4050) ? i : m), 0), n0 = b.Lc[k0];
        const ok = (a) => a && a.length && a.every((v) => v !== null && Number.isFinite(v));
        series.push({ x: b.lam, y: b.Lc.map((v) => v / n0), label: 'continuum only', color: PALETTE[3], dash: '5 3' });
        series.push({ x: b.lam, y: b.Lnu.map((v) => v / n0), label: 'full resolution', color: PALETTE[0], width: 1 });
        if (ok(b.L100)) series.push({ x: b.lam, y: b.L100.map((v) => v / n0), label: 'R = 100', color: PALETTE[1], width: 2 });
        if (ok(b.L30)) series.push({ x: b.lam, y: b.L30.map((v) => v / n0), label: 'R = 30', color: PALETTE[2], width: 2 });
        txt = `Break ⟨f_ν⟩(4000–4100)/⟨f_ν⟩(3620–3720): continuum alone ${fmt(b.brk.cont, 3)}; with the series ${fmt(b.brk.full, 3)}` +
          (ok(b.L100) ? ` (R=100 ${fmt(b.brk.R100, 3)}, R=30 ${fmt(b.brk.R30, 3)}).` : ' (break windows only, no smoothing).');
      }
      ctl.push(select(bar2, { label: 'run', options: D.balmer.map((x) => ({ value: x.tag, label: x.tag })), value: S.tbal || 'n40_b12', onChange: (v) => { S.tbal = v; draw(); } }));
      xlabel = 'λ [Å] (vacuum)'; xlim = [3550, 4150];
    } else if (S.tview === 'axisym') {
      const ln = S.tline || 'Halpha';
      D.axisym.filter((x) => x.line === ln).forEach((x, i) => series.push({ x: x.v, y: x.F, label: `θ_o = ${x.theta}°, i = ${x.incl}° (EW ${fmt(x.EW, 3)} Å)`, color: PALETTE[i % PALETTE.length] }));
      ctl.push(select(bar2, { label: 'line', options: [{ value: 'Halpha', label: 'Hα' }, { value: 'Hbeta', label: 'Hβ' }], value: ln, onChange: (v) => { S.tline = v; draw(); } }));
      txt = 'θ_o: half opening angle of the polar outflow cones (0 = pure infall, 90 = pure outflow); i: angle of the line of sight to the axis.';
      xlim = [-2500, 2500];
    } else if (S.tview === 'cool') {
      const q = S.tcool || 'sed';
      (D.cool || []).forEach((c, i) => {
        if (q === 'sed') series.push({ x: c.lam, y: c.nuLnu, label: c.label, color: PALETTE[i % PALETTE.length], width: c.name.endsWith('nlte') ? 1.5 : 1, dash: c.name.endsWith('nlte') ? null : '5 3' });
        else if (c.xH.length) {
          const k = c.r.map((r, j) => j).filter((j) => c.r[j] >= 0.95);
          series.push({ x: k.map((j) => c.r[j]), y: k.map((j) => c.xH[j]), label: c.label, color: PALETTE[i % PALETTE.length] });
        }
      });
      ctl.push(select(bar2, { label: 'plot', options: [{ value: 'sed', label: 'spectrum νL_ν [erg/s] vs λ [Å]' }, { value: 'x', label: 'ionized fraction x_H vs r / R_star' }], value: q, onChange: (v) => { S.tcool = v; plot.resetView(); draw(); } }));
      if (q === 'sed') { xlabel = 'λ [Å] (vacuum)'; ylog = true; xlim = [900, 25000]; }
      else { xlabel = 'r / R (each star its own R)'; ylog = true; }
      const ylim = q === 'sed' ? [1e40, 3e44] : [1e-3, 1.2];
      plot.setAxes({ xlabel, ylabel: q === 'sed' ? 'νL_ν [erg/s]' : 'x_H', xlog: true, ylog, xlim, ylim, title: '' });
      plot.setSeries(series); plot.draw();
      note.innerHTML = 'Dashed: the stars alone. Solid: fed, with non-LTE hydrogen in the infall (lines included where computed). The C30 non-LTE solution has no Unsöld–Lucy correction and loses 41 % of TITAN\'s luminosity (see the book).';
      onChange(S);
      return;
    } else {
      D.he.forEach((x, i) => series.push({ x: x.v, y: x.F, label: `${x.name}: EW ${fmt(x.EW, 3)} Å, N(2³S) ${fmt(x.N3, 2)} cm⁻²`, color: PALETTE[i] }));
      xlim = [-2500, 2500];
      txt = 'EW &lt; 0: net emission (positive = absorption). N(2³S): radial column of the metastable level.';
    }
    plot.setAxes({ xlabel, ylabel: VIEWS[S.tview], xlog: false, ylog, xlim, title: '' });
    plot.setSeries(series);
    plot.draw();
    note.innerHTML = txt;
    onChange(S);
  };
  select(bar, { label: 'show', options: Object.keys(VIEWS).map((k) => ({ value: k, label: VIEWS[k] })), value: S.tview, onChange: (v) => { S.tview = v; plot.resetView(); draw(); } });
  draw();
  return { draw };
}
