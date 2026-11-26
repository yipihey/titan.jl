// landmarks.js - Interactive Linked Explorer for Titan.jl Landmark Astrophysics Reproductions (Waves 1-4).
// Click any point on one graph to highlight the exact corresponding state across all linked graphs.
import { loadData, h, caption, fmt } from '../core.js';
import { LinePlot, HeatmapPlot, PALETTE } from '../viz/plot.js';
import { select } from '../viz/ui.js';

const SUITES = [
  { value: 'ensman',    label: 'Wave 1a: Ensman (1994) Radiative Shocks (SUBC / SUPC)' },
  { value: 'omukai',    label: 'Wave 1b: Omukai (1998/2001) Collapse & Ripamonti-Abel H₂ Lines' },
  { value: 'whalen',    label: 'Wave 2a: Whalen+04 / Kitayama+04 / Alvarez+06 Photoevaporation' },
  { value: 'schaerer',  label: 'Wave 2b: Schaerer (2002) Pop III Grid & Raiter+10 High-n_e' },
  { value: 'sms',       label: 'Wave 3a: Woods+17 / Haemmerlé+18 SMS GR Collapse & Live 1D GR' },
  { value: 'quasistar', label: 'Wave 3b: Begelman+08 Quasi-Stars & Freese+16 Dark Stars (JWST)' },
  { value: 'pisn',      label: 'Wave 4: Heger+02 Pair Instability, Kasen+11 PISN & Woosley+07 PPISN' },
];

const REFS = {
  ensman:    'Ensman (1994, ApJ 424, 275); Stahler, Shu & Taam (1980, ApJ 241, 637). Click any spatial zone to highlight T_gas, T_rad, ρ, and velocity u simultaneously.',
  omukai:    'Omukai & Nishi (1998, ApJ 508, 141); Omukai (2001, ApJ 546, 635); Ripamonti & Abel (2004, MNRAS 348, 1019). Click any density n_H to link thermal track, H₂ fraction, and ExoMol line luminosities.',
  whalen:    'Whalen, Abel & Norman (2004, ApJ 610, 14); Kitayama+04 (ApJ 613, 631); Alvarez+06 (ApJ 639, 621). Click any time to inspect f_esc, f_bound, core evacuation, and nebular line fading across halo masses.',
  schaerer:  'Schaerer (2002, A&A 382, 28); Raiter, Schaerer & Fosbury (2010, A&A 523, A64). Click any Pop III mass M_* or density n_e to compare Titan.jl against Schaerer Table 3 & Raiter+10.',
  sms:       'Woods+17 (ApJ 842, L6); Umeda+16 (ApJ 830, L34); Haemmerlé+18 (MNRAS 474, 2757); Chandrasekhar (1964, ApJ 140, 417). Click any accretion rate or collapse time to compare Full 1PN/TOV vs Truncated KEPLER GR.',
  quasistar: 'Begelman+08 (MNRAS 387, 1649); Ball+11 (MNRAS 414, 2751); Spolyar, Freese & Gondolo (2008, PRL 100, 051101); Freese+16 (RoPP 79, 066902). Click any object on the HR locus or wavelength on the JWST SED.',
  pisn:      'Heger & Woosley (2002, ApJ 567, 532); Kasen, Woosley & Heger (2011, ApJ 734, 102); Woosley, Blinnikov & Heger (2007, Nature 450, 390). Click the Γ₁ pair valley or lightcurves to inspect shock breakout and PPISN shell collision.',
};

function setTable(tbl, headers, rows, hlIdx = -1) {
  tbl.innerHTML = '<tr>' + headers.map((x) => `<th>${x}</th>`).join('') + '</tr>' +
    rows.map((r, i) => `<tr class="${i === hlIdx ? 'is-hl' : ''}">` + r.map((c) => `<td>${c}</td>`).join('') + '</tr>').join('');
}

export async function mountLandmarks(host, S, onChange = () => {}, opts = {}) {
  const suiteInit = opts.suite || S.lmsuite || 'ensman';
  S.lmsuite = SUITES.some((x) => x.value === suiteInit) ? suiteInit : 'ensman';

  const bar = h('div', 'qd-controls', host);
  const selbar = h('div', 'qd-selbar', host);
  const stage = h('div', null, host);
  const tbox = h('div', 'qd-panel', host);
  const tbl = h('table', 'qd-table', tbox);
  const cap = h('p', 'qd-note', host);

  let plots = [];
  const clearStage = () => { plots.forEach((p) => p.destroy && p.destroy()); plots = []; stage.innerHTML = ''; };

  async function render() {
    clearStage();
    const k = S.lmsuite;
    cap.textContent = REFS[k] || '';
    selbar.innerHTML = '<span><strong>Click any point</strong> on a graph below to highlight the corresponding physical state across all linked plots.</span>';

    if (k === 'ensman') {
      const D = (await loadData('lm_shocks_collapse.json')).ensman;
      const g = h('div', 'qd-grid2', stage);
      const b1 = h('div', 'qd-panel', g), b2 = h('div', 'qd-panel', g);
      const p1 = new LinePlot(b1, { width: 520, height: 340, selectMode: 'x' });
      const p2 = new LinePlot(b2, { width: 520, height: 340, selectMode: 'x' });
      plots = [p1, p2];
      const x1 = D.subc.r_cm.map((v) => v / 1e10), x2 = D.supc.r_cm.map((v) => v / 1e10);
      p1.setAxes({ xlabel: 'x [10¹⁰ cm]', ylabel: 'Temperature [K]', ylog: true, title: 'Subcritical Shock (SUBC, 6 km/s)' });
      p1.setSeries([
        { x: x1, y: D.subc.T_gas, label: 'T_gas (Zel’dovich spike = 847 K)', color: '#D55E00', width: 2 },
        { x: x1, y: D.subc.T_rad, label: 'T_rad (precursor)', color: '#0072B2', dash: '5 3', width: 1.8 },
      ]);
      p2.setAxes({ xlabel: 'x [10¹⁰ cm]', ylabel: 'Temperature [K]', ylog: true, title: 'Supercritical Shock (SUPC, 20 km/s)' });
      p2.setSeries([
        { x: x2, y: D.supc.T_gas, label: 'T_gas (spike = 6485 K)', color: '#D55E00', width: 2 },
        { x: x2, y: D.supc.T_rad, label: 'T_rad (T_pre/T_post = 1.039)', color: '#0072B2', dash: '5 3', width: 1.8 },
      ]);
      const sync = (which, idx) => {
        const d = which === 'subc' ? D.subc : D.supc;
        const i = Math.max(0, Math.min(d.r_cm.length - 1, idx));
        const xv = d.r_cm[i] / 1e10;
        const p = which === 'subc' ? p1 : p2;
        p.setHighlights([
          { x: xv, y: d.T_gas[i], color: '#D55E00', label: `T_gas=${fmt(d.T_gas[i], 4)} K` },
          { x: xv, y: d.T_rad[i], color: '#0072B2', vline: false },
        ]);
        selbar.innerHTML = `<strong>${which.toUpperCase()} zone ${d.zone[i]}</strong> ` +
          `<span>x = ${fmt(xv, 4)}×10¹⁰ cm</span> <span>T_gas = ${fmt(d.T_gas[i], 4)} K</span> ` +
          `<span>T_rad = ${fmt(d.T_rad[i], 4)} K</span> <span>ρ = ${fmt(d.rho[i], 3)} g/cm³</span> <span>u = ${fmt(d.u_kms[i], 3)} km/s</span>`;
      };
      p1.onSelect((hit) => sync('subc', hit.pointIndex));
      p2.onSelect((hit) => sync('supc', hit.pointIndex));
      p1.draw(); p2.draw();
      setTable(tbl, ['Benchmark', 'Published Target (Ensman 1994)', 'Titan.jl Moving-Mesh RHD', 'Status'], [
        ['SUBC (6 km/s)', 'T_spike ≈ 810–850 K; T_pre ≪ T_post', 'T_spike = 846.9 K; T_pre/T_post = 0.278 (3.74× jump)', 'PASS'],
        ['SUPC (20 km/s)', 'T_spike ≈ 6000–6500 K; T_pre ≈ T_post', 'T_spike = 6485 K; T_pre/T_post = 1.039 (6.48× jump, 6.7 s)', 'PASS'],
      ]);
    } else if (k === 'omukai') {
      const D = (await loadData('lm_shocks_collapse.json')).omukai;
      const g = h('div', 'qd-grid2', stage);
      const b1 = h('div', 'qd-panel', g), b2 = h('div', 'qd-panel', g);
      const p1 = new LinePlot(b1, { width: 520, height: 340, selectMode: 'x' });
      const p2 = new LinePlot(b2, { width: 520, height: 340, selectMode: 'x' });
      plots = [p1, p2];
      const J0 = D.runs.J0, HL = D.h2lines;
      p1.setAxes({ xlabel: 'Hydrogen Density n_H [cm⁻³]', ylabel: 'Temperature T [K]', xlog: true, ylog: true, title: 'Omukai (1998/2001): 22-Decade Collapse & J_crit', legendPos: 'tl' });
      p1.setSeries([
        { x: J0.nH, y: J0.T, label: 'J₂₁ = 0 (Pop III minihalo)', color: '#0072B2', width: 2.2 },
        { x: D.runs.J30.nH, y: D.runs.J30.T, label: 'J₂₁ = 30 (sub-critical)', color: '#009E73', dash: '5 3' },
        { x: D.runs.J1e3.nH, y: D.runs.J1e3.T, label: 'J₂₁ = 10³ (T* = 10⁴ K, DCBH)', color: '#E69F00', width: 1.8 },
        { x: D.runs.J1e4.nH, y: D.runs.J1e4.T, label: 'J₂₁ = 10⁴ (atomic cooling)', color: '#D55E00', width: 2.1 },
      ]);
      p2.setAxes({ xlabel: 'Hydrogen Density n_H [cm⁻³]', ylabel: 'Emergent Luminosity [L☉]', xlog: true, ylog: true, ylim: [1e-4, 1e3], title: 'Ripamonti & Abel (2004): ExoMol H₂ Lines vs CIE', legendPos: 'tl' });
      p2.setSeries([
        { x: HL.nH, y: HL.L_H2tot_Lsun, label: 'Total ExoMol H₂ (183 lines)', color: '#0072B2', width: 2.2 },
        { x: HL.nH, y: HL.L_S1_Lsun, label: '0–0 S(1) 17.0 μm', color: '#009E73', dash: '5 3' },
        { x: HL.nH, y: HL.L_S3_Lsun, label: '0–0 S(3) 9.66 μm', color: '#CC79A7', dash: '3 3' },
        { x: HL.nH, y: HL.L_v1S1_Lsun, label: '1–0 S(1) 2.12 μm', color: '#E69F00', dash: '2 2' },
        { x: HL.nH, y: HL.L_cont_Lsun, label: 'Borysow CIA/CIE continuum', color: '#D55E00', width: 2.0 },
      ]);
      const syncN = (nVal) => {
        const nearest = (arr, v) => { let bi = 0, bd = Infinity; for (let i = 0; i < arr.length; i++) { const d = Math.abs(Math.log10(arr[i]) - Math.log10(v)); if (d < bd) { bd = d; bi = i; } } return bi; };
        const i0 = nearest(J0.nH, nVal), i1 = nearest(D.runs.J1e4.nH, nVal), ih = nearest(HL.nH, nVal);
        p1.setHighlights([
          { x: J0.nH[i0], y: J0.T[i0], color: '#0072B2', label: `J=0: ${fmt(J0.T[i0], 3)} K` },
          { x: D.runs.J1e4.nH[i1], y: D.runs.J1e4.T[i1], color: '#D55E00', vline: false },
        ]);
        p2.setHighlights([
          { x: HL.nH[ih], y: HL.L_H2tot_Lsun[ih], color: '#0072B2', label: `H₂: ${fmt(HL.L_H2tot_Lsun[ih], 3)} L☉` },
          { x: HL.nH[ih], y: Math.max(HL.L_cont_Lsun[ih] || 0, 1e-9), color: '#D55E00', vline: false },
        ]);
        selbar.innerHTML = `<strong>n_H = ${fmt(J0.nH[i0], 3)} cm⁻³</strong> ` +
          `<span>T(J=0) = ${fmt(J0.T[i0], 4)} K</span> <span>x_H₂ = ${fmt(J0.xH2[i0], 3)}</span> ` +
          `<span>T(J=10⁴) = ${fmt(D.runs.J1e4.T[i1], 4)} K</span> <span>L_H₂ = ${fmt(HL.L_H2tot_Lsun[ih], 3)} L☉</span> <span>L_CIE = ${fmt(HL.L_cont_Lsun[ih], 3)} L☉</span>`;
      };
      p1.onSelect((hit) => syncN(hit.x));
      p2.onSelect((hit) => syncN(hit.x));
      p1.draw(); p2.draw(); syncN(1e10);
      setTable(tbl, ['Metric', 'Published Target', 'Titan.jl Result'], [
        ['Loitering minimum (J=0)', 'T_min ≈ 200–250 K at n ~ 10⁴ cm⁻³ (Omukai & Nishi 1998)', 'T_min = 236.3 K at 6.9×10³ cm⁻³; max x_H₂ = 0.500'],
        ['LW J_crit bifurcation', 'Atomic-cooling jump at J > J_crit (Omukai 2001)', '66× hotter at n_H = 10⁴ cm⁻³ (J=10⁴ vs J=0)'],
        ['ExoMol H₂ vs CIE crossover', 'H₂ peaks ~10¹⁰; CIE dominates 10¹³–10¹⁶ (Ripamonti & Abel 2004)', 'Peak L_H₂ ≈ 180 L☉ @ 10¹⁰; CIE overtakes 10¹³–10¹⁸ cm⁻³'],
      ]);
    } else if (k === 'whalen') {
      const D = (await loadData('lm_feedback_popiii.json')).whalen;
      const g = h('div', 'qd-grid2', stage);
      const b1 = h('div', 'qd-panel', g), b2 = h('div', 'qd-panel', g);
      const p1 = new LinePlot(b1, { width: 520, height: 340, selectMode: 'x' });
      const p2 = new LinePlot(b2, { width: 520, height: 340, selectMode: 'x' });
      plots = [p1, p2];
      const h3 = D.halos['3e5'], h6 = D.halos['1e6'], h5 = D.halos['5e6'];
      p1.setAxes({ xlabel: 'Time t [Myr]', ylabel: 'f_esc (solid) & f_bound (dashed)', xlog: true, xlim: [0.005, 2.6], ylim: [-0.02, 1.06], title: 'Minihalo Photoevaporation (3×10⁵–5×10⁶ M☉)', legendPos: 'bl' });
      p1.setSeries([
        { x: h3.t_Myr, y: h3.f_esc_H, label: '3×10⁵ M☉ f_esc (0.981)', color: '#009E73', width: 2.1 },
        { x: h3.t_Myr, y: h3.f_bound, label: '3×10⁵ M☉ f_bound (0.397)', color: '#009E73', dash: '5 3', width: 1.6 },
        { x: h6.t_Myr, y: h6.f_esc_H, label: '10⁶ M☉ f_esc', color: '#0072B2', width: 2.1 },
        { x: h6.t_Myr, y: h6.f_bound, label: '10⁶ M☉ f_bound', color: '#0072B2', dash: '5 3', width: 1.6 },
        { x: h5.t_Myr, y: h5.f_bound, label: '5×10⁶ M☉ f_bound (bound)', color: '#D55E00', dash: '4 2', width: 1.8 },
      ]);
      p2.setAxes({ xlabel: 'Time t [Myr]', ylabel: 'Nebular Line Luminosity [L☉] (10⁶ M☉ Halo)', xlog: true, ylog: true, xlim: [0.005, 2.6], ylim: [50, 3e5], title: 'JWST Atomic Lines During Evacuation (ModelAtom)', legendPos: 'bl' });
      p2.setSeries([
        { x: h6.t_Myr, y: h6.L_Lya, label: 'Lyα 1216 Å', color: '#D55E00', width: 2.1 },
        { x: h6.t_Myr, y: h6.L_1640, label: 'He II λ1640', color: '#CC79A7', width: 2.0 },
        { x: h6.t_Myr, y: h6.L_Hb, label: 'Hβ λ4861', color: '#0072B2', width: 1.9 },
        { x: h6.t_Myr, y: h6.L_10830, label: 'He I λ10830 (fades 8×)', color: '#009E73', width: 2.0 },
      ]);
      const syncT = (tVal) => {
        const nearest = (arr, v) => { let bi = 0, bd = Infinity; for (let i = 0; i < arr.length; i++) { const d = Math.abs(arr[i] - v); if (d < bd) { bd = d; bi = i; } } return bi; };
        const i3 = nearest(h3.t_Myr, tVal), i6 = nearest(h6.t_Myr, tVal);
        p1.setHighlights([
          { x: h3.t_Myr[i3], y: h3.f_esc_H[i3], color: '#009E73', label: `f_esc=${fmt(h3.f_esc_H[i3], 3)}` },
          { x: h6.t_Myr[i6], y: h6.f_bound[i6], color: '#0072B2', vline: false },
        ]);
        p2.setHighlights([
          { x: h6.t_Myr[i6], y: h6.L_1640[i6], color: '#CC79A7', label: `1640: ${fmt(h6.L_1640[i6], 3)}` },
          { x: h6.t_Myr[i6], y: h6.L_10830[i6], color: '#009E73', vline: false },
        ]);
        selbar.innerHTML = `<strong>t = ${fmt(h6.t_Myr[i6], 3)} Myr</strong> ` +
          `<span>10⁶ M☉ core n_H = ${fmt(h6.nH_core[i6], 3)} cm⁻³</span> <span>R_HII = ${fmt(h6.R_HII_pc[i6], 3)} pc</span> ` +
          `<span>v_sh = ${fmt(h6.u_max_kms[i6], 3)} km/s</span> <span>He I 10830 = ${fmt(h6.L_10830[i6], 3)} L☉</span>`;
      };
      p1.onSelect((hit) => syncT(hit.x));
      p2.onSelect((hit) => syncT(hit.x));
      p1.draw(); p2.draw(); syncT(0.5);
      setTable(tbl, ['Minihalo Mass', 'Shock v_peak', 'Final Core n_c', 'f_esc(2.5 Myr)', 'f_bound(2.5 Myr)'], [
        ['3×10⁵ M☉', '27.0 km/s', '0.16 cm⁻³', '0.981', '0.397 (unbound — Whalen+04 / Alvarez+06)'],
        ['1×10⁶ M☉', '19.0 km/s', '1.6 cm⁻³ (8,000× drop)', '0.88', '0.78 (transition — Kitayama+04)'],
        ['5×10⁶ M☉', '12.5 km/s', 'trapped core', '0.12', '0.99 (gravitationally bound)'],
      ]);
    } else if (k === 'schaerer') {
      const D = await loadData('lm_feedback_popiii.json');
      const sc = D.schaerer, rt = D.raiter;
      const g = h('div', 'qd-grid2', stage);
      const b1 = h('div', 'qd-panel', g), b2 = h('div', 'qd-panel', g);
      const p1 = new LinePlot(b1, { width: 520, height: 340, selectMode: 'x' });
      const p2 = new LinePlot(b2, { width: 520, height: 340, selectMode: 'x' });
      plots = [p1, p2];
      p1.setAxes({ xlabel: 'Pop III ZAMS Mass M_* [M☉]', ylabel: 'T_eff [kK]', xlog: true, title: 'Schaerer (2002) Pop III ZAMS T_eff (5–500 M☉)', legendPos: 'br' });
      p1.setSeries([
        { x: sc.M_Msun, y: sc.Teff.map((v) => v / 1e3), label: 'Titan.jl 4D shooting (2.6% RMS)', color: '#0072B2', points: true, width: 2.2 },
        { x: sc.M_Msun, y: sc.Teff_S02.map((v) => v / 1e3), label: 'Schaerer (2002) Table 3', color: '#D55E00', dash: '5 3', points: true, hollow: true, width: 1.8 },
      ]);
      p2.setAxes({ xlabel: 'Electron Density n_e [cm⁻³]', ylabel: 'Line Ratio / Hβ', xlog: true, ylog: true, title: 'Raiter+10 High-n_e Diagnostics (l-resolved ModelAtom)', legendPos: 'tl' });
      p2.setSeries([
        { x: rt.ne, y: rt.HeI10830_Hb, label: 'He I λ10830 / Hβ (5× collisional rise)', color: '#009E73', points: true, width: 2.2 },
        { x: rt.ne, y: rt.Ha_Hb, label: 'Hα / Hβ', color: '#D55E00', width: 1.8 },
        { x: rt.ne, y: rt.HeII1640_Hb, label: 'He II λ1640 / Hβ', color: '#CC79A7', width: 1.9 },
        { x: rt.ne, y: rt.L_2g_Lsun.map((v, i) => v / Math.max(rt.L_Lya_Lsun[i], 1e-9) * 10), label: '10 × (L_2γ / L_Lyα) (quenched >10⁴)', color: '#0072B2', dash: '5 3', width: 1.9 },
      ]);
      const syncM = (i) => {
        p1.setHighlights([{ x: sc.M_Msun[i], y: sc.Teff[i] / 1e3, color: '#0072B2', label: `${sc.M_Msun[i]} M☉: ${fmt(sc.Teff[i]/1e3, 3)} kK` }]);
        selbar.innerHTML = `<strong>Pop III ZAMS M = ${sc.M_Msun[i]} M☉</strong> ` +
          `<span>T_eff = ${fmt(sc.Teff[i]/1e3, 4)} kK (S02: ${fmt(sc.Teff_S02[i]/1e3, 4)} kK)</span> ` +
          `<span>log Q_H = ${fmt(sc.logQH[i], 4)}</span> <span>log Q_HeII = ${fmt(sc.logQHe2[i], 4)} (S02: ${fmt(sc.logQHe2_S02[i], 4)})</span> ` +
          `<span>EW(1640) = ${fmt(sc.EW_1640_A[i], 3)} Å</span>`;
        setTable(tbl, ['M [M☉]', 'T_eff [kK] (Titan / S02)', 'log L/L☉', 'log Q_H', 'log Q_HeII (Titan / S02)', 'EW(He II 1640) [Å]'],
          sc.M_Msun.map((m, j) => [m, `${fmt(sc.Teff[j]/1e3, 3)} / ${fmt(sc.Teff_S02[j]/1e3, 3)}`, fmt(sc.logL[j], 4), fmt(sc.logQH[j], 4), `${fmt(sc.logQHe2[j], 4)} / ${fmt(sc.logQHe2_S02[j], 4)}`, fmt(sc.EW_1640_A[j], 3)]), i);
      };
      const syncNe = (i) => {
        p2.setHighlights([
          { x: rt.ne[i], y: rt.HeI10830_Hb[i], color: '#009E73', label: `10830/Hβ=${fmt(rt.HeI10830_Hb[i], 3)}` },
          { x: rt.ne[i], y: rt.L_2g_Lsun[i] / Math.max(rt.L_Lya_Lsun[i], 1e-9) * 10, color: '#0072B2', vline: false },
        ]);
        selbar.innerHTML = `<strong>Raiter+10 Nebula n_e = ${fmt(rt.ne[i], 3)} cm⁻³</strong> ` +
          `<span>He I 10830/Hβ = ${fmt(rt.HeI10830_Hb[i], 3)}</span> <span>He II 1640/Hβ = ${fmt(rt.HeII1640_Hb[i], 3)}</span> ` +
          `<span>L_2γ/L_Lyα = ${fmt(rt.L_2g_Lsun[i]/rt.L_Lya_Lsun[i], 3)}</span> <span>EW(1640) = ${fmt(rt.EW_1640_A[i], 3)} Å</span>`;
      };
      p1.onSelect((hit) => syncM(hit.pointIndex));
      p2.onSelect((hit) => syncNe(hit.pointIndex));
      p1.draw(); p2.draw(); syncM(sc.M_Msun.indexOf(200) >= 0 ? sc.M_Msun.indexOf(200) : sc.M_Msun.length - 2);
    } else if (k === 'sms') {
      const D = (await loadData('lm_exotics_pisn.json')).sms;
      const mc = D.mcoll, gr = D.gr;
      const g = h('div', 'qd-grid2', stage);
      const b1 = h('div', 'qd-panel', g), b2 = h('div', 'qd-panel', g);
      const p1 = new LinePlot(b1, { width: 520, height: 340, selectMode: 'x' });
      const p2 = new LinePlot(b2, { width: 520, height: 340, selectMode: 'x' });
      plots = [p1, p2];
      p1.setAxes({ xlabel: 'Accretion Rate Ṁ [M☉/yr]', ylabel: 'GR Collapse Mass M_coll [10⁵ M☉]', xlog: true, ylim: [0, 3.8], title: 'SMS GR Collapse: Full 1PN vs Truncated KEPLER 1PN', legendPos: 'tl' });
      p1.setSeries([
        { x: mc.Mdot_Msunyr, y: mc.Mcoll_KEP_Msun.map((v) => v / 1e5), label: 'Titan.jl Trunc-1PN (κ=0.7445)', color: '#E69F00', points: true, width: 2.2 },
        { x: mc.Mdot_Msunyr, y: mc.Mcoll_W17_Msun.map((v) => v / 1e5), label: 'Woods+17 KEPLER ref', color: '#0072B2', dash: '5 3', points: true, hollow: true, width: 1.8 },
        { x: mc.Mdot_Msunyr, y: mc.Mcoll_Msun.map((v) => v / 1e5), label: 'Titan.jl Full 1PN/TOV (κ=1.1245)', color: '#D55E00', points: true, width: 2.3 },
        { x: mc.Mdot_Msunyr, y: mc.Mcoll_H18_Msun.map((v) => v / 1e5), label: 'Haemmerlé+18 GENEC ref', color: '#CC79A7', dash: '2 3', points: true, hollow: true, width: 1.8 },
      ]);
      const gc = gr.GR_supercrit_0p7Rc, gs = gr.GR_subcrit_1p3Rc, gn = gr.Newt_0p7Rc;
      p2.setAxes({ xlabel: 'Time t / t_dyn', ylabel: 'Central Compression ρ_c / ρ_c,0', ylog: true, ylim: [0.7, 250], title: 'Live 1D RHD+GR Collapse (2.5×10⁵ M☉ SMS Core)', legendPos: 'tl' });
      p2.setSeries([
        { x: gc.t_tdyn, y: gc.rho_rat, label: 'GR (0.72 R_crit): COLLAPSE (141×, −0.163c)', color: '#D55E00', width: 2.3 },
        { x: gs.t_tdyn, y: gs.rho_rat, label: 'GR (1.35 R_crit): STABLE', color: '#0072B2', dash: '5 3', width: 1.9 },
        { x: gn.t_tdyn, y: gn.rho_rat, label: 'Newtonian (0.72 R_crit): STABLE', color: '#009E73', dash: '2 3', width: 1.8 },
      ]);
      const syncMc = (i) => {
        p1.setHighlights([
          { x: mc.Mdot_Msunyr[i], y: mc.Mcoll_KEP_Msun[i] / 1e5, color: '#E69F00', label: `Trunc: ${fmt(mc.Mcoll_KEP_Msun[i]/1e5, 3)}×10⁵` },
          { x: mc.Mdot_Msunyr[i], y: mc.Mcoll_Msun[i] / 1e5, color: '#D55E00', vline: false },
        ]);
        selbar.innerHTML = `<strong>Ṁ = ${mc.Mdot_Msunyr[i]} M☉/yr</strong> ` +
          `<span>Full 1PN = ${fmt(mc.Mcoll_Msun[i]/1e5, 3)}×10⁵ M☉ (H18: ${fmt(mc.Mcoll_H18_Msun[i]/1e5, 3)})</span> ` +
          `<span>Trunc 1PN = ${fmt(mc.Mcoll_KEP_Msun[i]/1e5, 3)}×10⁵ M☉ (W17: ${fmt(mc.Mcoll_W17_Msun[i]/1e5, 3)})</span> ` +
          `<span>T_c = ${fmt(mc.Tc_K[i]/1e8, 3)}×10⁸ K</span> <span>X_c = ${fmt(mc.X_c[i], 2)}</span>`;
        setTable(tbl, ['Ṁ [M☉/yr]', 'Full 1PN (Titan)', 'Haemmerlé+18', 'Trunc 1PN (Titan)', 'Woods+17 KEPLER', 'T_c [10⁸ K]'],
          mc.Mdot_Msunyr.map((md, j) => [md, Math.round(mc.Mcoll_Msun[j]).toLocaleString(), Math.round(mc.Mcoll_H18_Msun[j]).toLocaleString(), Math.round(mc.Mcoll_KEP_Msun[j]).toLocaleString(), Math.round(mc.Mcoll_W17_Msun[j]).toLocaleString(), fmt(mc.Tc_K[j]/1e8, 3)]), i);
      };
      p1.onSelect((hit) => syncMc(hit.pointIndex));
      p2.onSelect((hit) => {
        const i = hit.pointIndex;
        if (i < gc.t_tdyn.length) {
          p2.setHighlights([{ x: gc.t_tdyn[i], y: gc.rho_rat[i], color: '#D55E00', label: `${fmt(gc.rho_rat[i], 3)}×` }]);
          selbar.innerHTML = `<strong>Live 1D RHD+GR (0.72 R_crit)</strong> <span>t = ${fmt(gc.t_tdyn[i], 3)} t_dyn</span> <span>ρ_c/ρ_c0 = ${fmt(gc.rho_rat[i], 4)}×</span> <span>u_min/c = ${fmt(gc.umin_c[i], 3)}</span>`;
        }
      });
      p1.draw(); p2.draw(); syncMc(4);
    } else if (k === 'quasistar') {
      const D = await loadData('lm_exotics_pisn.json');
      const qs = D.quasistar.seq, ds = D.darkstar, sed = D.seds;
      const g = h('div', 'qd-grid2', stage);
      const b1 = h('div', 'qd-panel', g), b2 = h('div', 'qd-panel', g);
      const p1 = new LinePlot(b1, { width: 520, height: 340 });
      const p2 = new LinePlot(b2, { width: 520, height: 340, selectMode: 'x' });
      plots = [p1, p2];
      p1.setAxes({ xlabel: 'Effective Temperature T_eff [K] (hot ← → cool)', ylabel: 'Luminosity L [L☉]', xlog: true, ylog: true, xrev: true, xlim: [3200, 1.4e5], ylim: [1e5, 1e11], title: 'Quasi-Star (LRD) & Supermassive Dark Star HR Locus', legendPos: 'bl' });
      p1.setMarkers([{ x: 4000, label: 'Begelman floor', color: '#888' }]);
      p1.setSeries([
        { x: qs.Teff_K, y: qs.L_Lsun, label: 'Zero-Z Quasi-Stars (10⁴→10⁶ M☉)', color: '#D55E00', points: true, width: 2.2 },
        { x: ds.Teff_K, y: ds.L_Lsun, label: 'Supermassive Dark Stars (10²→10⁶ M☉)', color: '#6f42c1', points: true, width: 2.2 },
        { x: [95500], y: [1.3e6], label: '100 M☉ Pop III ZAMS', color: '#0072B2', line: false, points: true, radius: 5 },
      ]);
      p2.setAxes({ xlabel: 'Rest-Frame Wavelength λ_rest [μm]', ylabel: 'Spectral Luminosity λL_λ [L☉]', xlog: true, ylog: true, xlim: [0.08, 4.5], ylim: [1e6, 1e10], title: 'JWST SED Discriminant (10⁵ M☉ Exotic Objects)', legendPos: 'br' });
      p2.setBands([{ x0: 0.09, x1: 0.45, color: '#56B4E9', opacity: 0.22, label: 'JWST NIRSpec (z=10)' }]);
      p2.setSeries([
        { x: sed.lam_um, y: sed.lLl_QS_1e5, label: '10⁵ M☉ Quasi-Star (4.9 kK, LRD)', color: '#D55E00', width: 2.2 },
        { x: sed.lam_um, y: sed.lLl_SMS_1e5, label: '10⁵ M☉ Accreting SMS (5.2 kK)', color: '#E69F00', dash: '5 3', width: 1.9 },
        { x: sed.lam_um, y: sed.lLl_DS_1e5, label: '10⁵ M☉ Dark Star (17.0 kK)', color: '#6f42c1', dash: '6 2 2 2', width: 2.2 },
        { x: sed.lam_um, y: sed.lLl_ZAMS_1e5, label: '10⁵ M☉ Pop III ZAMS (95.5 kK)', color: '#0072B2', dash: '2 3', width: 1.9 },
      ]);
      p1.onSelect((hit) => {
        const i = hit.pointIndex;
        if (hit.seriesIndex === 0) {
          p1.setHighlights([{ x: qs.Teff_K[i], y: qs.L_Lsun[i], color: '#D55E00', label: `${fmt(qs.Menv_Msun[i], 2)} M☉ QS` }]);
          selbar.innerHTML = `<strong>Zero-Z Quasi-Star M_env = ${qs.Menv_Msun[i].toLocaleString()} M☉ (M_BH = ${qs.Mbh_Msun[i].toLocaleString()} M☉)</strong> ` +
            `<span>T_eff = ${fmt(qs.Teff_K[i], 4)} K</span> <span>R_ph = ${fmt(qs.R_AU[i], 3)} AU</span> <span>L = ${fmt(qs.L_Lsun[i], 3)} L☉</span>`;
        } else if (hit.seriesIndex === 1) {
          p1.setHighlights([{ x: ds.Teff_K[i], y: ds.L_Lsun[i], color: '#6f42c1', label: `${fmt(ds.M_Msun[i], 2)} M☉ DS` }]);
          selbar.innerHTML = `<strong>Supermassive Dark Star M = ${ds.M_Msun[i].toLocaleString()} M☉ (n_eff = ${fmt(ds.n_eff[i], 3)})</strong> ` +
            `<span>T_eff = ${fmt(ds.Teff_K[i], 4)} K</span> <span>R = ${fmt(ds.R_AU[i], 3)} AU</span> <span>T_c = ${fmt(ds.Tc_K[i], 3)} K (fusion-free!)</span>`;
        }
      });
      p2.onSelect((hit) => {
        const i = hit.pointIndex, lam = sed.lam_um[i];
        p2.setHighlights([
          { x: lam, y: sed.lLl_QS_1e5[i], color: '#D55E00' },
          { x: lam, y: sed.lLl_DS_1e5[i], color: '#6f42c1', vline: false },
        ]);
        selbar.innerHTML = `<strong>λ_rest = ${fmt(lam, 3)} μm (λ_obs at z=10 = ${fmt(lam * 11, 3)} μm)</strong> ` +
          `<span>QS = ${fmt(sed.lLl_QS_1e5[i], 3)} L☉</span> <span>DS = ${fmt(sed.lLl_DS_1e5[i], 3)} L☉</span> <span>ZAMS = ${fmt(sed.lLl_ZAMS_1e5[i], 3)} L☉</span>`;
      });
      p1.draw(); p2.draw();
      setTable(tbl, ['10⁵ M☉ Primordial Object', 'T_eff [K]', 'Radius', 'L_bol [L☉]', 'JWST z=10 NIRSpec Signature'], [
        ['Begelman+08/Ball+11 Quasi-Star', '4,886 K', '396 AU', '3.73×10⁹', 'Steep red rest-UV slope (Little Red Dot)'],
        ['Hosokawa+13 Accreting SMS', '5,200 K', '350 AU', '3.50×10⁹', 'Cool red supergiant continuum'],
        ['Freese+16 Dark Star (100 GeV)', '17,020 K', '32.3 AU', '3.65×10⁹', 'Bright rest-frame UV peak (~0.18 μm), fusion-free'],
        ['Pop III ZAMS Cluster (1000×100 M☉)', '95,500 K', '4.1 R☉/star', '1.30×10⁹', 'Ultra-blue Rayleigh–Jeans FUV + He II λ1640'],
      ]);
    } else if (k === 'pisn') {
      const D = (await loadData('lm_exotics_pisn.json')).pisn;
      const pm = D.pair_map, blc = D.breakout_lc, plc = D.ppisn_lc;
      const g = h('div', 'qd-grid2', stage);
      const b1 = h('div', 'qd-panel', g), b2 = h('div', 'qd-panel', g);
      const hm = new HeatmapPlot(b1, { width: 520, height: 340, colormap: 'viridis' });
      const p2 = new LinePlot(b2, { width: 520, height: 340, selectMode: 'x' });
      plots = [hm, p2];
      hm.setAxes({ xlabel: 'log₁₀ ρ [g cm⁻³]', ylabel: 'log₁₀ T [K]', title: 'Heger & Woosley (2002): Exact e± Pair-Instability Valley (Γ₁ < 4/3)' });
      hm.setGrid({ xs: pm.xs, ys: pm.ys, z: pm.gamma1_pairs, zmin: 1.22, zmax: 1.42, contour: 4 / 3 });
      hm.setOverlays([
        { x: [3.2, 4.1, 4.9, 5.7, 6.2], y: [8.6, 8.85, 9.05, 9.28, 9.42], color: '#ffffff', label: '100 M☉ He (PISN)' },
        { x: [3.7, 4.6, 5.5, 6.3], y: [8.55, 8.78, 9.02, 9.22], color: '#ffd43b', dash: [5, 3], label: '50 M☉ He (PPISN)' },
      ]);
      p2.setAxes({ xlabel: 'Time t [days]', ylabel: 'Bolometric Luminosity L [erg s⁻¹]', xlog: true, ylog: true, xlim: [0.05, 120], ylim: [1e41, 5e45], title: 'Kasen+11 PISN Breakout & Woosley+07 PPISN Shell Collision', legendPos: 'tr' });
      p2.setSeries([
        { x: blc.t_day, y: blc.L_ergs, label: '200 M☉ PISN Shock Breakout (peak 0.78 d)', color: '#0072B2', width: 2.1 },
        { x: plc.t_day, y: plc.L_ergs, label: 'PPISN Shell Collision (1.77×10⁴⁴ @ 40.5 d)', color: '#D55E00', width: 2.2 },
      ]);
      hm.onSelect((hit) => {
        hm.setHighlight({ x: hit.x, y: hit.y });
        selbar.innerHTML = `<strong>C/O Core State</strong> <span>log₁₀ ρ = ${fmt(hit.x, 3)}</span> <span>log₁₀ T = ${fmt(hit.y, 3)}</span> ` +
          `<span>Γ₁ (with e± pairs) = ${fmt(hit.z, 4)}</span> <span>${hit.z < 4 / 3 ? '⚡ UNSTABLE (Γ₁ < 4/3)' : 'stable (Γ₁ ≥ 4/3)'}</span>`;
      });
      p2.onSelect((hit) => {
        const arr = hit.seriesIndex === 0 ? blc : plc;
        const i = hit.pointIndex;
        p2.setHighlights([{ x: arr.t_day[i], y: arr.L_ergs[i], color: hit.seriesIndex === 0 ? '#0072B2' : '#D55E00', label: `${fmt(arr.L_ergs[i], 3)} erg/s` }]);
        selbar.innerHTML = `<strong>${hit.seriesIndex === 0 ? 'Kasen+11 PISN Breakout' : 'Woosley+07 PPISN Collision'}</strong> ` +
          `<span>t = ${fmt(arr.t_day[i], 3)} d</span> <span>L = ${fmt(arr.L_ergs[i], 3)} erg/s</span> <span>T_ph = ${fmt(arr.Tph_K[i], 4)} K</span>`;
      });
      hm.draw(); p2.draw();
      setTable(tbl, ['Transient Benchmark', 'Published Target', 'Titan.jl Result'], [
        ['Heger & Woosley (2002) Pair Valley', 'Γ₁ < 4/3 at T ~ 10⁹ K, ρ ~ 10⁴–10⁵ g/cm³', 'Exact StellarEOS: min Γ₁ = 1.221 @ 1.18×10⁹ K (605 unstable pts)'],
        ['Kasen+11 200 M☉ RSG PISN Breakout', 'Luminous UV breakout flash (< 1 d), v ~ 2×10⁴ km/s', 't_peak = 0.78 d, v_max = 2.4×10⁴ km/s (HLLC + LP-FLD)'],
        ['Woosley+07 PPISN Shell Collision', 'SLSN L_peak ~ 10⁴⁴ erg/s, E_rad > 10⁵⁰ erg', 'L_peak = 1.77×10⁴⁴ erg/s @ 40.5 d, E_rad = 4.18×10⁵⁰ erg'],
      ]);
    }
    onChange(S);
  }

  select(bar, { label: 'Landmark suite', options: SUITES, value: S.lmsuite, onChange: (v) => { S.lmsuite = v; render(); } });
  await render();
  return { draw: render };
}
