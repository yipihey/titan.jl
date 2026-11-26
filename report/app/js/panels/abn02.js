// abn02.js - Interactive Linked Explorer for Abel, Bryan & Norman (2002, Science 295, 93) End-to-End Evolution:
// Stage 1 (abn02_red -> abn02_final overlay @ 1e12 -> Core Birth @ 3e20 cm⁻³) -> Stage 2 (3D Disk vs Strict 1D) -> Stage 3 (2.2 Myr Photoevaporation & JWST Lines).
import { loadData, h, fmt } from '../core.js';
import { LinePlot } from '../viz/plot.js';

export async function mountAbn02(host, S, onChange = () => {}) {
  const D = await loadData('lm_abn02.json');
  const selbar = h('div', 'qd-selbar', host);
  selbar.innerHTML = '<span><strong>Click any point</strong> on Stage 1, 2, or 3 below to highlight the corresponding shell, protostellar mass, or photoevaporation epoch across linked graphs.</span>';

  const g = h('div', 'qd-grid2', host);
  const b1 = h('div', 'qd-panel', g), b2 = h('div', 'qd-panel', g);
  const b3 = h('div', 'qd-panel', g), b4 = h('div', 'qd-panel', g);
  const p1 = new LinePlot(b1, { width: 520, height: 320, selectMode: 'x' });
  const p2 = new LinePlot(b2, { width: 520, height: 320, selectMode: 'x' });
  const p3 = new LinePlot(b3, { width: 520, height: 320, selectMode: 'x' });
  const p4 = new LinePlot(b4, { width: 520, height: 320, selectMode: 'x' });

  const red = D.red, fin = D.fin, s12 = D.at_1e12, sc = D.core_birth, pr = D.protostar, pp = D.photoevap_polar, p1d = D.photoevap_1d;

  // (1) Stage 1 Radial Density Collapse n(r)
  p1.setAxes({ xlabel: 'Radius r [pc]', ylabel: 'Hydrogen Density n_H [cm⁻³]', xlog: true, ylog: true, xlim: [1e-8, 30], ylim: [10, 3e21], title: '(1) Stage 1: ABN02 Collapse → Core Birth', legendPos: 'bl' });
  p1.setSeries([
    { x: red.n_r.x.map((v) => 10 ** v), y: red.n_r.y.map((v) => 10 ** v), label: 'ABN02 Fig. 2 red IC (9 Myr)', color: '#888', dash: '2 3', width: 1.7, selectable: false },
    { x: fin.n_r.x.map((v) => 10 ** v), y: fin.n_r.y.map((v) => 10 ** v), label: 'ABN02 Fig. 2 final (Enzo 3D)', color: '#222', dash: '5 3', width: 2.0, selectable: false },
    { x: s12.r_pc, y: s12.nH, label: 'Titan.jl @ n_c = 10¹² (Δ=0.18 dex)', color: '#0072B2', width: 2.2 },
    { x: sc.r_pc, y: sc.nH, label: 'Titan.jl Core Birth (3×10²⁰ cm⁻³)', color: '#D55E00', width: 2.2 },
  ]);

  // (2) Stage 1 Thermal & Chemical Validation T(M_enc)
  p2.setAxes({ xlabel: 'Enclosed Gas Mass M_enc [M☉]', ylabel: 'Temperature T [K]', xlog: true, ylog: true, xlim: [1e-3, 5e3], ylim: [100, 4e4], title: '(2) Stage 1: T(M) vs ABN02 Fig. 2 (Δ=0.20 dex)', legendPos: 'tr' });
  p2.setSeries([
    { x: fin.T_M.x.map((v) => 10 ** v), y: fin.T_M.y.map((v) => 10 ** v), label: 'ABN02 final (optically thin H₂)', color: '#222', dash: '5 3', width: 2.0, selectable: false },
    { x: s12.m_Msun, y: s12.T, label: 'Titan.jl @ 10¹² (ExoMol Sobolev)', color: '#0072B2', width: 2.2 },
    { x: sc.m_Msun, y: sc.T, label: 'Titan.jl Core Birth (shock 12.5 kK)', color: '#D55E00', dash: '6 2 2 2', width: 2.0 },
  ]);

  // (3) Stage 2 Protostellar Radius: 3D Disk vs Strict 1D Spherical
  p3.setAxes({ xlabel: 'Protostellar Mass M_* [M☉]', ylabel: 'Stellar Radius R_* [R☉]', xlog: true, ylog: true, xlim: [0.025, 300], ylim: [0.5, 300], title: '(3) Stage 2: 3D Disk (106 M☉) vs Strict 1D (>220 M☉)', legendPos: 'tl' });
  p3.setSeries([
    { x: pr.M_Msun, y: pr.R_disk_Rsun, label: '3D Disk-Regulated (v_rot/v_Kep=0.4 → ZAMS)', color: '#009E73', width: 2.3 },
    { x: pr.M_Msun, y: pr.R_1d_Rsun, label: 'Strict 1D Spherical (bloated → PISN/BH)', color: '#D55E00', dash: '5 3', width: 2.0 },
    { x: pr.M_Msun, y: pr.Rzams_Rsun, label: 'Pop III ZAMS R_ZAMS', color: '#666', dash: '2 3', width: 1.4, selectable: false },
  ]);
  p3.setMarkers([
    { x: 105.7, label: 'Disk breakout 106 M☉', color: '#009E73' },
    { x: 237.0, label: '1D ZAMS 237 M☉', color: '#D55E00' },
  ]);

  // (4) Stage 3 Full-Lifetime JWST Atomic Line Lightcurves (Log-Time Axis!)
  p4.setAxes({ xlabel: 'Post-Breakout Lifetime t [Myr]', ylabel: 'Emergent Line Luminosity [L☉]', xlog: true, ylog: true, xlim: [4.5e-3, 2.35], ylim: [0.2, 2e6], title: '(4) Stage 3: ABN02 Photoevaporation JWST Lines', legendPos: 'bl' });
  p4.setSeries([
    { x: pp.t_Myr, y: pp.L_Lya, label: 'Lyα (Disk 106 M☉)', color: '#D55E00', width: 2.2 },
    { x: p1d.t_Myr, y: p1d.L_Lya, label: 'Lyα (Strict 1D 200 M☉)', color: '#D55E00', dash: '5 3', width: 1.5 },
    { x: pp.t_Myr, y: pp.L_1640, label: 'He II λ1640 (Disk)', color: '#CC79A7', width: 2.1 },
    { x: pp.t_Myr, y: pp.L_Hb, label: 'Hβ λ4861 (Disk)', color: '#0072B2', width: 1.9 },
    { x: pp.t_Myr, y: pp.L_10830, label: 'He I λ10830 (4-decade He III cliff)', color: '#009E73', width: 2.1 },
  ]);

  // Cross-plot linking!
  const syncShell = (idx) => {
    const i = Math.max(0, Math.min(s12.m_Msun.length - 1, idx));
    p1.setHighlights([
      { x: s12.r_pc[i], y: s12.nH[i], color: '#0072B2', label: `M=${fmt(s12.m_Msun[i], 2)} M☉` },
      { x: sc.r_pc[i], y: sc.nH[i], color: '#D55E00', vline: false },
    ]);
    p2.setHighlights([
      { x: s12.m_Msun[i], y: s12.T[i], color: '#0072B2', label: `T=${fmt(s12.T[i], 3)} K` },
      { x: sc.m_Msun[i], y: sc.T[i], color: '#D55E00', vline: false },
    ]);
    // Also highlight matching mass on Stage 2 if in [0.025, 300] Msun!
    const m = s12.m_Msun[i];
    if (m >= 0.025 && m <= 300) {
      let bi = 0, bd = Infinity;
      for (let j = 0; j < pr.M_Msun.length; j++) { const d = Math.abs(Math.log10(pr.M_Msun[j]) - Math.log10(m)); if (d < bd) { bd = d; bi = j; } }
      p3.setHighlights([{ x: pr.M_Msun[bi], y: pr.R_disk_Rsun[bi], color: '#009E73' }, { x: pr.M_Msun[bi], y: pr.R_1d_Rsun[bi], color: '#D55E00', vline: false }]);
    }
    selbar.innerHTML = `<strong>Stage 1 Shell M_enc = ${fmt(m, 3)} M☉</strong> ` +
      `<span>@10¹²: r = ${fmt(s12.r_pc[i], 3)} pc, n_H = ${fmt(s12.nH[i], 3)} cm⁻³, T = ${fmt(s12.T[i], 4)} K, f_H₂ = ${fmt(s12.fH2[i], 3)}</span> ` +
      `<span>| Core birth: Ṁ_1D = ${fmt(sc.Mdot_inst_Msunyr[i], 3)} M☉/yr</span>`;
  };
  const syncProto = (idx) => {
    const i = Math.max(0, Math.min(pr.M_Msun.length - 1, idx));
    p3.setHighlights([
      { x: pr.M_Msun[i], y: pr.R_disk_Rsun[i], color: '#009E73', label: `Disk: ${fmt(pr.R_disk_Rsun[i], 3)} R☉` },
      { x: pr.M_Msun[i], y: pr.R_1d_Rsun[i], color: '#D55E00', vline: false },
    ]);
    selbar.innerHTML = `<strong>Stage 2 Protostar M_* = ${fmt(pr.M_Msun[i], 3)} M☉</strong> ` +
      `<span>3D Disk: Ṁ = ${fmt(pr.Mdot_disk[i], 3)}, R = ${fmt(pr.R_disk_Rsun[i], 3)} R☉, T_eff = ${fmt(pr.Teff_disk_K[i]/1e3, 3)} kK</span> ` +
      `<span>| Strict 1D: Ṁ = ${fmt(pr.Mdot_1d[i], 3)}, R = ${fmt(pr.R_1d_Rsun[i], 3)} R☉</span>`;
  };
  const syncEvap = (idx) => {
    const i = Math.max(0, Math.min(pp.t_Myr.length - 1, idx));
    p4.setHighlights([
      { x: pp.t_Myr[i], y: pp.L_Lya[i], color: '#D55E00', label: `t=${fmt(pp.t_Myr[i], 3)} Myr` },
      { x: pp.t_Myr[i], y: pp.L_1640[i], color: '#CC79A7', vline: false },
      { x: pp.t_Myr[i], y: pp.L_10830[i], color: '#009E73', vline: false },
    ]);
    selbar.innerHTML = `<strong>Stage 3 Photoevaporation t = ${fmt(pp.t_Myr[i], 3)} Myr</strong> ` +
      `<span>R_HII = ${fmt(pp.R_HII_pc[i], 3)} pc</span> <span>n_core = ${fmt(pp.nH_core[i], 3)} cm⁻³</span> ` +
      `<span>f_esc = ${fmt(pp.f_esc_H[i], 3)}</span> <span>f_bound = ${fmt(pp.f_bound[i], 3)}</span> ` +
      `<span>He II 1640 = ${fmt(pp.L_1640[i], 3)} L☉</span> <span>He I 10830 = ${fmt(pp.L_10830[i], 3)} L☉</span>`;
  };

  p1.onSelect((hit) => syncShell(hit.pointIndex));
  p2.onSelect((hit) => syncShell(hit.pointIndex));
  p3.onSelect((hit) => syncProto(hit.pointIndex));
  p4.onSelect((hit) => syncEvap(hit.pointIndex));
  p1.draw(); p2.draw(); p3.draw(); p4.draw();
  syncShell(22);

  const pnote = h('p', 'qd-note', host, 'Exact Abel, Bryan & Norman (2002, Science 295, 93) collapse from abn02_red_fig2.csv (M = 7,137 M☉, exact dm_z19_r halo), validated directly against abn02_final_fig2.csv (0.18 dex RMS in n_H, 0.20 dex in T), continued past ABN02’s optical-depth halt to hydrostatic core birth (3×10²⁰ cm⁻³), then evolved through ZAMS and 2.2 Myr live 1D RHD minihalo photoevaporation.');
  onChange(S);
  return { draw: () => { p1.draw(); p2.draw(); p3.draw(); p4.draw(); } };
}
