// popiii.js - Linked Explorer for Titan.jl Pop III Protostar RHD Runs (HOY12 Pulsations, Episodic Accretion,
// Shock Detachment) & Live 1D RHD H II / He III / LW Breakout with Atomic Line Lightcurves.
import { loadData, h, fmt } from '../core.js';
import { LinePlot } from '../viz/plot.js';
import { select, slider } from '../viz/ui.js';

export async function mountPopIII(host, S, onChange = () => {}, opts = {}) {
  const D = await loadData('popiii_runs.json');
  const keys = [...Object.keys(D.runs), ...(D.hii_breakout ? ['hii_breakout'] : [])];
  const initRun = opts.run || S.prun || keys[0];
  S.prun = keys.includes(initRun) ? initRun : keys[0];

  const bar = h('div', 'qd-controls', host);
  const selbar = h('div', 'qd-selbar', host);
  const stage = h('div', null, host);
  const cap = h('p', 'qd-note', host);

  let plots = [];
  const clear = () => { plots.forEach((p) => p.destroy && p.destroy()); plots = []; stage.innerHTML = ''; };

  function render() {
    clear();
    const rk = S.prun;
    if (rk === 'hii_breakout') {
      const H = D.hii_breakout;
      cap.textContent = 'Live 1D spherical RHD H II / He III / LW breakout around an accreting 45 M☉ Pop III star (Option A 35-group IonRay + l-resolved ModelAtom). Click any time on either plot to link front radii, infall reversal, and atomic lines.';
      const g = h('div', 'qd-grid2', stage);
      const b1 = h('div', 'qd-panel', g), b2 = h('div', 'qd-panel', g);
      const p1 = new LinePlot(b1, { width: 520, height: 340, selectMode: 'x' });
      const p2 = new LinePlot(b2, { width: 520, height: 340, selectMode: 'x' });
      plots = [p1, p2];
      p1.setAxes({ xlabel: 'Time t [yr]', ylabel: 'Front / Shock Radius [pc]', ylog: true, ylim: [1e-3, 1.0], title: 'Dynamical H II / He III / LW Breakout & Shock', legendPos: 'tl' });
      p1.setSeries([
        { x: H.t_yr, y: H.R_LW_pc, label: 'LW H₂ dissociation front', color: '#E69F00', dash: '5 3', width: 1.9 },
        { x: H.t_yr, y: H.R_HII_pc, label: 'H II ionization front', color: '#D55E00', width: 2.3 },
        { x: H.t_yr, y: H.R_HeIII_pc, label: 'He III double-ionization front', color: '#CC79A7', width: 2.0 },
        { x: H.t_yr, y: H.R_shock_pc, label: 'D-type champagne shock', color: '#0072B2', dash: '2 3', width: 1.9 },
      ]);
      p2.setAxes({ xlabel: 'Time t [yr]', ylabel: 'Atomic Line Luminosity [L☉]', ylog: true, ylim: [10, 2e6], title: 'Emergent Atomic Lines (l-resolved ModelAtom)', legendPos: 'br' });
      p2.setSeries([
        { x: H.t_yr, y: H.L_Lya_Lsun, label: 'Lyα 1216 Å', color: '#D55E00', width: 2.2 },
        { x: H.t_yr, y: H.L_2g_Lsun, label: 'H 2-photon continuum', color: '#E69F00', dash: '5 3', width: 1.7 },
        { x: H.t_yr, y: H.L_HeII1640_Lsun, label: 'He II λ1640', color: '#CC79A7', width: 2.1 },
        { x: H.t_yr, y: H.L_Hb_Lsun, label: 'Hβ λ4861', color: '#0072B2', width: 1.9 },
        { x: H.t_yr, y: H.L_HeI10830_Lsun, label: 'He I λ10830', color: '#009E73', width: 2.0 },
      ]);
      const sync = (i) => {
        const k = Math.max(0, Math.min(H.t_yr.length - 1, i));
        p1.setHighlights([
          { x: H.t_yr[k], y: H.R_HII_pc[k], color: '#D55E00', label: `R_HII=${fmt(H.R_HII_pc[k], 3)} pc` },
          { x: H.t_yr[k], y: H.R_HeIII_pc[k], color: '#CC79A7', vline: false },
        ]);
        p2.setHighlights([
          { x: H.t_yr[k], y: H.L_Lya_Lsun[k], color: '#D55E00', label: `Lyα=${fmt(H.L_Lya_Lsun[k], 3)} L☉` },
          { x: H.t_yr[k], y: H.L_HeII1640_Lsun[k], color: '#CC79A7', vline: false },
        ]);
        selbar.innerHTML = `<strong>t = ${fmt(H.t_yr[k], 4)} yr</strong> ` +
          `<span>T_eff = ${fmt(H.Teff[k]/1e3, 3)} kK</span> <span>R_HII = ${fmt(H.R_HII_pc[k], 3)} pc</span> ` +
          `<span>u_max = ${fmt(H.u_max_kms[k], 3)} km/s</span> <span>He II 1640 = ${fmt(H.L_HeII1640_Lsun[k], 3)} L☉</span> ` +
          `<span>He I 10830 = ${fmt(H.L_HeI10830_Lsun[k], 3)} L☉</span>`;
      };
      p1.onSelect((hit) => sync(hit.pointIndex));
      p2.onSelect((hit) => sync(hit.pointIndex));
      p1.draw(); p2.draw(); sync(Math.floor(H.t_yr.length * 0.75));
      onChange(S);
      return;
    }

    const R = D.runs[rk], ts = R.ts, snaps = R.snaps || [];
    cap.textContent = R.desc + '. Click any time on the left graph to highlight that epoch and load the closest live radial interior profile on the right!';
    const g = h('div', 'qd-grid2', stage);
    const b1 = h('div', 'qd-panel', g), b2 = h('div', 'qd-panel', g);
    const p1 = new LinePlot(b1, { width: 520, height: 340, selectMode: 'x' });
    const p2 = new LinePlot(b2, { width: 520, height: 340, selectMode: 'x' });
    plots = [p1, p2];
    p1.setAxes({ xlabel: 'Time t [yr]', ylabel: 'Radius [R☉]', ylog: true, title: 'Shock (R_*) & Photosphere (R_ph) Evolution', legendPos: 'tl' });
    p1.setSeries([
      { x: ts.t_yr, y: ts.R_star, label: 'Accretion shock R_* [R☉]', color: '#D55E00', width: 2.2 },
      { x: ts.t_yr, y: ts.R_ph, label: 'Photosphere R_ph [R☉]', color: '#0072B2', dash: '5 3', width: 1.9 },
    ]);

    const drawSnap = (tsIdx) => {
      const i = Math.max(0, Math.min(ts.t_yr.length - 1, tsIdx));
      p1.setHighlights([
        { x: ts.t_yr[i], y: ts.R_star[i], color: '#D55E00', label: `R_*=${fmt(ts.R_star[i], 3)} R☉` },
        { x: ts.t_yr[i], y: ts.R_ph[i], color: '#0072B2', vline: false },
      ]);
      const st = ts.step[i];
      if (snaps.length) {
        let best = snaps[0], bd = Infinity;
        for (const sn of snaps) { const d = Math.abs(sn.step - st); if (d < bd) { bd = d; best = sn; } }
        const pr = best.prof;
        p2.setAxes({ xlabel: 'Radius r [R☉]', ylabel: 'Temperature T [K] & Density (scaled)', xlog: true, ylog: true, title: `Radial Interior Profile (step ${best.step})`, legendPos: 'bl' });
        const sers = [];
        if (pr.T) sers.push({ x: pr.r_Rsun, y: pr.T, label: 'Temperature T [K]', color: '#D55E00', width: 2.2 });
        if (pr.rho) sers.push({ x: pr.r_Rsun, y: pr.rho.map((v) => v * 1e12), label: '10¹² × ρ [g/cm³]', color: '#0072B2', dash: '5 3', width: 1.8 });
        if (pr.Gamma_Edd) sers.push({ x: pr.r_Rsun, y: pr.Gamma_Edd.map((v) => Math.max(v, 1e-4) * 1e4), label: '10⁴ × Γ_Edd (super-Edd > 10⁴)', color: '#009E73', dash: '2 3', width: 1.7 });
        p2.setSeries(sers);
        p2.setMarkers([{ x: ts.R_star[i], label: 'R_*', color: '#D55E00' }]);
        p2.draw();
      } else {
        p2.setAxes({ xlabel: 'Time t [yr]', ylabel: 'Luminosity L_ph & L_Edd [L☉]', ylog: true, title: 'Bolometric Luminosity vs Eddington Limit', legendPos: 'tl' });
        p2.setSeries([
          { x: ts.t_yr, y: ts.L_ph, label: 'L_ph [L☉]', color: '#0072B2', width: 2.1 },
          { x: ts.t_yr, y: ts.L_Edd, label: 'L_Edd [L☉]', color: '#D55E00', dash: '5 3', width: 1.8 },
        ]);
        p2.setHighlights([{ x: ts.t_yr[i], y: ts.L_ph[i], color: '#0072B2' }]);
        p2.draw();
      }
      selbar.innerHTML = `<strong>t = ${fmt(ts.t_yr[i], 4)} yr (step ${st})</strong> ` +
        `<span>M_* = ${fmt(ts.M_star[i], 4)} M☉</span> <span>R_* = ${fmt(ts.R_star[i], 4)} R☉</span> ` +
        `<span>T_eff = ${fmt(ts.Teff[i], 4)} K</span> <span>max(dlnρ/dlnr) = ${fmt(ts.max_dlnrho_dlnr[i], 3)}</span> ` +
        `<span>max Γ_Edd = ${fmt(ts.max_Gamma_Edd[i], 3)}</span>`;
    };

    p1.onSelect((hit) => drawSnap(hit.pointIndex));
    p1.draw();
    drawSnap(Math.floor(ts.t_yr.length * 0.8));
    onChange(S);
  }

  const optsList = keys.map((k) => ({
    value: k,
    label: k === 'hii_breakout' ? 'Live 1D RHD H II / He III / LW Breakout (Option A, 45 M☉)' : D.runs[k].desc.slice(0, 72),
  }));
  select(bar, { label: 'Pop III RHD run', options: optsList, value: S.prun, onChange: (v) => { S.prun = v; render(); } });
  render();
  return { draw: render };
}
