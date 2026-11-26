// linespec.js - Linked Primordial Spectral Features & Atomic/Molecular Line Profile Explorer.
// Click any spectral line or wavelength on the SED (or click a line pill) to zoom the linked
// high-resolution velocity-space line profile (Doppler + Voigt/Stark + P-Cygni/double-peaked Lyα)
// and highlight its density diagnostics from Raiter+10 / ModelAtom!
import { loadData, h, fmt } from '../core.js';
import { LinePlot } from '../viz/plot.js';
import { select, slider } from '../viz/ui.js';

const hP = 6.62607015e-27, kB = 1.380649e-16, cc = 2.99792458e10, SIGMA = 5.670374419e-5;

const LINES = [
  { id: 'Lya',    name: 'Lyα 1216 Å',       lam_um: 0.121567, species: 'H I (2p→1s)',   kind: 'resonant', color: '#D55E00', ew0: 260, m_amu: 1.0 },
  { id: 'He1640', name: 'He II λ1640',      lam_um: 0.164042, species: 'He II (3→2)',   kind: 'recomb',   color: '#CC79A7', ew0: 42,  m_amu: 4.0 },
  { id: 'Hb',     name: 'Hβ λ4861',         lam_um: 0.486133, species: 'H I (4→2)',     kind: 'balmer',   color: '#0072B2', ew0: 145, m_amu: 1.0 },
  { id: 'Ha',     name: 'Hα λ6563',         lam_um: 0.656280, species: 'H I (3→2)',     kind: 'balmer',   color: '#e03131', ew0: 410, m_amu: 1.0 },
  { id: 'He10830',name: 'He I λ10830',      lam_um: 1.08303,  species: 'He I (2³P→2³S)',kind: 'triplet',  color: '#009E73', ew0: 95,  m_amu: 4.0 },
  { id: 'Paa',    name: 'Paα 1.875 μm',     lam_um: 1.87510,  species: 'H I (4→3)',     kind: 'recomb',   color: '#56B4E9', ew0: 55,  m_amu: 1.0 },
  { id: 'H2_10S1',name: 'H₂ 1–0 S(1) 2.12 μm', lam_um: 2.1218, species: 'H₂ ro-vib',   kind: 'mol',      color: '#E69F00', ew0: 18,  m_amu: 2.0 },
];

// Fast pseudo-Voigt profile (Thompson-Cox-Hastings) for thermal Doppler + Stark/radiative damping
function voigt(x, a) {
  const fG = 2.35482, fL = 2 * a;
  const f = Math.pow(fG ** 5 + 2.69269 * fG ** 4 * fL + 2.42843 * fG ** 3 * fL ** 2 + 4.47163 * fG ** 2 * fL ** 3 + 0.07842 * fG * fL ** 4 + fL ** 5, 0.2);
  const eta = 1.36603 * (fL / f) - 0.47719 * (fL / f) ** 2 + 0.11116 * (fL / f) ** 3;
  const g = Math.exp(-x * x), l = 1 / (1 + (x / Math.max(a, 1e-4)) ** 2);
  return (1 - eta) * g + eta * l;
}

function planckLam(lam_um, Teff, L_Lsun) {
  const lam = lam_um * 1e-4;
  const x = hP * cc / (lam * kB * Teff);
  if (x > 250) return 1e-30;
  const Blam = 2 * hP * cc * cc / (lam ** 5) / Math.expm1(x);
  return L_Lsun * (Math.PI * lam * Blam / (SIGMA * Teff ** 4));
}

export async function mountLineSpec(host, S, onChange = () => {}) {
  const [fb, ex] = await Promise.all([loadData('lm_feedback_popiii.json'), loadData('lm_exotics_pisn.json')]);
  const sc = fb.schaerer, rt = fb.raiter, seds = ex.seds;

  let curLine = LINES.find((l) => l.id === S.lsline) || LINES[0];
  let z = S.lsz ?? 10;
  let Tgas = 2.0e4, vBulk = 25.0, logNe = 3.0, logN = 18.5;

  const bar = h('div', 'qd-controls', host);
  const bar2 = h('div', 'qd-controls', host);
  const pillRow = h('div', 'qd-controls', host);
  const selbar = h('div', 'qd-selbar', host);

  const topBox = h('div', 'qd-panel', host);
  const pSed = new LinePlot(topBox, { width: 920, height: 320, selectMode: 'x' });
  const botGrid = h('div', 'qd-grid2', host);
  const bProf = h('div', 'qd-panel', botGrid), bDiag = h('div', 'qd-panel', botGrid);
  const pProf = new LinePlot(bProf, { width: 450, height: 310, selectMode: 'x' });
  const pDiag = new LinePlot(bDiag, { width: 450, height: 310, selectMode: 'x' });

  const pillBtns = new Map();
  h('span', 'qd-label', pillRow, 'Highlight line:');
  for (const ln of LINES) {
    const b = h('button', 'qd-pill', pillRow, ln.name);
    b.type = 'button';
    b.addEventListener('click', () => { curLine = ln; S.lsline = ln.id; drawAll(); });
    pillBtns.set(ln.id, b);
  }

  function nearestNeIdx() {
    const ne = 10 ** logNe;
    let bi = 0, bd = Infinity;
    for (let i = 0; i < rt.ne.length; i++) {
      const d = Math.abs(Math.log10(rt.ne[i]) - logNe);
      if (d < bd) { bd = d; bi = i; }
    }
    return bi;
  }

  function drawAll() {
    for (const [id, b] of pillBtns) b.classList.toggle('is-active', id === curLine.id);
    const iNe = nearestNeIdx();
    const ne = 10 ** logNe;
    const f2g = rt.L_2g_Lsun[iNe] / Math.max(rt.L_2g_Lsun[0], 1e-9);
    const fac10830 = rt.HeI10830_Hb[iNe] / Math.max(rt.HeI10830_Hb[0], 1e-9);

    // 1. Composite SED (rest-frame or observed at redshift z) with 2γ nebular continuum + Gaussian emission spikes!
    const lamRest = [];
    for (let i = 0; i < 320; i++) lamRest.push(0.08 * Math.pow(4.5 / 0.08, i / 319));
    // Insert fine sampling right around each line so emission lines render crisply
    for (const ln of LINES) {
      for (const dv of [-0.015, -0.004, 0, 0.004, 0.015]) lamRest.push(ln.lam_um * (1 + dv));
    }
    lamRest.sort((a, b) => a - b);

    const scaleX = (l) => (z > 0 ? l * (1 + z) : l);
    const xGrid = lamRest.map(scaleX);

    // 10^5 Msun Pop III ZAMS + nebular 2γ + recombination lines
    const yPop3 = lamRest.map((l) => {
      let y = planckLam(l, 95500, 1.3e9);
      // H 2-photon continuum peaks at ~0.145 μm (quenched smoothly at high n_e per Raiter+10)
      if (l > 0.1216) {
        const y2g = f2g * 3.8e8 * Math.exp(-0.5 * ((Math.log(l / 0.16)) / 0.55) ** 2);
        y += y2g;
      }
      for (const ln of LINES) {
        if (ln.kind === 'mol') continue;
        const sig = ln.lam_um * 0.0035;
        const amp = (ln.id === 'He10830' ? fac10830 : 1.0) * (ln.ew0 / 120) * 6.5e8;
        y += amp * Math.exp(-0.5 * ((l - ln.lam_um) / sig) ** 2);
      }
      return y;
    });
    const yQS = lamRest.map((l) => planckLam(l, 4886, 3.73e9));
    const yDS = lamRest.map((l) => planckLam(l, 17020, 3.65e9));

    const xLab = z > 0 ? `Observed Wavelength λ_obs [μm] (z = ${fmt(z, 3)})` : 'Rest-Frame Wavelength λ_rest [μm]';
    pSed.setAxes({ xlabel: xLab, ylabel: 'Spectral Luminosity λL_λ [L☉]', xlog: true, ylog: true, xlim: [scaleX(0.085), scaleX(3.5)], ylim: [1e6, 8e9], title: 'Composite Primordial SEDs + Nebular 2γ & Atomic Emission Lines', legendPos: 'br' });
    pSed.setBands(z > 0 ? [
      { x0: 0.6, x1: 5.3, color: '#56B4E9', opacity: 0.18, label: 'JWST NIRSpec (0.6–5.3 μm)' },
      { x0: 5.3, x1: 28.0, color: '#E69F00', opacity: 0.13, label: 'JWST MIRI' },
    ] : [{ x0: 0.09, x1: 0.45, color: '#56B4E9', opacity: 0.18, label: 'Rest-UV at z=10' }]);
    pSed.setMarkers(LINES.map((ln) => ({ x: scaleX(ln.lam_um), label: ln.id === curLine.id ? ln.name : '', color: ln.id === curLine.id ? ln.color : 'rgba(130,130,130,0.45)', dash: ln.id === curLine.id ? '4 2' : '2 4' })));
    pSed.setSeries([
      { x: xGrid, y: yPop3, label: `Pop III ZAMS + Nebula (n_e=${fmt(ne, 2)})`, color: '#0072B2', width: 2.0 },
      { x: xGrid, y: yDS, label: '10⁵ M☉ Supermassive Dark Star (17.0 kK)', color: '#6f42c1', dash: '6 2 2 2', width: 2.0 },
      { x: xGrid, y: yQS, label: '10⁵ M☉ Quasi-Star (4.9 kK, LRD)', color: '#D55E00', dash: '5 3', width: 2.0 },
    ]);
    pSed.draw();
    const xCur = scaleX(curLine.lam_um);
    pSed.setHighlights([{ x: xCur, y: 1.5e9, color: curLine.color, label: curLine.name }]);

    // 2. High-resolution velocity-space line profile (v in [-450, +450] km/s)
    const vTh = Math.sqrt(2 * kB * Tgas / (curLine.m_amu * 1.6726e-24)) / 1e5; // km/s
    const vTurb = 12.0;
    const bD = Math.sqrt(vTh * vTh + vTurb * vTurb);
    // Stark + natural damping parameter a (scales with n_e for Balmer/He lines)
    const aDamp = (curLine.kind === 'balmer' ? 0.03 + 0.15 * Math.pow(ne / 1e6, 0.67) : curLine.kind === 'resonant' ? 0.04 : 0.015);
    const vs = [], fEm = [], fAbs = [];
    const tau0 = Math.pow(10, logN - 17.0) * (curLine.kind === 'resonant' ? 30.0 : curLine.kind === 'triplet' ? 3.0 : 0.4);
    for (let i = 0; i <= 240; i++) {
      const v = -450 + 900 * (i / 240);
      vs.push(v);
      let prof;
      if (curLine.kind === 'resonant') {
        // Neufeld / Verhamme double-peaked resonant Lyα scattering profile with red/blue asymmetry from vBulk!
        const vPeak = Math.min(260, 28 * Math.pow(Math.max(tau0, 1), 0.33));
        const red = voigt((v - (vPeak + 0.6 * vBulk)) / (bD * 1.6), aDamp * 2);
        const blue = Math.exp(-0.012 * Math.max(vBulk, 0)) * voigt((v + (vPeak - 0.5 * vBulk)) / (bD * 1.5), aDamp * 2);
        prof = 1.0 + 3.2 * (red + 0.72 * blue);
      } else if (curLine.kind === 'triplet' || curLine.kind === 'balmer') {
        // Recombination emission + P-Cygni blue-shifted absorption trough when expanding (vBulk > 0) or inverse P-Cygni when infalling (vBulk < 0)!
        const em = voigt((v - 0.25 * vBulk) / bD, aDamp);
        const tauAbs = Math.min(3.5, 0.25 * Math.sqrt(tau0)) * voigt((v + vBulk) / (bD * 0.9), aDamp);
        prof = (1.0 + (curLine.id === 'He10830' ? 2.2 * fac10830 : 2.6) * em) * Math.exp(-0.55 * tauAbs);
      } else {
        prof = 1.0 + 2.4 * voigt((v - vBulk) / bD, aDamp);
      }
      fEm.push(prof);
      fAbs.push(1.0 + 2.4 * voigt(v / bD, aDamp));
    }
    pProf.setAxes({ xlabel: 'Velocity Offset Δv [km s⁻¹]', ylabel: 'Normalized Flux F_λ / F_cont', xlim: [-420, 420], title: `Line Profile: ${curLine.name} (${curLine.species})`, legendPos: 'tr' });
    pProf.setMarkers([{ x: 0, label: 'v = 0', color: '#888' }]);
    pProf.setHMarkers([{ y: 1.0, label: 'continuum', color: '#888', dash: '2 3' }]);
    pProf.setSeries([
      { x: vs, y: fEm, label: `Dynamic profile (v=${fmt(vBulk, 3)} km/s, b=${fmt(bD, 3)})`, color: curLine.color, width: 2.3 },
      { x: vs, y: fAbs, label: 'Pure static Voigt reference', color: '#888', dash: '5 3', width: 1.4, selectable: false },
    ]);
    pProf.draw();

    // 3. Linked Raiter+10 Density Diagnostics Panel
    pDiag.setAxes({ xlabel: 'Electron Density n_e [cm⁻³]', ylabel: 'Line Ratio / Hβ & 2γ Factor', xlog: true, ylog: true, title: 'Raiter+10 Density Diagnostics (Click n_e)', legendPos: 'tl' });
    pDiag.setSeries([
      { x: rt.ne, y: rt.HeI10830_Hb, label: 'He I λ10830 / Hβ', color: '#009E73', points: true, width: 2.1 },
      { x: rt.ne, y: rt.HeII1640_Hb, label: 'He II λ1640 / Hβ', color: '#CC79A7', width: 1.9 },
      { x: rt.ne, y: rt.Ha_Hb, label: 'Hα / Hβ', color: '#D55E00', width: 1.8 },
      { x: rt.ne, y: rt.L_2g_Lsun.map((v) => v / Math.max(rt.L_2g_Lsun[0], 1e-9)), label: '2γ continuum survival f_2γ', color: '#0072B2', dash: '5 3', width: 1.9 },
    ]);
    pDiag.draw();
    pDiag.setHighlights([
      { x: rt.ne[iNe], y: rt.HeI10830_Hb[iNe], color: '#009E73', label: `n_e=${fmt(rt.ne[iNe], 2)}` },
      { x: rt.ne[iNe], y: f2g, color: '#0072B2', vline: false },
    ]);

    selbar.innerHTML = `<strong>Selected Line: ${curLine.name}</strong> ` +
      `<span>λ_rest = ${curLine.lam_um} μm</span> <span>λ_obs(z=${fmt(z, 3)}) = ${fmt(xCur, 4)} μm</span> ` +
      `<span>b_D = ${fmt(bD, 3)} km/s</span> <span>n_e = ${fmt(ne, 2)} cm⁻³</span> ` +
      `<span>2γ survival = ${fmt(f2g * 100, 3)}%</span> <span>He I 10830/Hβ = ${fmt(rt.HeI10830_Hb[iNe], 3)}</span>`;
    onChange(S);
  }

  pSed.onSelect((hit) => {
    const lamR = z > 0 ? hit.x / (1 + z) : hit.x;
    let best = LINES[0], bd = Infinity;
    for (const ln of LINES) { const d = Math.abs(Math.log10(ln.lam_um) - Math.log10(lamR)); if (d < bd) { bd = d; best = ln; } }
    curLine = best; S.lsline = best.id; drawAll();
  });
  pDiag.onSelect((hit) => {
    logNe = Math.log10(hit.x);
    sNe.set(logNe);
    drawAll();
  });

  slider(bar, { label: 'Redshift z', min: 0, max: 20, step: 0.5, value: z, onInput: (v) => { z = v; S.lsz = v; drawAll(); } });
  const sNe = slider(bar, { label: 'log₁₀ n_e [cm⁻³]', min: 2, max: 8, step: 0.25, value: logNe, onInput: (v) => { logNe = v; drawAll(); } });
  slider(bar2, { label: 'T_gas [K]', min: 3000, max: 50000, step: 1000, value: Tgas, onInput: (v) => { Tgas = v; drawAll(); } });
  slider(bar2, { label: 'Bulk v [km/s] (+out / −in)', min: -120, max: 120, step: 5, value: vBulk, onInput: (v) => { vBulk = v; drawAll(); } });
  slider(bar2, { label: 'log₁₀ N_col [cm⁻²]', min: 16, max: 21, step: 0.25, value: logN, onInput: (v) => { logN = v; drawAll(); } });

  drawAll();
  return { draw: drawAll };
}
