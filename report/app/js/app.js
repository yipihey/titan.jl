// app.js - the standalone Titan.jl interactive laboratory: grouped tabs for Landmark Atlas,
// Spectral & Line Profile Studio, and the Quasi-Star WASM + TITAN RHD lab.
import { makeDefaults, encodeHash, decodeHash } from './state.js';
import { mountLandmarks } from './panels/landmarks.js';
import { mountAbn02 } from './panels/abn02.js';
import { mountPopIII } from './panels/popiii.js';
import { mountLineSpec } from './panels/linespec.js';
import { mountMicro } from './panels/micro.js';
import { mountEnvelope } from './panels/envelope.js';
import { mountShoot } from './panels/shoot.js';
import { mountBranches } from './panels/branches.js';
import { mountBench } from './panels/bench.js';
import { mountQstar } from './panels/qstar.js';
import { mountSpec } from './panels/spec.js';
import { copyText, flash } from './viz/download.js';

const GROUPS = [
  { name: 'Landmarks & Pop III', panels: ['landmarks', 'abn02', 'popiii', 'linespec'] },
  { name: 'Quasi-Star & Core Lab', panels: ['micro', 'env', 'shoot', 'branch', 'bench', 'qstar', 'spec'] },
];

const LABELS = {
  landmarks: 'Landmark Atlas (Waves 1–4)',
  abn02: 'ABN02 End-to-End',
  popiii: 'Pop III RHD & Breakout',
  linespec: 'Spectra & Line Profiles',
  micro: 'Microphysics (WASM)',
  env: 'Steady envelope',
  shoot: 'Shooting scan',
  branch: 'Branches & fold',
  bench: 'TITAN shocks',
  qstar: 'Quasi-star in TITAN',
  spec: 'Quasi-star Non-LTE',
};

const EXPLAIN = {
  landmarks: 'Quantitative reproductions of 17 landmark astrophysics papers across shocks, primordial collapse, photoevaporation, Pop III spectra, SMS GR collapse, quasi-stars, dark stars, and PISN/PPISN. Click any point on one graph to highlight linked plots.',
  abn02: 'End-to-end evolution of the exact Abel, Bryan & Norman (2002) minihalo from abn02_red (8.3×10⁵ cm⁻³) through abn02_final overlay (10¹²), hydrostatic core birth (3×10²⁰), protostar/ZAMS (3D disk vs 1D), and 2.2 Myr photoevaporation.',
  popiii: 'Live 1D RHD Pop III protostellar accretion, κ-mechanism pulsations, episodic burst cycles, shock detachment, and 35-group Option A H II / He III / LW breakout with emergent atomic lines.',
  linespec: 'Interactive primordial spectral studio: compare Pop III ZAMS + 2γ nebular continuum, Supermassive Dark Stars, and Quasi-Stars (LRDs) at any redshift z, click any line to inspect its velocity-space profile, and link to Raiter+10 density diagnostics.',
  micro: 'Saha ionization of H/He, Rosseland and Planck means, ∇_ad and Γ₁ along temperature at fixed density (WASM Rust core).',
  env: 'A steady quasi-star envelope solved live in the browser via WebAssembly: all roots of the shooting problem and full interior profiles.',
  shoot: 'The shooting function over trial radii; steady stars are its crossings of 1. Raise Ṁ to watch the two roots merge at the fold.',
  branch: 'All steady quasi-star solutions as a function of the infall rate, published and recomputed here.',
  bench: 'The radiating-shock and shock-tube benchmarks of TITAN, run by the Julia port against the original 1995 Fortran code.',
  qstar: 'The quasi-star outer shell in Julia TITAN: hydrostatic tests and adaptive-grid infall runs.',
  spec: 'Emergent LTE & non-LTE continuum + hydrogen line spectra of TITAN quasi-star infall snapshots.',
};

const MOUNT = {
  landmarks: mountLandmarks, abn02: mountAbn02, popiii: mountPopIII, linespec: mountLineSpec,
  micro: mountMicro, env: mountEnvelope, shoot: mountShoot, branch: mountBranches,
  bench: mountBench, qstar: mountQstar, spec: mountSpec,
};

const S = decodeHash(location.hash);
const host = document.getElementById('lab');
const bar = document.getElementById('lab-bar');

for (const grp of GROUPS) {
  const box = document.createElement('span');
  box.className = 'qd-tabgroup';
  const lbl = document.createElement('span');
  lbl.className = 'qd-tabgroup-label';
  lbl.textContent = grp.name;
  box.appendChild(lbl);
  for (const p of grp.panels) {
    const b = document.createElement('button');
    b.type = 'button';
    b.className = 'qd-btn' + (p === S.panel ? ' is-active' : '');
    b.textContent = LABELS[p];
    b.addEventListener('click', () => { const T = makeDefaults(); T.panel = p; location.hash = encodeHash(T); });
    box.appendChild(b);
  }
  bar.appendChild(box);
}

const btn = document.createElement('button');
btn.type = 'button'; btn.className = 'qd-btn'; btn.textContent = 'copy link';
btn.style.marginLeft = 'auto';
btn.title = 'Copy permalink encoding this exact experiment and controls';
btn.addEventListener('click', async () => { await copyText(location.href); flash(btn, 'copied'); });
bar.appendChild(btn);

document.getElementById('lab-explain').textContent = EXPLAIN[S.panel] || '';
window.addEventListener('hashchange', () => { if (location.hash.replace(/^#/, '') !== encodeHash(S)) location.reload(); });
MOUNT[S.panel](host, S, (S) => { history.replaceState(null, '', '#' + encodeHash(S)); })
  .catch((e) => { host.textContent = 'failed: ' + e; console.error(e); });
