// app.js - the standalone laboratory: one experiment at a time, all state in the URL hash, a copy-link button.
import { makeDefaults, encodeHash, decodeHash, PANELS } from './state.js';
import { mountMicro } from './panels/micro.js';
import { mountEnvelope } from './panels/envelope.js';
import { mountShoot } from './panels/shoot.js';
import { mountBranches } from './panels/branches.js';
import { mountBench } from './panels/bench.js';
import { mountQstar } from './panels/qstar.js';
import { copyText, flash } from './viz/download.js';

const LABELS = {
  micro: 'Microphysics', env: 'Steady envelope', shoot: 'Shooting scan', branch: 'Branches and fold',
  bench: 'TITAN benchmarks', qstar: 'Quasi-star in TITAN',
};
const EXPLAIN = {
  micro: 'Saha ionization of H/He, the Rosseland and Planck means, ∇_ad and Γ₁ along temperature at fixed density (book §2).',
  env: 'A steady quasi-star envelope solved in the browser: all roots of the shooting problem, then the profile of one (book §3).',
  shoot: 'The shooting function over trial radii; steady stars are its crossings of 1. Raise Ṁ to watch the two roots merge at the fold (book §4).',
  branch: 'All steady solutions as a function of the infall rate, published and recomputed here (book §4).',
  bench: 'The radiating-shock and shock-tube benchmarks of TITAN, run by the Julia port, against the original Fortran (book §6).',
  qstar: 'The quasi-star outer shell in the Julia TITAN: the hydrostatic test and the first infall runs (book §7).',
};
const MOUNT = { micro: mountMicro, env: mountEnvelope, shoot: mountShoot, branch: mountBranches, bench: mountBench, qstar: mountQstar };

const S = decodeHash(location.hash);
const host = document.getElementById('lab');
const bar = document.getElementById('lab-bar');
const tabs = document.createElement('span'); tabs.className = 'qd-row'; bar.appendChild(tabs);
for (const p of PANELS) {
  const b = document.createElement('button'); b.type = 'button'; b.className = 'qd-btn' + (p === S.panel ? ' is-active' : '');
  b.textContent = LABELS[p];
  b.addEventListener('click', () => { const T = makeDefaults(); T.panel = p; location.hash = encodeHash(T); });
  tabs.appendChild(b);
}
const btn = document.createElement('button'); btn.type = 'button'; btn.className = 'qd-btn'; btn.textContent = 'copy link';
btn.title = 'Copy the address of this exact experiment (the part after # encodes every control)';
btn.addEventListener('click', async () => { await copyText(location.href); flash(btn, 'copied'); });
bar.appendChild(btn);
document.getElementById('lab-explain').textContent = EXPLAIN[S.panel];
window.addEventListener('hashchange', () => { if (location.hash.replace(/^#/, '') !== encodeHash(S)) location.reload(); });
MOUNT[S.panel](host, S, (S) => { history.replaceState(null, '', '#' + encodeHash(S)); })
  .catch((e) => { host.textContent = 'failed: ' + e; console.error(e); });
