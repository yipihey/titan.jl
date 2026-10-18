// widgets.js - chapter widgets: <div class="qd-widget" data-widget="NAME" [data-run="…"]></div> plus this module script.
// Each widget starts from the default state (state.js); the "open in the lab" buttons carry their own hash.
import { makeDefaults } from './state.js';
import { mountMicro } from './panels/micro.js';
import { mountEnvelope } from './panels/envelope.js';
import { mountShoot } from './panels/shoot.js';
import { mountBranches } from './panels/branches.js';
import { mountBench } from './panels/bench.js';
import { mountQstar } from './panels/qstar.js';
import { mountSpec } from './panels/spec.js';
import { mountAnalytic } from './panels/analytic.js';

const MOUNT = {
  micro: (el, S) => mountMicro(el, S),
  envelope: (el, S) => mountEnvelope(el, S),
  shoot: (el, S) => { S.mdot = 5; return mountShoot(el, S); },
  branches: (el, S) => mountBranches(el, S),
  bench: (el, S) => mountBench(el, S),
  qstar: (el, S) => mountQstar(el, S, () => {}, { run: el.dataset.run }),
  analytic: (el, S) => mountAnalytic(el, S, () => {}, { view: el.dataset.view }),
  spec: (el, S) => { if (el.dataset.view) S.sview = el.dataset.view; if (el.dataset.compare) S.spec2 = el.dataset.compare; if (el.dataset.xax) S.sxax = el.dataset.xax; return mountSpec(el, S, () => {}, { spec: el.dataset.spec }); },
};
for (const el of document.querySelectorAll('.qd-widget')) {
  const m = MOUNT[el.dataset.widget];
  if (!m) continue;
  if (el.dataset.mounted) continue;
  el.dataset.mounted = '1';
  m(el, makeDefaults()).catch((e) => { el.textContent = 'widget failed: ' + e; console.error(e); });
}
