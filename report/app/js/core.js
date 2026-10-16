// core.js - the WebAssembly core (Rust crate qs_core) with its opacity table, loaded once per page.
import init, * as qs from '../pkg/qs_core.js';

let ready = null;
/** Resolves to the module's exports once the WebAssembly and the opacity table are loaded. */
export function loadCore() {
  if (!ready) {
    ready = (async () => {
      await init();
      const res = await fetch(new URL('../data/kaptable.bin', import.meta.url));
      if (!res.ok) throw new Error('opacity table: HTTP ' + res.status);
      qs.load_table(new Uint8Array(await res.arrayBuffer()));
      return qs;
    })();
  }
  return ready;
}

const cache = new Map();
/** JSON data file from app/data (cached). */
export function loadData(name) {
  if (!cache.has(name)) {
    cache.set(name, fetch(new URL(`../data/${name}`, import.meta.url)).then((r) => {
      if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`);
      return r.json();
    }));
  }
  return cache.get(name);
}

/** Model constants shared by the panels. */
export const MODEL = { M: 1e6, logRmin: 14.0, logRmax: 16.6, dlogR: 0.025 };

/** Let the browser paint a status message before a long synchronous computation. */
export const yieldToPaint = () => new Promise((r) => requestAnimationFrame(() => setTimeout(r, 0)));

/** Element helper. */
export function h(tag, cls, parent, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  if (parent) parent.appendChild(e);
  return e;
}

/** Caption with a back-link to the section of the book that explains the panel. */
export function caption(parent, text, href, label) {
  const p = h('p', 'qd-note', parent, text + ' ');
  if (href) { const a = h('a', null, p, label || '§ book'); a.href = href; }
  return p;
}

/** Number formatting for read-outs. */
export const fmt = (v, d = 3) => (!Number.isFinite(v) ? '–' : (Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-2 && v !== 0)) ? v.toExponential(d - 1) : String(+v.toPrecision(d)));
