// state.js - the complete control state of the lab and its URL-hash encoding.
// S uses readable names; the hash uses the short keys of SCHEMA. Only values that differ from the defaults are written,
// so a bare URL is the default experiment and any hash is a permalink, a preset and a regression fixture at once.
// Lists and vectors are joined with '_' (never commas).

export const PANELS = ['micro', 'env', 'shoot', 'branch', 'bench', 'qstar', 'spec'];
const SPECTRA = ['star_t0', 'A220_late_ext', 'A220_late_ext_nlte', 'A220_late_ext10_nlte', 'A220_late_ion', 'A220_late_hotion', 'A220_late_noext', 'A220_corona_ext'];

export function makeDefaults() {
  return {
    panel: 'env',     // which experiment the lab shows
    lrho: -12,        // micro: log10 ρ [g cm⁻³]
    fbh: 0.1,         // M_BH / M
    alpha: 0.8,       // α_MLT
    eps: 1e-3,        // ε_in = r_in / R
    mdot: 0,          // Ṁ [M_sun/yr] (env, shoot)
    root: -1,         // env: which steady root (index in increasing R; −1 = the largest)
    bench: 'SUPC',    // bench: SOD | SUBC | SUPC
    bvar: 'T',        // bench variable: rho | T | Trad | u
    run: 'QSTAR',     // qstar: R2 | S6 | A220 | A440 | QSTAR | E1 | A1
    qvar: 'T',        // qstar variable: T | rho | u | Lconv | Lrad | dr
    ti: 1,            // time slider position in [0, 1] (bench, qstar)
    qhist: 'radii',   // qstar history plot: radii | mdot | teff | tpk | lum
    spec: 'A220_late_ext',  // spec: the spectrum shown (lib/spectra/NAME.jls, exported to spectra.json)
    spec2: 'star_t0', // spec: comparison spectrum, or 'none'
    sview: 'sed',     // spec: sed | ratio | T | xH | Halpha | Hbeta
    sxax: 'lam',      // spec: x axis lam | E
    sbb: 0,           // spec: 1 = overlay blackbodies
  };
}

// [stateKey, hashKey, type, range-or-options]
const SCHEMA = [
  ['panel', 'p', 'enum', PANELS],
  ['lrho', 'lr', 'num', [-20, -4]],
  ['fbh', 'fbh', 'num', [1e-3, 0.5]],
  ['alpha', 'al', 'num', [0.1, 3]],
  ['eps', 'eps', 'num', [1e-5, 0.1]],
  ['mdot', 'md', 'num', [0, 100]],
  ['root', 'rt', 'int', [-1, 10]],
  ['bench', 'b', 'enum', ['SOD', 'SUBC', 'SUPC']],
  ['bvar', 'bv', 'enum', ['rho', 'T', 'Trad', 'u']],
  ['run', 'run', 'enum', ['R2', 'S6', 'A220', 'A440', 'QSTAR', 'E1', 'A1']],
  ['qvar', 'v', 'enum', ['T', 'rho', 'u', 'Lconv', 'Lrad', 'dr']],
  ['ti', 'ti', 'num', [0, 1]],
  ['qhist', 'qh', 'enum', ['radii', 'mdot', 'teff', 'tpk', 'lum']],
  ['spec', 'sp', 'enum', SPECTRA],
  ['spec2', 'sp2', 'enum', ['none', ...SPECTRA]],
  ['sview', 'sv', 'enum', ['sed', 'ratio', 'T', 'xH', 'full', 'Halpha', 'Hbeta']],
  ['sxax', 'sx', 'enum', ['lam', 'E']],
  ['sbb', 'bb', 'int', [0, 1]],
];
export const HASH_KEYS = SCHEMA.map((e) => e[1]);

const num = (v) => String(+Number(v).toPrecision(6));
const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));

/** Hash string (without '#') for state S. */
export function encodeHash(S) {
  const D = makeDefaults();
  const p = [];
  for (const [sk, hk, type] of SCHEMA) {
    const v = S[sk];
    if (JSON.stringify(v) === JSON.stringify(D[sk])) continue;
    p.push(`${hk}=${type === 'enum' ? v : type === 'int' ? String(Math.round(v)) : num(v)}`);
  }
  return p.join('&');
}

/** State from a hash string (unknown or invalid entries are ignored). */
export function decodeHash(hash) {
  const S = makeDefaults();
  for (const kv of (hash || '').replace(/^#/, '').split('&')) {
    const [hk, raw] = kv.split('=');
    const ent = SCHEMA.find((e) => e[1] === hk);
    if (!ent || raw === undefined) continue;
    const [sk, , type, range] = ent;
    if (type === 'enum') { if (range.includes(raw)) S[sk] = raw; continue; }
    const v = type === 'int' ? parseInt(raw, 10) : parseFloat(raw);
    if (Number.isFinite(v)) S[sk] = range ? clamp(v, range[0], range[1]) : v;
  }
  return S;
}
