// state.js - the complete control state of the Titan.jl laboratory and its URL-hash encoding.
// S uses readable names; the hash uses the short keys of SCHEMA. Only values that differ from the defaults are written,
// so a bare URL is the default experiment and any hash is a permalink, a preset and a regression fixture at once.
// Lists and vectors are joined with '_' (never commas).

export const PANELS = ['landmarks', 'abn02', 'popiii', 'linespec', 'micro', 'env', 'shoot', 'branch', 'bench', 'qstar', 'spec'];
const SPECTRA = ['star_t0', 'A220_late_ext', 'W30_late_ext', 'A220_late_ext_nlte', 'W30_late_ext_nlte', 'A220_late_ext10_nlte', 'A220_late_ion', 'A220_late_hotion', 'A220_late_noext', 'A220_corona_ext'];
const LM_SUITES = ['ensman', 'omukai', 'whalen', 'schaerer', 'sms', 'quasistar', 'pisn'];
const POPIII_RUNS = ['hoy12_md1e1_full', 'hoy12_md3e1_full', 'episodic_cycle', 'episodic_short', 'hoy12_md3e2', 'hii_breakout'];
const LINE_IDS = ['Lya', 'He1640', 'Hb', 'Ha', 'He10830', 'Paa', 'H2_10S1'];

export function makeDefaults() {
  return {
    panel: 'landmarks', // default panel: Landmark Reproductions Explorer
    lmsuite: 'ensman',  // landmarks: ensman | omukai | whalen | schaerer | sms | quasistar | pisn
    prun: 'hoy12_md1e1_full', // popiii: Pop III RHD run or hii_breakout
    lsline: 'Lya',      // linespec: highlighted spectral line
    lsz: 10,            // linespec: redshift z
    lrho: -12,          // micro: log10 ρ [g cm⁻³]
    fbh: 0.1,           // M_BH / M
    alpha: 0.8,         // α_MLT
    eps: 1e-3,          // ε_in = r_in / R
    mdot: 0,            // Ṁ [M_sun/yr] (env, shoot)
    root: -1,           // env: which steady root (index in increasing R; −1 = the largest)
    bench: 'SUPC',      // bench: SOD | SUBC | SUPC
    bvar: 'T',          // bench variable: rho | T | Trad | u
    run: 'QSTAR',       // qstar: R2 | S6 | A220 | A440 | W30 | QSTAR | E1 | A1
    qvar: 'T',          // qstar variable: T | rho | u | Lconv | Lrad | dr
    ti: 1,              // time slider position in [0, 1] (bench, qstar)
    qhist: 'radii',     // qstar history plot: radii | mdot | teff | tpk | lum
    spec: 'A220_late_ext',
    spec2: 'star_t0',
    sview: 'sed',
    sxax: 'lam',
    sbb: 0,
    aview: 'x',
    amdot: 3,
    afuv: 1,
    atrap: 1,
    alog: 0,
    tview: 'mc',
    tmcall: 0,
    tbal: 'n40_b12',
    tline: 'Halpha',
    tcool: 'sed',
  };
}

// [stateKey, hashKey, type, range-or-options]
const SCHEMA = [
  ['panel', 'p', 'enum', PANELS],
  ['lmsuite', 'lms', 'enum', LM_SUITES],
  ['prun', 'pr', 'enum', POPIII_RUNS],
  ['lsline', 'lsl', 'enum', LINE_IDS],
  ['lsz', 'lsz', 'num', [0, 20]],
  ['lrho', 'lr', 'num', [-20, -4]],
  ['fbh', 'fbh', 'num', [1e-3, 0.5]],
  ['alpha', 'al', 'num', [0.1, 3]],
  ['eps', 'eps', 'num', [1e-5, 0.1]],
  ['mdot', 'md', 'num', [0, 100]],
  ['root', 'rt', 'int', [-1, 10]],
  ['bench', 'b', 'enum', ['SOD', 'SUBC', 'SUPC']],
  ['bvar', 'bv', 'enum', ['rho', 'T', 'Trad', 'u']],
  ['run', 'run', 'enum', ['R2', 'S6', 'A220', 'A440', 'W30', 'QSTAR', 'E1', 'A1']],
  ['qvar', 'v', 'enum', ['T', 'rho', 'u', 'Lconv', 'Lrad', 'dr']],
  ['ti', 'ti', 'num', [0, 1]],
  ['qhist', 'qh', 'enum', ['radii', 'mdot', 'teff', 'tpk', 'lum']],
  ['spec', 'sp', 'enum', SPECTRA],
  ['spec2', 'sp2', 'enum', ['none', ...SPECTRA]],
  ['sview', 'sv', 'enum', ['sed', 'ratio', 'T', 'xH', 'full', 'Halpha', 'Hbeta']],
  ['sxax', 'sx', 'enum', ['lam', 'E']],
  ['sbb', 'bb', 'int', [0, 1]],
  ['aview', 'av', 'enum', ['x', 'dL', 'SB', 'fuv']],
  ['amdot', 'amd', 'num', [0.1, 30]],
  ['afuv', 'afu', 'num', [1e-3, 10]],
  ['atrap', 'atr', 'int', [0, 1]],
  ['alog', 'alg', 'int', [0, 1]],
  ['tview', 'tv', 'enum', ['mc', 'mcincl', 'balmer', 'axisym', 'he', 'cool']],
  ['tmcall', 'tma', 'int', [0, 1]],
  ['tbal', 'tb', 'enum', ['n40_b12', 'n40_lte', 'n80_b12']],
  ['tline', 'tl', 'enum', ['Halpha', 'Hbeta']],
  ['tcool', 'tc', 'enum', ['sed', 'x']],
];
export const HASH_KEYS = SCHEMA.map((e) => e[1]);

const num = (v) => String(+Number(v).toPrecision(6));
const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));

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
