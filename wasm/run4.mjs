import fs from 'fs';
const { instance } = await WebAssembly.instantiate(fs.readFileSync('out/probe.wasm'));
const E = instance.exports, nT = 941, nR = 521;
const buf = new Float64Array(fs.readFileSync('out/kap64.bin').buffer.slice(0));
E.table_init(nT, nR, 3.3, 0.0050000000000003375, -12.0, 0.025000000000000355);
for (let i = 0; i < nT; i++) for (let j = 0; j < nR; j++) E.table_set(i + 1, j + 1, buf[i * nR + j], buf[nT * nR + i * nR + j]);
const tryf = (n, f) => { try { console.log(n, f()); } catch (e) { console.log(n, 'TRAP', e.message); } };
tryf('thermo_nad', () => E.thermo_nad(10, 7000));
tryf('mlt', () => E.mlt_nabla(1e4, 7.5e4, 7.5e15, 1.9e39));
tryf('Pb', () => E.bc_Pb(8e15, 0));
tryf('tau', () => E.bc_tau(8e15, 1));
tryf('integ', () => E.integ(8e15, E.bc_Pb(8e15, 0), 7300));
tryf('rhs_dr', () => E.rhs_dr(Math.log(8), 8e15, 1.98e39, Math.log(7300)));
tryf('rhs_rho', () => E.rhs_rho(Math.log(8), 8e15, 1.98e39, Math.log(7300)));
tryf('mlt_full', () => E.mlt_full(1e7, 2.6e5, 3.4e15, 1.6e39));
tryf('integd', () => E.integd(8e15, E.bc_Pb(8e15, 0), 7300));
console.log('MARK step/stage/phase/accepted =', [1, 2, 3, 4].map((i) => E.mark(i)));
for (const [f, md] of [[0.1, 0.0], [0.1, 5.0], [0.01, 10.0]]) {
  const t0 = performance.now(); const n = E.all_roots2(f, md, 0.8, 1e-3);
  const r = []; for (let k = 1; k <= n; k++) r.push(E.root2(k));
  console.log(`wasm all_roots2(${f}, ${md}): n=${n}`, r.map(String).join(' '), `${((performance.now() - t0) / 1000).toFixed(2)} s`);
}
