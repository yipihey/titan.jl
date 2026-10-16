import fs from 'fs';
const { instance } = await WebAssembly.instantiate(fs.readFileSync('out/probe.wasm'));
const E = instance.exports, nT = 941, nR = 521;
const buf = new Float64Array(fs.readFileSync('out/kap64.bin').buffer.slice(0));
E.table_init(nT, nR, 3.3, 0.0050000000000003375, -12.0, 0.025000000000000355);
for (let i = 0; i < nT; i++) for (let j = 0; j < nR; j++) E.table_set(i + 1, j + 1, buf[i * nR + j], buf[nT * nR + i * nR + j]);
const rows = [[8.26785524971e15, 8639.089107546477], [8.26785524972e15, 8639.089107535654], [8.26785524973e15, 8639.089107524862], [8.26785524974e15, 8639.089107514072], [8e15, 8936.425607290605], [3e15, 21114.728502780086]];
for (const [R, Tb] of rows) { const Pb = E.bc_Pb(R, 5); console.log(R, Pb, E.integd(R, Pb, Tb)); }
