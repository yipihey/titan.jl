
"use strict";
const $ = id => document.getElementById(id);
const tok = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const RUNS = {};
const ST = { run: "base64", i: 0, cursor: null, playing: false, timer: null, xlim: { logr: [null, null], logm: [null, null] } };
// run key -> [name, colour token]; a run keeps its colour whichever run is selected
const RUNDEF = { nohd64: ["nohd64 (no HD)", "--c7"], full6: ["full6", "--c6"], base64: ["base64 baseline", "--c8"], f128: ["f128 (full6 on 128 zones)", "--c6"], c146: ["c146 (146 zones)", "--run-b"], c128: ["c128 production", "--run-b"], full4: ["full4", "--c4"], full3: ["full3", "--c5"], fx3: ["fx3 (286 zones)", "--c7"], full2: ["full2", "--c3"], ref: ["reference", "--run-a"] };
const col = k => tok(RUNDEF[k][1]), nm = k => RUNDEF[k][0];

const PANELS = [
  { t: "Temperature (log K)", c: [["logT", "gas"], ["logTrad", "radiation", true]] },
  { t: "Density (log g cm⁻³)", c: [["logrho"]] },
  { t: "Velocity (km s⁻¹)", c: [["u"]], ref: 0, cur: true },
  { t: "Luminosity (L☉)", c: [["L"]], log: true },
  { t: "Ionised, molecular and H⁻ fraction", c: [["xHp", "x(H⁺)"], ["xH2", "x(H₂)", true], ["logxHm", "x(H⁻)", 2, true]], log: true, fix: [-7, 0.3] },
  { t: "Grey continuum optical depth (log τ_R)", c: [["logtau"]], ref: Math.log10(2 / 3) },
  { t: "Eddington ratio Γ", c: [["Gamma"]], log: true, ref: 1 },
  { t: "F / cE", c: [["FcE"]], log: true, ref: 1 },
  { t: "Inflow rate Ṁ (M☉ yr⁻¹)", c: [["Mdot"]], log: true },
  // nuclei, whatever their chemical state; colour is the nucleus (5th entry), the overlay run is dashed in the same colours
  { t: "H, D and He nuclei (mass fraction)", c: [["logXH", "¹H", false, true, "--run-a"], ["logXD", "D", false, true, "--run-b"], ["logXHe3", "³He", false, true, "--c3"], ["logXHe4", "⁴He", false, true, "--c4"]],
    log: true, fix: [-12, 0.3], leg: true },
  { t: "C, N and O (mass fraction)", c: [["logXC", "C", false, true, "--c5"], ["logXN", "N", false, true, "--c6"], ["logXO", "O", false, true, "--c7"]],
    log: true, cur: true, span: 12, floor: -30, leg: true, none: "No C, N or O above a mass fraction of 1e-30 in the range shown." },
];
// burning processes: [scalar, legend label, colour token (the fuel's colour in the abundance panels), short label]
const NUC = [["L_pp", "hydrogen (pp chain to ³He)", "--run-a", "pp"], ["L_D", "deuterium", "--run-b", "D"], ["L_He3", "³He", "--c3", "³He"],
  ["L_3a", "helium (3α, α captures)", "--c4", "3α"], ["L_CNO", "CNO cycle", "--c5", "CNO"]];
const MARKNAME = { D: "Deuterium burning", pp: "Hydrogen burning (pp chain)", CNO: "CNO cycle" };
const TRACK = [["M_star", "Protostar mass (M☉)", true], ["R_star", "Protostar radius (R☉)", true], ["r_ph", "Photosphere radius (R☉)", true],
  ["T_eff_lum", "T_eff from L and r_ph (K)", false], ["T_c", "Central temperature (K)", true], ["logn_c", "log central density (cm⁻³)", false],
  ["dt_s", "Time step (s)", true], ["u_min", "Fastest infall (km s⁻¹)", false], ["L_ph", "Luminosity at photosphere (L☉)", true]];

async function loadRun(key) {
  try {
    const meta = await (await fetch(key + "/frames.json")).json();
    const b64 = (await (await fetch(key + "/frames.b64.txt")).text()).trim(), bin = atob(b64), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    const buf = u8.buffer;
    const names = meta.fields.map(f => f[0]), nf = names.length;
    const frames = meta.frames.map(fr => {
      const a = new Float32Array(buf, fr.off, fr.n * nf), f = {};
      names.forEach((k, j) => f[k] = a.subarray(j * fr.n, (j + 1) * fr.n));
      const s = {}; meta.scalars.forEach((k, j) => s[k] = fr.s[j] === null ? NaN : fr.s[j]);
      return { f, s, n: fr.n };
    });
    const t0 = frames.length ? frames[0].s.t_yr : 0;
    frames.forEach(fr => fr.s.t = fr.s.t_yr - t0);
    RUNS[key] = { frames, ranges: {}, nuc: meta.nuc || null, t0 };
  } catch (e) { RUNS[key] = { frames: [], ranges: {}, err: String(e) }; }
  // spectra (optional): per spectrum ν L_ν at the edge, at the photosphere, r(τ_ν = 2/3), H₂ line luminosities at both places
  try {
    const meta = await (await fetch(key + "/spec.json")).json();
    const bin = atob((await (await fetch(key + "/spec.b64.txt")).text()).trim()), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    const n = meta.nnu, nl = meta.nl, per = 3 * n + 2 * nl, a = new Float32Array(u8.buffer);
    const S = meta.spectra.map((row, j) => { const s = {}; meta.scalars.forEach((k, i) => s[k] = row[i] === null ? NaN : row[i]);
      const o = j * per; return { s, Le: a.subarray(o, o + n), Lp: a.subarray(o + n, o + 2 * n), rt: a.subarray(o + 2 * n, o + 3 * n),
        le: a.subarray(o + 3 * n, o + 3 * n + nl), lp: a.subarray(o + 3 * n + nl, o + per) } });
    RUNS[key].spec = { lam: meta.lam.map(Math.log10), lamline: meta.lamline.map(Math.log10), S };
  } catch (e) { RUNS[key].spec = null }
}

function finite(v) { return Number.isFinite(v) }
function pct(arr, p) { if (!arr.length) return NaN; const s = arr.slice().sort((a, b) => a - b); return s[Math.min(s.length - 1, Math.max(0, Math.round(p / 100 * (s.length - 1))))] }
// x range of the run and y range of every panel; with win = [lo, hi] the y ranges are those of the zones inside that x window
function rangesFor(key, xk, win) {
  const R = RUNS[key], ck = win ? xk + "|" + win[0] + "|" + win[1] : xk; if (R.ranges[ck]) return R.ranges[ck];
  if (win) { const old = Object.keys(R.ranges).filter(k => k.includes("|")); if (old.length > 40) old.forEach(k => delete R.ranges[k]) }
  const xs = [];
  const step = Math.max(1, Math.floor(R.frames.length / 120));
  if (!win) {
    for (let i = 0; i < R.frames.length; i += step) for (const k of [xk, xk + "_o"]) if (R.frames[i].f[k]) for (const v of R.frames[i].f[k]) if (finite(v)) xs.push(v);
    const last = R.frames[R.frames.length - 1].f[xk + "_o"]; if (last) xs.push(last[last.length - 1]);
  }
  const out = { x: win || [pct(xs, 0), pct(xs, 100)], y: [] };
  for (const P of PANELS) {
    if (P.fix) { out.y.push(P.fix); continue }
    const vs = [];
    const take = fr => { const xc = fr.f[xk]; for (const [k, , , isl] of P.c) { const a = fr.f[k]; if (!a) continue;
      for (let j = 0; j < a.length; j++) { const v = a[j]; if (!finite(v) || (win && !(xc[j] >= win[0] && xc[j] <= win[1]))) continue;
        if (P.log && !isl) { if (v > 0) vs.push(Math.log10(v)) } else vs.push(v) } } };
    for (let i = 0; i < R.frames.length; i += step) take(R.frames[i]);
    if (!vs.length) { out.y.push(P.log ? [-3, 0] : [0, 1]); continue }
    if (P.full) { vs.length = 0; for (const fr of R.frames) take(fr) }
    let lo = pct(vs, P.full ? 0 : P.log ? 1 : 0.5), hi = pct(vs, P.full ? 100 : 99.5), pad = 0.05 * (hi - lo || 1);
    out.y.push(P.log ? [lo - 0.3, hi + 0.3] : [lo - pad, hi + pad]);
  }
  R.ranges[xk] = out; return out;
}

function niceTicks(lo, hi, n) {
  const span = hi - lo; if (!(span > 0)) return [lo];
  const raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag;
  const st = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag, t = [];
  for (let v = Math.ceil(lo / st) * st; v <= hi + 1e-9 * st; v += st) t.push(Math.abs(v) < st * 1e-9 ? 0 : v);
  return t;
}
function fmt(v) {
  if (!finite(v)) return "–"; const a = Math.abs(v);
  if (a !== 0 && (a >= 1e5 || a < 1e-3)) return v.toExponential(2);
  return a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : v.toPrecision(3);
}
function setup(cv) {
  const r = cv.getBoundingClientRect(), d = window.devicePixelRatio || 1;
  const W = Math.max(10, Math.round(r.width * d)), H = Math.max(10, Math.round(r.height * d));
  if (cv.width !== W || cv.height !== H) { cv.width = W; cv.height = H }
  const g = cv.getContext("2d"); g.setTransform(d, 0, 0, d, 0, 0); g.clearRect(0, 0, r.width, r.height);
  return { g, w: r.width, h: r.height };
}
// axes in data space; ylog: values are log10 already, ticks labelled as powers of ten
function axes(g, w, h, xl, yl, o) {
  const m = { l: 46, r: 8, t: 6, b: o.xlabel ? 30 : 20 };
  const X = v => m.l + (v - xl[0]) / (xl[1] - xl[0]) * (w - m.l - m.r);
  const Y = v => h - m.b - (v - yl[0]) / (yl[1] - yl[0]) * (h - m.t - m.b);
  g.font = "11px " + tok("--mono"); g.lineWidth = 1;
  g.strokeStyle = tok("--grid"); g.fillStyle = tok("--muted");
  const dec = []; if (o.ylog && yl[1] - yl[0] <= 9) for (let v = Math.ceil(yl[0] - 1e-9); v <= yl[1] + 1e-9; v++) dec.push(v);      // every decade
  const yt = dec.length >= 2 ? dec : o.ylog ? niceTicks(yl[0], yl[1], 4).filter(v => Math.abs(v - Math.round(v)) < 1e-6) : niceTicks(yl[0], yl[1], 4);
  const ytt = yt.length >= 2 ? yt : niceTicks(yl[0], yl[1], 4);
  g.textAlign = "right"; g.textBaseline = "middle";
  for (const v of ytt) { const y = Y(v); g.beginPath(); g.moveTo(m.l, y); g.lineTo(w - m.r, y); g.stroke();
    g.fillText(o.ylog ? (Math.abs(v - Math.round(v)) < 1e-6 ? "1e" + Math.round(v) : (10 ** v).toPrecision(2)) : fmt(v), m.l - 5, y) }
  g.textAlign = "center"; g.textBaseline = "top";
  for (const v of niceTicks(xl[0], xl[1], Math.max(3, Math.floor((w - 60) / 70))).filter(v => !o.xint || Math.abs(v - Math.round(v)) < 1e-6)) {
    const x = X(v); g.beginPath(); g.moveTo(x, m.t); g.lineTo(x, h - m.b); g.stroke();
    g.textAlign = x > w - m.r - 18 ? "right" : "center"; g.fillText(o.xfmt ? o.xfmt(v) : fmt(v), x > w - m.r - 18 ? w - 1 : x, h - m.b + 4) }
  g.textAlign = "center";
  if (o.xlabel) { g.fillStyle = tok("--faint"); g.fillText(o.xlabel, (m.l + w - m.r) / 2, h - 13) }
  X.w = w - m.r;
  return { X, Y, m };
}
function line(g, xs, ys, X, Y, col, dash, wid) {
  g.save(); g.strokeStyle = col; g.lineWidth = wid || 1.7; g.lineJoin = "round"; if (dash) g.setLineDash(dash === 2 ? [1.5, 3] : [5, 4]);
  g.beginPath(); let on = false;
  for (let i = 0; i < xs.length; i++) { const x = xs[i], y = ys[i];
    if (!finite(x) || !finite(y)) { on = false; continue } const px = X(x), py = Y(y); on ? g.lineTo(px, py) : g.moveTo(px, py); on = true }
  g.stroke(); g.restore();
}

// a zone's value drawn across the zone, from its inner to its outer interface (the first zone from its centre)
function steps(fr, xk, ys) {
  const xc = fr.f[xk], xo = fr.f[xk + "_o"];
  if (!xo) return [Array.from(xc), ys];
  const X = [], Y = [];
  for (let k = 0; k < ys.length; k++) { X.push(k === 0 ? xc[0] : xo[k - 1], xo[k]); Y.push(ys[k], ys[k]) }
  return [X, Y];
}
function zoneAt(fr, xk, x) {
  const xc = fr.f[xk], xo = fr.f[xk + "_o"];
  if (xo) { for (let k = 0; k < xo.length; k++) if (x <= xo[k]) return k; return xo.length - 1 }
  let bi = -1, bd = Infinity; for (let i = 0; i < xc.length; i++) { const d = Math.abs(xc[i] - x); if (d < bd) { bd = d; bi = i } } return bi;
}
// the run drawn dashed over the selected one: the chosen run, or (nearest other run) the first other run with a frame at this stage
function overlayRun(s) {
  const v = $("ovr").value; if (v === "none") return null;
  if (v !== "auto") return v === ST.run ? null : v;
  const c = Object.keys(RUNDEF).filter(k => k !== ST.run);
  return c.find(k => match(k, s) >= 0) || c[0];
}
function onTrack(key) { const v = $("ovr").value; return key === ST.run || v === "auto" || v === key }
function match(key, s) {
  const F = RUNS[key] ? RUNS[key].frames : []; if (!F.length) return -1;
  const useM = finite(s.M_star) && s.M_star >= 0.5 && s.logn_c > 23;
  let best = -1, bd = Infinity;
  F.forEach((fr, j) => { const q = fr.s;
    const d = useM ? (finite(q.M_star) && q.logn_c > 23 ? Math.abs(Math.log(q.M_star / s.M_star)) : Infinity) : Math.abs(q.logn_c - s.logn_c);
    if (d < bd) { bd = d; best = j } });
  return bd < (useM ? 0.1 : 0.3) ? best : -1;
}

// the x range of the profiles: the limits typed in; an empty box means the limit of the run's full range
function xwin(xk, auto) {
  const [a, b] = ST.xlim[xk]; if (a === null && b === null) return auto;
  const lo = a === null ? auto[0] : a, hi = b === null ? auto[1] : b;
  return hi > lo ? [lo, hi] : auto;
}
function syncX() {
  const xk = $("xax").value, R = RUNS[ST.run], auto = R && R.frames.length ? rangesFor(ST.run, xk).x : null;
  $("xlab").textContent = (xk === "logr" ? "log r (R☉)" : "log m (M☉)") + " from";
  [["xmin", 0], ["xmax", 1]].forEach(([id, j]) => { const v = ST.xlim[xk][j]; $(id).value = v === null ? "" : v; $(id).placeholder = auto && finite(auto[j]) ? auto[j].toFixed(2) : "" });
}
function buildProfiles() {
  const box = $("profiles"); box.innerHTML = "";
  PANELS.forEach((P, k) => {
    const c = document.createElement("div"); c.className = "card prof";
    c.innerHTML = `<div class="hd"><h3>${P.t}</h3><div class="val" id="pv${k}"></div></div>` +
      (P.leg ? `<div class="key leg">` + P.c.map(([, lab, , , ct]) => `<span><i class="sw" style="background:var(${ct})"></i>${lab}</span>`).join("") + `</div>` : "") + `<canvas id="pc${k}"></canvas>`;
    box.appendChild(c);
    const cv = c.querySelector("canvas");
    const mv = e => { const r = cv.getBoundingClientRect(); ST.cursor = { k, px: e.clientX - r.left }; drawProfiles() };
    cv.addEventListener("pointermove", mv); cv.addEventListener("pointerdown", mv);
    cv.addEventListener("pointerleave", () => { ST.cursor = null; drawProfiles() });
  });
}
function drawProfiles() {
  const R = RUNS[ST.run]; if (!R || !R.frames.length) return;
  const xk = $("xax").value, auto = rangesFor(ST.run, xk).x, xl = xwin(xk, auto), zoom = xl !== auto;
  const rg = zoom ? rangesFor(ST.run, xk, xl) : rangesFor(ST.run, xk), fr = R.frames[ST.i];
  const inx = (F, j) => !zoom || (F.f[xk][j] >= xl[0] && F.f[xk][j] <= xl[1]);
  const msg = (g, t) => { g.fillStyle = tok("--muted"); g.font = "12.5px " + tok("--sans"); g.textAlign = "left"; g.textBaseline = "top"; g.fillText(t, 12, 14) };
  const ok = overlayRun(fr.s), oj = ok ? match(ok, fr.s) : -1, ofr = oj >= 0 ? RUNS[ok].frames[oj] : null;
  const ca = col(ST.run), cb = ok ? col(ok) : ca;
  let cx = null;
  PANELS.forEach((P, k) => {
    const cv = $("pc" + k), { g, w, h } = setup(cv);
    let yl = rg.y[k];
    if (!P.c.some(([key]) => fr.f[key])) { msg(g, "Not recorded for this run."); return }
    if (P.cur) {      // range of the frame shown (and of the overlay frame) inside the x range shown, not of the whole run
      let lo = Infinity, hi = -Infinity;
      for (const F of [fr, ofr]) if (F) for (const [key] of P.c) { const a = F.f[key]; if (!a) continue;
        for (let j = 0; j < a.length; j++) { const v = a[j]; if (finite(v) && inx(F, j) && !(v < P.floor)) { if (v < lo) lo = v; if (v > hi) hi = v } } }
      if (P.none && !(lo <= hi)) { msg(g, P.none); return }
      if (P.ref !== undefined) { lo = Math.min(lo, P.ref); hi = Math.max(hi, P.ref) }
      if (P.span && hi - lo > P.span) lo = hi - P.span;
      if (lo <= hi) { const pad = 0.06 * (hi - lo || 1); yl = [lo - pad, hi + pad] }
    }
    const xsp = xl[1] - xl[0], xdg = Math.min(6, Math.max(2, Math.ceil(-Math.log10(xsp)) + 1));      // enough digits for a narrow range
    const A = axes(g, w, h, xl, yl, { ylog: P.log, xlabel: xk === "logr" ? "log r (R☉)" : "log m (M☉)", xfmt: xsp < 1 ? v => v.toFixed(xdg) : null });
    g.save(); g.beginPath(); g.rect(A.m.l, A.m.t, w - A.m.l - A.m.r, h - A.m.t - A.m.b); g.clip();
    if (P.ref !== undefined) { g.strokeStyle = tok("--faint"); g.setLineDash([2, 3]); g.beginPath(); const y = A.Y(P.log ? Math.log10(P.ref) : P.ref); g.moveTo(A.m.l, y); g.lineTo(w, y); g.stroke(); g.setLineDash([]) }
    const tr = v => P.log ? (v > 0 ? Math.log10(v) : NaN) : v, id = v => v;      // a field stored as log10 is drawn as it is
    if (ofr) P.c.forEach(([key, , alt, isl, ct]) => { if (!ofr.f[key]) return; const [xs, ys] = steps(ofr, xk, Array.from(ofr.f[key], isl ? id : tr)); line(g, xs, ys, A.X, A.Y, ct ? tok(ct) : cb, alt === 2 ? 2 : true, alt ? 1.2 : 1.4) });
    P.c.forEach(([key, , alt, isl, ct]) => { if (!fr.f[key]) return; const [xs, ys] = steps(fr, xk, Array.from(fr.f[key], isl ? id : tr)); line(g, xs, ys, A.X, A.Y, ct ? tok(ct) : ca, alt, 1.8) });
    g.restore();
    if (ST.cursor) {
      const xv = xl[0] + (ST.cursor.px - A.m.l) / (w - A.m.l - A.m.r) * (xl[1] - xl[0]); cx = xv;
      const px = A.X(xv); if (px >= A.m.l && px <= w - A.m.r) { g.strokeStyle = tok("--muted"); g.lineWidth = 1; g.beginPath(); g.moveTo(px, A.m.t); g.lineTo(px, h - A.m.b); g.stroke() }
    }
  });
  PANELS.forEach((P, k) => {
    const el = $("pv" + k);
    if (cx === null || !P.c.some(([key]) => fr.f[key])) { el.textContent = ""; return }
    const bi = zoneAt(fr, xk, cx);
    const xlab = xk === "logr" ? "r " + fmt(10 ** fr.f.logr[bi]) + " R☉" : "m " + fmt(10 ** fr.f.logm[bi]) + " M☉";      // the zone's centre
    el.textContent = bi < 0 ? "" : xlab + "  ·  " + P.c.map(([key, lab, , isl]) => (lab ? lab + " " : "") + (fr.f[key] ? fmt(isl ? 10 ** fr.f[key][bi] : fr.f[key][bi]) : "–")).join("  ·  ");
  });
  const kk = $("key");
  kk.innerHTML = `<span><i class="sw" style="background:${ca}"></i>${nm(ST.run)} frame ${ST.i + 1} of ${R.frames.length}</span>` +
    (ok ? (ofr ? `<span><i class="sw" style="background:${cb};height:0;border-top:2px dashed ${cb}"></i>${nm(ok)} at step ${ofr.s.step} (matched)</span>` : `<span>${nm(ok)}: no frame at this stage yet</span>`) : "") +
    `<span>Second curve in a panel (radiation T, x(H₂)) is drawn dashed in the same colour, x(H⁻) dotted. In the two abundance panels the colour is the nucleus (legend in the panel) and the overlay run is dashed.</span>`;
}

function chips() {
  const s = RUNS[ST.run].frames[ST.i].s, c = [];
  const add = (lab, v) => c.push(`<div class="chip"><span>${lab}</span><b>${v}</b></div>`);
  add("step", s.step); add("t", s.t_yr.toExponential(6) + " yr"); add("Δt", fmt(s.dt_s) + " s"); add("log n_c", s.logn_c.toFixed(2));
  add("T_c", fmt(s.T_c) + " K");
  if (s.logn_c > 23 && finite(s.M_star)) { add("M*", s.M_star.toFixed(3) + " M☉"); add("R*", fmt(s.R_star) + " R☉") }
  if (finite(s.r_ph)) { add("r_ph", fmt(s.r_ph) + " R☉"); add("T_eff", fmt(s.T_eff_lum) + " K") }
  add("zones", s.nz);
  $("chips").innerHTML = c.join("");
}

// the spectrum nearest in step to the frame shown
function specFor(key, step) {
  const sp = RUNS[key] && RUNS[key].spec; if (!sp || !sp.S.length) return null;
  let b = null, bd = Infinity; for (const q of sp.S) { const d = Math.abs(q.s.step - step); if (d < bd) { bd = d; b = q } } return b;
}
function drawSpec() {
  const cv = $("spec"), cv2 = $("spect"), R = RUNS[ST.run], fr = R.frames[ST.i], sp = R.spec, q = specFor(ST.run, fr.s.step), ca = col(ST.run);
  if (!q) {
    for (const c of [cv, cv2]) { const { g } = setup(c); g.fillStyle = tok("--muted"); g.font = "13px " + tok("--sans"); g.fillText("No spectra for this run yet.", 16, 24) }
    $("specchips").innerHTML = ""; $("speckey").innerHTML = ""; $("spv").textContent = ""; $("spv2").textContent = ""; return }
  const lg = a => Array.from(a, v => v > 0 ? Math.log10(v) : NaN);
  const ye = lg(q.Le), yp = lg(q.Lp), yli = lg(q.le), yr = lg(q.rt), s = q.s;
  let bb = null;        // blackbody of the run's T_eff at the run's photosphere radius
  if (finite(s.T_eff) && finite(s.r_ph)) { const A4 = 4 * Math.PI * (s.r_ph * 6.957e10) ** 2;
    bb = sp.lam.map(l => { const nu = 2.99792e10 / (10 ** l * 1e-4), x = 6.62607e-27 * nu / (1.380649e-16 * s.T_eff);
      const v = A4 * Math.PI * 2 * 6.62607e-27 * nu ** 4 / 2.99792e10 ** 2 / Math.expm1(x) / 3.83e33; return v > 0 && finite(v) ? Math.log10(v) : NaN }) }
  const all = ye.concat(yp, yli).filter(finite), top = all.length ? Math.ceil(Math.max(...all) + 0.3) : 0;
  const xl = [-1.35, 3.05], yl = [top - 9, top], lamfmt = v => { const p = 10 ** v; return p >= 1 ? String(Math.round(p)) : p.toPrecision(1) };
  let cl = null;
  {
    const { g, w, h } = setup(cv), A = axes(g, w, h, xl, yl, { ylog: true, xint: true, xfmt: lamfmt });
    g.save(); g.beginPath(); g.rect(A.m.l, A.m.t, w - A.m.l - A.m.r, h - A.m.t - A.m.b); g.clip();
    g.strokeStyle = tok("--c3"); g.lineWidth = 1; g.beginPath();
    sp.lamline.forEach((l, j) => { if (!finite(yli[j]) || yli[j] < yl[0]) return; const x = Math.round(A.X(l)) + 0.5; g.moveTo(x, A.Y(yl[0])); g.lineTo(x, A.Y(yli[j])) }); g.stroke();
    if (bb) { g.save(); g.setLineDash([2, 3]); line(g, sp.lam, bb, A.X, A.Y, tok("--faint"), false, 1.3); g.restore() }
    line(g, sp.lam, yp, A.X, A.Y, ca, true, 1.5); line(g, sp.lam, ye, A.X, A.Y, ca, false, 1.9);
    g.restore();
    if (ST.scur !== null && ST.scur !== undefined) { const xv = xl[0] + (ST.scur - A.m.l) / (w - A.m.l - A.m.r) * (xl[1] - xl[0]);
      if (xv >= xl[0] && xv <= xl[1]) { cl = xv; const px = A.X(xv); g.strokeStyle = tok("--muted"); g.lineWidth = 1; g.beginPath(); g.moveTo(px, A.m.t); g.lineTo(px, h - A.m.b); g.stroke() } }
  }
  {
    const { g, w, h } = setup(cv2), fin = yr.filter(finite), rp = s.r_ph > 0 ? Math.log10(s.r_ph) : NaN;
    if (!fin.length) { g.fillStyle = tok("--muted"); g.font = "13px " + tok("--sans"); g.fillText("Transparent at every wavelength: τ_λ stays below 2/3 down to the centre.", 16, 24) }
    else {
      if (finite(rp)) fin.push(rp);
      const lo = Math.min(...fin), hi = Math.max(...fin), pad = Math.max(0.15, 0.08 * (hi - lo)), A = axes(g, w, h, xl, [lo - pad, hi + pad], { ylog: true, xint: true, xfmt: lamfmt });
      g.save(); g.beginPath(); g.rect(A.m.l, A.m.t, w - A.m.l - A.m.r, h - A.m.t - A.m.b); g.clip();
      if (finite(rp)) { g.strokeStyle = tok("--faint"); g.setLineDash([2, 3]); g.beginPath(); g.moveTo(A.m.l, A.Y(rp)); g.lineTo(w, A.Y(rp)); g.stroke(); g.setLineDash([]) }
      line(g, sp.lam, yr, A.X, A.Y, ca, false, 1.9); g.restore();
      if (cl !== null) { const px = A.X(cl); g.strokeStyle = tok("--muted"); g.lineWidth = 1; g.beginPath(); g.moveTo(px, A.m.t); g.lineTo(px, h - A.m.b); g.stroke() }
    }
  }
  if (cl !== null) {
    let bi = 0, bd = Infinity; sp.lam.forEach((l, j) => { const d = Math.abs(l - cl); if (d < bd) { bd = d; bi = j } });
    let li = -1, ld = 0.012; sp.lamline.forEach((l, j) => { const d = Math.abs(l - cl); if (d < ld && q.le[j] > 0) { ld = d; li = j } });
    $("spv").textContent = "λ " + fmt(10 ** sp.lam[bi]) + " μm  ·  edge " + fmt(q.Le[bi]) + "  ·  photosphere " + fmt(q.Lp[bi]) + (li >= 0 ? "  ·  H₂ line " + (10 ** sp.lamline[li]).toFixed(2) + " μm " + fmt(q.le[li]) : "");
    $("spv2").textContent = "λ " + fmt(10 ** sp.lam[bi]) + " μm  ·  r " + fmt(q.rt[bi]) + " R☉";
  } else { $("spv").textContent = ""; $("spv2").textContent = "" }
  const c = [], add = (lab, v) => c.push(`<div class="chip"><span>${lab}</span><b>${v}</b></div>`), rat = (a, b) => finite(a) && finite(b) && b > 0 ? " (ratio " + (a / b).toFixed(2) + ")" : "";
  add("spectrum of step", s.step + (s.step !== fr.s.step ? " (frame: " + fr.s.step + ")" : ""));
  add("edge ∫L_ν dν", fmt(s.Lc_edge) + " L☉"); add("run's L at edge", s.Lrun_edge === 0 ? "none: no continuum transport yet" : fmt(s.Lrun_edge) + " L☉" + rat(s.Lc_edge, s.Lrun_edge));
  if (s.kph > 0) { add("photosphere ∫L_ν dν", fmt(s.Lc_ph) + " L☉"); add("run's L there", fmt(s.Lrun_ph) + " L☉" + rat(s.Lc_ph, s.Lrun_ph)); add("T_eff", fmt(s.T_eff) + " K") }
  else add("photosphere", "none: τ_R < 2/3 to the centre");
  add("H₂ lines at edge", fmt(s.Ll_edge) + " L☉"); add("run's L_lines", fmt(s.Ll_run) + " L☉");
  if (finite(s.lscale) && Math.abs(s.lscale - 1) > 0.02) add("line sticks scaled by", s.lscale.toFixed(2));
  $("specchips").innerHTML = c.join("");
  $("speckey").innerHTML = `<span><i class="sw" style="background:${ca}"></i>continuum at the mesh edge</span>` +
    `<span><i class="sw" style="background:${ca};height:0;border-top:2px dashed ${ca}"></i>continuum through the photosphere</span>` +
    `<span><i class="sw" style="background:${tok("--c3")};width:2px;height:12px"></i>H₂ lines at the edge (luminosity per line)</span>` +
    `<span><i class="sw" style="background:${tok("--faint")};height:0;border-top:2px dotted ${tok("--faint")}"></i>blackbody of T_eff at r_ph · grey photosphere</span>`;
}
for (const id of ["spec", "spect"]) {
  const mv = e => { const r = $(id).getBoundingClientRect(); ST.scur = e.clientX - r.left; drawSpec() };
  $(id).addEventListener("pointermove", mv); $(id).addEventListener("pointerdown", mv);
  $(id).addEventListener("pointerleave", () => { ST.scur = null; drawSpec() });
}

function drawTrack() {
  const cv = $("track"), { g, w, h } = setup(cv);
  const [yk, , ylog] = TRACK.find(t => t[0] === $("ty").value), xk = $("tx").value;
  const series = [];
  for (const key of ["ref", "full2", "full3", "full4", "full6", "fx3", "c146", "f128", "c128", "nohd64", "base64"]) {
    const R = RUNS[key]; if (!R || !R.frames.length) continue;
    if (!onTrack(key)) continue;
    const pts = R.frames.map(fr => { let x = fr.s[xk], y = fr.s[yk];
      if ((yk === "M_star" || yk === "R_star" || xk === "M_star") && !(fr.s.logn_c > 23)) { if (yk === "M_star" || yk === "R_star") y = NaN; if (xk === "M_star") x = NaN }
      if (ylog) y = y > 0 ? Math.log10(y) : NaN;
      if (xk === "M_star") x = x > 0 ? Math.log10(x) : NaN;
      return [x, y] });
    series.push({ key, pts });
  }
  const xs = [], ys = []; series.forEach(S => S.pts.forEach(([x, y]) => { if (finite(x) && finite(y)) { xs.push(x); ys.push(y) } }));
  if (!xs.length) { g.fillStyle = tok("--muted"); g.font = "13px " + tok("--sans"); g.fillText("No values of this quantity yet.", 16, 24); cv._map = null; return }
  let xl = [Math.min(...xs), Math.max(...xs)], yl = yk === "u_min" ? [Math.min(...ys), Math.max(...ys)] : [pct(ys, 0.5), pct(ys, 99.5)];
  if (xl[1] === xl[0]) xl = [xl[0] - 1, xl[1] + 1]; const p = 0.05 * (yl[1] - yl[0] || 1); yl = [yl[0] - p, yl[1] + p];
  const A = axes(g, w, h, xl, yl, { ylog, xlabel: xk === "t" ? "years since the run's first frame" : xk === "step" ? "step" : "protostar mass (M☉)", xfmt: xk === "M_star" ? v => fmt(10 ** v) : null });
  g.save(); g.beginPath(); g.rect(A.m.l, A.m.t, w - A.m.l - A.m.r, h - A.m.t - A.m.b); g.clip();
  for (const S of series) line(g, S.pts.map(p => p[0]), S.pts.map(p => p[1]), A.X, A.Y, col(S.key), S.key !== ST.run, 1.6);
  const cur = series.find(S => S.key === ST.run);
  if (cur) { const [x, y] = cur.pts[ST.i];
    if (finite(x)) { g.strokeStyle = tok("--fg"); g.lineWidth = 1; g.beginPath(); g.moveTo(A.X(x), A.m.t); g.lineTo(A.X(x), h - A.m.b); g.stroke() }
    if (finite(x) && finite(y)) { g.fillStyle = col(ST.run); g.strokeStyle = tok("--panel"); g.lineWidth = 2; g.beginPath(); g.arc(A.X(x), A.Y(y), 4.5, 0, 7); g.fill(); g.stroke() } }
  g.restore();
  drawMarks(g, A, h, xk);
  cv._map = { A, pts: cur ? cur.pts : [] };
  // readout: the selected run's point nearest the pointer, or the frame shown
  let ri = cur ? ST.i : -1;
  if (cur && ST.tcur !== null && ST.tcur !== undefined) { let bd = Infinity; cur.pts.forEach(([x, y], j) => { if (!finite(x) || !finite(y)) return; const d = Math.abs(A.X(x) - ST.tcur); if (d < bd) { bd = d; ri = j } });
    const [x, y] = cur.pts[ri]; if (finite(x) && finite(y)) { g.strokeStyle = tok("--muted"); g.lineWidth = 1; g.beginPath(); g.moveTo(A.X(x), A.m.t); g.lineTo(A.X(x), h - A.m.b); g.stroke();
      g.strokeStyle = col(ST.run); g.lineWidth = 1.5; g.beginPath(); g.arc(A.X(x), A.Y(y), 5, 0, 7); g.stroke() } }
  $("trh").textContent = TRACK.find(t => t[0] === yk)[1];
  if (ri >= 0) { const q = R0(ri), xv = xk === "t" ? "t " + fmt(q.t) + " yr" : xk === "step" ? "step " + q.step : "M* " + fmt(q.M_star) + " M☉";
    $("trv").textContent = xv + "  ·  " + fmt(q[yk]) + (xk !== "step" ? "  ·  step " + q.step : "") } else $("trv").textContent = "";
}
function R0(i) { return RUNS[ST.run].frames[i].s }

// onset of D, pp and CNO burning in the selected run, as labelled dashed lines on a chart whose x axis is xk
function drawMarks(g, A, h, xk, bottom) {
  const R = RUNS[ST.run]; if (!R || !R.nuc) return;
  g.save(); g.font = "11px " + tok("--mono"); g.textBaseline = "top"; let row = 0;
  for (const m of R.nuc.marks) {
    if (m.state !== "at") continue;
    const x = xk === "t" ? m.t_yr - R.t0 : xk === "step" ? m.step : m.M_star > 0 ? Math.log10(m.M_star) : NaN; if (!finite(x)) continue;
    const px = A.X(x); g.strokeStyle = tok("--muted"); g.lineWidth = 1; g.setLineDash([4, 3]); g.beginPath(); g.moveTo(px, A.m.t); g.lineTo(px, h - A.m.b); g.stroke(); g.setLineDash([]);
    const lab = m.key + " onset", tw = g.measureText(lab).width, left = px + 4 + tw > A.X.w;
    g.fillStyle = tok("--muted"); g.textAlign = left ? "right" : "left"; g.fillText(lab, px + (left ? -4 : 4), bottom ? h - A.m.b - 15 - 13 * row++ : A.m.t + 2 + 13 * row++);
  }
  g.restore();
}
function drawNuc() {
  const cv = $("nuc"), { g, w, h } = setup(cv), R = RUNS[ST.run], xk = $("tx").value;
  const note = t => { g.fillStyle = tok("--muted"); g.font = "13px " + tok("--sans"); g.fillText(t, 16, 24); cv._map = null; $("nuv").textContent = "" };
  $("nuckey").innerHTML = NUC.map(([, lab, c]) => `<span><i class="sw" style="background:${tok(c)}"></i>${lab}</span>`).join("") +
    `<span><i class="sw" style="background:${tok("--faint")}"></i>luminosity at the photosphere (running median)</span>` +
    `<span><i class="sw" style="background:none;height:0;border-top:2px dotted ${tok("--faint")}"></i>0.1 % of it: the onset level (dotted)</span>`;
  if (!R.nuc) { $("nucchips").innerHTML = ""; note("The nuclear abundances of this run were not recorded."); return }
  const lg = v => v > 0 ? Math.log10(v) : NaN;
  const xs = R.frames.map(fr => { const s = fr.s; return xk === "M_star" ? (s.logn_c > 23 ? lg(s.M_star) : NaN) : s[xk] });
  const ser = NUC.map(([k]) => R.frames.map(fr => lg(fr.s[k]))), ref = R.frames.map(fr => lg(fr.s.L_ph_med)), thr = ref.map(v => v + Math.log10(R.nuc.thr));
  // status of the three marked processes
  $("nucchips").innerHTML = R.nuc.marks.map(m => {
    const M = finite(m.M_star) && m.M_star > 0 ? m.M_star.toFixed(2) + " M☉" : null, p = 100 * m.frac;
    const f = finite(p) && p > 0 ? (p >= 0.01 ? p.toPrecision(2) : p.toExponential(1)) + " % of L in the latest snapshot" : "no energy from it yet";
    const v = m.state === "at" ? "onset at " + (M ? M + ", " : "") + "step " + m.step : m.state === "before" ? "already on in the first recorded snapshot" + (M ? " (" + M + ")" : "") + "; now " + f : "not yet, " + f;
    return `<div class="chip"><span>${MARKNAME[m.key] || m.key}</span><b>${v}</b></div>` }).join("");
  const fx = xs.filter(finite), fy = ser.flat().concat(ref, thr).filter(finite);
  if (!fx.length || !fy.length) { note(xk === "M_star" ? "No protostar in this run yet: choose time or step for the track's x axis." : "No luminosities to show yet."); return }
  let xl = [Math.min(...fx), Math.max(...fx)]; if (xl[1] === xl[0]) xl = [xl[0] - 1, xl[1] + 1];
  // by default down to six decades below the onset level; "All decades" shows everything
  const ft = thr.filter(finite), hi = Math.ceil(Math.max(...fy) + 0.2), all = $("nucall").checked || !ft.length;
  const lo = all ? Math.max(Math.floor(Math.min(...fy)), hi - 40) : Math.floor(Math.min(...ft)) - 6;
  const A = axes(g, w, h, xl, [lo, hi], { ylog: true, xlabel: xk === "t" ? "years since the run's first frame" : xk === "step" ? "step" : "protostar mass (M☉)",
    xfmt: xk === "M_star" ? v => { const m = 10 ** v, sp = xl[1] - xl[0]; return sp < 0.3 && m >= 1 ? m.toFixed(sp < 0.03 ? 2 : 1) : fmt(m) } : null });
  g.save(); g.beginPath(); g.rect(A.m.l, A.m.t, w - A.m.l - A.m.r, h - A.m.t - A.m.b); g.clip();
  line(g, xs, ref, A.X, A.Y, tok("--faint"), false, 1.3); line(g, xs, thr, A.X, A.Y, tok("--faint"), 2, 1.5);
  NUC.forEach(([, , c], j) => line(g, xs, ser[j], A.X, A.Y, tok(c), false, 1.9));
  const x0 = xs[ST.i]; if (finite(x0)) { g.strokeStyle = tok("--fg"); g.lineWidth = 1; g.beginPath(); g.moveTo(A.X(x0), A.m.t); g.lineTo(A.X(x0), h - A.m.b); g.stroke() }
  g.restore();
  drawMarks(g, A, h, xk, true);
  cv._map = { A, xs };
  let ri = ST.i;
  if (ST.ncur !== null && ST.ncur !== undefined) { let bd = Infinity; xs.forEach((x, j) => { if (!finite(x)) return; const d = Math.abs(A.X(x) - ST.ncur); if (d < bd) { bd = d; ri = j } });
    const x = xs[ri]; if (finite(x)) { g.strokeStyle = tok("--muted"); g.lineWidth = 1; g.beginPath(); g.moveTo(A.X(x), A.m.t); g.lineTo(A.X(x), h - A.m.b); g.stroke();
      NUC.forEach(([, , c], j) => { const y = ser[j][ri]; if (!finite(y) || y < lo) return; g.fillStyle = tok(c); g.strokeStyle = tok("--panel"); g.lineWidth = 2; g.beginPath(); g.arc(A.X(x), A.Y(y), 4, 0, 7); g.fill(); g.stroke() }) } }
  const q = R.frames[ri].s, xv = xk === "t" ? "t " + fmt(q.t) + " yr" : xk === "step" ? "step " + q.step : "M* " + fmt(q.M_star) + " M☉";
  $("nuv").textContent = xv + "  ·  " + NUC.map(([k, , , sh]) => sh + " " + (q[k] > 0 ? fmt(q[k]) : "0")).join("  ·  ") + "  ·  L " + fmt(q.L_ph_med) + (xk !== "step" ? "  ·  step " + q.step : "");
}
{ const mv = e => { const r = $("nuc").getBoundingClientRect(); ST.ncur = e.clientX - r.left; drawNuc() };
  $("nuc").addEventListener("pointermove", mv); $("nuc").addEventListener("pointerleave", () => { ST.ncur = null; drawNuc() });
  $("nuc").addEventListener("click", e => { const m = $("nuc")._map; if (!m) return; const px = e.clientX - $("nuc").getBoundingClientRect().left;
    let bi = -1, bd = Infinity; m.xs.forEach((x, j) => { if (!finite(x)) return; const d = Math.abs(m.A.X(x) - px); if (d < bd) { bd = d; bi = j } }); if (bi >= 0) seek(bi) }) }
{ const mv = e => { const r = $("track").getBoundingClientRect(); ST.tcur = e.clientX - r.left; drawTrack() };
  $("track").addEventListener("pointermove", mv); $("track").addEventListener("pointerleave", () => { ST.tcur = null; drawTrack() }) }
$("track").addEventListener("click", e => {
  const m = $("track")._map; if (!m) return; const r = $("track").getBoundingClientRect(), px = e.clientX - r.left;
  let bi = -1, bd = Infinity; m.pts.forEach(([x], j) => { if (!finite(x)) return; const d = Math.abs(m.A.X(x) - px); if (d < bd) { bd = d; bi = j } });
  if (bi >= 0) seek(bi);
});

let PERF = []; const PERFS = {};
function parsePerf(txt) {
  const L = txt.trim().split("\n"), hd = L[0].split(",");
  return L.slice(1).map(l => { const v = l.split(","), o = {}; hd.forEach((k, j) => o[k] = k === "dtlim" ? v[j] : +v[j]); return o }).filter(o => o.nsteps > 0);
}
function phaseOf(o) {
  const n = o.log_nc;
  if (n < 23.35) return n < 8 ? "Collapse: log n_c < 8" : n < 16 ? "Collapse: log n_c 8–16" : n < 21 ? "Collapse: log n_c 16–21" : "Collapse: log n_c 21 to remap";
  const M = o.mstar; return !(M > 0) ? "Accretion: no protostar found" : M < 10 ? "Accretion: M* < 10 M☉" : M < 30 ? "Accretion: M* 10–30 M☉" : M < 100 ? "Accretion: M* 30–100 M☉" : "Accretion: M* > 100 M☉";
}
const PCATS = [["t_step", "Hydro and radiation", "--run-a"], ["t_chem", "Chemistry", "--run-b"], ["t_lines", "H₂ lines", "--c3"], ["t_other", "Output and the rest", "--c4"]];
function drawPerf() {
  const cv = $("perf"), cv2 = $("perf2");
  { const pk = PERFS[ST.run] && PERFS[ST.run].length ? ST.run : "base64"; PERF = PERFS[pk] || []; $("perfh").textContent = "Where the time goes (" + nm(pk) + ")" }
  $("perfkey").innerHTML = PCATS.map(([, lab, c]) => `<span><i class="sw" style="background:${tok(c)};height:8px;width:10px"></i>${lab}</span>`).join("") +
    `<span><i class="sw" style="background:${tok("--c5")}"></i>GC share</span>`;
  if (!PERF.length) { const { g } = setup(cv); g.fillStyle = tok("--muted"); g.font = "13px " + tok("--sans"); g.fillText("No cost records yet.", 16, 24); setup(cv2); return }
  {
    const { g, w, h } = setup(cv);
    const xs = PERF.map(o => o.step), tot = PERF.map(o => (o.t_step + o.t_chem + o.t_lines + o.t_other) / o.nsteps);
    const xl = [Math.min(...xs) - 50, Math.max(...xs) + 50], yl = [0, Math.max(pct(tot, 98) * 1.1, 1e-3)];
    const A = axes(g, w, h, xl, yl, { xlabel: "step  ·  wall seconds per step" });
    const bw = Math.max(1, (A.X(xl[0] + 100) - A.X(xl[0])) - 1);
    g.save(); g.beginPath(); g.rect(A.m.l, A.m.t, w - A.m.l - A.m.r, h - A.m.t - A.m.b); g.clip();
    PERF.forEach(o => { let y0 = 0; const x = A.X(o.step - o.nsteps / 2) - bw / 2;
      for (const [k, , c] of PCATS) { const v = o[k] / o.nsteps; g.fillStyle = tok(c); const yt = A.Y(y0 + v), yb = A.Y(y0); g.fillRect(x, yt, bw, Math.max(0, yb - yt - (bw > 4 ? 1 : 0))); y0 += v } });
    g.restore();
  }
  {
    const { g, w, h } = setup(cv2);
    const xs = PERF.map(o => o.step), div = [], gcs = [];
    PERF.forEach((o, j) => { const p = j > 0 ? PERF[j - 1] : null;
      div.push(p && p.pid === o.pid ? (o.ndiv - p.ndiv) / o.nsteps : NaN);
      gcs.push(o.gc_s / Math.max(1e-9, o.t_step + o.t_chem + o.t_lines + o.t_other)) });
    const yl = [0, Math.max(1, ...div.filter(finite)) * 1.1];
    const A = axes(g, w, h, [Math.min(...xs) - 50, Math.max(...xs) + 50], yl, { xlabel: "step  ·  Newton divergences per step (bars), GC share of wall time (line, same axis)" });
    const bw = Math.max(1, (A.X(Math.min(...xs) + 100) - A.X(Math.min(...xs))) - 1);
    g.fillStyle = tok("--faint"); PERF.forEach((o, j) => { if (!finite(div[j])) return; const yt = A.Y(div[j]); g.fillRect(A.X(o.step - o.nsteps / 2) - bw / 2, yt, bw, A.Y(0) - yt) });
    line(g, xs, gcs, A.X, A.Y, tok("--c5"), false, 1.8);
  }
  const ph = new Map();
  PERF.forEach((o, j) => { const k = phaseOf(o), p = j > 0 ? PERF[j - 1] : null;
    const r = ph.get(k) || { n: 0, step: 0, chem: 0, lines: 0, oth: 0, gc: 0, alloc: 0, div: 0, first: o.step, last: o.step };
    r.n += o.nsteps; r.step += o.t_step; r.chem += o.t_chem; r.lines += o.t_lines; r.oth += o.t_other; r.gc += o.gc_s; r.alloc += o.alloc_GB;
    if (p && p.pid === o.pid) r.div += o.ndiv - p.ndiv; r.last = o.step; ph.set(k, r) });
  let html = "<thead><tr><th>Phase</th><th>steps</th><th>wall h</th><th>s / step</th><th>hydro+rad</th><th>chemistry</th><th>H₂ lines</th><th>other</th><th>GC</th><th>alloc GB / step</th><th>divergences / step</th></tr></thead><tbody>";
  for (const [k, r] of ph) { const T = r.step + r.chem + r.lines + r.oth, P = v => (100 * v / T).toFixed(0) + " %";
    html += `<tr><td>${k}</td><td class="mono">${r.n}</td><td class="mono">${(T / 3600).toFixed(2)}</td><td class="mono">${(T / r.n).toFixed(3)}</td><td class="mono">${P(r.step)}</td><td class="mono">${P(r.chem)}</td><td class="mono">${P(r.lines)}</td><td class="mono">${P(r.oth)}</td><td class="mono">${P(r.gc)}</td><td class="mono">${(r.alloc / r.n).toFixed(2)}</td><td class="mono">${(r.div / r.n).toFixed(2)}</td></tr>` }
  $("phases").innerHTML = html + "</tbody>";
}

function seek(i) { const R = RUNS[ST.run]; if (!R.frames.length) return; ST.i = Math.max(0, Math.min(R.frames.length - 1, i)); $("frame").value = ST.i; drawAll() }
function drawAll() { const R = RUNS[ST.run]; if (!R || !R.frames.length) { $("chips").innerHTML = `<div class="empty">No frames in this run yet${R && R.err ? " (" + R.err + ")" : ""}.</div>`; return } chips(); drawProfiles(); drawSpec(); drawTrack(); drawNuc(); drawPerf() }
function selectRun() { ST.run = $("run").value; const R = RUNS[ST.run]; $("frame").max = Math.max(0, R.frames.length - 1); syncX(); seek(R.frames.length - 1) }

$("run").addEventListener("change", selectRun);
$("ovr").addEventListener("change", drawAll);
$("xax").addEventListener("change", () => { syncX(); drawAll() });
[["xmin", 0], ["xmax", 1]].forEach(([id, j]) => $(id).addEventListener("input", e => {
  const t = e.target.value.trim(), v = t === "" ? null : Number(t); ST.xlim[$("xax").value][j] = v !== null && finite(v) ? v : null; drawProfiles() }));
$("xauto").addEventListener("click", () => { ST.xlim[$("xax").value] = [null, null]; syncX(); drawProfiles() });
$("ty").addEventListener("change", drawTrack); $("tx").addEventListener("change", () => { drawTrack(); drawNuc() }); $("nucall").addEventListener("change", drawNuc);
$("frame").addEventListener("input", e => seek(+e.target.value));
$("prev").addEventListener("click", () => seek(ST.i - 1)); $("next").addEventListener("click", () => seek(ST.i + 1));
$("play").addEventListener("click", () => {
  ST.playing = !ST.playing; $("play").textContent = ST.playing ? "Pause" : "Play";
  clearInterval(ST.timer);
  if (ST.playing) { if (ST.i >= RUNS[ST.run].frames.length - 1) seek(0);
    ST.timer = setInterval(() => { if (ST.i >= RUNS[ST.run].frames.length - 1) { if ($("loop").checked) seek(0); else $("play").click(); return } seek(ST.i + 1) }, 110) }
});
document.addEventListener("keydown", e => { if (e.target.tagName === "SELECT" || e.target.type === "number") return; if (e.key === "ArrowLeft") seek(ST.i - 1); if (e.key === "ArrowRight") seek(ST.i + 1) });
let rz; window.addEventListener("resize", () => { clearTimeout(rz); rz = setTimeout(drawAll, 120) });
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", drawAll);
new MutationObserver(drawAll).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });

TRACK.forEach(([k, lab]) => { const o = document.createElement("option"); o.value = k; o.textContent = lab; $("ty").appendChild(o) });
buildProfiles();
(async () => {
  await Promise.all(Object.keys(RUNDEF).map(loadRun));
  for (const k of ["base64", "nohd64", "c128", "f128", "full4", "full3", "full6", "c146", "fx3", "full2"]) { try { PERFS[k] = parsePerf(await (await fetch(k + "/perf.csv")).text()) } catch (e) { PERFS[k] = [] } }
  if (!RUNS.base64.frames.length) $("run").value = RUNS.c128.frames.length ? "c128" : "full6";
  selectRun();
})();
