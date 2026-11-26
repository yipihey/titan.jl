// plot.js - crisp, theme-aware SVG line/scatter & 2D heatmap plots with interactive zoom/pan,
// per-axis lin/log & min-max controls, SVG/PNG/CSV export, and linked cross-plot highlighting.
//
// Axes, ticks and text use `currentColor`, so a plot inherits the text colour
// of its container and works on light and dark pages. Series use an
// Okabe-Ito categorical palette.

import { COLORMAPS } from './colormaps.js';
import { linearTicks, logTicks, formatLinear, formatLog } from './ticks.js';
import { downloadText, downloadBlob, slugify, flash } from './download.js';

const NS = 'http://www.w3.org/2000/svg';

/** Okabe-Ito colour-blind-friendly palette (black replaced by mid grey for dark pages). */
export const PALETTE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#F0E442', '#999999'];

let uid = 0;

function el(tag, attrs, parent, text) {
  const e = document.createElementNS(NS, tag);
  if (attrs) for (const k in attrs) if (attrs[k] !== undefined && attrs[k] !== null) e.setAttribute(k, attrs[k]);
  if (text !== undefined) e.textContent = text;
  if (parent) parent.appendChild(e);
  return e;
}

const fmt = (x) => (Math.round(x * 100) / 100).toString();

export class LinePlot {
  /**
   * @param {HTMLElement} container
   * @param {{width?:number, height?:number, margin?:{l?:number,r?:number,t?:number,b?:number}, interactive?:boolean, selectMode?:string}} [opts]
   */
  constructor(container, opts = {}) {
    this.container = container;
    this.opts = { width: 520, height: 340, selectMode: 'nearest', ...opts };
    this.userMargin = opts.margin || {};
    this.axes = { xlog: false, ylog: false, xrev: false };
    this.series = [];
    this.markers = [];
    this.hmarkers = [];
    this.bands = [];
    this.highlights = [];
    this._onSelect = null;
    this._hlLayer = null;
    this.id = ++uid;
    const { width: W, height: H } = this.opts;
    this.svg = el('svg', {
      viewBox: `0 0 ${W} ${H}`, role: 'img', class: 'qd-plot',
      style: `display:block;width:100%;height:auto;max-width:${W}px;overflow:visible;color:inherit;font:12px/1.2 system-ui,-apple-system,"Segoe UI",sans-serif`,
    });
    container.appendChild(this.svg);
    this.view = { xlim: null, ylim: null, xlog: null, ylog: null };
    this._lastDefaults = null;
    this.interactive = opts.interactive !== false;
    if (this.interactive) this._buildTools();
  }

  /** @param {{xlog?:boolean, ylog?:boolean, xrev?:boolean, xlabel?:string, ylabel?:string, xlim?:number[], ylim?:number[], title?:string, legendPos?:string}} axes */
  setAxes(axes) { this.axes = { ...this.axes, ...axes }; }

  _effAxes() {
    const A = this.axes, v = this.view;
    const d = { xlog: !!A.xlog, ylog: !!A.ylog };
    if (this._lastDefaults) {
      if (this._lastDefaults.xlog !== d.xlog) { v.xlog = null; v.xlim = null; }
      if (this._lastDefaults.ylog !== d.ylog) { v.ylog = null; v.ylim = null; }
    }
    this._lastDefaults = d;
    return { ...A, xlog: v.xlog ?? d.xlog, ylog: v.ylog ?? d.ylog, xlim: v.xlim ?? A.xlim, ylim: v.ylim ?? A.ylim };
  }

  resetView() { this.view = { xlim: null, ylim: null, xlog: null, ylog: null }; this.draw(); }

  /** Register click-to-select callback for linked cross-plot highlighting. */
  onSelect(fn) { this._onSelect = fn; return this; }

  /** Set horizontal reference lines [{y, label, color, dash}]. */
  setHMarkers(hmarkers) { this.hmarkers = hmarkers || []; }

  /** Set shaded background bands [{x0, x1, y0, y1, color, opacity, label}]. */
  setBands(bands) { this.bands = bands || []; }

  /** Set highlighted points [{x, y, color, label, vline, hline}] and fast-update overlay. */
  setHighlights(pts) {
    this.highlights = pts || [];
    this._renderHighlights();
  }

  // --- interaction ----------------------------------------------------------

  _buildTools() {
    const doc = document;
    const bar = this.tools = doc.createElement('div');
    bar.className = 'qd-plot-tools';
    bar.style.maxWidth = `${this.opts.width}px`;
    bar.title = 'Click a point to highlight across linked plots. Wheel: zoom (over an axis: that axis only). Drag: pan. Double-click: reset.';
    const mk = (tag, cls, parent, text) => { const e = doc.createElement(tag); if (cls) e.className = cls; if (text !== undefined) e.textContent = text; parent.appendChild(e); return e; };
    const axisTools = (key, label) => {
      const g = mk('span', 'qd-pt-axis', bar);
      mk('span', 'qd-pt-label', g, label);
      const b = mk('button', 'qd-pt-btn', g, 'lin'); b.type = 'button'; b.title = `${label} axis: linear / logarithmic`;
      b.addEventListener('click', () => this._toggleLog(key));
      const lo = mk('input', 'qd-pt-in', g); lo.type = 'text'; lo.inputMode = 'decimal'; lo.title = `${label} minimum`;
      mk('span', 'qd-pt-dash', g, '–');
      const hi = mk('input', 'qd-pt-in', g); hi.type = 'text'; hi.inputMode = 'decimal'; hi.title = `${label} maximum`;
      const apply = () => this._setRange(key, parseFloat(lo.value), parseFloat(hi.value));
      for (const inp of [lo, hi]) { inp.addEventListener('change', apply); inp.addEventListener('keydown', (ev) => { if (ev.key === 'Enter') { apply(); inp.blur(); } }); }
      return { b, lo, hi };
    };
    this._tx = axisTools('x', 'x');
    this._ty = axisTools('y', 'y');
    const r = mk('button', 'qd-pt-btn qd-pt-reset', bar, 'reset'); r.type = 'button'; r.title = 'Reset zoom, ranges and scales';
    r.addEventListener('click', () => this.resetView());
    const ex = mk('span', 'qd-pt-export', bar);
    const bs = mk('button', 'qd-pt-btn', ex, 'svg'); bs.type = 'button'; bs.title = 'Download this plot as a standalone SVG file';
    bs.addEventListener('click', () => { this.downloadSVG(); flash(bs, 'saved'); });
    const bp = mk('button', 'qd-pt-btn', ex, 'png'); bp.type = 'button'; bp.title = 'Download this plot as a high-resolution PNG image';
    bp.addEventListener('click', () => { this.downloadPNG(); flash(bp, 'saved'); });
    const bc = mk('button', 'qd-pt-btn', ex, 'csv'); bc.type = 'button'; bc.title = 'Download the plotted series as a CSV file (full precision)';
    bc.addEventListener('click', () => { this.downloadCSV(); flash(bc, 'saved'); });
    this._readout = mk('span', 'qd-pt-readout', bar, '');
    this.container.appendChild(bar);
    this._bindEvents();
  }

  _syncTools() {
    if (!this.tools || !this._scales) return;
    const A = this._effAxes(), S = this._scales;
    const f = (v) => Number.isFinite(v) ? String(+v.toPrecision(4)) : '';
    this._tx.b.textContent = A.xlog ? 'log' : 'lin'; this._ty.b.textContent = A.ylog ? 'log' : 'lin';
    if (document.activeElement !== this._tx.lo) this._tx.lo.value = f(S.x0);
    if (document.activeElement !== this._tx.hi) this._tx.hi.value = f(S.x1);
    if (document.activeElement !== this._ty.lo) this._ty.lo.value = f(S.y0);
    if (document.activeElement !== this._ty.hi) this._ty.hi.value = f(S.y1);
    const changed = this.view.xlim || this.view.ylim || this.view.xlog !== null || this.view.ylog !== null;
    this.tools.classList.toggle('is-changed', !!changed);
  }

  _toggleLog(key) {
    const A = this._effAxes();
    const cur = key === 'x' ? A.xlog : A.ylog;
    this.view[key + 'log'] = !cur;
    const lim = this.view[key + 'lim'];
    if (!cur && lim && lim[0] <= 0) this.view[key + 'lim'] = null;
    this.draw();
  }

  _setRange(key, lo, hi) {
    if (!Number.isFinite(lo) || !Number.isFinite(hi) || !(hi > lo)) { this._syncTools(); return; }
    const A = this._effAxes();
    if ((key === 'x' ? A.xlog : A.ylog) && lo <= 0) { this._syncTools(); return; }
    this.view[key + 'lim'] = [lo, hi];
    this.draw();
  }

  _pos(ev) {
    const r = this.svg.getBoundingClientRect();
    const k = this.opts.width / Math.max(1, r.width);
    return { px: (ev.clientX - r.left) * k, py: (ev.clientY - r.top) * k };
  }

  _region(px, py) {
    const S = this._scales; if (!S) return null;
    const { m, pw, ph } = S;
    const inX = px >= m.l && px <= m.l + pw, inY = py >= m.t && py <= m.t + ph;
    if (inX && inY) return 'plot';
    if (inX && py > m.t + ph && py <= this.opts.height) return 'x';
    if (inY && px < m.l && px >= 0) return 'y';
    return null;
  }

  _t(key) {
    const A = this._effAxes(), log = key === 'x' ? A.xlog : A.ylog;
    return log ? { f: Math.log10, g: (u) => Math.pow(10, u) } : { f: (v) => v, g: (u) => u };
  }

  _lim(key) { const S = this._scales, t = this._t(key); return key === 'x' ? [t.f(S.x0), t.f(S.x1)] : [t.f(S.y0), t.f(S.y1)]; }

  _pxToX(px) {
    const S = this._scales, A = this._effAxes();
    const u = (px - S.m.l) / S.pw;
    const frac = A.xrev ? (1 - u) : u;
    const [a, b] = this._lim('x');
    return this._t('x').g(a + frac * (b - a));
  }

  _pyToY(py) {
    const S = this._scales;
    const [a, b] = this._lim('y');
    return this._t('y').g(a + (S.m.t + S.ph - py) / S.ph * (b - a));
  }

  _zoomAxis(key, center, factor) {
    const t = this._t(key), [a, b] = this._lim(key);
    const c = t.f(center);
    const na = c - (c - a) / factor, nb = c + (b - c) / factor;
    if (!(nb > na) || !Number.isFinite(na) || !Number.isFinite(nb)) return;
    if (Math.abs(nb - na) < 1e-9 * (Math.abs(c) + 1e-300)) return;
    this.view[key + 'lim'] = [t.g(na), t.g(nb)];
  }

  _panAxis(key, frac) {
    const t = this._t(key), [a, b] = this._lim(key);
    const d = (b - a) * frac;
    this.view[key + 'lim'] = [t.g(a + d), t.g(b + d)];
  }

  _pickNearest(px, py) {
    const S = this._scales, A = this._effAxes();
    if (!S || !this.series.length) return null;
    const xOnly = this.opts.selectMode === 'x';
    let best = null, bestD = Infinity;
    this.series.forEach((s, si) => {
      if (!s || !s.x || !s.y || s.selectable === false) return;
      const n = Math.min(s.x.length, s.y.length);
      for (let i = 0; i < n; i++) {
        const xv = s.x[i], yv = s.y[i];
        if (!Number.isFinite(xv) || !Number.isFinite(yv) || (A.xlog && xv <= 0) || (A.ylog && yv <= 0)) continue;
        const sx = S.sx(xv), sy = S.sy(yv);
        if (sx < S.m.l - 2 || sx > S.m.l + S.pw + 2 || sy < S.m.t - 2 || sy > S.m.t + S.ph + 2) continue;
        const d = xOnly ? Math.abs(sx - px) : Math.hypot(sx - px, sy - py);
        if (d < bestD) { bestD = d; best = { seriesIndex: si, pointIndex: i, x: xv, y: yv, series: s, dist: d }; }
      }
    });
    return best;
  }

  _bindEvents() {
    const svg = this.svg;
    svg.style.touchAction = 'pan-y';
    svg.addEventListener('wheel', (ev) => {
      const { px, py } = this._pos(ev);
      const reg = this._region(px, py);
      if (!reg || !this._scales) return;
      ev.preventDefault();
      const factor = Math.exp(-Math.sign(ev.deltaY) * Math.min(1, Math.abs(ev.deltaY) / 100) * 0.25);
      const xc = this._pxToX(px), yc = this._pyToY(py);
      if (reg === 'plot' || reg === 'x') this._zoomAxis('x', xc, factor);
      if (reg === 'plot' || reg === 'y') this._zoomAxis('y', yc, factor);
      this.draw();
    }, { passive: false });
    let drag = null;
    svg.addEventListener('pointerdown', (ev) => {
      if (ev.button !== 0) return;
      const { px, py } = this._pos(ev);
      const reg = this._region(px, py);
      if (!reg) return;
      drag = { reg, px, py, startPx: px, startPy: py, moved: false };
      svg.setPointerCapture(ev.pointerId);
    });
    svg.addEventListener('pointermove', (ev) => {
      const { px, py } = this._pos(ev);
      if (drag) {
        const S = this._scales, A = this._effAxes();
        const dx = px - drag.px, dy = py - drag.py;
        if (!drag.moved && Math.hypot(px - drag.startPx, py - drag.startPy) < 4) return;
        drag.moved = true;
        if (drag.reg === 'plot' || drag.reg === 'x') this._panAxis('x', (A.xrev ? 1 : -1) * dx / S.pw);
        if (drag.reg === 'plot' || drag.reg === 'y') this._panAxis('y', dy / S.ph);
        drag.px = px; drag.py = py;
        svg.style.cursor = 'grabbing';
        this.draw();
        return;
      }
      const reg = this._region(px, py);
      svg.style.cursor = reg === 'plot' ? (this._onSelect ? 'pointer' : 'crosshair') : reg ? 'ew-resize' : '';
      if (reg === 'y') svg.style.cursor = 'ns-resize';
      if (this._readout) {
        if (reg === 'plot') {
          const x = this._pxToX(px), y = this._pyToY(py);
          this._readout.textContent = `x = ${+x.toPrecision(4)}, y = ${+y.toPrecision(4)}`;
        } else this._readout.textContent = '';
      }
    });
    const end = (ev) => {
      if (!drag) return;
      const wasClick = !drag.moved && drag.reg === 'plot';
      const { px, py } = drag;
      drag = null; svg.style.cursor = '';
      try { svg.releasePointerCapture(ev.pointerId); } catch (e) { /* ignore */ }
      if (wasClick && this._onSelect) {
        const hit = this._pickNearest(px, py);
        if (hit) this._onSelect(hit);
      }
    };
    svg.addEventListener('pointerup', end);
    svg.addEventListener('pointercancel', end);
    svg.addEventListener('pointerleave', () => { if (!drag && this._readout) this._readout.textContent = ''; });
    svg.addEventListener('dblclick', (ev) => { ev.preventDefault(); this.resetView(); });
  }

  setSeries(series) { this.series = series; }
  setMarkers(markers) { this.markers = markers; }

  toSVGString() { return new XMLSerializer().serializeToString(this.svg); }

  toStandaloneSVG() {
    const W = this.opts.width, H = this.opts.height;
    const c = this.svg.cloneNode(true);
    c.setAttribute('xmlns', NS);
    c.setAttribute('width', W); c.setAttribute('height', H);
    c.setAttribute('style', 'color:#222;overflow:visible;font:12px/1.2 system-ui,-apple-system,"Segoe UI",sans-serif');
    const bg = document.createElementNS(NS, 'rect');
    for (const [k, v] of Object.entries({ x: 0, y: 0, width: W, height: H, fill: '#fff' })) bg.setAttribute(k, v);
    c.insertBefore(bg, c.firstChild);
    let out = new XMLSerializer().serializeToString(c);
    if (!/^<svg[^>]*\sxmlns=/.test(out)) out = out.replace(/^<svg/, `<svg xmlns="${NS}"`);
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + out + '\n';
  }

  _fileStem() { return slugify(this.axes.title || this.axes.ylabel || 'plot'); }

  downloadSVG() { return downloadText(this.toStandaloneSVG(), this._fileStem() + '.svg', 'image/svg+xml'); }

  /** Export crisp 2× high-DPI PNG image rendered from the standalone SVG. */
  downloadPNG(scale = 2) {
    const W = this.opts.width, H = this.opts.height;
    const svgStr = this.toStandaloneSVG();
    const blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const img = new Image();
    img.onload = () => {
      const cv = document.createElement('canvas');
      cv.width = Math.round(W * scale); cv.height = Math.round(H * scale);
      const ctx = cv.getContext('2d');
      ctx.fillStyle = '#ffffff';
      ctx.fillRect(0, 0, cv.width, cv.height);
      ctx.scale(scale, scale);
      ctx.drawImage(img, 0, 0, W, H);
      URL.revokeObjectURL(url);
      cv.toBlob((b) => { if (b) downloadBlob(b, this._fileStem() + '.png'); }, 'image/png');
    };
    img.onerror = () => { URL.revokeObjectURL(url); };
    img.src = url;
    return this._fileStem() + '.png';
  }

  toCSV() {
    const S = this.series.filter((s) => s && s.x && s.y);
    const name = (s, k) => (s.label ? String(s.label) : `series${k + 1}`);
    const q = (t) => (/[",\n\r]/.test(t) ? '"' + t.replace(/"/g, '""') + '"' : t);
    const num = (v) => String(v);
    const n0 = S.length ? Math.min(S[0].x.length, S[0].y.length) : 0;
    const shared = S.length > 0 && S.every((s) => {
      if (Math.min(s.x.length, s.y.length) !== n0 || s.x.length !== S[0].x.length) return false;
      for (let i = 0; i < n0; i++) if (s.x[i] !== S[0].x[i] && !(s.x[i] !== s.x[i] && S[0].x[i] !== S[0].x[i])) return false;
      return true;
    });
    const L = [];
    if (shared) {
      L.push(['x', ...S.map(name)].map(q).join(','));
      for (let i = 0; i < n0; i++) L.push([num(S[0].x[i]), ...S.map((s) => num(s.y[i]))].join(','));
    } else {
      L.push('series,x,y');
      S.forEach((s, k) => { const nm = q(name(s, k)), n = Math.min(s.x.length, s.y.length); for (let i = 0; i < n; i++) L.push(`${nm},${num(s.x[i])},${num(s.y[i])}`); });
    }
    return L.join('\n') + '\n';
  }

  downloadCSV() { return downloadText(this.toCSV(), this._fileStem() + '.csv', 'text/csv'); }

  destroy() { this.svg.remove(); if (this.tools) this.tools.remove(); }

  _drawSeries() { return this.series; }
  _showYTicks() { return true; }
  _defaultMargin() { return {}; }
  _decorate() {}

  _dataRange(series, key, log) {
    let lo = Infinity, hi = -Infinity;
    for (const s of series) {
      const a = s[key];
      for (let i = 0; i < a.length; i++) {
        const v = a[i];
        if (!Number.isFinite(v) || (log && v <= 0)) continue;
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
    }
    return [lo, hi];
  }

  _limits(series, key, log, user, pad) {
    let [lo, hi] = this._dataRange(series, key, log);
    if (!(lo <= hi)) { lo = log ? 1 : 0; hi = log ? 10 : 1; }
    if (lo === hi) { if (log) { lo /= 2; hi *= 2; } else { const d = Math.abs(lo) || 1; lo -= d / 2; hi += d / 2; } }
    if (log) { const e = Math.pow(hi / lo, pad); lo /= e; hi *= e; }
    else { const d = (hi - lo) * pad; lo -= d; hi += d; }
    if (user && user.length === 2) {
      if (Number.isFinite(user[0]) && !(log && user[0] <= 0)) lo = user[0];
      if (Number.isFinite(user[1]) && !(log && user[1] <= 0)) hi = user[1];
    }
    return [lo, hi];
  }

  _renderHighlights() {
    if (!this._hlLayer || !this._scales) return;
    const g = this._hlLayer, S = this._scales, A = this._effAxes();
    while (g.firstChild) g.removeChild(g.firstChild);
    for (const pt of this.highlights) {
      if (!pt) continue;
      const col = pt.color || '#e03131';
      const hasX = Number.isFinite(pt.x) && !(A.xlog && pt.x <= 0);
      const hasY = Number.isFinite(pt.y) && !(A.ylog && pt.y <= 0);
      const px = hasX ? S.sx(pt.x) : NaN;
      const py = hasY ? S.sy(pt.y) : NaN;
      if (hasX && px >= S.m.l && px <= S.m.l + S.pw && pt.vline !== false) {
        el('line', { x1: fmt(px), x2: fmt(px), y1: S.m.t, y2: S.m.t + S.ph, stroke: col, 'stroke-width': 1.2, 'stroke-dasharray': '4 3', 'stroke-opacity': 0.75 }, g);
      }
      if (hasY && py >= S.m.t && py <= S.m.t + S.ph && pt.hline) {
        el('line', { x1: S.m.l, x2: S.m.l + S.pw, y1: fmt(py), y2: fmt(py), stroke: col, 'stroke-width': 1.1, 'stroke-dasharray': '4 3', 'stroke-opacity': 0.65 }, g);
      }
      if (hasX && hasY && px >= S.m.l && px <= S.m.l + S.pw && py >= S.m.t && py <= S.m.t + S.ph) {
        el('circle', { cx: fmt(px), cy: fmt(py), r: 7, fill: col, 'fill-opacity': 0.22, stroke: col, 'stroke-width': 1.6 }, g);
        el('circle', { cx: fmt(px), cy: fmt(py), r: 3.2, fill: col, stroke: '#fff', 'stroke-width': 1.2 }, g);
        if (pt.label) {
          const right = px > S.m.l + S.pw * 0.65;
          el('text', { x: fmt(px + (right ? -10 : 10)), y: fmt(Math.max(S.m.t + 14, py - 8)), 'text-anchor': right ? 'end' : 'start', fill: col, 'font-size': 11, 'font-weight': 600 }, g, pt.label);
        }
      }
    }
  }

  draw() {
    const A = this._effAxes(), svg = this.svg, W = this.opts.width, H = this.opts.height;
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    svg.setAttribute('aria-label', A.title || (A.ylabel && A.xlabel ? `${A.ylabel} versus ${A.xlabel}` : 'plot'));
    const series = this._drawSeries();
    const [x0, x1] = this._limits(series, 'x', A.xlog, A.xlim, 0);
    const [y0, y1] = this._limits(series, 'y', A.ylog, A.ylim, A.ylog ? 0.03 : 0.05);

    const xt = A.xlog ? this._logTicks(x0, x1) : this._linTicks(x0, x1, Math.max(3, Math.floor(this.opts.width / 90)));
    const yt = A.ylog ? this._logTicks(y0, y1) : this._linTicks(y0, y1, Math.max(3, Math.floor(this.opts.height / 55)));
    const maxYLabel = this._showYTicks() ? Math.max(0, ...yt.labels.map((s) => s.length)) : 0;
    const dm = this._defaultMargin();
    const m = {
      l: this.userMargin.l ?? (dm.l ?? (14 * (A.ylabel ? 1 : 0) + 14 + maxYLabel * 6.6 + 6)),
      r: this.userMargin.r ?? (dm.r ?? 16),
      t: this.userMargin.t ?? (A.title ? 30 : 12),
      b: this.userMargin.b ?? (A.xlabel ? 46 : 30),
    };
    const pw = W - m.l - m.r, ph = H - m.t - m.b;
    const lx = A.xlog ? [Math.log10(x0), Math.log10(x1)] : [x0, x1];
    const ly = A.ylog ? [Math.log10(y0), Math.log10(y1)] : [y0, y1];
    const sx = (v) => {
      const f = ((A.xlog ? Math.log10(v) : v) - lx[0]) / (lx[1] - lx[0]);
      return m.l + (A.xrev ? (1 - f) : f) * pw;
    };
    const sy = (v) => m.t + ph - ((A.ylog ? Math.log10(v) : v) - ly[0]) / (ly[1] - ly[0]) * ph;
    this._scales = { sx, sy, m, pw, ph, x0, x1, y0, y1 };

    const clipId = `qd-clip-${this.id}`;
    const defs = el('defs', null, svg);
    el('rect', { x: m.l, y: m.t, width: pw, height: ph }, el('clipPath', { id: clipId }, defs));

    // shaded bands (behind grid & series)
    const plot = el('g', { 'clip-path': `url(#${clipId})` }, svg);
    for (const bd of this.bands) {
      let rx = m.l, rw = pw, ry = m.t, rh = ph;
      if (bd.x0 !== undefined && bd.x1 !== undefined) {
        if (A.xlog && (bd.x0 <= 0 || bd.x1 <= 0)) continue;
        const a = sx(bd.x0), b = sx(bd.x1);
        rx = Math.min(a, b); rw = Math.abs(b - a);
      }
      if (bd.y0 !== undefined && bd.y1 !== undefined) {
        if (A.ylog && (bd.y0 <= 0 || bd.y1 <= 0)) continue;
        const a = sy(bd.y0), b = sy(bd.y1);
        ry = Math.min(a, b); rh = Math.abs(b - a);
      }
      el('rect', { x: fmt(rx), y: fmt(ry), width: fmt(rw), height: fmt(rh), fill: bd.color || '#56B4E9', 'fill-opacity': bd.opacity ?? 0.16 }, plot);
      if (bd.label) el('text', { x: fmt(rx + 6), y: fmt(ry + 14), fill: bd.color || '#0072B2', 'font-size': 10.5, 'font-weight': 600, opacity: 0.9 }, plot, bd.label);
    }

    // grid + ticks
    const grid = el('g', { stroke: 'currentColor', 'stroke-opacity': 0.12, 'stroke-width': 1, fill: 'none' }, svg);
    const ticks = el('g', { stroke: 'currentColor', 'stroke-width': 1 }, svg);
    const labels = el('g', { fill: 'currentColor' }, svg);
    xt.major.forEach((v, k) => {
      const x = Math.round(sx(v)) + 0.5;
      if (x < m.l - 0.5 || x > m.l + pw + 0.5) return;
      el('line', { x1: x, x2: x, y1: m.t, y2: m.t + ph }, grid);
      el('line', { x1: x, x2: x, y1: m.t + ph, y2: m.t + ph + 5 }, ticks);
      el('text', { x, y: m.t + ph + 18, 'text-anchor': 'middle', stroke: 'none' }, labels, xt.labels[k]);
    });
    xt.minor.forEach((v) => { const x = Math.round(sx(v)) + 0.5; if (x >= m.l && x <= m.l + pw) el('line', { x1: x, x2: x, y1: m.t + ph, y2: m.t + ph + 3, 'stroke-opacity': 0.6 }, ticks); });
    if (this._showYTicks()) {
      yt.major.forEach((v, k) => {
        const y = Math.round(sy(v)) + 0.5;
        if (y < m.t - 0.5 || y > m.t + ph + 0.5) return;
        el('line', { x1: m.l, x2: m.l + pw, y1: y, y2: y }, grid);
        el('line', { x1: m.l - 5, x2: m.l, y1: y, y2: y }, ticks);
        el('text', { x: m.l - 8, y: y + 4, 'text-anchor': 'end', stroke: 'none' }, labels, yt.labels[k]);
      });
      yt.minor.forEach((v) => { const y = Math.round(sy(v)) + 0.5; if (y >= m.t && y <= m.t + ph) el('line', { x1: m.l - 3, x2: m.l, y1: y, y2: y, 'stroke-opacity': 0.6 }, ticks); });
    }

    // frame
    el('rect', { x: m.l + 0.5, y: m.t + 0.5, width: pw, height: ph, fill: 'none', stroke: 'currentColor', 'stroke-width': 1 }, svg);

    // series layer (re-append above grid)
    const seriesLayer = el('g', { 'clip-path': `url(#${clipId})` }, svg);
    this._decorate(seriesLayer, { sx, sy, m, pw, ph });
    series.forEach((s, k) => {
      const color = s.color || PALETTE[k % PALETTE.length];
      const op = s.opacity ?? 1;
      const n = Math.min(s.x.length, s.y.length);
      const ok = (i) => {
        const x = s.x[i], y = s.y[i];
        return Number.isFinite(x) && Number.isFinite(y) && !(A.xlog && x <= 0) && !(A.ylog && y <= 0);
      };
      if (s.line !== false) {
        let d = '', pen = false;
        for (let i = 0; i < n; i++) {
          if (!ok(i)) { pen = false; continue; }
          d += (pen ? 'L' : 'M') + fmt(sx(s.x[i])) + ' ' + fmt(sy(s.y[i]));
          pen = true;
        }
        if (d) el('path', { d, fill: 'none', stroke: color, 'stroke-width': s.width ?? 1.6, 'stroke-dasharray': s.dash, 'stroke-opacity': op, 'stroke-linejoin': 'round', 'stroke-linecap': 'round' }, seriesLayer);
      }
      if (s.points) {
        const g = el('g', s.hollow ? { fill: 'none', stroke: color, 'stroke-width': 1.5, 'stroke-opacity': op } : { fill: color, 'fill-opacity': op }, seriesLayer);
        for (let i = 0; i < n; i++) if (ok(i)) el('circle', { cx: fmt(sx(s.x[i])), cy: fmt(sy(s.y[i])), r: s.radius ?? 2.6 }, g);
      }
    });

    // vertical markers
    this.markers.forEach((mk, k) => {
      if (!Number.isFinite(mk.x) || (A.xlog && mk.x <= 0)) return;
      const x = sx(mk.x);
      if (x < m.l || x > m.l + pw) return;
      const color = mk.color || 'currentColor';
      el('line', { x1: x, x2: x, y1: m.t, y2: m.t + ph, stroke: color, 'stroke-dasharray': mk.dash || '5 3', 'stroke-width': 1.2, 'stroke-opacity': 0.8 }, svg);
      if (mk.label) {
        const t = el('text', { x: x + 4, y: m.t + 12 + 12 * (k % 2), fill: color, 'font-size': 11 }, svg, mk.label);
        if (x > m.l + pw - 80) { t.setAttribute('x', x - 4); t.setAttribute('text-anchor', 'end'); }
      }
    });

    // horizontal markers
    this.hmarkers.forEach((mk) => {
      if (!Number.isFinite(mk.y) || (A.ylog && mk.y <= 0)) return;
      const y = sy(mk.y);
      if (y < m.t || y > m.t + ph) return;
      const color = mk.color || 'currentColor';
      el('line', { x1: m.l, x2: m.l + pw, y1: y, y2: y, stroke: color, 'stroke-dasharray': mk.dash || '4 3', 'stroke-width': 1.1, 'stroke-opacity': 0.75 }, svg);
      if (mk.label) el('text', { x: m.l + 5, y: y - 4, fill: color, 'font-size': 10.5 }, svg, mk.label);
    });

    // highlight overlay layer
    this._hlLayer = el('g', { class: 'qd-hl-layer' }, svg);
    this._renderHighlights();

    // titles / axis labels
    if (A.title) el('text', { x: m.l + pw / 2, y: 18, 'text-anchor': 'middle', 'font-size': 13.5, 'font-weight': 600, fill: 'currentColor' }, svg, A.title);
    if (A.xlabel) el('text', { x: m.l + pw / 2, y: H - 6, 'text-anchor': 'middle', fill: 'currentColor' }, svg, A.xlabel);
    if (A.ylabel) el('text', { transform: `translate(14 ${m.t + ph / 2}) rotate(-90)`, 'text-anchor': 'middle', fill: 'currentColor' }, svg, A.ylabel);

    this._legend(series, m, pw, ph, A.legendPos);
    this._syncTools();
  }

  _legend(series, m, pw, ph, pos = 'tr') {
    const entries = series.map((s, k) => ({ s, k })).filter((e) => e.s.label);
    if (!entries.length) return;
    const w = Math.max(...entries.map((e) => e.s.label.length)) * 6.3 + 34;
    const hgt = entries.length * 15 + 8;
    const g = el('g', { 'font-size': 11.2 }, this.svg);
    const x0 = pos.includes('l') ? m.l + 8 : m.l + pw - w - 6;
    const y0 = pos.includes('b') ? m.t + ph - hgt - 6 : m.t + 6;
    el('rect', { x: x0 - 4, y: y0, width: w + 6, height: hgt, rx: 4, fill: 'var(--qd-legend-bg, rgba(128,128,128,0.08))', stroke: 'currentColor', 'stroke-opacity': 0.14 }, g);
    entries.forEach((e, row) => {
      const y = y0 + 14 + row * 15, color = e.s.color || PALETTE[e.k % PALETTE.length];
      if (e.s.line !== false) el('line', { x1: x0, x2: x0 + 18, y1: y - 4, y2: y - 4, stroke: color, 'stroke-width': e.s.width ?? 1.6, 'stroke-dasharray': e.s.dash, 'stroke-opacity': e.s.opacity ?? 1 }, g);
      if (e.s.points) el('circle', e.s.hollow ? { cx: x0 + 9, cy: y - 4, r: 2.5, fill: 'none', stroke: color, 'stroke-width': 1.5 } : { cx: x0 + 9, cy: y - 4, r: 2.5, fill: color }, g);
      el('text', { x: x0 + 24, y, fill: 'currentColor' }, g, e.s.label);
    });
  }

  _linTicks(lo, hi, target) {
    const t = linearTicks(lo, hi, target);
    return { major: t.ticks, minor: [], labels: t.ticks.map((v) => formatLinear(v, t.step)) };
  }

  _logTicks(lo, hi) {
    const t = logTicks(lo, hi);
    return { major: t.major, minor: t.minor, labels: t.major.map(formatLog) };
  }
}

/**
 * Interactive 2D scalar heatmap (e.g. e± pair-instability Γ₁(ρ,T) valley or Kippenhahn map)
 * with click-to-select, PNG download, and optional overlay contours/tracks.
 */
export class HeatmapPlot {
  constructor(container, opts = {}) {
    this.container = container;
    this.opts = { width: 520, height: 340, colormap: 'viridis', ...opts };
    this.axes = { xlabel: '', ylabel: '', clabel: '', title: '' };
    this.grid = null;
    this.overlays = [];
    this.highlight = null;
    this._onSelect = null;
    this.wrap = document.createElement('div');
    this.wrap.style.maxWidth = `${this.opts.width}px`;
    this.canvas = document.createElement('canvas');
    this.canvas.width = this.opts.width * 2;
    this.canvas.height = this.opts.height * 2;
    this.canvas.style.cssText = 'display:block;width:100%;height:auto;cursor:crosshair;border-radius:4px';
    this.wrap.appendChild(this.canvas);
    const bar = this.tools = document.createElement('div');
    bar.className = 'qd-plot-tools';
    const bp = document.createElement('button'); bp.type = 'button'; bp.className = 'qd-pt-btn'; bp.textContent = 'png';
    bp.title = 'Download heatmap as high-resolution PNG';
    bp.addEventListener('click', () => { this.downloadPNG(); flash(bp, 'saved'); });
    bar.appendChild(bp);
    this._readout = document.createElement('span'); this._readout.className = 'qd-pt-readout';
    bar.appendChild(this._readout);
    this.wrap.appendChild(bar);
    container.appendChild(this.wrap);
    this._bindEvents();
  }

  setAxes(axes) { this.axes = { ...this.axes, ...axes }; }
  setGrid(grid) { this.grid = grid; } // { xs, ys, z (ny*nx), zmin, zmax, contour }
  setOverlays(overlays) { this.overlays = overlays || []; }
  onSelect(fn) { this._onSelect = fn; return this; }
  setHighlight(pt) { this.highlight = pt; this.draw(); }

  _geom() {
    const W = this.opts.width, H = this.opts.height;
    const m = { l: 56, r: 62, t: this.axes.title ? 28 : 12, b: 42 };
    return { W, H, m, pw: W - m.l - m.r, ph: H - m.t - m.b };
  }

  _bindEvents() {
    const pos = (ev) => {
      const r = this.canvas.getBoundingClientRect();
      return { px: (ev.clientX - r.left) * this.opts.width / Math.max(1, r.width), py: (ev.clientY - r.top) * this.opts.height / Math.max(1, r.height) };
    };
    this.canvas.addEventListener('pointermove', (ev) => {
      if (!this.grid) return;
      const { px, py } = pos(ev), { m, pw, ph } = this._geom(), G = this.grid;
      if (px < m.l || px > m.l + pw || py < m.t || py > m.t + ph) { this._readout.textContent = ''; return; }
      const fx = (px - m.l) / pw, fy = 1 - (py - m.t) / ph;
      const ix = Math.min(G.xs.length - 1, Math.max(0, Math.round(fx * (G.xs.length - 1))));
      const iy = Math.min(G.ys.length - 1, Math.max(0, Math.round(fy * (G.ys.length - 1))));
      const zv = G.z[iy * G.xs.length + ix];
      this._readout.textContent = `x=${+G.xs[ix].toPrecision(4)}, y=${+G.ys[iy].toPrecision(4)}, z=${Number.isFinite(zv) ? +zv.toPrecision(4) : '–'}`;
    });
    this.canvas.addEventListener('click', (ev) => {
      if (!this.grid || !this._onSelect) return;
      const { px, py } = pos(ev), { m, pw, ph } = this._geom(), G = this.grid;
      if (px < m.l || px > m.l + pw || py < m.t || py > m.t + ph) return;
      const fx = (px - m.l) / pw, fy = 1 - (py - m.t) / ph;
      const ix = Math.min(G.xs.length - 1, Math.max(0, Math.round(fx * (G.xs.length - 1))));
      const iy = Math.min(G.ys.length - 1, Math.max(0, Math.round(fy * (G.ys.length - 1))));
      this._onSelect({ ix, iy, x: G.xs[ix], y: G.ys[iy], z: G.z[iy * G.xs.length + ix] });
    });
  }

  downloadPNG() {
    const name = slugify(this.axes.title || 'heatmap') + '.png';
    this.canvas.toBlob((b) => { if (b) downloadBlob(b, name); }, 'image/png');
    return name;
  }

  draw() {
    if (!this.grid) return;
    const G = this.grid, { W, H, m, pw, ph } = this._geom();
    const ctx = this.canvas.getContext('2d');
    const dark = document.documentElement.classList.contains('quarto-dark') || document.documentElement.dataset.bsTheme === 'dark';
    const fg = dark ? '#e6e8eb' : '#1f2328', bg = dark ? '#181b20' : '#ffffff';
    ctx.save(); ctx.setTransform(2, 0, 0, 2, 0, 0);
    ctx.fillStyle = bg; ctx.fillRect(0, 0, W, H);

    const nx = G.xs.length, ny = G.ys.length;
    const x0 = G.xs[0], x1 = G.xs[nx - 1], y0 = G.ys[0], y1 = G.ys[ny - 1];
    const zmin = G.zmin ?? Math.min(...G.z), zmax = G.zmax ?? Math.max(...G.z);
    const cm = COLORMAPS[this.opts.colormap] || COLORMAPS.viridis;
    const cw = pw / Math.max(1, nx - 1), ch = ph / Math.max(1, ny - 1);
    for (let iy = 0; iy < ny; iy++) {
      const py = m.t + ph - ((iy + 0.5) / ny) * ph;
      for (let ix = 0; ix < nx; ix++) {
        const v = G.z[iy * nx + ix];
        const t = Math.max(0, Math.min(1, (v - zmin) / Math.max(1e-30, zmax - zmin)));
        const c = cm(t);
        ctx.fillStyle = `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`;
        ctx.fillRect(m.l + (ix / nx) * pw, m.t + ph - ((iy + 1) / ny) * ph, pw / nx + 0.6, ph / ny + 0.6);
      }
    }
    const sx = (v) => m.l + (v - x0) / (x1 - x0) * pw;
    const sy = (v) => m.t + ph - (v - y0) / (y1 - y0) * ph;

    // Optional threshold contour highlight (e.g. Γ₁ = 4/3)
    if (Number.isFinite(G.contour)) {
      ctx.fillStyle = 'rgba(255,255,255,0.85)';
      for (let iy = 0; iy < ny - 1; iy++) {
        for (let ix = 0; ix < nx - 1; ix++) {
          const a = G.z[iy * nx + ix] < G.contour, b = G.z[iy * nx + ix + 1] < G.contour, c = G.z[(iy + 1) * nx + ix] < G.contour;
          if (a !== b || a !== c) ctx.fillRect(sx(G.xs[ix]) - 1, sy(G.ys[iy]) - 1, 2.2, 2.2);
        }
      }
    }

    // Overlays
    ctx.save(); ctx.beginPath(); ctx.rect(m.l, m.t, pw, ph); ctx.clip();
    for (const ov of this.overlays) {
      ctx.strokeStyle = ov.color || '#ffffff'; ctx.lineWidth = ov.width || 2;
      ctx.setLineDash(ov.dash || []);
      ctx.beginPath();
      for (let i = 0; i < ov.x.length; i++) {
        const px = sx(ov.x[i]), py = sy(ov.y[i]);
        if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      }
      ctx.stroke();
      ctx.setLineDash([]);
      if (ov.label && ov.x.length) {
        const k = Math.floor(ov.x.length * 0.72);
        ctx.font = '600 10.5px system-ui,sans-serif'; ctx.fillStyle = ov.color || '#ffffff';
        ctx.fillText(ov.label, sx(ov.x[k]) + 4, sy(ov.y[k]) - 4);
      }
    }
    if (this.highlight) {
      const hx = sx(this.highlight.x), hy = sy(this.highlight.y);
      ctx.strokeStyle = '#ff6b6b'; ctx.lineWidth = 1.4; ctx.setLineDash([4, 3]);
      ctx.beginPath(); ctx.moveTo(hx, m.t); ctx.lineTo(hx, m.t + ph); ctx.moveTo(m.l, hy); ctx.lineTo(m.l + pw, hy); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = '#ff6b6b'; ctx.beginPath(); ctx.arc(hx, hy, 4.5, 0, Math.PI * 2); ctx.fill();
      ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5; ctx.stroke();
    }
    ctx.restore();

    // Colorbar
    const cbx = m.l + pw + 10, cbw = 12;
    for (let j = 0; j < ph; j++) {
      const c = cm(1 - j / Math.max(1, ph - 1));
      ctx.fillStyle = `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`;
      ctx.fillRect(cbx, m.t + j, cbw, 1.5);
    }
    ctx.strokeStyle = fg; ctx.lineWidth = 0.8; ctx.strokeRect(cbx, m.t, cbw, ph);
    ctx.fillStyle = fg; ctx.font = '10px system-ui,sans-serif'; ctx.textAlign = 'left';
    for (let k = 0; k <= 4; k++) {
      const v = zmin + (zmax - zmin) * k / 4;
      ctx.fillText(String(+v.toPrecision(3)), cbx + cbw + 4, m.t + ph - (k / 4) * ph + 3);
    }

    // Axes & ticks
    ctx.strokeStyle = fg; ctx.lineWidth = 1; ctx.strokeRect(m.l + 0.5, m.t + 0.5, pw, ph);
    ctx.fillStyle = fg; ctx.font = '11px system-ui,sans-serif';
    const xt = linearTicks(x0, x1, 6);
    ctx.textAlign = 'center';
    for (const v of xt.ticks) {
      const x = sx(v);
      ctx.beginPath(); ctx.moveTo(x, m.t + ph); ctx.lineTo(x, m.t + ph + 4); ctx.stroke();
      ctx.fillText(formatLinear(v, xt.step), x, m.t + ph + 16);
    }
    const yt = linearTicks(y0, y1, 5);
    ctx.textAlign = 'right';
    for (const v of yt.ticks) {
      const y = sy(v);
      ctx.beginPath(); ctx.moveTo(m.l - 4, y); ctx.lineTo(m.l, y); ctx.stroke();
      ctx.fillText(formatLinear(v, yt.step), m.l - 7, y + 3.5);
    }
    ctx.textAlign = 'center';
    if (this.axes.title) { ctx.font = '600 13px system-ui,sans-serif'; ctx.fillText(this.axes.title, m.l + pw / 2, 17); }
    ctx.font = '11.5px system-ui,sans-serif';
    if (this.axes.xlabel) ctx.fillText(this.axes.xlabel, m.l + pw / 2, H - 6);
    if (this.axes.ylabel) { ctx.save(); ctx.translate(14, m.t + ph / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(this.axes.ylabel, 0, 0); ctx.restore(); }
    ctx.restore();
  }
}

export class WaterfallPlot extends LinePlot {
  constructor(container, opts = {}) {
    super(container, opts);
    this.offset = opts.offset;
    this.colormap = opts.colormap || 'viridis';
  }
  _defaultMargin() { return { l: this.axes.ylabel ? 28 : 12, r: 60 }; }
  _showYTicks() { return false; }
  _amplitude() {
    let amp = 0;
    for (const s of this.series) {
      const [lo, hi] = this._dataRange([s], 'y', false);
      if (hi - lo > amp) amp = hi - lo;
    }
    return amp || 1;
  }
  _drawSeries() {
    const K = this.series.length;
    const off = this.offset ?? 0.6 * this._amplitude();
    this._off = off;
    const cm = COLORMAPS[this.colormap] || COLORMAPS.viridis;
    return this.series.map((s, k) => {
      const c = cm(K > 1 ? 0.15 + 0.7 * k / (K - 1) : 0.5);
      const y = new Float64Array(s.y.length);
      for (let i = 0; i < y.length; i++) y[i] = s.y[i] + k * off;
      return { ...s, y, base: k * off, label: undefined, _label: s.label, color: s.color || `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})` };
    });
  }
  _decorate(plot, { sx, sy, m, pw }) {
    const series = this._drawSeries();
    series.forEach((s) => {
      const y = sy(s.base);
      el('line', { x1: m.l, x2: m.l + pw, y1: y, y2: y, stroke: 'currentColor', 'stroke-opacity': 0.12, 'stroke-dasharray': '2 3' }, plot);
    });
    const g = el('g', { fill: 'currentColor', 'font-size': 11 }, this.svg);
    series.forEach((s) => {
      if (s._label) el('text', { x: m.l + pw + 6, y: sy(s.base) + 4, fill: s.color }, g, s._label);
    });
  }
}
