// 2D mechanism animation player.
//
// Consumes an AnimationSpec (see backend/animation.py) and renders it as an
// animated, ball-and-stick SVG. The spec is a list of keyframes; this player
// interpolates between consecutive frames:
//   - an atom id present in both frames slides from one position to the other;
//   - an id in only one frame fades in or out.
// That one rule renders the atom-tracked SN2 inversion and the atom-mapped
// morph mechanisms alike.
//
// Electron behaviour is the point of the visualization, so it is drawn
// explicitly rather than implied by lines: lone pairs (counted from formal
// charge) sit on each atom, curly arrows carry flowing electron-pair dots from
// source to destination, and forming/breaking bonds show electrons streaming
// in or out. None of this is in the spec — it is derived deterministically
// from the atoms and bonds, so every reaction gets it for free.
//
// Interactivity is deliberately limited to the *camera, playback, and
// inspection* — scrub, step between phases, zoom/pan, hover for atom info.
// Every visible state is a deterministic interpolation along the fixed
// mechanism path, so the viewer can explore but can never drag the chemistry
// into an unrealistic configuration.

const SVG_NS = "http://www.w3.org/2000/svg";

const ELEMENT_COLORS = {
  H: "#e8e8e8", C: "#9aa0a6", N: "#5577ff", O: "#ff4d4d", F: "#7ee6a8",
  Cl: "#3cc94f", Br: "#c1622f", I: "#a64dff", S: "#ffd633", P: "#ff9f40",
  Na: "#8866ff", K: "#8f4bd6", Ag: "#c0c0c8", Mg: "#3fd6a0", Ca: "#4fd0e0",
};
const ELEMENT_RADII = {
  H: 9, C: 15, N: 14, O: 14, F: 13, Cl: 17, Br: 18, I: 20, S: 16, P: 16,
};
const ELEMENT_NAMES = {
  H: "Hydrogen", C: "Carbon", N: "Nitrogen", O: "Oxygen", F: "Fluorine",
  Cl: "Chlorine", Br: "Bromine", I: "Iodine", S: "Sulfur", P: "Phosphorus",
  Na: "Sodium", K: "Potassium", Ag: "Silver", Mg: "Magnesium", Ca: "Calcium",
};
// Group (valence) electron counts, for lone-pair bookkeeping.
const GROUP_VE = {
  H: 1, B: 3, C: 4, N: 5, O: 6, F: 7, Si: 4, P: 5, S: 6, Cl: 7,
  Br: 7, I: 7, Na: 1, K: 1, Mg: 2, Ca: 2, Ag: 1, Zn: 2,
};
const ROLE_COLORS = {
  nucleophile: "#2fd0ff", electrophile: "#ff9f43", leaving_group: "#ff5e5e",
  acid_proton: "#ffd633", base_site: "#2fd0ff",
};
const ROLE_LABELS = {
  nucleophile: "nucleophile", electrophile: "electrophile",
  leaving_group: "leaving group", acid_proton: "acidic proton",
  base_site: "basic site",
};
// One visual language for "electron": this cyan is used for lone pairs and for
// every flowing electron dot, so the eye tracks charge as it moves.
const ELECTRON_COLOR = "#8fe8ff";

const colorFor = (el) => ELEMENT_COLORS[el] || "#ff80c0";
const radiusFor = (el) => ELEMENT_RADII[el] || 14;
const lerp = (a, b, u) => a + (b - a) * u;
const smoothstep = (u) => u * u * (3 - 2 * u);
const bondKey = (a, b) => (a < b ? `${a}|${b}` : `${b}|${a}`);

// Nonbonding electron pairs on an atom, from formal charge: FC = V - B - N,
// so N (nonbonding electrons) = V - B - FC. Returns whole lone pairs.
function lonePairCount(el, charge, orderSum) {
  const v = GROUP_VE[el];
  if (v === undefined) return 0;
  const nonbonding = v - orderSum - (charge || 0);
  if (nonbonding < 1) return 0;
  return Math.min(4, Math.floor(nonbonding / 2));
}

// Place `count` lone pairs in the widest angular gaps left by the bonded
// neighbours, so they sit where real lone pairs would — pointing away from
// bonds (e.g. a nucleophile's pair aims straight at the carbon it attacks).
function lonePairAngles(atom, neighbors, count) {
  const occupied = neighbors.map((n) => Math.atan2(n.y - atom.y, n.x - atom.x));
  const out = [];
  for (let k = 0; k < count; k++) {
    let bestAngle = (k / Math.max(count, 1)) * 2 * Math.PI;
    if (occupied.length) {
      const sorted = [...occupied].sort((a, b) => a - b);
      let bestGap = -1;
      for (let i = 0; i < sorted.length; i++) {
        const a = sorted[i];
        const b = i === sorted.length - 1 ? sorted[0] + 2 * Math.PI : sorted[i + 1];
        const gap = b - a;
        if (gap > bestGap) { bestGap = gap; bestAngle = a + gap / 2; }
      }
    }
    out.push(bestAngle);
    occupied.push(bestAngle);
  }
  return out;
}

function el(tag, attrs = {}) {
  const node = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  return node;
}

// Shade a hex color: amt > 0 lightens toward white, amt < 0 darkens toward
// black. Used to fake a lit sphere (bright highlight, darker rim).
function shade(hex, amt) {
  const n = parseInt(hex.slice(1), 16);
  const r = (n >> 16) & 255, g = (n >> 8) & 255, b = n & 255;
  const target = amt >= 0 ? 255 : 0;
  const a = Math.abs(amt);
  const mix = (c) => Math.round(c + (target - c) * a);
  return `rgb(${mix(r)},${mix(g)},${mix(b)})`;
}

export class MechanismPlayer {
  constructor(svg, { caption, equation, tooltip } = {}) {
    this.svg = svg;
    this.captionEl = caption;
    this.equationEl = equation;
    this.tooltipEl = tooltip;
    this.spec = null;
    this.frames = [];
    this.progress = 0;
    this.playing = false;
    this.speed = 1;
    this.loop = false;
    this._raf = null;
    this._last = null;
    this._clock = 0; // wall-clock seconds, for pulse/flow effects
    this._base = null; // base viewBox [x,y,w,h]
    this.view = { scale: 1, tx: 0, ty: 0 };
    this._installCamera();
  }

  load(spec) {
    this.spec = spec;
    this.frames = [...spec.frames].sort((a, b) => a.t - b.t);
    this.progress = 0;
    this.resetView();
    if (this.equationEl) this.equationEl.textContent = spec.equation_names;
    this._fitViewBox();
    this.render();
    this.onLoad?.(this.frames);
  }

  _fitViewBox() {
    let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
    for (const f of this.frames) {
      for (const a of f.atoms) {
        const r = radiusFor(a.el) + 8;
        minx = Math.min(minx, a.x - r); miny = Math.min(miny, a.y - r);
        maxx = Math.max(maxx, a.x + r); maxy = Math.max(maxy, a.y + r);
      }
    }
    const pad = 48;
    this._base = [minx - pad, miny - pad, maxx - minx + 2 * pad, maxy - miny + 2 * pad];
    this._applyView();
  }

  _applyView() {
    if (!this._base) return;
    const [x, y, w, h] = this._base;
    const s = this.view.scale;
    const vw = w / s, vh = h / s;
    // center of base, shifted by pan (pan is in base units)
    const cx = x + w / 2 + this.view.tx;
    const cy = y + h / 2 + this.view.ty;
    this.svg.setAttribute("viewBox", `${cx - vw / 2} ${cy - vh / 2} ${vw} ${vh}`);
  }

  // ---- scene resolution ----
  _scene() {
    const frames = this.frames;
    if (frames.length === 1) return this._frameScene(frames[0], frames[0], 0);
    const p = this.progress;
    let i = 0;
    while (i < frames.length - 2 && p > frames[i + 1].t) i++;
    const A = frames[i], B = frames[i + 1];
    const span = Math.max(B.t - A.t, 1e-6);
    const u = smoothstep(Math.min(Math.max((p - A.t) / span, 0), 1));
    return this._frameScene(A, B, u);
  }

  _frameScene(A, B, u) {
    const aAtoms = new Map(A.atoms.map((a) => [a.id, a]));
    const bAtoms = new Map(B.atoms.map((a) => [a.id, a]));
    const atoms = new Map();
    for (const id of new Set([...aAtoms.keys(), ...bAtoms.keys()])) {
      const a = aAtoms.get(id), b = bAtoms.get(id);
      if (a && b) atoms.set(id, { x: lerp(a.x, b.x, u), y: lerp(a.y, b.y, u), el: b.el, charge: b.charge, role: b.role || a.role, opacity: 1 });
      else if (a) atoms.set(id, { x: a.x, y: a.y, el: a.el, charge: a.charge, role: a.role, opacity: 1 - u });
      else atoms.set(id, { x: b.x, y: b.y, el: b.el, charge: b.charge, role: b.role, opacity: u });
    }

    const aBonds = new Map(A.bonds.map((bd) => [bondKey(bd.a, bd.b), bd]));
    const bBonds = new Map(B.bonds.map((bd) => [bondKey(bd.a, bd.b), bd]));
    const bonds = [];
    for (const key of new Set([...aBonds.keys(), ...bBonds.keys()])) {
      const a = aBonds.get(key), b = bBonds.get(key);
      const ref = b || a;
      let opacity, state;
      if (a && b) { opacity = 1; state = b.state; }
      else if (a) { opacity = 1 - u; state = a.state === "normal" ? "breaking" : a.state; }
      else { opacity = u; state = b.state === "normal" ? "forming" : b.state; }
      bonds.push({ a: ref.a, b: ref.b, order: ref.order, state, opacity });
    }

    const arrows = [];
    for (const ar of A.arrows || []) arrows.push({ ...ar, opacity: 1 - u });
    for (const ar of B.arrows || []) arrows.push({ ...ar, opacity: u });
    return { atoms, bonds, arrows };
  }

  render() {
    if (!this.spec) return;
    const { atoms, bonds, arrows } = this._scene();
    while (this.svg.firstChild) this.svg.removeChild(this.svg.firstChild);

    const usedEls = new Set([...atoms.values()].map((a) => a.el));
    this.svg.appendChild(this._defs(usedEls));

    // Per-atom bonding tally (for lone-pair counts) and neighbour list (for
    // lone-pair placement). Only count bonds that are substantially present.
    const adj = new Map();
    for (const id of atoms.keys()) adj.set(id, { orderSum: 0, neighbors: [] });
    for (const bd of bonds) {
      if (bd.opacity < 0.35) continue;
      const o = Math.max(1, Math.round(bd.order));
      const ea = adj.get(bd.a), eb = adj.get(bd.b);
      if (ea && atoms.get(bd.b)) { ea.orderSum += o; ea.neighbors.push(atoms.get(bd.b)); }
      if (eb && atoms.get(bd.a)) { eb.orderSum += o; eb.neighbors.push(atoms.get(bd.a)); }
    }

    const root = el("g");
    this.svg.appendChild(root);
    const gBonds = el("g"), gArrows = el("g"), gAtoms = el("g"), gElectrons = el("g");
    root.append(gBonds, gArrows, gAtoms, gElectrons);

    for (const bd of bonds) this._drawBond(gBonds, bd, atoms, gElectrons);
    for (const ar of arrows) this._drawArrow(gArrows, ar, atoms, gElectrons);
    for (const [id, a] of atoms) this._drawAtom(gAtoms, id, a, adj.get(id), gElectrons);

    let active = this.frames[0], idx = 0;
    this.frames.forEach((f, i) => { if (this.progress >= f.t - 1e-9) { active = f; idx = i; } });
    if (this.captionEl) this.captionEl.textContent = active.caption;
    this.onPhase?.(idx, active);
  }

  _defs(usedEls) {
    const defs = el("defs");
    for (const e of usedEls) {
      const base = colorFor(e);
      const grad = el("radialGradient", { id: `grad-${e}`, cx: "35%", cy: "30%", r: "75%" });
      grad.appendChild(el("stop", { offset: "0%", "stop-color": shade(base, 0.55) }));
      grad.appendChild(el("stop", { offset: "55%", "stop-color": base }));
      grad.appendChild(el("stop", { offset: "100%", "stop-color": shade(base, -0.32) }));
      defs.appendChild(grad);
    }
    const glow = el("filter", { id: "glow", x: "-60%", y: "-60%", width: "220%", height: "220%" });
    glow.appendChild(el("feGaussianBlur", { stdDeviation: "3.2", result: "b" }));
    const merge = el("feMerge");
    merge.appendChild(el("feMergeNode", { in: "b" }));
    merge.appendChild(el("feMergeNode", { in: "SourceGraphic" }));
    glow.appendChild(merge);
    defs.appendChild(glow);

    // Softer, hotter glow specifically for electron dots.
    const eglow = el("filter", { id: "eglow", x: "-150%", y: "-150%", width: "400%", height: "400%" });
    eglow.appendChild(el("feGaussianBlur", { stdDeviation: "2.4", result: "b" }));
    const emerge = el("feMerge");
    emerge.appendChild(el("feMergeNode", { in: "b" }));
    emerge.appendChild(el("feMergeNode", { in: "SourceGraphic" }));
    eglow.appendChild(emerge);
    defs.appendChild(eglow);

    const marker = el("marker", { id: "arrowhead", markerWidth: 9, markerHeight: 9, refX: 6, refY: 3, orient: "auto", markerUnits: "strokeWidth" });
    marker.appendChild(el("path", { d: "M0,0 L6,3 L0,6 Z", fill: "#ffd24d" }));
    defs.appendChild(marker);
    return defs;
  }

  // A small glowing electron dot at (x, y).
  _electron(g, x, y, opacity = 1, r = 2.6) {
    g.appendChild(el("circle", {
      cx: x, cy: y, r, fill: ELECTRON_COLOR, opacity,
      filter: "url(#eglow)",
    }));
  }

  _drawBond(g, bd, atoms, gElectrons) {
    const A = atoms.get(bd.a), B = atoms.get(bd.b);
    if (!A || !B) return;
    const op = bd.opacity * Math.min(A.opacity, B.opacity);
    const dx = B.x - A.x, dy = B.y - A.y;
    const len = Math.hypot(dx, dy) || 1;
    const px = -dy / len, py = dx / len;
    const ux = dx / len, uy = dy / len;
    const order = Math.round(bd.order);

    let stroke = "#b9c0c9", dash = "none", width = 5.5, dashoff = 0;
    if (bd.state === "forming") { stroke = "#3ce06a"; dash = "7 7"; dashoff = -this._clock * 22; }
    else if (bd.state === "breaking") { stroke = "#ff6363"; dash = "7 7"; dashoff = this._clock * 22; }

    const offsets = order >= 3 ? [-5, 0, 5] : order === 2 ? [-3.4, 3.4] : [0];
    for (const off of offsets) {
      g.appendChild(el("line", {
        x1: A.x + px * off, y1: A.y + py * off, x2: B.x + px * off, y2: B.y + py * off,
        stroke, "stroke-width": width, "stroke-linecap": "round",
        "stroke-dasharray": dash, "stroke-dashoffset": dashoff, opacity: op,
      }));
    }

    // Electrons streaming along a changing bond: into the bond as it forms,
    // out toward the leaving atom as it breaks. Margin keeps dots off the
    // atom spheres.
    if (bd.state === "forming" || bd.state === "breaking") {
      const m = 0.22;
      const flow = (this._clock * 0.55) % 1; // 0..1 loop
      const s = bd.state === "forming" ? 1 - flow : flow; // forming: B->A (inward), breaking: A->B (outward)
      const t = m + s * (1 - 2 * m);
      const x = A.x + ux * len * t, y = A.y + uy * len * t;
      this._electron(gElectrons, x, y, op * 0.95);
      const t2 = m + ((s + 0.5) % 1) * (1 - 2 * m);
      this._electron(gElectrons, A.x + ux * len * t2, A.y + uy * len * t2, op * 0.5, 2.1);
    }
  }

  _drawAtom(g, id, a, bonding, gElectrons) {
    if (a.opacity <= 0.01) return;
    const r = radiusFor(a.el);
    const grp = el("g", { opacity: a.opacity, cursor: "pointer" });

    if (a.role && a.role !== "normal" && ROLE_COLORS[a.role]) {
      const pulse = 0.5 + 0.5 * Math.sin(this._clock * 3.2);
      grp.appendChild(el("circle", {
        cx: a.x, cy: a.y, r: r + 6 + pulse * 3, fill: "none",
        stroke: ROLE_COLORS[a.role], "stroke-width": 3,
        opacity: 0.55 + 0.4 * pulse, filter: "url(#glow)",
      }));
    }
    grp.appendChild(el("circle", { cx: a.x + 1.4, cy: a.y + 2.2, r, fill: "#000", opacity: 0.28 }));
    grp.appendChild(el("circle", {
      cx: a.x, cy: a.y, r, fill: `url(#grad-${a.el})`,
      stroke: "#0c0f16", "stroke-width": 1.2,
    }));
    const label = el("text", {
      x: a.x, y: a.y, fill: a.el === "H" || a.el === "C" ? "#11151c" : "#ffffff",
      "font-size": r * 1.1, "font-weight": 700, "text-anchor": "middle",
      "dominant-baseline": "central", "font-family": "system-ui, sans-serif",
      "pointer-events": "none",
    });
    label.textContent = a.el;
    grp.appendChild(label);
    if (a.charge) {
      const sign = a.charge > 0 ? "+" : "−";
      const mag = Math.abs(a.charge) === 1 ? "" : Math.abs(a.charge);
      const chg = el("text", {
        x: a.x + r * 0.92, y: a.y - r * 0.82, fill: "#ffe066",
        "font-size": r * 0.95, "font-weight": 700, "text-anchor": "middle",
        "font-family": "system-ui, sans-serif", "pointer-events": "none",
      });
      chg.textContent = mag + sign;
      grp.appendChild(chg);
    }
    if (this.tooltipEl) {
      grp.addEventListener("pointerenter", (e) => this._showTip(e, a));
      grp.addEventListener("pointermove", (e) => this._showTip(e, a));
      grp.addEventListener("pointerleave", () => this._hideTip());
    }
    g.appendChild(grp);

    // Lone pairs — the electrons sitting on the atom. Drawn into the shared
    // electron layer (above bonds) so the cyan reads as one substance whether
    // it is resting on an atom or flowing along an arrow.
    const pairs = lonePairCount(a.el, a.charge, bonding ? bonding.orderSum : 0);
    if (pairs > 0 && a.opacity > 0.2) {
      const angles = lonePairAngles(a, bonding ? bonding.neighbors : [], pairs);
      const dist = r + 8;
      // Anions/nucleophiles breathe slightly to draw the eye to the reactive pair.
      const live = a.charge < 0 || a.role === "nucleophile" || a.role === "base_site";
      const breath = live ? 0.85 + 0.15 * Math.sin(this._clock * 3.0) : 0.8;
      for (const ang of angles) {
        const bx = a.x + Math.cos(ang) * dist, by = a.y + Math.sin(ang) * dist;
        const ppx = -Math.sin(ang), ppy = Math.cos(ang);
        for (const o of [-3.1, 3.1]) {
          this._electron(gElectrons, bx + ppx * o, by + ppy * o, a.opacity * breath, 2.3);
        }
      }
    }
  }

  _drawArrow(g, ar, atoms, gElectrons) {
    const A = atoms.get(ar.src), B = atoms.get(ar.dst);
    if (!A || !B || ar.opacity <= 0.02) return;
    const mx = (A.x + B.x) / 2, my = (A.y + B.y) / 2;
    const dx = B.x - A.x, dy = B.y - A.y;
    const len = Math.hypot(dx, dy) || 1;
    const px = -dy / len, py = dx / len;
    const cx = mx + px * 26, cy = my + py * 26;
    g.appendChild(el("path", {
      d: `M ${A.x} ${A.y} Q ${cx} ${cy} ${B.x} ${B.y}`, fill: "none",
      stroke: "#ffd24d", "stroke-width": 3, opacity: ar.opacity,
      "stroke-linecap": "round", "marker-end": "url(#arrowhead)", filter: "url(#glow)",
    }));

    // An electron pair flowing along the arrow, tail -> head, so the curly
    // arrow shows electrons actually moving rather than just a path.
    const bez = (s) => ({
      x: (1 - s) * (1 - s) * A.x + 2 * (1 - s) * s * cx + s * s * B.x,
      y: (1 - s) * (1 - s) * A.y + 2 * (1 - s) * s * cy + s * s * B.y,
    });
    const head = (this._clock * 0.6) % 1;
    for (const lag of [0, 0.07]) {
      const s = head - lag;
      if (s < 0 || s > 1) continue;
      const p = bez(s);
      this._electron(gElectrons, p.x, p.y, ar.opacity, 2.5);
    }
  }

  // ---- atom inspection (read-only) ----
  _showTip(evt, a) {
    const t = this.tooltipEl;
    const role = a.role && a.role !== "normal" ? ` · ${ROLE_LABELS[a.role] || a.role}` : "";
    const charge = a.charge ? ` (${a.charge > 0 ? "+" : ""}${a.charge})` : "";
    t.innerHTML = `<b>${ELEMENT_NAMES[a.el] || a.el}</b> ${a.el}${charge}${role}`;
    t.style.display = "block";
    const rect = this.svg.getBoundingClientRect();
    t.style.left = `${evt.clientX - rect.left + 14}px`;
    t.style.top = `${evt.clientY - rect.top + 12}px`;
  }
  _hideTip() { if (this.tooltipEl) this.tooltipEl.style.display = "none"; }

  // ---- camera (pan / zoom only; never edits chemistry) ----
  _installCamera() {
    let dragging = false, lx = 0, ly = 0;
    this.svg.addEventListener("wheel", (e) => {
      e.preventDefault();
      const factor = e.deltaY < 0 ? 1.12 : 1 / 1.12;
      this.zoomBy(factor);
    }, { passive: false });
    this.svg.addEventListener("pointerdown", (e) => {
      if (e.target.closest("g[cursor]")) return; // let atoms get their hover
      dragging = true; lx = e.clientX; ly = e.clientY;
      this.svg.setPointerCapture(e.pointerId); this.svg.style.cursor = "grabbing";
    });
    this.svg.addEventListener("pointermove", (e) => {
      if (!dragging || !this._base) return;
      const scale = (this._base[2] / this.svg.clientWidth) / this.view.scale;
      this.view.tx -= (e.clientX - lx) * scale;
      this.view.ty -= (e.clientY - ly) * scale;
      lx = e.clientX; ly = e.clientY;
      this._clampPan(); this._applyView();
    });
    const end = () => { dragging = false; this.svg.style.cursor = "grab"; };
    this.svg.addEventListener("pointerup", end);
    this.svg.addEventListener("pointercancel", end);
    this.svg.style.cursor = "grab";
  }
  zoomBy(factor) {
    this.view.scale = Math.min(Math.max(this.view.scale * factor, 0.6), 4);
    this._clampPan(); this._applyView();
  }
  _clampPan() {
    if (!this._base) return;
    const [, , w, h] = this._base;
    const maxX = w * 0.5, maxY = h * 0.5;
    this.view.tx = Math.min(Math.max(this.view.tx, -maxX), maxX);
    this.view.ty = Math.min(Math.max(this.view.ty, -maxY), maxY);
  }
  resetView() { this.view = { scale: 1, tx: 0, ty: 0 }; this._applyView(); }

  // ---- playback ----
  play() {
    if (this.playing) return;
    if (this.progress >= 1) this.progress = 0;
    this.playing = true;
    this._last = performance.now();
    const step = (now) => {
      if (!this.playing) return;
      const dt = (now - this._last) / 1000;
      this._last = now;
      this._clock += dt;
      const dur = (this.spec.duration_ms || 9000) / 1000;
      this.progress += (dt * this.speed) / dur;
      if (this.progress >= 1) {
        this.progress = this.loop ? this.progress - 1 : 1;
        if (!this.loop) { this.playing = false; this.render(); this.onProgress?.(1); this.onEnd?.(); return; }
      }
      this.render();
      this.onProgress?.(this.progress);
      this._raf = requestAnimationFrame(step);
    };
    this._raf = requestAnimationFrame(step);
  }
  pause() { this.playing = false; if (this._raf) cancelAnimationFrame(this._raf); }
  toggle() { this.playing ? this.pause() : this.play(); }
  restart() { this.progress = 0; this.render(); this.play(); }
  seek(p) { this.pause(); this.progress = Math.min(Math.max(p, 0), 1); this._clock += 0.016; this.render(); this.onProgress?.(this.progress); }
  setSpeed(s) { this.speed = s; }
  setLoop(on) { this.loop = on; }

  stepPhase(dir) {
    const ts = this.frames.map((f) => f.t);
    const eps = 1e-4;
    let target;
    if (dir > 0) target = ts.find((t) => t > this.progress + eps);
    else target = [...ts].reverse().find((t) => t < this.progress - eps);
    if (target === undefined) target = dir > 0 ? 1 : 0;
    this.seek(target);
  }
}
