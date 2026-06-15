// Phase 4 - 2D mechanism animation player.
//
// Consumes an AnimationSpec (see backend/animation.py) and renders it as an
// animated SVG. The spec is a list of keyframes; this player interpolates
// between consecutive frames:
//   - an atom id present in both frames slides from one position to the other;
//   - an id in only one frame fades in or out.
// That single rule renders both the atom-tracked SN2 inversion and the
// cross-fade mechanisms (acid-base, esterification, ...) with no per-type code.

const SVG_NS = "http://www.w3.org/2000/svg";

// CPK-ish palette tuned for a dark background.
const ELEMENT_COLORS = {
  H: "#e8e8e8", C: "#9aa0a6", N: "#5577ff", O: "#ff4d4d", F: "#7ee6a8",
  Cl: "#3cc94f", Br: "#c1622f", I: "#a64dff", S: "#ffd633", P: "#ff9f40",
  Na: "#8866ff", K: "#8f4bd6", Ag: "#c0c0c8", Mg: "#3fd6a0", Ca: "#4fd0e0",
};
const ELEMENT_RADII = {
  H: 9, C: 15, N: 14, O: 14, F: 13, Cl: 17, Br: 18, I: 20, S: 16, P: 16,
};
const ROLE_COLORS = {
  nucleophile: "#2fd0ff",
  electrophile: "#ff9f43",
  leaving_group: "#ff5e5e",
  acid_proton: "#ffd633",
  base_site: "#2fd0ff",
};

const colorFor = (el) => ELEMENT_COLORS[el] || "#ff80c0";
const radiusFor = (el) => ELEMENT_RADII[el] || 14;
const lerp = (a, b, u) => a + (b - a) * u;
const bondKey = (a, b) => (a < b ? `${a}|${b}` : `${b}|${a}`);

function el(tag, attrs = {}) {
  const node = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  return node;
}

export class MechanismPlayer {
  constructor(svg, captionEl, equationEl) {
    this.svg = svg;
    this.captionEl = captionEl;
    this.equationEl = equationEl;
    this.spec = null;
    this.frames = [];
    this.progress = 0;
    this.playing = false;
    this.speed = 1;
    this._raf = null;
    this._last = null;
  }

  load(spec) {
    this.spec = spec;
    this.frames = [...spec.frames].sort((a, b) => a.t - b.t);
    this.progress = 0;
    if (this.equationEl) this.equationEl.textContent = spec.equation_names;
    this._fitViewBox();
    this.render();
  }

  _fitViewBox() {
    let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
    for (const f of this.frames) {
      for (const a of f.atoms) {
        const r = radiusFor(a.el);
        minx = Math.min(minx, a.x - r);
        miny = Math.min(miny, a.y - r);
        maxx = Math.max(maxx, a.x + r);
        maxy = Math.max(maxy, a.y + r);
      }
    }
    const pad = 48;
    minx -= pad; miny -= pad; maxx += pad; maxy += pad;
    this.svg.setAttribute("viewBox", `${minx} ${miny} ${maxx - minx} ${maxy - miny}`);
  }

  // Resolve the interpolated scene at the current progress.
  _scene() {
    const frames = this.frames;
    if (frames.length === 1) return this._frameScene(frames[0], frames[0], 0);
    const p = this.progress;
    let i = 0;
    while (i < frames.length - 2 && p > frames[i + 1].t) i++;
    const A = frames[i], B = frames[i + 1];
    const span = Math.max(B.t - A.t, 1e-6);
    const u = Math.min(Math.max((p - A.t) / span, 0), 1);
    return this._frameScene(A, B, u);
  }

  _frameScene(A, B, u) {
    const aAtoms = new Map(A.atoms.map((a) => [a.id, a]));
    const bAtoms = new Map(B.atoms.map((a) => [a.id, a]));
    const atoms = new Map(); // id -> {x,y,el,charge,role,opacity}
    for (const id of new Set([...aAtoms.keys(), ...bAtoms.keys()])) {
      const a = aAtoms.get(id), b = bAtoms.get(id);
      if (a && b) {
        atoms.set(id, {
          x: lerp(a.x, b.x, u), y: lerp(a.y, b.y, u),
          el: b.el, charge: b.charge, role: b.role || a.role, opacity: 1,
        });
      } else if (a) {
        atoms.set(id, { x: a.x, y: a.y, el: a.el, charge: a.charge, role: a.role, opacity: 1 - u });
      } else {
        atoms.set(id, { x: b.x, y: b.y, el: b.el, charge: b.charge, role: b.role, opacity: u });
      }
    }

    const aBonds = new Map(A.bonds.map((bd) => [bondKey(bd.a, bd.b), bd]));
    const bBonds = new Map(B.bonds.map((bd) => [bondKey(bd.a, bd.b), bd]));
    const bonds = [];
    for (const key of new Set([...aBonds.keys(), ...bBonds.keys()])) {
      const a = aBonds.get(key), b = bBonds.get(key);
      const ref = b || a;
      let opacity, state;
      // The destination frame is the source of truth: a bond that finished
      // forming (forming -> normal) should read as a normal bond once reached.
      if (a && b) { opacity = 1; state = b.state; }
      else if (a) { opacity = 1 - u; state = a.state === "normal" ? "breaking" : a.state; }
      else { opacity = u; state = b.state === "normal" ? "forming" : b.state; }
      bonds.push({ a: ref.a, b: ref.b, order: ref.order, state, opacity });
    }

    // Arrows belong to a keyframe; ramp A's down and B's up across the segment.
    const arrows = [];
    for (const ar of A.arrows || []) arrows.push({ ...ar, opacity: 1 - u });
    for (const ar of B.arrows || []) arrows.push({ ...ar, opacity: u });

    return { atoms, bonds, arrows };
  }

  render() {
    const { atoms, bonds, arrows } = this._scene();
    while (this.svg.firstChild) this.svg.removeChild(this.svg.firstChild);

    const gBonds = el("g");
    const gAtoms = el("g");
    const gArrows = el("g");
    this.svg.append(gBonds, gArrows, gAtoms);

    for (const bd of bonds) this._drawBond(gBonds, bd, atoms);
    for (const ar of arrows) this._drawArrow(gArrows, ar, atoms);
    for (const [, a] of atoms) this._drawAtom(gAtoms, a);

    // Caption: the latest keyframe we've reached.
    let active = this.frames[0];
    for (const f of this.frames) if (this.progress >= f.t - 1e-9) active = f;
    if (this.captionEl) this.captionEl.textContent = active.caption;
  }

  _drawBond(g, bd, atoms) {
    const A = atoms.get(bd.a), B = atoms.get(bd.b);
    if (!A || !B) return;
    const op = bd.opacity * Math.min(A.opacity, B.opacity) + 0.0;
    const dx = B.x - A.x, dy = B.y - A.y;
    const len = Math.hypot(dx, dy) || 1;
    const px = -dy / len, py = dx / len; // unit perpendicular
    const order = Math.round(bd.order);

    let stroke = "#b9c0c9", dash = "none", width = 5;
    if (bd.state === "forming") { stroke = "#3ce06a"; dash = "6 7"; }
    else if (bd.state === "breaking") { stroke = "#ff6363"; dash = "6 7"; }

    const offsets = order >= 3 ? [-5, 0, 5] : order === 2 ? [-3.2, 3.2] : [0];
    for (const off of offsets) {
      g.appendChild(el("line", {
        x1: A.x + px * off, y1: A.y + py * off,
        x2: B.x + px * off, y2: B.y + py * off,
        stroke, "stroke-width": width, "stroke-linecap": "round",
        "stroke-dasharray": dash, opacity: op,
      }));
    }
  }

  _drawAtom(g, a) {
    if (a.opacity <= 0.01) return;
    const r = radiusFor(a.el);
    const grp = el("g", { opacity: a.opacity });
    if (a.role && a.role !== "normal" && ROLE_COLORS[a.role]) {
      grp.appendChild(el("circle", {
        cx: a.x, cy: a.y, r: r + 7,
        fill: "none", stroke: ROLE_COLORS[a.role], "stroke-width": 3, opacity: 0.85,
      }));
    }
    grp.appendChild(el("circle", {
      cx: a.x, cy: a.y, r, fill: colorFor(a.el),
      stroke: "#0c0f16", "stroke-width": 1.5,
    }));
    const label = el("text", {
      x: a.x, y: a.y, fill: a.el === "H" || a.el === "C" ? "#11151c" : "#ffffff",
      "font-size": r * 1.15, "font-weight": 700,
      "text-anchor": "middle", "dominant-baseline": "central",
      "font-family": "system-ui, sans-serif",
    });
    label.textContent = a.el;
    grp.appendChild(label);
    if (a.charge) {
      const sign = a.charge > 0 ? "+" : "−";
      const mag = Math.abs(a.charge) === 1 ? "" : Math.abs(a.charge);
      const chg = el("text", {
        x: a.x + r * 0.9, y: a.y - r * 0.8, fill: "#ffe066",
        "font-size": r * 0.95, "font-weight": 700, "text-anchor": "middle",
        "font-family": "system-ui, sans-serif",
      });
      chg.textContent = mag + sign;
      grp.appendChild(chg);
    }
    g.appendChild(grp);
  }

  _drawArrow(g, ar, atoms) {
    const A = atoms.get(ar.src), B = atoms.get(ar.dst);
    if (!A || !B || ar.opacity <= 0.02) return;
    const mx = (A.x + B.x) / 2, my = (A.y + B.y) / 2;
    const dx = B.x - A.x, dy = B.y - A.y;
    const len = Math.hypot(dx, dy) || 1;
    const px = -dy / len, py = dx / len;
    const bow = 26;
    const cx = mx + px * bow, cy = my + py * bow;
    const path = el("path", {
      d: `M ${A.x} ${A.y} Q ${cx} ${cy} ${B.x} ${B.y}`,
      fill: "none", stroke: "#ffd24d", "stroke-width": 3,
      opacity: ar.opacity, "marker-end": "url(#arrowhead)",
    });
    g.appendChild(this._ensureArrowhead());
    g.appendChild(path);
  }

  _ensureArrowhead() {
    const defs = el("defs");
    const marker = el("marker", {
      id: "arrowhead", markerWidth: 8, markerHeight: 8,
      refX: 6, refY: 3, orient: "auto", markerUnits: "strokeWidth",
    });
    const poly = el("path", { d: "M0,0 L6,3 L0,6 Z", fill: "#ffd24d" });
    marker.appendChild(poly);
    defs.appendChild(marker);
    return defs;
  }

  // ---- playback ----
  play() {
    if (this.playing) return;
    this.playing = true;
    this._last = performance.now();
    const step = (now) => {
      if (!this.playing) return;
      const dt = (now - this._last) / 1000;
      this._last = now;
      const dur = (this.spec.duration_ms || 9000) / 1000;
      this.progress += (dt * this.speed) / dur;
      if (this.progress >= 1) { this.progress = 1; this.playing = false; this.render(); this.onEnd?.(); return; }
      this.render();
      this.onProgress?.(this.progress);
      this._raf = requestAnimationFrame(step);
    };
    this._raf = requestAnimationFrame(step);
  }

  pause() { this.playing = false; if (this._raf) cancelAnimationFrame(this._raf); }
  toggle() { this.playing ? this.pause() : (this.progress >= 1 ? this.restart() : this.play()); }
  restart() { this.progress = 0; this.render(); this.play(); }
  seek(p) { this.pause(); this.progress = Math.min(Math.max(p, 0), 1); this.render(); }
  setSpeed(s) { this.speed = s; }
}
