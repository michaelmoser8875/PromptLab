// Reaction-coordinate energy diagram — doubles as the timeline scrubber.
//
// The horizontal axis is the reaction coordinate (progress 0..1); the curve is
// a schematic energy profile (reactants -> transition state -> products). The
// draggable marker rides the curve, so the only thing a user can do is move
// along the real mechanism path — they cannot scrub to an off-pathway, unphysical
// state. Keyframe phases are marked and labeled.
//
// Energies are schematic (a teaching cue for "uphill to the transition state,
// downhill to products"), not computed thermochemistry.

const SVG_NS = "http://www.w3.org/2000/svg";
const W = 1000, H = 250;
const PAD_L = 54, PAD_R = 28, PAD_T = 26, PAD_B = 52;
const X0 = PAD_L, X1 = W - PAD_R, Y0 = H - PAD_B, Y1 = PAD_T;

const PHASE_NAMES = {
  reactants: "Reactants", approach: "Approach", transition: "Transition state",
  products: "Products", reaction: "Reaction", encounter: "Encounter",
};
const phaseName = (id) => PHASE_NAMES[id] || id.charAt(0).toUpperCase() + id.slice(1);

function el(tag, attrs = {}, text) {
  const n = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
  if (text != null) n.textContent = text;
  return n;
}

function energyOf(frame, i, n) {
  if (i === 0) return 0.40;
  if (i === n - 1) return 0.16;
  if (/transition|reaction/.test(frame.id)) return 0.96;
  return 0.66;
}

const xOf = (t) => X0 + t * (X1 - X0);
const yOf = (e) => Y0 - e * (Y0 - Y1);

export class EnergyDiagram {
  constructor(svg, { onSeek } = {}) {
    this.svg = svg;
    this.onSeek = onSeek;
    this.frames = [];
    this.points = []; // {t, e}
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    svg.setAttribute("preserveAspectRatio", "none");
    this._installDrag();
  }

  // Schematic energy as a function of progress, smooth through the keyframes.
  _evalE(t) {
    const p = this.points;
    if (p.length === 0) return 0.4;
    if (t <= p[0].t) return p[0].e;
    if (t >= p[p.length - 1].t) return p[p.length - 1].e;
    let i = 0;
    while (i < p.length - 1 && t > p[i + 1].t) i++;
    const p0 = p[Math.max(i - 1, 0)], p1 = p[i], p2 = p[i + 1], p3 = p[Math.min(i + 2, p.length - 1)];
    const local = (t - p1.t) / Math.max(p2.t - p1.t, 1e-6);
    const t2 = local * local, t3 = t2 * local;
    // Catmull-Rom (uniform) on the energy values.
    return 0.5 * ((2 * p1.e) + (-p0.e + p2.e) * local +
      (2 * p0.e - 5 * p1.e + 4 * p2.e - p3.e) * t2 +
      (-p0.e + 3 * p1.e - 3 * p2.e + p3.e) * t3);
  }

  setFrames(frames) {
    this.frames = frames;
    const n = frames.length;
    this.points = frames.map((f, i) => ({ t: f.t, e: energyOf(f, i, n) }));
    this._render(0);
  }

  setProgress(p) { this._render(p); }

  _render(progress) {
    const svg = this.svg;
    while (svg.firstChild) svg.removeChild(svg.firstChild);

    const defs = el("defs");
    const grad = el("linearGradient", { id: "ediag-fill", x1: "0", y1: "0", x2: "0", y2: "1" });
    grad.appendChild(el("stop", { offset: "0%", "stop-color": "#2f6dff", "stop-opacity": "0.35" }));
    grad.appendChild(el("stop", { offset: "100%", "stop-color": "#2f6dff", "stop-opacity": "0.02" }));
    defs.appendChild(grad);
    svg.appendChild(defs);

    // axes
    svg.appendChild(el("line", { x1: X0, y1: Y1 - 6, x2: X0, y2: Y0, stroke: "#2a3344", "stroke-width": 1.5 }));
    svg.appendChild(el("line", { x1: X0, y1: Y0, x2: X1 + 6, y2: Y0, stroke: "#2a3344", "stroke-width": 1.5 }));
    svg.appendChild(el("text", { x: 16, y: (Y0 + Y1) / 2, fill: "#8b94a3", "font-size": 15, "font-family": "system-ui, sans-serif", "text-anchor": "middle", transform: `rotate(-90 16 ${(Y0 + Y1) / 2})` }, "Energy"));
    svg.appendChild(el("text", { x: (X0 + X1) / 2, y: H - 12, fill: "#8b94a3", "font-size": 15, "font-family": "system-ui, sans-serif", "text-anchor": "middle" }, "Reaction coordinate →"));

    if (!this.points.length) return;

    // sample the smooth curve
    const N = 200;
    const pts = [];
    for (let k = 0; k <= N; k++) {
      const t = k / N;
      pts.push([xOf(t), yOf(this._evalE(t))]);
    }
    const dCurve = "M " + pts.map(([x, y]) => `${x.toFixed(1)} ${y.toFixed(1)}`).join(" L ");
    const dArea = dCurve + ` L ${X1} ${Y0} L ${X0} ${Y0} Z`;
    svg.appendChild(el("path", { d: dArea, fill: "url(#ediag-fill)" }));
    svg.appendChild(el("path", { d: dCurve, fill: "none", stroke: "#6ea0ff", "stroke-width": 3, "stroke-linejoin": "round" }));

    // keyframe phase markers + labels
    this.frames.forEach((f, i) => {
      const x = xOf(f.t), y = yOf(this.points[i].e);
      svg.appendChild(el("line", { x1: x, y1: y, x2: x, y2: Y0, stroke: "#2a3344", "stroke-width": 1, "stroke-dasharray": "3 4" }));
      svg.appendChild(el("circle", { cx: x, cy: y, r: 4.5, fill: "#0c0f16", stroke: "#6ea0ff", "stroke-width": 2 }));
      const anchor = i === 0 ? "start" : i === this.frames.length - 1 ? "end" : "middle";
      svg.appendChild(el("text", { x, y: Y0 + 22, fill: "#aeb8c8", "font-size": 14, "font-family": "system-ui, sans-serif", "text-anchor": anchor }, phaseName(f.id)));
    });

    // progress marker
    const mx = xOf(progress), my = yOf(this._evalE(progress));
    svg.appendChild(el("line", { x1: mx, y1: Y1 - 6, x2: mx, y2: Y0, stroke: "#ffd24d", "stroke-width": 1.5, opacity: 0.5 }));
    svg.appendChild(el("circle", { cx: mx, cy: my, r: 9, fill: "#ffd24d", opacity: 0.18 }));
    svg.appendChild(el("circle", { cx: mx, cy: my, r: 6, fill: "#ffd24d", stroke: "#0c0f16", "stroke-width": 2 }));
  }

  _installDrag() {
    const toProgress = (e) => {
      const rect = this.svg.getBoundingClientRect();
      const fx = (e.clientX - rect.left) / rect.width * W; // into viewBox units
      return Math.min(Math.max((fx - X0) / (X1 - X0), 0), 1);
    };
    let dragging = false;
    const move = (e) => { if (dragging) this.onSeek?.(toProgress(e)); };
    this.svg.style.cursor = "ew-resize";
    this.svg.addEventListener("pointerdown", (e) => {
      dragging = true; this.svg.setPointerCapture(e.pointerId); this.onSeek?.(toProgress(e));
    });
    this.svg.addEventListener("pointermove", move);
    const end = () => { dragging = false; };
    this.svg.addEventListener("pointerup", end);
    this.svg.addEventListener("pointercancel", end);
  }
}
