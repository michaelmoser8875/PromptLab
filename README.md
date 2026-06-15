# PromptLab

> Natural-language molecular reaction simulator — describe a reaction in plain
> English and watch it animate.

PromptLab takes a plain-English description of a chemical reaction, resolves the
molecules, balances the equation, classifies the reaction type, and animates the
mechanism. It is honest about its scope: **mechanism-templated animation, not
first-principles molecular dynamics.**

## Pipeline

| Phase | Module | What it does |
|------|--------|--------------|
| 1 | [`backend/nl_parser.py`](backend/nl_parser.py) | Natural language → named compounds → canonical SMILES (PubChem, LLM fallback, RDKit-validated) |
| 2 | [`backend/equation_builder.py`](backend/equation_builder.py) | Reactants → predicted products → **atom-balanced** equation (exact rational nullspace) |
| 3 | [`backend/reaction_classifier.py`](backend/reaction_classifier.py) | Equation → reaction type (LLM proposes, RDKit SMARTS verifies) |
| 4 | [`backend/animation.py`](backend/animation.py) + [`backend/templates/`](backend/templates/) + [`frontend/`](frontend/) | Equation + type → **2D mechanism animation** |

Supported reaction types: **SN2, acid-base, esterification, combustion,
precipitation.** Anything else classifies as `unsupported` and is not animated.

## Setup

```bash
pip install -r requirements.txt
export ANTHROPIC_API_KEY=...   # needed for phases 1–3 (parsing/classification)
```

Phase 4 geometry is fully deterministic (RDKit only) and needs **no API key**.

## Phase 4 — 2D mechanism animations

The animation builder turns a balanced equation + reaction-type label into an
**AnimationSpec**: a JSON document of keyframes (`approach → transition state →
products`) that the frontend interpolates. Geometry comes from RDKit's 2D
coordinates of the *actual* molecules, so the templates are **swappable** —
change the reactants and the animation re-choreographs itself with no code
change (e.g. bromomethane → 1-chlorobutane, hydroxide → cyanide).

Each frame is a full scene snapshot: atoms (id, element, position, charge,
role) and bonds (order, state: `normal`/`forming`/`breaking`). The renderer
transitions between consecutive frames — an atom id present in both is
interpolated, an id in only one fades. That single rule drives both the
atom-tracked **SN2** inversion (backside attack → planar TS → Walden inversion)
and the cross-fade mechanisms for the other types.

### Generate specs and view the animations

```bash
# 1. Build example specs (offline, no API key) into frontend/specs/
python -m backend.gen_examples

# 2. Serve the frontend and open the printed URL
cd frontend && python -m http.server 8000
#    → http://localhost:8000/

# Or build a spec straight from a prompt (uses the API for phases 1–3):
python -m backend.animation "react bromomethane with hydroxide" \
    --out=frontend/specs/custom.json
```

The viewer has a reaction picker, play/pause/scrub, and speed control.
Deep-link to a frozen frame with `?spec=<file>&t=<0..1>` or auto-play with
`?play`.

Roles are color-coded: nucleophile/base (cyan), electrophile (orange), leaving
group (red), acidic proton (yellow); bonds forming (green dashed) and breaking
(red dashed).

### Tests

```bash
python tests/test_phase4.py     # offline: spec validity + SN2 mechanism + swap-invariance
python tests/test_phase3.py     # live (needs ANTHROPIC_API_KEY): classification accuracy
```

## Repo layout

```
backend/
  nl_parser.py          # Phase 1
  equation_builder.py   # Phase 2
  reaction_classifier.py# Phase 3
  animation.py          # Phase 4 spec model + RDKit 2D layout + dispatch
  templates/            # Phase 4 mechanism choreographers (sn2 + generic engine)
  gen_examples.py       # emit example specs offline
frontend/
  index.html            # 2D animation viewer
  src/animations/2d/player.js
  specs/                # generated AnimationSpec JSON
tests/
  reaction_prompts.json # curated prompt set
  test_phase{1,2,3,4}.py
```

## Roadmap

Phase 5 reuses the Phase 4 template metadata (which atoms move, which bonds
change) to drive 3D ball-and-stick animations.
