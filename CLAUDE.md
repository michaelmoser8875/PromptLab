# PromptLab

> Natural-language molecular reaction simulator — describe a reaction in plain English and watch it animate.

## Overview

PromptLab is an agentic chemistry playground that takes a natural-language description of a chemical reaction, resolves the participating molecules, classifies the reaction type, and produces an animated visualization of the mechanism. Reaction types are implemented as parameterized templates so that any compatible set of molecules can be swapped into a known mechanism.

## Goals

- Accept plain-English input (e.g., "what happens when methanol reacts with acetic acid?")
- Resolve named compounds to canonical molecular representations (SMILES)
- Construct balanced reaction equations
- Identify the reaction class (SN2, esterification, acid-base, etc.)
- Animate the mechanism in 2D, then 3D
- Be honest about scope: this is mechanism-templated animation, not first-principles MD

## Non-Goals

- Quantum-chemical accuracy
- Predicting unknown reaction outcomes
- Replacing real molecular dynamics packages
- Supporting every possible reaction — a curated set of common templates is enough

## Tech Stack (proposed)

- **Language:** Python for backend logic, JavaScript/React for the frontend
- **Cheminformatics:** [RDKit](https://www.rdkit.org/) for SMILES parsing, functional-group detection, 2D/3D coordinate generation
- **Name resolution:** [PubChem REST API](https://pubchem.ncbi.nlm.nih.gov/rest/pug/) for compound name → SMILES lookup
- **LLM:** Claude (via Anthropic API) for natural-language parsing and reaction classification
- **2D rendering:** SVG + a JS animation library (Framer Motion, GSAP, or hand-rolled with `requestAnimationFrame`)
- **3D rendering:** [3Dmol.js](https://3dmol.csb.pitt.edu/) or [NGL Viewer](https://nglviewer.org/)

## Project Phases

### Phase 1 — Natural Language → SMILES Identifier

Build the front-of-pipeline component that turns a user's prompt into structured molecular data.

**Deliverables**
- A function/agent that accepts a natural-language string
- Extracts named compounds from the prompt
- Resolves each name to a canonical SMILES string
- Returns a structured object: `{ reactants: [...], products: [...] (if mentioned), context: "..." }`

**Approach**
- Use an LLM to extract compound names from free text
- Resolve names via PubChem (`/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON`)
- Fall back to LLM-generated SMILES for compounds PubChem doesn't return
- Validate every SMILES with RDKit before passing downstream

**Done when**
- 20+ test prompts covering common molecules all resolve correctly
- Invalid/ambiguous inputs return a helpful error rather than nonsense

---

### Phase 2 — Reaction Equation Builder

Take the resolved reactants and produce a balanced chemical equation.

**Deliverables**
- A function that accepts a list of reactant SMILES and (optionally) hints about reaction context
- Returns a balanced equation: reactants → products with stoichiometric coefficients
- Equation is represented as structured data (not just a string) so later phases can consume it

**Approach**
- Use the LLM to propose products given the reactants
- Use RDKit to verify atom balance (count atoms on both sides; adjust coefficients)
- For simple cases (combustion, acid-base, precipitation), prefer deterministic balancing
- For mechanism-driven cases (SN2, esterification), trust the templated product structure

**Done when**
- Every reaction the system supports produces an atom-balanced equation
- Equations render cleanly as both structured data and human-readable strings

---

### Phase 3 — Reaction Type Identifier

Classify the reaction into one of the supported templates.

**Deliverables**
- A classifier that maps `(reactants, products)` → reaction type label
- Initial supported set: SN2, acid-base, esterification, combustion, precipitation (adjust as needed)
- Returns `null` or "unsupported" cleanly when nothing matches

**Approach**
- Two-pass classification: LLM proposes a label, then a rule-based check (using RDKit functional-group SMARTS) verifies the structural fingerprint matches the template
- Rule-based verification prevents the LLM from hallucinating reaction types
- Each supported reaction type has an explicit set of structural requirements documented in code

**Done when**
- 90%+ accuracy on a hand-built test set of 30–50 reaction prompts
- Unsupported reactions fail gracefully with a useful message

---

### Phase 4 — 2D Mechanism Animations (Swappable Templates)

Build the visual core: per-reaction-type animation templates that accept arbitrary compatible molecules.

**Deliverables**
- One animation template per supported reaction type
- Each template is parameterized: it takes the resolved molecules and renders them through the canonical mechanism (approach → transition state → products)
- Atoms and bonds are drawn from the actual SMILES, not hardcoded per example

**Approach**
- Render 2D structures from SMILES using RDKit's coordinate generation
- Export atom positions to SVG and animate via keyframes
- Each template specifies: which atoms are the "actors," which bonds break/form, and a canonical motion path
- Suggested first template: **SN2** — backside attack and Walden inversion are visually striking and instructive

**Done when**
- All supported reaction types have working 2D animations
- Swapping reactants (e.g., a different alkyl halide in SN2) produces a correct animation without code changes

---

### Phase 5 — 3D Animations

Upgrade the visualization to 3D ball-and-stick rendering.

**Deliverables**
- 3D animated version of each Phase 4 template
- Camera framing that highlights the mechanistically important atoms
- Ability to pause, rewind, and rotate the scene

**Approach**
- Use RDKit to generate 3D conformers (ETKDG)
- Render with 3Dmol.js or NGL Viewer
- Animate by interpolating between keyframe conformers (reactant geometry → transition-state geometry → product geometry)
- Reuse the Phase 4 template metadata (which atoms move, which bonds change) — the choreography is shared between 2D and 3D

**Done when**
- All supported reaction types have a 3D animation
- The user can toggle between 2D and 3D views of the same reaction

---

## Stretch Ideas (Post-Phase-5)

- Voice input
- Side-by-side comparison of two reactions
- "Quiz mode": show an animation and ask the user to name the reaction type
- Energy diagram overlay
- Export animations as GIF/MP4
- Expand reaction template library

## Repo Layout (suggested)

```
promptlab/
├── README.md
├── CLAUDE.md
├── backend/
│   ├── nl_parser.py        # Phase 1
│   ├── equation_builder.py # Phase 2
│   ├── reaction_classifier.py       # Phase 3
│   └── templates/          # Phase 4 & 5 mechanism specs
├── frontend/
│   ├── src/
│   │   ├── animations/2d/  # Phase 4
│   │   └── animations/3d/  # Phase 5
└── tests/
    └── reaction_prompts.json
```

## Working Notes for Claude

- Build one phase at a time and verify it independently before moving on
- Every phase should have a small, hand-curated test set checked into the repo
- Prefer deterministic logic (RDKit, rules) over LLM calls wherever both work — LLM is for interpretation, not arithmetic
- When in doubt about a reaction's mechanism, ask the user rather than guessing