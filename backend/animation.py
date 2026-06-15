"""Phase 4 - 2D mechanism animation builder.

Turns a balanced equation + reaction-type label into a *swappable* animation
spec: a JSON document of keyframes that a renderer interpolates to animate the
mechanism (approach -> transition state -> products).

Design
------
- Geometry is deterministic. RDKit generates 2D coordinates for each molecule
  (no LLM, no hardcoded per-example positions). Swapping one alkyl halide for
  another in an SN2 reaction therefore changes the picture but not the code.
- The choreography is described per reaction *type* by a template (see
  backend/templates/). Each template emits keyframes; this module owns the
  shared data model and the RDKit layout primitives the templates build on.
- The same template metadata - which atoms are actors, which bonds break/form -
  is intended to drive the Phase 5 3D animations too, so it lives in structured
  fields rather than baked into pixels.

A *frame* is a full snapshot of the scene: a list of atoms (each with a stable
id, element, position, role) and bonds (each with an order and a state:
normal / forming / breaking). The renderer transitions between consecutive
frames: an atom id present in both frames is interpolated; an id in only one
fades in or out. This single rule covers both "tracked" mechanisms (SN2, where
actor atoms keep their ids and slide through the inversion) and "cross-fade"
mechanisms (reactant scene morphing into a distinct product scene).
"""
from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass, field
from typing import Optional

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem

from backend.equation_builder import BalancedEquation, build_equation
from backend.nl_parser import parse_reaction

RDLogger.DisableLog("rdApp.*")

# RDKit lays bonds out at ~1.5 units; scale into a comfortable SVG pixel space.
SCALE = 42.0

# Roles drive highlighting in the renderer (and, later, 3D camera framing).
ROLE_NORMAL = "normal"
ROLE_NUCLEOPHILE = "nucleophile"
ROLE_ELECTROPHILE = "electrophile"
ROLE_LEAVING_GROUP = "leaving_group"
ROLE_ACID_PROTON = "acid_proton"
ROLE_BASE_SITE = "base_site"

# Bond states.
BOND_NORMAL = "normal"
BOND_FORMING = "forming"
BOND_BREAKING = "breaking"


@dataclass
class Atom:
    id: str
    el: str
    x: float
    y: float
    charge: int = 0
    role: str = ROLE_NORMAL

    def moved(self, dx: float, dy: float, new_id: Optional[str] = None) -> "Atom":
        return Atom(
            id=new_id if new_id is not None else self.id,
            el=self.el,
            x=self.x + dx,
            y=self.y + dy,
            charge=self.charge,
            role=self.role,
        )


@dataclass
class Bond:
    a: str  # atom id
    b: str  # atom id
    order: float = 1.0
    state: str = BOND_NORMAL


@dataclass
class Arrow:
    """A curved electron-pushing arrow, drawn from one atom toward another."""

    src: str  # atom id (tail)
    dst: str  # atom id (head)
    kind: str = "electron_pair"  # or "single_electron"


@dataclass
class Frame:
    id: str
    t: float  # normalized time in [0, 1]
    caption: str
    atoms: list[Atom] = field(default_factory=list)
    bonds: list[Bond] = field(default_factory=list)
    arrows: list[Arrow] = field(default_factory=list)


@dataclass
class AnimationSpec:
    reaction_type: str
    equation_smiles: str
    equation_names: str
    frames: list[Frame]
    duration_ms: int = 9000

    def to_dict(self) -> dict:
        return {
            "reaction_type": self.reaction_type,
            "equation_smiles": self.equation_smiles,
            "equation_names": self.equation_names,
            "duration_ms": self.duration_ms,
            "frames": [
                {
                    "id": f.id,
                    "t": f.t,
                    "caption": f.caption,
                    "atoms": [asdict(a) for a in f.atoms],
                    "bonds": [asdict(b) for b in f.bonds],
                    "arrows": [asdict(ar) for ar in f.arrows],
                }
                for f in self.frames
            ],
        }


# ---------------------------------------------------------------------------
# RDKit 2D layout primitives
# ---------------------------------------------------------------------------


@dataclass
class MolLayout:
    """A single molecule laid out in 2D, in scaled screen coordinates."""

    smiles: str
    atoms: list[Atom]  # ids are local indices as strings: "0", "1", ...
    bonds: list[Bond]  # a/b reference those local ids

    def bounds(self) -> tuple[float, float, float, float]:
        xs = [a.x for a in self.atoms]
        ys = [a.y for a in self.atoms]
        return min(xs), min(ys), max(xs), max(ys)

    def width(self) -> float:
        minx, _, maxx, _ = self.bounds()
        return maxx - minx

    def height(self) -> float:
        _, miny, _, maxy = self.bounds()
        return maxy - miny

    def centroid(self) -> tuple[float, float]:
        n = len(self.atoms)
        return (sum(a.x for a in self.atoms) / n, sum(a.y for a in self.atoms) / n)


def layout_molecule(smiles: str, with_h: bool = True) -> MolLayout:
    """RDKit 2D coordinates for `smiles`, scaled into screen space (y points down).

    Hydrogens are added by default so proton transfers and O-H / N-H actors are
    visible; pass with_h=False for a heavy-atom skeleton.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot lay out invalid SMILES: {smiles!r}")
    if with_h:
        mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()

    atoms: list[Atom] = []
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        atoms.append(
            Atom(
                id=str(a.GetIdx()),
                el=a.GetSymbol(),
                x=p.x * SCALE,
                y=-p.y * SCALE,  # flip to screen coordinates
                charge=a.GetFormalCharge(),
            )
        )
    bonds: list[Bond] = []
    for b in mol.GetBonds():
        bonds.append(
            Bond(
                a=str(b.GetBeginAtomIdx()),
                b=str(b.GetEndAtomIdx()),
                order=b.GetBondTypeAsDouble(),
            )
        )
    return MolLayout(smiles=smiles, atoms=atoms, bonds=bonds)


def translate(layout: MolLayout, dx: float, dy: float, prefix: str) -> MolLayout:
    """Return a copy shifted by (dx, dy) with every atom id namespaced by prefix."""
    atoms = [a.moved(dx, dy, new_id=f"{prefix}{a.id}") for a in layout.atoms]
    bonds = [
        Bond(a=f"{prefix}{b.a}", b=f"{prefix}{b.b}", order=b.order, state=b.state)
        for b in layout.bonds
    ]
    return MolLayout(smiles=layout.smiles, atoms=atoms, bonds=bonds)


def recenter(layout: MolLayout) -> MolLayout:
    """Shift a layout so its centroid sits at the origin (ids unchanged)."""
    cx, cy = layout.centroid()
    atoms = [a.moved(-cx, -cy) for a in layout.atoms]
    return MolLayout(smiles=layout.smiles, atoms=atoms, bonds=list(layout.bonds))


# ---------------------------------------------------------------------------
# Template dispatch
# ---------------------------------------------------------------------------


def build_animation(eq: BalancedEquation, reaction_type: str) -> AnimationSpec:
    """Build an animation spec for a balanced equation of the given type.

    Falls back to the generic approach/converge/cross-fade template for any
    type without a dedicated choreographer, so the call always succeeds for a
    supported, balanced equation.
    """
    # Imported here to avoid a circular import (templates import this module).
    from backend.templates import get_template

    template = get_template(reaction_type)
    return template(eq, reaction_type)


def _main() -> None:
    import os

    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    out_flag = next((a for a in sys.argv[1:] if a.startswith("--out=")), None)
    prompt = " ".join(args).strip()
    if not prompt:
        prompt = "react bromomethane with hydroxide"
        print(f"(no prompt given; using default: {prompt!r})", file=sys.stderr)

    from backend.reaction_classifier import classify

    parsed = parse_reaction(prompt)
    eq = build_equation(parsed)
    classification = classify(eq)
    if classification.reaction_type == "unsupported":
        print(
            f"Cannot animate: reaction classified as unsupported "
            f"({classification.reason})",
            file=sys.stderr,
        )
        sys.exit(1)

    spec = build_animation(eq, classification.reaction_type)
    payload = json.dumps(spec.to_dict(), indent=2)
    if out_flag:
        path = out_flag.split("=", 1)[1]
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(payload)
        print(f"wrote {path}")
    else:
        print(payload)


if __name__ == "__main__":
    _main()
