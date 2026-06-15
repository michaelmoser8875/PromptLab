"""Shared building blocks for Phase 4 mechanism templates.

Two things live here:

- Layout helpers (arrange a set of molecules into a tidy row, parse the
  namespaced atom ids the renderer interpolates on, tag actor atoms via SMARTS).
- `generic_template`: the default approach -> converge -> cross-fade
  choreography. Reactant molecules drift together, then the scene cross-fades
  into a freshly laid-out product scene. It is correct for any balanced
  equation and needs no per-example tuning, which is what makes templates
  swappable: change the molecules, not the code.

A reaction type with a more instructive mechanism (e.g. SN2's backside attack)
ships its own template instead and only borrows the helpers here.
"""
from __future__ import annotations

import re
from typing import Callable, Optional

from rdkit import Chem

from backend.animation import (
    AnimationSpec,
    Atom,
    Bond,
    Frame,
    MolLayout,
    layout_molecule,
    recenter,
    translate,
)
from backend.equation_builder import BalancedEquation, Species

# role_fn(species, which) -> {(species_index, local_atom_id): role}
RoleFn = Callable[[list[Species], str], dict]

_ID_RE = re.compile(r"^(?P<prefix>[a-z]+)(?P<sp>\d+)_(?P<local>\d+)$")


def match_atom_indices(smiles: str, smarts: str) -> set[str]:
    """Local atom ids (AddHs indexing, as strings) in `smiles` matching `smarts`.

    Indexing matches layout_molecule (which adds hydrogens), so the returned
    ids line up with the atoms the renderer draws.
    """
    mol = Chem.MolFromSmiles(smiles)
    pat = Chem.MolFromSmarts(smarts)
    if mol is None or pat is None:
        return set()
    mol = Chem.AddHs(mol)
    out: set[str] = set()
    for match in mol.GetSubstructMatches(pat):
        for idx in match:
            out.add(str(idx))
    return out


def arrange_row(
    layouts: list[MolLayout], gap: float, prefix: str
) -> tuple[list[Atom], list[Bond]]:
    """Place molecules left-to-right with `gap` between bounding boxes.

    Atom ids become ``f"{prefix}{species_index}_{local_id}"`` so the renderer
    can tell which scene (reactant vs product) an atom belongs to and which
    molecule within it.
    """
    atoms: list[Atom] = []
    bonds: list[Bond] = []
    x = 0.0
    for i, layout in enumerate(layouts):
        centered = recenter(layout)
        minx, _, maxx, _ = centered.bounds()
        width = maxx - minx
        dx = x - minx
        placed = translate(centered, dx, 0.0, prefix=f"{prefix}{i}_")
        atoms.extend(placed.atoms)
        bonds.extend(placed.bonds)
        x += width + gap

    # Recenter the whole row on the origin.
    if atoms:
        cx = (min(a.x for a in atoms) + max(a.x for a in atoms)) / 2
        cy = (min(a.y for a in atoms) + max(a.y for a in atoms)) / 2
        atoms = [a.moved(-cx, -cy) for a in atoms]
    return atoms, bonds


def apply_roles(atoms: list[Atom], role_map: dict) -> None:
    """Mutate atoms in place, setting roles from a {(species, local): role} map."""
    for a in atoms:
        m = _ID_RE.match(a.id)
        if not m:
            continue
        key = (int(m.group("sp")), m.group("local"))
        if key in role_map:
            a.role = role_map[key]


def _layouts_for(species: list[Species]) -> list[MolLayout]:
    return [layout_molecule(s.smiles) for s in species]


def generic_template(
    eq: BalancedEquation,
    reaction_type: str,
    caption_reactants: str,
    caption_converge: str,
    caption_products: str,
    role_fn: Optional[RoleFn] = None,
) -> AnimationSpec:
    """Approach -> converge -> cross-fade choreography for a balanced equation."""
    react_layouts = _layouts_for(eq.reactants)
    prod_layouts = _layouts_for(eq.products)

    # Frame 1: reactants spread apart.
    r_atoms_spread, r_bonds_spread = arrange_row(react_layouts, gap=110.0, prefix="r")
    # Frame 2: same molecules (same ids) drawn closer -> they slide together.
    r_atoms_near, r_bonds_near = arrange_row(react_layouts, gap=34.0, prefix="r")
    # Frame 3: products, freshly laid out -> cross-fade from the reactant scene.
    p_atoms, p_bonds = arrange_row(prod_layouts, gap=70.0, prefix="p")

    if role_fn is not None:
        r_roles = role_fn(eq.reactants, "reactants")
        apply_roles(r_atoms_spread, r_roles)
        apply_roles(r_atoms_near, r_roles)
        p_roles = role_fn(eq.products, "products")
        apply_roles(p_atoms, p_roles)

    frames = [
        Frame(
            id="reactants",
            t=0.0,
            caption=caption_reactants,
            atoms=r_atoms_spread,
            bonds=r_bonds_spread,
        ),
        Frame(
            id="encounter",
            t=0.45,
            caption=caption_converge,
            atoms=r_atoms_near,
            bonds=r_bonds_near,
        ),
        Frame(
            id="products",
            t=1.0,
            caption=caption_products,
            atoms=p_atoms,
            bonds=p_bonds,
        ),
    ]
    return AnimationSpec(
        reaction_type=reaction_type,
        equation_smiles=eq.to_string(),
        equation_names=eq.to_string(use_names=True),
        frames=frames,
    )
