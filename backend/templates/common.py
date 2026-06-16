"""Shared building blocks for mechanism templates.

The default choreographer here is an **atom-mapped morph**: reactant atoms and
product atoms are matched to a shared id space, so the renderer slides each atom
from where it starts to where it ends up while bonds that disappear break and
bonds that appear form. That gives every reaction real on-screen dynamics
(motion + bond changes) rather than a cross-fade, using only the balanced
equation and RDKit geometry - so swapping reactants re-choreographs the
animation with no code change.

Stoichiometric coefficients are honored by drawing each species the right
number of times, which also makes the per-element atom counts on the two sides
match exactly so the mapping is complete.

A reaction with a more instructive, specific mechanism (e.g. SN2's backside
attack) ships its own template and only borrows the layout helpers here.
"""
from __future__ import annotations

import math
from collections import defaultdict
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

# role_fn(smiles, which) -> {local_atom_id: role}, evaluated per molecule.
RoleFn = Callable[[str, str], dict]

# Don't draw an absurd number of copies for huge coefficients.
_MAX_COPIES = 8


def match_atom_indices(smiles: str, smarts: str) -> set[str]:
    """Local atom ids (AddHs indexing, as strings) in `smiles` matching `smarts`."""
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
    """Place molecules left-to-right with `gap` between bounding boxes, centered.

    Atom ids become ``f"{prefix}{slot_index}_{local_id}"``.
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

    if atoms:
        cx = (min(a.x for a in atoms) + max(a.x for a in atoms)) / 2
        cy = (min(a.y for a in atoms) + max(a.y for a in atoms)) / 2
        atoms = [a.moved(-cx, -cy) for a in atoms]
    return atoms, bonds


def _orient_horizontal(layout: MolLayout) -> MolLayout:
    """Rotate a molecule so its longest axis is horizontal (tidier rows, less
    vertical sprawl from molecules RDKit happens to draw tall)."""
    atoms = layout.atoms
    if len(atoms) < 2:
        return layout
    cx = sum(a.x for a in atoms) / len(atoms)
    cy = sum(a.y for a in atoms) / len(atoms)
    sxx = sum((a.x - cx) ** 2 for a in atoms)
    syy = sum((a.y - cy) ** 2 for a in atoms)
    sxy = sum((a.x - cx) * (a.y - cy) for a in atoms)
    theta = 0.5 * math.atan2(2 * sxy, sxx - syy)  # major-axis angle
    cs, sn = math.cos(-theta), math.sin(-theta)
    rotated = []
    for a in atoms:
        dx, dy = a.x - cx, a.y - cy
        rotated.append(Atom(
            id=a.id, el=a.el, charge=a.charge, role=a.role,
            x=cx + dx * cs - dy * sn, y=cy + dx * sn + dy * cs,
        ))
    return MolLayout(smiles=layout.smiles, atoms=rotated, bonds=list(layout.bonds))


def _fragments(smiles: str) -> list[str]:
    """Split a species into its disconnected fragments (e.g. a salt's ions).

    RDKit packs disconnected components poorly, so each ion/molecule is laid
    out on its own. This also makes ionic mechanisms read as ions migrating to
    new partners.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [smiles]
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    return [Chem.MolToSmiles(f) for f in frags] or [smiles]


def _expand(species: list[Species]) -> tuple[list[MolLayout], list[str]]:
    """Lay out each species `coefficient` times, split into fragments."""
    layouts: list[MolLayout] = []
    smiles_of: list[str] = []
    for s in species:
        copies = min(max(s.coefficient, 1), _MAX_COPIES)
        for _ in range(copies):
            for frag in _fragments(s.smiles):
                layouts.append(_orient_horizontal(layout_molecule(frag)))
                smiles_of.append(frag)
    return layouts, smiles_of


def _apply_roles(
    atoms: list[Atom], smiles_of: list[str], role_fn: Optional[RoleFn], which: str
) -> None:
    if role_fn is None:
        return
    cache: dict[str, dict] = {}
    for a in atoms:
        slot_local = a.id.split("_", 1)
        if len(slot_local) != 2:
            continue
        slot = int(slot_local[0].lstrip("rp"))
        local = slot_local[1]
        smi = smiles_of[slot]
        if smi not in cache:
            cache[smi] = role_fn(smi, which)
        role = cache[smi].get(local)
        if role:
            a.role = role


def _atom_map(r_atoms: list[Atom], p_atoms: list[Atom]) -> dict[str, str]:
    """Match reactant atom ids to product atom ids, element by element.

    Greedy nearest-pair within each element bucket - cheap, and since both
    scenes are centered on the origin it keeps motion local and orderly.
    """
    r_by: dict[str, list[Atom]] = defaultdict(list)
    p_by: dict[str, list[Atom]] = defaultdict(list)
    for a in r_atoms:
        r_by[a.el].append(a)
    for a in p_atoms:
        p_by[a.el].append(a)

    mapping: dict[str, str] = {}
    for el, rs in r_by.items():
        ps = p_by.get(el, [])
        candidates = []
        for ra in rs:
            for pa in ps:
                d = (ra.x - pa.x) ** 2 + (ra.y - pa.y) ** 2
                candidates.append((d, ra.id, pa.id))
        candidates.sort()
        used_r: set[str] = set()
        used_p: set[str] = set()
        for _, rid, pid in candidates:
            if rid in used_r or pid in used_p:
                continue
            mapping[rid] = pid
            used_r.add(rid)
            used_p.add(pid)
    return mapping


def generic_template(
    eq: BalancedEquation,
    reaction_type: str,
    caption_reactants: str,
    caption_react: str,
    caption_products: str,
    role_fn: Optional[RoleFn] = None,
) -> AnimationSpec:
    """Atom-mapped morph: atoms slide from reactant to product layout while
    bonds break and form."""
    react_layouts, r_smiles = _expand(eq.reactants)
    prod_layouts, p_smiles = _expand(eq.products)

    r_atoms, r_bonds = arrange_row(react_layouts, gap=78.0, prefix="r")
    p_atoms, p_bonds = arrange_row(prod_layouts, gap=70.0, prefix="p")
    _apply_roles(r_atoms, r_smiles, role_fn, "reactants")
    _apply_roles(p_atoms, p_smiles, role_fn, "products")

    mapping = _atom_map(r_atoms, p_atoms)
    # Shared ids: every matched pair gets one stable id ("s{n}"); the few
    # unmatched atoms (only if counts somehow differ) keep a scene-local id and
    # simply fade.
    shared_of_r: dict[str, str] = {}
    shared_of_p: dict[str, str] = {}
    for n, (rid, pid) in enumerate(mapping.items()):
        sid = f"s{n}"
        shared_of_r[rid] = sid
        shared_of_p[pid] = sid

    r_by_id = {a.id: a for a in r_atoms}
    p_by_id = {a.id: a for a in p_atoms}

    def sid_r(rid: str) -> str:
        return shared_of_r.get(rid, "u_" + rid)

    def sid_p(pid: str) -> str:
        return shared_of_p.get(pid, "u_" + pid)

    # Reactant-frame atoms (start positions, reactant charges/roles).
    react_frame_atoms = [
        Atom(id=sid_r(a.id), el=a.el, x=a.x, y=a.y, charge=a.charge, role=a.role)
        for a in r_atoms
    ]
    # Product-frame atoms (end positions, product charges/roles).
    prod_frame_atoms = [
        Atom(id=sid_p(a.id), el=a.el, x=a.x, y=a.y, charge=a.charge, role=a.role)
        for a in p_atoms
    ]
    # Transition-frame atoms: matched atoms at the midpoint of their travel.
    trans_atoms: list[Atom] = []
    for rid, pid in mapping.items():
        ra, pa = r_by_id[rid], p_by_id[pid]
        trans_atoms.append(Atom(
            id=shared_of_r[rid], el=ra.el,
            x=(ra.x + pa.x) / 2, y=(ra.y + pa.y) / 2,
            charge=pa.charge, role=ra.role if ra.role != "normal" else pa.role,
        ))

    react_bonds = [Bond(a=sid_r(b.a), b=sid_r(b.b), order=b.order) for b in r_bonds]
    prod_bonds = [Bond(a=sid_p(b.a), b=sid_p(b.b), order=b.order) for b in p_bonds]
    # In the transition frame, show both sets at once: old bonds breaking, new
    # bonds forming, survivors normal.
    prod_keys = {frozenset((b.a, b.b)) for b in prod_bonds}
    react_keys = {frozenset((b.a, b.b)) for b in react_bonds}
    trans_bonds: list[Bond] = []
    for b in react_bonds:
        state = "normal" if frozenset((b.a, b.b)) in prod_keys else "breaking"
        trans_bonds.append(Bond(a=b.a, b=b.b, order=b.order, state=state))
    for b in prod_bonds:
        if frozenset((b.a, b.b)) not in react_keys:
            trans_bonds.append(Bond(a=b.a, b=b.b, order=b.order, state="forming"))

    frames = [
        Frame(id="reactants", t=0.0, caption=caption_reactants,
              atoms=react_frame_atoms, bonds=react_bonds),
        Frame(id="reaction", t=0.5, caption=caption_react,
              atoms=trans_atoms, bonds=trans_bonds),
        Frame(id="products", t=1.0, caption=caption_products,
              atoms=prod_frame_atoms, bonds=prod_bonds),
    ]
    return AnimationSpec(
        reaction_type=reaction_type,
        equation_smiles=eq.to_string(),
        equation_names=eq.to_string(use_names=True),
        frames=frames,
    )
