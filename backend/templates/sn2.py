"""SN2 mechanism template - backside attack with Walden inversion.

The showcase template. Unlike the generic cross-fade, SN2 atoms keep stable
ids across every frame, so the renderer interpolates them smoothly: the
nucleophile drives in from the face opposite the leaving group, the three
spectator substituents flatten through a planar transition state and then
invert (the umbrella flip), and the leaving group departs as a free ion.

Geometry is deterministic and reads off the actual molecules:
- the electrophilic carbon and its leaving halide are found structurally;
- each spectator substituent keeps its own RDKit-generated internal geometry
  but is fanned around the carbon and inverted as a rigid subtree.

So swapping bromomethane for 1-chlorobutane (or hydroxide for cyanide) changes
the picture without touching this file.
"""
from __future__ import annotations

import math
from collections import deque
from dataclasses import replace
from typing import Optional

from rdkit import Chem

from backend.animation import (
    SCALE,
    AnimationSpec,
    Arrow,
    Atom,
    Bond,
    Frame,
    ROLE_ELECTROPHILE,
    ROLE_LEAVING_GROUP,
    ROLE_NUCLEOPHILE,
    BOND_BREAKING,
    BOND_FORMING,
    layout_molecule,
)
from backend.equation_builder import BalancedEquation, Species

_HALOGENS = {"F", "Cl", "Br", "I"}
_BOND = 1.5 * SCALE  # one display bond length

# Nucleophile carbon distance along the attack axis (-x), per phase.
_NU_FAR = -3.4 * SCALE
_NU_NEAR = -1.9 * SCALE
_NU_BONDED = -_BOND

# Leaving-group carbon distance along +x, per phase.
_LG_START = _BOND
_LG_TS = 1.7 * SCALE
_LG_GONE = 3.4 * SCALE


def _addh_mol(smiles: str) -> Optional[Chem.Mol]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.AddHs(mol)


def _find_center(smiles: str) -> Optional[tuple[int, int]]:
    """Return (carbon_idx, halide_idx) for the SN2 reacting center, or None."""
    mol = _addh_mol(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _HALOGENS:
            continue
        nbrs = atom.GetNeighbors()
        if len(nbrs) != 1:
            continue
        c = nbrs[0]
        if c.GetSymbol() == "C" and c.GetTotalDegree() == 4:
            return c.GetIdx(), atom.GetIdx()
    return None


def _find_nucleophile_atom(smiles: str) -> Optional[int]:
    """Index of the attacking atom: an anionic heteroatom, else an amine N."""
    mol = _addh_mol(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() < 0 and atom.GetSymbol() in ("O", "S", "N", "C"):
            return atom.GetIdx()
    pat = Chem.MolFromSmarts("[NX3;H1,H2,H3;!$(N=*);!$(NC=O)]")
    if pat is not None:
        m = mol.GetSubstructMatches(pat)
        if m:
            return m[0][0]
    return None


def _split_reactants(
    eq: BalancedEquation,
) -> Optional[tuple[Species, int, int, Species, int]]:
    """Identify (alkyl_halide, c_idx, lg_idx, nucleophile, nu_idx)."""
    halide: Optional[tuple[Species, int, int]] = None
    for s in eq.reactants:
        center = _find_center(s.smiles)
        if center is not None:
            halide = (s, center[0], center[1])
            break
    if halide is None:
        return None
    halide_species = halide[0]

    for s in eq.reactants:
        if s is halide_species:
            continue
        nu = _find_nucleophile_atom(s.smiles)
        if nu is not None:
            return halide_species, halide[1], halide[2], s, nu
    return None


def _subtree(mol: Chem.Mol, start: int, blocked: int) -> set[int]:
    """Atoms reachable from `start` without passing through `blocked`."""
    seen = {start, blocked}
    out = {start}
    q = deque([start])
    while q:
        cur = q.popleft()
        for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
            j = nbr.GetIdx()
            if j not in seen:
                seen.add(j)
                out.add(j)
                q.append(j)
    return out


def _build(eq: BalancedEquation, reaction_type: str) -> AnimationSpec:
    split = _split_reactants(eq)
    if split is None:
        raise ValueError("SN2 template: could not locate reacting center")
    halide_sp, c_idx, lg_idx, nu_sp, nu_idx = split

    e_layout = layout_molecule(halide_sp.smiles)
    n_layout = layout_molecule(nu_sp.smiles)
    e_pos = {int(a.id): (a.x, a.y) for a in e_layout.atoms}
    e_meta = {int(a.id): a for a in e_layout.atoms}
    n_meta = {int(a.id): a for a in n_layout.atoms}

    mol = _addh_mol(halide_sp.smiles)
    assert mol is not None

    # Spectator substituents on the carbon (everything bonded but the halide).
    subs = [
        nbr.GetIdx()
        for nbr in mol.GetAtomWithIdx(c_idx).GetNeighbors()
        if nbr.GetIdx() != lg_idx
    ]
    subtrees = {s: _subtree(mol, s, c_idx) for s in subs}

    cx, cy = e_pos[c_idx]

    def fan_positions(angles: list[float]) -> dict[int, tuple[float, float]]:
        """Target position of each substituent's first atom at the given angles."""
        out: dict[int, tuple[float, float]] = {}
        for s, ang in zip(subs, angles):
            out[s] = (_BOND * math.cos(ang), _BOND * math.sin(ang))
        return out

    n = len(subs)
    # Left-leaning umbrella (reactant), planar vertical (TS), right-leaning (product).
    spread = math.radians(58)
    left = [math.radians(180) + (i - (n - 1) / 2) * spread / max(n, 1) for i in range(n)]
    vert = [math.radians(90) + (i - (n - 1) / 2) * math.radians(150) / max(n, 1) for i in range(n)]
    right = [math.pi - a for a in left]

    def electrophile_atoms(
        fan: dict[int, tuple[float, float]], lg_x: float
    ) -> list[Atom]:
        atoms: list[Atom] = []
        # carbon at origin
        c = e_meta[c_idx]
        atoms.append(Atom(id="e" + str(c_idx), el=c.el, x=0.0, y=0.0,
                          charge=c.charge, role=ROLE_ELECTROPHILE))
        # leaving group along +x
        lg = e_meta[lg_idx]
        atoms.append(Atom(id="e" + str(lg_idx), el=lg.el, x=lg_x, y=0.0,
                          charge=lg.charge, role=ROLE_LEAVING_GROUP))
        # each substituent subtree, rigidly shifted so its root hits the fan target
        for s in subs:
            root_now = e_pos[s]
            target = fan[s]
            shift = (target[0] - (root_now[0] - cx), target[1] - (root_now[1] - cy))
            for idx in subtrees[s]:
                px, py = e_pos[idx]
                meta = e_meta[idx]
                atoms.append(Atom(
                    id="e" + str(idx),
                    el=meta.el,
                    x=(px - cx) + shift[0],
                    y=(py - cy) + shift[1],
                    charge=meta.charge,
                ))
        return atoms

    def nucleophile_atoms(nu_x: float, bonded: bool) -> list[Atom]:
        # orient so the attacking atom is at origin and the bulk points -x (left)
        nx, ny = n_meta[nu_idx].x, n_meta[nu_idx].y
        shifted = [(replace(a), a.x - nx, a.y - ny) for a in n_layout.atoms]
        cxn = sum(dx for _, dx, _ in shifted) / len(shifted)
        mirror = cxn > 0
        atoms: list[Atom] = []
        for a, dx, dy in shifted:
            x = -dx if mirror else dx
            atoms.append(Atom(
                id="n" + a.id,
                el=a.el,
                x=nu_x + x,
                y=dy,
                charge=a.charge,
                role=ROLE_NUCLEOPHILE if int(a.id) == nu_idx else "normal",
            ))
        return atoms

    # Static intramolecular bonds (ids namespaced to match the atoms above).
    def e_bonds(state_clg: str) -> list[Bond]:
        out: list[Bond] = []
        for b in e_layout.bonds:
            ia, ib = int(b.a), int(b.b)
            is_clg = {ia, ib} == {c_idx, lg_idx}
            out.append(Bond(a="e" + b.a, b="e" + b.b, order=b.order,
                            state=state_clg if is_clg else "normal"))
        return out

    n_bonds = [Bond(a="n" + b.a, b="n" + b.b, order=b.order) for b in n_layout.bonds]
    nu_c_bond = Bond(a="n" + str(nu_idx), b="e" + str(c_idx), order=1.0)

    def frame(
        fid: str, t: float, caption: str, fan: list[float], nu_x: float, lg_x: float,
        clg_state: str, nu_state: Optional[str], arrows: Optional[list[Arrow]] = None,
    ) -> Frame:
        atoms = electrophile_atoms(fan_positions(fan), lg_x) + nucleophile_atoms(
            nu_x, bonded=nu_state == "normal"
        )
        bonds = e_bonds(clg_state) + n_bonds
        if nu_state is not None:
            bonds.append(replace(nu_c_bond, state=nu_state))
        return Frame(id=fid, t=t, caption=caption, atoms=atoms, bonds=bonds,
                     arrows=arrows or [])

    nu_label = nu_sp.name
    lg_el = e_meta[lg_idx].el
    arrows_ts = [
        Arrow(src="n" + str(nu_idx), dst="e" + str(c_idx)),
        Arrow(src="e" + str(c_idx), dst="e" + str(lg_idx)),
    ]

    frames = [
        frame("reactants", 0.0,
              f"{nu_label} approaches the carbon from the side opposite the "
              f"{lg_el} leaving group.",
              left, _NU_FAR, _LG_START, "normal", None),
        frame("approach", 0.4,
              "Backside attack: the nucleophile's lone pair lines up with the "
              "C–{0} bond.".format(lg_el),
              left, _NU_NEAR, _LG_START, "normal", "forming"),
        frame("transition", 0.68,
              "Transition state: the C–Nu bond half-forms as the "
              f"C–{lg_el} bond half-breaks; the carbon goes planar.",
              vert, _NU_NEAR * 0.7, _LG_TS, BOND_BREAKING, BOND_FORMING, arrows_ts),
        frame("products", 1.0,
              f"Inversion complete: the {lg_el} departs as an ion and the "
              "substituents have flipped through the carbon (Walden inversion).",
              right, _NU_BONDED, _LG_GONE, BOND_BREAKING, "normal"),
    ]

    return AnimationSpec(
        reaction_type=reaction_type,
        equation_smiles=eq.to_string(),
        equation_names=eq.to_string(use_names=True),
        frames=frames,
        duration_ms=10000,
    )


def sn2_template(eq: BalancedEquation, reaction_type: str) -> AnimationSpec:
    return _build(eq, reaction_type)
