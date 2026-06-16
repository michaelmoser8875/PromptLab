"""Phase 4 mechanism templates.

`get_template(reaction_type)` returns a callable ``(BalancedEquation, str) ->
AnimationSpec``. SN2 has a dedicated, atom-tracked choreographer; the other
supported types use the generic approach/converge/cross-fade engine, configured
with type-specific captions and (where it helps) actor-atom highlighting.

Every template is parameterized purely by the molecules in the equation, so
swapping reactants produces a new animation without code changes.
"""
from __future__ import annotations

from typing import Callable

from rdkit import Chem

from backend.animation import (
    AnimationSpec,
    ROLE_ACID_PROTON,
    ROLE_BASE_SITE,
    ROLE_ELECTROPHILE,
    ROLE_LEAVING_GROUP,
    ROLE_NUCLEOPHILE,
)
from backend.equation_builder import BalancedEquation, Species
from backend.templates.common import generic_template
from backend.templates.sn2 import sn2_template

Template = Callable[[BalancedEquation, str], AnimationSpec]

_HALOGENS = {"F", "Cl", "Br", "I"}


# ---------------------------------------------------------------------------
# Role helpers (best-effort actor highlighting for the generic engine)
# ---------------------------------------------------------------------------


def _smarts_matches(smiles: str, smarts: str) -> list[tuple[int, ...]]:
    mol = Chem.MolFromSmiles(smiles)
    pat = Chem.MolFromSmarts(smarts)
    if mol is None or pat is None:
        return []
    mol = Chem.AddHs(mol)
    return list(mol.GetSubstructMatches(pat))


def _ester_roles(smiles: str, which: str) -> dict:
    roles: dict = {}
    if which == "reactants":
        for m in _smarts_matches(smiles, "[CX3](=[OX1])[OX2H1]"):
            roles[str(m[0])] = ROLE_ELECTROPHILE   # carbonyl carbon
            roles[str(m[2])] = ROLE_LEAVING_GROUP   # hydroxyl O (-> water)
        for m in _smarts_matches(smiles, "[#6;X4][OX2H1]"):
            roles[str(m[1])] = ROLE_NUCLEOPHILE     # alcohol oxygen
    else:  # products
        for m in _smarts_matches(smiles, "[CX3](=[OX1])[OX2][#6]"):
            roles[str(m[0])] = ROLE_ELECTROPHILE    # ester carbonyl
            roles[str(m[2])] = ROLE_NUCLEOPHILE     # ester bridging O
        for m in _smarts_matches(smiles, "[OX2H2]"):
            roles[str(m[0])] = ROLE_LEAVING_GROUP   # water oxygen
    return roles


def _acid_base_roles(smiles: str, which: str) -> dict:
    """Mark acidic protons and basic sites by structure (both sides)."""
    roles: dict = {}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return roles
    mol = Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        idx = str(atom.GetIdx())
        # basic site: a lone-pair donor (anion or neutral amine N)
        if atom.GetFormalCharge() < 0 and sym in ("O", "N", "S"):
            roles[idx] = ROLE_BASE_SITE
        elif sym == "N" and atom.GetTotalNumHs() >= 1 and atom.GetFormalCharge() == 0:
            if not any(n.GetSymbol() == "C" and any(
                b.GetBondTypeAsDouble() == 2.0 for b in n.GetBonds()
            ) for n in atom.GetNeighbors()):
                roles[idx] = ROLE_BASE_SITE
        # acidic proton: H on a halogen, or H on an O that hangs off an oxoacid
        elif sym == "H":
            heavy = atom.GetNeighbors()[0] if atom.GetDegree() == 1 else None
            if heavy is None:
                continue
            if heavy.GetSymbol() in _HALOGENS:
                roles[idx] = ROLE_ACID_PROTON
            elif heavy.GetSymbol() == "O":
                if any(n.GetSymbol() in ("C", "S", "N", "P")
                       for n in heavy.GetNeighbors()):
                    roles[idx] = ROLE_ACID_PROTON
    return roles


# ---------------------------------------------------------------------------
# Type templates
# ---------------------------------------------------------------------------


def _acid_base(eq: BalancedEquation, rtype: str) -> AnimationSpec:
    return generic_template(
        eq, rtype,
        "An acid and a base meet in solution.",
        "A proton transfers from the acid to the base.",
        "The products: a salt and/or water.",
        role_fn=_acid_base_roles,
    )


def _esterification(eq: BalancedEquation, rtype: str) -> AnimationSpec:
    return generic_template(
        eq, rtype,
        "A carboxylic acid and an alcohol come together (acid-catalyzed).",
        "The alcohol oxygen attacks the carbonyl; the acid's –OH leaves as water.",
        "An ester forms, releasing a molecule of water.",
        role_fn=_ester_roles,
    )


def _combustion(eq: BalancedEquation, rtype: str) -> AnimationSpec:
    return generic_template(
        eq, rtype,
        "Fuel and oxygen are brought together.",
        "Bonds break and oxygen inserts as the fuel is oxidized.",
        "Carbon dioxide and water are released, along with heat.",
    )


def _precipitation(eq: BalancedEquation, rtype: str) -> AnimationSpec:
    return generic_template(
        eq, rtype,
        "Two soluble salts are dissolved together.",
        "The dissolved ions swap partners.",
        "An insoluble salt precipitates out of solution.",
    )


def _generic_default(eq: BalancedEquation, rtype: str) -> AnimationSpec:
    return generic_template(
        eq, rtype,
        "The reactants are brought together.",
        "They react.",
        "The products form.",
    )


_REGISTRY: dict[str, Template] = {
    "sn2": sn2_template,
    "acid-base": _acid_base,
    "esterification": _esterification,
    "combustion": _combustion,
    "precipitation": _precipitation,
}


def get_template(reaction_type: str) -> Template:
    return _REGISTRY.get(reaction_type, _generic_default)
