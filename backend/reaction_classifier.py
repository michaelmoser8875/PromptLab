"""Phase 3 - Reaction type classifier.

Takes a BalancedEquation (from Phase 2) and classifies it into one of a
small set of supported reaction templates. Two-pass design:

  1. Claude proposes a label (one of SUPPORTED_TYPES, or "unsupported").
  2. A rule-based check using RDKit SMARTS verifies that the structural
     fingerprint of the equation matches the proposed label. If the rule
     check fails, the result is downgraded to "unsupported" - the LLM
     cannot override structural evidence.

The supported set is intentionally small: each entry has explicit
structural requirements documented alongside its verifier. Reactions that
don't fit any template fail gracefully with reaction_type="unsupported".
"""
from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass
from typing import Callable, Optional

import anthropic
from rdkit import Chem
from rdkit import RDLogger

from backend.equation_builder import BalancedEquation, Species, build_equation
from backend.nl_parser import parse_reaction

RDLogger.DisableLog("rdApp.*")

MODEL = "claude-opus-4-7"

UNSUPPORTED = "unsupported"
SUPPORTED_TYPES = (
    "sn2",
    "acid-base",
    "esterification",
    "combustion",
    "precipitation",
)


@dataclass
class Classification:
    reaction_type: str  # one of SUPPORTED_TYPES, or UNSUPPORTED
    proposed_by_llm: str
    verified: bool
    reason: str

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# SMARTS helpers
# ---------------------------------------------------------------------------

_PATTERNS: dict[str, str] = {
    # SN2
    "alkyl_halide": "[CX4;!$(C=*)][F,Cl,Br,I]",
    "anion_nucleophile": "[O-,S-,N-,C-;!$(*=*)]",
    "neutral_amine_nucleophile": "[NX3;H2,H3;!$(N=*);!$(NC=O);!$(N[!#1;!#6])]",
    "free_halide_ion": "[F-,Cl-,Br-,I-]",
    # Esterification
    "carboxylic_acid": "[CX3](=[OX1])[OX2H1]",
    "alcohol": "[#6;X4][OX2H1]",
    "ester": "[CX3](=[OX1])[OX2][#6]",
    # Combustion / common
    "molecular_oxygen": "[OX1]=[OX1]",
    "carbon_dioxide": "O=C=O",
    "water": "[OX2H2]",
    # Acid-base
    # Note: no X<n> constraint here - in RDKit SMARTS, X counts non-implicit
    # neighbors, which is 0 for free ions like [OH-] and [O-2].
    "hydroxide": "[OH-]",
    "oxide_ion": "[O-2]",
    "carbonate": "[CX3](=[OX1])([OX1-])[OX1-]",
    "bicarbonate": "[CX3](=[OX1])([OX2H1])[OX1-]",
    "ammonia_or_primary_amine": "[NX3;H2,H3;!$(N=*);!$(NC=O)]",
    "hydrohalic_acid": "[F,Cl,Br,I;X1;H1]",
    "sulfuric_or_sulfonic_oh": "[#16X4](=[OX1])(=[OX1])[OX2H1]",
    "phosphoric_oh": "[#15X4](=[OX1])[OX2H1]",
    "nitric_oh_zwitterion": "[NX3+](=[OX1])([OX1-])[OX2H1]",
}

_compiled: dict[str, Chem.Mol] = {}


def _pattern(name: str) -> Chem.Mol:
    if name not in _compiled:
        pat = Chem.MolFromSmarts(_PATTERNS[name])
        if pat is None:
            raise RuntimeError(f"Bad SMARTS for {name!r}: {_PATTERNS[name]!r}")
        _compiled[name] = pat
    return _compiled[name]


def _mol(smiles: str) -> Optional[Chem.Mol]:
    return Chem.MolFromSmiles(smiles)


def _matches(species: list[Species], pattern_name: str) -> bool:
    pat = _pattern(pattern_name)
    for s in species:
        m = _mol(s.smiles)
        if m is not None and m.HasSubstructMatch(pat):
            return True
    return False


def _matches_any(species: list[Species], pattern_names: list[str]) -> bool:
    return any(_matches(species, name) for name in pattern_names)


# Periodic-table metals (everything not in this nonmetal/metalloid set). Used
# to spot metal oxides and ionic salts.
_NONMETALS = {
    "H", "He",
    "C", "N", "O", "F", "Ne",
    "P", "S", "Cl", "Ar",
    "Se", "Br", "Kr",
    "I", "Xe",
    "Rn",
    "B", "Si", "Ge", "As", "Sb", "Te", "At",
}


def _has_metal_atom(smiles: str) -> bool:
    m = _mol(smiles)
    if m is None:
        return False
    return any(a.GetSymbol() not in _NONMETALS for a in m.GetAtoms())


def _has_ion(smiles: str) -> bool:
    m = _mol(smiles)
    if m is None:
        return False
    return any(a.GetFormalCharge() != 0 for a in m.GetAtoms())


def _has_metal_oxide_product(products: list[Species]) -> bool:
    """A product containing both a metal and an oxygen (covers MgO, Fe2O3, etc.)."""
    for s in products:
        m = _mol(s.smiles)
        if m is None:
            continue
        syms = {a.GetSymbol() for a in m.GetAtoms()}
        if "O" in syms and any(sym not in _NONMETALS for sym in syms):
            return True
    return False


# ---------------------------------------------------------------------------
# Per-type verifiers
# ---------------------------------------------------------------------------


def _verify_sn2(eq: BalancedEquation) -> tuple[bool, str]:
    if not _matches(eq.reactants, "alkyl_halide"):
        return False, "no alkyl halide reactant"
    nu = _matches(eq.reactants, "anion_nucleophile") or _matches(
        eq.reactants, "neutral_amine_nucleophile"
    )
    if not nu:
        return False, "no nucleophile reactant"
    # The halide should have been displaced - either the alkyl halide is
    # gone from the products, or a free halide ion is present.
    if _matches(eq.products, "alkyl_halide") and not _matches(
        eq.products, "free_halide_ion"
    ):
        return False, "alkyl halide persists with no free halide ion in products"
    return True, "alkyl halide + nucleophile -> substitution product"


def _verify_esterification(eq: BalancedEquation) -> tuple[bool, str]:
    if not _matches(eq.reactants, "carboxylic_acid"):
        return False, "no carboxylic acid reactant"
    if not _matches(eq.reactants, "alcohol"):
        return False, "no alcohol reactant"
    if not _matches(eq.products, "ester"):
        return False, "no ester product"
    if not _matches(eq.products, "water"):
        return False, "no water product"
    return True, "acid + alcohol -> ester + water"


def _verify_combustion(eq: BalancedEquation) -> tuple[bool, str]:
    if not _matches(eq.reactants, "molecular_oxygen"):
        return False, "no O2 reactant"
    has_oxide_product = (
        _matches(eq.products, "carbon_dioxide")
        or _matches(eq.products, "water")
        or _has_metal_oxide_product(eq.products)
    )
    if not has_oxide_product:
        return False, "no oxide / CO2 / water product"
    return True, "fuel + O2 -> oxide products"


_ACID_PATTERNS = [
    "carboxylic_acid",
    "hydrohalic_acid",
    "sulfuric_or_sulfonic_oh",
    "phosphoric_oh",
    "nitric_oh_zwitterion",
]

_BASE_PATTERNS = [
    "hydroxide",
    "oxide_ion",
    "carbonate",
    "bicarbonate",
    "ammonia_or_primary_amine",
]


def _verify_acid_base(eq: BalancedEquation) -> tuple[bool, str]:
    has_acid = _matches_any(eq.reactants, _ACID_PATTERNS)
    has_base = _matches_any(eq.reactants, _BASE_PATTERNS)
    if not has_acid:
        return False, "no acidic reactant"
    if not has_base:
        return False, "no basic reactant"
    has_salt_or_water = (
        _matches(eq.products, "water")
        or any(_has_ion(s.smiles) for s in eq.products)
    )
    if not has_salt_or_water:
        return False, "no water or ionic salt product"
    return True, "acid + base -> salt and/or water"


def _verify_precipitation(eq: BalancedEquation) -> tuple[bool, str]:
    if len(eq.reactants) < 2:
        return False, "precipitation needs two ionic reactants"
    if not all(_has_ion(s.smiles) for s in eq.reactants):
        return False, "not all reactants are ionic"
    if not any(_has_ion(s.smiles) for s in eq.products):
        return False, "no ionic product"
    # Classic double-displacement signature: product ions must be a
    # rearrangement of reactant ions (no new element appearing on the right).
    def ion_symbols(species: list[Species]) -> set[str]:
        out: set[str] = set()
        for s in species:
            m = _mol(s.smiles)
            if m is None:
                continue
            for a in m.GetAtoms():
                if a.GetFormalCharge() != 0:
                    out.add(a.GetSymbol())
        return out

    if not ion_symbols(eq.products).issubset(ion_symbols(eq.reactants)):
        return False, "product ions are not a rearrangement of reactant ions"
    return True, "ionic exchange between two salts"


_VERIFIERS: dict[str, Callable[[BalancedEquation], tuple[bool, str]]] = {
    "sn2": _verify_sn2,
    "acid-base": _verify_acid_base,
    "esterification": _verify_esterification,
    "combustion": _verify_combustion,
    "precipitation": _verify_precipitation,
}


# ---------------------------------------------------------------------------
# LLM proposal
# ---------------------------------------------------------------------------

_PROPOSAL_SYSTEM = """\
You classify chemical reactions into one of these labels:

- sn2:            bimolecular nucleophilic substitution at a saturated carbon
                  (e.g. CH3Br + OH- -> CH3OH + Br-)
- acid-base:      Bronsted acid-base proton transfer
                  (e.g. NaOH + HCl -> NaCl + H2O, NH3 + HCl -> NH4Cl)
- esterification: carboxylic acid + alcohol -> ester + water
- combustion:     fuel + O2 -> oxide products
- precipitation:  aqueous double-displacement that yields an insoluble salt
                  (e.g. AgNO3 + NaCl -> AgCl(s) + NaNO3)
- unsupported:    anything else (addition, elimination, EAS, redox between
                  metals, decomposition, dehydration, hydration, etc.)

You will be given a balanced chemical equation (with compound names and SMILES).
Return the single best label. If the reaction does not cleanly fit one of the
supported types, return "unsupported". Do not guess - "unsupported" is the
correct answer when nothing fits.
"""

_PROPOSAL_SCHEMA = {
    "type": "object",
    "properties": {
        "reaction_type": {
            "type": "string",
            "enum": list(SUPPORTED_TYPES) + [UNSUPPORTED],
        },
    },
    "required": ["reaction_type"],
    "additionalProperties": False,
}


def _propose_label(eq: BalancedEquation, client: anthropic.Anthropic) -> str:
    lines = ["Equation:", eq.to_string(use_names=True), eq.to_string()]
    if eq.context:
        lines.append(f"Context hint from user: {eq.context}")
    user_msg = "\n".join(lines) + "\n\nClassify this reaction."
    response = client.messages.create(
        model=MODEL,
        max_tokens=300,
        system=_PROPOSAL_SYSTEM,
        messages=[{"role": "user", "content": user_msg}],
        output_config={
            "format": {
                "type": "json_schema",
                "schema": _PROPOSAL_SCHEMA,
            },
        },
    )
    text = next(b.text for b in response.content if b.type == "text")
    return json.loads(text)["reaction_type"]


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


def classify(
    eq: BalancedEquation,
    client: Optional[anthropic.Anthropic] = None,
) -> Classification:
    """Classify a balanced equation into a supported reaction type or UNSUPPORTED.

    LLM proposes; SMARTS verifies. A proposal that fails verification is
    downgraded to UNSUPPORTED so we never hallucinate a structural match.
    """
    if client is None:
        client = anthropic.Anthropic()

    proposed = _propose_label(eq, client)
    if proposed == UNSUPPORTED:
        return Classification(
            reaction_type=UNSUPPORTED,
            proposed_by_llm=UNSUPPORTED,
            verified=False,
            reason="LLM declined to assign a supported label",
        )

    verifier = _VERIFIERS.get(proposed)
    if verifier is None:
        return Classification(
            reaction_type=UNSUPPORTED,
            proposed_by_llm=proposed,
            verified=False,
            reason=f"no verifier registered for {proposed!r}",
        )

    ok, why = verifier(eq)
    if not ok:
        return Classification(
            reaction_type=UNSUPPORTED,
            proposed_by_llm=proposed,
            verified=False,
            reason=f"verification failed for {proposed!r}: {why}",
        )
    return Classification(
        reaction_type=proposed,
        proposed_by_llm=proposed,
        verified=True,
        reason=why,
    )


def _main() -> None:
    prompt = " ".join(sys.argv[1:]).strip()
    if not prompt:
        prompt = "react bromomethane with hydroxide"
        print(f"(no prompt given; using default: {prompt!r})", file=sys.stderr)
    parsed = parse_reaction(prompt)
    eq = build_equation(parsed)
    print(eq.to_string(use_names=True))
    result = classify(eq)
    print(json.dumps(result.to_dict(), indent=2))


if __name__ == "__main__":
    _main()
