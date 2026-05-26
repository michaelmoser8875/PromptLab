"""Phase 2 - Reaction equation builder.

Takes the structured output of Phase 1 (reactants, optional user-named
products, context) and produces an atom-balanced chemical equation:

    BalancedEquation(
        reactants=[Species(name, smiles, coefficient), ...],
        products=[Species(name, smiles, coefficient), ...],
        context="...",
    )

When the user did not name the products, Claude proposes them from the
reactants + context hint. All SMILES are canonicalized through RDKit.
Stoichiometric coefficients are then determined exactly by computing the
null space of the element-balance matrix in rational arithmetic, so the
balance is correct or the call raises BalanceError - never silently wrong.
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import asdict, dataclass, field
from fractions import Fraction
from typing import Optional

import anthropic
from rdkit import Chem
from rdkit import RDLogger

from backend.nl_parser import (
    Compound,
    ParsedReaction,
    ResolutionError,
    _canonicalize,
    parse_reaction,
)

RDLogger.DisableLog("rdApp.*")

MODEL = "claude-opus-4-7"


@dataclass
class Species:
    name: str
    smiles: str
    coefficient: int = 1


@dataclass
class BalancedEquation:
    reactants: list[Species]
    products: list[Species]
    catalysts: list[Compound] = field(default_factory=list)
    context: str = ""

    def to_dict(self) -> dict:
        return {
            "reactants": [asdict(s) for s in self.reactants],
            "products": [asdict(s) for s in self.products],
            "catalysts": [asdict(c) for c in self.catalysts],
            "context": self.context,
        }

    def to_string(self, use_names: bool = False) -> str:
        def render(species: list[Species]) -> str:
            parts = []
            for s in species:
                label = s.name if use_names else s.smiles
                prefix = "" if s.coefficient == 1 else f"{s.coefficient} "
                parts.append(f"{prefix}{label}")
            return " + ".join(parts)

        arrow = "->"
        if self.catalysts:
            cat = ", ".join(
                (c.name if use_names else c.smiles) for c in self.catalysts
            )
            arrow = f"-[{cat}]->"
        return f"{render(self.reactants)} {arrow} {render(self.products)}"


class BalanceError(Exception):
    """Raised when an equation cannot be balanced with positive integer coefficients."""


# ---------------------------------------------------------------------------
# Product proposal (LLM)
# ---------------------------------------------------------------------------

_PRODUCT_SYSTEM = """\
You predict the products of a chemical reaction given its reactants.

You will be told the reactant compounds (name + SMILES) and a short context
hint about the reaction type. Return the expected products as JSON.

For each product, give:
- name: common English name (e.g. "water", "methyl acetate", "sodium chloride")
- smiles: a valid SMILES string for that product

Rules:
- Return only major products. Skip catalysts, solvents, and trace byproducts.
- If a product is an ionic salt (e.g. NaCl, NaOH), give the neutral salt
  SMILES with separated ions (e.g. "[Na+].[Cl-]"), not just one ion.
- Do not include stoichiometric coefficients - one entry per distinct product.
- All atoms present in the reactants must end up somewhere in the products.
- If you cannot predict the products with confidence, return an empty list.
"""

_PRODUCT_SCHEMA = {
    "type": "object",
    "properties": {
        "products": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "name": {"type": "string"},
                    "smiles": {"type": "string"},
                },
                "required": ["name", "smiles"],
                "additionalProperties": False,
            },
        },
    },
    "required": ["products"],
    "additionalProperties": False,
}


def propose_products(
    reactants: list[Compound],
    context: str,
    client: anthropic.Anthropic,
) -> list[Compound]:
    """Ask Claude for the products of a reaction. Returns canonical Compounds."""
    reactant_lines = "\n".join(f"- {c.name}: {c.smiles}" for c in reactants)
    ctx_line = f"\nContext: {context}" if context else ""
    user_msg = f"Reactants:\n{reactant_lines}{ctx_line}\n\nList the products."

    response = client.messages.create(
        model=MODEL,
        max_tokens=2000,
        system=_PRODUCT_SYSTEM,
        messages=[{"role": "user", "content": user_msg}],
        output_config={
            "format": {
                "type": "json_schema",
                "schema": _PRODUCT_SCHEMA,
            },
        },
    )
    text = next(b.text for b in response.content if b.type == "text")
    data = json.loads(text)

    products: list[Compound] = []
    for item in data["products"]:
        canonical = _canonicalize(item["smiles"])
        if canonical is None:
            raise ResolutionError(
                f"Claude proposed invalid SMILES for product {item['name']!r}: "
                f"{item['smiles']!r}"
            )
        products.append(Compound(name=item["name"], smiles=canonical, source="llm"))
    return products


# ---------------------------------------------------------------------------
# Atom counting (RDKit)
# ---------------------------------------------------------------------------


def _count_atoms(smiles: str) -> dict[str, int]:
    """Element counts for a SMILES, including implicit hydrogens."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise BalanceError(f"Cannot count atoms in invalid SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    counts: dict[str, int] = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        counts[sym] = counts.get(sym, 0) + 1
    return counts


def _net_charge(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise BalanceError(f"Cannot compute charge for invalid SMILES: {smiles!r}")
    return Chem.GetFormalCharge(mol)


# ---------------------------------------------------------------------------
# Balancing via rational nullspace
# ---------------------------------------------------------------------------


def _rref(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    """Reduced row echelon form (in place) using exact rationals."""
    if not matrix:
        return matrix
    rows = len(matrix)
    cols = len(matrix[0])
    pivot_row = 0
    for col in range(cols):
        pivot = None
        for r in range(pivot_row, rows):
            if matrix[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pv = matrix[pivot_row][col]
        matrix[pivot_row] = [x / pv for x in matrix[pivot_row]]
        for r in range(rows):
            if r == pivot_row or matrix[r][col] == 0:
                continue
            factor = matrix[r][col]
            matrix[r] = [
                matrix[r][i] - factor * matrix[pivot_row][i] for i in range(cols)
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return matrix


def _nullspace(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    """Basis for the null space of `matrix` (rational entries)."""
    if not matrix:
        return []
    cols = len(matrix[0])
    rref = _rref([row[:] for row in matrix])
    pivot_cols: list[int] = []
    for row in rref:
        for c, v in enumerate(row):
            if v != 0:
                pivot_cols.append(c)
                break
    free_cols = [c for c in range(cols) if c not in pivot_cols]

    basis: list[list[Fraction]] = []
    for free in free_cols:
        vec = [Fraction(0)] * cols
        vec[free] = Fraction(1)
        for r, p in enumerate(pivot_cols):
            vec[p] = -rref[r][free]
        basis.append(vec)
    return basis


def _scale_to_positive_ints(vec: list[Fraction]) -> Optional[list[int]]:
    """Scale a rational vector to its smallest positive-integer multiple."""
    nonzero = [v for v in vec if v != 0]
    if not nonzero:
        return None
    if all(v < 0 for v in nonzero):
        vec = [-v for v in vec]
        nonzero = [-v for v in nonzero]
    if not all(v >= 0 for v in vec):
        return None  # mixed signs - not a physical solution

    denom_lcm = 1
    for v in vec:
        denom_lcm = denom_lcm * v.denominator // math.gcd(denom_lcm, v.denominator)
    ints = [int(v * denom_lcm) for v in vec]
    g = 0
    for x in ints:
        g = math.gcd(g, abs(x))
    if g == 0:
        return None
    return [x // g for x in ints]


def _balance(
    reactants: list[Compound], products: list[Compound]
) -> tuple[list[int], list[int]]:
    """Find positive integer coefficients for reactants and products."""
    if not reactants or not products:
        raise BalanceError("Need at least one reactant and one product to balance.")

    species = reactants + products
    atom_counts = [_count_atoms(s.smiles) for s in species]
    elements = sorted({el for ac in atom_counts for el in ac})

    rows: list[list[Fraction]] = []
    for el in elements:
        row: list[Fraction] = []
        for i, ac in enumerate(atom_counts):
            sign = 1 if i < len(reactants) else -1
            row.append(Fraction(sign * ac.get(el, 0)))
        rows.append(row)

    charge_row: list[Fraction] = []
    for i, s in enumerate(species):
        sign = 1 if i < len(reactants) else -1
        charge_row.append(Fraction(sign * _net_charge(s.smiles)))
    if any(v != 0 for v in charge_row):
        rows.append(charge_row)

    basis = _nullspace(rows)
    if not basis:
        raise BalanceError(
            "Overdetermined system - no nontrivial coefficient solution exists."
        )
    if len(basis) > 1:
        raise BalanceError(
            f"Underdetermined system ({len(basis)} free parameters) - "
            "the proposed product set is ambiguous or incomplete."
        )
    coeffs = _scale_to_positive_ints(basis[0])
    if coeffs is None:
        raise BalanceError(
            "No positive integer coefficient assignment balances the equation."
        )
    return coeffs[: len(reactants)], coeffs[len(reactants):]


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


def build_equation(
    parsed: ParsedReaction,
    client: Optional[anthropic.Anthropic] = None,
) -> BalancedEquation:
    """Predict products (if missing) and produce a balanced equation."""
    if not parsed.reactants:
        raise ValueError("No reactants to balance.")

    products = list(parsed.products)
    if not products:
        if client is None:
            client = anthropic.Anthropic()
        products = propose_products(parsed.reactants, parsed.context, client)
    if not products:
        raise BalanceError("No products predicted for this reaction.")

    r_coeffs, p_coeffs = _balance(parsed.reactants, products)

    # A reactant whose balanced coefficient is 0 is not consumed - it is
    # acting as a catalyst (or spectator / solvent). Pull it out of the
    # main equation and report it separately. Same for products that drop
    # to 0 (the LLM proposed an extraneous species).
    kept_reactants = [
        (c, k) for c, k in zip(parsed.reactants, r_coeffs) if k != 0
    ]
    catalysts = [c for c, k in zip(parsed.reactants, r_coeffs) if k == 0]
    kept_products = [(p, k) for p, k in zip(products, p_coeffs) if k != 0]
    if not kept_reactants:
        raise BalanceError("Every reactant dropped out during balancing.")
    if not kept_products:
        raise BalanceError("All products dropped out during balancing.")

    return BalancedEquation(
        reactants=[
            Species(name=c.name, smiles=c.smiles, coefficient=k)
            for c, k in kept_reactants
        ],
        products=[
            Species(name=p.name, smiles=p.smiles, coefficient=k)
            for p, k in kept_products
        ],
        catalysts=catalysts,
        context=parsed.context,
    )


def _main() -> None:
    prompt = " ".join(sys.argv[1:]).strip()
    if not prompt:
        prompt = "burn methane in oxygen"
        print(f"(no prompt given; using default: {prompt!r})", file=sys.stderr)
    parsed = parse_reaction(prompt)
    eq = build_equation(parsed)
    print(json.dumps(eq.to_dict(), indent=2))
    print()
    print(eq.to_string())
    print(eq.to_string(use_names=True))


if __name__ == "__main__":
    _main()
