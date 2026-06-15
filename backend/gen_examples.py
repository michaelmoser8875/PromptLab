"""Generate example animation specs for the frontend, offline (no API calls).

Builds a balanced equation by hand for one representative reaction of each
supported type (plus an SN2 variant, to demonstrate that swapping reactants
re-choreographs the animation with no code change), runs it through the Phase 4
template engine, and writes the specs to frontend/specs/.

    python -m backend.gen_examples
"""
from __future__ import annotations

import json
from pathlib import Path

from backend.animation import build_animation
from backend.equation_builder import BalancedEquation, Species

REPO_ROOT = Path(__file__).resolve().parent.parent
SPECS_DIR = REPO_ROOT / "frontend" / "specs"


def _eq(reactants, products) -> BalancedEquation:
    return BalancedEquation(
        reactants=[Species(n, s, c) for n, s, c in reactants],
        products=[Species(n, s, c) for n, s, c in products],
    )


# (filename, reaction_type, equation)
EXAMPLES = [
    (
        "sn2_bromomethane_hydroxide",
        "sn2",
        _eq(
            [("bromomethane", "CBr", 1), ("hydroxide", "[OH-]", 1)],
            [("methanol", "CO", 1), ("bromide", "[Br-]", 1)],
        ),
    ),
    (
        "sn2_chlorobutane_cyanide",
        "sn2",
        _eq(
            [("1-chlorobutane", "CCCCCl", 1), ("cyanide", "[C-]#N", 1)],
            [("pentanenitrile", "CCCCC#N", 1), ("chloride", "[Cl-]", 1)],
        ),
    ),
    (
        "acid_base_naoh_hcl",
        "acid-base",
        _eq(
            [("sodium hydroxide", "[Na+].[OH-]", 1), ("hydrochloric acid", "Cl", 1)],
            [("sodium chloride", "[Na+].[Cl-]", 1), ("water", "O", 1)],
        ),
    ),
    (
        "esterification_acetic_methanol",
        "esterification",
        _eq(
            [("acetic acid", "CC(=O)O", 1), ("methanol", "CO", 1)],
            [("methyl acetate", "CC(=O)OC", 1), ("water", "O", 1)],
        ),
    ),
    (
        "combustion_methane",
        "combustion",
        _eq(
            [("methane", "C", 1), ("oxygen", "O=O", 2)],
            [("carbon dioxide", "O=C=O", 1), ("water", "O", 2)],
        ),
    ),
    (
        "precipitation_agno3_nacl",
        "precipitation",
        _eq(
            [
                ("silver nitrate", "[Ag+].[O-][N+](=O)[O-]", 1),
                ("sodium chloride", "[Na+].[Cl-]", 1),
            ],
            [
                ("silver chloride", "[Ag+].[Cl-]", 1),
                ("sodium nitrate", "[Na+].[O-][N+](=O)[O-]", 1),
            ],
        ),
    ),
]


def main() -> None:
    SPECS_DIR.mkdir(parents=True, exist_ok=True)
    index = []
    for name, rtype, eq in EXAMPLES:
        spec = build_animation(eq, rtype)
        path = SPECS_DIR / f"{name}.json"
        path.write_text(json.dumps(spec.to_dict(), indent=2), encoding="utf-8")
        index.append({"file": f"{name}.json", "label": spec.equation_names,
                      "type": rtype})
        print(f"wrote {path.relative_to(REPO_ROOT)}  ({len(spec.frames)} frames)")
    (SPECS_DIR / "index.json").write_text(
        json.dumps({"specs": index}, indent=2), encoding="utf-8"
    )
    print(f"wrote {(SPECS_DIR / 'index.json').relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
