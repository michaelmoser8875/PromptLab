"""Phase 4 - animation spec validation (offline, no API key needed).

The Phase 4 engine is fully deterministic (RDKit geometry + structural rules),
so this test builds balanced equations by hand and checks the emitted specs:

  - every supported type produces a structurally valid spec;
  - SN2 tracks its actor atoms across frames and contains the breaking/forming
    bond events and electron-pushing arrows of a backside attack;
  - swapping reactants (different alkyl halide + nucleophile) yields a
    structurally analogous SN2 spec with no code change - the swappability the
    phase is built around.

Run:  python tests/test_phase4.py   (exit 0 on success)
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from backend.animation import build_animation  # noqa: E402
from backend.equation_builder import BalancedEquation, Species  # noqa: E402
from backend.gen_examples import EXAMPLES  # noqa: E402

_VALID_BOND_STATES = {"normal", "forming", "breaking"}
_failures: list[str] = []


def check(cond: bool, msg: str) -> None:
    if not cond:
        _failures.append(msg)


def eq(reactants, products) -> BalancedEquation:
    return BalancedEquation(
        reactants=[Species(n, s, 1) for n, s in reactants],
        products=[Species(n, s, 1) for n, s in products],
    )


def validate_structure(name: str, spec) -> None:
    d = spec.to_dict()
    frames = d["frames"]
    check(len(frames) >= 3, f"{name}: expected >=3 frames, got {len(frames)}")
    ts = [f["t"] for f in frames]
    check(ts == sorted(ts), f"{name}: frame times not sorted: {ts}")
    check(abs(ts[0]) < 1e-9, f"{name}: first frame t!=0 ({ts[0]})")
    check(abs(ts[-1] - 1.0) < 1e-9, f"{name}: last frame t!=1 ({ts[-1]})")

    for f in frames:
        ids = [a["id"] for a in f["atoms"]]
        check(len(ids) == len(set(ids)), f"{name}/{f['id']}: duplicate atom ids")
        idset = set(ids)
        for b in f["bonds"]:
            check(b["a"] in idset and b["b"] in idset,
                  f"{name}/{f['id']}: bond references missing atom")
            check(b["state"] in _VALID_BOND_STATES,
                  f"{name}/{f['id']}: bad bond state {b['state']!r}")
        for a in f["atoms"]:
            check(isinstance(a["el"], str) and a["el"].isalpha(),
                  f"{name}/{f['id']}: bad element {a['el']!r}")
        for ar in f["arrows"]:
            check(ar["src"] in idset and ar["dst"] in idset,
                  f"{name}/{f['id']}: arrow references missing atom")


def roles_in(spec) -> set[str]:
    return {a["role"] for f in spec.to_dict()["frames"] for a in f["atoms"]}


def bond_states_in(spec) -> set[str]:
    return {b["state"] for f in spec.to_dict()["frames"] for b in f["bonds"]}


def sn2_signature(spec) -> tuple:
    """A structure-only fingerprint that should match across swapped SN2 inputs."""
    frames = spec.to_dict()["frames"]
    frame_ids = tuple(f["id"] for f in frames)
    role_set = frozenset(r for r in roles_in(spec) if r != "normal")
    state_set = frozenset(s for s in bond_states_in(spec) if s != "normal")
    # actor atoms keep stable ids across all frames -> same id set every frame
    idsets = [frozenset(a["id"] for a in f["atoms"]) for f in frames]
    tracked = all(s == idsets[0] for s in idsets)
    return frame_ids, role_set, state_set, tracked


def main() -> int:
    # 1. Every curated example produces a valid spec.
    for name, rtype, equation in EXAMPLES:
        spec = build_animation(equation, rtype)
        validate_structure(name, spec)
        check(spec.reaction_type == rtype, f"{name}: type mismatch")

    # 2. SN2 mechanism specifics.
    sn2 = build_animation(
        eq([("bromomethane", "CBr"), ("hydroxide", "[OH-]")],
           [("methanol", "CO"), ("bromide", "[Br-]")]),
        "sn2",
    )
    r = roles_in(sn2)
    for role in ("nucleophile", "electrophile", "leaving_group"):
        check(role in r, f"sn2: missing role {role}")
    check("forming" in bond_states_in(sn2), "sn2: no forming bond")
    check("breaking" in bond_states_in(sn2), "sn2: no breaking bond")
    has_arrows = any(f["arrows"] for f in sn2.to_dict()["frames"])
    check(has_arrows, "sn2: no electron-pushing arrows")

    # 3. Swap-invariance: a different alkyl halide + nucleophile, same template.
    sn2_b = build_animation(
        eq([("1-iodopropane", "CCCI"), ("methoxide", "[O-]C")],
           [("methyl propyl ether", "CCCOC"), ("iodide", "[I-]")]),
        "sn2",
    )
    validate_structure("sn2_swap", sn2_b)
    check(sn2_signature(sn2) == sn2_signature(sn2_b),
          "sn2: swapped reactants produced a different mechanism signature")

    # 4. Generic engine produces a tracked morph, not a cross-fade: reactant and
    #    product frames share atom ids (atoms slide) and bonds break and form.
    comb = build_animation(
        eq([("methane", "C"), ("oxygen", "O=O")],
           [("carbon dioxide", "O=C=O"), ("water", "O")]),
        "combustion",
    )
    f = comb.to_dict()["frames"]
    react_ids = {a["id"] for a in f[0]["atoms"]}
    prod_ids = {a["id"] for a in f[-1]["atoms"]}
    shared = react_ids & prod_ids
    check(len(shared) >= min(len(react_ids), len(prod_ids)) // 2,
          "combustion: expected shared (tracked) atom ids across frames")
    states = bond_states_in(comb)
    check("breaking" in states and "forming" in states,
          "combustion: expected bonds to break and form during the morph")

    total = len(EXAMPLES) + 3
    if _failures:
        print(f"FAIL — {len(_failures)} check(s) failed:")
        for m in _failures:
            print(f"  - {m}")
        return 1
    print(f"PASS — {total} spec groups validated "
          f"({len(EXAMPLES)} curated examples + SN2 mechanism + swap + tracked morph)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
