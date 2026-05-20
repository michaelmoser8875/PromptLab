"""Phase 1 smoke runner.

Runs `parse_reaction` against every prompt in `reaction_prompts.json`,
checking that each one resolves without exceptions and that every returned
SMILES round-trips through RDKit.

Requires ANTHROPIC_API_KEY in the environment. Without it, the runner prints
a warning and exits 0 so it doesn't fail in CI that lacks the secret.

Usage:
    python tests/test_phase1.py
"""
from __future__ import annotations

import json
import os
import sys
import time
import traceback
from pathlib import Path

# Allow `python tests/test_phase1.py` from anywhere.
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from rdkit import Chem  # noqa: E402

from backend.nl_parser import ResolutionError, parse_reaction  # noqa: E402


PROMPTS_FILE = Path(__file__).with_name("reaction_prompts.json")


def smiles_roundtrips(smiles: str) -> bool:
    return Chem.MolFromSmiles(smiles) is not None


def run_one(prompt: str) -> tuple[bool, str]:
    try:
        result = parse_reaction(prompt)
    except ResolutionError as e:
        return False, f"resolution error: {e}"
    except Exception as e:  # pragma: no cover - unexpected failures
        return False, f"unexpected {type(e).__name__}: {e}"

    if not result.reactants:
        return False, "no reactants extracted"

    for c in result.reactants + result.products:
        if not smiles_roundtrips(c.smiles):
            return False, f"invalid SMILES for {c.name!r}: {c.smiles!r}"

    reactant_summary = ", ".join(
        f"{c.name}={c.smiles} ({c.source})" for c in result.reactants
    )
    return True, reactant_summary


def main() -> int:
    if not os.environ.get("ANTHROPIC_API_KEY"):
        print(
            "ANTHROPIC_API_KEY is not set - skipping live test run.",
            file=sys.stderr,
        )
        return 0

    cases = json.loads(PROMPTS_FILE.read_text(encoding="utf-8"))["prompts"]
    passed = 0
    failures: list[tuple[str, str]] = []

    for i, case in enumerate(cases, start=1):
        prompt = case["prompt"]
        start = time.monotonic()
        ok, info = run_one(prompt)
        elapsed = time.monotonic() - start
        marker = "PASS" if ok else "FAIL"
        print(f"[{i:2d}/{len(cases)}] {marker} ({elapsed:5.1f}s) {prompt}")
        if ok:
            print(f"        -> {info}")
            passed += 1
        else:
            print(f"        -> {info}")
            failures.append((prompt, info))

    print()
    print(f"Result: {passed}/{len(cases)} passed")
    if failures:
        print("Failures:")
        for prompt, info in failures:
            print(f"  - {prompt}")
            print(f"      {info}")

    # Phase 1 "done when" is >=20 of the 22 prompts pass.
    return 0 if passed >= 20 else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        sys.exit(130)
    except Exception:
        traceback.print_exc()
        sys.exit(2)
