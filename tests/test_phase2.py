"""Phase 2 smoke runner.

For every prompt in `reaction_prompts.json`, runs the full pipeline
(parse_reaction -> build_equation) and verifies:

  - the equation builds without raising
  - element counts balance on both sides (with stoichiometry)
  - net charge balances on both sides
  - the human-readable string renders cleanly

Requires ANTHROPIC_API_KEY. Without it, prints a warning and exits 0 so it
doesn't fail in CI that lacks the secret.

Usage:
    python tests/test_phase2.py
"""
from __future__ import annotations

import json
import os
import sys
import time
import traceback
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from backend.equation_builder import (  # noqa: E402
    BalanceError,
    BalancedEquation,
    _count_atoms,
    _net_charge,
    build_equation,
)
from backend.nl_parser import ResolutionError, parse_reaction  # noqa: E402


PROMPTS_FILE = Path(__file__).with_name("reaction_prompts.json")


def verify_balance(eq: BalancedEquation) -> tuple[bool, str]:
    """Check that atom and charge totals match on both sides."""
    left: dict[str, int] = {}
    right: dict[str, int] = {}
    for s in eq.reactants:
        for el, n in _count_atoms(s.smiles).items():
            left[el] = left.get(el, 0) + s.coefficient * n
    for s in eq.products:
        for el, n in _count_atoms(s.smiles).items():
            right[el] = right.get(el, 0) + s.coefficient * n
    if left != right:
        return False, f"atom mismatch: left={left} right={right}"

    lq = sum(s.coefficient * _net_charge(s.smiles) for s in eq.reactants)
    rq = sum(s.coefficient * _net_charge(s.smiles) for s in eq.products)
    if lq != rq:
        return False, f"charge mismatch: left={lq} right={rq}"
    return True, "balanced"


def run_one(prompt: str) -> tuple[bool, str]:
    try:
        parsed = parse_reaction(prompt)
    except ResolutionError as e:
        return False, f"phase-1 resolution error: {e}"
    except Exception as e:
        return False, f"phase-1 unexpected {type(e).__name__}: {e}"

    try:
        eq = build_equation(parsed)
    except (BalanceError, ResolutionError) as e:
        return False, f"phase-2 {type(e).__name__}: {e}"
    except Exception as e:
        return False, f"phase-2 unexpected {type(e).__name__}: {e}"

    ok, why = verify_balance(eq)
    if not ok:
        return False, why

    rendered = eq.to_string(use_names=True)
    if not rendered or "->" not in rendered:
        return False, f"render failed: {rendered!r}"
    return True, rendered


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
        print(f"        -> {info}")
        if ok:
            passed += 1
        else:
            failures.append((prompt, info))

    print()
    print(f"Result: {passed}/{len(cases)} balanced")
    if failures:
        print("Failures:")
        for prompt, info in failures:
            print(f"  - {prompt}")
            print(f"      {info}")

    # Phase 2 bar: most prompts in the curated set should balance. A handful
    # of unusual cases (e.g. iron rust, dehydration with H2SO4 as catalyst)
    # are expected to be flaky because the LLM has to predict products in
    # contexts where multiple plausible answers exist.
    return 0 if passed >= 18 else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        sys.exit(130)
    except Exception:
        traceback.print_exc()
        sys.exit(2)
