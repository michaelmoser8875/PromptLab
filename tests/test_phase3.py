"""Phase 3 smoke runner.

For every prompt in `reaction_prompts.json` that carries an
`expected_classification` field, runs the full pipeline
(parse_reaction -> build_equation -> classify) and verifies the returned
label matches.

Requires ANTHROPIC_API_KEY. Without it, prints a warning and exits 0 so it
doesn't fail in CI that lacks the secret.

Usage:
    python tests/test_phase3.py
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

from backend.equation_builder import BalanceError, build_equation  # noqa: E402
from backend.nl_parser import ResolutionError, parse_reaction  # noqa: E402
from backend.reaction_classifier import (  # noqa: E402
    SUPPORTED_TYPES,
    UNSUPPORTED,
    Classification,
    classify,
)

PROMPTS_FILE = Path(__file__).with_name("reaction_prompts.json")


def run_one(prompt: str) -> tuple[bool, str, Classification | None]:
    try:
        parsed = parse_reaction(prompt)
    except ResolutionError as e:
        return False, f"phase-1 resolution error: {e}", None
    except Exception as e:
        return False, f"phase-1 unexpected {type(e).__name__}: {e}", None

    try:
        eq = build_equation(parsed)
    except (BalanceError, ResolutionError) as e:
        return False, f"phase-2 {type(e).__name__}: {e}", None
    except Exception as e:
        return False, f"phase-2 unexpected {type(e).__name__}: {e}", None

    try:
        result = classify(eq)
    except Exception as e:
        return False, f"phase-3 unexpected {type(e).__name__}: {e}", None
    return True, eq.to_string(use_names=True), result


def main() -> int:
    if not os.environ.get("ANTHROPIC_API_KEY"):
        print(
            "ANTHROPIC_API_KEY is not set - skipping live test run.",
            file=sys.stderr,
        )
        return 0

    cases = json.loads(PROMPTS_FILE.read_text(encoding="utf-8"))["prompts"]
    cases = [c for c in cases if "expected_classification" in c]

    valid_labels = set(SUPPORTED_TYPES) | {UNSUPPORTED}
    for c in cases:
        if c["expected_classification"] not in valid_labels:
            print(
                f"FATAL: prompt {c['prompt']!r} has invalid "
                f"expected_classification {c['expected_classification']!r}",
                file=sys.stderr,
            )
            return 2

    passed = 0
    failures: list[tuple[str, str, str]] = []  # prompt, expected, info

    for i, case in enumerate(cases, start=1):
        prompt = case["prompt"]
        expected = case["expected_classification"]
        start = time.monotonic()
        ok, eq_str, result = run_one(prompt)
        elapsed = time.monotonic() - start

        if not ok:
            print(f"[{i:2d}/{len(cases)}] FAIL ({elapsed:5.1f}s) {prompt}")
            print(f"        -> {eq_str}")
            failures.append((prompt, expected, eq_str))
            continue

        assert result is not None
        got = result.reaction_type
        match = got == expected
        marker = "PASS" if match else "FAIL"
        print(f"[{i:2d}/{len(cases)}] {marker} ({elapsed:5.1f}s) {prompt}")
        print(f"        eq:       {eq_str}")
        print(
            f"        expected: {expected}    "
            f"got: {got}    proposed: {result.proposed_by_llm}    "
            f"verified: {result.verified}"
        )
        print(f"        reason:   {result.reason}")
        if match:
            passed += 1
        else:
            failures.append(
                (
                    prompt,
                    expected,
                    f"got={got} proposed={result.proposed_by_llm} "
                    f"verified={result.verified} reason={result.reason}",
                )
            )

    total = len(cases)
    print()
    print(f"Result: {passed}/{total} classified correctly")
    if failures:
        print("Failures:")
        for prompt, expected, info in failures:
            print(f"  - {prompt}")
            print(f"      expected={expected}    {info}")

    # Phase 3 "done when" is >=90% accuracy on the curated set.
    threshold = int(round(0.90 * total))
    return 0 if passed >= threshold else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        sys.exit(130)
    except Exception:
        traceback.print_exc()
        sys.exit(2)
