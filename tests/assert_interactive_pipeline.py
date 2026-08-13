#!/usr/bin/env python3

import argparse
import re
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the guided CLI pipeline")
    parser.add_argument("--executable", required=True)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--require-avx2", action="store_true")
    arguments = parser.parse_args()

    completed = subprocess.run(
        [arguments.executable, "--progress"],
        input=arguments.input.read_text(encoding="utf-8"),
        text=True,
        capture_output=True,
        timeout=30,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"pipeline exited with {completed.returncode}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )

    expected_order = [
        "Stage 1/7 - Black-hole inputs",
        "Stage 2/7 - Kerr integration controls",
        "Stage 3/7 - Algebraic Kerr reservoir",
        "Stage 4/7 - Particle and event inputs",
        "Stage 5/7 - Candidate parameter grid",
        "Stage 6/7 - Exhaustive phase-map evaluation",
        "Bounded Penrose phase-space map",
        "Stage 7/7 - Selected Penrose event",
        "Idealized equatorial Penrose event",
    ]
    position = -1
    for marker in expected_order:
        next_position = completed.stdout.find(marker, position + 1)
        if next_position < 0:
            raise AssertionError(f"missing or out-of-order pipeline marker: {marker}")
        position = next_position

    if "Choose an action" in completed.stdout:
        raise AssertionError("the legacy free-order action menu is still present")
    if "Dijkstra Penrose search" in completed.stdout:
        raise AssertionError("the guided pipeline still uses scalar Dijkstra evaluation")
    if arguments.require_avx2 and "completed / avx2-four-lane-single-thread" not in completed.stdout:
        raise AssertionError("guided phase map did not dispatch its four-state batch through AVX2")
    if not re.search(r"Progress: 4/4 \(100\.0%\) \[#{24}\]", completed.stderr):
        raise AssertionError(f"missing final percentage progress update:\n{completed.stderr}")

    print("Validated sequential input stages and percentage progress")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
