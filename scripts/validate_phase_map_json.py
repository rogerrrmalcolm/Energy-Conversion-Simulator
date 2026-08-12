#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate an AVX2 phase-map JSON result")
    parser.add_argument("result", type=Path)
    arguments = parser.parse_args()

    # PowerShell may add a UTF-8 BOM when redirecting output on Windows.
    with arguments.result.open(encoding="utf-8-sig") as stream:
        document = json.load(stream)

    if document.get("schema_version") != 1 or document.get("command") != "map":
        raise ValueError("expected a version-1 phase-map result")
    if document.get("status") != "completed" or document.get("complete") is not True:
        raise ValueError("phase map did not complete")
    if document.get("execution_backend") != "avx2-four-lane-single-thread":
        raise ValueError("phase map did not report the AVX2 backend")

    diagnostics = document.get("diagnostics", {})
    if int(diagnostics.get("avx2_four_lane_batches", 0)) < 1:
        raise ValueError("phase map did not dispatch an AVX2 batch")
    if int(diagnostics.get("nodes_evaluated", 0)) < 1:
        raise ValueError("phase map did not evaluate any nodes")

    print(
        "Validated phase map: "
        f"{diagnostics['nodes_evaluated']} nodes, "
        f"{diagnostics['avx2_four_lane_batches']} AVX2 batches"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
