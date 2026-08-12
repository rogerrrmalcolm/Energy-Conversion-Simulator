#!/usr/bin/env python3

import argparse
import json
import re
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the CLI and validate its JSON contract")
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--expected-command", required=True)
    parser.add_argument("--expected-status")
    parser.add_argument("--require-avx2", action="store_true")
    parser.add_argument("--require-progress", action="store_true")
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    options = parser.parse_args()

    arguments = options.arguments
    if arguments and arguments[0] == "--":
        arguments = arguments[1:]
    completed = subprocess.run(
        [
            str(options.executable),
            *arguments,
            "--format",
            "json",
            "--progress" if options.require_progress else "--no-progress",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    if options.require_progress:
        if not re.search(r"Progress: \d+/\d+ \(100\.0%\)", completed.stderr):
            raise ValueError(f"JSON command did not report final progress: {completed.stderr}")
    elif completed.stderr:
        raise ValueError(f"JSON command unexpectedly wrote to stderr: {completed.stderr}")

    document = json.loads(completed.stdout)
    if document.get("schema_version") != 1:
        raise ValueError("JSON output does not use schema version 1")
    if document.get("command") != options.expected_command:
        raise ValueError(
            f"expected command {options.expected_command!r}, got {document.get('command')!r}"
        )
    if options.expected_status and document.get("status") != options.expected_status:
        raise ValueError(
            f"expected status {options.expected_status!r}, got {document.get('status')!r}"
        )
    if options.require_avx2:
        if document.get("execution_backend") != "avx2-four-lane-single-thread":
            raise ValueError("JSON output did not select the AVX2 backend")
        if document.get("diagnostics", {}).get("avx2_four_lane_batches", 0) < 1:
            raise ValueError("JSON output did not report an AVX2 batch")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
