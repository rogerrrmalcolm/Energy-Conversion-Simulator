#!/usr/bin/env python3

import argparse
import json
import math
import os
import statistics
import subprocess
from dataclasses import dataclass
from pathlib import Path


STATUS_KEYS = (
    "outside_ergosphere",
    "physics_invalid",
    "captured_or_non_escaping",
    "escaping_without_target",
    "integration_failed",
    "goal_feasible",
)


@dataclass(frozen=True)
class SimulationResult:
    backend: str
    status: str
    nodes: int
    counts: tuple[int, ...]
    parameters: tuple[float, float, float]
    extracted_energy: float
    efficiency_percent: float
    maximum_residual: float
    captured_termination: str
    escaping_termination: str
    avx2_batches: int
    elapsed_ms: float


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare scalar and AVX2 phase-map runs")
    parser.add_argument("--scalar", type=Path, required=True)
    parser.add_argument("--avx2", type=Path, required=True)
    parser.add_argument("--scenario", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--relative-tolerance", type=float, default=1e-9)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-10)
    parser.add_argument("--require-speedup", action="store_true")
    return parser.parse_args()


def parse_output(output: str) -> SimulationResult:
    document = json.loads(output)
    if document.get("command") != "map":
        raise ValueError("simulation JSON is not a phase-map result")
    diagnostics = document["diagnostics"]
    candidate = document.get("best_candidate")
    if candidate is None:
        raise ValueError("phase-map JSON does not contain a best candidate")
    split = candidate["split"]
    counts = tuple(int(diagnostics["status_counts"][key]) for key in STATUS_KEYS)
    return SimulationResult(
        backend=str(document["execution_backend"]),
        status=str(document["status"]),
        nodes=int(diagnostics["nodes_evaluated"]),
        counts=counts,
        parameters=(
            float(split["split_radius_over_m"]),
            float(split["incoming_lz_over_m_m"]),
            float(split["split_angle_rad"]),
        ),
        extracted_energy=float(candidate["net_extracted_energy"]),
        efficiency_percent=100.0 * float(candidate["penrose_efficiency"]),
        maximum_residual=float(candidate["maximum_normalized_residual"]),
        captured_termination=str(candidate["captured_termination"]),
        escaping_termination=str(candidate["escaping_termination"]),
        avx2_batches=int(diagnostics["avx2_four_lane_batches"]),
        elapsed_ms=float(diagnostics["elapsed_ms"]),
    )


def execute(executable: Path, scenario: Path) -> SimulationResult:
    completed = subprocess.run(
        [
            str(executable.resolve()),
            "map",
            str(scenario.resolve()),
            "--format",
            "json",
            "--no-progress",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return parse_output(completed.stdout)


def compare_results(
    scalar: SimulationResult,
    avx2: SimulationResult,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> list[str]:
    errors = []
    for field in ("status", "nodes", "counts", "captured_termination", "escaping_termination"):
        if getattr(scalar, field) != getattr(avx2, field):
            errors.append(f"{field} differs: {getattr(scalar, field)} != {getattr(avx2, field)}")
    scalar_values = scalar.parameters + (
        scalar.extracted_energy,
        scalar.efficiency_percent,
        scalar.maximum_residual,
    )
    avx2_values = avx2.parameters + (
        avx2.extracted_energy,
        avx2.efficiency_percent,
        avx2.maximum_residual,
    )
    for index, (left, right) in enumerate(zip(scalar_values, avx2_values)):
        if not math.isclose(left, right, rel_tol=relative_tolerance, abs_tol=absolute_tolerance):
            errors.append(f"numeric result {index} differs: {left} != {right}")
    if scalar.avx2_batches != 0:
        errors.append("scalar executable unexpectedly reported AVX2 batches")
    if avx2.avx2_batches == 0:
        errors.append("AVX2 executable did not dispatch an AVX2 batch")
    return errors


def timed_runs(
    scalar: Path, avx2: Path, scenario: Path, repetitions: int
) -> tuple[list[SimulationResult], list[SimulationResult]]:
    execute(scalar, scenario)
    execute(avx2, scenario)
    results = {"scalar": [], "avx2": []}
    executables = {"scalar": scalar, "avx2": avx2}
    for repetition in range(repetitions):
        order = ("scalar", "avx2") if repetition % 2 == 0 else ("avx2", "scalar")
        for backend in order:
            results[backend].append(execute(executables[backend], scenario))
    return results["scalar"], results["avx2"]


def report(
    scalar_runs: list[SimulationResult], avx2_runs: list[SimulationResult], errors: list[str]
) -> float:
    scalar_ms = statistics.median(run.elapsed_ms for run in scalar_runs)
    avx2_ms = statistics.median(run.elapsed_ms for run in avx2_runs)
    speedup = scalar_ms / avx2_ms
    correctness = "PASS" if not errors else "FAIL"
    lines = [
        "## Scalar versus AVX2",
        "",
        f"Correctness: **{correctness}**",
        "",
        "| Backend | Median simulation time |",
        "| --- | ---: |",
        f"| Scalar batch4 | {scalar_ms:.3f} ms |",
        f"| AVX2 batch4 | {avx2_ms:.3f} ms |",
        "",
        f"Measured AVX2 speedup: **{speedup:.3f}x**",
    ]
    summary = "\n".join(lines)
    print(summary)
    if summary_path := os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(summary_path, "a", encoding="utf-8") as stream:
            stream.write(summary + "\n")
    return speedup


def main() -> int:
    arguments = parse_arguments()
    if arguments.repetitions < 1:
        raise ValueError("repetitions must be positive")
    scalar_runs, avx2_runs = timed_runs(
        arguments.scalar, arguments.avx2, arguments.scenario, arguments.repetitions
    )
    errors = compare_results(
        scalar_runs[0],
        avx2_runs[0],
        arguments.relative_tolerance,
        arguments.absolute_tolerance,
    )
    speedup = report(scalar_runs, avx2_runs, errors)
    for error in errors:
        print(f"ERROR: {error}")
    if errors:
        return 1
    if speedup <= 1.0:
        print("::warning::AVX2 was not faster on this shared runner")
        return 1 if arguments.require_speedup else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
