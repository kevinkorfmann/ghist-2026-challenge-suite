from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_growth_bgs.spectrum import fit_growth_spectrum
from ghist_growth_bgs.simulation import (
    SIMULATION_SCALING_SCHEME,
    require_physical_scaling,
)


def load_group(files: list[Path]) -> list[tuple[dict, np.ndarray]]:
    return [
        (
            json.loads(path.read_text()),
            np.load(path.with_suffix(".sfs.npz"))["spectrum"],
        )
        for path in files
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stress-dir", type=Path, required=True)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=64)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    stress_files = sorted(
        args.stress_dir.glob("task_*.json"),
        key=lambda path: int(path.stem.rsplit("_", 1)[1]),
    )
    if len(stress_files) != 12:
        raise ValueError(f"expected 12 stress fits, found {len(stress_files)}")
    groups = {
        "scale_5|burnin_0.1": load_group(stress_files[0:4]),
        "scale_20|burnin_0.1": load_group(stress_files[4:8]),
        "scale_10|burnin_0.5": load_group(stress_files[8:12]),
        "scale_10|burnin_0.1": load_group(
            [args.baseline_dir / f"task_{index}.json" for index in range(16, 20)]
        ),
    }
    results = {}
    truths = [
        row["truth_in_analysis_units"] for group in groups.values() for row, _ in group
    ]
    if any(truth != truths[0] for truth in truths[1:]):
        raise ValueError(
            "stress comparison requires identical physical demographic truth"
        )
    truth = truths[0]
    for index, (label, group) in enumerate(groups.items()):
        require_physical_scaling([row for row, _ in group])
        spectrum = np.sum([value for _, value in group], axis=0)
        exposure = sum(row["spectrum"]["physical_exposure_bp"] for row, _ in group)
        mutation_rate = group[0][0]["truth_in_analysis_units"]["mutation_rate"]
        fit = fit_growth_spectrum(
            spectrum,
            sequence_length=exposure,
            mutation_rate=mutation_rate,
            projection=len(spectrum) - 1,
            restarts=args.restarts,
            workers=args.workers,
            maxiter=600,
            seed=7596000 + 1009 * index,
            convention="exponential",
        )
        results[label] = {
            "fit": fit,
            "relative_error": {
                "generations": fit["generations"] / truth["generations"] - 1,
                "post_growth_scaling": fit["post_growth_scaling"]
                / truth["post_growth_scaling"]
                - 1,
            },
            "chromosome_arms": sorted(
                row["generator"]["chromosome"] for row, _ in group
            ),
            "physical_exposure_bp": exposure,
            "accepted_sites": sum(
                row["spectrum"]["accepted_sites"] for row, _ in group
            ),
        }
        print(json.dumps({label: results[label]}, indent=2), flush=True)
    payload = {
        "simulation_scaling_scheme": SIMULATION_SCALING_SCHEME,
        "truth": truth,
        "results": results,
        "maximum_absolute_relative_error": {
            key: max(abs(row["relative_error"][key]) for row in results.values())
            for key in ("generations", "post_growth_scaling")
        },
        "freeze_gate_threshold": {
            "generations": 0.20,
            "post_growth_scaling": 0.05,
        },
        "freeze_gate_passed": all(
            max(abs(row["relative_error"][key]) for row in results.values()) < threshold
            for key, threshold in {
                "generations": 0.20,
                "post_growth_scaling": 0.05,
            }.items()
        ),
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
