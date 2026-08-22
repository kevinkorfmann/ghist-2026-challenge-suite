from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_growth_bgs.spectrum import fit_growth_spectrum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=64)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=7593991)
    args = parser.parse_args()
    files = sorted(args.input_dir.glob("task_*.json"))
    if len(files) != 36:
        raise ValueError(f"expected 36 recovery fits, found {len(files)}")

    grouped: dict[tuple[float, float], list[tuple[dict, np.ndarray]]] = {}
    for path in files:
        row = json.loads(path.read_text())
        truth = row["truth_in_analysis_units"]
        key = (float(truth["generations"]), float(truth["post_growth_scaling"]))
        spectrum = np.load(path.with_suffix(".sfs.npz"))["spectrum"]
        grouped.setdefault(key, []).append((row, spectrum))
    if any(len(group) != 4 for group in grouped.values()):
        raise ValueError("each truth point must have four chromosome-arm replicates")

    results = {}
    for index, ((generations, scaling), group) in enumerate(sorted(grouped.items())):
        combined_spectrum = np.sum([spectrum for _, spectrum in group], axis=0)
        exposure = sum(row["spectrum"]["physical_exposure_bp"] for row, _ in group)
        mutation_rate = group[0][0]["truth_in_analysis_units"]["mutation_rate"]
        exponential = fit_growth_spectrum(
            combined_spectrum,
            sequence_length=exposure,
            mutation_rate=mutation_rate,
            projection=len(combined_spectrum) - 1,
            restarts=args.restarts,
            workers=args.workers,
            maxiter=600,
            seed=args.seed + 1009 * index,
            convention="exponential",
        )
        sudden = fit_growth_spectrum(
            combined_spectrum,
            sequence_length=exposure,
            mutation_rate=mutation_rate,
            projection=len(combined_spectrum) - 1,
            restarts=args.restarts,
            workers=args.workers,
            maxiter=600,
            seed=args.seed + 1009 * index + 1,
            convention="sudden",
        )
        key = f"generations_{generations:g}|scaling_{scaling:g}"
        results[key] = {
            "truth": {
                "generations": generations,
                "post_growth_scaling": scaling,
            },
            "chromosome_arms": sorted(
                row["generator"]["chromosome"] for row, _ in group
            ),
            "physical_exposure_bp": exposure,
            "accepted_sites": sum(row["spectrum"]["accepted_sites"] for row, _ in group),
            "exponential": exponential,
            "sudden": sudden,
            "relative_error": {
                "generations": exponential["generations"] / generations - 1,
                "post_growth_scaling": (
                    exponential["post_growth_scaling"] / scaling - 1
                ),
            },
            "log_likelihood_exponential_minus_sudden": (
                exponential["multinomial_log_likelihood"]
                - sudden["multinomial_log_likelihood"]
            ),
        }
        print(json.dumps({key: results[key]}, indent=2), flush=True)
    payload = {
        "source": str(args.input_dir),
        "truth_points": len(results),
        "replicates_per_truth": 4,
        "results": results,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
