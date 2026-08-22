from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

from ghist_secondary_contact.features import inverse_transform_targets, transform_targets


def _point_rrmse(prediction: np.ndarray, truth: np.ndarray, point_id: np.ndarray) -> float:
    estimates, targets = [], []
    for identifier in np.unique(point_id):
        selected = point_id == identifier
        estimates.append(prediction[selected].mean(axis=0))
        targets.append(truth[selected][0])
    relative = (np.asarray(estimates) - np.asarray(targets)) / np.asarray(targets)
    return float(np.sqrt(np.mean(np.sum(relative**2, axis=1))))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluations", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--domain-filter", required=True)
    parser.add_argument("--relative-selection-cutoff", type=float, default=1.25)
    parser.add_argument("--score-power", type=float, default=2.0)
    args = parser.parse_args()
    if args.relative_selection_cutoff < 1 or args.score_power < 0:
        raise ValueError("invalid selection cutoff or score power")
    paths = [Path(path) for path in sorted(glob.glob(args.evaluations))]
    candidates = []
    for path in paths:
        report_path = path.with_name(path.name.replace(".evaluation.npz", ".json"))
        report = json.loads(report_path.read_text())
        if report.get("domain_filter") == args.domain_filter:
            candidates.append((path, report))
    if not candidates:
        raise ValueError("no evaluation archives satisfy the domain filter")
    scores = np.asarray([report["selection_rrmse"] for _, report in candidates])
    best = float(scores.min())
    admitted = scores <= best * args.relative_selection_cutoff
    selected = [candidate for candidate, keep in zip(candidates, admitted) if keep]
    scores = scores[admitted]
    weights = np.power(np.maximum(scores, 1e-8), -args.score_power)
    weights /= weights.sum()
    archives = [np.load(path) for path, _ in selected]
    reference = archives[0]
    for archive in archives[1:]:
        for key in ("indices", "truth", "point_id"):
            if not np.array_equal(archive[key], reference[key]):
                raise ValueError(f"evaluation archives disagree on {key}")
    transformed = np.asarray(
        [transform_targets(archive["prediction"]) for archive in archives]
    )
    ensemble = inverse_transform_targets(
        np.average(transformed, axis=0, weights=weights)
    )
    truth = reference["truth"]
    point_id = reference["point_id"]
    result = {
        "policy": "selection_gated_inverse_score_weighting",
        "domain_filter": args.domain_filter,
        "relative_selection_cutoff": args.relative_selection_cutoff,
        "inverse_score_power": args.score_power,
        "best_selection_rrmse": best,
        "candidate_models": len(candidates),
        "admitted_models": [
            {
                "evaluation": str(path),
                "selection_rrmse": report["selection_rrmse"],
                "untouched_evaluation_rrmse": report[
                    "untouched_evaluation_rrmse"
                ],
                "weight": float(weight),
            }
            for (path, report), weight in zip(selected, weights)
        ],
        "ensemble_untouched_evaluation_rrmse": _point_rrmse(
            ensemble, truth, point_id
        ),
        "ensemble_median_absolute_relative_error": np.median(
            np.abs((ensemble - truth) / truth), axis=0
        ).tolist(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
