from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
from ghist_secondary_contact.contract import require_spectrum_contract
import torch

from ghist_secondary_contact.features import inverse_transform_targets
from ghist_secondary_contact.model import SecondaryContactNet
from ghist_secondary_contact.spectrum import transform_spectrum_counts
from ghist_secondary_contact.vcf_fast import scan_secondary_contact_vcf


OUTPUT_NAMES = (
    "mainland_population_size",
    "island_population_size",
    "generations_since_split",
    "generations_since_migration",
    "migration_rate",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--checkpoints", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--domain-filter", type=str)
    parser.add_argument("--relative-selection-cutoff", type=float, default=1.25)
    parser.add_argument("--score-power", type=float, default=2.0)
    args = parser.parse_args()
    checkpoints = [Path(path) for path in sorted(glob.glob(args.checkpoints))]
    if not checkpoints:
        raise FileNotFoundError(f"no checkpoints match {args.checkpoints!r}")
    if args.relative_selection_cutoff < 1 or args.score_power < 0:
        raise ValueError("invalid selection cutoff or score power")
    candidates = []
    for path in checkpoints:
        report_path = path.with_suffix(".json")
        report = json.loads(report_path.read_text()) if report_path.exists() else {}
        if (
            args.domain_filter is not None
            and report.get("domain_filter") != args.domain_filter
        ):
            continue
        candidates.append((path, report))
    if not candidates:
        raise ValueError("no checkpoints satisfy the requested domain filter")
    selection_scores = np.asarray(
        [
            float(
                report.get(
                    "selection_rrmse",
                    torch.load(path, map_location="cpu", weights_only=False)[
                        "selection_rrmse"
                    ],
                )
            )
            for path, report in candidates
        ]
    )
    best_selection = float(selection_scores.min())
    admitted = selection_scores <= best_selection * args.relative_selection_cutoff
    candidates = [candidate for candidate, keep in zip(candidates, admitted) if keep]
    selection_scores = selection_scores[admitted]
    weights = np.power(np.maximum(selection_scores, 1e-8), -args.score_power)
    weights /= weights.sum()
    spectrum, metadata = scan_secondary_contact_vcf(args.vcf)
    predictions = []
    uncertainty = []
    records = []
    for (path, report), weight in zip(candidates, weights):
        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
        require_spectrum_contract(checkpoint)
        spectrum_transform = checkpoint.get("spectrum_transform", "log_l2")
        transformed_spectrum = transform_spectrum_counts(spectrum, spectrum_transform)
        configuration = checkpoint["model_configuration"]
        model = SecondaryContactNet(**configuration)
        model.load_state_dict(checkpoint["model"])
        model.eval()
        normalization = checkpoint["normalization"]
        window_mean = np.asarray(normalization["window_mean"])
        window_scale = np.asarray(normalization["window_scale"])
        global_mean = np.asarray(normalization["global_mean"])
        global_scale = np.asarray(normalization["global_scale"])
        windows = (np.zeros((1, 1, 1)) - window_mean) / window_scale
        global_summary = (
            np.asarray([[np.log1p(metadata["records"])]]) - global_mean
        ) / global_scale
        with torch.no_grad():
            mean, log_variance = model(
                torch.as_tensor(windows, dtype=torch.float32),
                torch.as_tensor(transformed_spectrum[None, None], dtype=torch.float32),
                torch.as_tensor(global_summary, dtype=torch.float32),
            )
        target_mean = np.asarray(normalization["target_mean"])
        target_scale = np.asarray(normalization["target_scale"])
        transformed = mean.numpy()[0] * target_scale + target_mean
        predictions.append(inverse_transform_targets(transformed))
        uncertainty.append(np.sqrt(np.exp(log_variance.numpy()[0])) * target_scale)
        records.append(
            {
                "checkpoint": str(path),
                "best_epoch": checkpoint["epoch"],
                "selection_rrmse": checkpoint["selection_rrmse"],
                "domain_filter": report.get("domain_filter"),
                "spectrum_transform": spectrum_transform,
                "ensemble_weight": float(weight),
                "prediction": dict(zip(OUTPUT_NAMES, predictions[-1].tolist())),
                "transformed_standard_deviation": uncertainty[-1].tolist(),
            }
        )
    transformed_predictions = np.asarray(
        [
            np.asarray(
                (
                    np.log(value[0]),
                    np.log(value[1]),
                    np.log(value[2]),
                    np.log(value[3] / value[2]) - np.log1p(-value[3] / value[2]),
                    np.log(value[4]),
                )
            )
            for value in predictions
        ]
    )
    ensemble = inverse_transform_targets(
        np.average(transformed_predictions, axis=0, weights=weights)
    )
    result = {
        "policy": "selection_gated_local_prediction_only",
        "selection_rule": {
            "domain_filter": args.domain_filter,
            "relative_cutoff": args.relative_selection_cutoff,
            "inverse_score_power": args.score_power,
            "best_selection_rrmse": best_selection,
            "candidate_checkpoints": len(checkpoints),
            "admitted_checkpoints": len(candidates),
        },
        "observed_records": metadata["records"],
        "models": records,
        "transformed_coordinate_ensemble": dict(zip(OUTPUT_NAMES, ensemble.tolist())),
        "model_spread_transformed": np.sqrt(
            np.average(
                (
                    transformed_predictions
                    - np.average(transformed_predictions, axis=0, weights=weights)
                )
                ** 2,
                axis=0,
                weights=weights,
            )
        ).tolist(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
