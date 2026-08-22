from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import torch

from ghist_secondary_contact.features import inverse_transform_targets
from ghist_secondary_contact.model import SecondaryContactNet
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
    args = parser.parse_args()
    checkpoints = [Path(path) for path in sorted(glob.glob(args.checkpoints))]
    if not checkpoints:
        raise FileNotFoundError(f"no checkpoints match {args.checkpoints!r}")
    spectrum, metadata = scan_secondary_contact_vcf(args.vcf)
    transformed_spectrum = np.log1p(spectrum.astype(np.float64))
    transformed_spectrum /= np.linalg.norm(transformed_spectrum) + 1e-12
    predictions = []
    uncertainty = []
    records = []
    for path in checkpoints:
        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
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
    ensemble = inverse_transform_targets(transformed_predictions.mean(axis=0))
    result = {
        "policy": "local_prediction_only",
        "observed_records": metadata["records"],
        "models": records,
        "transformed_coordinate_ensemble": dict(zip(OUTPUT_NAMES, ensemble.tolist())),
        "model_spread_transformed": transformed_predictions.std(axis=0).tolist(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
