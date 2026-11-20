from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np

from ghist_secondary_contact.features import inverse_transform_targets
from ghist_secondary_contact.vcf_fast import scan_secondary_contact_vcf
from ghist_secondary_contact.contract import require_spectrum_contract
from ghist_secondary_contact.spectrum import transform_spectrum_counts


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
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    checkpoint = joblib.load(args.model)
    require_spectrum_contract(checkpoint)
    spectrum, metadata = scan_secondary_contact_vcf(args.vcf)
    transformed_spectrum = transform_spectrum_counts(
        spectrum, checkpoint["spectrum_transform"]
    )
    features = np.concatenate(
        (transformed_spectrum.ravel(), [np.log1p(metadata["records"])]), axis=0
    )[None].astype(np.float32)
    reduced = checkpoint["pca"].transform(features)
    transformed = (
        checkpoint["model"].predict(reduced)[0] * checkpoint["target_scale"]
        + checkpoint["target_mean"]
    )
    prediction = inverse_transform_targets(transformed)
    result = {
        "policy": "local_prediction_only",
        "model": str(args.model),
        "selection_report": checkpoint["report"],
        "observed_records": metadata["records"],
        "prediction": dict(zip(OUTPUT_NAMES, prediction.tolist())),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
