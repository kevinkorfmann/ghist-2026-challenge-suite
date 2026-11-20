"""One fixed Task12 histogram-v2 fit on original calibration simulations only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import joblib
import numpy as np
from sklearn.base import clone
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.utils.class_weight import compute_sample_weight

from ghist_relatedness.validation import validate_simulation_splits
from scripts.fit_human_empirical_background import EmpiricalRegime, augment_empirical
from scripts.fit_human_background_corrected import corrected_feature_names
from scripts.fit_human_family_relationships import LABELS, matrix
from scripts.validate_human_background_correction import read_pedsim_with_lengths
from scripts.task12_background_correction import METADATA_KEY, correct_records_v2
from scripts.task12_histogram_protocol import (
    load_frozen_protocol,
    require_protocol_digest,
    sha256,
    validate_pairs,
    validate_protocol,
    write_json,
)


def build_training_records(paths, settings):
    rng = np.random.default_rng(settings["training_seed"])
    loaded = [read_pedsim_with_lengths(path, settings["genome_cm"]) for path in paths]
    output, regimes = [], []
    for index, (clean, lengths) in enumerate(loaded):
        validate_pairs(clean)
        reference = [
            record
            for other, (rows, _) in enumerate(loaded)
            if other != index
            for record in rows
        ]
        output.extend(clean)
        for _ in range(settings["augmentations"]):
            regime = EmpiricalRegime(
                float(rng.uniform(0.25, 1.75)),
                float(rng.uniform(0.1, 0.7)),
                float(rng.uniform(0.6, 1.5)),
            )
            augmented = augment_empirical(
                clean,
                lengths,
                settings["background_profile"],
                regime,
                rng,
                settings["genome_cm"],
            )
            corrected = correct_records_v2(
                augmented,
                reference,
                minimum_cm=settings["minimum_cm"],
                genome_cm=settings["genome_cm"],
            )
            if validate_pairs(corrected) != validate_pairs(clean):
                raise ValueError("Training augmentation changed target pairs")
            output.extend(corrected)
            regimes.append(vars(regime))
    return output, regimes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    snapshot = Path(__file__).resolve().parents[1]
    protocol, protocol_digest = load_frozen_protocol(args.protocol)
    validate_protocol(protocol, args.root, snapshot)
    require_protocol_digest(args.protocol, protocol_digest)
    settings = protocol["training"]
    if args.out_dir.resolve() != (snapshot / settings["artifact_directory"]).resolve():
        raise ValueError("Candidate destination differs from frozen protocol")
    args.out_dir.mkdir(parents=True, exist_ok=False)
    historical = joblib.load(
        args.root / protocol["confirmation"]["comparators"][0]["path"]
    )
    if (
        historical["corrected_model_weight"] != 1
        or historical["feature_names"] != corrected_feature_names(5)
        or historical["background_profile"] != settings["background_profile"]
    ):
        raise ValueError("Historical model feature/profile/blend contract differs")
    previous = historical["corrected_payload"]
    if previous["kind"] != "direct" or not isinstance(
        previous["model"], HistGradientBoostingClassifier
    ):
        raise ValueError("Expected fixed historical HGB architecture")
    estimator = clone(previous["model"]).set_params(early_stopping=False)
    if estimator.get_params() != protocol["fixed_hyperparameters"]:
        raise ValueError("Fixed Task12 HGB hyperparameters differ")
    paths = [args.root / name for name in settings["calibration"]]
    splits = validate_simulation_splits({"calibration": paths})
    started = time.monotonic()
    records, regimes = build_training_records(paths, settings)
    if (
        len(records) != protocol["expected_training_rows"]
        or regimes != protocol["training_regimes_reference"]
    ):
        raise ValueError(
            "Calibration augmentation differs from the fixed original recipe"
        )
    x = matrix(records, historical["feature_names"])
    y = np.asarray([row["label"] for row in records])
    if not np.isfinite(x).all() or set(y) != set(LABELS):
        raise ValueError("Invalid training features or missing Task12 labels")
    require_protocol_digest(args.protocol, protocol_digest)
    estimator.fit(x, y, sample_weight=compute_sample_weight("balanced", y))
    require_protocol_digest(args.protocol, protocol_digest)
    candidate = {
        "feature_names": list(historical["feature_names"]),
        "background_profile": list(historical["background_profile"]),
        "corrected_model_weight": 1.0,
        "corrected_payload": {"kind": "direct", "model": estimator},
        METADATA_KEY: settings["metadata"],
        "training_provenance": {
            "status": "CANDIDATE_NOT_PROMOTED",
            "protocol_sha256": protocol_digest,
            "calibration_splits": splits,
            "training_seed": settings["training_seed"],
            "historical_model_sha256": protocol["confirmation"]["comparators"][0][
                "sha256"
            ],
            "inherited_settings": [
                "feature_names",
                "background_profile",
                "fixed architecture and hyperparameters",
            ],
            "historical_metadata_keys_not_reused": sorted(
                set(historical)
                - {
                    "feature_names",
                    "background_profile",
                    "corrected_model_weight",
                    "corrected_payload",
                }
            ),
            "historical_selection_or_validation_scores_reused": False,
        },
    }
    model_path = args.out_dir / "model.joblib"
    joblib.dump(candidate, model_path)
    validate_protocol(protocol, args.root, snapshot)
    require_protocol_digest(args.protocol, protocol_digest)
    report = {
        "status": "CANDIDATE_NOT_PROMOTED",
        "protocol_sha256": protocol_digest,
        "model_sha256": sha256(model_path),
        "runtime_manifest_sha256": protocol["runtime_manifest_sha256"],
        "source_splits": splits,
        "training_pairs": len(records),
        "training_regimes": regimes,
        "hyperparameters": estimator.get_params(),
        METADATA_KEY: settings["metadata"],
        "seconds": time.monotonic() - started,
        "component_only": True,
        "holdouts_generated": False,
        "public_answers_used": False,
        "submission_count": 0,
    }
    require_protocol_digest(args.protocol, protocol_digest)
    write_json(args.out_dir / "training.json", report)
    require_protocol_digest(args.protocol, protocol_digest)
    write_json(
        args.out_dir / "candidate_lock.json",
        {
            "protocol_sha256": protocol_digest,
            "model_sha256": sha256(model_path),
            "training_sha256": sha256(args.out_dir / "training.json"),
            "source_sha256": protocol["source_sha256"],
            "runtime_manifest_sha256": protocol["runtime_manifest_sha256"],
            "candidate_promoted": False,
        },
    )
    print(
        json.dumps(
            {
                k: report[k]
                for k in ("status", "model_sha256", "training_pairs", "seconds")
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
