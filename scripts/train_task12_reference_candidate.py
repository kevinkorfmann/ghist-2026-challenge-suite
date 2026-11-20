"""Refit the fixed historical Task 12 architecture with independent references."""
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

from ghist_relatedness.task9_release import sha256
from ghist_relatedness.validation import validate_simulation_splits
from scripts.fit_human_empirical_background import training_set
from scripts.fit_human_family_relationships import matrix


def validate_protocol(protocol, root):
    snapshot = Path(__file__).resolve().parents[1]
    for name, digest in protocol["source_sha256"].items():
        if sha256(snapshot / name) != digest:
            raise ValueError(f"Locked Task12 source changed: {name}")
    for role in ("historical_model", "baseline_model"):
        if sha256(root / protocol[role]) != protocol[f"{role}_sha256"]:
            raise ValueError(f"Locked Task12 model changed: {role}")
    for name, digest in protocol["input_sha256"].items():
        if sha256(root / name) != digest:
            raise ValueError(f"Locked Task12 input changed: {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    validate_protocol(protocol, args.root)
    historical_path = args.root / protocol["historical_model"]
    historical = joblib.load(historical_path)
    previous = historical["corrected_payload"]
    if previous["kind"] != "direct" or not isinstance(previous["model"], HistGradientBoostingClassifier):
        raise ValueError("Expected the fixed historical direct histogram-gradient model")
    if historical["corrected_model_weight"] != 1.0:
        raise ValueError("Protocol assumes the recorded historical blend weight 1.0")
    paths = [args.root / name for name in protocol["calibration"]]
    splits = validate_simulation_splits({"calibration": paths})
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if (args.out_dir / "model.joblib").exists():
        raise FileExistsError("Use a separate versioned candidate directory")
    started = time.monotonic()
    records, regimes = training_set(
        paths, historical["background_profile"], protocol["augmentations"],
        protocol["training_seed"], protocol["minimum_cm"], protocol["genome_cm"],
    )
    x = matrix(records, historical["feature_names"])
    y = np.asarray([row["label"] for row in records])
    estimator = clone(previous["model"]).set_params(early_stopping=False)
    if estimator.max_iter != protocol["max_iter"]:
        raise ValueError("Historical architecture no longer matches fixed protocol")
    estimator.fit(x, y, sample_weight=compute_sample_weight("balanced", y))
    candidate = {**historical, "corrected_payload": {"kind": "direct", "model": estimator}}
    destination = args.out_dir / "model.joblib"
    joblib.dump(candidate, destination)
    report = {
        "status": "CANDIDATE_NOT_PROMOTED", "component_only": True,
        "protocol_sha256": sha256(args.protocol), "historical_model_sha256": sha256(historical_path),
        "model_sha256": sha256(destination), "source_splits": splits,
        "background_profile": historical["background_profile"],
        "profile_provenance": "Reuse the frozen unlabeled empirical bin profile; no challenge inputs read",
        "hyperparameters": estimator.get_params(), "training_pairs": len(records),
        "training_regimes": regimes, "seconds": time.monotonic() - started,
        "independent_correction_reference": True,
        "correction_reference_policy": "All other calibration replicates, leaving out the augmented replicate",
        "architecture_or_weight_search": False, "public_answers_used": False,
        "final_data_runs": False, "submission_count": 0,
    }
    (args.out_dir / "training.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
