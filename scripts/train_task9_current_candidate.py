"""Refit the historical forest on the current label contract without selection."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier

from ghist_relatedness.ibd import matrix, read_sex_map
from ghist_relatedness.labels import IBD_LABEL_MEIOSES
from ghist_relatedness.validation import validate_simulation_splits
from scripts.fit_task9_ibd_classifier import load_simulations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-model", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, nargs="+", required=True)
    parser.add_argument("--sex-map", type=Path, required=True)
    parser.add_argument("--genome-cm", type=float, default=3350.933402)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    provenance = validate_simulation_splits({"calibration": args.calibration})
    baseline = joblib.load(args.baseline_model)
    if not isinstance(baseline["model"], RandomForestClassifier):
        raise ValueError("Expected the frozen Task9 random forest")
    records = load_simulations(
        args.calibration, args.genome_cm, read_sex_map(args.sex_map)
    )
    x, features = matrix(records, baseline["features"])
    labels = np.asarray([row["label"] for row in records])
    model = clone(baseline["model"]).set_params(n_jobs=args.workers)
    model.fit(x, labels)
    if set(model.classes_) != set(IBD_LABEL_MEIOSES):
        raise ValueError("Calibration does not cover every current Task9 class")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump({"model": model, "features": features}, args.out_dir / "model.joblib")
    with args.baseline_model.open("rb") as handle:
        baseline_hash = hashlib.file_digest(handle, "sha256").hexdigest()
    report = {
        "status": "CANDIDATE_NOT_FROZEN",
        "selection": "No search; preserve frozen forest hyperparameters and remove retired training labels",
        "training_pairs": len(records),
        "labels": list(model.classes_),
        "source_splits": provenance,
        "baseline_sha256": baseline_hash,
        "public_answers_used": False,
        "new_evaluation_used_for_selection": False,
        "hyperparameters": model.get_params(),
        "submission_count": 0,
    }
    (args.out_dir / "training.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
