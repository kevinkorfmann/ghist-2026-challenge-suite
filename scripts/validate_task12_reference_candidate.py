"""Fixed Task 12 candidate versus historical empirical model on new pedigrees."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import joblib
import numpy as np

from ghist_relatedness.task9_release import sha256
from ghist_relatedness.validation import validate_simulation_splits
from ghist_submission.validate import validate_submission
from scripts.fit_human_background_corrected import correct_records, observed_predictions
from scripts.fit_human_empirical_background import EmpiricalRegime, augment_empirical
from scripts.fit_human_family_relationships import MEIOSES
from scripts.validate_human_background_correction import read_pedsim_with_lengths
from scripts.train_task12_reference_candidate import validate_protocol


def score(predictions, records):
    truth = {(row["indID1"], row["indID2"]): row["label"] for row in records}
    correct, degree, confidences = [], [], []
    by_class = {}
    brier, nll = [], []
    for row in predictions:
        label = truth[(row["indID1"], row["indID2"])]
        hit = row["relation"] == label
        correct.append(hit)
        degree.append(MEIOSES[row["relation"]] == MEIOSES[label])
        confidences.append(row["confidence"])
        by_class.setdefault(label, []).append(hit)
        probabilities = row["probabilities"]
        values = np.asarray(list(probabilities.values()), dtype=float)
        if (not np.isfinite(values).all() or (values < 0).any()
                or not np.isclose(values.sum(), 1.0) or label not in probabilities):
            raise ValueError("Invalid relationship probabilities")
        brier.append(sum((value - (name == label)) ** 2 for name, value in probabilities.items()))
        nll.append(-np.log(max(probabilities[label], 1e-15)))
    correct, confidences = np.asarray(correct), np.asarray(confidences)
    ece = 0.0
    for index in range(10):
        selected = (confidences >= index / 10) & (confidences < (index + 1) / 10 if index < 9 else confidences <= 1)
        if selected.any():
            ece += selected.mean() * abs(correct[selected].mean() - confidences[selected].mean())
    return {"pairs": len(records), "accuracy": float(correct.mean()),
            "balanced_accuracy": float(np.mean([np.mean(v) for v in by_class.values()])),
            "meioses_accuracy": float(np.mean(degree)), "ece_10_equal_width_bins": float(ece),
            "multiclass_brier_sum": float(np.mean(brier)), "negative_log_likelihood": float(np.mean(nll)),
            "recall_by_label": {k: float(np.mean(v)) for k, v in by_class.items()},
            "confidence_is_uncalibrated": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--holdouts", type=Path, nargs=4, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    validate_protocol(protocol, args.root)
    expected_names = [f"holdout_{seed}.seg" for seed in protocol["holdout_seeds"]]
    if [path.name for path in args.holdouts] != expected_names:
        raise ValueError("Holdout inputs do not match the frozen seed order")
    training = json.loads((args.candidate_dir / "training.json").read_text())
    candidate_path = args.candidate_dir / "model.joblib"
    if training["protocol_sha256"] != sha256(args.protocol) or training["model_sha256"] != sha256(candidate_path):
        raise ValueError("Candidate or protocol changed after fitting")
    # The reference intentionally belongs to calibration, never to holdout.
    calibration = [args.root / name for name in protocol["calibration"]]
    split_manifest = validate_simulation_splits({"calibration": calibration, "holdout": args.holdouts})
    reference, _ = read_pedsim_with_lengths(calibration[0], protocol["genome_cm"])
    models = {"historical": joblib.load(args.root / protocol["historical_model"]),
              "candidate": joblib.load(candidate_path)}
    baseline = joblib.load(args.root / protocol["baseline_model"])
    profile = models["historical"]["background_profile"]
    if profile != models["candidate"]["background_profile"]:
        raise ValueError("Profile changed between fixed comparison models")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    report = {"component_only": True, "protocol_sha256": sha256(args.protocol),
              "model_sha256": {"historical": sha256(args.root / protocol["historical_model"]), "candidate": sha256(candidate_path)},
              "split_manifest": split_manifest, "runs": [], "regimes_share_pedigrees": True,
              "independent_correction_reference": True, "full_vcf_pipeline_validated": False,
              "public_answers_used": False, "candidate_promoted": False, "submission_count": 0}
    for index, path in enumerate(args.holdouts):
        records, lengths = read_pedsim_with_lengths(path, protocol["genome_cm"])
        for regime_index, settings in enumerate(protocol["evaluation_regimes"]):
            regime = EmpiricalRegime(**settings)
            augmented = augment_empirical(records, lengths, profile, regime,
                np.random.default_rng(protocol["augmentation_seed"] + 100 * index + regime_index), protocol["genome_cm"])
            corrected = correct_records(augmented, reference, protocol["minimum_cm"], protocol["genome_cm"])
            result = {"seed": protocol["holdout_seeds"][index], "regime": settings, "models": {}}
            expected_pairs = {(row["indID1"], row["indID2"]) for row in records}
            for name, model in models.items():
                predictions = observed_predictions(model["corrected_payload"], baseline["model"],
                    baseline["feature_names"], corrected, model["feature_names"], model["corrected_model_weight"])
                output = args.out_dir / f"seed{index}_regime{regime_index}_{name}.csv"
                with output.open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=("indID1", "indID2", "relation"))
                    writer.writeheader()
                    writer.writerows({k: row[k] for k in writer.fieldnames} for row in predictions)
                result["models"][name] = {**score(predictions, corrected),
                    "validation": validate_submission(output, 12, "testing", expected_pairs),
                    "output_sha256": sha256(output)}
            report["runs"].append(result)
    gains = [row["models"]["candidate"]["accuracy"] - row["models"]["historical"]["accuracy"] for row in report["runs"]]
    degree_gains = [row["models"]["candidate"]["meioses_accuracy"] - row["models"]["historical"]["meioses_accuracy"] for row in report["runs"]]
    per_seed = [float(np.mean([gain for gain, row in zip(gains, report["runs"], strict=True) if row["seed"] == seed])) for seed in protocol["holdout_seeds"]]
    report["mean_exact_gain_by_seed"] = per_seed
    report["mean_exact_gain"] = float(np.mean(gains))
    report["mean_meioses_gain"] = float(np.mean(degree_gains))
    report["gate"] = {"mean_exact_gain_at_least_0.02": np.mean(gains) >= .02,
                      "positive_mean_gain_in_at_least_80pct_seeds": np.mean(np.asarray(per_seed) > 0) >= .80,
                      "worst_exact_regression_at_most_0.01": min(gains) >= -.01,
                      "mean_meioses_regression_at_most_0.005": np.mean(degree_gains) >= -.005}
    report["gate"] = {key: bool(value) for key, value in report["gate"].items()}
    report["gate"]["passed"] = all(report["gate"].values())
    (args.out_dir / "holdout.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "runs"}, indent=2))


if __name__ == "__main__":
    main()
