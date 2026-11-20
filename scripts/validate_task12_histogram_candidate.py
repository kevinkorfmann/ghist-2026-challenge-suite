"""Generate fresh Task12 confirmation only after a fixed candidate is locked."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import time

import joblib
import numpy as np

from ghist_relatedness.validation import validate_simulation_splits
from ghist_submission.validate import validate_submission
from scripts.fit_human_background_corrected import observed_predictions
from scripts.fit_human_empirical_background import EmpiricalRegime, augment_empirical
from scripts.validate_human_background_correction import read_pedsim_with_lengths
from scripts.task12_background_correction import (
    FULL_HISTOGRAM,
    METADATA_KEY,
    correct_for_model,
)
from scripts.task12_histogram_protocol import (
    compare_runs,
    load_frozen_protocol,
    require_protocol_digest,
    score_predictions,
    sha256,
    validate_pairs,
    validate_protocol,
    write_json,
)


def validate_candidate_lock(candidate_dir, protocol_path, protocol, protocol_digest):
    candidate_dir = Path(candidate_dir)
    require_protocol_digest(protocol_path, protocol_digest)
    lock = json.loads((candidate_dir / "candidate_lock.json").read_text())
    checks = {
        "protocol_sha256": protocol_digest,
        "model_sha256": sha256(candidate_dir / "model.joblib"),
        "training_sha256": sha256(candidate_dir / "training.json"),
        "source_sha256": protocol["source_sha256"],
        "runtime_manifest_sha256": protocol["runtime_manifest_sha256"],
    }
    if any(lock.get(key) != value for key, value in checks.items()):
        raise ValueError("Task12 candidate lock differs from fitted artifact/protocol")
    return lock


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    snapshot = Path(__file__).resolve().parents[1]
    protocol, protocol_digest = load_frozen_protocol(args.protocol)
    validate_protocol(protocol, args.root, snapshot)
    lock = validate_candidate_lock(
        args.candidate_dir, args.protocol, protocol, protocol_digest
    )
    args.out_dir.mkdir(parents=True, exist_ok=False)
    settings, training = protocol["confirmation"], protocol["training"]
    started = time.monotonic()
    models = {
        item["name"]: joblib.load(args.root / item["path"])
        for item in settings["comparators"][:2]
    }
    models["full_histogram_v2"] = joblib.load(args.candidate_dir / "model.joblib")
    if (
        models["full_histogram_v2"].get(METADATA_KEY) != training["metadata"]
        or training["metadata"]["version"] != FULL_HISTOGRAM
    ):
        raise ValueError("Fitted candidate correction metadata differs")
    for name, model in models.items():
        if (
            model["background_profile"] != training["background_profile"]
            or model["corrected_model_weight"] != 1
        ):
            raise ValueError(f"Model profile/blend mismatch: {name}")
        if name != "full_histogram_v2" and METADATA_KEY in model:
            raise ValueError(
                "Historical comparator unexpectedly changed its correction metadata"
            )
    baseline = joblib.load(args.root / protocol["baseline_model"])
    calibration = [args.root / name for name in training["calibration"]]
    reference, _ = read_pedsim_with_lengths(calibration[0], training["genome_cm"])
    validate_pairs(reference)
    report = {
        "protocol_sha256": protocol_digest,
        "candidate_lock_sha256": sha256(args.candidate_dir / "candidate_lock.json"),
        "candidate_model_sha256": lock["model_sha256"],
        "runs": [],
        "generated": [],
        "component_only": True,
        "full_pipeline_validated": False,
        "candidate_promoted": False,
        "regimes_share_pedigrees": True,
        "public_answers_used": False,
        "submission_count": 0,
    }
    generated = []
    for index, seed in enumerate(settings["seeds"]):
        # No generation is permitted before the unchanged fit lock is verified.
        validate_candidate_lock(
            args.candidate_dir, args.protocol, protocol, protocol_digest
        )
        prefix = args.out_dir / f"holdout_{seed}"
        if list(args.out_dir.glob(prefix.name + ".*")):
            raise FileExistsError("Confirmation seed output already exists")
        command = [
            str(args.root / "external/ped-sim/ped-sim"),
            "-d",
            str(snapshot / settings["template"]),
            "-m",
            str(args.root / "external/ped-sim-map/refined_mf.simmap"),
            "--intf",
            str(args.root / "external/ped-sim/interfere/nu_p_campbell.tsv"),
            "--seed",
            str(seed),
            "-o",
            str(prefix),
        ]
        with prefix.with_suffix(".stdout").open("x") as handle:
            subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
        require_protocol_digest(args.protocol, protocol_digest)
        segment, log = prefix.with_suffix(".seg"), prefix.with_suffix(".log")
        found = re.findall(r"Random seed:\s*(\d+)", log.read_text())
        if found != [str(seed)]:
            raise ValueError("Generated simulation seed provenance differs")
        generated.append(segment)
        split_manifest = validate_simulation_splits(
            {"calibration": calibration, "confirmation": generated}
        )
        records, lengths = read_pedsim_with_lengths(segment, training["genome_cm"])
        expected_pairs = validate_pairs(records)
        report["generated"].append(
            {
                "seed": seed,
                "command": command,
                "segment_sha256": sha256(segment),
                "log_sha256": sha256(log),
                "expected_template_pairs": 7000,
                "segment_present_pairs": len(records),
                "zero_ibd_omitted_pairs": 7000 - len(records),
            }
        )
        for regime_index, regime_settings in enumerate(settings["regimes"]):
            augmentation_seed = (
                settings["augmentation_seed"] + 100 * index + regime_index
            )
            augmented = augment_empirical(
                records,
                lengths,
                training["background_profile"],
                EmpiricalRegime(**regime_settings),
                np.random.default_rng(augmentation_seed),
                training["genome_cm"],
            )
            if validate_pairs(augmented) != expected_pairs:
                raise ValueError("Confirmation augmentation changed target pairs")
            result = {
                "seed": seed,
                "regime_index": regime_index,
                "regime": regime_settings,
                "augmentation_seed": augmentation_seed,
                "models": {},
            }
            for name, model in models.items():
                corrected = correct_for_model(
                    augmented,
                    reference,
                    model,
                    minimum_cm=training["minimum_cm"],
                    genome_cm=training["genome_cm"],
                )
                predictions = observed_predictions(
                    model["corrected_payload"],
                    baseline["model"],
                    baseline["feature_names"],
                    corrected,
                    model["feature_names"],
                    model["corrected_model_weight"],
                )
                metrics = score_predictions(predictions, corrected)
                output = args.out_dir / f"seed{index}_regime{regime_index}_{name}.csv"
                with output.open("x", newline="") as handle:
                    writer = csv.DictWriter(
                        handle, fieldnames=("indID1", "indID2", "relation")
                    )
                    writer.writeheader()
                    writer.writerows(
                        {k: row[k] for k in writer.fieldnames} for row in predictions
                    )
                result["models"][name] = {
                    **metrics,
                    "csv_sha256": sha256(output),
                    "validation": validate_submission(
                        output, 12, "testing", expected_pairs
                    ),
                }
            write_json(args.out_dir / f"seed{index}_regime{regime_index}.json", result)
            report["runs"].append(result)
        print(
            json.dumps({"completed_seed": seed, "segment_pairs": len(records)}),
            flush=True,
        )
    report["split_manifest"] = split_manifest
    report["comparisons"] = {
        name: compare_runs(report["runs"], name, settings["gate"])
        for name in models
        if name != "full_histogram_v2"
    }
    metrics = (
        "accuracy",
        "balanced_accuracy",
        "meioses_accuracy",
        "template_accuracy_counting_unpredicted_as_wrong",
        "negative_log_likelihood",
        "multiclass_brier_sum",
        "ece_10_equal_width_bins",
    )
    report["aggregate"] = {
        name: {
            metric: float(np.mean([r["models"][name][metric] for r in report["runs"]]))
            for metric in metrics
        }
        for name in models
    }
    report["seconds"] = time.monotonic() - started
    validate_candidate_lock(
        args.candidate_dir, args.protocol, protocol, protocol_digest
    )
    validate_protocol(protocol, args.root, snapshot)
    require_protocol_digest(args.protocol, protocol_digest)
    for item, path in zip(report["generated"], generated, strict=True):
        if sha256(path) != item["segment_sha256"]:
            raise ValueError("Generated confirmation input changed during evaluation")
    require_protocol_digest(args.protocol, protocol_digest)
    write_json(args.out_dir / "confirmation.json", report)
    print(
        json.dumps(
            {
                "aggregate": report["aggregate"],
                "comparisons": report["comparisons"],
                "seconds": report["seconds"],
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
