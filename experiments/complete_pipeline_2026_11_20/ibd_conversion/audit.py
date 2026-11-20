"""Compare converters on preserved synthetic Task11 cohorts; never refit a model."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import joblib
import numpy as np


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def pair_summaries(module, ibd: Path, map_path: Path, threshold: float):
    segments, _, id_to_sample = module.convert_segments(ibd, map_path, threshold)
    summaries = defaultdict(lambda: {"ibd1_cm": 0.0, "ibd2_cm": 0.0, "segments": 0})
    for left, right, _, _, _, full, length in segments:
        pair = tuple(sorted((id_to_sample[left], id_to_sample[right])))
        summary = summaries[pair]
        summary["ibd2_cm" if full else "ibd1_cm"] += length
        summary["segments"] += 1
    return summaries


def main():
    # Imported historical modules are evidence, including their source folders.
    sys.dont_write_bytecode = True
    parser = argparse.ArgumentParser()
    parser.add_argument("--historical-root", type=Path, required=True)
    parser.add_argument("--corrected", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    source = args.historical_root / "source"
    sys.path[:0] = [str(source), str(source / "src"), str(source / "scripts")]
    baseline = load_module(
        "baseline_converter", source / "scripts/reconstruct_bonsai_pedigree.py"
    )
    corrected = load_module("corrected_converter", args.corrected)
    forest = load_module(
        "frozen_forest_features", source / "scripts/fit_task11_dog_relationships.py"
    )
    result = {
        "scope": "Read-only reanalysis of existing synthetic smoke; no independent validation or refitting",
        "submission_ready": False,
        "submission_count": 0,
        "job_id": os.environ.get("SLURM_JOB_ID"),
        "converter_sha256": {
            "baseline": digest(Path(baseline.__file__)),
            "corrected": digest(args.corrected),
        },
        "audit_sha256": digest(Path(__file__)),
        "cohorts": [],
    }
    rows = []
    feature_rows = []
    for seed in (2026112711, 2026112811):
        cohort = args.historical_root / "results" / f"seed_{seed}"
        protocol = json.loads((cohort / "protocol.json").read_text())
        truth = json.loads((cohort / "truth.json").read_text())
        assert (
            result["converter_sha256"]["baseline"]
            == protocol["sources"]["scripts/reconstruct_bonsai_pedigree.py"]
        )
        assert (
            digest(Path(forest.__file__))
            == protocol["sources"]["scripts/fit_task11_dog_relationships.py"]
        )
        map_path = Path(protocol["provenance"]["map"]["path"])
        model_path = Path(protocol["provenance"]["model"]["path"])
        assert digest(map_path) == protocol["provenance"]["map"]["sha256"]
        assert digest(model_path) == protocol["provenance"]["model"]["sha256"]
        watched = sorted(path for path in cohort.rglob("*") if path.is_file()) + [
            map_path,
            model_path,
            Path(baseline.__file__),
            Path(forest.__file__),
        ]
        before = {str(path): digest(path) for path in watched}
        sample_family = {
            sample: index
            for index, family in enumerate(truth["families"])
            for sample in family
        }
        cohort_summary = {"seed": seed, "conversion": [], "masked_features": []}
        for threshold in (0.0, 10.0):
            old = pair_summaries(
                baseline, cohort / "segments.ibd.gz", map_path, threshold
            )
            new = pair_summaries(
                corrected, cohort / "segments.ibd.gz", map_path, threshold
            )
            threshold_rows = []
            for pair in sorted(set(old) | set(new)):
                family = sample_family[pair[0]]
                stratum = (
                    truth["family_mating_types"][family]
                    if sample_family[pair[1]] == family
                    else "cross_family"
                )
                row = {
                    "seed": seed,
                    "threshold_cm": threshold,
                    "sample1": pair[0],
                    "sample2": pair[1],
                    "parent_mating": stratum,
                }
                for key in ("ibd1_cm", "ibd2_cm", "segments"):
                    row[f"old_{key}"] = old[pair][key]
                    row[f"new_{key}"] = new[pair][key]
                    row[f"delta_{key}"] = new[pair][key] - old[pair][key]
                threshold_rows.append(row)
            rows.extend(threshold_rows)
            by_stratum = {}
            for stratum in sorted({row["parent_mating"] for row in threshold_rows}):
                selected = [
                    row for row in threshold_rows if row["parent_mating"] == stratum
                ]
                by_stratum[stratum] = {
                    "pairs": len(selected),
                    "pairs_with_ibd2_change": sum(
                        abs(row["delta_ibd2_cm"]) > 1e-8 for row in selected
                    ),
                    **{
                        key: sum(row[key] for row in selected)
                        for key in (
                            "old_ibd1_cm",
                            "new_ibd1_cm",
                            "old_ibd2_cm",
                            "new_ibd2_cm",
                            "delta_ibd2_cm",
                        )
                    },
                }
            cohort_summary["conversion"].append(
                {
                    "threshold_cm": threshold,
                    "pairs": len(threshold_rows),
                    "pairs_with_ibd2_change": sum(
                        abs(row["delta_ibd2_cm"]) > 1e-8 for row in threshold_rows
                    ),
                    "old_ibd2_cm": sum(row["old_ibd2_cm"] for row in threshold_rows),
                    "new_ibd2_cm": sum(row["new_ibd2_cm"] for row in threshold_rows),
                    "max_pair_ibd2_decrease_cm": max(
                        -row["delta_ibd2_cm"] for row in threshold_rows
                    ),
                    "max_union_length_difference_cm": max(
                        abs(row["delta_ibd1_cm"] + row["delta_ibd2_cm"])
                        for row in threshold_rows
                    ),
                    "by_parent_mating": by_stratum,
                }
            )
        model = joblib.load(model_path)
        model.n_jobs = 2
        for hbd_threshold in (2, 5, 10):
            matrices = []
            for converter in (baseline, corrected):
                forest.convert_segments = converter.convert_segments
                forest.read_map = converter.read_map
                matrix, pairs, family_indices = forest.read_observed(
                    cohort / "segments.ibd.gz",
                    cohort / "segments.hbd.gz",
                    map_path,
                    cohort / "clusters.json",
                    hbd_threshold,
                )
                matrices.append(matrix)
            old_x, new_x = matrices
            previous = json.loads(
                (cohort / f"pair_predictions_hbd{hbd_threshold}.json").read_text()
            )
            previous_by_pair = {(row["indID1"], row["indID2"]): row for row in previous}
            preserved_x = np.asarray(
                [
                    [
                        previous_by_pair[pair]["features"][name]
                        for name in forest.feature_names()
                    ]
                    for pair in pairs
                ]
            )
            np.testing.assert_allclose(old_x, preserved_x, atol=1e-8, rtol=1e-12)
            old_prob, new_prob = model.predict_proba(old_x), model.predict_proba(new_x)
            preserved_prob = np.asarray(
                [
                    [
                        previous_by_pair[pair]["probabilities"][label]
                        for label in model.classes_
                    ]
                    for pair in pairs
                ]
            )
            np.testing.assert_allclose(old_prob, preserved_prob, atol=1e-12, rtol=1e-12)
            old_pred = model.classes_[old_prob.argmax(axis=1)]
            new_pred = model.classes_[new_prob.argmax(axis=1)]
            delta = new_x - old_x
            changed = ~np.isclose(old_x, new_x, atol=1e-8, rtol=1e-12)
            for index, pair in enumerate(pairs):
                feature_rows.append(
                    {
                        "seed": seed,
                        "hbd_threshold_cm": hbd_threshold,
                        "sample1": pair[0],
                        "sample2": pair[1],
                        "changed_features": int(changed[index].sum()),
                        "max_absolute_feature_change": float(
                            np.max(np.abs(delta[index]))
                        ),
                        "delta_total_ibd1_cm": float(delta[index, 2]),
                        "delta_total_ibd2_cm": float(delta[index, 3]),
                        "old_forest_label": str(old_pred[index]),
                        "new_forest_label": str(new_pred[index]),
                        "max_absolute_probability_change": float(
                            np.max(np.abs(new_prob[index] - old_prob[index]))
                        ),
                    }
                )
            cohort_summary["masked_features"].append(
                {
                    "hbd_threshold_cm": hbd_threshold,
                    "pairs": len(pairs),
                    "baseline_features_and_probabilities_reproduced": True,
                    "pairs_with_feature_change": int(np.any(changed, axis=1).sum()),
                    "changed_feature_names": [
                        name
                        for name, flag in zip(
                            forest.feature_names(), np.any(changed, axis=0), strict=True
                        )
                        if flag
                    ],
                    "max_absolute_feature_change": float(np.max(np.abs(delta))),
                    "forest_label_changes": int(np.sum(old_pred != new_pred)),
                    "old_correct_siblings": int(np.sum(old_pred == "siblings")),
                    "new_correct_siblings": int(np.sum(new_pred == "siblings")),
                    "max_absolute_probability_change": float(
                        np.max(np.abs(new_prob - old_prob))
                    ),
                }
            )
        after = {str(path): digest(path) for path in watched}
        assert before == after, (
            "Historical input or result changed during read-only audit"
        )
        cohort_summary["historical_files_unchanged"] = True
        cohort_summary["historical_sha256"] = before
        result["cohorts"].append(cohort_summary)
    for filename, data in (
        ("conversion_pairs.csv", rows),
        ("masked_feature_pairs.csv", feature_rows),
    ):
        with (args.out / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    (args.out / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "job_id": result["job_id"],
                "cohorts": len(result["cohorts"]),
                "submission_ready": False,
            }
        )
    )


if __name__ == "__main__":
    main()
