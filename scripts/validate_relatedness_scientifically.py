"""Score frozen relatedness estimators on independently generated Ped-sim data.

This is an exact-IBD component audit, not VCF, detector or pedigree validation.
No fitting or model selection is performed. Confidence intervals resample whole
simulated pedigrees, keeping all pairs from each pedigree together.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path

import joblib
import numpy as np

from ghist_relatedness.ibd import (
    pair_features,
    read_segment_rows,
    read_sex_map,
    matrix as ibd_matrix,
)
from ghist_relatedness.labels import IBD_LABEL_MEIOSES
from scripts.predict_task9_current import accepted_probabilities
from scripts.fit_human_family_relationships import LABELS, MEIOSES, matrix, read_pedsim
from scripts.fit_human_background_corrected import align_probabilities, correct_records
from scripts.fit_human_background_robust import predict_probabilities
from scripts.fit_human_empirical_background import augment_empirical, EmpiricalRegime
from scripts.validate_human_background_correction import read_pedsim_with_lengths
from scripts.fit_task11_dog_relationships import read_pedsim as read_dog
from scripts.consensus_task11_dog_pedigrees import MEIOSES as DOG_MEIOSES


def source_hash(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def pedigree(sample: str) -> str:
    return re.split(r"_g\d", sample, maxsplit=1)[0]


def probability_metrics(classes, probability, truth, groups, meioses, seed=20261120):
    classes, truth, groups = np.asarray(classes), np.asarray(truth), np.asarray(groups)
    probability = np.asarray(probability, dtype=float)
    if probability.shape != (len(truth), len(classes)) or not np.all(
        np.isfinite(probability)
    ):
        raise ValueError("Invalid probability matrix")
    if np.any(probability < 0) or not np.allclose(probability.sum(axis=1), 1):
        raise ValueError("Probabilities must sum to one")
    lookup = {str(label): index for index, label in enumerate(classes)}
    indices = np.asarray([lookup[str(label)] for label in truth])
    prediction = classes[np.argmax(probability, axis=1)]
    correct = prediction == truth
    confidence = probability.max(axis=1)
    bins = []
    for index in range(10):
        selected = (confidence >= index / 10) & (
            (confidence < (index + 1) / 10) if index < 9 else (confidence <= 1)
        )
        if selected.any():
            bins.append(
                {
                    "lower": index / 10,
                    "upper": (index + 1) / 10,
                    "pairs": int(selected.sum()),
                    "mean_confidence": float(confidence[selected].mean()),
                    "accuracy": float(correct[selected].mean()),
                }
            )
    _, inverse = np.unique(groups, return_inverse=True)
    sizes = np.bincount(inverse)
    successes = np.bincount(inverse, weights=correct)
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(1000):
        chosen = rng.integers(0, len(sizes), size=len(sizes))
        estimates.append(float(successes[chosen].sum() / sizes[chosen].sum()))
    target = np.eye(len(classes))[indices]
    class_recall = {
        str(label): float(correct[truth == label].mean()) for label in np.unique(truth)
    }
    return {
        "pairs": len(truth),
        "independent_pedigrees": len(sizes),
        "accuracy": float(correct.mean()),
        "balanced_accuracy": float(np.mean(list(class_recall.values()))),
        "meioses_accuracy": float(
            np.mean(
                [
                    meioses[str(a)] == meioses[str(b)]
                    for a, b in zip(prediction, truth, strict=True)
                ]
            )
        ),
        "pedigree_bootstrap_accuracy_95pct": np.quantile(
            estimates, [0.025, 0.975]
        ).tolist(),
        "multiclass_brier_sum": float(
            np.mean(np.sum((probability - target) ** 2, axis=1))
        ),
        "negative_log_likelihood": float(
            -np.log(
                np.maximum(probability[np.arange(len(truth)), indices], 1e-15)
            ).mean()
        ),
        "ece_10_equal_width_bins": float(
            sum(
                row["pairs"] * abs(row["mean_confidence"] - row["accuracy"])
                for row in bins
            )
            / len(truth)
        ),
        "reliability_bins": bins,
        "recall_by_label": class_recall,
        "confidence_is_uncalibrated": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--segments", type=Path, required=True)
    parser.add_argument("--kind", choices=("task9", "human", "dog"), required=True)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--task9-model", type=Path)
    args = parser.parse_args()
    root = args.root
    result = {
        "component_only": True,
        "kind": args.kind,
        "scenario": args.scenario,
        "seed": args.seed,
        "segment_sha256": source_hash(args.segments),
        "public_answers_used": False,
        "models_refitted": False,
        "runs": {},
        "limitations": [
            "Exact segments bypass genotype error, phasing and IBD detection.",
            "Pairs with no IBD are absent from Ped-sim segment files.",
            "Pedigree templates match calibration; this is seed and crossover-model generalization, not novel topology validation.",
        ],
    }
    if args.kind == "task9":
        sex_map = read_sex_map(root / "external/ped-sim-map/refined_mf.simmap")
        model_path = args.task9_model or (
            root / "artifacts/task9_ibd_segments/model_v3_shape/task9_model.joblib"
        )
        payload = joblib.load(model_path)
        payload["model"].n_jobs = 2
        original = read_segment_rows(args.segments)
        expected = {
            "siblings": 0.5,
            "grandparent_grandchild": 0.25,
            "avuncular": 0.25,
            "double_cousin": 0.25,
            "full_1st_cousin": 0.125,
            "half_1st_cousin": 0.0625,
            "full_2nd_cousin": 0.03125,
            "full_3rd_cousin": 0.0078125,
            "full_4th_cousin": 0.001953125,
            "full_5th_cousin": 0.00048828125,
        }
        for tag, scale in [
            ("unchanged", 1.0),
            ("map_length_minus10pct", 0.9),
            ("map_length_plus10pct", 1.1),
        ]:
            rows = copy.deepcopy(original)
            if scale != 1:
                for row in rows:
                    for key in ("length", "start_cM", "end_cM"):
                        row[key] = str(float(row[key]) * scale)
                used_map = {
                    k: (v[0], v[1] * scale, v[2] * scale) for k, v in sex_map.items()
                }
            else:
                used_map = sex_map
            records = pair_features(rows, 3350.933402, used_map)
            records = [
                r
                for r in records
                if pedigree(str(r["indID1"])).rstrip("0123456789") in IBD_LABEL_MEIOSES
            ]
            truth = [pedigree(str(r["indID1"])).rstrip("0123456789") for r in records]
            x, _ = ibd_matrix(records, payload["features"])
            raw_probability = payload["model"].predict_proba(x)
            accepted = np.asarray(
                [str(label) in IBD_LABEL_MEIOSES for label in payload["model"].classes_]
            )
            scorable = raw_probability[:, accepted].sum(axis=1) > 0
            rejected = int(np.sum(~scorable))
            original_count = len(records)
            records = [r for r, keep in zip(records, scorable, strict=True) if keep]
            truth = [t for t, keep in zip(truth, scorable, strict=True) if keep]
            classes, prob = accepted_probabilities(
                payload["model"].classes_, raw_probability[scorable]
            )
            metric = probability_metrics(
                classes,
                prob,
                truth,
                [pedigree(str(r["indID1"])) for r in records],
                IBD_LABEL_MEIOSES,
            )
            metric["zero_accepted_probability_pairs"] = rejected
            metric["eligible_pairs_before_probability_failure"] = original_count
            metric["strict_production_predictor_would_abort"] = rejected > 0
            metric["accuracy_counting_probability_failures_wrong"] = (
                metric["accuracy"] * len(records) / original_count
            )
            levels = np.unique(list(expected.values()))
            estimated = np.asarray(
                [
                    levels[
                        np.argmin(
                            np.abs(
                                np.log2(levels)
                                - np.log2(max(float(r["relatedness"]), 1e-12))
                            )
                        )
                    ]
                    for r in records
                ]
            )
            metric["mendelian_expected_sharing_class_accuracy"] = float(
                np.mean(estimated == np.asarray([expected[t] for t in truth]))
            )
            metric["classifier_expected_sharing_class_accuracy"] = float(
                np.mean(
                    [
                        expected[str(c)] == expected[t]
                        for c, t in zip(
                            classes[np.argmax(prob, axis=1)], truth, strict=True
                        )
                    ]
                )
            )
            result["runs"][tag] = metric
        result["expected_pairs_before_zero_ibd_censoring"] = 11000
        result["model_sha256"] = source_hash(model_path)
    elif args.kind == "dog":
        path = (
            root
            / "artifacts/task11_multifamily_inbreeding/pair_model_hbd5/model.joblib"
        )
        model = joblib.load(path)
        model.n_jobs = 2
        x, y = read_dog(args.segments)
        # read_dog sorts pairs; recover that exact grouping for the bootstrap.
        rows = read_segment_rows(args.segments)
        pairs = sorted({tuple(sorted((r["indID1"], r["indID2"]))) for r in rows})
        result["runs"]["unchanged"] = probability_metrics(
            model.classes_,
            model.predict_proba(x),
            y,
            [pedigree(p[0]) for p in pairs],
            DOG_MEIOSES,
        )
        result["expected_pairs_before_zero_ibd_censoring"] = 4000
        result["model_sha256"] = source_hash(path)
        result["limitations"].append(
            "Dog founders are unrelated; observed inbreeding and HBD-mask fragmentation are not simulated."
        )
    else:
        baseline_path = (
            root / "artifacts/task10_large_family/pair_model_min5/model.joblib"
        )
        clean_path = (
            root / "artifacts/task10_large_family/clean_spectrum_v1/model.joblib"
        )
        empirical_path = (
            root
            / "artifacts/task12_large_family_background/empirical_background_v1/model.joblib"
        )
        baseline, clean, empirical = [
            joblib.load(p) for p in (baseline_path, clean_path, empirical_path)
        ]
        baseline["model"].n_jobs = 2
        records, lengths = read_pedsim_with_lengths(args.segments, 3346.298892538)
        truth = [r["label"] for r in records]
        groups = [pedigree(str(r["indID1"])) for r in records]
        classes = np.asarray(LABELS)
        pbase = align_probabilities(
            baseline["model"].classes_,
            baseline["model"].predict_proba(matrix(records, baseline["feature_names"])),
            classes,
        )
        pclean = align_probabilities(
            clean["candidate_model"].classes_,
            clean["candidate_model"].predict_proba(
                matrix(records, clean["candidate_feature_names"])
            ),
            classes,
        )
        pclean = (
            float(clean["candidate_weight"]) * pclean
            + (1 - float(clean["candidate_weight"])) * pbase
        )
        result["runs"]["task10_clean"] = probability_metrics(
            classes, pclean, truth, groups, MEIOSES
        )
        # Reference is calibration replicate 0, never this fresh evaluation file.
        ref_path = root / "artifacts/task10_human_calibration/replicate_0.seg"
        reference = read_pedsim(ref_path, 3346.298892538)
        result["correction_reference_sha256"] = source_hash(ref_path)
        for intensity in (0.0, 0.6, 1.4):
            augmented = augment_empirical(
                records,
                lengths,
                empirical["background_profile"],
                EmpiricalRegime(intensity, 0.4, 1.1),
                np.random.default_rng(args.seed + int(intensity * 100)),
                3346.298892538,
            )
            corrected = correct_records(augmented, reference, 5, 3346.298892538)
            cls, prob = predict_probabilities(
                empirical["corrected_payload"],
                matrix(corrected, empirical["feature_names"]),
            )
            weight = float(empirical["corrected_model_weight"])
            if weight != 1:
                prob = weight * prob + (1 - weight) * align_probabilities(
                    baseline["model"].classes_,
                    baseline["model"].predict_proba(
                        matrix(corrected, baseline["feature_names"])
                    ),
                    cls,
                )
            result["runs"][f"task12_independent_reference_intensity_{intensity}"] = (
                probability_metrics(cls, prob, truth, groups, MEIOSES)
            )
        result["expected_pairs_before_zero_ibd_censoring"] = 7000
        result["model_sha256"] = {
            str(p.relative_to(root)): source_hash(p)
            for p in (baseline_path, clean_path, empirical_path)
        }
        result["limitations"].append(
            "Background is independently seeded segment addition with the already frozen unlabeled cohort bin profile; it does not reproduce correlated IBD from background demography."
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: {
                    name: value
                    for name, value in v.items()
                    if name
                    in (
                        "accuracy",
                        "balanced_accuracy",
                        "meioses_accuracy",
                        "pairs",
                        "ece_10_equal_width_bins",
                    )
                }
                for k, v in result["runs"].items()
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
