from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

import joblib
import numpy as np

try:
    from scripts.fit_human_background_corrected import (
        correct_records,
        corrected_feature_names,
        evaluate_blend,
        observed_predictions,
    )
    from scripts.fit_human_background_robust import (
        fit_direct_models,
        fit_hierarchical_model,
    )
    from scripts.fit_human_family_relationships import (
        BIN_EDGES,
        matrix,
        read_observed,
        summarize_raw_lengths,
    )
    from scripts.validate_human_background_correction import read_pedsim_with_lengths
except ModuleNotFoundError:
    from fit_human_background_corrected import (
        correct_records,
        corrected_feature_names,
        evaluate_blend,
        observed_predictions,
    )
    from fit_human_background_robust import (
        fit_direct_models,
        fit_hierarchical_model,
    )
    from fit_human_family_relationships import (
        BIN_EDGES,
        matrix,
        read_observed,
        summarize_raw_lengths,
    )
    from validate_human_background_correction import read_pedsim_with_lengths


@dataclass(frozen=True)
class EmpiricalRegime:
    intensity: float
    pair_heterogeneity: float
    within_bin_shape: float = 1.0


SELECTION_REGIMES = (
    EmpiricalRegime(0.0, 0.0),
    EmpiricalRegime(0.50, 0.18, 0.8),
    EmpiricalRegime(0.75, 0.28, 1.2),
    EmpiricalRegime(1.00, 0.38, 1.0),
    EmpiricalRegime(1.25, 0.48, 0.7),
    EmpiricalRegime(1.50, 0.58, 1.3),
)

DEVELOPMENT_REGIMES = (
    EmpiricalRegime(0.0, 0.0),
    EmpiricalRegime(0.40, 0.22, 1.3),
    EmpiricalRegime(0.65, 0.32, 0.7),
    EmpiricalRegime(0.90, 0.42, 1.1),
    EmpiricalRegime(1.15, 0.52, 0.9),
    EmpiricalRegime(1.40, 0.62, 1.2),
)

FINAL_REGIMES = (
    EmpiricalRegime(0.55, 0.25, 0.75),
    EmpiricalRegime(0.80, 0.35, 1.25),
    EmpiricalRegime(1.00, 0.45, 0.85),
    EmpiricalRegime(1.20, 0.55, 1.15),
    EmpiricalRegime(1.45, 0.65, 1.0),
)


def estimate_bin_profile(
    observed: list[dict[str, object]], reference: list[dict[str, object]]
) -> list[float]:
    rates = []
    for index in range(len(BIN_EDGES) - 1):
        key = f"length_bin_{index}"
        rates.append(
            max(
                0.0,
                float(np.median([float(row[key]) for row in observed]))
                - float(np.median([float(row[key]) for row in reference])),
            )
        )
    return rates


def sample_bin_lengths(
    count: int,
    lower: float,
    upper: float,
    shape: float,
    rng: np.random.Generator,
) -> np.ndarray:
    if count == 0:
        return np.asarray([], dtype=float)
    if np.isinf(upper):
        return lower + rng.exponential(max(10.0, lower / 2.0), size=count)
    fractions = rng.beta(shape, shape, size=count)
    return lower + fractions * (upper - lower)


def augment_empirical(
    records: list[dict[str, object]],
    source_lengths: dict[tuple[str, str], list[float]],
    profile: list[float],
    regime: EmpiricalRegime,
    rng: np.random.Generator,
    genome_cm: float,
) -> list[dict[str, object]]:
    output = []
    for record in records:
        if regime.pair_heterogeneity > 0:
            effect = float(
                rng.lognormal(
                    -0.5 * regime.pair_heterogeneity**2,
                    regime.pair_heterogeneity,
                )
            )
        else:
            effect = 1.0
        background = []
        for rate, lower, upper in zip(
            profile, BIN_EDGES[:-1], BIN_EDGES[1:], strict=True
        ):
            count = int(rng.poisson(rate * regime.intensity * effect))
            background.extend(
                sample_bin_lengths(
                    count,
                    float(lower),
                    float(upper),
                    regime.within_bin_shape,
                    rng,
                ).tolist()
            )
        pair = tuple(sorted((str(record["indID1"]), str(record["indID2"]))))
        row: dict[str, object] = {
            "indID1": record["indID1"],
            "indID2": record["indID2"],
            "label": record["label"],
        }
        lengths = source_lengths[pair] + background
        row.update(summarize_raw_lengths([("0", value) for value in lengths], genome_cm))
        output.append(row)
    return output


def training_set(
    paths: list[Path],
    profile: list[float],
    augmentations: int,
    seed: int,
    minimum_cm: int,
    genome_cm: float,
):
    rng = np.random.default_rng(seed)
    output = []
    regimes = []
    for path in paths:
        reference, lengths = read_pedsim_with_lengths(path, genome_cm)
        output.extend(reference)
        for _ in range(augmentations):
            regime = EmpiricalRegime(
                intensity=float(rng.uniform(0.25, 1.75)),
                pair_heterogeneity=float(rng.uniform(0.1, 0.7)),
                within_bin_shape=float(rng.uniform(0.6, 1.5)),
            )
            # A small long-segment perturbation prevents dependence on a single
            # empirical tail while the bin profile controls the dominant short IBD.
            augmented = augment_empirical(
                reference, lengths, profile, regime, rng, genome_cm
            )
            output.extend(
                correct_records(augmented, reference, minimum_cm, genome_cm)
            )
            regimes.append(asdict(regime))
    return output, regimes


def regime_datasets(
    path: Path,
    profile: list[float],
    regimes: tuple[EmpiricalRegime, ...],
    seed: int,
    minimum_cm: int,
    genome_cm: float,
):
    reference, lengths = read_pedsim_with_lengths(path, genome_cm)
    output = []
    for index, regime in enumerate(regimes):
        augmented = augment_empirical(
            reference,
            lengths,
            profile,
            regime,
            np.random.default_rng(seed + index),
            genome_cm,
        )
        output.append(
            (
                regime,
                correct_records(augmented, reference, minimum_cm, genome_cm),
            )
        )
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-model", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, nargs="+", required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--development-evaluation", type=Path, required=True)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--background-reference", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--minimum-cm", type=int, default=5)
    parser.add_argument("--augmentations", type=int, default=4)
    parser.add_argument("--genome-cm", type=float, default=3346.298892538)
    parser.add_argument("--seed", type=int, default=122226)
    args = parser.parse_args()

    observed = read_observed(args.observed)
    background_reference = read_observed(args.background_reference)
    profile = estimate_bin_profile(observed, background_reference)
    baseline = joblib.load(args.baseline_model)
    names = corrected_feature_names(args.minimum_cm)
    train, training_regimes = training_set(
        args.calibration,
        profile,
        args.augmentations,
        args.seed,
        args.minimum_cm,
        args.genome_cm,
    )
    train_x = matrix(train, names)
    train_y = np.asarray([str(record["label"]) for record in train])
    candidates = fit_direct_models(train_x, train_y, args.seed)
    candidates.append(fit_hierarchical_model(train_x, train_y, args.seed + 100))
    selection = regime_datasets(
        args.selection,
        profile,
        SELECTION_REGIMES,
        args.seed + 1000,
        args.minimum_cm,
        args.genome_cm,
    )
    results = []
    fitted = {}
    for name, payload in candidates:
        fitted[name] = payload
        for weight in (0.25, 0.50, 0.75, 1.0):
            result = evaluate_blend(
                payload,
                baseline["model"],
                baseline["feature_names"],
                selection,
                names,
                weight,
            )
            result.update({"name": name, "corrected_model_weight": weight})
            results.append(result)
    results.sort(
        key=lambda row: (
            row["mean_exact_accuracy"],
            row["worst_exact_accuracy"],
            row["mean_meioses_accuracy"],
        ),
        reverse=True,
    )
    winner = results[0]
    selected_payload = fitted[str(winner["name"])]
    weight = float(winner["corrected_model_weight"])
    development = regime_datasets(
        args.development_evaluation,
        profile,
        DEVELOPMENT_REGIMES,
        args.seed + 2000,
        args.minimum_cm,
        args.genome_cm,
    )
    development_result = evaluate_blend(
        selected_payload,
        baseline["model"],
        baseline["feature_names"],
        development,
        names,
        weight,
    )
    corrected_observed = correct_records(
        observed, background_reference, args.minimum_cm, args.genome_cm
    )
    predictions = observed_predictions(
        selected_payload,
        baseline["model"],
        baseline["feature_names"],
        corrected_observed,
        names,
        weight,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "corrected_payload": selected_payload,
            "feature_names": names,
            "corrected_model_weight": weight,
            "background_profile": profile,
        },
        args.out_dir / "model.joblib",
    )
    (args.out_dir / "predictions.json").write_text(json.dumps(predictions, indent=2) + "\n")
    with (args.out_dir / "candidate.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("indID1", "indID2", "relation"))
        writer.writeheader()
        writer.writerows(
            {key: row[key] for key in writer.fieldnames} for row in predictions
        )
    summary = {
        "state": "DEVELOPMENT",
        "design": "unlabeled empirical background-bin simulation followed by corrected-domain classification",
        "background_profile": profile,
        "background_profile_bin_edges": list(BIN_EDGES[:-1]) + ["infinity"],
        "training_pairs": len(train),
        "training_regimes": training_regimes,
        "selection": results,
        "selected_model": winner["name"],
        "corrected_model_weight": weight,
        "development_evaluation": development_result,
        "observed_pairs": len(predictions),
        "observed_label_counts": dict(Counter(row["relation"] for row in predictions)),
        "final_holdout_evaluated": False,
        "submission_count": 0,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {
                "background_profile": profile,
                "selected_model": summary["selected_model"],
                "corrected_model_weight": weight,
                "selection_winner": winner,
                "development_evaluation": development_result,
                "observed_label_counts": summary["observed_label_counts"],
                "submission_count": 0,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
