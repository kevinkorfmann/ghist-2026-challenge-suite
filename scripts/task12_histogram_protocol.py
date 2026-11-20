"""Provenance and scoring contracts for the fixed Task12 histogram experiment."""

from __future__ import annotations

from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import re
import sys

import numpy as np

from scripts.fit_human_family_relationships import LABELS, MEIOSES, relationship_label


PACKAGES = ("numpy", "scipy", "scikit-learn", "joblib", "threadpoolctl", "PyYAML")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_frozen_protocol(path):
    """Bind the parsed configuration to the exact bytes read once at entry."""
    data = Path(path).read_bytes()
    return json.loads(data), hashlib.sha256(data).hexdigest()


def require_protocol_digest(path, expected):
    if sha256(path) != expected:
        raise ValueError("Task12 protocol changed after it was parsed")


def validate_seed_registry(protocol, root):
    registry = protocol["seed_registry"]
    expected_logs = {
        str(Path(name).with_suffix(".log"))
        for name in protocol["training"]["calibration"]
    }
    records = registry["calibration_logs"]
    if (
        len(records) != len(expected_logs)
        or {row["path"] for row in records} != expected_logs
    ):
        raise ValueError("Calibration seed logs do not match calibration inputs")
    calibration_seeds = []
    for row in records:
        path = Path(root) / row["path"]
        if (
            sha256(path) != row["sha256"]
            or protocol["input_sha256"].get(row["path"]) != row["sha256"]
        ):
            raise ValueError("Calibration seed-log provenance changed")
        seeds = [
            int(value)
            for value in re.findall(r"Random seed:\s*(\d+)", path.read_text())
        ]
        if seeds != [row["seed"]]:
            raise ValueError("Calibration simulation seed differs from registry")
        calibration_seeds.extend(seeds)
    settings = protocol["confirmation"]
    confirmation = settings["seeds"]
    augmentation = [
        settings["augmentation_seed"] + 100 * i + j
        for i in range(len(confirmation))
        for j in range(len(settings["regimes"]))
    ]
    if (
        registry["confirmation_seeds"] != confirmation
        or registry["augmentation_seeds"] != augmentation
    ):
        raise ValueError("Confirmation random streams differ from seed registry")
    streams = calibration_seeds + confirmation + augmentation
    if any(
        isinstance(seed, bool) or not isinstance(seed, int) or seed < 0
        for seed in streams
    ) or len(set(streams)) != len(streams):
        raise ValueError(
            "Calibration, confirmation and augmentation seeds must be distinct"
        )
    if set(confirmation + augmentation) & set(
        protocol["prior_confirmation_seeds_excluded"]
    ):
        raise ValueError("Confirmation random streams reuse a prior confirmation seed")


def write_json(path, value):
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"Refusing to replace frozen evidence: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    serialized = json.dumps(value, indent=2, allow_nan=False) + "\n"
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("x") as handle:
        handle.write(serialized)
    temporary.replace(path)


def runtime_manifest():
    files, versions = {}, {}
    for name in PACKAGES:
        distribution = importlib.metadata.distribution(name)
        versions[name] = distribution.version
        for entry in distribution.files or []:
            path = Path(distribution.locate_file(entry)).resolve()
            if path.is_file() and path.suffix != ".pyc":
                files[str(path)] = sha256(path)
    executable = Path(sys.executable).resolve()
    return {
        "python": platform.python_version(),
        "system": platform.system(),
        "machine": platform.machine(),
        "executable": str(executable),
        "executable_sha256": sha256(executable),
        "versions": versions,
        "package_files_sha256": dict(sorted(files.items())),
    }


def validate_protocol(protocol, root, snapshot):
    root, snapshot = Path(root), Path(snapshot)
    validate_seed_registry(protocol, root)
    for name, expected in protocol["source_sha256"].items():
        if sha256(snapshot / name) != expected:
            raise ValueError(f"Frozen Task12 source changed: {name}")
    for name, expected in protocol["input_sha256"].items():
        if sha256(root / name) != expected:
            raise ValueError(f"Frozen Task12 input/model changed: {name}")
    runtime_path = snapshot / protocol["runtime_manifest"]
    if sha256(runtime_path) != protocol["runtime_manifest_sha256"]:
        raise ValueError("Frozen Task12 runtime manifest changed")
    expected = json.loads(runtime_path.read_text())
    if runtime_manifest() != expected:
        raise ValueError("Task12 runtime differs from its frozen manifest")


def validate_pairs(records):
    pairs = set()
    for row in records:
        pair = (row["indID1"], row["indID2"])
        if pair[0] >= pair[1] or pair in pairs:
            raise ValueError(
                "Task12 simulated pairs must be unique, canonical and non-self"
            )
        if row["label"] not in LABELS or any(
            relationship_label(sample) != row["label"] for sample in pair
        ):
            raise ValueError("Task12 simulated relationship labels disagree")
        pairs.add(pair)
    if not pairs:
        raise ValueError("Empty Task12 simulation")
    return pairs


def score_predictions(predictions, records, expected_per_label=500):
    truth = {(row["indID1"], row["indID2"]): row["label"] for row in records}
    if len(truth) != len(records):
        raise ValueError("Duplicate truth pairs")
    seen, confidence, correct, degree, brier, nll = set(), [], [], [], [], []
    confusion = {label: Counter() for label in LABELS}
    for row in predictions:
        pair = (row["indID1"], row["indID2"])
        if pair not in truth or pair in seen:
            raise ValueError("Prediction pair coverage/uniqueness mismatch")
        seen.add(pair)
        label = truth[pair]
        probabilities = row["probabilities"]
        if set(probabilities) != set(LABELS):
            raise ValueError(
                "Probability labels do not match the complete Task12 domain"
            )
        values = np.asarray([probabilities[name] for name in LABELS], dtype=float)
        if (
            not np.isfinite(values).all()
            or (values < 0).any()
            or not np.isclose(values.sum(), 1, rtol=1e-8, atol=1e-10)
        ):
            raise ValueError("Invalid Task12 probability distribution")
        # Score the emitted decision, accepting any exact tied maximum rather
        # than imposing a second tie-break order inside the measurement code.
        predicted = row["relation"]
        if (
            predicted not in probabilities
            or probabilities[predicted] != values.max()
            or not np.isclose(
                row["confidence"], probabilities[predicted], rtol=1e-8, atol=1e-10
            )
        ):
            raise ValueError("Prediction label/confidence differs from probabilities")
        confidence.append(float(probabilities[predicted]))
        correct.append(predicted == label)
        degree.append(MEIOSES[predicted] == MEIOSES[label])
        brier.append(
            sum((probabilities[name] - (name == label)) ** 2 for name in LABELS)
        )
        nll.append(-np.log(max(probabilities[label], 1e-15)))
        confusion[label][predicted] += 1
    if seen != set(truth) or not seen:
        raise ValueError("Prediction pair coverage mismatch")
    counts = Counter(truth.values())
    if set(counts) != set(LABELS) or any(
        count > expected_per_label for count in counts.values()
    ):
        raise ValueError("Unexpected template-label denominator")
    correct, confidence = np.asarray(correct), np.asarray(confidence)
    ece = 0.0
    for index in range(10):
        selected = (confidence >= index / 10) & (
            (confidence < (index + 1) / 10) if index < 9 else (confidence <= 1)
        )
        if selected.any():
            ece += selected.mean() * abs(
                correct[selected].mean() - confidence[selected].mean()
            )
    expected_total = expected_per_label * len(LABELS)
    recall = {label: confusion[label][label] / counts[label] for label in LABELS}
    return {
        "segment_present_pairs": len(records),
        "expected_template_pairs": expected_total,
        "segment_present_by_label": dict(counts),
        "zero_ibd_omitted_by_label": {
            label: expected_per_label - counts[label] for label in LABELS
        },
        "accuracy": float(correct.mean()),
        "balanced_accuracy": float(np.mean(list(recall.values()))),
        "meioses_accuracy": float(np.mean(degree)),
        "recall_by_label": recall,
        "template_accuracy_counting_unpredicted_as_wrong": float(
            correct.sum() / expected_total
        ),
        "template_meioses_accuracy_counting_unpredicted_as_wrong": float(
            sum(degree) / expected_total
        ),
        "template_recall_counting_unpredicted_as_wrong": {
            label: confusion[label][label] / expected_per_label for label in LABELS
        },
        "confusion": {
            label: {predicted: confusion[label][predicted] for predicted in LABELS}
            for label in LABELS
        },
        "negative_log_likelihood": float(np.mean(nll)),
        "multiclass_brier_sum": float(np.mean(brier)),
        "ece_10_equal_width_bins": float(ece),
        "confidence_is_uncalibrated": True,
    }


def compare_runs(runs, comparator, settings):
    gains = [
        r["models"]["full_histogram_v2"]["accuracy"]
        - r["models"][comparator]["accuracy"]
        for r in runs
    ]
    degrees = [
        r["models"]["full_histogram_v2"]["meioses_accuracy"]
        - r["models"][comparator]["meioses_accuracy"]
        for r in runs
    ]
    seeds = sorted({row["seed"] for row in runs})
    by_seed, count_gains = {}, {}
    for seed in seeds:
        paired = [row for row in runs if row["seed"] == seed]
        denominators = {
            row["models"][name]["segment_present_pairs"]
            for row in paired
            for name in ("full_histogram_v2", comparator)
        }
        if len(denominators) != 1 or next(iter(denominators)) <= 0:
            raise ValueError("Paired Task12 seed/regime denominators must agree")
        denominator = next(iter(denominators))
        # Each regime reuses the same seed's pairs. Sum the integer diagonal
        # counts before dividing so a mathematical tie cannot become positive
        # through subtraction/averaging of rounded floating accuracies.
        difference = sum(
            row["models"]["full_histogram_v2"]["confusion"][label][label]
            - row["models"][comparator]["confusion"][label][label]
            for row in paired
            for label in LABELS
        )
        count_gains[str(seed)] = difference
        by_seed[str(seed)] = difference / (denominator * len(paired))
    gate = {
        "mean_exact_gain": np.mean(gains) >= settings["mean_exact_gain_min"],
        "positive_seed_fraction": sum(value > 0 for value in count_gains.values())
        / len(seeds)
        >= settings["positive_mean_gain_fraction_by_seed_min"],
        "worst_regime_exact_gain": min(gains)
        >= settings["worst_regime_exact_gain_min"],
        "mean_meioses_gain": np.mean(degrees) >= settings["mean_meioses_gain_min"],
    }
    gate = {name: bool(value) for name, value in gate.items()}
    return {
        "mean_exact_gain": float(np.mean(gains)),
        "mean_exact_gain_by_seed": by_seed,
        "exact_correct_count_gain_by_seed": count_gains,
        "minimum_regime_exact_gain": min(gains),
        "mean_meioses_gain": float(np.mean(degrees)),
        "gate": {**gate, "passed": all(gate.values())},
    }
