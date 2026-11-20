"""Apply a frozen Task 9 model under the current accepted-label contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import tempfile

import joblib
import numpy as np

from ghist_relatedness.ibd import matrix, pair_features, read_segment_rows, read_sex_map
from ghist_relatedness.labels import IBD_LABEL_MEIOSES
from ghist_relatedness.task9_release import input_pairs, load_release, validate_payload
from ghist_submission.validate import validate_submission


def accepted_probabilities(classes, probabilities):
    classes = np.asarray(classes)
    probabilities = np.asarray(probabilities, dtype=float)
    if (
        classes.ndim != 1
        or len(set(classes)) != len(classes)
        or probabilities.ndim != 2
        or probabilities.shape[1] != len(classes)
        or not np.all(np.isfinite(probabilities))
        or np.any(probabilities < 0)
    ):
        raise ValueError("Invalid Task 9 model probability matrix")
    keep = np.asarray([str(label) in IBD_LABEL_MEIOSES for label in classes])
    values = probabilities[:, keep]
    total = values.sum(axis=1, keepdims=True)
    if (
        not values.shape[1]
        or np.any(total <= 0)
        or not np.all(np.isfinite(values))
        or np.any(values < 0)
    ):
        raise ValueError(
            "Model assigned invalid or zero probability to accepted Task 9 labels"
        )
    return classes[keep], values / total


def file_hash(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--model", type=Path, help="Explicit unversioned/legacy model")
    selection.add_argument("--release", type=Path, help="Checksum-pinned release manifest")
    parser.add_argument("--artifact-root", type=Path, default=Path("."))
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--sex-map", type=Path, required=True)
    parser.add_argument("--dataset", choices=("testing", "final"), required=True)
    parser.add_argument("--evaluation", type=Path)
    parser.add_argument("--genome-cm", type=float)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--pairs", type=Path, help="Optional independent expected-pair CSV")
    parser.add_argument("--input-provenance", type=Path, help="JSON with entity, version and sha256")
    parser.add_argument("--synthetic-input", action="store_true")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.workers < 1:
        raise ValueError("workers must be positive")
    manifest = None
    if args.release:
        manifest, args.model, args.genome_cm = load_release(
            args.release, args.artifact_root, args.sex_map, args.genome_cm
        )
    elif args.genome_cm is None:
        args.genome_cm = 3350.933402
    provenance = None
    if args.input_provenance:
        provenance = json.loads(args.input_provenance.read_text())
        expected_entity = "syn75202632" if args.dataset == "testing" else "syn75202631"
        if (
            provenance.get("entity") != expected_entity
            or type(provenance.get("version")) is not int
            or provenance["version"] < 1
            or provenance.get("sha256") != file_hash(args.observed)
        ):
            raise ValueError("Input provenance does not match the requested Task 9 dataset")
    if args.dataset == "final" and not args.synthetic_input and provenance is None:
        raise ValueError("Final inputs require recorded Synapse entity/version/checksum provenance")
    payload = joblib.load(args.model)
    model, features = validate_payload(payload, manifest)
    if hasattr(model, "n_jobs"):
        model.n_jobs = args.workers
    sex_map = read_sex_map(args.sex_map)
    segments = read_segment_rows(args.observed)
    expected_pairs = input_pairs(segments)
    if args.pairs:
        with args.pairs.open(newline="") as handle:
            supplied_pairs = {
                tuple(sorted((row["indID1"], row["indID2"])))
                for row in csv.DictReader(handle)
            }
        if supplied_pairs != expected_pairs:
            raise ValueError("Expected pairs differ from the segment input; zero-IBD/unlisted pairs are unsupported")
    rows = pair_features(segments, args.genome_cm, sex_map)
    x, _ = matrix(rows, features)
    raw = model.predict_proba(x)
    labels, probabilities = accepted_probabilities(model.classes_, raw)
    predictions = labels[np.argmax(probabilities, axis=1)]
    original = model.classes_[np.argmax(raw, axis=1)]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output = args.out_dir / f"task9_{args.dataset}_candidate.csv"
    with tempfile.NamedTemporaryFile(mode="w", newline="", dir=args.out_dir, delete=False) as handle:
        temporary_output = Path(handle.name)
        writer = csv.writer(handle)
        writer.writerow(("indID1", "indID2", "meioses_count", "relation"))
        for row, label in zip(rows, predictions, strict=True):
            writer.writerow(
                (row["indID1"], row["indID2"], IBD_LABEL_MEIOSES[str(label)], label)
            )
    try:
        validation = validate_submission(temporary_output, 9, args.dataset, expected_pairs)
    except Exception:
        temporary_output.unlink(missing_ok=True)
        raise
    report = {
        "dataset": args.dataset,
        "synthetic_input": args.synthetic_input,
        "release_id": manifest["release_id"] if manifest else None,
        "policy": "Use explicitly selected fixed model; renormalize accepted labels; no fitting or answer-based selection",
        "pair_domain": "unordered_pairs_present_in_segment_csv",
        "input_provenance": provenance,
        "changed_predictions": int(np.sum(original != predictions)),
        "validation": validation,
        "source_sha256": {
            str(p): file_hash(p) for p in (args.model, args.observed, args.sex_map)
        },
        "output_sha256": file_hash(temporary_output),
        "submission_count": 0,
    }
    if args.evaluation:
        from scripts.fit_task9_ibd_classifier import load_simulations

        evaluation = load_simulations([args.evaluation], args.genome_cm, sex_map)
        evaluation_x, _ = matrix(evaluation, features)
        evaluation_labels, probability = accepted_probabilities(
            model.classes_, model.predict_proba(evaluation_x)
        )
        predicted = evaluation_labels[np.argmax(probability, axis=1)]
        truth = np.asarray([row["label"] for row in evaluation])
        report["evaluation_on_current_label_set"] = {
            "pairs": len(truth),
            "accuracy": float(np.mean(predicted == truth)),
            "meioses_accuracy": float(
                np.mean(
                    [
                        IBD_LABEL_MEIOSES[str(a)] == IBD_LABEL_MEIOSES[str(b)]
                        for a, b in zip(predicted, truth, strict=True)
                    ]
                )
            ),
            "source": str(args.evaluation),
            "sha256": file_hash(args.evaluation),
        }
    temporary_output.replace(output)
    with tempfile.NamedTemporaryFile(mode="w", dir=args.out_dir, delete=False) as handle:
        temporary_report = Path(handle.name)
        handle.write(json.dumps(report, indent=2) + "\n")
    temporary_report.replace(args.out_dir / "summary.json")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
