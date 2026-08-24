from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODEL_ROOT = ROOT / "artifacts/task9_ibd_segments/model_v3_shape"
SUMMARY = MODEL_ROOT / "summary.json"
RAW_PREDICTIONS = MODEL_ROOT / "task9_testing_predictions.csv"
CANDIDATE = MODEL_ROOT / "task9_candidate.csv"
FREEZE = ROOT / "experiments/ibd_segments/task9_freeze.json"
FIGURE_ROOT = ROOT / "figures/task9_ibd_segments/task9_confusion"

MEIOSES = {
    "siblings": 2,
    "half_sibling": 2,
    "grandparent_grandchild": 2,
    "avuncular": 3,
    "double_cousin": 4,
    "full_1st_cousin": 4,
    "half_1st_cousin": 4,
    "full_2nd_cousin": 6,
    "half_2nd_cousin": 6,
    "full_3rd_cousin": 8,
    "full_4th_cousin": 10,
    "full_5th_cousin": 12,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_candidate() -> int:
    with RAW_PREDICTIONS.open() as source:
        rows = list(csv.DictReader(source))
    with CANDIDATE.open("w", newline="") as destination:
        writer = csv.writer(destination)
        writer.writerow(("indID1", "indID2", "meioses_count", "relation"))
        for row in rows:
            relation = row.get("relation", row.get("relationship"))
            if relation not in MEIOSES:
                raise ValueError(f"Unaccepted relation: {relation}")
            writer.writerow((row["indID1"], row["indID2"], MEIOSES[relation], relation))
    return len(rows)


def meiosis_accuracy(labels: list[str], confusion: np.ndarray) -> float:
    correct = sum(
        confusion[left, right]
        for left, label1 in enumerate(labels)
        for right, label2 in enumerate(labels)
        if MEIOSES[label1] == MEIOSES[label2]
    )
    return float(correct / np.sum(confusion))


def make_figure(labels: list[str], confusion: np.ndarray) -> None:
    normalized = confusion / confusion.sum(axis=1, keepdims=True)
    abbreviations = {
        "avuncular": "AV",
        "double_cousin": "DC",
        "full_1st_cousin": "1C",
        "full_2nd_cousin": "2C",
        "full_3rd_cousin": "3C",
        "full_4th_cousin": "4C",
        "full_5th_cousin": "5C",
        "grandparent_grandchild": "GP",
        "half_1st_cousin": "H1C",
        "half_2nd_cousin": "H2C",
        "half_sibling": "HS",
        "siblings": "SIB",
    }
    ticks = [abbreviations[label] for label in labels]
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 8,
            "font.weight": "normal",
            "axes.labelweight": "normal",
            "axes.titleweight": "normal",
            "pdf.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    image = ax.imshow(normalized, cmap="cividis", vmin=0, vmax=1, aspect="equal")
    ax.set_xticks(np.arange(len(labels)), ticks, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(labels)), ticks)
    ax.set_xlabel("Predicted relation")
    ax.set_ylabel("Simulated relation")
    ax.set_title("Task 9 untouched relationship-label evaluation")
    for row in range(normalized.shape[0]):
        for column in range(normalized.shape[1]):
            value = normalized[row, column]
            if value >= 0.10:
                ax.text(
                    column,
                    row,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    color="white" if value < 0.55 else "black",
                    fontsize=6.5,
                    fontweight="normal",
                )
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Row-normalized fraction")
    fig.tight_layout()
    FIGURE_ROOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_ROOT.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(FIGURE_ROOT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    summary = json.loads(SUMMARY.read_text())
    pairs = write_candidate()
    locked = summary["locked_evaluation"]
    labels = locked["labels"]
    confusion = np.asarray(locked["confusion_matrix"], dtype=float)
    meiosis_score = meiosis_accuracy(labels, confusion)
    make_figure(labels, confusion)
    with CANDIDATE.open() as handle:
        candidate_rows = list(csv.DictReader(handle))
    pair_keys = {(row["indID1"], row["indID2"]) for row in candidate_rows}
    gates = {
        "fresh_evaluation_has_at_least_10000_pairs": locked["pairs"] >= 10_000,
        "exact_label_accuracy_at_least_80pct": locked["accuracy"] >= 0.80,
        "balanced_accuracy_at_least_75pct": locked["balanced_accuracy"] >= 0.75,
        "meiosis_accuracy_at_least_84pct": meiosis_score >= 0.84,
        "testing_pair_count_matches_input": pairs == summary["observed_pairs"],
        "testing_pairs_unique": len(pair_keys) == pairs,
        "official_schema_exact": list(candidate_rows[0])
        == ["indID1", "indID2", "meioses_count", "relation"],
        "all_testing_labels_accepted": all(
            row["relation"] in MEIOSES for row in candidate_rows
        ),
    }
    gates["passed"] = all(gates.values())
    if not gates["passed"]:
        raise RuntimeError(f"Task 9 robustness gate failed: {gates}")
    figure_paths = (FIGURE_ROOT.with_suffix(".pdf"), FIGURE_ROOT.with_suffix(".png"))
    freeze = {
        "challenge": "Relatedness from identical-by-descent segments",
        "status": "method frozen; no submission",
        "submission_count": 0,
        "testing_prediction": {
            "pairs": pairs,
            "path": str(CANDIDATE.relative_to(ROOT)),
            "label_counts": summary["observed_label_counts"],
        },
        "selected_method": {
            "model": summary["selected_model"],
            "features": len(summary["feature_names"]),
            "description": "sex-map-aware totals, ordered segment shape, binned IBD1/IBD2 spectra, and chromosome sharing",
        },
        "validation": {
            "protocol": summary["protocol"],
            "validation": summary["validation_candidates"][summary["selected_model"]],
            "untouched_evaluation": {
                **locked,
                "meiosis_count_accuracy": meiosis_score,
            },
            "testing_confidence": summary["observed_confidence"],
            "robustness_gate": gates,
        },
        "data_provenance": {
            "testing_csv_sha256": "f1afa39793c32efc7da67129807f18f9de1f3c09e4a7f81e897a2cb3e1c55746",
            "final_csv_sealed_sha256": "9e5d2280397625b0db6ec68ef5eaf882df603736b79429eec9aed6eba7c5b851",
        },
        "source_artifacts": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (SUMMARY, CANDIDATE)
        },
        "figures": {
            str(path.relative_to(ROOT)): sha256(path) for path in figure_paths
        },
        "betty_job_ids": ["7824555", "7824935"],
    }
    FREEZE.parent.mkdir(parents=True, exist_ok=True)
    FREEZE.write_text(json.dumps(freeze, indent=2) + "\n")
    print(json.dumps(freeze, indent=2))


if __name__ == "__main__":
    main()
