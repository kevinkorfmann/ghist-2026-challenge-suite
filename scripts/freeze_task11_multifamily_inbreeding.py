from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from ghist_relatedness.labels import DOG_LABELS
from ghist_submission.validate import (
    vcf_samples,
    within_family_pairs,
    validate_submission,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_ROOT = ROOT / "artifacts/task11_multifamily_inbreeding"
CLUSTERS = ARTIFACT_ROOT / "family_clusters/summary.json"
PAIR_MODEL = ARTIFACT_ROOT / "pair_model_hbd5/summary.json"
CONSENSUS = ARTIFACT_ROOT / "consensus/summary.json"
CANDIDATE = ARTIFACT_ROOT / "consensus/task11_candidate.csv"
FREEZE = ROOT / "experiments/multifamily_inbreeding/task11_freeze.json"
FIGURE_ROOT = ROOT / "figures/task11_multifamily_inbreeding/task11_validation"
TESTING = (
    ROOT
    / "data/raw/relatedness/multifamily_inbreeding/testing/relatedness3.testing.vcf.gz"
)
FINAL = (
    ROOT
    / "data/holdout_final/relatedness/multifamily_inbreeding/relatedness3.final.vcf.gz"
)
FINAL_SHA256 = "5bbdd1a850f65308ed2376a159cbe63c7540e76eab9e91bb77b75c127a2fc553"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    clusters = json.loads(CLUSTERS.read_text())
    pair_model = json.loads(PAIR_MODEL.read_text())
    consensus = json.loads(CONSENSUS.read_text())
    with CANDIDATE.open() as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        header = reader.fieldnames

    pairs = [tuple(sorted((row["indID1"], row["indID2"]))) for row in rows]
    accepted = set(DOG_LABELS.values())
    samples = vcf_samples(TESTING)
    expected = within_family_pairs(clusters["families"], samples)
    coverage = validate_submission(CANDIDATE, 11, "testing", expected)
    gates = {
        "all_input_samples_partitioned": sum(
            len(family) for family in clusters["families"]
        )
        == len(samples),
        "family_partition_stable": clusters["minimum_adjusted_rand"] == 1.0
        and clusters["stability_runs"] >= 16,
        "fresh_pair_evaluation_above_60_percent": pair_model["fresh_evaluation"][
            "balanced_accuracy"
        ]
        >= 0.60,
        "official_relation_only_schema": header == ["indID1", "indID2", "relation"],
        "all_within_family_pairs_present": coverage["input_pair_coverage_verified"],
        "pairs_unique": len(set(pairs)) == len(pairs),
        "all_labels_accepted": {row["relation"] for row in rows}.issubset(accepted),
        "hbd_threshold_stability": consensus["hbd_sensitivity_stability"] >= 0.98,
        "fallback_hbd_threshold_stability": consensus[
            "fallback_hbd_sensitivity_stability"
        ]
        >= 0.93,
        "sibling_components_consistent": consensus["sibling_component_conflicts"] == 0,
    }
    gates["passed"] = all(gates.values())
    if not gates["passed"]:
        raise RuntimeError(f"Task 11 robustness gate failed: {gates}")

    labels = pair_model["labels"]
    confusion = np.asarray(pair_model["fresh_evaluation"]["confusion"], dtype=float)
    normalized = np.divide(
        confusion,
        confusion.sum(axis=1, keepdims=True),
        out=np.zeros_like(confusion),
        where=confusion.sum(axis=1, keepdims=True) > 0,
    )
    short = ["sib", "half-sib", "av", "half-av", "av-1R", "av-2R", "1C", "half-1C"]
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 7,
            "font.weight": "normal",
            "axes.labelweight": "normal",
            "axes.titleweight": "normal",
            "pdf.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(
        1, 2, figsize=(7.2, 3.1), gridspec_kw={"width_ratios": [1.15, 0.85]}
    )
    image = axes[0].imshow(normalized, vmin=0, vmax=1, cmap="Blues")
    for row in range(normalized.shape[0]):
        for column in range(normalized.shape[1]):
            value = normalized[row, column]
            if value >= 0.10:
                axes[0].text(
                    column,
                    row,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    color="white" if value > 0.55 else "black",
                    fontsize=5.5,
                    fontweight="normal",
                )
    axes[0].set_xticks(range(len(short)), short, rotation=45, ha="right")
    axes[0].set_yticks(range(len(short)), short)
    axes[0].set_xlabel("Predicted")
    axes[0].set_ylabel("Simulated truth")
    axes[0].set_title("Untouched pairwise evaluation")
    fig.colorbar(image, ax=axes[0], fraction=0.046, pad=0.04, label="Row fraction")

    counts = [consensus["label_counts"][label] for label in labels]
    axes[1].barh(range(len(short)), counts, color="#0072B2")
    axes[1].set_yticks(range(len(short)), short)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Frozen pair calls")
    axes[1].set_title("Eight-family consensus")
    axes[1].spines[["top", "right"]].set_visible(False)
    axes[1].grid(axis="x", color="#D9D9D9", linewidth=0.5)
    fig.suptitle("Task 11 dog relatedness validation", fontweight="normal")
    fig.tight_layout()
    FIGURE_ROOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_ROOT.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(FIGURE_ROOT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    freeze = {
        "challenge": "Relatedness in a multifamily cohort with close inbreeding",
        "status": "method frozen; no submission",
        "submission_count": 0,
        "testing_prediction": {
            "pairs": len(rows),
            "families": len(clusters["families"]),
            "label_counts": consensus["label_counts"],
            "path": str(CANDIDATE.relative_to(ROOT)),
            "schema": header,
        },
        "selected_method": {
            "family_partition": "stable graph clustering over long-segment sharing",
            "pair_model": "dog-map Ped-sim forest using IBD1/IBD2 totals and segment-length shape",
            "inbreeding_control": "pairwise union-HBD masking at 2, 5, and 10 cM",
            "topology_consensus": "accepted pedigree tuples, meiosis-constrained pair fallbacks, and full-sibling component consistency",
        },
        "validation": {
            "input_pair_coverage": coverage,
            "fresh_pair_accuracy": pair_model["fresh_evaluation"]["accuracy"],
            "fresh_pair_balanced_accuracy": pair_model["fresh_evaluation"][
                "balanced_accuracy"
            ],
            "fresh_pair_pairs": pair_model["evaluation_pairs"],
            "minimum_family_partition_ari": clusters["minimum_adjusted_rand"],
            "family_partition_runs": clusters["stability_runs"],
            "accepted_topology_pairs": consensus["accepted_topology_pairs"],
            "pair_model_fallbacks": consensus["pair_model_fallbacks"],
            "hbd_threshold_stability": consensus["hbd_sensitivity_stability"],
            "fallback_hbd_threshold_stability": consensus[
                "fallback_hbd_sensitivity_stability"
            ],
            "degree_two_relation_vector_refinements": consensus[
                "degree_two_relation_vector_refinements"
            ],
            "full_sibling_component_propagations": consensus[
                "full_sibling_relation_vector_propagations"
            ],
            "robustness_gate": gates,
        },
        "data_provenance": {
            "testing_vcf_sha256": sha256(TESTING),
            "final_vcf_sealed_sha256": sha256(FINAL)
            if FINAL.exists()
            else FINAL_SHA256,
        },
        "source_artifacts": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (CLUSTERS, PAIR_MODEL, CONSENSUS, CANDIDATE)
        },
        "figures": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                FIGURE_ROOT.with_suffix(".pdf"),
                FIGURE_ROOT.with_suffix(".png"),
            )
        },
        "betty_job_ids": ["7824956", "7824985", "7824994", "7825033"],
    }
    FREEZE.parent.mkdir(parents=True, exist_ok=True)
    FREEZE.write_text(json.dumps(freeze, indent=2) + "\n")
    print(json.dumps(freeze, indent=2))


if __name__ == "__main__":
    main()
