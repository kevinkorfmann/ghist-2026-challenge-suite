from __future__ import annotations

import argparse
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

if __package__:
    from scripts.consensus_task11_dog_pedigrees import METHOD_VERSION
else:
    from consensus_task11_dog_pedigrees import METHOD_VERSION


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
RELEASE_SOURCES = (
    "scripts/consensus_task11_dog_pedigrees.py",
    "scripts/freeze_task11_multifamily_inbreeding.py",
    "scripts/prepare_hapibd_vcf.py",
    "scripts/summarize_hapibd.py",
    "scripts/cluster_task11_families.py",
    "scripts/reconstruct_bonsai_pedigree.py",
    "scripts/reconstruct_task11_pedigrees.py",
    "scripts/fit_task11_dog_relationships.py",
    "src/ghist_relatedness/labels.py",
    "src/ghist_relatedness/ibd.py",
    "src/ghist_relatedness/validation.py",
)
RELEASE_ARTIFACTS = (
    "artifacts/task11_multifamily_inbreeding/consensus/task11_candidate.csv",
    "artifacts/task11_multifamily_inbreeding/consensus/summary.json",
    "artifacts/task11_multifamily_inbreeding/family_clusters/summary.json",
    "artifacts/task11_multifamily_inbreeding/pair_model_hbd5/summary.json",
    *(
        f"artifacts/task11_multifamily_inbreeding/pair_model_hbd{threshold}/{name}"
        for threshold in (2, 5, 10)
        for name in ("model.joblib", "predictions.json")
    ),
)
RELEASE_CHECKS = (
    "independent_complete_pipeline",
    "disjoint_from_training_and_development",
    "varied_pedigree_topologies",
    "matched_dog_founder_model",
    "permitted_missingness_validated",
    "unconditional_pair_coverage_scored",
    "predeclared_acceptance_criteria_passed",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require_release_validation(evidence_path: Path | None, root: Path = ROOT) -> dict:
    """Require independent whole-task evidence before opening legacy inputs.

    Hashes bind the attestation to reviewed protocol/results, source, models and
    candidate. This verifies provenance and declared gates, not scientific truth.
    The sibling-only development cohorts cannot qualify as confirmation.
    """
    if evidence_path is None:
        raise RuntimeError(
            "Task 11 is not release-ready: --validation-evidence is required. "
            "Historical pair accuracy and topology-conditioned stability do not "
            "validate the complete pipeline."
        )
    evidence = json.loads(evidence_path.read_text())
    if evidence.get("schema_version") != "task11-complete-validation-v1":
        raise ValueError("Invalid Task 11 complete-validation evidence schema")
    if evidence.get("method_version") != METHOD_VERSION:
        raise ValueError(
            "Validation evidence does not match the selected Task 11 method"
        )
    checks = evidence.get("checks", {})
    if any(checks.get(check) is not True for check in RELEASE_CHECKS):
        raise ValueError("Task 11 complete-validation checks did not all pass")
    if set(evidence.get("relationship_labels_evaluated", [])) != set(
        DOG_LABELS.values()
    ):
        raise ValueError(
            "Task 11 evidence must evaluate all eight accepted relationships"
        )
    cohort_ids = evidence.get("confirmation_cohort_ids", [])
    if (
        not cohort_ids
        or any(not isinstance(value, str) or not value for value in cohort_ids)
        or len(set(cohort_ids)) != len(cohort_ids)
        or any(
            seed in value
            for value in cohort_ids
            for seed in ("2026112711", "2026112811")
        )
    ):
        raise ValueError(
            "Task 11 confirmation cohorts must be identified, unique and new"
        )
    for key, required in (
        ("source_sha256", RELEASE_SOURCES),
        ("artifact_sha256", RELEASE_ARTIFACTS),
    ):
        recorded = evidence.get(key, {})
        if not set(required).issubset(recorded):
            raise ValueError(f"Task 11 evidence is missing required {key} bindings")
        for relative in required:
            if sha256(root / relative) != recorded[relative]:
                raise ValueError(f"Task 11 validation hash mismatch: {relative}")
    records = {
        key: evidence.get(key, {})
        for key in ("validation_protocol", "validation_report")
    }
    records.update(
        {
            key: evidence.get("execution_assets", {}).get(key, {})
            for key in ("dog_recombination_map", "hapibd_jar", "inference_config")
        }
    )
    for key, record in records.items():
        relative = Path(record.get("path", ""))
        if (
            not str(record.get("path", ""))
            or relative.is_absolute()
            or ".." in relative.parts
        ):
            raise ValueError(f"Task 11 evidence needs a repository-relative {key} path")
        if sha256(root / relative) != record.get("sha256"):
            raise ValueError(f"Task 11 {key} hash mismatch")
    return {
        "path": str(evidence_path),
        "sha256": sha256(evidence_path),
        "method_version": evidence["method_version"],
        "confirmation_cohort_ids": cohort_ids,
        "checks": checks,
        "provenance_verified": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation-evidence", type=Path)
    args = parser.parse_args()
    # Must precede legacy artifact reads, genotype inspection and any writes.
    complete_validation = require_release_validation(args.validation_evidence)
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
        "fallback_hbd_threshold_stability": (
            consensus["fallback_hbd_sensitivity_pairs"] == 0
            or (
                consensus["fallback_hbd_sensitivity_stability"] is not None
                and consensus["fallback_hbd_sensitivity_stability"] >= 0.93
            )
        ),
        "sibling_components_consistent": consensus["sibling_component_conflicts"] == 0,
        "all_final_labels_have_pair_probability_support": consensus[
            "final_labels_with_zero_pair_probability"
        ]
        == 0,
        "audited_consensus_method": consensus["method_version"] == METHOD_VERSION,
        "independent_complete_pipeline_validation": complete_validation[
            "provenance_verified"
        ],
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

    counts = [consensus["label_counts"].get(label, 0) for label in labels]
    axes[1].barh(range(len(short)), counts, color="#0072B2")
    axes[1].set_yticks(range(len(short)), short)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Frozen pair calls")
    axes[1].set_title(f"{len(clusters['families'])}-family consensus")
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
            "method_version": METHOD_VERSION,
            "family_partition": "stable graph clustering over long-segment sharing",
            "pair_model": "dog-map Ped-sim forest using IBD1/IBD2 totals and segment-length shape",
            "inbreeding_control": "pairwise union-HBD masking at 2, 5, and 10 cM",
            "topology_consensus": "accepted pedigree tuples, meiosis-constrained pair fallbacks, and full-sibling component consistency",
        },
        "validation": {
            "complete_pipeline_evidence": complete_validation,
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
