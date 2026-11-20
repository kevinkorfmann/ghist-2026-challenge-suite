from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from ghist_single_sweep.recovery import validate_blind_records


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "artifacts/task7_single_sweep_bgs/audit_summary.json"
LOCALIZATION = ROOT / "artifacts/task7_single_sweep_bgs/localization/summary.json"
RAISD = ROOT / "artifacts/task7_single_sweep_bgs/raisd_summary.json"
FREEZE = ROOT / "experiments/single_sweep_bgs/task7_freeze.json"
BED = ROOT / "artifacts/task7_single_sweep_bgs/task7_candidate.bed"
FIGURE_ROOT = ROOT / "figures/task7_single_sweep_bgs/task7_localization"
HALF_WIDTH = 80_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_matched_recovery(recovery: dict) -> None:
    validate_blind_records(recovery.get("records", []))
    if recovery.get("inference_method") != "median_multiscale_raisd_gap_ld":
        raise ValueError(
            "Task 7 needs recovery of the exact ten-estimate inference method"
        )
    for record in recovery["records"]:
        if (
            record.get("background_selection") is not True
            or record.get("demography") != "bottleneck_growth"
        ):
            raise ValueError(
                "Task 7 requires independently simulated bottleneck/growth plus BGS recovery"
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Historical testing export; requires new matched BGS recovery"
    )
    parser.add_argument("--matched-recovery", type=Path, required=True)
    args = parser.parse_args()
    recovery_path = args.matched_recovery
    audit = json.loads(AUDIT.read_text())
    localization = json.loads(LOCALIZATION.read_text())
    raisd = json.loads(RAISD.read_text())
    recovery = json.loads(recovery_path.read_text())
    validate_matched_recovery(recovery)

    evidence: list[tuple[str, str, float]] = []
    for scale in ("25000", "50000", "100000", "250000", "500000"):
        evidence.append(
            (
                f"Composite {int(scale) // 1000} kb",
                "Composite scan",
                float(audit["scales"][scale]["top_windows"][0]["center"]),
            )
        )
    evidence.append(("RAiSD", "Independent scan", float(raisd["maximum"]["position"])))
    evidence.append(
        (
            "Largest variant gap",
            "Fine localization",
            float(localization["largest_variant_gap"]["midpoint"]),
        )
    )
    for scale in ("50000", "100000", "250000"):
        evidence.append(
            (
                f"Split-flank linkage {int(scale) // 1000} kb",
                "Fine localization",
                float(
                    localization["flank_ld_scans"][scale]["maximum"][
                        "candidate_position"
                    ]
                ),
            )
        )
    positions = np.asarray([position for _, _, position in evidence])
    center = int(np.median(positions))
    start = center - HALF_WIDTH
    end = center + HALF_WIDTH
    distant_mu = max(
        peak["mu"]
        for peak in raisd["separated_peaks_100kb"]
        if abs(peak["position"] - center) > 2_000_000
    )
    signal_ratio = float(raisd["maximum"]["mu"] / distant_mu)
    recovery_observed_like = recovery["untouched_evaluation"]["observed_like_consensus"]
    gates = {
        "exact_scale_hard_sweep_recovery_passed": bool(
            recovery["robustness_gate"]["passed"]
        ),
        "all_observed_localizers_inside_interval": bool(
            np.all(np.abs(positions - center) <= HALF_WIDTH)
        ),
        "localizer_span_at_most_160kb": bool(np.ptp(positions) <= 160_000),
        "independent_peak_at_least_fourfold_background": signal_ratio >= 4.0,
        "recovery_observed_like_max_error_inside_interval": bool(
            recovery_observed_like["quantiles"]["1.0"] <= HALF_WIDTH
        ),
    }
    gates["passed"] = all(gates.values())
    if not gates["passed"]:
        raise RuntimeError(f"Task 7 robustness gate failed: {gates}")

    BED.parent.mkdir(parents=True, exist_ok=True)
    BED.write_text(f"7\t{start}\t{end}\n")

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
    colors = {
        "Composite scan": "#0072B2",
        "Independent scan": "#D55E00",
        "Fine localization": "#009E73",
    }
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    ax.axvspan(start / 1e6, end / 1e6, color="#56B4E9", alpha=0.18, linewidth=0)
    ax.axvline(center / 1e6, color="black", linewidth=1.0)
    for index, (label, group, position) in enumerate(evidence):
        ax.scatter(
            position / 1e6, index, color=colors[group], marker="o", s=27, zorder=3
        )
    ax.set_yticks(np.arange(len(evidence)), [label for label, _, _ in evidence])
    ax.invert_yaxis()
    ax.set_xlabel("Chromosome 7 position (Mb)")
    ax.set_title("Task 7 sweep localization evidence")
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.6)
    ax.spines[["top", "right", "left"]].set_visible(False)
    handles = [
        plt.Line2D([], [], marker="o", linestyle="none", color=color, label=group)
        for group, color in colors.items()
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right")
    fig.tight_layout()
    FIGURE_ROOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_ROOT.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(FIGURE_ROOT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    freeze = {
        "challenge": "Single sweep with background selection",
        "status": "method frozen; no submission",
        "submission_count": 0,
        "testing_prediction": {
            "chromosome": "7",
            "start": start,
            "end": end,
            "center": center,
        },
        "selected_method": {
            "localizer": "median of ten multiscale, independent, gap, and split-flank estimates",
            "interval_policy": "80 kb half-width covering every accepted observed localizer",
        },
        "validation": {
            "localizer_count": len(evidence),
            "localizer_span_bp": float(np.ptp(positions)),
            "independent_peak_to_background_ratio": signal_ratio,
            "matched_bgs_recovery": {
                "replicates": recovery["replicates"],
                "untouched_observed_like_replicates": recovery[
                    "evaluation_consensus_replicates"
                ],
                "median_absolute_error": recovery_observed_like[
                    "median_absolute_error"
                ],
                "maximum_absolute_error": recovery_observed_like["quantiles"]["1.0"],
                "coverage_within_50kb": recovery_observed_like["within_50kb"],
            },
            "candidate_hit_score_if_truth_inside": 1.0 / (1.0 + HALF_WIDTH / 100_000),
            "robustness_gate": gates,
        },
        "observed_estimates": {label: position for label, _, position in evidence},
        "data_provenance": {
            "testing_vcf_sha256": audit["files"]["vcf"]["sha256"],
            "accessibility_bed_sha256": audit["files"]["accessible_bed"]["sha256"],
            "accepted_biallelic_sites": audit["accepted_biallelic_sites"],
            "haplotypes": audit["haplotypes"],
        },
        "source_artifacts": {
            str(path): sha256(path)
            for path in (AUDIT, LOCALIZATION, RAISD, recovery_path, BED)
        },
        "figures": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                FIGURE_ROOT.with_suffix(".pdf"),
                FIGURE_ROOT.with_suffix(".png"),
            )
        },
        "betty_job_ids": ["7824764", "7824777", "7824787"],
    }
    FREEZE.parent.mkdir(parents=True, exist_ok=True)
    FREEZE.write_text(json.dumps(freeze, indent=2) + "\n")
    print(json.dumps(freeze, indent=2))


if __name__ == "__main__":
    main()
