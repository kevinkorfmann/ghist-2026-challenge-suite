from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from select_background_sweep_candidates import read_raisd, separated_peaks


ROOT = Path(__file__).resolve().parents[1]
PROVISIONAL = ROOT / "artifacts/background_sweeps/provisional_candidates.json"
RAISD_REPORT = next(
    (ROOT / "artifacts/task8_multiple_sweeps_bgs/raisd").glob("RAiSD_Report.*")
)
EXONS = ROOT / "data/reference/sweeps/havana104_chr7_exons.bed"
RECOVERY = ROOT / "artifacts/task5_single_sweep/selection_exact_recovery_summary.json"
CANDIDATE = ROOT / "artifacts/background_sweeps/task8_candidate.bed"
FREEZE = ROOT / "experiments/multiple_sweeps_bgs/task8_freeze.json"
FIGURE_ROOT = ROOT / "figures/task8_multiple_sweeps_bgs/task8_peak_threshold"
HALF_WIDTH = 60_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_exons() -> list[tuple[int, int]]:
    intervals = []
    with EXONS.open() as handle:
        for line in handle:
            _, start, end = line.split()
            intervals.append((int(start), int(end)))
    return intervals


def main() -> None:
    provisional = json.loads(PROVISIONAL.read_text())
    recovery = json.loads(RECOVERY.read_text())
    candidates = sorted(provisional["task8"], key=lambda item: item["exonic_center"])
    centers = np.asarray([item["exonic_center"] for item in candidates], dtype=int)
    exons = read_exons()
    center_is_exonic = [
        any(start <= center < end for start, end in exons) for center in centers
    ]

    report = read_raisd(RAISD_REPORT)
    peak_indices = separated_peaks(report, 1_000_000)[:15]
    peak_mu = report[peak_indices, 6]
    threshold = float(provisional["calibration"]["task8_detection_threshold"])
    selected_count = int(np.sum(peak_mu >= threshold))
    weakest_selected = float(np.min(peak_mu[peak_mu >= threshold]))
    strongest_rejected = float(np.max(peak_mu[peak_mu < threshold]))
    observed_like = recovery["untouched_evaluation"]["observed_like_consensus"]
    gates = {
        "task7_background_null_calibrated": provisional["calibration"][
            "task7_background_maximum_mu_outside_2mb"
        ]
        > 0,
        "exactly_six_peaks_cross_threshold": selected_count == 6,
        "weakest_selected_exceeds_next_peak": weakest_selected > strongest_rejected,
        "all_selected_peaks_at_least_twice_background": all(
            item["raisd_to_task7_null_ratio"] >= 2.0 for item in candidates
        ),
        "all_consensus_centers_are_exonic": all(center_is_exonic),
        "all_exon_snaps_at_most_20kb": all(
            item["exon_snap_distance"] <= 20_000 for item in candidates
        ),
        "selected_centers_separated_by_at_least_1mb": bool(
            np.all(np.diff(centers) >= 1_000_000)
        ),
        "transferred_recovery_passed": bool(recovery["robustness_gate"]["passed"]),
        "transferred_observed_like_max_error_inside_interval": bool(
            observed_like["quantiles"]["1.0"] <= HALF_WIDTH
        ),
    }
    gates["passed"] = all(gates.values())
    if not gates["passed"]:
        raise RuntimeError(f"Task 8 robustness gate failed: {gates}")

    CANDIDATE.parent.mkdir(parents=True, exist_ok=True)
    with CANDIDATE.open("w") as handle:
        for center in centers:
            handle.write(f"7\t{center - HALF_WIDTH}\t{center + HALF_WIDTH}\n")

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
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    positions_mb = report[peak_indices, 0] / 1e6
    selected = peak_mu >= threshold
    ax.scatter(
        positions_mb[~selected],
        peak_mu[~selected],
        color="#999999",
        marker="o",
        s=24,
        label="Below threshold",
    )
    ax.scatter(
        positions_mb[selected],
        peak_mu[selected],
        color="#D55E00",
        marker="o",
        s=34,
        label="Frozen candidate",
        zorder=3,
    )
    ax.axhline(
        threshold,
        color="#0072B2",
        linestyle="--",
        linewidth=1.0,
        label="Twice Task 7 background",
    )
    for number, index in enumerate(np.flatnonzero(selected), start=1):
        ax.annotate(
            str(number),
            (positions_mb[index], peak_mu[index]),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            fontsize=7,
            fontweight="normal",
        )
    ax.set_xlabel("Chromosome 7 position (Mb)")
    ax.set_ylabel("RAiSD μ")
    ax.set_title("Task 8 separated-peak detection threshold")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    FIGURE_ROOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_ROOT.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(FIGURE_ROOT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    freeze = {
        "challenge": "Multiple sweeps with background selection",
        "status": "method frozen; no submission",
        "submission_count": 0,
        "testing_prediction": {
            "regions": [
                {
                    "chromosome": "7",
                    "start": int(center - HALF_WIDTH),
                    "end": int(center + HALF_WIDTH),
                    "center": int(center),
                }
                for center in centers
            ],
            "path": str(CANDIDATE.relative_to(ROOT)),
        },
        "selected_method": {
            "detector": "one-megabase-separated RAiSD peaks exceeding twice the Task 7 background maximum",
            "localizer": "median multiscale consensus snapped to the nearest Havana release 104 exon",
            "interval_policy": "60 kb half-width transferred from exact-scale hard-sweep recovery",
        },
        "validation": {
            "threshold_mu": threshold,
            "selected_peaks": selected_count,
            "weakest_selected_mu": weakest_selected,
            "strongest_rejected_mu": strongest_rejected,
            "weakest_to_strongest_rejected_ratio": weakest_selected
            / strongest_rejected,
            "candidate_to_background_ratios": [
                item["raisd_to_task7_null_ratio"] for item in candidates
            ],
            "exon_snap_distances": [item["exon_snap_distance"] for item in candidates],
            "transferred_untouched_observed_like": {
                "replicates": recovery["evaluation_consensus_replicates"],
                "median_absolute_error": observed_like["median_absolute_error"],
                "maximum_absolute_error": observed_like["quantiles"]["1.0"],
                "coverage_within_50kb": observed_like["within_50kb"],
            },
            "candidate_hit_score_if_each_truth_inside": 1.0
            / (1.0 + HALF_WIDTH / 100_000),
            "robustness_gate": gates,
        },
        "data_provenance": {
            "testing_vcf_sha256": "9b686087b7fc39cc7dab2f844c43af87b95932b9adebff982cf0e964da9db6e0",
            "final_vcf_sealed_sha256": "cec3b5172e27f54354934ade9437913f97a9e6a5c578aed691c993ed02050714",
            "exon_annotation_sha256": sha256(EXONS),
        },
        "source_artifacts": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (PROVISIONAL, RAISD_REPORT, RECOVERY, CANDIDATE)
        },
        "figures": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (FIGURE_ROOT.with_suffix(".pdf"), FIGURE_ROOT.with_suffix(".png"))
        },
        "betty_job_ids": ["7824764", "7824777", "7824787"],
    }
    FREEZE.parent.mkdir(parents=True, exist_ok=True)
    FREEZE.write_text(json.dumps(freeze, indent=2) + "\n")
    print(json.dumps(freeze, indent=2))


if __name__ == "__main__":
    main()
