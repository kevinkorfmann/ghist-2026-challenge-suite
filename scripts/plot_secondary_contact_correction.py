from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


PARAMETER_LABELS = (
    "Mainland size",
    "Island size",
    "Split time",
    "Contact time",
    "Migration rate",
)
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
MAGENTA = "#CC79A7"
CONDITION_LABELS = {
    "nominal": "Nominal",
    "recombination": "Recombination",
    "mutation": "Mutation",
    "local_size": "Local size",
    "migration_asymmetry": "Asymmetric\nmigration",
    "polarization": "Polarization",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--stress-transfer", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, required=True)
    args = parser.parse_args()
    calibration = json.loads(args.calibration.read_text())
    transfer = json.loads(args.stress_transfer.read_text())
    evaluation = calibration["untouched_evaluation"]
    raw_scores = np.asarray(evaluation["raw"]["replicate_scores"], dtype=float)
    corrected_scores = np.asarray(
        evaluation["corrected"]["replicate_scores"], dtype=float
    )
    if raw_scores.shape != corrected_scores.shape:
        raise ValueError("raw and corrected evaluation scores must be paired")

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "font.weight": "normal",
            "axes.titleweight": "normal",
            "axes.labelweight": "normal",
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure, axes = plt.subplots(1, 3, figsize=(12.0, 3.45))

    # a: untouched evaluation scores preserve parameter-point grouping.
    for raw, corrected in zip(raw_scores, corrected_scores):
        axes[0].plot((0, 1), (raw, corrected), color="#A8A8A8", alpha=0.65, linewidth=0.65)
    axes[0].scatter(
        np.zeros(len(raw_scores)), raw_scores, color=BLUE, marker="o", s=22, edgecolor="white", linewidth=0.4, zorder=3
    )
    axes[0].scatter(
        np.ones(len(corrected_scores)), corrected_scores, color=ORANGE, marker="s", s=22, edgecolor="white", linewidth=0.4, zorder=3
    )
    axes[0].set_xticks((0, 1), ("Raw", "Corrected"))
    axes[0].set_xlim(-0.35, 1.35)
    axes[0].set_ylabel("Competition RRMSE")
    axes[0].set_title(f"Untouched varying-demography evaluation (n={len(raw_scores)})")

    # b: natural-scale error for each scored parameter.
    raw_errors = np.asarray(
        list(evaluation["raw"]["median_absolute_relative_error"].values())
    )
    corrected_errors = np.asarray(
        list(evaluation["corrected"]["median_absolute_relative_error"].values())
    )
    positions = np.arange(len(PARAMETER_LABELS))
    width = 0.36
    axes[1].bar(
        positions - width / 2,
        raw_errors * 100,
        width,
        color=BLUE,
        alpha=0.72,
        label="Raw",
    )
    axes[1].bar(
        positions + width / 2,
        corrected_errors * 100,
        width,
        color=ORANGE,
        alpha=0.78,
        label="Corrected",
    )
    axes[1].set_xticks(positions, PARAMETER_LABELS, rotation=35, ha="right")
    axes[1].set_ylabel("Median absolute error (%)")
    axes[1].set_title("Untouched parameter error")
    axes[1].legend(frameon=False, loc="upper right")

    # c: score change under nominal and deliberately perturbed generators.
    conditions = transfer["conditions"]
    labels = list(conditions)
    for index, label in enumerate(labels):
        differences = np.asarray(
            conditions[label]["corrected_minus_raw"]["replicate_differences"],
            dtype=float,
        )
        improves = differences.mean() < 0
        offsets = np.linspace(-0.13, 0.13, len(differences))
        axes[2].scatter(
            index + offsets,
            differences,
            color=GREEN if improves else MAGENTA,
            marker="o" if improves else "^",
            s=14,
            alpha=0.68,
            linewidth=0,
        )
        axes[2].scatter(
            index,
            differences.mean(),
            color="#222222",
            marker="D",
            s=30,
            zorder=4,
        )
    axes[2].axhline(0, color="#777777", linewidth=0.8, linestyle="--")
    axes[2].set_xticks(
        np.arange(len(labels)),
        [CONDITION_LABELS.get(label, label.replace("_", " ").title()) for label in labels],
        rotation=24,
        ha="right",
    )
    axes[2].set_ylabel("Corrected − raw RRMSE")
    axes[2].set_title("Transfer under generator perturbations")

    for panel, axis in zip(("a", "b", "c"), axes):
        axis.text(-0.14, 1.04, panel, transform=axis.transAxes, fontsize=10, va="bottom")
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.55)
    figure.tight_layout(w_pad=2.0)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out_prefix.with_suffix(".png"), dpi=360, bbox_inches="tight")
    figure.savefig(args.out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
