from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--candidate-summary", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    validation = json.loads(args.validation.read_text())
    candidate = json.loads(args.candidate_summary.read_text())
    colors = {"exact": "#0072B2", "meioses": "#D55E00"}

    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    axis = axes[0]
    positions = {
        ("exact", "baseline"): 0.0,
        ("exact", "candidate"): 1.0,
        ("meioses", "baseline"): 2.5,
        ("meioses", "candidate"): 3.5,
    }
    for metric in ("exact", "meioses"):
        key = f"{metric}_accuracy"
        baseline = np.asarray([row[key] for row in validation["baseline_runs"]])
        candidate_values = np.asarray(
            [row[key] for row in validation["candidate_runs"]]
        )
        left = positions[(metric, "baseline")]
        right = positions[(metric, "candidate")]
        for baseline_value, candidate_value in zip(
            baseline, candidate_values, strict=True
        ):
            axis.plot(
                (left, right),
                (baseline_value, candidate_value),
                color="#999999",
                alpha=0.32,
                linewidth=0.65,
                zorder=1,
            )
        axis.scatter(
            np.full(baseline.size, left),
            baseline,
            color=colors[metric],
            facecolors="white",
            linewidths=0.7,
            s=15,
            alpha=0.7,
            zorder=2,
        )
        axis.scatter(
            np.full(candidate_values.size, right),
            candidate_values,
            color=colors[metric],
            s=15,
            alpha=0.7,
            zorder=2,
        )
        axis.plot(
            (left - 0.17, left + 0.17),
            (baseline.mean(), baseline.mean()),
            color="#111111",
            linewidth=2.0,
            zorder=3,
        )
        axis.plot(
            (right - 0.17, right + 0.17),
            (candidate_values.mean(), candidate_values.mean()),
            color="#111111",
            linewidth=2.0,
            zorder=3,
        )
    axis.set_xticks(
        [0.0, 1.0, 2.5, 3.5],
        ["Previous", "Frozen", "Previous", "Frozen"],
    )
    axis.text(0.5, -0.14, "Exact label", ha="center", transform=axis.get_xaxis_transform())
    axis.text(3.0, -0.14, "Meiosis degree", ha="center", transform=axis.get_xaxis_transform())
    axis.set_ylim(0.28, 0.62)
    axis.set_ylabel("Held-out accuracy")
    axis.set_title("A  Empirical-background holdouts")
    axis.spines[["top", "right"]].set_visible(False)

    axis = axes[1]
    label_counts = candidate["label_counts"]
    ordered = sorted(label_counts, key=label_counts.get)
    labels = [label.replace("full_", "").replace("_", " ") for label in ordered]
    values = [label_counts[label] for label in ordered]
    axis.barh(np.arange(len(ordered)), values, color="#009E73")
    axis.set_yticks(np.arange(len(ordered)), labels, fontsize=7.5)
    axis.set_xlabel("Predicted pairs")
    axis.set_title("B  Frozen testing prediction")
    axis.spines[["top", "right"]].set_visible(False)
    axis.text(
        0.98,
        0.02,
        f"{candidate['pairs']:,} pairs; {candidate['individuals']} individuals",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="#444444",
    )

    figure.tight_layout()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out_dir / "task12_validation.pdf", bbox_inches="tight")
    figure.savefig(args.out_dir / "task12_validation.png", dpi=220, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
