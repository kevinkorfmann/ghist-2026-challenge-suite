from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-summaries", type=Path, nargs="+", required=True)
    parser.add_argument("--candidate-summary", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    pair_summaries = [json.loads(path.read_text()) for path in args.pair_summaries]
    pair_summaries.sort(key=lambda item: item["minimum_cm"])
    candidate = json.loads(args.candidate_summary.read_text())
    thresholds = [item["minimum_cm"] for item in pair_summaries]
    colors = {
        "accuracy": "#0072B2",
        "balanced_accuracy": "#009E73",
        "meioses_accuracy": "#D55E00",
    }
    labels = {
        "accuracy": "Exact relation",
        "balanced_accuracy": "Balanced exact",
        "meioses_accuracy": "Meiosis degree",
    }

    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    axis = axes[0]
    for metric in colors:
        values = [item["fresh_evaluation"][metric] for item in pair_summaries]
        axis.plot(
            thresholds,
            values,
            marker="o",
            linewidth=2,
            color=colors[metric],
            label=labels[metric],
        )
    axis.axvline(5, color="#555555", linestyle=":", linewidth=1)
    axis.set_xticks(thresholds)
    axis.set_ylim(0.35, 0.72)
    axis.set_xlabel("Minimum IBD segment (cM)")
    axis.set_ylabel("Untouched held-out accuracy")
    axis.set_title("A  Pair-model threshold selection")
    axis.legend(frameon=False, fontsize=8)
    axis.spines[["top", "right"]].set_visible(False)

    axis = axes[1]
    counts = candidate["label_counts"]
    ordered = sorted(counts, key=counts.get)
    names = [label.replace("full_", "").replace("_", " ") for label in ordered]
    values = [counts[label] for label in ordered]
    axis.barh(np.arange(len(ordered)), values, color="#56B4E9")
    axis.set_yticks(np.arange(len(ordered)), names, fontsize=7.5)
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
    figure.savefig(args.out_dir / "task10_validation.pdf", bbox_inches="tight")
    figure.savefig(args.out_dir / "task10_validation.png", dpi=220, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
