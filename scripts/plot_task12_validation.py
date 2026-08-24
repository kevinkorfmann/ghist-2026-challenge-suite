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
    entries = {
        (entry["quantile"], entry["trim_shortest_shape"]): entry
        for entry in validation["aggregates"]
    }
    quantiles = (0.25, 0.5, 0.75)
    x = np.arange(len(quantiles))
    width = 0.18
    colors = {"exact": "#0072B2", "meioses": "#D55E00"}

    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    axis = axes[0]
    for offset, (metric, trim, hatch) in enumerate(
        (
            ("exact", False, ""),
            ("exact", True, "///"),
            ("meioses", False, ""),
            ("meioses", True, "///"),
        )
    ):
        key = f"mean_{metric}_accuracy"
        values = [entries[(quantile, trim)][key] for quantile in quantiles]
        positions = x + (offset - 1.5) * width
        axis.bar(
            positions,
            values,
            width,
            color=colors[metric],
            alpha=0.9 if trim else 0.48,
            hatch=hatch,
            edgecolor="white" if not trim else colors[metric],
            linewidth=0.7,
            label=f"{metric.capitalize()}, {'trimmed' if trim else 'untrimmed'}",
        )
    axis.axhline(
        validation["baseline_exact_accuracy"],
        color=colors["exact"],
        linestyle=":",
        linewidth=1.2,
    )
    axis.axhline(
        validation["baseline_meioses_accuracy"],
        color=colors["meioses"],
        linestyle=":",
        linewidth=1.2,
    )
    axis.set_xticks(x, ["Q25", "Q50", "Q75"])
    axis.set_ylim(0.3, 0.65)
    axis.set_ylabel("Held-out accuracy")
    axis.set_xlabel("Background correction")
    axis.set_title("A  Structured-background stress test")
    axis.legend(frameon=False, fontsize=8, ncol=2, loc="lower left")
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
