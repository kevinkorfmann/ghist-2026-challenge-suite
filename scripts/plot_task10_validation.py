from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--holdout-summaries", type=Path, nargs="+", required=True)
    parser.add_argument("--candidate-summary", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    holdouts = [json.loads(path.read_text()) for path in args.holdout_summaries]
    runs = [run for summary in holdouts for run in summary["runs"]]
    candidate = json.loads(args.candidate_summary.read_text())
    exact = 100 * np.asarray([run["exact_difference"] for run in runs])
    degree = 100 * np.asarray([run["meioses_difference"] for run in runs])

    plt.rcParams.update(
        {"font.family": "sans-serif", "font.size": 9, "axes.titleweight": "normal"}
    )

    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    axis = axes[0]
    positions = np.arange(1, len(runs) + 1)
    axis.axhline(0, color="#666666", linewidth=1, linestyle=":")
    axis.plot(positions, exact, marker="o", color="#0072B2", linewidth=1.8, label="Exact relation")
    axis.plot(positions, degree, marker="s", color="#D55E00", linewidth=1.8, label="Meiosis degree")
    axis.set_xticks(positions)
    axis.set_xlabel("Untouched simulation replicate")
    axis.set_ylabel("Paired accuracy gain (percentage points)")
    axis.set_title("A  Independent holdout gains")
    axis.legend(frameon=False, fontsize=8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.text(
        0.98,
        0.04,
        f"Mean exact gain: {exact.mean():+.2f} points\n{int(np.sum(exact > 0))}/{len(exact)} replicates improved",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="#333333",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
    )

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
