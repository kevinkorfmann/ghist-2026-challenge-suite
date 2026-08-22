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
COLORS = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00")
MARKERS = ("o", "s", "^", "D", "v")


def labelled_path(value: str) -> tuple[str, Path]:
    label, separator, path = value.partition("=")
    if not separator or not label or not path:
        raise argparse.ArgumentTypeError("expected LABEL=PATH")
    return label, Path(path)


def pairwise_selected_difference(source: dict) -> tuple[str, np.ndarray, float, list[float]]:
    selected, contender = source["ranking"][:2]
    direct = f"{selected}_minus_{contender}"
    reverse = f"{contender}_minus_{selected}"
    if direct in source["paired_comparisons"]:
        record = source["paired_comparisons"][direct]
        sign = 1.0
    elif reverse in source["paired_comparisons"]:
        record = source["paired_comparisons"][reverse]
        sign = -1.0
    else:
        raise KeyError(f"missing paired comparison for {selected!r} and {contender!r}")
    differences = sign * np.asarray(
        list(record["replicate_differences"].values()), dtype=np.float64
    )
    mean = sign * float(record["mean_score_difference"])
    interval = sorted(
        sign * np.asarray(record["mean_difference_bootstrap_interval_95"], dtype=float)
    )
    return f"{selected} − {contender}", differences, mean, interval


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison", action="append", type=labelled_path, required=True)
    parser.add_argument("--out-prefix", type=Path, required=True)
    args = parser.parse_args()
    comparisons = [
        (label, json.loads(path.read_text()), str(path))
        for label, path in args.comparison
    ]
    if any(not source["selection_valid"] for _, source, _ in comparisons):
        raise ValueError("all resolution comparisons must pass boundary checks")

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
    figure, axes = plt.subplots(1, 3, figsize=(12.2, 3.45))

    # a: individual fixed-center scores, retaining their matched-replicate identity.
    label, fixed, _ = comparisons[0]
    resolutions = list(fixed["summaries"])
    score_vectors = [
        np.asarray(list(fixed["summaries"][name]["replicate_scores"].values()))
        for name in resolutions
    ]
    boxes = axes[0].boxplot(
        score_vectors,
        positions=np.arange(len(resolutions)),
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 1.2},
        whiskerprops={"color": "#555555", "linewidth": 0.8},
        capprops={"color": "#555555", "linewidth": 0.8},
    )
    common_replicates = set.intersection(
        *(
            set(fixed["summaries"][name]["replicate_scores"])
            for name in resolutions
        )
    )
    for replicate in sorted(common_replicates, key=int):
        axes[0].plot(
            np.arange(len(resolutions)),
            [
                fixed["summaries"][name]["replicate_scores"][replicate]
                for name in resolutions
            ],
            color="#A7A7A7",
            linewidth=0.45,
            alpha=0.48,
            zorder=1,
        )
    for index, (box, values) in enumerate(zip(boxes["boxes"], score_vectors)):
        box.set(facecolor=COLORS[index % len(COLORS)], alpha=0.42, edgecolor="#333333")
        offsets = np.linspace(-0.16, 0.16, len(values))
        axes[0].scatter(
            index + offsets,
            values,
            s=12,
            marker=MARKERS[index % len(MARKERS)],
            color=COLORS[index % len(COLORS)],
            edgecolor="white",
            linewidth=0.35,
            zorder=3,
        )
    axes[0].set_xticks(np.arange(len(resolutions)), resolutions)
    axes[0].set_ylabel("Competition RRMSE")
    axes[0].set_title(f"Matched chromosomes: {label}")

    # b: score-aligned paired contrasts with bootstrap intervals for the mean.
    contrast_labels = []
    for index, (comparison_label, source, _) in enumerate(comparisons):
        contrast, differences, mean, interval = pairwise_selected_difference(source)
        contrast_labels.append(f"{comparison_label}\n{contrast}\nn={len(differences)}")
        offsets = np.linspace(-0.13, 0.13, len(differences))
        axes[1].scatter(
            index + offsets,
            differences,
            s=11,
            marker=MARKERS[index % len(MARKERS)],
            color=COLORS[index % len(COLORS)],
            alpha=0.72,
            linewidth=0,
        )
        axes[1].errorbar(
            index,
            mean,
            yerr=[[mean - interval[0]], [interval[1] - mean]],
            marker="D",
            markersize=5,
            color="#222222",
            capsize=3,
            linewidth=1.1,
            zorder=4,
        )
    axes[1].axhline(0, color="#777777", linewidth=0.8, linestyle="--")
    axes[1].set_xticks(np.arange(len(contrast_labels)), contrast_labels)
    axes[1].set_ylabel("Selected − contender RRMSE")
    axes[1].set_title("Paired risk difference (mean and 95% interval)")

    # c: parameter-resolved natural-scale error for the broadest comparison.
    _, broad, _ = comparisons[-1]
    heatmap_labels = list(broad["summaries"])
    matrix = np.asarray(
        [
            list(broad["summaries"][name]["median_absolute_relative_error"].values())
            for name in heatmap_labels
        ]
    )
    image = axes[2].imshow(matrix * 100.0, cmap="cividis", aspect="auto")
    axes[2].set_xticks(np.arange(len(PARAMETER_LABELS)), PARAMETER_LABELS, rotation=38, ha="right")
    axes[2].set_yticks(np.arange(len(heatmap_labels)), heatmap_labels)
    axes[2].set_title("Median absolute parameter error (%)")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix[row, column] * 100.0
            color = "white" if value < 0.55 * np.nanmax(matrix * 100.0) else "black"
            axes[2].text(column, row, f"{value:.1f}", ha="center", va="center", color=color, fontsize=7)
    colorbar = figure.colorbar(image, ax=axes[2], fraction=0.045, pad=0.03)
    colorbar.set_label("Error (%)")

    for panel, axis in zip(("a", "b", "c"), axes):
        axis.text(-0.14, 1.04, panel, transform=axis.transAxes, fontsize=10, va="bottom")
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.55)
    axes[2].grid(False)
    figure.tight_layout(w_pad=2.0)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out_prefix.with_suffix(".png"), dpi=360, bbox_inches="tight")
    figure.savefig(args.out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
