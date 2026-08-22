from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
GRAY = "#6B6B6B"
BLACK = "#111111"


def style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
            "font.size": 8,
            "font.weight": "normal",
            "axes.labelsize": 8.5,
            "axes.labelweight": "normal",
            "axes.titleweight": "normal",
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.7,
            "lines.linewidth": 1.3,
        }
    )


def clean_axis(axis) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(width=0.7, length=3)
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.45, alpha=0.65)
    axis.set_axisbelow(True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--distance-fit", type=Path, required=True)
    parser.add_argument("--middle-summary", type=Path, required=True)
    parser.add_argument("--audit-summary", type=Path, required=True)
    parser.add_argument("--bootstrap-dir", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, required=True)
    args = parser.parse_args()
    style()

    distance = json.loads(args.distance_fit.read_text())["results"]
    middle = json.loads(args.middle_summary.read_text())["strata"][
        "noncoding_1kb_100kb"
    ]["by_projection"]
    exposure = json.loads(args.audit_summary.read_text())[
        "physical_bp_by_distance_stratum"
    ]
    bootstrap_rows = []
    for path in sorted(args.bootstrap_dir.glob("task_*.json")):
        bootstrap_rows.extend(json.loads(path.read_text())["results"])
    if len(bootstrap_rows) != 128:
        raise ValueError(f"expected 128 bootstrap replicates, found {len(bootstrap_rows)}")

    keys = [
        "distance_coding_overlap",
        "distance_noncoding_1bp_1kb",
        "distance_noncoding_1kb_10kb",
        "distance_noncoding_10kb_100kb",
        "distance_noncoding_gt100kb",
    ]
    labels = ["Coding", "1 bp–1 kb", "1–10 kb", "10–100 kb", ">100 kb"]
    exposure_keys = [
        "coding_overlap",
        "noncoding_1bp_1kb",
        "noncoding_1kb_10kb",
        "noncoding_10kb_100kb",
        "noncoding_gt100kb",
    ]
    colors = [GRAY, ORANGE, GREEN, GREEN, BLUE]
    generations = np.asarray([distance[key]["generations"] for key in keys])
    scaling = np.asarray([distance[key]["post_growth_scaling"] for key in keys])
    x = np.arange(len(keys))

    figure, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
    ax_a, ax_b, ax_c, ax_d = axes.flat
    for axis in axes.flat:
        clean_axis(axis)

    ax_a.scatter(x, generations, c=colors, s=34, edgecolor=BLACK, linewidth=0.45, zorder=3)
    ax_a.axhline(5000, color=BLACK, linestyle="--", linewidth=0.9, label="Provisional value")
    ax_a.set_ylabel("Growth start (generations)")
    ax_a.set_xticks(x, labels, rotation=24, ha="right")
    ax_a.set_ylim(2300, 7500)

    ax_b.scatter(x, scaling, c=colors, s=34, edgecolor=BLACK, linewidth=0.45, zorder=3)
    ax_b.axhline(2.0, color=BLACK, linestyle="--", linewidth=0.9)
    ax_b.set_ylabel("Post-growth size ratio")
    ax_b.set_xticks(x, labels, rotation=24, ha="right")
    ax_b.set_ylim(1.75, 4.45)
    for index, key in enumerate(exposure_keys):
        ax_b.text(
            index,
            4.30,
            f"{exposure[key] / 1e6:.1f} Mb",
            ha="center",
            va="top",
            fontsize=6.5,
            color=GRAY,
        )

    projections = np.asarray(sorted(int(key) for key in middle))
    middle_generation = np.asarray(
        [middle[str(projection)]["generations"] for projection in projections]
    )
    middle_scaling = np.asarray(
        [middle[str(projection)]["post_growth_scaling"] for projection in projections]
    )
    ax_c.plot(
        projections,
        100 * (middle_generation / 5000 - 1),
        color=BLUE,
        marker="o",
        label="Growth start",
    )
    ax_c.plot(
        projections,
        100 * (middle_scaling / 2 - 1),
        color=ORANGE,
        marker="s",
        linestyle="--",
        label="Size ratio",
    )
    ax_c.axhline(0, color=BLACK, linewidth=0.8)
    ax_c.set_xlabel("Projected haplotypes")
    ax_c.set_ylabel("Difference from provisional value (%)")
    ax_c.set_xticks(projections)
    ax_c.legend(frameon=False, loc="upper left")

    boot_generation = np.asarray([row["generations"] for row in bootstrap_rows]) / 5000
    boot_scaling = np.asarray(
        [row["post_growth_scaling"] for row in bootstrap_rows]
    ) / 2
    rng = np.random.default_rng(640742)
    for index, (values, color) in enumerate(
        ((boot_generation, BLUE), (boot_scaling, ORANGE))
    ):
        jitter = rng.uniform(-0.11, 0.11, len(values))
        ax_d.scatter(
            np.full(len(values), index) + jitter,
            values,
            color=color,
            alpha=0.24,
            s=9,
            linewidth=0,
        )
        low, median, high = np.quantile(values, (0.025, 0.5, 0.975))
        ax_d.plot([index, index], [low, high], color=BLACK, linewidth=1.1)
        ax_d.scatter([index], [median], color="white", edgecolor=BLACK, s=30, zorder=4)
    ax_d.axhline(1, color=BLACK, linestyle="--", linewidth=0.9)
    ax_d.set_xticks([0, 1], ["Growth start", "Size ratio"])
    ax_d.set_ylabel("Bootstrap estimate / provisional value")
    ax_d.set_xlim(-0.45, 1.45)
    ax_d.text(
        0.98,
        0.03,
        "n = 128; 5 Mb blocks",
        transform=ax_d.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.8,
        color=GRAY,
    )

    for label, axis in zip("ABCD", axes.flat):
        axis.text(
            -0.14,
            1.04,
            label,
            transform=axis.transAxes,
            ha="left",
            va="bottom",
            fontsize=9,
            fontweight="normal",
        )
    figure.subplots_adjust(left=0.10, right=0.985, bottom=0.13, top=0.97, wspace=0.34, hspace=0.44)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(args.out_prefix.with_suffix(".png"), dpi=400, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
