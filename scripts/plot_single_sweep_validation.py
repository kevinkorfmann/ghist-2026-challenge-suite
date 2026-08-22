from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


BLUE = "#2B6CB0"
ORANGE = "#DD6B20"
GREEN = "#2F855A"
PURPLE = "#805AD5"
RED = "#C53030"
GRAY = "#52606D"
PALE = "#E8F1F8"


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def ensemble_error(summary: dict) -> np.ndarray:
    records = summary["records"]
    truth = np.asarray([record["truth"]["position"] for record in records], dtype=float)
    estimate = np.asarray(
        [
            np.median(
                [
                    record["estimates"]["largest_gap_midpoint"],
                    record["estimates"]["composite"],
                    record["estimates"]["flank_ld_100kb"],
                ]
            )
            for record in records
        ]
    )
    indices = np.asarray([record["seed"] - 75_222_000 for record in records])
    selection = np.isin((indices * 17 + 3) % 5, [0, 1])
    return np.abs(estimate[~selection] - truth[~selection])


def ecdf(errors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(errors) / 1_000.0
    y = np.arange(1, x.size + 1) / x.size
    return np.r_[0, x], np.r_[0, y]


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.12, 1.04, label, transform=axis.transAxes, fontsize=12, va="top")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan", type=Path, required=True)
    parser.add_argument("--localization", type=Path, required=True)
    parser.add_argument("--flank-ld", type=Path, required=True)
    parser.add_argument("--raisd", type=Path, required=True)
    parser.add_argument("--generic-recovery", type=Path, required=True)
    parser.add_argument("--chromosome-recovery", type=Path, required=True)
    parser.add_argument("--scale-stress", type=Path, required=True)
    parser.add_argument("--interval", type=Path, required=True)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--png", type=Path, required=True)
    args = parser.parse_args()

    scan = np.genfromtxt(args.scan, delimiter=",", names=True)
    ld = np.genfromtxt(args.flank_ld, delimiter=",", names=True)
    localization = load(args.localization)
    raisd = load(args.raisd)
    generic = load(args.generic_recovery)
    chromosome = load(args.chromosome_recovery)
    scale = load(args.scale_stress)
    interval = load(args.interval)

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.weight": "normal",
            "axes.titleweight": "normal",
            "axes.labelweight": "normal",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(3, 2, figsize=(11.6, 11.4))
    fig.subplots_adjust(left=0.09, right=0.98, top=0.94, bottom=0.07, hspace=0.42, wspace=0.27)
    start, end = interval["interval_selection"]["interval"]
    center = interval["interval_selection"]["center"]

    # a: chromosome-wide scan.
    ax = axes[0, 0]
    valid = np.isfinite(scan["composite_score"])
    ax.plot(scan["center"][valid] / 1e6, scan["composite_score"][valid], color=BLUE, lw=0.8)
    ax.axvspan(start / 1e6, end / 1e6, color=ORANGE, alpha=0.25, lw=0)
    ax.set(xlabel="Chromosome 7 position (Mb)", ylabel="Composite sweep score", title="One dominant chromosome-wide signal")
    ax.set_xlim(0, np.nanmax(scan["end"]) / 1e6)
    panel_label(ax, "a")

    # b: local signatures on directly interpretable scales.
    ax = axes[0, 1]
    zoom = (scan["center"] >= 37_700_000) & (scan["center"] <= 38_900_000)
    background = np.nanmedian(scan["pi"])
    pi_depletion = np.clip(1.0 - scan["pi"] / background, -0.15, 1.0)
    negative_d = np.clip(-scan["tajima_d"] / 3.0, -0.15, 1.0)
    x = scan["center"][zoom] / 1e6
    ax.plot(x, pi_depletion[zoom], color=BLUE, lw=1.5, label="Diversity depletion")
    ax.plot(x, scan["h12"][zoom], color=ORANGE, lw=1.5, label="H12")
    ax.plot(x, negative_d[zoom], color=GREEN, lw=1.5, label="−Tajima's D / 3")
    ax.axvspan(start / 1e6, end / 1e6, color=GRAY, alpha=0.12, lw=0)
    ax.set(xlabel="Position (Mb)", ylabel="Scaled signature", title="Independent signatures converge locally")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    panel_label(ax, "b")

    # c: fine estimators and the scoring-aligned interval.
    ax = axes[1, 0]
    observed = interval["observed_estimates"]
    labels = ["Split-flank linkage", "RAiSD", "Gap midpoint", "Composite", "Locked center"]
    values = [
        observed["flank_ld_100kb"],
        observed["raisd"],
        observed["largest_gap_midpoint"],
        observed["composite_100kb"],
        center,
    ]
    colors = [GREEN, PURPLE, BLUE, ORANGE, RED]
    ax.axvspan(start / 1e6, end / 1e6, color=PALE, lw=0)
    for index, (value, color) in enumerate(zip(values, colors, strict=True)):
        ax.scatter(value / 1e6, index, color=color, s=42, edgecolor="white", linewidth=0.5, zorder=3)
    ax.axvline(center / 1e6, color=RED, lw=1.0, ls="--")
    ax.set_yticks(range(len(labels)), labels)
    ax.set(xlabel="Position (Mb)", title="Fine localization and candidate interval")
    ax.set_xlim(38.20, 38.30)
    panel_label(ax, "c")

    # d: untouched localization recovery.
    ax = axes[1, 1]
    for errors, label, color in [
        (ensemble_error(generic), "Varying-map background", BLUE),
        (ensemble_error(chromosome), "Chromosome-7 map", ORANGE),
    ]:
        ex, ey = ecdf(errors)
        ax.step(ex, ey, where="post", color=color, lw=1.8, label=f"{label} (n={errors.size})")
    ax.axvline(interval["interval_selection"]["final_half_width"] / 1_000, color=RED, ls="--", lw=1.1, label="35 kb candidate half-width")
    ax.set(xlabel="Absolute localization error (kb)", ylabel="Cumulative fraction", title="Untouched recovery after method lock")
    ax.set_xlim(0, 150)
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    panel_label(ax, "d")

    # e: exact and accelerated forward simulations.
    ax = axes[2, 0]
    scales = np.asarray(sorted(int(key) for key in scale["by_scale"]))
    median = np.asarray([scale["by_scale"][str(key)]["median_all_three"]["median_absolute_error"] for key in scales]) / 1_000
    q90 = np.asarray([scale["by_scale"][str(key)]["median_all_three"]["quantiles"]["0.9"] for key in scales]) / 1_000
    ax.plot(scales, median, marker="o", color=BLUE, lw=1.5, label="Median")
    ax.plot(scales, q90, marker="s", color=ORANGE, lw=1.5, label="90th percentile")
    ax.axhline(100, color=GRAY, lw=0.9, ls=":")
    ax.set_xticks(scales)
    ax.set(xlabel="Forward-simulation scale factor", ylabel="Absolute error (kb)", title="Localization survives simulation scaling")
    ax.legend(frameon=False, fontsize=8)
    panel_label(ax, "e")

    # f: interval length is chosen using the published challenge objective.
    ax = axes[2, 1]
    curve = interval["interval_selection"]["selection_scores"]
    widths = np.asarray(sorted(int(key) for key in curve)) / 1_000
    scores = np.asarray([curve[str(int(width * 1_000))] for width in widths])
    ax.plot(widths, scores, color=BLUE, lw=1.6, marker="o", ms=3)
    ax.axvline(interval["interval_selection"]["selection_optimal_half_width"] / 1_000, color=GREEN, lw=1.1, ls="--", label="Calibration optimum")
    ax.axvline(interval["interval_selection"]["final_half_width"] / 1_000, color=RED, lw=1.1, ls="--", label="Consensus-guarded width")
    ax.set(xlabel="Interval half-width (kb)", ylabel="Mean challenge score", title="Scoring-aligned interval selection")
    ax.set_xlim(0, 205)
    ax.legend(frameon=False, fontsize=8)
    panel_label(ax, "f")

    fig.suptitle("Challenge 5 | single hard-sweep localization", fontsize=16, y=0.985, fontweight="normal")
    fig.text(0.5, 0.018, "Testing chromosome 7 only • final chromosome 9 remains sealed • no submission", ha="center", color=GRAY, fontsize=9)
    args.pdf.parent.mkdir(parents=True, exist_ok=True)
    args.png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.pdf, dpi=300)
    fig.savefig(args.png, dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
