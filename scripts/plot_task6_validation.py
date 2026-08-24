from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--localization", type=Path, required=True)
    parser.add_argument("--raisd", type=Path, required=True)
    parser.add_argument("--transfer", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    localization = json.loads(args.localization.read_text())
    raisd = json.loads(args.raisd.read_text())
    transfer = json.loads(args.transfer.read_text())
    center = float(localization["largest_variant_gap"]["midpoint"])
    half_width = int(transfer["selected_half_width"])
    points = [
        ("Largest variant gap", center),
        ("RAiSD maximum", float(raisd["maximum"]["position"])),
        (
            "Flank LD, 50 kb",
            float(localization["flank_ld_scans"]["50000"]["maximum"]["candidate_position"]),
        ),
        (
            "Flank LD, 100 kb",
            float(localization["flank_ld_scans"]["100000"]["maximum"]["candidate_position"]),
        ),
    ]

    figure, axes = plt.subplots(1, 2, figsize=(11.2, 4.5))
    axis = axes[0]
    axis.axvspan(
        (center - half_width) / 1e6,
        (center + half_width) / 1e6,
        color="#F0E442",
        alpha=0.35,
        label="Frozen interval",
    )
    for index, (label, position) in enumerate(points):
        axis.scatter(position / 1e6, index, s=55, color="#0072B2", zorder=3)
    axis.set_yticks(np.arange(len(points)), [label for label, _ in points])
    axis.set_xlim(16.50, 16.64)
    axis.set_xlabel("Chromosome 7 position (Mb)")
    axis.set_title("A  Independent observed localizers")
    axis.legend(frameon=False, fontsize=8, loc="lower right")
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="y", length=0)

    all_metrics = transfer["all_regime_untouched_evaluation"]["metrics"]
    signal = transfer["signal_matched_untouched_evaluation"]
    frequency = transfer["high_frequency_untouched_evaluation"]
    names = ["All regimes", "Signal matched", "High frequency"]
    medians = [
        all_metrics["median_absolute_error"] / 1000,
        signal["median_absolute_error"] / 1000,
        frequency["median_absolute_error"] / 1000,
    ]
    q90 = [
        all_metrics["quantiles"]["0.9"] / 1000,
        signal["q90_absolute_error"] / 1000,
        frequency["q90_absolute_error"] / 1000,
    ]
    x = np.arange(len(names))
    width = 0.34
    axis = axes[1]
    axis.bar(x - width / 2, medians, width, color="#009E73", label="Median")
    axis.bar(x + width / 2, q90, width, color="#D55E00", label="90th percentile")
    axis.axhline(250, color="#555555", linestyle=":", linewidth=1, label="250 kb gate")
    axis.set_xticks(x, names, rotation=15, ha="right")
    axis.set_ylabel("Untouched localization error (kb)")
    axis.set_title("B  Recovery across sweep regimes")
    axis.legend(frameon=False, fontsize=8)
    axis.spines[["top", "right"]].set_visible(False)

    figure.tight_layout()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out_dir / "task6_validation.pdf", bbox_inches="tight")
    figure.savefig(args.out_dir / "task6_validation.png", dpi=220, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
