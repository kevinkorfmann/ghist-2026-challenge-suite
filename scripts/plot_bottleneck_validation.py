from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    "nominal": "#0072B2",
    "rate": "#E69F00",
    "observation": "#009E73",
    "combined": "#CC79A7",
    "reference": "#4D4D4D",
}
SCENARIO_LABELS = {
    0: "Nominal",
    1: "Rate variation",
    2: "Observation / diversity",
    3: "Combined",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.weight": "normal",
            "axes.titleweight": "normal",
            "axes.labelweight": "normal",
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def inverse_targets(values: np.ndarray) -> np.ndarray:
    result = np.empty_like(values, dtype=np.float64)
    result[..., 0] = np.exp(values[..., 0])
    result[..., 1] = 1 / (1 + np.exp(-values[..., 1]))
    return result


def plot_validation(input_path: Path, output_stem: Path) -> None:
    configure_style()
    data = np.load(input_path)
    truth = inverse_targets(data["target_latent"])
    deep = inverse_targets(data["deep_mean_latent"])
    reference = inverse_targets(data["reference_latent"])
    scenarios = data["nuisance_scenario"].astype(int)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), constrained_layout=True)
    parameter_specs = (
        (0, "Decline time (generations)"),
        (1, "Post-decline fraction"),
    )
    for ax, (column, label) in zip(axes[0], parameter_specs):
        for scenario in sorted(np.unique(scenarios)):
            mask = scenarios == scenario
            ax.scatter(
                truth[mask, column],
                deep[mask, column],
                s=9,
                alpha=0.55,
                linewidths=0,
                color=list(COLORS.values())[scenario],
                label=SCENARIO_LABELS.get(scenario, f"Scenario {scenario}"),
            )
        limits = np.nanpercentile(
            np.concatenate((truth[:, column], deep[:, column])), [0.5, 99.5]
        )
        ax.plot(limits, limits, color="#262626", linewidth=0.8, linestyle="--")
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.set_xlabel(f"True {label.lower()}")
        ax.set_ylabel(f"Estimated {label.lower()}")
        if column == 0:
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.legend(frameon=False, loc="upper left")

    relative_time_error = np.abs(deep[:, 0] - truth[:, 0]) / truth[:, 0]
    scenario_ids = sorted(np.unique(scenarios))
    positions = np.arange(len(scenario_ids))
    ax = axes[1, 0]
    time_groups = [relative_time_error[scenarios == scenario] for scenario in scenario_ids]
    boxes = ax.boxplot(
        time_groups,
        positions=positions,
        widths=0.58,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#111111", "linewidth": 0.8},
        boxprops={"linewidth": 0.7},
        whiskerprops={"linewidth": 0.7},
        capprops={"linewidth": 0.7},
    )
    for patch, scenario in zip(boxes["boxes"], scenario_ids):
        patch.set_facecolor(list(COLORS.values())[scenario])
        patch.set_alpha(0.75)
    ax.set_xticks(positions, [SCENARIO_LABELS[i] for i in scenario_ids], rotation=24, ha="right")
    ax.set_ylabel("Absolute relative time error")

    ax = axes[1, 1]
    width = 0.34
    deep_rmse = []
    reference_rmse = []
    for scenario in scenario_ids:
        mask = scenarios == scenario
        deep_rmse.append(np.sqrt(np.mean(np.square(deep[mask, 1] - truth[mask, 1]))))
        reference_rmse.append(
            np.sqrt(np.mean(np.square(reference[mask, 1] - truth[mask, 1])))
        )
    ax.bar(positions - width / 2, deep_rmse, width, color="#0072B2", label="Context ensemble")
    ax.bar(
        positions + width / 2,
        reference_rmse,
        width,
        color=COLORS["reference"],
        hatch="///",
        label="Permutation-invariant reference",
    )
    ax.set_xticks(positions, [SCENARIO_LABELS[i] for i in scenario_ids], rotation=24, ha="right")
    ax.set_ylabel("Fraction RMSE")
    ax.legend(frameon=False)

    for label, ax in zip("abcd", axes.flat):
        ax.text(-0.13, 1.04, label, transform=ax.transAxes, fontsize=9, fontweight="normal")

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=400, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-stem", required=True, type=Path)
    args = parser.parse_args()
    plot_validation(args.input, args.output_stem)


if __name__ == "__main__":
    main()
