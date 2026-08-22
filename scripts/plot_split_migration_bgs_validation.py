from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "experiments" / "split_migration_bgs"
BOOTSTRAP = ROOT / "artifacts" / "task4_split_migration_bgs" / "bootstrap_1mb"
OUTPUT = ROOT / "figures" / "task4_split_migration_bgs"
PARAMETERS = (
    "AFR_population_size",
    "EUR_population_size",
    "generations_since_split",
    "migration_rate",
)
LABELS = ("AFR size", "EUR size", "Split time", "Migration")
COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7")
CANDIDATE = np.asarray((16_000.0, 24_000.0, 4_500.0, 0.00015))


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def vector(row: dict) -> np.ndarray:
    return np.asarray([row[key] for key in PARAMETERS], dtype=np.float64)


def relative_percent(values: np.ndarray, reference: np.ndarray) -> np.ndarray:
    return 100 * (values / reference - 1)


def panel_label(axis, label: str) -> None:
    axis.text(
        -0.14,
        1.06,
        label,
        transform=axis.transAxes,
        fontsize=12,
        fontweight="normal",
        va="top",
    )


def main() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
        }
    )
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 7.4), constrained_layout=True)

    distance_files = (
        ("Coding", "distance_coding_p50_40.json"),
        ("1 bp–1 kb", "distance_near_p50_40.json"),
        ("1–100 kb", "baseline_p50_40.json"),
        ("All noncoding", "distance_all_noncoding_p50_40.json"),
    )
    distance_values = np.asarray(
        [vector(load(EVIDENCE / filename)["fit"]) for _, filename in distance_files]
    )
    axis = axes[0, 0]
    x = np.arange(len(distance_files))
    for index, (label, color) in enumerate(zip(LABELS, COLORS)):
        axis.plot(
            x,
            relative_percent(distance_values[:, index], CANDIDATE[index]),
            marker="o",
            linewidth=1.5,
            markersize=4.5,
            color=color,
            label=label,
        )
    axis.axhline(0, color="#444444", linewidth=0.9, linestyle="--")
    axis.set_xticks(x, [label for label, _ in distance_files], rotation=18, ha="right")
    axis.set_ylabel("Difference from candidate (%)")
    axis.set_title("Distance from coding sequence")
    axis.legend(ncol=2, frameon=False, loc="lower left")
    panel_label(axis, "a")

    baseline = load(EVIDENCE / "baseline_summary.json")["by_projection"]
    projection_keys = ("50x40", "70x50", "94x68")
    projection_values = np.asarray([vector(baseline[key]) for key in projection_keys])
    axis = axes[0, 1]
    x = np.arange(len(projection_keys))
    for index, (label, color) in enumerate(zip(LABELS, COLORS)):
        axis.plot(
            x,
            relative_percent(projection_values[:, index], CANDIDATE[index]),
            marker="o",
            linewidth=1.5,
            markersize=4.5,
            color=color,
            label=label,
        )
    axis.axhline(0, color="#444444", linewidth=0.9, linestyle="--")
    axis.set_xticks(x, ("50 × 40", "70 × 50", "94 × 68"))
    axis.set_ylabel("Difference from candidate (%)")
    axis.set_title("Projection stability (1–100 kb)")
    panel_label(axis, "b")

    recovery_files = [
        EVIDENCE / "selection_candidate_aggregate.json",
        EVIDENCE / "selection_near_observed_aggregate.json",
        *(EVIDENCE / "selection_grid_aggregate" / f"point_{i}.json" for i in range(4)),
    ]
    recoveries = [load(path) for path in recovery_files]
    recovery_errors = np.asarray(
        [relative_percent(vector(row["fit"]), vector(row["truth"])) for row in recoveries]
    )
    axis = axes[1, 0]
    x = np.arange(len(recoveries))
    for index, (label, color) in enumerate(zip(LABELS, COLORS)):
        axis.plot(
            x,
            recovery_errors[:, index],
            marker="o",
            linewidth=1.2,
            markersize=4.2,
            color=color,
            label=label,
        )
    axis.axhline(0, color="#444444", linewidth=0.9, linestyle="--")
    axis.set_xticks(x, ("Larger", "Near", "Grid 1", "Grid 2", "Grid 3", "Grid 4"))
    axis.set_ylabel("Recovery error (%)")
    axis.set_title("Six held-out truth settings")
    axis.set_ylim(-16, 16)
    panel_label(axis, "c")

    bootstrap_files = sorted(
        BOOTSTRAP.glob("task_*.json"),
        key=lambda path: int(path.stem.rsplit("_", 1)[1]),
    )
    if len(bootstrap_files) != 64:
        raise ValueError(f"expected 64 bootstrap results, found {len(bootstrap_files)}")
    bootstrap_values = np.asarray([vector(load(path)["fit"]) for path in bootstrap_files])
    bootstrap_relative = relative_percent(bootstrap_values, CANDIDATE)
    axis = axes[1, 1]
    rng = np.random.default_rng(640750)
    for index, (label, color) in enumerate(zip(LABELS, COLORS)):
        jitter = rng.uniform(-0.12, 0.12, len(bootstrap_relative))
        axis.scatter(
            index + jitter,
            bootstrap_relative[:, index],
            s=12,
            alpha=0.42,
            color=color,
            edgecolors="none",
        )
        low, median, high = np.quantile(
            bootstrap_relative[:, index], (0.025, 0.5, 0.975)
        )
        axis.plot([index, index], [low, high], color="#222222", linewidth=1.3)
        axis.plot(
            [index - 0.13, index + 0.13],
            [median, median],
            color="#222222",
            linewidth=2.0,
        )
    axis.axhline(0, color="#444444", linewidth=0.9, linestyle="--")
    axis.set_xticks(np.arange(4), LABELS, rotation=18, ha="right")
    axis.set_ylabel("Difference from candidate (%)")
    axis.set_title("64 chromosome-stratified 1 Mb bootstraps")
    panel_label(axis, "d")

    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.55, alpha=0.7)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT / "task4_validation.pdf", bbox_inches="tight")
    figure.savefig(OUTPUT / "task4_validation.png", dpi=400, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
