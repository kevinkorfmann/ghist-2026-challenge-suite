from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    "moments": "#0072B2",
    "dadi": "#E69F00",
    "joint": "#009E73",
    "ld": "#56B4E9",
    "learned": "#CC79A7",
    "score": "#D55E00",
    "diagnostic": "#6B6B6B",
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
            "axes.titlesize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def natural(latent: np.ndarray) -> np.ndarray:
    return np.column_stack((np.exp(latent[:, 0]), 1 / (1 + np.exp(-latent[:, 1]))))


def read_estimate(path: Path) -> tuple[float, float]:
    data = json.loads(path.read_text())
    return float(data["generations"]), float(data["post_decline_fraction"])


def strip_svg_metadata(path: Path) -> None:
    text = path.read_text()
    start = text.find("<metadata>")
    stop = text.find("</metadata>")
    if start >= 0 and stop >= start:
        path.write_text(text[:start] + text[stop + len("</metadata>") :])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--moments-validation", type=Path, required=True)
    parser.add_argument("--dadi-validation", type=Path, required=True)
    parser.add_argument("--learned-validation", type=Path, required=True)
    parser.add_argument("--split-record", type=Path, required=True)
    parser.add_argument("--moments", type=Path, required=True)
    parser.add_argument("--dadi", type=Path, required=True)
    parser.add_argument("--joint", type=Path, required=True)
    parser.add_argument("--ld-only", type=Path, required=True)
    parser.add_argument("--learned-observed", type=Path, required=True)
    parser.add_argument("--score-gate", type=Path, required=True)
    parser.add_argument("--diagnostic-observed", type=Path, required=True)
    parser.add_argument("--bootstrap", type=Path, required=True)
    parser.add_argument("--output-stem", type=Path, required=True)
    args = parser.parse_args()

    configure_style()
    moments_validation = np.load(args.moments_validation)
    dadi_validation = np.load(args.dadi_validation)
    learned_validation = np.load(args.learned_validation)
    split = json.loads(args.split_record.read_text())["split_cells"]
    moments_index = {int(cell): i for i, cell in enumerate(moments_validation["cell_id"])}
    dadi_index = {int(cell): i for i, cell in enumerate(dadi_validation["cell_id"])}
    learned_index = {int(cell): i for i, cell in enumerate(learned_validation["cell_id"])}
    cells = [
        int(cell)
        for cell in split["evaluation"]
        if int(cell) in moments_index and int(cell) in dadi_index and int(cell) in learned_index
    ]
    mi = np.asarray([moments_index[cell] for cell in cells])
    di = np.asarray([dadi_index[cell] for cell in cells])
    li = np.asarray([learned_index[cell] for cell in cells])
    truth = natural(moments_validation["target_latent"][mi])
    recovery = {
        "moments": natural(moments_validation["raw_latent"][mi]),
        "dadi": natural(dadi_validation["raw_latent"][di]),
        "learned": natural(learned_validation["summary_calibrated_latent"][li]),
    }

    score_gate = json.loads(args.score_gate.read_text())
    score_prediction = score_gate["observed"]["selected_candidate"]
    observed = {
        "moments": read_estimate(args.moments),
        "dadi": read_estimate(args.dadi),
        "joint": read_estimate(args.joint),
        "ld": read_estimate(args.ld_only),
        "learned": read_estimate(args.learned_observed),
        "score": (float(score_prediction[0]), float(score_prediction[1])),
        "diagnostic": read_estimate(args.diagnostic_observed),
    }
    bootstrap = json.loads(args.bootstrap.read_text())
    generation_interval = np.asarray(bootstrap["interval_95"][0], dtype=float)
    fraction_interval = np.asarray(bootstrap["interval_95"][1], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.55), constrained_layout=True)
    labels = {"moments": "moments SFS", "dadi": "dadi SFS", "learned": "Context ensemble"}
    markers = {"moments": "o", "dadi": "s", "learned": "^"}
    for column, ax in enumerate(axes[0]):
        for method in ("moments", "dadi", "learned"):
            ax.scatter(
                truth[:, column],
                recovery[method][:, column],
                s=18,
                marker=markers[method],
                facecolors="none" if method == "learned" else COLORS[method],
                edgecolors=COLORS[method],
                linewidths=0.7,
                alpha=0.78,
                label=labels[method],
            )
        values = np.concatenate([truth[:, column], *[value[:, column] for value in recovery.values()]])
        lower, upper = np.quantile(values, [0.01, 0.99])
        ax.plot([lower, upper], [lower, upper], color="#303030", linestyle="--", linewidth=0.8)
        ax.set_xlim(lower, upper)
        ax.set_ylim(lower, upper)
        if column == 0:
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("True decline time (generations)")
            ax.set_ylabel("Estimated decline time (generations)")
            ax.legend(frameon=False, loc="upper left")
        else:
            ax.set_xlabel("True post-decline fraction")
            ax.set_ylabel("Estimated post-decline fraction")

    method_order = ["moments", "dadi", "joint", "ld", "learned", "score", "diagnostic"]
    method_labels = [
        "moments SFS (primary)",
        "dadi SFS",
        "Joint SFS + LD",
        "LD only",
        "Context ensemble",
        "Score calibrator",
        "FastCXT diagnostic",
    ]
    y = np.arange(len(method_order))[::-1]
    for column, ax in enumerate(axes[1]):
        values = np.asarray([observed[method][column] for method in method_order])
        for position, method, value in zip(y, method_order, values):
            if method == "diagnostic":
                ax.scatter(
                    value,
                    position,
                    s=32,
                    marker="x",
                    color=COLORS[method],
                    linewidths=1.2,
                    zorder=3,
                )
            else:
                ax.scatter(
                    value,
                    position,
                    s=38 if method == "moments" else 28,
                    marker="o",
                    facecolors="none" if method in ("learned", "score") else COLORS[method],
                    edgecolors=COLORS[method],
                    linewidths=1.0,
                    zorder=3,
                )
        interval = generation_interval if column == 0 else fraction_interval
        primary = observed["moments"][column]
        ax.errorbar(
            primary,
            y[0],
            xerr=np.asarray([[primary - interval[0]], [interval[1] - primary]]),
            fmt="none",
            ecolor=COLORS["moments"],
            elinewidth=1.0,
            capsize=2.5,
            zorder=2,
        )
        ax.axvspan(interval[0], interval[1], color=COLORS["moments"], alpha=0.07, linewidth=0)
        ax.set_yticks(y, method_labels)
        ax.grid(axis="x", color="#D9D9D9", linewidth=0.45)
        ax.set_axisbelow(True)
        if column == 0:
            ax.set_xlabel("Observed decline time (generations)")
        else:
            ax.set_xlabel("Observed post-decline fraction")

    for label, ax in zip("abcd", axes.flat):
        ax.text(-0.13, 1.04, label, transform=ax.transAxes, fontsize=9, fontweight="normal")

    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(args.output_stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(args.output_stem.with_suffix(".png"), dpi=450, bbox_inches="tight")
    strip_svg_metadata(args.output_stem.with_suffix(".svg"))
    plt.close(fig)


if __name__ == "__main__":
    main()
