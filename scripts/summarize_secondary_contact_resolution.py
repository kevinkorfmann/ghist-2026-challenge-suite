from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PARAMETERS = (
    ("population1_effective_size", "Mainland size"),
    ("population2_effective_size", "Island size"),
    ("generations_since_split", "Split time"),
    ("generations_since_contact", "Contact time"),
    ("migration_rate", "Migration rate"),
)
BOOTSTRAP_KEYS = {
    "population1_effective_size": "mainland_population_size",
    "population2_effective_size": "island_population_size",
    "generations_since_split": "generations_since_split",
    "generations_since_contact": "generations_since_migration",
    "migration_rate": "migration_rate",
}


def labelled_path(value: str) -> tuple[str, Path]:
    label, separator, path = value.partition("=")
    if not separator or not label or not path:
        raise argparse.ArgumentTypeError("expected LABEL=PATH")
    return label, Path(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fit", action="append", type=labelled_path, required=True)
    parser.add_argument("--bootstrap", action="append", type=labelled_path, default=[])
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, required=True)
    args = parser.parse_args()

    fits = []
    for label, path in args.fit:
        source = json.loads(path.read_text())
        fits.append(
            {
                "label": label,
                "path": str(path),
                "projection": source.get("projection"),
                "folded": source["folded"],
                "integration_dt_fac": source.get("integration_dt_fac", 0.02),
                "adaptive_integration": source.get("adaptive_integration", False),
                "values": {key: source[key] for key, _ in PARAMETERS},
            }
        )
    bootstraps = []
    for label, path in args.bootstrap:
        source = json.loads(path.read_text())
        bootstraps.append(
            {
                "label": label,
                "path": str(path),
                "replicates": source["replicates"],
                "competition_rrmse": source["competition_rrmse"],
                "values": {
                    key: source["observed_bias_corrected"][BOOTSTRAP_KEYS[key]]
                    for key, _ in PARAMETERS
                },
                "intervals": {
                    key: source["observed_sampling_interval_95"][BOOTSTRAP_KEYS[key]]
                    for key, _ in PARAMETERS
                },
            }
        )
    value_matrix = np.asarray(
        [[fit["values"][key] for key, _ in PARAMETERS] for fit in fits]
    )
    result = {
        "fits": fits,
        "bootstraps": bootstraps,
        "fit_relative_range": dict(
            zip(
                (key for key, _ in PARAMETERS),
                ((value_matrix.max(axis=0) - value_matrix.min(axis=0))
                 / np.median(value_matrix, axis=0)).tolist(),
            )
        ),
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2) + "\n")

    labels = [fit["label"] for fit in fits]
    positions = np.arange(len(labels))
    figure, axes = plt.subplots(1, len(PARAMETERS), figsize=(15.5, 3.8))
    for axis, (key, title) in zip(axes, PARAMETERS):
        values = np.asarray([fit["values"][key] for fit in fits])
        axis.plot(positions, values, marker="o", color="#315B7D", linewidth=1.7)
        for offset, bootstrap in enumerate(bootstraps):
            center = bootstrap["values"][key]
            lower, upper = bootstrap["intervals"][key]
            x = len(labels) + 0.65 * offset
            axis.errorbar(
                x,
                center,
                yerr=[[center - lower], [upper - center]],
                marker="D",
                color="#B05A3C",
                capsize=3,
            )
        axis.set_title(title, fontsize=10)
        axis.grid(axis="y", alpha=0.22)
        axis.tick_params(axis="y", labelsize=8)
        if key == "migration_rate":
            axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, -3))
    all_labels = labels + [item["label"] for item in bootstraps]
    all_positions = list(positions) + [len(labels) + 0.65 * i for i in range(len(bootstraps))]
    for axis in axes:
        axis.set_xticks(all_positions, all_labels, rotation=35, ha="right", fontsize=8)
    figure.suptitle("Secondary Contact resolution and finite-site calibration", fontsize=12)
    figure.tight_layout()
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out_prefix.with_suffix(".png"), dpi=240, bbox_inches="tight")
    figure.savefig(args.out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
