from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit, logit
from scipy.stats import qmc


PARAMETER_NAMES = (
    "ancestral_population_size",
    "mainland_population_size",
    "island_population_size",
    "generations_since_split",
    "generations_since_contact",
    "migration_rate",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--points", type=int, default=24)
    parser.add_argument("--replicates", type=int, default=2)
    parser.add_argument(
        "--center",
        type=float,
        nargs=6,
        default=(32_500, 21_000, 76_000, 4_600, 1_550, 0.0002),
    )
    parser.add_argument(
        "--half-width-factors",
        type=float,
        nargs=5,
        default=(1.30, 1.30, 1.35, 1.30, 1.50),
        metavar=("ANCESTRAL", "MAINLAND", "ISLAND", "SPLIT", "MIGRATION"),
    )
    parser.add_argument("--contact-logit-half-width", type=float, default=0.7)
    parser.add_argument("--seed", type=int, default=20260823)
    args = parser.parse_args()
    if args.points < 10 or args.replicates < 1:
        parser.error("at least ten points and one replicate are required")
    center = np.asarray(args.center, dtype=np.float64)
    if np.any(center <= 0) or center[4] >= center[3]:
        parser.error("the center must have positive values and contact younger than split")
    half_widths = np.asarray(args.half_width_factors, dtype=np.float64)
    if np.any(half_widths <= 1) or args.contact_logit_half_width <= 0:
        parser.error("all design half-widths must be greater than one or zero")

    unit = qmc.LatinHypercube(6, seed=args.seed).random(args.points)
    values = np.empty((args.points, 6), dtype=np.float64)
    logged_columns = (0, 1, 2, 3, 5)
    values[:, logged_columns] = np.exp(
        np.log(center[list(logged_columns)])
        + (2.0 * unit[:, logged_columns] - 1.0) * np.log(half_widths)
    )
    center_fraction = center[4] / center[3]
    fractions = expit(
        logit(center_fraction)
        + (2.0 * unit[:, 4] - 1.0) * args.contact_logit_half_width
    )
    values[:, 4] = values[:, 3] * fractions
    # Anchor the design at its nominal center and keep that point in calibration.
    values[0] = center

    n_calibration = round(0.60 * args.points)
    n_selection = round(0.20 * args.points)
    order = np.random.default_rng(args.seed + 1).permutation(np.arange(1, args.points))
    calibration = {0, *order[: n_calibration - 1].tolist()}
    selection = set(order[n_calibration - 1 : n_calibration - 1 + n_selection].tolist())
    rows = []
    for point in range(args.points):
        if point in calibration:
            partition = "calibration"
        elif point in selection:
            partition = "selection"
        else:
            partition = "evaluation"
        for replicate in range(args.replicates):
            rows.append(
                {
                    "task_index": len(rows),
                    "point": point,
                    "replicate": replicate,
                    "partition": partition,
                    "truth": values[point].tolist(),
                }
            )

    result = {
        "parameter_names": PARAMETER_NAMES,
        "points": args.points,
        "replicates_per_point": args.replicates,
        "center": center.tolist(),
        "half_width_factors": half_widths.tolist(),
        "contact_logit_half_width": args.contact_logit_half_width,
        "seed": args.seed,
        "partition_point_counts": {
            name: len({row["point"] for row in rows if row["partition"] == name})
            for name in ("calibration", "selection", "evaluation")
        },
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2))


if __name__ == "__main__":
    main()
