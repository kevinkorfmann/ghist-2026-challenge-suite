from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ghist_secondary_contact.coalescent import fit_secondary_contact_interval_counts
from ghist_secondary_contact.demography import SecondaryContactParameters


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--initial", type=float, nargs=6, required=True, metavar=(
        "ANCESTRAL", "MAINLAND", "ISLAND", "SPLIT", "CONTACT", "MIGRATION"
    ))
    parser.add_argument("--seed", type=int, default=20260822)
    parser.add_argument("--local-only", action="store_true")
    args = parser.parse_args()
    data = np.load(args.features)
    if "polegon_global_interval_counts" in data:
        counts = data["polegon_global_interval_counts"]
    elif "total_interval_counts" in data:
        counts = data["total_interval_counts"]
    elif "interval_counts" in data:
        counts = data["interval_counts"].sum(axis=0)
    else:
        raise KeyError("feature archive does not contain interval counts")
    ancestral, mainland, island, split, contact, migration = args.initial
    initial = SecondaryContactParameters(
        ancestral_size=ancestral,
        population1_size=mainland,
        population2_size=island,
        split_time=split,
        contact_time=contact,
        migration_rate=migration,
    )
    result = fit_secondary_contact_interval_counts(
        counts,
        initial=initial,
        global_search=not args.local_only,
        seed=args.seed,
    )
    result["source_features"] = str(args.features)
    result["policy"] = "local_fit_only"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        key: value for key, value in result.items()
        if key not in {
            "empirical_interval_probabilities",
            "predicted_interval_probabilities",
            "predicted_hazards",
        }
    }, indent=2))


if __name__ == "__main__":
    main()
