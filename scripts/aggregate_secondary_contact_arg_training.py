from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ghist_secondary_contact.contract import (
    SPECTRUM_CONTRACT,
    require_spectrum_contract,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cells", type=int, required=True)
    args = parser.parse_args()
    records = []
    missing = []
    for cell in range(args.cells):
        root = args.root / f"cell_{cell:05d}"
        simulation_path = root / "simulated.features.npz"
        singer_path = root / "singer_hazards.npz"
        polegon_path = root / "polegon_hazards.npz"
        if not all(
            path.exists() for path in (simulation_path, singer_path, polegon_path)
        ):
            missing.append(cell)
            continue
        simulation = np.load(simulation_path)
        require_spectrum_contract(simulation)
        singer = np.load(singer_path)
        polegon = np.load(polegon_path)
        records.append(
            {
                "cell_id": int(simulation["cell_id"]),
                "point_id": int(simulation["point_id"]),
                "partition": int(simulation["partition"]),
                "design_domain": int(simulation["design_domain"]),
                "truth": simulation["truth"],
                "ancestral_size": float(simulation["ancestral_size"]),
                "nuisance": simulation["nuisance"],
                "truth_hazards": simulation["truth_hazards"],
                "global_summary": np.concatenate(
                    (
                        singer["mean"].reshape(-1),
                        singer["standard_deviation"].reshape(-1),
                        polegon["mean"].reshape(-1),
                    )
                ),
            }
        )
    if missing:
        raise FileNotFoundError(f"missing ARG training cells: {missing[:20]}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        spectrum_contract=SPECTRUM_CONTRACT,
        spectrum_transform="log_l2",
        cell_id=np.asarray([record["cell_id"] for record in records]),
        point_id=np.asarray([record["point_id"] for record in records]),
        partition=np.asarray(
            [record["partition"] for record in records], dtype=np.int8
        ),
        design_domain=np.asarray(
            [record["design_domain"] for record in records], dtype=np.int8
        ),
        truth=np.asarray([record["truth"] for record in records], dtype=np.float64),
        ancestral_size=np.asarray(
            [record["ancestral_size"] for record in records], dtype=np.float64
        ),
        nuisance=np.asarray(
            [record["nuisance"] for record in records], dtype=np.float64
        ),
        truth_hazards=np.asarray(
            [record["truth_hazards"] for record in records], dtype=np.float32
        ),
        windows=np.zeros((args.cells, 1, 1), dtype=np.float32),
        joint_sfs=np.zeros((args.cells, 87, 49), dtype=np.float32),
        global_summary=np.asarray(
            [record["global_summary"] for record in records], dtype=np.float32
        ),
    )
    report = {
        "cells": args.cells,
        "parameter_points": len({record["point_id"] for record in records}),
        "global_features": len(records[0]["global_summary"]),
        "partitions": {
            str(value): int(sum(record["partition"] == value for record in records))
            for value in (0, 1, 2)
        },
        "domains": {
            str(value): int(sum(record["design_domain"] == value for record in records))
            for value in (0, 1, 2)
        },
    }
    args.out.with_suffix(args.out.suffix + ".json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
