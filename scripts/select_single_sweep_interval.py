from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


EVALUATOR_WEIGHT = 1e-5
EVALUATOR_COMMIT = "43b95b482ab0b3a39659226f51088f1c65ba39c9"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def competition_score(errors: np.ndarray, half_width: float) -> np.ndarray:
    hit = errors <= half_width
    weighted_false_positive = 2.0 * half_width * EVALUATOR_WEIGHT
    return hit.astype(float) / (1.0 + weighted_false_positive / 2.0)


def ensemble_errors(summary: dict, seed_origin: int) -> tuple[np.ndarray, np.ndarray]:
    if summary["held_out_method_gate"]["selected"] != "median_all_three":
        raise ValueError("recovery campaign did not select the locked ensemble")
    records = summary["records"]
    truth = np.asarray([record["truth"]["position"] for record in records], dtype=float)
    estimates = np.asarray(
        [
            np.median(
                [
                    record["estimates"]["largest_gap_midpoint"],
                    record["estimates"]["composite"],
                    record["estimates"]["flank_ld_100kb"],
                ]
            )
            for record in records
        ],
        dtype=float,
    )
    task_indices = np.asarray([record["seed"] - seed_origin for record in records])
    selection = np.isin((task_indices * 17 + 3) % 5, [0, 1])
    return np.abs(estimates - truth), selection


def choose_interval(
    summaries: list[dict],
    observed_estimates: dict[str, float],
    half_widths: np.ndarray,
    seed_origin: int = 75_222_000,
    safety_padding: int = 2_000,
) -> dict:
    domains = [ensemble_errors(summary, seed_origin) for summary in summaries]
    selection_scores = {
        str(int(width)): float(
            np.mean(
                [
                    competition_score(errors[selection], width).mean()
                    for errors, selection in domains
                ]
            )
        )
        for width in half_widths
    }
    selected_width = int(max(selection_scores, key=selection_scores.get))
    center = float(
        np.median(
            [
                observed_estimates["largest_gap_midpoint"],
                observed_estimates["composite_100kb"],
                observed_estimates["flank_ld_100kb"],
            ]
        )
    )
    maximum_deviation = max(abs(value - center) for value in observed_estimates.values())
    consensus_guard = int(math.ceil(maximum_deviation / 1_000.0) * 1_000 + safety_padding)
    final_width = max(selected_width, consensus_guard)
    evaluation = []
    for errors, selection in domains:
        held_out = errors[~selection]
        evaluation.append(
            {
                "replicates": int(held_out.size),
                "coverage": float(np.mean(held_out <= final_width)),
                "mean_competition_score": float(
                    competition_score(held_out, final_width).mean()
                ),
            }
        )
    return {
        "center": int(round(center)),
        "selection_optimal_half_width": selected_width,
        "maximum_observed_estimator_deviation": maximum_deviation,
        "consensus_guard_half_width": consensus_guard,
        "final_half_width": final_width,
        "interval": [int(round(center - final_width)), int(round(center + final_width))],
        "selection_scores": selection_scores,
        "held_out_evaluation_after_lock": evaluation,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generic-recovery", type=Path, required=True)
    parser.add_argument("--chromosome-recovery", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--localization", type=Path, required=True)
    parser.add_argument("--raisd", type=Path, required=True)
    parser.add_argument("--exact-recovery", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    args = parser.parse_args()

    paths = [
        args.generic_recovery,
        args.chromosome_recovery,
        args.audit,
        args.localization,
        args.raisd,
        args.exact_recovery,
    ]
    generic, chromosome, audit, localization, raisd, exact = [
        json.loads(path.read_text()) for path in paths
    ]
    observed = {
        "largest_gap_midpoint": float(localization["largest_variant_gap"]["midpoint"]),
        "composite_100kb": float(audit["scales"]["100000"]["top_windows"][0]["center"]),
        "flank_ld_100kb": float(
            localization["flank_ld_scans"]["100000"]["maximum"]["candidate_position"]
        ),
        "raisd": float(raisd["maximum"]["position"]),
    }
    half_widths = np.asarray(
        [5_000, 10_000, 15_000, 20_000, 25_000, 30_000, 40_000, 50_000, 60_000,
         75_000, 100_000, 125_000, 150_000, 200_000],
        dtype=float,
    )
    interval = choose_interval([generic, chromosome], observed, half_widths)
    interval["pre_exact_half_width"] = interval["final_half_width"]
    interval["exact_calibration_half_width"] = exact["interval_width_selection"]["selected_half_width"]
    interval["final_half_width"] = max(
        interval["final_half_width"], interval["exact_calibration_half_width"]
    )
    interval["interval"] = [
        interval["center"] - interval["final_half_width"],
        interval["center"] + interval["final_half_width"],
    ]
    result = {
        "state": "CANDIDATE; NOT FROZEN",
        "chromosome": audit["contig"],
        "selected_method": "median of gap, 100 kb composite, and 100 kb split-flank linkage",
        "observed_estimates": observed,
        "interval_selection": interval,
        "evaluator_alignment": {
            "formula": "hit / (1 + interval_length * 1e-5 / 2)",
            "source_repository": "https://github.com/tjstruck/GHIST_2025_Sweeps",
            "source_commit": EVALUATOR_COMMIT,
        },
        "source_artifacts": {
            str(path): sha256(path) for path in paths
        },
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    args.bed.parent.mkdir(parents=True, exist_ok=True)
    args.bed.write_text(
        f"{audit['contig']}\t{interval['interval'][0]}\t{interval['interval'][1]}\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
