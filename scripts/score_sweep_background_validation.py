"""Score a prespecified, segment-scale independent sweep benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ghist_single_sweep.multiple import calibrate_maximum, detect_regions
from ghist_single_sweep.peaks import read_raisd


def region_metrics(intervals, positions):
    covered = [any(start <= p < end for start, end in intervals) for p in positions]
    true_regions = [
        any(start <= p < end for p in positions) for start, end in intervals
    ]
    return {
        "true_sweeps": len(positions),
        "covered_sweeps": sum(covered),
        "detected_regions": len(intervals),
        "true_regions": sum(true_regions),
        "false_regions": len(intervals) - sum(true_regions),
        "covered_bases": sum(end - start for start, end in intervals),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--coding-bed", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    design = json.loads(args.design.read_text())
    records = {}
    for part, seeds in (
        ("calibration", design["calibration_null_seeds"]),
        ("evaluation", design["evaluation_seeds"]),
    ):
        records[part] = []
        for index, seed in enumerate(seeds):
            path = args.input_dir / f"{part}_{index}.json"
            record = json.loads(path.read_text())
            if record["seed"] != seed:
                raise ValueError("Simulation seed does not match prespecified design")
            if part == "calibration" and record["truth_positions"]:
                raise ValueError("Null calibration contains a selected sweep")
            if (
                part == "evaluation"
                and len(record["truth_positions"]) != design["evaluation_counts"][index]
            ):
                raise ValueError("Sweep count does not match prespecified design")
            record["report"] = read_raisd(path.with_suffix(".raisd.txt"))
            record["record_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            record["report_sha256"] = hashlib.sha256(
                path.with_suffix(".raisd.txt").read_bytes()
            ).hexdigest()
            records[part].append(record)
    context = records["calibration"][0]["generator"]
    if any(
        record["generator"] != context for group in records.values() for record in group
    ):
        raise ValueError("Null and evaluation generator settings differ")
    exons = []
    for line in args.coding_bed.read_text().splitlines():
        fields = line.split()
        if fields and not line.startswith("#") and fields[0].removeprefix("chr") == "7":
            a, b = (
                max(context["left"], int(fields[1])),
                min(context["right"], int(fields[2])),
            )
            if a < b:
                exons.append((a, b))
    maxima = [float(record["report"][:, 6].max()) for record in records["calibration"]]
    settings = dict(design["detector"])
    threshold = calibrate_maximum(maxima, settings.pop("alpha"))
    evaluations = []
    for record in records["evaluation"]:
        prediction = detect_regions(
            record["report"],
            threshold=threshold,
            chromosome_length=context["right"],
            exons=exons,
            **settings,
        )
        prediction["intervals"] = [
            [max(context["left"], start), min(context["right"], end)]
            for start, end in prediction["intervals"]
        ]
        scores = region_metrics(prediction["intervals"], record["truth_positions"])
        evaluations.append(
            {
                "seed": record["seed"],
                "truth_positions": record["truth_positions"],
                "intervals": prediction["intervals"],
                **scores,
                "single_sweep_errors": record["single_sweep_errors"],
                "raisd_top_peak_error": abs(
                    record["raisd_peaks"][0]["position"] - record["truth_positions"][0]
                )
                if len(record["truth_positions"]) == 1
                else None,
                "record_sha256": record["record_sha256"],
                "report_sha256": record["report_sha256"],
            }
        )
    strata = {}
    for count in (0, 1, 2):
        group = [record for record in evaluations if record["true_sweeps"] == count]
        regions = sum(record["detected_regions"] for record in group)
        total_truth = sum(record["true_sweeps"] for record in group)
        strata[str(count)] = {
            "replicates": len(group),
            "replicates_with_any_detection": sum(
                record["detected_regions"] > 0 for record in group
            ),
            "covered_sweeps": sum(record["covered_sweeps"] for record in group),
            "true_sweeps": total_truth,
            "true_regions": sum(record["true_regions"] for record in group),
            "detected_regions": regions,
            "recall": sum(record["covered_sweeps"] for record in group) / total_truth
            if total_truth
            else None,
            "region_precision": sum(record["true_regions"] for record in group)
            / regions
            if regions
            else None,
            "mean_covered_bases": float(
                np.mean([record["covered_bases"] for record in group])
            ),
        }
        if count == 1:
            methods = list(group[0]["single_sweep_errors"]) + ["raisd"]
            strata[str(count)]["unconditional_localization"] = {}
            for method in methods:
                errors = np.asarray(
                    [
                        record["raisd_top_peak_error"]
                        if method == "raisd"
                        else record["single_sweep_errors"][method]
                        for record in group
                    ]
                )
                strata[str(count)]["unconditional_localization"][method] = {
                    "median_absolute_error": float(np.median(errors)),
                    "q90_absolute_error": float(np.quantile(errors, 0.9)),
                    "within_80kb": int(np.sum(errors <= 80000)),
                }
    payload = {
        "status": "segment-scale scientific diagnostic; not submission calibration",
        "design": design,
        "generator": context,
        "null_maxima": maxima,
        "threshold": threshold,
        "strata": strata,
        "evaluation_records": evaluations,
        "calibration_record_hashes": [
            r["record_sha256"] for r in records["calibration"]
        ],
        "calibration_report_hashes": [
            r["report_sha256"] for r in records["calibration"]
        ],
        "independent_unit": "one independently simulated 1 Mb segment",
        "public_testing_answers_used": False,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(strata, indent=2))


if __name__ == "__main__":
    main()
