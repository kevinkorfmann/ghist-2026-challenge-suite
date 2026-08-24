from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def read_raisd(path: Path) -> np.ndarray:
    return np.loadtxt(path, comments="/")


def separated_peaks(report: np.ndarray, separation: int) -> list[int]:
    order = np.argsort(report[:, 6])[::-1]
    selected: list[int] = []
    for index in order:
        position = int(report[index, 0])
        if all(abs(position - int(report[other, 0])) >= separation for other in selected):
            selected.append(int(index))
    return selected


def nearest_exonic_position(position: int, intervals: list[tuple[int, int]]) -> int:
    best_position = position
    best_distance = float("inf")
    for start, end in intervals:
        candidate = min(max(position, start), end - 1)
        distance = abs(candidate - position)
        if distance < best_distance:
            best_position = candidate
            best_distance = distance
    return best_position


def local_audit_centers(directory: Path, position: int) -> list[int]:
    centers: list[int] = []
    for width in (25_000, 50_000, 100_000):
        scan = np.genfromtxt(directory / f"scan_{width}.csv", delimiter=",", names=True)
        mask = np.abs(scan["center"] - position) <= 500_000
        local = np.flatnonzero(mask & np.isfinite(scan["composite_score"]))
        if local.size:
            index = local[np.argmax(scan["composite_score"][local])]
            centers.append(int(scan["center"][index]))
    return centers


def localization_centers(path: Path) -> tuple[list[int], dict[str, object]]:
    summary = json.loads(path.read_text())
    centers = [int(summary["largest_variant_gap"]["midpoint"])]
    for width in ("50000", "100000", "250000"):
        centers.append(int(summary["flank_ld_scans"][width]["maximum"]["candidate_position"]))
    return centers, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task7-dir", type=Path, required=True)
    parser.add_argument("--task8-dir", type=Path, required=True)
    parser.add_argument("--exons", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    args = parser.parse_args()

    task7_report_path = next((args.task7_dir / "raisd").glob("RAiSD_Report.*"))
    task8_report_path = next((args.task8_dir / "raisd").glob("RAiSD_Report.*"))
    task7_report = read_raisd(task7_report_path)
    task8_report = read_raisd(task8_report_path)
    task7_primary = int(task7_report[np.argmax(task7_report[:, 6]), 0])
    task7_null = task7_report[np.abs(task7_report[:, 0] - task7_primary) > 2_000_000]
    task7_null_maximum = float(np.max(task7_null[:, 6]))
    threshold = 2.0 * task7_null_maximum

    task7_evidence = [task7_primary]
    task7_evidence.extend(local_audit_centers(args.task7_dir / "baseline_scans", task7_primary))
    local_centers, task7_localization = localization_centers(args.task7_dir / "localization" / "summary.json")
    task7_evidence.extend(local_centers)
    task7_center = int(np.median(task7_evidence))

    exons = []
    with args.exons.open() as handle:
        for line in handle:
            _, start, end = line.split()
            exons.append((int(start), int(end)))

    candidates = []
    task8_localizations = []
    for localization_path in sorted((args.task8_dir / "localization").glob("region_*/summary.json")):
        localization = json.loads(localization_path.read_text())
        midpoint = float(np.mean(localization["region"]))
        task8_localizations.append((midpoint, localization_path))
    for peak_index in separated_peaks(task8_report, 1_000_000):
        peak = task8_report[peak_index]
        if peak[6] < threshold:
            continue
        position = int(peak[0])
        _, localization_path = min(
            task8_localizations, key=lambda item: abs(item[0] - position)
        )
        local_centers, localization = localization_centers(localization_path)
        evidence = [position]
        evidence.extend(local_audit_centers(args.task8_dir / "baseline_scans", position))
        evidence.extend(local_centers)
        consensus = int(np.median(evidence))
        exonic_center = nearest_exonic_position(consensus, exons)
        candidates.append(
            {
                "raisd_position": position,
                "raisd_mu": float(peak[6]),
                "raisd_to_task7_null_ratio": float(peak[6] / task7_null_maximum),
                "evidence_centers": evidence,
                "consensus_center": consensus,
                "exonic_center": exonic_center,
                "exon_snap_distance": abs(exonic_center - consensus),
                "localization": localization,
            }
        )

    task7_half_width = 80_000
    task8_half_width = 60_000
    summary = {
        "state": "PROVISIONAL",
        "calibration": {
            "task7_primary_position": task7_primary,
            "task7_background_maximum_mu_outside_2mb": task7_null_maximum,
            "task8_detection_threshold": threshold,
            "rule": "one-megabase-separated RAiSD peaks exceeding twice the Task 7 BGS null maximum",
        },
        "task7": {
            "center": task7_center,
            "start": task7_center - task7_half_width,
            "end": task7_center + task7_half_width,
            "evidence_centers": task7_evidence,
            "localization": task7_localization,
        },
        "task8": candidates,
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    args.bed.parent.mkdir(parents=True, exist_ok=True)
    with args.bed.open("w") as handle:
        for candidate in sorted(candidates, key=lambda item: int(item["exonic_center"])):
            center = int(candidate["exonic_center"])
            handle.write(f"chr7\t{center - task8_half_width}\t{center + task8_half_width}\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
