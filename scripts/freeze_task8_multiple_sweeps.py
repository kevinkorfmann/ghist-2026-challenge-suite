"""Freeze only a candidate with matched, independently evaluated recovery.

The archived Task 7 null / Task 5 recovery transfer is insufficient for Task 8.
This entry point deliberately requires new explicit evidence and output paths.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ghist_submission.validate import validate_sweeps


def validated_candidate(candidate, recovery, candidate_sha256):
    if recovery.get("candidate_sha256") != candidate_sha256:
        raise ValueError("Recovery must identify the exact candidate configuration")
    if (
        recovery.get("public_testing_answers_used") is not False
        or recovery.get("independent_evaluation") is not True
    ):
        raise ValueError(
            "Recovery must use independent evaluation, without public testing answers"
        )
    if (
        recovery.get("demography") != "bottleneck_growth"
        or recovery.get("background_selection") is not True
        or recovery.get("multiple_sweeps") is not True
    ):
        raise ValueError(
            "Task 8 needs matched multiple-sweep, bottleneck/growth, and BGS recovery"
        )
    # Thresholds must be chosen before evaluating new replicates. No gate is
    # based on a fixed number of peaks in the observed challenge dataset.
    for metric in ("precision", "recall", "interval_coverage"):
        value = recovery.get("evaluation", {}).get(metric)
        minimum = recovery.get("prespecified_minimum", {}).get(metric)
        if (
            not isinstance(value, (int, float))
            or not isinstance(minimum, (int, float))
            or not np.isfinite([value, minimum]).all()
            or not 0 < minimum <= value <= 1
        ):
            raise ValueError(
                f"Matched recovery fails its prespecified {metric} requirement"
            )
    if recovery.get("evaluation_replicates", 0) < 2:
        raise ValueError("Independent recovery needs multiple evaluation replicates")
    if not candidate.get("intervals"):
        raise ValueError(
            "No supported sweep candidates; refusing to fabricate a region"
        )
    return candidate["intervals"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--recovery", type=Path, required=True)
    parser.add_argument("--dataset", choices=("testing", "final"), required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    candidate = json.loads(args.candidate.read_text())
    recovery = json.loads(args.recovery.read_text())
    digest = hashlib.sha256(args.candidate.read_bytes()).hexdigest()
    intervals = validated_candidate(candidate, recovery, digest)
    chromosome = "7" if args.dataset == "testing" else "9"
    if candidate.get("chromosome") != chromosome:
        raise ValueError("Candidate chromosome does not match the requested dataset")
    # Validate a temporary file before replacing any previous candidate.
    import tempfile

    args.bed.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", dir=args.bed.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        for start, end in intervals:
            handle.write(f"chr{chromosome}\t{start}\t{end}\n")
    try:
        validation = validate_sweeps(temporary, 8, args.dataset)
        temporary.replace(args.bed)
    finally:
        temporary.unlink(missing_ok=True)
    result = {
        "challenge": 8,
        "status": "method frozen; no submission",
        "submission_count": 0,
        "dataset": args.dataset,
        "prediction": validation,
        "candidate_sha256": digest,
        "recovery": recovery,
        "recovery_sha256": hashlib.sha256(args.recovery.read_bytes()).hexdigest(),
        "bed_sha256": hashlib.sha256(args.bed.read_bytes()).hexdigest(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
