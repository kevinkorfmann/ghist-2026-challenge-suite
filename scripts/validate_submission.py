from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import yaml

from ghist_submission.validate import all_pairs, validate_submission


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate locally; never upload or submit"
    )
    parser.add_argument("file", type=Path)
    parser.add_argument("--challenge", type=int, choices=range(1, 13), required=True)
    parser.add_argument("--dataset", choices=("testing", "final"), required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--samples",
        type=Path,
        help="One sample ID per line; require every unordered pair",
    )
    group.add_argument(
        "--pairs",
        type=Path,
        help="CSV with indID1/indID2 columns listing expected pairs",
    )
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    expected = None
    try:
        if args.samples:
            expected = all_pairs(args.samples.read_text().splitlines())
        if args.pairs:
            with args.pairs.open(newline="") as handle:
                expected = {
                    (row["indID1"], row["indID2"]) for row in csv.DictReader(handle)
                }
        report = validate_submission(args.file, args.challenge, args.dataset, expected)
    except (ValueError, KeyError, OSError, yaml.YAMLError) as error:
        report = {"format_valid": False, "error": str(error)}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["format_valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
