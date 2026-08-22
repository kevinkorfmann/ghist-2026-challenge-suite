from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--info", type=Path, required=True)
    parser.add_argument("--site-report", type=Path, required=True)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    with args.report.open() as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0] == "//" or len(fields) < 7:
                continue
            rows.append(
                {
                    "position": int(fields[0]),
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                    "variation": float(fields[3]),
                    "frequency_spectrum": float(fields[4]),
                    "linkage_disequilibrium": float(fields[5]),
                    "mu": float(fields[6]),
                }
            )
    rows.sort(key=lambda row: row["mu"], reverse=True)
    selected = []
    for row in rows:
        if all(abs(row["position"] - prior["position"]) >= 100_000 for prior in selected):
            selected.append(row)
        if len(selected) == 20:
            break
    result = {
        "state": "COMPLETED",
        "version": "2.9",
        "window_sites": 50,
        "evaluated_locations": len(rows),
        "maximum": rows[0],
        "separated_peaks_100kb": selected,
        "job_id": args.job_id,
        "source_artifacts": {
            "report": {"path": str(args.report), "sha256": sha256(args.report)},
            "info": {"path": str(args.info), "sha256": sha256(args.info)},
            "site_report": {
                "path": str(args.site_report),
                "sha256": sha256(args.site_report),
            },
        },
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

