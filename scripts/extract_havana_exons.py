from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--chromosome", default="7")
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    intervals: list[tuple[int, int]] = []
    with gzip.open(args.gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[0] != args.chromosome or fields[2] != "exon":
                continue
            if fields[1] != "havana" and 'gene_source "havana"' not in fields[8]:
                continue
            intervals.append((int(fields[3]) - 1, int(fields[4])))
    intervals.sort()
    merged: list[list[int]] = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    args.bed.parent.mkdir(parents=True, exist_ok=True)
    with args.bed.open("w") as handle:
        for start, end in merged:
            handle.write(f"chr{args.chromosome}\t{start}\t{end}\n")
    summary = {
        "state": "COMPLETED",
        "source": {"path": str(args.gtf), "sha256": sha256(args.gtf)},
        "chromosome": args.chromosome,
        "havana_exon_records": len(intervals),
        "merged_intervals": len(merged),
        "merged_bases": sum(end - start for start, end in merged),
        "bed": str(args.bed),
        "submission_count": 0,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
