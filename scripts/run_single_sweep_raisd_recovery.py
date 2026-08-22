from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

import tskit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trees", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--raisd", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    truth_record = json.loads(args.truth.read_text())
    tree_sequence = tskit.load(args.trees)
    with tempfile.TemporaryDirectory(prefix="single-sweep-raisd-") as directory:
        workdir = Path(directory)
        vcf_path = workdir / "simulation.vcf"
        with vcf_path.open("w") as handle:
            tree_sequence.write_vcf(
                handle,
                contig_id="simulated",
                allow_position_zero=True,
            )
        subprocess.run(
            [
                str(args.raisd.resolve()),
                "-n",
                "recovery",
                "-I",
                str(vcf_path),
                "-w",
                "50",
                "-R",
                "-f",
            ],
            cwd=workdir,
            env=os.environ.copy(),
            check=True,
            stdout=subprocess.DEVNULL,
        )
        maximum = None
        with (workdir / "RAiSD_Report.recovery").open() as handle:
            for line in handle:
                fields = line.split()
                if not fields or fields[0] == "//" or len(fields) < 7:
                    continue
                row = {
                    "position": int(fields[0]),
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                    "variation": float(fields[3]),
                    "frequency_spectrum": float(fields[4]),
                    "linkage_disequilibrium": float(fields[5]),
                    "mu": float(fields[6]),
                }
                if maximum is None or row["mu"] > maximum["mu"]:
                    maximum = row
    if maximum is None:
        raise ValueError("reference report contained no evaluated locations")
    true_position = float(truth_record["truth"]["position"])
    result = {
        "state": "COMPLETED",
        "truth_position": true_position,
        "maximum": maximum,
        "absolute_error": abs(maximum["position"] - true_position),
        "source_trees": str(args.trees),
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

