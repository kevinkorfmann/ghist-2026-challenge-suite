"""Run the developmental likelihood pipeline on explicitly supplied BAMs.

This command does not retrieve challenge data or write a submission. Estimates
remain marked development-only until independent scientific release gates pass.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ghist_lowcoverage.pipeline import TASKS, run_pipeline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", choices=TASKS, required=True)
    parser.add_argument("--bam-dir", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--angsd", type=Path, required=True)
    parser.add_argument("--realsfs", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--restarts", type=int, default=8)
    parser.add_argument("--gl-model", type=int, choices=[1, 2], default=1)
    args = parser.parse_args()
    if args.threads < 1 or args.restarts < 2:
        parser.error("positive threads and at least two optimizer restarts required")
    result = run_pipeline(args.bam_dir, args.reference, args.task, args.output, args.angsd, args.realsfs,
                          threads=args.threads, restarts=args.restarts, genotype_likelihood_model=args.gl_model)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
