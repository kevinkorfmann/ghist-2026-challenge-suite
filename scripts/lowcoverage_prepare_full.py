"""Generate a fixed independent full-length low-coverage validation cohort."""
import argparse
from pathlib import Path

from ghist_lowcoverage.pipeline import write_json
from ghist_lowcoverage.remap import remap_synthetic_bams
from scripts.lowcoverage_validate import DESIGNS, simulate_bams


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', choices=DESIGNS, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bwa', type=Path, required=True)
    parser.add_argument('--seed-offset', type=int, default=10000)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    simulation = simulate_bams(args.task, args.output / 'synthetic', 10_000_000, seed_offset=args.seed_offset)
    remap = remap_synthetic_bams(args.output / 'synthetic', args.output / 'remapped', args.bwa, args.threads)
    write_json(args.output / 'prepared.json', {'simulation': simulation, 'remap': remap, 'inference_performed': False})


if __name__ == '__main__':
    main()
