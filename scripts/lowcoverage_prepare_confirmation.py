"""Prepare one prespecified confirmation cohort without demographic inference."""
import argparse
import json
from pathlib import Path

from ghist_lowcoverage.pipeline import sha256, write_json
from ghist_lowcoverage.remap import remap_synthetic_bams
from scripts.lowcoverage_validate import simulate_bams


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--task', choices=['bottleneck', 'split_isolation'], required=True)
    parser.add_argument('--regime', type=int, choices=range(3), required=True)
    parser.add_argument('--replicate', type=int, choices=[0, 1], required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bwa', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    design = dict(protocol['confirmation_designs'][args.task][args.regime])
    design['seed'] = design.pop('seeds')[args.replicate]
    args.output.mkdir(parents=True, exist_ok=False)
    simulation = simulate_bams(args.task, args.output / 'synthetic', 10_000_000,
        design_override=design, depth_sigma=.25 if args.replicate == 0 else .5,
        quality=30 if args.replicate == 0 else 20)
    remapped = remap_synthetic_bams(args.output / 'synthetic', args.output / 'remapped', args.bwa, args.threads)
    write_json(args.output / 'prepared.json', {'task': args.task, 'regime': args.regime,
        'replicate': args.replicate, 'protocol_sha256': sha256(args.protocol),
        'source_sha256': sha256(Path(__file__)), 'simulation': simulation, 'remapped': remapped,
        'inference_performed': False, 'testing_data_access': False, 'final_data_access': False})


if __name__ == '__main__':
    main()
