"""Run a separately registered window-EM comparison on synthetic data only."""
import argparse
import json
from pathlib import Path

from ghist_lowcoverage.pipeline import write_json, sha256
from ghist_lowcoverage.winsfs import make_safs, run_window_comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--winsfs', type=Path, required=True)
    parser.add_argument('--angsd', type=Path)
    parser.add_argument('--existing-inference', type=Path)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--max-epochs', type=int, default=20)
    args = parser.parse_args()
    simulation_path = args.case / 'synthetic/simulation.json'
    simulation = json.loads(simulation_path.read_text())
    if simulation['testing_data_access'] or simulation['final_data_access']:
        raise ValueError('independent synthetic data only')
    args.output.mkdir(parents=True, exist_ok=False)
    if args.existing_inference:
        audit = json.loads((args.existing_inference/'inputs.json').read_text())
        for sample in audit['samples'].values():
            if sha256(Path(sample['path'])) != sample['sha256']:
                raise ValueError('input BAM changed')
        if sha256(Path(audit['reference'])) != audit['reference_sha256']:
            raise ValueError('input reference changed')
        indices = [args.existing_inference/f'population_{i}.saf.idx' for i in range(len(audit['population_order']))]
    else:
        if args.angsd is None:
            parser.error('--angsd required when making new SAFs')
        directory = args.case / ('remapped' if (args.case/'remapped').is_dir() else 'synthetic')
        audit, indices = make_safs(directory, directory/'reference.fasta', simulation['task'], args.output/'safs',
                                  args.angsd, threads=args.threads, expected_length=simulation['length'])
    comparison = run_window_comparison(audit, indices, args.output/'window_comparison', args.winsfs,
                                      threads=args.threads, max_epochs=args.max_epochs)
    errors = {label: {k: estimate[k]/value-1 for k,value in simulation['truth'].items()}
              for label, estimate in comparison['estimates'].items()}
    write_json(args.output/'scored.json', {'truth': simulation['truth'], 'relative_errors': errors,
        'comparison': comparison, 'simulation_sha256': sha256(simulation_path),
        'source_sha256': sha256(Path(__file__)),
        'input_likelihood_protocol_sha256': sha256((args.existing_inference or args.output/'safs')/'protocol.json'),
        'testing_data_access': False, 'final_data_access': False,
        'scientifically_ready': False})


if __name__ == '__main__':
    main()
