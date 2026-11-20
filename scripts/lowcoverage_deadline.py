"""Run the frozen Task14 stabilization option on independent synthetic data."""
import argparse
import json
from pathlib import Path
import shutil

from ghist_lowcoverage.pipeline import sha256
from ghist_lowcoverage.stabilized import run_stabilized_comparison
from ghist_lowcoverage.submission import format_submission_draft
from ghist_lowcoverage.winsfs import make_safs
from scripts.lowcoverage_deadline_guard import check_clock, input_binding, verify, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--freeze-sha', required=True)
    parser.add_argument('--mode', choices=['development', 'confirmation'], required=True)
    parser.add_argument('--regime', type=int, choices=range(3), required=True)
    parser.add_argument('--replicate', type=int, choices=[0, 1], required=True)
    args = parser.parse_args()
    root = args.root
    assert root.name == 'lowcoverage-20261120'
    clock_record = check_clock(root, args.mode)
    verify(root, args.freeze_sha)
    case, before = input_binding(root, args.mode, args.regime, args.replicate)
    stem = f'split_r{args.regime}_p{args.replicate}_v1'
    output = root / f'window_extended_{args.mode}_{stem}'
    output.mkdir(parents=True, exist_ok=False)
    record = {'freeze_sha256': args.freeze_sha, 'mode': args.mode, 'regime': args.regime,
        'replicate': args.replicate, 'input_binding': before, 'runtime_source_tools_formatter_verified': True,
        'wrapper_sha256': sha256(Path(__file__)), 'clock_guard': clock_record}
    if args.mode == 'confirmation':
        prepared = json.loads((case / 'prepare_binding.json').read_text())
        assert prepared['input_binding'] == before and prepared['freeze_sha256'] == args.freeze_sha
    write(output / 'initial_binding.json', record)
    if args.mode == 'development':
        previous = root / f'window_confirmation_{stem}'
        old_binding = json.loads((previous / 'execution_binding.json').read_text())
        assert old_binding['input_binding'] == before
        audit = json.loads((previous / 'safs/inputs.json').read_text())
        assert audit == before['audit']
        old_window = json.loads((previous / 'window_comparison/protocol.json').read_text())
        for row in old_window['saf_files']:
            assert sha256(Path(row['path'])) == row['sha256']
        indices = [previous / f'safs/population_{i}.saf.idx' for i in range(2)]
        (output / 'safs').mkdir()
        for name in ['inputs.json', 'protocol.json']:
            shutil.copyfile(previous / 'safs' / name, output / 'safs' / name)
    else:
        audit, indices = make_safs(case / 'remapped', case / 'remapped/reference.fasta', 'split_isolation',
            output / 'safs', root / 'external/angsd/angsd', threads=4, expected_length=10_000_000)
    comparison = run_stabilized_comparison(audit, indices, output / 'window_comparison',
        root / 'external/winsfs-release-v0.7.0/winsfs', threads=4, max_epochs=500)
    verify(root, args.freeze_sha)
    _, after = input_binding(root, args.mode, args.regime, args.replicate)
    assert before == after
    for row in json.loads((output / 'window_comparison/protocol.json').read_text())['saf_files']:
        assert sha256(Path(row['path'])) == row['sha256']
    simulation = json.loads((case / 'synthetic/simulation.json').read_text())
    errors = {label: {key: estimate[key] / value - 1 for key, value in simulation['truth'].items()}
              for label, estimate in comparison['estimates'].items()}
    write(output / 'scored.json', {'truth': simulation['truth'], 'relative_errors': errors,
        'comparison': comparison, 'simulation_sha256': before['simulation_sha256'],
        'freeze_sha256': args.freeze_sha, 'testing_data_access': False, 'final_data_access': False,
        'scientifically_ready': False})
    draft = format_submission_draft('split_isolation', comparison['estimates']['seed1'],
        root / 'official_templates/split_isolation_low_coverage.submission.yaml')
    import hashlib
    import yaml
    assert yaml.safe_load(draft)['parameters'] == comparison['estimates']['seed1']
    record.update({'source_runtime_input_pre_post_verified': True, 'synthetic_formatter_schema_passed': True,
        'synthetic_formatted_string_sha256': hashlib.sha256(draft.encode()).hexdigest(),
        'initial_binding_sha256': sha256(output / 'initial_binding.json'),
        'output_bindings': {name: sha256(output / name) for name in [
            'scored.json', 'safs/protocol.json', 'safs/inputs.json', 'window_comparison/protocol.json']}})
    write(output / 'execution_binding.json', record)


if __name__ == '__main__':
    main()
