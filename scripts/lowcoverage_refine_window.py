"""Bounded full-batch Task13 development diagnostic from three window starts."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from ghist_lowcoverage.pipeline import audit_inputs, fit_spectrum, realsfs_diagnostics, sha256, write_json
from ghist_lowcoverage.winsfs import read_winsfs


def initialization_counts(values):
    values = np.asarray(values, dtype=float)
    if values.shape != (41,) or np.any(~np.isfinite(values)) or np.any(values < 0) or np.any(values[21:] != 0):
        raise ValueError('Task13 globally folded 41-bin count spectrum required')
    result = values.copy()
    result[:21] = np.maximum(result[:21], 1e-8)
    return result


def input_binding(root):
    audit = json.loads((root / 'window_full_bottleneck_v7/safs/inputs.json').read_text())
    directory = Path(audit['reference']).parent
    assert audit_inputs(directory, Path(audit['reference']), 'bottleneck', 10_000_000) == audit
    source = root / 'source_v7/source_manifest.json'
    assert sha256(source) == '43a608f4bde3334fc8239f48b4dd4834ee3a87ee547acb15b379bfc4838e21f4'
    for name, expected in json.loads(source.read_text()).items():
        assert sha256(root / 'source_v7' / name) == expected
    window = json.loads((root / 'window_full_bottleneck_v7/window_comparison/protocol.json').read_text())
    for row in window['saf_files']:
        assert sha256(Path(row['path'])) == row['sha256']
    runtime_path = root / 'runtime_task14_confirmation_v1.json'
    assert sha256(runtime_path) == '5cff417b7e3c0d28e386943adf6f14767825a8a41c0c1406023172f95b16d48c'
    runtime = json.loads(runtime_path.read_text())
    for name, expected in runtime['files'].items():
        assert sha256(Path(name)) == expected
    assert sha256(root / 'external/angsd/misc/realSFS') == '3292b27066892d303a2d2a0fa8ca4941a21198ed24ce8e2ea409c48c2e7e0dd2'
    return {'input_audit': audit, 'saf_files': window['saf_files'], 'source_manifest_sha256': sha256(source),
            'runtime_manifest_sha256': sha256(runtime_path), 'script_sha256': sha256(Path(__file__)),
            'realSFS_sha256': sha256(root / 'external/angsd/misc/realSFS')}


def comparison(left, right, first, second):
    return {'frequency_L1': float(np.abs(left / left.sum() - right / right.sum()).sum()),
            'segregating_relative_L1': float(np.abs(left[1:] - right[1:]).sum() / left[1:].sum()),
            'target_relative_change': {key: second['targets'][key] / value - 1 for key, value in first['targets'].items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--label', choices=['seed1', 'seed2', 'continued'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.root.name != 'lowcoverage-20261120':
        raise ValueError('isolated synthetic experiment root required')
    registered = args.root / 'realsfs_refinement_development_v1.json'
    assert sha256(registered) == 'd25f94eae94fd168c6bb845f295ceb6eba2b97d7af10e16ad669c0f9667c7120'
    before = input_binding(args.root)
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output / 'initial_binding.json', before)
    starting_path = args.root / f'window_full_bottleneck_v7/window_comparison/{args.label}.folded.sfs'
    original = read_winsfs(starting_path, [40])
    starting = initialization_counts(original)
    initial = args.output / 'initialization.sfs'
    np.savetxt(initial, starting[None, :], fmt='%.17g')
    protocol = {'registered_protocol_sha256': sha256(registered), 'starting_sfs_sha256': sha256(starting_path),
                'initialization_sha256': sha256(initial), 'floored_bins': np.flatnonzero(starting != original).tolist(),
                'initialization_count_added': float((starting - original).sum()),
                'script_sha256': sha256(Path(__file__)), 'commands': [], 'scientifically_ready': False}
    spectra, estimates = {}, {}
    for phase in ['refined', 'continued']:
        result_path = args.output / f'{phase}.sfs'
        start = initial if phase == 'refined' else args.output / 'refined.sfs'
        command = [str(args.root / 'external/angsd/misc/realSFS'),
            str(args.root / 'window_full_bottleneck_v7/safs/population_0.saf.idx'),
            '-fold', '1', '-P', '4', '-maxIter', '250', '-tole', '1e-8', '-seed', '20261120', '-sfs', str(start)]
        protocol['commands'].append(command)
        write_json(args.output / 'protocol.json', protocol)
        log_path = args.output / f'{phase}.log'
        with result_path.open('w') as result, log_path.open('w') as log:
            subprocess.run(command, stdout=result, stderr=log, check=True)
        values = np.loadtxt(result_path)
        if values.shape != (41,) or not 0 < values.sum() <= 10_000_001:
            raise ValueError('invalid retained exposure or dimensions')
        estimate = fit_spectrum(values, [40], 'bottleneck')
        estimate['realsfs_diagnostics'] = realsfs_diagnostics(log_path.read_text())
        estimate['callable_fraction'] = float(values.sum() / 10_000_000)
        write_json(args.output / f'{phase}.estimate.json', estimate)
        spectra[phase], estimates[phase] = values, estimate
    after = input_binding(args.root)
    assert before == after, 'source or inputs changed during execution'
    simulation_path = args.root / 'full_bottleneck/synthetic/simulation.json'
    simulation = json.loads(simulation_path.read_text())
    assert not simulation['testing_data_access'] and not simulation['final_data_access']
    scores = {phase: {key: estimate['targets'][key] / value - 1 for key, value in simulation['truth'].items()}
              for phase, estimate in estimates.items()}
    write_json(args.output / 'scored.json', {'starting_label': args.label, 'truth': simulation['truth'],
        'estimates': {phase: estimate['targets'] for phase, estimate in estimates.items()},
        'relative_errors': scores, 'continuation_comparison': comparison(spectra['refined'], spectra['continued'],
            estimates['refined'], estimates['continued']), 'source_and_input_pre_post_verification_passed': True,
        'initial_binding_sha256': sha256(args.output / 'initial_binding.json'),
        'simulation_sha256': sha256(simulation_path), 'protocol_sha256': sha256(args.output / 'protocol.json'),
        'scientifically_ready': False, 'scope': 'Previously inspected Task13 development only; v7 failures remain retained.'})


if __name__ == '__main__':
    main()
