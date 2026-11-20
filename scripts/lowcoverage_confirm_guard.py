"""Freeze and verify the complete Task14 synthetic confirmation execution."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import sys
import sysconfig


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def write(path, value):
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def freeze(root):
    manifest = root / 'source_v7/source_manifest.json'
    assert sha(manifest) == '43a608f4bde3334fc8239f48b4dd4834ee3a87ee547acb15b379bfc4838e21f4'
    environment = json.loads((root / 'environment_v7.json').read_text())
    runtime = {'python': sys.version, 'packages': {}, 'files': {}}
    assert runtime['python'] == environment['python']
    for name, version in environment['packages'].items():
        distribution = metadata.distribution(name)
        assert distribution.version == version, name
        runtime['packages'][name] = version
        for relative in distribution.files or []:
            path = Path(distribution.locate_file(relative)).resolve()
            if path.is_file() and '__pycache__' not in path.parts and path.suffix not in ('.pyc', '.pyo'):
                runtime['files'][str(path)] = sha(path)
    for path in [Path(sys.executable).resolve(),
                 Path(sysconfig.get_config_var('LIBDIR')) / sysconfig.get_config_var('LDLIBRARY')]:
        if path.is_file():
            runtime['files'][str(path)] = sha(path)
    runtime_path = root / 'runtime_task14_confirmation_v1.json'
    write(runtime_path, runtime)
    names = ['source_v7/source_manifest.json', 'validation_acceptance_v1.json',
             'winsfs_full_development_v1.json', 'runtime_task14_confirmation_v1.json',
             'external/angsd/angsd', 'external/bwa/bwa', 'external/winsfs-release-v0.7.0/winsfs',
             'official_templates/split_isolation_low_coverage.submission.yaml',
             'official_templates/provenance.json', 'lowcoverage_confirm_guard_v1.py',
             'prepare_task14_confirmation_v1.sbatch', 'infer_task14_confirmation_v1.sbatch']
    plan = json.loads((root / 'validation_acceptance_v1.json').read_text())
    record = {
        'frozen_at_utc': datetime.now(timezone.utc).isoformat(), 'task': 'split_isolation',
        'scope': 'Six prespecified independent synthetic Task14 cohorts only; no challenge input or Task13 confirmation seeds.',
        'bindings': {name: sha(root / name) for name in names},
        'source_files': json.loads(manifest.read_text()),
        'confirmation_designs': plan['confirmation_designs']['split_isolation'],
        'noise_profiles': {'replicate0': {'quality': 30, 'depth_sigma': .25},
                           'replicate1': {'quality': 20, 'depth_sigma': .5}},
        'method': 'Unchanged source_v7 GL2/SAF1 minMapQ30/minQ20; winsfs0.7.0 blocks500/window100/tolerance1e-4/max_epochs50; seeds761001/761002 and initialized continuation; global folding; v5 eight starts each L-BFGS-B and Powell.',
        'primary_estimate': 'seed1; seed2 and continuation are numerical checks, not extra independent cohorts.',
        'numerical_gates': {'callable_fraction_min': .99, 'each_run_stop_before_epoch': 50,
            'all_three_optimizer_agreements': True, 'target_relative_change_max_exclusive': .01,
            'full_frequency_L1_max_exclusive': .001, 'segregating_relative_L1_max_exclusive': .05},
        'scientific_gates': {'all_targets_all_cohorts_absolute_relative_error_max': .30,
            'median_absolute_relative_error_max': .20,
            'median_definition': 'All 18 primary seed1 target errors across six cohorts, with each noise stratum additionally reported separately.',
            'extra_conservative_check': 'Also require <=30% recovery for every target in alternate-seed/continued fits.'},
        'failure_rule': 'Preserve every failure/timeout and all six cohorts. No refits, seed replacement, tolerance change, or dropping Q20/high-variation strata under this confirmation label.',
        'historical_gates': 'Earlier realSFS criteria and Task13 failures remain unpassed; this is separately registered Task14 window-estimator confirmation.',
        'budget': 'At most2 concurrent4CPU jobs,5632M/CPU; six preparation jobs capped30min each then six inference jobs capped45min each. Worst scheduler-free bound3h45m; expected30-50min.',
        'release_rule': 'No scientific readiness until all six cohorts, all numerical/scientific gates, provenance checks, and independent audit pass; scope limited to tested model/noise domain.',
        'scientifically_ready': False,
    }
    target = root / 'task14_confirmation_freeze_v1.json'
    write(target, record)
    print(json.dumps({'freeze': str(target), 'sha256': sha(target), 'runtime_files': len(runtime['files'])}))


def verify(root, frozen_sha, regime, replicate, phase):
    from ghist_lowcoverage.pipeline import audit_inputs

    frozen_path = root / 'task14_confirmation_freeze_v1.json'
    assert sha(frozen_path) == frozen_sha
    frozen = json.loads(frozen_path.read_text())
    for name, expected in frozen['bindings'].items():
        assert sha(root / name) == expected, name
    for name, expected in frozen['source_files'].items():
        assert sha(root / 'source_v7' / name) == expected, name
    runtime = json.loads((root / 'runtime_task14_confirmation_v1.json').read_text())
    for name, expected in runtime['files'].items():
        assert sha(Path(name)) == expected, name
    assert sys.version == runtime['python']
    assert all(metadata.version(name) == expected for name, expected in runtime['packages'].items())
    stem = f'split_r{regime}_p{replicate}_v1'
    case = root / f'full_confirmation_{stem}'
    output = root / f'window_confirmation_{stem}'
    record = {'freeze_sha256': frozen_sha, 'runtime_source_tools_formatter_verified': True,
              'regime': regime, 'replicate': replicate, 'phase': phase,
              'verified_at_utc': datetime.now(timezone.utc).isoformat()}
    if phase != 'prepare_initial':
        simulation = json.loads((case / 'synthetic/simulation.json').read_text())
        design = dict(frozen['confirmation_designs'][regime])
        design['seed'] = design.pop('seeds')[replicate]
        assert simulation['design'] == design and simulation['length'] == 10_000_000
        noise = frozen['noise_profiles'][f'replicate{replicate}']
        assert simulation['quality'] == noise['quality']
        assert simulation['individual_depth_lognormal_sigma'] == noise['depth_sigma']
        assert not simulation['testing_data_access'] and not simulation['final_data_access']
        directory = case / 'remapped'
        audit = audit_inputs(directory, directory / 'reference.fasta', 'split_isolation', 10_000_000)
        indices = [Path(str(directory / 'reference.fasta') + '.fai')]
        indices += [Path(row['path'] + '.bai') for row in audit['samples'].values()]
        inputs = {'audit': audit, 'index_sha256': {str(path): sha(path) for path in indices},
                  'simulation_sha256': sha(case / 'synthetic/simulation.json')}
        record['input_binding'] = inputs
        if phase.startswith('infer'):
            initial = json.loads((case / 'prepare_binding.json').read_text())
            assert initial['input_binding'] == inputs
            assert initial['freeze_sha256'] == frozen_sha
        if phase == 'infer_final':
            import yaml
            from ghist_lowcoverage.submission import format_submission_draft

            assert json.loads((output / 'safs/inputs.json').read_text()) == audit
            targets = json.loads((output / 'scored.json').read_text())['comparison']['estimates']['seed1']
            draft = format_submission_draft('split_isolation', targets,
                root / 'official_templates/split_isolation_low_coverage.submission.yaml')
            assert yaml.safe_load(draft)['parameters'] == targets
            record['synthetic_formatter_schema_passed'] = True
            record['synthetic_formatted_string_sha256'] = hashlib.sha256(draft.encode()).hexdigest()
            record['output_bindings'] = {name: sha(output / name) for name in [
                'scored.json', 'safs/protocol.json', 'window_comparison/protocol.json']}
    destinations = {'prepare_initial': root / f'full_confirmation_{stem}.initial.json',
                    'prepare_final': case / 'prepare_binding.json',
                    'infer_initial': root / f'window_confirmation_{stem}.initial.json',
                    'infer_final': output / 'execution_binding.json'}
    write(destinations[phase], record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--freeze', action='store_true')
    parser.add_argument('--freeze-sha')
    parser.add_argument('--regime', type=int, choices=range(3))
    parser.add_argument('--replicate', type=int, choices=[0, 1])
    parser.add_argument('--phase', choices=['prepare_initial', 'prepare_final', 'infer_initial', 'infer_final'])
    args = parser.parse_args()
    if args.root.name != 'lowcoverage-20261120':
        raise ValueError('isolated synthetic experiment root required')
    if args.freeze:
        freeze(args.root)
    else:
        assert args.freeze_sha and args.phase and args.regime is not None and args.replicate is not None
        verify(args.root, args.freeze_sha, args.regime, args.replicate, args.phase)


if __name__ == '__main__':
    main()
