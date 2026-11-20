"""Freeze/check the new Task14 stabilization policy; preserve old failures."""
import argparse
from datetime import datetime, timezone
import importlib.metadata as metadata
import json
from pathlib import Path
import sys

from ghist_lowcoverage.pipeline import audit_inputs, sha256


PROTOCOL_SHA = '851464f4ef60c675a2b59e90ec9176b080e4fc6da2997721b50a53b215d98c79'


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def verify(root, frozen_sha):
    frozen_path = root / 'task14_stabilization_freeze_v1.json'
    assert sha256(frozen_path) == frozen_sha
    frozen = json.loads(frozen_path.read_text())
    for name, expected in frozen['bindings'].items():
        assert sha256(root / name) == expected, name
    for name, expected in frozen['source_files'].items():
        assert sha256(root / 'source_v8' / name) == expected, name
    runtime = json.loads((root / 'runtime_task14_confirmation_v1.json').read_text())
    for name, expected in runtime['files'].items():
        assert sha256(Path(name)) == expected, name
    assert sys.version == runtime['python']
    assert all(metadata.version(name) == version for name, version in runtime['packages'].items())
    import ghist_lowcoverage.pipeline as pipeline
    assert Path(pipeline.__file__).resolve() == (root / 'source_v8/src/ghist_lowcoverage/pipeline.py').resolve()
    assert sha256(root / 'task14_stabilization_protocol_v1.json') == PROTOCOL_SHA
    return frozen


def input_binding(root, mode, regime, replicate):
    assert mode in ['development', 'confirmation']
    assert regime in range(3) and replicate in [0, 1]
    if mode == 'development':
        assert regime == 0
        case = root / f'full_confirmation_split_r0_p{replicate}_v1'
        plans = json.loads((root / 'task14_confirmation_freeze_v1.json').read_text())['confirmation_designs']
    else:
        case = root / f'full_stabilized_confirmation_split_r{regime}_p{replicate}_v1'
        plans = json.loads((root / 'task14_stabilization_protocol_v1.json').read_text())['confirmation_designs']['split_isolation']
    simulation_path = case / 'synthetic/simulation.json'
    simulation = json.loads(simulation_path.read_text())
    expected = dict(plans[regime])
    expected['seed'] = expected.pop('seeds')[replicate]
    assert simulation['design'] == expected and simulation['length'] == 10_000_000
    assert simulation['quality'] == (30 if replicate == 0 else 20)
    assert simulation['individual_depth_lognormal_sigma'] == (.25 if replicate == 0 else .5)
    assert not simulation['testing_data_access'] and not simulation['final_data_access']
    directory = case / 'remapped'
    audit = audit_inputs(directory, directory / 'reference.fasta', 'split_isolation', 10_000_000)
    indices = [Path(str(directory / 'reference.fasta') + '.fai')]
    indices += [Path(row['path'] + '.bai') for row in audit['samples'].values()]
    return case, {'audit': audit, 'index_sha256': {str(path): sha256(path) for path in indices},
                  'simulation_sha256': sha256(simulation_path)}


def freeze(root, manifest_sha):
    manifest = root / 'source_v8/source_manifest.json'
    assert sha256(manifest) == manifest_sha
    assert sha256(root / 'task14_stabilization_protocol_v1.json') == PROTOCOL_SHA
    assert sha256(root / 'runtime_task14_confirmation_v1.json') == '5cff417b7e3c0d28e386943adf6f14767825a8a41c0c1406023172f95b16d48c'
    assert sha256(root / 'task14_confirmation_freeze_v1.json') == '1f61090830f5d33cd62eb117488d32f906e670d8819add5f30edb1b26b0b11bc'
    old = json.loads((root / 'task14_confirmation_freeze_v1.json').read_text())
    new_sources = json.loads(manifest.read_text())
    for name, expected in old['source_files'].items():
        assert new_sources[name] == expected and sha256(root / 'source_v8' / name) == expected
    for name in ['external/angsd/angsd', 'external/bwa/bwa', 'external/winsfs-release-v0.7.0/winsfs',
                 'official_templates/split_isolation_low_coverage.submission.yaml', 'official_templates/provenance.json']:
        assert sha256(root / name) == old['bindings'][name]
    names = ['source_v8/source_manifest.json', 'task14_stabilization_protocol_v1.json',
        'runtime_task14_confirmation_v1.json', 'task14_confirmation_freeze_v1.json',
        'external/angsd/angsd', 'external/bwa/bwa', 'external/winsfs-release-v0.7.0/winsfs',
        'official_templates/split_isolation_low_coverage.submission.yaml', 'official_templates/provenance.json',
        'develop_task14_stabilization_v1.sbatch', 'prepare_task14_stabilized_confirmation_v1.sbatch',
        'infer_task14_stabilized_confirmation_v1.sbatch']
    record = {'frozen_at_utc': datetime.now(timezone.utc).isoformat(),
        'bindings': {name: sha256(root / name) for name in names},
        'source_files': json.loads(manifest.read_text()), 'protocol_sha256': PROTOCOL_SHA,
        'scope': 'One prospectively fixed stabilization option; development then conditional fresh confirmation.',
        'previous_confirmation_status': 'Failed; all original six outcomes remain retained.',
        'scientifically_ready': False}
    path = root / 'task14_stabilization_freeze_v1.json'
    write(path, record)
    verify(root, sha256(path))
    print(json.dumps({'freeze_sha256': sha256(path), 'source_files': len(record['source_files'])}))


def prepare_phase(root, frozen_sha, regime, replicate, phase):
    verify(root, frozen_sha)
    development = json.loads((root / 'task14_stabilized_development_gate_audit_v1.json').read_text())
    assert development['freeze_sha256'] == frozen_sha and development['all_registered_gates_passed']
    stem = f'split_r{regime}_p{replicate}_v1'
    record = {'freeze_sha256': frozen_sha, 'regime': regime, 'replicate': replicate, 'phase': phase,
        'runtime_source_tools_formatter_verified': True, 'verified_at_utc': datetime.now(timezone.utc).isoformat(),
        'development_gate_report_sha256': sha256(root / 'task14_stabilized_development_gate_audit_v1.json')}
    if phase == 'prepare_initial':
        write(root / f'full_stabilized_confirmation_{stem}.initial.json', record)
    else:
        case, inputs = input_binding(root, 'confirmation', regime, replicate)
        record['input_binding'] = inputs
        write(case / 'prepare_binding.json', record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--freeze', action='store_true')
    parser.add_argument('--source-manifest-sha')
    parser.add_argument('--freeze-sha')
    parser.add_argument('--regime', type=int, choices=range(3))
    parser.add_argument('--replicate', type=int, choices=[0, 1])
    parser.add_argument('--phase', choices=['prepare_initial', 'prepare_final'])
    args = parser.parse_args()
    assert args.root.name == 'lowcoverage-20261120'
    if args.freeze:
        freeze(args.root, args.source_manifest_sha)
    else:
        prepare_phase(args.root, args.freeze_sha, args.regime, args.replicate, args.phase)


if __name__ == '__main__':
    main()
