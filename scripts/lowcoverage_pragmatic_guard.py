"""Bind the separate pragmatic Task14 development policy without promoting old failures."""
import argparse
from datetime import datetime, timezone
import importlib.metadata as metadata
import json
from pathlib import Path
import sys

from ghist_lowcoverage.pipeline import audit_inputs, sha256


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def check_clock(root, stage):
    assert stage in ['development', 'final']
    plan = json.loads((root / 'task14_pragmatic_protocol_v1.json').read_text())
    now = datetime.now(timezone.utc)
    deadline = datetime.fromisoformat(plan['stage_latest_start_utc'][stage].replace('Z', '+00:00'))
    if now > deadline:
        raise RuntimeError(f'{stage} actual start exceeds registered deadline: {now.isoformat()} > {deadline.isoformat()}')
    return {'stage': stage, 'actual_start_utc': now.isoformat(), 'latest_start_utc': deadline.isoformat()}


def verify(root, frozen_sha):
    path = root / 'task14_pragmatic_freeze_v1.json'
    assert sha256(path) == frozen_sha
    frozen = json.loads(path.read_text())
    for name, expected in frozen['bindings'].items():
        assert sha256(root / name) == expected, name
    for name, expected in frozen['source_files'].items():
        assert sha256(root / 'source_v12' / name) == expected, name
    runtime = json.loads((root / 'runtime_task14_confirmation_v1.json').read_text())
    for name, expected in runtime['files'].items():
        assert sha256(Path(name)) == expected, name
    assert sys.version == runtime['python']
    assert all(metadata.version(name) == version for name, version in runtime['packages'].items())
    import ghist_lowcoverage.pipeline as pipeline
    import ghist_lowcoverage.stabilized as stabilized
    import ghist_lowcoverage.submission as submission
    import ghist_lowcoverage.winsfs as winsfs
    for module in [pipeline, stabilized, submission, winsfs]:
        assert Path(module.__file__).resolve() == (root / 'source_v12/src/ghist_lowcoverage' / Path(module.__file__).name).resolve()
    return frozen


def input_binding(root, mode, regime, replicate):
    assert mode == 'development' and regime == 0 and replicate in [0, 1]
    case = root / f'full_confirmation_split_r0_p{replicate}_v1'
    plans = json.loads((root / 'task14_confirmation_freeze_v1.json').read_text())['confirmation_designs']
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


def freeze(root, manifest_sha, protocol_sha):
    manifest = root / 'source_v12/source_manifest.json'
    assert sha256(manifest) == manifest_sha
    assert sha256(root / 'task14_pragmatic_protocol_v1.json') == protocol_sha
    assert sha256(root / 'task14_extension_freeze_v3.json') == 'ed4dc55bd3729d04a8dfa09dc0fabb16b24b36f7301cb24870f5b85d34e89047'
    old = json.loads((root / 'task14_extension_freeze_v3.json').read_text())
    sources = json.loads(manifest.read_text())
    for name, expected in old['source_files'].items():
        assert sources[name] == expected and sha256(root / 'source_v12' / name) == expected
    assert not list(root.glob('full_extended_confirmation*'))
    assert not list(root.glob('full_stabilized_confirmation*'))
    assert not list(root.glob('window_extended_development*'))
    names = ['source_v12/source_manifest.json', 'task14_pragmatic_protocol_v1.json',
        'task14_pragmatic_authorization_v1.json', 'task14_extension_strict_path_closed_v1.json',
        'runtime_task14_confirmation_v1.json', 'task14_confirmation_freeze_v1.json',
        'task14_confirmation_gate_audit_v1.json', 'task14_stabilization_freeze_v1.json',
        'task14_stabilized_development_gate_audit_v1.json', 'task14_stabilization_v1_no_fresh_dispatch.json',
        'task14_extension_freeze_v1.json', 'task14_extension_freeze_v2.json', 'task14_extension_freeze_v3.json',
        'task14_extension_v1_dispatch_failure.json', 'task14_extension_v2_dispatch_failure.json',
        'external/angsd/angsd', 'external/bwa/bwa', 'external/winsfs-release-v0.7.0/winsfs',
        'official_templates/split_isolation_low_coverage.submission.yaml', 'official_templates/provenance.json',
        'develop_task14_pragmatic_v1.sbatch']
    for name in ['external/angsd/angsd', 'external/bwa/bwa', 'external/winsfs-release-v0.7.0/winsfs',
                 'official_templates/split_isolation_low_coverage.submission.yaml', 'official_templates/provenance.json']:
        assert sha256(root / name) == old['bindings'][name]
    record = {'frozen_at_utc': datetime.now(timezone.utc).isoformat(),
        'bindings': {name: sha256(root / name) for name in names}, 'source_files': sources,
        'protocol_sha256': protocol_sha, 'scope': 'Separate pragmatic development, both existing Q30/Q20 cases, no fresh generation.',
        'strict_qualification': False, 'scientifically_ready': False,
        'previous_failures': 'All original six confirmation, cap200 stopping and strict cap500 clock failures remain retained.'}
    path = root / 'task14_pragmatic_freeze_v1.json'
    write(path, record)
    verify(root, sha256(path))
    print(json.dumps({'freeze_sha256': sha256(path), 'source_files': len(sources)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--source-manifest-sha', required=True)
    parser.add_argument('--protocol-sha', required=True)
    args = parser.parse_args()
    assert args.root.name == 'lowcoverage-20261120'
    freeze(args.root, args.source_manifest_sha, args.protocol_sha)


if __name__ == '__main__':
    main()
