"""Prepare a new frozen pair from archived script bytes; never invoke SLiM."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REFERENCE_HASHES = {
    'reference/protocol.json': '52fb8147aba90613e08823c84dc1cacfa2f440f0b4b162e1db1cebfaec9f4f33',
    'reference/interval20.slim': 'f591fdc135e8c9c17076527475e41423d826c346b0e143ce34c9de40b7d7e26e',
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    for name, expected in REFERENCE_HASHES.items():
        if digest(HERE / name) != expected:
            raise ValueError(f'Archived reference differs: {name}')
    old = json.loads((HERE / 'reference/protocol.json').read_text())
    runtime = json.loads((HERE / 'runtime_reference.json').read_text())
    left = (HERE / 'reference/interval20.slim').read_text()
    needle = 'initializeTreeSeq(simplificationInterval=20, timeUnit="generations");'
    replacement = 'initializeTreeSeq(simplificationInterval=100, timeUnit="generations");'
    if left.count(needle) != 1:
        raise ValueError('Expected exactly one fixed interval20 declaration')
    scripts = {'interval20.slim': left, 'interval100.slim': left.replace(needle, replacement)}
    protocol_path = HERE / 'protocol.json'
    if protocol_path.exists():
        raise FileExistsError(protocol_path)
    for name, text in scripts.items():
        path = HERE / name
        if path.exists() and path.read_text() != text:
            raise ValueError(f'Existing prepared script differs: {name}')
        if not path.exists():
            path.write_text(text)
    names = ['prepare.py', 'run.py', 'run.sbatch', 'runtime_reference.json',
             *REFERENCE_HASHES, *scripts]
    protocol = {
        'schema': 'paired_slim_fixed_interval20_vs100_software_fixture_v1',
        'seed': old['seed'], 'biology': old['biology'],
        'arms': {'interval20': 'explicit fixed simplificationInterval20 with default NULL ratio',
                 'interval100': 'explicit fixed simplificationInterval100 with default NULL ratio'},
        'comparisons': old['comparisons'], 'resource_limits': old['resource_limits'],
        'script_difference': 'Exactly one initializeTreeSeq fixed-interval argument; all other bytes identical.',
        'source_fixture_protocol_sha256': REFERENCE_HASHES['reference/protocol.json'],
        'source_sha256': {name: digest(HERE / name) for name in names},
        'engine_path': old['engine_path'], 'engine_sha256': old['engine_sha256'],
        'slim_path': old['slim_path'], 'slim_sha256': old['slim_sha256'],
        'versions': runtime['versions'], 'python_version': runtime['python']['version'],
        'python_executable_sha256': runtime['python_executable_sha256'],
        'scientifically_validated': False,
        'scope': 'Bookkeeping equivalence software fixture only; no full-chromosome resource, accuracy or recovery claim.',
    }
    protocol_path.write_text(json.dumps(protocol, indent=2) + '\n')
    print(json.dumps({'protocol_sha256': digest(protocol_path), 'source_sha256': protocol['source_sha256']}))


if __name__ == '__main__':
    main()
