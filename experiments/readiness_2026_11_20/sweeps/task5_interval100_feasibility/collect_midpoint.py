"""Capture complete append-only monitor prefixes; do not inspect accuracy outputs."""
import datetime
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
EXPECTED_PROTOCOL = '9004cf2349334767c1bb9777d3ad9ef7063b961876e8dbd92ac7222dc1901307'


def prefix(path):
    raw = Path(path).read_bytes()
    last = raw.rfind(b'\n')
    if last < 0:
        raise ValueError('No complete monitor records')
    return raw[:last + 1]


def main():
    protocol_bytes = (HERE / 'protocol.json').read_bytes()
    if hashlib.sha256(protocol_bytes).hexdigest() != EXPECTED_PROTOCOL:
        raise ValueError('Unexpected frozen protocol')
    protocol = json.loads(protocol_bytes)
    diagnostics = Path(protocol['output_directory']) / 'diagnostics'
    resources = prefix(diagnostics / 'resources.jsonl')
    phases = prefix(diagnostics / 'phase.tsv')
    rows = [json.loads(line) for line in resources.splitlines()]
    scheduler = subprocess.run(['squeue', '-h', '-j', '10857763', '-o', '%T|%M|%l'],
                               check=True, capture_output=True, text=True).stdout.strip()
    if not scheduler.startswith('RUNNING|'):
        raise ValueError('Midpoint requires a running allocation')
    report = {
        'schema': 'task5_interval100_midpoint_observation_v1',
        'observed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'job_id': '10857763', 'state': 'RUNNING_PARTIAL', 'scheduler_snapshot': scheduler,
        'protocol_sha256': EXPECTED_PROTOCOL, 'resource_observations': len(rows),
        'resources_prefix_bytes': len(resources),
        'resources_prefix_sha256': hashlib.sha256(resources).hexdigest(),
        'phase_prefix_bytes': len(phases), 'phase_prefix_sha256': hashlib.sha256(phases).hexdigest(),
        'last_resource': rows[-1], 'last_phase': phases.decode().splitlines()[-1],
        'monitor_failure': (diagnostics / 'monitor_failure.json').exists(),
        'interpretation': 'Partial resource observation only; no full-run memory bound, accuracy or recovery claim. Original 45-minute cap remains unchanged.',
    }
    with (HERE / 'midpoint_observation.json').open('x') as handle:
        handle.write(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
