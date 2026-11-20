"""Read-only, posthoc resource/source collection; never read prediction metrics."""
import datetime
import hashlib
import importlib.metadata
import importlib.util
import json
import math
from pathlib import Path
import platform
import re

HERE = Path(__file__).resolve().parent
EXPECTED_PROTOCOL = '9004cf2349334767c1bb9777d3ad9ef7063b961876e8dbd92ac7222dc1901307'


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    protocol = json.loads((HERE / 'protocol.json').read_text())
    locks = {HERE / name: sha for name, sha in protocol['source_sha256'].items()}
    locks.update({Path(name): sha for name, sha in protocol['external_sha256'].items()})
    locks[HERE / 'protocol.json'] = EXPECTED_PROTOCOL
    mismatches = [str(path) for path, sha in locks.items() if digest(path) != sha]
    versions = {name: importlib.metadata.version(name) for name in protocol['versions']}
    if versions != protocol['versions'] or platform.python_version() != protocol['python']['version']:
        mismatches.append('runtime_versions')
    output = Path(protocol['output_directory'])
    diagnostics = output / 'diagnostics'
    resource_bytes = (diagnostics / 'resources.jsonl').read_bytes()
    rows = [json.loads(line) for line in resource_bytes.splitlines()]
    groups = [next(values for name, values in row['cgroups'].items()
                   if name.endswith('job_10857763')) for row in rows]
    timestamps = [row['unix_time'] for row in rows]
    gaps = [right - left for left, right in zip(timestamps, timestamps[1:])]
    if not all(math.isfinite(value) for value in timestamps) or any(gap < 0 for gap in gaps):
        raise ValueError('Nonfinite or nonmonotone resource timestamps')
    phase_bytes = (diagnostics / 'phase.tsv').read_bytes()
    lines = phase_bytes.decode().splitlines()
    if lines[0] != 'tick\tdiploid_population_size\tsaves\trestores':
        raise ValueError('Unexpected phase schema')
    phases = [list(map(int, line.split('\t'))) for line in lines[1:]]
    if any(len(row) != 4 or min(row) < 0 for row in phases):
        raise ValueError('Invalid phase record')
    midpoint_path = HERE / 'midpoint_observation.json'
    midpoint_checks = {}
    if midpoint_path.exists():
        midpoint = json.loads(midpoint_path.read_text())
        midpoint_checks = {
            'resources': hashlib.sha256(resource_bytes[:midpoint['resources_prefix_bytes']]).hexdigest()
            == midpoint['resources_prefix_sha256'],
            'phase': hashlib.sha256(phase_bytes[:midpoint['phase_prefix_bytes']]).hexdigest()
            == midpoint['phase_prefix_sha256'],
        }
    spec = importlib.util.spec_from_file_location('frozen_adapter', HERE / 'adapter.py')
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    original = (diagnostics / 'generated_original.slim').read_text()
    transformed = adapter.transform_script(original, diagnostics / 'raw_slim.trees', diagnostics / 'phase.tsv')
    script_matches = transformed == (diagnostics / 'generated_interval100.slim').read_text()
    times = [float(value) for value in re.search(r'defineConstant\("_T", c\(([^)]+)\)', original)[1].split(',')]
    sizes = [float(value) for value in re.search(r'_N = array\(c\(\s*//[^\n]*\n\s*c\(([^)]+)\)', original)[1].split(',')]
    rates = [float(value) for value in re.search(r'defineConstant\("growth_rates", array\(c\(\s*//[^\n]*\n\s*c\(([^)]+)\)', original)[1].split(',')]
    constants = {name: float(re.search(r'defineConstant\("' + name + r'", ([^)]+)\)', original)[1])
                 for name in ('Q', 'burn_in', 'generation_time')}
    if constants != {'Q': 1, 'burn_in': 1, 'generation_time': 30}:
        raise ValueError('Frozen physical schedule assumptions differ from generated source')
    end_tick = int(max(times) / constants['generation_time'] + 1 + round(sizes[0]))
    boundaries = [end_tick - int(value / constants['generation_time']) for value in times]
    schedule = []
    for index in range(1, len(times)):
        start, end = boundaries[index - 1] + 1, boundaries[index] - 1
        rate = rates[index] * (end - start + 1) / (end - start)
        schedule.append({'phase': index, 'start_tick': start, 'end_tick': end,
                         'initial_target_diploids': sizes[index],
                         'approximate_last_scheduled_target': round(sizes[index] * math.exp(rate * (end - boundaries[index - 1])))})
    report = {
        'schema': 'task5_interval100_terminal_resource_collection_v1',
        'collected_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'collector_sha256': digest(Path(__file__)), 'protocol_sha256': digest(HERE / 'protocol.json'),
        'frozen_source_count': len(protocol['source_sha256']), 'external_lock_count': len(protocol['external_sha256']),
        'post_run_mismatches': mismatches, 'runtime_versions': versions,
        'generated_script_exact_transformation': script_matches, 'midpoint_prefix_checks': midpoint_checks,
        'midpoint_available': midpoint_path.exists(),
        'adapter_status': json.loads((diagnostics / 'adapter_status.json').read_text()),
        'resource_trace': {'records': len(rows), 'first_unix_time': timestamps[0], 'last_unix_time': timestamps[-1],
                           'duration_seconds': timestamps[-1] - timestamps[0], 'maximum_sample_gap_seconds': max(gaps),
                           'stages': sorted(set(row['stage'] for row in rows)),
                           'last_usage_bytes': int(groups[-1]['memory.usage_in_bytes']),
                           'maximum_recorded_usage_bytes': max(int(group['memory.usage_in_bytes']) for group in groups),
                           'maximum_recorded_high_water_bytes': max(int(group['memory.max_usage_in_bytes']) for group in groups),
                           'limits_bytes': sorted(set(int(group['memory.limit_in_bytes']) for group in groups)),
                           'maximum_failcnt': max(int(group['memory.failcnt']) for group in groups),
                           'oom_controls': sorted(set(group['memory.oom_control'] for group in groups)),
                           'monitor_failure_file': (diagnostics / 'monitor_failure.json').exists()},
        'phase_trace': {'records': len(phases), 'first': phases[0], 'last': phases[-1],
                        'maximum_logged_tick': max(row[0] for row in phases),
                        'schedule_boundaries_reached': [tick for tick in boundaries
                                                        if any(row[0] >= tick for row in phases)],
                        'source_end_tick_reached': any(row[0] >= end_tick for row in phases),
                        'maximum_logged_diploids': max(row[1] for row in phases),
                        'maximum_saves': max(row[2] for row in phases), 'maximum_restores': max(row[3] for row in phases)},
        'source_schedule': {'constants': constants, 'end_tick': end_tick, 'epoch_boundaries': boundaries, 'catalog_sizes': sizes,
                            'growth_phases': schedule, 'interpretation': 'Source-derived targets only; observed phase coverage must be assessed separately from the phase trace.'},
        'output_files': {str(path.relative_to(output)): {'bytes': path.stat().st_size, 'sha256': digest(path)}
                         for path in sorted(output.rglob('*')) if path.is_file()},
        'scope': 'Resource/format/source evidence only. No prediction or accuracy content read; byte hashes bind retained outputs.',
        'scientifically_validated': False,
    }
    (HERE / 'terminal_collection.json').write_text(json.dumps(report, indent=2) + '\n')
    if mismatches or not script_matches or not all(midpoint_checks.values()):
        raise ValueError('Terminal source/evidence verification failed; report retained')
    print(json.dumps({'source_checks': 'PASSED', 'resource_trace': report['resource_trace'],
                      'phase_trace': report['phase_trace']}))


if __name__ == '__main__':
    main()
