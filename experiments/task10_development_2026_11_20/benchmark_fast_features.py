"""Bounded equivalence and runtime checks, using training cohorts only."""
import argparse
import itertools
import json
import math
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import joblib
from fast_multiway_features import family_features

from scripts.develop_task10_detected_ibd import (
    genome_length,
    read_exact_phase,
    read_haplotype_segments,
    read_map,
    save_json,
    sha256,
    verify_inputs,
)


def compare(actual, expected):
    if len(actual) != len(expected):
        raise ValueError('Feature pair count differs')
    maximum = 0.
    for left, right in zip(actual, expected, strict=True):
        if set(left) != set(right):
            raise ValueError('Feature vocabulary differs')
        for key, value in right.items():
            if key in ['indID1', 'indID2']:
                if left[key] != value:
                    raise ValueError('Pair identity differs')
            else:
                if not math.isclose(left[key], value, rel_tol=1e-12, abs_tol=1e-12):
                    raise ValueError(f'Feature differs: {key}: {left[key]} vs {value}')
                maximum = max(maximum, abs(left[key] - value) / (1 + abs(value)))
    return maximum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--people', type=int)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Use Slurm')
    source = args.source_root
    root = source.parent
    protocol = json.loads(args.protocol.read_text())
    if sha256(source / 'source_manifest.json') != protocol['frozen_manifest_sha256']:
        raise ValueError('Frozen source manifest changed')
    verify_inputs(json.loads((source / 'source_manifest.json').read_text()))
    for name, digest in protocol['optimization_source_sha256'].items():
        if sha256(Path(__file__).parent / name) != digest:
            raise ValueError('Optimization source changed')
    map_path = source / 'assets/external/ped-sim-map/refined_average.plink.map'
    genetic_map, genome_cm = read_map(map_path), genome_length(map_path)
    model = joblib.load(source / 'models/phase_model.joblib')
    if args.people is not None:
        prior = root / 'scaling_preflight'
        samples = json.loads((prior / 'roster.json').read_text())[:args.people]
        printed = json.loads((prior / 'printed_roster.json').read_text())
        started = time.monotonic()
        segments = read_exact_phase(prior / 'profile.seg', dict.fromkeys(itertools.combinations(sorted(printed), 2)))
        selected = set(samples)
        segments = {pair: rows for pair, rows in segments.items() if set(pair) <= selected}
        loaded = time.monotonic()
        rows, groups, _ = family_features(samples, segments, genetic_map, genome_cm, model)
        if len(rows) != args.people * (args.people - 1) // 2:
            raise ValueError('Incomplete profile pair roster')
        save_json(args.output / f'profile_{args.people}.json', {'state': 'COMPLETED', 'people': args.people,
            'pairs': len(rows), 'groups': len(groups), 'feature_seconds': time.monotonic() - loaded,
            'total_seconds': time.monotonic() - started,
            'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
        return
    args.output.mkdir(parents=True, exist_ok=False)
    bound, checks = {}, []
    for index in protocol['training_cohort_indices']:
        directory = root / 'detected_panel_v2' / f'cohort_{index:03d}'
        result = json.loads((directory / 'result.json').read_text())
        if result['item']['stage'] != 'training':
            raise ValueError('Equivalence may only read declared training cohorts')
        for name in ['detected.ibd.gz', 'identifiers.json', 'features.json', 'inferred_groups.json']:
            path = directory / name
            bound[str(path)] = result['output_sha256'][str(path)]
        verify_inputs(bound)
        grouped, _ = read_haplotype_segments(directory / 'detected.ibd.gz')
        segments = {}
        for (a, b, chrom), rows in grouped.items():
            segments.setdefault((a, b), {})[chrom] = rows
        samples = sorted(json.loads((directory / 'identifiers.json').read_text()).values())
        started = time.monotonic()
        rows, groups, _ = family_features(samples, segments, genetic_map, genome_cm, model)
        seconds = time.monotonic() - started
        maximum = compare(rows, json.loads((directory / 'features.json').read_text()))
        if list(map(list, groups)) != json.loads((directory / 'inferred_groups.json').read_text()):
            raise ValueError('Inferred sibling groups differ')
        save_json(args.output / f'training_{index:03d}_features.json', rows)
        checks.append({'cohort_index': index, 'pairs': len(rows), 'all_features_equivalent': True,
            'maximum_scaled_error': maximum, 'optimized_seconds': seconds})
    prior = root / 'scaling_preflight'
    bound.update({str(prior / name): sha256(prior / name) for name in ['profile.seg', 'roster.json', 'printed_roster.json']})
    save_json(args.output / 'equivalence.json', {'checks': checks, 'input_sha256': bound,
        'protocol_sha256': sha256(args.protocol), 'tolerance': {'relative': 1e-12, 'absolute': 1e-12},
        'biological_qualification': False})
    profiles = []
    for people, limit in protocol['profiles']:
        command = ['timeout', '--signal=TERM', '--kill-after=10', str(limit), '/usr/bin/time', '-v',
            sys.executable, str(Path(__file__).resolve()), '--source-root', str(source),
            '--protocol', str(args.protocol), '--output', str(args.output), '--people', str(people)]
        with (args.output / f'profile_{people}.log').open('w') as log:
            result = subprocess.run(command, check=False, stdout=log, stderr=subprocess.STDOUT)
        path = args.output / f'profile_{people}.json'
        entry = {'people': people, 'limit_seconds': limit, 'exit_code': result.returncode}
        entry.update(json.loads(path.read_text()) if path.exists() else {'state': 'TIMEOUT' if result.returncode == 124 else 'FAILED'})
        profiles.append(entry)
        save_json(args.output / 'result.json', {'equivalence': checks, 'profiles': profiles,
            'protocol_sha256': sha256(args.protocol), 'biological_qualification': False})
        print(json.dumps(entry), flush=True)
    verify_inputs(bound)


if __name__ == '__main__':
    main()
