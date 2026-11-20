"""Bounded development-only runtime preflight; never a biological qualification."""
from __future__ import annotations

import argparse
import itertools
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import joblib
from multiway_candidate import family_features, whole_family_definition

from scripts.develop_task10_detected_ibd import (
    genome_length,
    read_exact_phase,
    read_map,
    save_json,
    sha256,
    verify_inputs,
)


def worker(args, protocol):
    output = args.out_dir
    samples = json.loads((output / 'roster.json').read_text())[:args.people]
    expected = dict.fromkeys(itertools.combinations(sorted(json.loads((output / 'printed_roster.json').read_text())), 2))
    started = time.monotonic()
    segments = read_exact_phase(output / 'profile.seg', expected)
    selected = set(samples)
    segments = {pair: rows for pair, rows in segments.items() if set(pair) <= selected}
    map_path = args.source_root / 'assets/external/ped-sim-map/refined_average.plink.map'
    model = joblib.load(args.source_root / 'models/phase_model.joblib')
    loaded = time.monotonic()
    rows, groups, _ = family_features(samples, segments, read_map(map_path), genome_length(map_path), model)
    result = {'state': 'COMPLETED', 'people': args.people, 'pairs': len(rows),
        'inferred_groups': len(groups), 'feature_seconds': time.monotonic() - loaded,
        'load_seconds': loaded - started, 'total_seconds': time.monotonic() - started,
        'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'feature_names': sorted(set(rows[0]) - {'indID1', 'indID2'}),
        'biological_qualification': False}
    if len(rows) != args.people * (args.people - 1) // 2:
        raise ValueError('Incomplete feature pair roster')
    save_json(output / f'profile_{args.people}.json', result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--people', type=int)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Runtime profiling requires Slurm')
    protocol = json.loads(args.protocol.read_text())
    if sha256(Path(__file__).resolve()) != protocol['profile_source_sha256']:
        raise ValueError('Profiling source changed after protocol freeze')
    if sha256(args.source_root / 'source_manifest.json') != protocol['source_manifest_sha256']:
        raise ValueError('Frozen candidate manifest changed')
    manifest = json.loads((args.source_root / 'source_manifest.json').read_text())
    verify_inputs(manifest)
    if args.people is not None:
        return worker(args, protocol)
    output = args.out_dir
    output.mkdir(parents=True, exist_ok=False)
    save_json(output / 'predeclared_protocol.json', protocol)
    depths = protocol['depths']
    definition, _, expected = whole_family_definition([depths], 'scale')
    (output / 'profile.def').write_text(definition)
    printed = sorted({sample for pair in expected for sample in pair})
    # Round-robin over terminal branches before adding their second siblings.
    ordered = [f'scale00001_g{depth + 1}-b{branch}-i{child}'
        for child in [1, 2] for branch, depth in enumerate(depths, 1)]
    if set(ordered) != set(printed):
        raise ValueError('Declared roster does not match generated pedigree')
    save_json(output / 'printed_roster.json', printed)
    save_json(output / 'roster.json', ordered[:229])
    assets = args.source_root / 'assets'
    command = [assets / 'external/ped-sim/ped-sim', '-d', output / 'profile.def',
        '-m', assets / 'external/ped-sim-map/refined_mf.simmap', '-o', output / 'profile',
        '--intf', assets / 'external/ped-sim/interfere/nu_p_campbell.tsv',
        '--seed', str(protocol['meiosis_seed']), '--seg_haps', '--fam']
    with (output / 'simulation.log').open('w') as log:
        subprocess.run(list(map(str, command)), check=True, stdout=log, stderr=subprocess.STDOUT, timeout=120)
    results = []
    for people, limit in protocol['profiles']:
        words = ['/usr/bin/time', '-v', sys.executable, str(Path(__file__).resolve()),
            '--source-root', str(args.source_root), '--protocol', str(args.protocol),
            '--out-dir', str(output), '--people', str(people)]
        started = time.monotonic()
        # External timeout kills the complete process group and retains GNU time's RSS report.
        with (output / f'profile_{people}.log').open('w') as log:
            result = subprocess.run(['timeout', '--signal=TERM', '--kill-after=10', str(limit), *words],
                stdout=log, stderr=subprocess.STDOUT, check=False)
        entry = {'people': people, 'limit_seconds': limit, 'exit_code': result.returncode,
            'wall_seconds': time.monotonic() - started}
        path = output / f'profile_{people}.json'
        entry.update(json.loads(path.read_text()) if path.exists() else {'state': 'TIMEOUT' if result.returncode == 124 else 'FAILED'})
        results.append(entry)
        save_json(output / 'result.json', {'profiles': results, 'biological_qualification': False,
            'source_manifest_sha256': protocol['source_manifest_sha256'], 'protocol_sha256': sha256(args.protocol)})
        print(json.dumps(entry), flush=True)
    verify_inputs(manifest)


if __name__ == '__main__':
    main()
