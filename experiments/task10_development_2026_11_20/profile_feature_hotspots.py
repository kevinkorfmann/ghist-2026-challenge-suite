"""One bounded CPU profile of frozen features on development-only synthetic IBD."""
import cProfile
import itertools
import json
import os
import pstats
import signal
import time
from pathlib import Path

import joblib
from multiway_candidate import family_features, inferred_groups

from scripts.develop_task10_detected_ibd import (
    genome_length,
    phase_pair_features,
    read_exact_phase,
    read_map,
    save_json,
    verify_inputs,
)


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Use Slurm')
    root = Path('/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120')
    source, prior = root / 'detected-panel-source', root / 'scaling_preflight'
    output = root / 'feature_hotspots'
    output.mkdir(exist_ok=False)
    verify_inputs(json.loads((source / 'source_manifest.json').read_text()))
    samples = json.loads((prior / 'roster.json').read_text())[:128]
    printed = json.loads((prior / 'printed_roster.json').read_text())
    segments = read_exact_phase(prior / 'profile.seg', dict.fromkeys(itertools.combinations(sorted(printed), 2)))
    samples_set = set(samples)
    segments = {pair: rows for pair, rows in segments.items() if set(pair) <= samples_set}
    map_path = source / 'assets/external/ped-sim-map/refined_average.plink.map'
    genetic_map, genome_cm = read_map(map_path), genome_length(map_path)
    model = joblib.load(source / 'models/phase_model.joblib')
    started = time.monotonic()
    rows = [{'indID1': a, 'indID2': b, **phase_pair_features(segments.get((a, b), {}), genetic_map, genome_cm)}
        for a, b in itertools.combinations(sorted(samples), 2)]
    groups, _ = inferred_groups(samples, rows, model)
    pair_seconds = time.monotonic() - started
    profiler = cProfile.Profile()

    def stop(_signal, _frame):
        raise TimeoutError('Fixed 60-second profiling interval completed')

    signal.signal(signal.SIGALRM, stop)
    signal.alarm(60)
    state = 'COMPLETED'
    try:
        profiler.runcall(family_features, samples, segments, genetic_map, genome_cm, model)
    except TimeoutError:
        state = 'BOUNDED_PARTIAL_PROFILE'
    finally:
        signal.alarm(0)
        profiler.dump_stats(output / 'features.prof')
        with (output / 'hotspots.txt').open('w') as handle:
            pstats.Stats(profiler, stream=handle).strip_dirs().sort_stats('cumulative').print_stats(25)
    result = {'state': state, 'people': len(samples), 'pairs': len(rows), 'inferred_groups': len(groups),
        'group_sizes': sorted(map(len, groups)), 'pair_features_and_group_seconds': pair_seconds,
        'profile_seconds_cap': 60, 'biological_qualification': False}
    save_json(output / 'result.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
