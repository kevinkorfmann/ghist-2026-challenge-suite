"""Report completed stabilization runtime intervals without changing the estimator."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from ghist_lowcoverage.pipeline import sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert args.root.name == 'lowcoverage-20261120' and not args.output.exists()
    report = {
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'reporter_sha256': sha256(Path(__file__)),
        'scope': 'Both separately registered cap500 pragmatic Task14 development cases.',
        'timing_limit': ('Intervals are derived from filesystem modification timestamps, not instrumented subprocess timers. '
                         'The first spectrum interval includes SAF hashing/startup after initial binding; later spectrum '
                         'intervals include command preparation after the previous fit. Fit intervals start after folding.'),
        'cohorts': {},
    }
    for replicate in range(2):
        case = args.root / f'window_pragmatic_development_split_r0_p{replicate}_v1'
        window = case / 'window_comparison'
        previous = case / 'initial_binding.json'
        executed = case / 'execution_binding.json'
        assert executed.is_file(), 'Only completed cases may be reported.'
        rows = {}
        for label in ['seed1', 'seed2', 'continued']:
            log, folded, estimate = [window / f'{label}.{suffix}' for suffix in ['log', 'folded.sfs', 'estimate.json']]
            assert all(path.is_file() for path in [previous, log, folded, estimate])
            text = log.read_text()
            epochs = [int(value) for value in re.findall(r'Finished epoch ([0-9]+)', text)]
            changes = re.findall(r'Current log-likelihood ([^,\n]+), Δ=([^ ]+) ([<>≤≥]=?) ([^\n]+)', text)
            assert epochs and changes
            parsed = json.loads(estimate.read_text())
            rows[label] = {
                'spectrum_interval_seconds': log.stat().st_mtime - previous.stat().st_mtime,
                'fit_interval_seconds': estimate.stat().st_mtime - folded.stat().st_mtime,
                'epochs': epochs[-1], 'last_normalized_likelihood_change': float(changes[-1][1]),
                'printed_stopping_comparison': list(changes[-1]),
                'stopped_before_cap': parsed['winsfs_early_stop_observed'],
                'artifacts': {str(path.relative_to(args.root)): {
                    'sha256': sha256(path), 'mtime_unix_seconds': path.stat().st_mtime}
                    for path in [previous, log, folded, estimate]},
            }
            previous = estimate
        report['cohorts'][f'replicate{replicate}'] = {
            'runs': rows,
            'binding_to_final_binding_seconds': executed.stat().st_mtime - (case / 'initial_binding.json').stat().st_mtime,
            'execution_binding_sha256': sha256(executed),
        }
    write_json(args.output, report)
    print(json.dumps({key: {label: {k: v for k, v in row.items() if k != 'artifacts'}
                           for label, row in value['runs'].items()} for key, value in report['cohorts'].items()}))


if __name__ == '__main__':
    main()
