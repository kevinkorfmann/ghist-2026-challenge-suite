from __future__ import annotations

import argparse
import json
from pathlib import Path

from ghist_inference.bootstrap import summarize_payloads
from ghist_inference.checkpoint import atomic_json, file_digest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--task-dir', type=Path, required=True)
    parser.add_argument('--block-summary', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('Preserve existing bootstrap summary; use a fresh output')
    files = sorted(args.task_dir.glob('task_*.json'))
    identities = {str(path): file_digest(path) for path in [*files, args.block_summary]}
    payloads = [json.loads(path.read_text()) for path in files]
    summary = summarize_payloads(payloads, 4)
    summary.update(source_artifacts=identities,
                   supplied_block_design=json.loads(args.block_summary.read_text()),
                   supplied_block_design_scope='Descriptive auxiliary record; resampling input is bound by the block archive hash in the producer contract.')
    for path in [*files, args.block_summary]:
        if file_digest(path) != identities[str(path)]:
            raise ValueError('Bootstrap source changed during aggregation')
    atomic_json(args.out, summary)
    if summary['state'] != 'COMPLETED':
        raise RuntimeError('Incomplete or failed bootstrap plan retained; interval output suppressed')
    print(json.dumps({'state': summary['state'], 'replicates': summary['replicates']}))


if __name__ == '__main__':
    main()
