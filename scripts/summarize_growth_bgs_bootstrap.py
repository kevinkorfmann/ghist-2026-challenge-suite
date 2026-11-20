from __future__ import annotations

import argparse
import json
from pathlib import Path

from ghist_inference.bootstrap import summarize_payloads
from ghist_inference.checkpoint import atomic_json, file_digest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--job', type=int, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('Preserve existing bootstrap summary; use a fresh output')
    files = sorted(args.input_dir.glob('task_*.json'))
    identities = {str(path): file_digest(path) for path in files}
    payloads = [json.loads(path.read_text()) for path in files]
    summary = summarize_payloads(payloads, 3)
    summary.update(job_id=args.job, task_files=len(files), source_artifacts=identities,
                   projection=summary['bootstrap_contract']['settings']['projection'])
    for path in files:
        if file_digest(path) != identities[str(path)]:
            raise ValueError('Bootstrap task changed during aggregation')
    atomic_json(args.out, summary)
    if summary['state'] != 'COMPLETED':
        raise RuntimeError('Incomplete or failed bootstrap plan retained; interval output suppressed')
    print(json.dumps({'state': summary['state'], 'replicates': summary['replicates']}))


if __name__ == '__main__':
    main()
