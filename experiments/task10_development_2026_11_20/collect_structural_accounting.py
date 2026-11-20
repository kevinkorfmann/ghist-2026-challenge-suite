"""Retain final Slurm accounting for every declared qualification cohort."""
import argparse
import json
import re
import subprocess
from pathlib import Path


def kib(value):
    if not value:
        return 0
    match = re.fullmatch(r'([0-9.]+)([KMGT]?)', value)
    if not match:
        raise ValueError(f'Unexpected Slurm RSS unit: {value}')
    number, suffix = match.groups()
    return float(number) * {'': 1 / 1024, 'K': 1, 'M': 1024, 'G': 1024**2, 'T': 1024**3}[suffix]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--stage', choices=['small', 'complete'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Accounting snapshot already exists')
    records = {}
    for index in range(8 if args.stage == 'small' else 10):
        path = args.root / f'cohort_{index:03d}/result.json'
        if not path.exists():
            continue
        result = json.loads(path.read_text())
        job = result['slurm_job_id']
        command = ['sacct', '-j', job, '--noheader', '--parsable2',
            '--format=JobIDRaw,State,ElapsedRaw,MaxRSS']
        text = subprocess.check_output(command, text=True)
        rows = [line.strip().split('|') for line in text.splitlines() if line.strip()]
        main_row = next(row for row in rows if row[0] == job)
        records[str(index)] = {'job_id': job, 'state': main_row[1],
            'elapsed_seconds': int(main_row[2]), 'max_rss_kib': max(kib(row[3]) for row in rows),
            'sacct_command': command, 'raw_sacct': text,
            'memory_interpretation': 'Reported Slurm MaxRSS across job steps, with scheduled memory limit enforced separately'}
    args.output.write_text(json.dumps(records, indent=2) + '\n')


if __name__ == '__main__':
    main()
