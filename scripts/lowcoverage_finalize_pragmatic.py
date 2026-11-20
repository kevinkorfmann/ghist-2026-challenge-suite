"""Export both completed pragmatic development cases; never access final inputs."""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys


def main():
    root = Path('/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120')
    assert root.is_dir()
    assert all((root / f'window_pragmatic_development_split_r0_p{replicate}_v1/execution_binding.json').is_file()
               for replicate in range(2)), 'Both cases must complete before aggregation.'
    environment = dict(os.environ, TZ='UTC', PYTHONPATH=f'{root}/source_v12/src:{root}/source_v12')
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', '10858309',
        '--format=JobID,State,ExitCode,Elapsed,Start,End,AllocCPUS'], env=environment, text=True)
    columns = ['job_id', 'state', 'exit_code', 'elapsed', 'start', 'end', 'allocated_cpus']
    rows = [dict(zip(columns, line.split('|'))) for line in raw.splitlines() if line.strip()]
    assert {row['job_id'] for row in rows} == {'10858309_0', '10858309_1'}
    assert all(row['state'] == 'COMPLETED' and row['exit_code'] == '0:0' for row in rows)
    cpu_seconds = 0
    for row in rows:
        hours, minutes, seconds = map(int, row['elapsed'].split(':'))
        cpu_seconds += int(row['allocated_cpus']) * (3600 * hours + 60 * minutes + seconds)
    with (root / 'task14_pragmatic_development_accounting_v1.json').open('x') as stream:
        json.dump({'recorded_at_utc': datetime.now(timezone.utc).isoformat(), 'requested_timezone': 'UTC',
            'rows': rows, 'allocated_cpu_hours': cpu_seconds / 3600, 'registered_cpu_hour_cap': 22 / 3,
            'process_completion_does_not_establish_scientific_qualification': True}, stream, indent=2)
    commands = [
        [sys.executable, str(root / 'source_v12/scripts/lowcoverage_assess_pragmatic.py'), '--root', str(root),
            '--freeze-sha', '7a98a0e446cf5b5221fef48bace8743b5a604f7cd807fdff0c26df3df001b1bb',
            '--mode', 'development', '--output', str(root / 'task14_pragmatic_development_gate_audit_v1.json')],
        [sys.executable, str(root / 'lowcoverage_assess_pragmatic_usefulness.py'), '--root', str(root),
            '--output', str(root / 'task14_pragmatic_usefulness_assessment_v1.json')],
        [sys.executable, str(root / 'lowcoverage_pragmatic_timing.py'), '--root', str(root),
            '--output', str(root / 'task14_pragmatic_development_timing_v1.json')],
        [sys.executable, str(root / 'lowcoverage_collect_completed.py'), '--root', str(root),
            '--output', str(root / 'evidence_checkpoint_v6.tar.gz')],
    ]
    for command in commands:
        subprocess.run(command, env=environment, check=True)


if __name__ == '__main__':
    main()
