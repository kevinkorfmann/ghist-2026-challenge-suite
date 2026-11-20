"""Read-only final sweep integrity audit; no accuracy or qualification claim."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


p = argparse.ArgumentParser()
p.add_argument('--run-dir', type=Path, required=True)
p.add_argument('--source-root', type=Path, required=True)
p.add_argument('--task', type=int, required=True)
p.add_argument('--input-sha256', required=True)
a = p.parse_args()
c = json.loads((a.run_dir / 'completion.json').read_text())
assert c['state'] == 'COMPLETED' and c['task'] == a.task
assert c['scientifically_validated'] is False and c['external_submissions'] == 0
assert c['inputs']['vcf']['sha256'] == a.input_sha256
assert c['contig'].removeprefix('chr') == '9' and c['contig_length'] == 138394717
assert c['samples_diploid'] == 25
checked = {}
for name, digest in c['artifact_sha256'].items():
    rel = Path(name)
    assert not rel.is_absolute() and '..' not in rel.parts
    assert sha(a.run_dir / rel) == digest, name
    checked[name] = digest
for name, info in c['inputs'].items():
    assert sha(info['path']) == info['sha256'], name
for name, digest in c['source_sha256'].items():
    assert sha(a.source_root / name) == digest, name
package = a.source_root / 'experiments/submission_package_2026_11_20/sweeps_best_effort'
assert sha(package / 'source_freeze.json') == c['source_freeze_sha256']
assert sha(package / 'run_final_sweep.py') == c['runner_sha256']
assert sha(package / 'method_protocol.json') == c['protocol_sha256']
freeze = json.loads((package / 'source_freeze.json').read_text())
assert freeze['source_sha256'] == c['source_sha256']
mask = []
for line in Path(c['inputs']['mask']['path']).read_text().splitlines():
    if line.strip() and not line.startswith('#'):
        chrom, left, right = line.split()[:3]
        if chrom.removeprefix('chr') == '9':
            mask.append((int(left), int(right)))
intervals = []
for line in (a.run_dir / 'deadline_submission.bed').read_text().splitlines():
    fields = line.split()
    assert len(fields) == 3 and fields[0] == '9'
    left, right = map(int, fields[1:])
    assert 0 <= left < right <= c['contig_length']
    assert any(lo <= left < right <= hi for lo, hi in mask)
    intervals.append([left, right])
assert intervals and (a.task == 8 or len(intervals) == 1)
assert intervals == sorted(intervals) == c['final_intervals']
assert all(x[1] <= y[0] for x, y in zip(intervals, intervals[1:]))
assert c['format_validation']['format_valid'] is True
assert c['format_validation']['challenge'] == a.task
assert c['format_validation']['dataset'] == 'final'
raw = json.loads((a.run_dir / 'raw/prediction.json').read_text())
centers = [raw['center']] if a.task in (5, 7) else [raw['position']] if a.task == 6 else [r['center'] for r in raw['candidates']]
assert centers == [r['raw_center'] for r in c['geometry_adjustments']]
half_width = {5: 60000, 6: 75000, 7: 80000, 8: 60000}[a.task]
for item in c['geometry_adjustments']:
    center = item['raw_center']
    adjusted = min((min(max(center, lo), hi - 1) for lo, hi in mask), key=lambda x: (abs(x-center), x))
    lo, hi = next((lo, hi) for lo, hi in mask if lo <= adjusted < hi)
    expected = [max(0, lo, adjusted-half_width), min(c['contig_length'], hi, adjusted+half_width)]
    assert item['final_center'] == adjusted and item['interval'] == expected
if a.task == 8:
    assert raw['task7_completion_sha256'] is not None
    assert raw['calibrated_null'] is False
print(json.dumps({
    'task': a.task, 'checked_utc': datetime.now(timezone.utc).isoformat(),
    'status': 'PASS', 'run_directory': str(a.run_dir),
    'candidate_sha256': checked['deadline_submission.bed'],
    'completion_sha256': sha(a.run_dir / 'completion.json'),
    'artifact_sha256': checked, 'final_intervals': intervals,
    'source_and_input_hashes_pass': True, 'geometry_reconstructed': True,
    'scientific_qualification': False, 'submission_performed': False,
}, indent=2))
