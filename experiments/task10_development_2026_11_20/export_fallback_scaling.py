"""Export both fresh cohorts with raw bindings; preserve absent outputs explicitly."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import tarfile

if not os.environ.get('SLURM_JOB_ID'):
    raise RuntimeError('Use the scheduled export wrapper')
root = Path('/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120')
out = root / 'structural_fallback_export_v1'
out.mkdir()
paths = list((root / 'structural_fallback_assessment_v1').glob('*.json'))
paths += [root / 'structural_fallback_accounting_v1.json',
    root / 'structural-launch-fallback-v1/launch_gate.json',
    root / 'structural-launch-fallback-v1/confidence_freeze.json',
    root / 'structural-fallback-source-v1/source_manifest.json',
    root / 'fallback-assessment-source-v1/source_manifest.json']
for index in [8, 9]:
    paths += [root / 'structural_fallback_scaling_v1' / f'cohort_{index:03d}' / name for name in [
        'result.json', 'truth.json', 'identifiers.json', 'structural_predictions.json',
        'pair_only_predictions.json', 'baseline_predictions.json', 'structural_tree.json',
        'missingness_contract/result.json', 'preregistered_design.json']]
    paths.append(root / f'fallback-large-10859039_{index}.out')
manifest, missing = {}, []
for original in paths:
    if not original.exists():
        missing.append(str(original))
        continue
    relative = original.relative_to(root)
    destination = out / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(original.read_bytes()).hexdigest()
    shutil.copy2(original, destination)
    if hashlib.sha256(destination.read_bytes()).hexdigest() != digest or hashlib.sha256(original.read_bytes()).hexdigest() != digest:
        raise ValueError('Evidence changed during export')
    manifest[str(relative)] = {'remote_path': str(original), 'sha256': digest}
(out / 'missing_files.json').write_text(json.dumps(missing, indent=2) + '\n')
(out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
with tarfile.open(root / 'structural_fallback_export_v1.tar.gz', 'w:gz') as archive:
    archive.add(out, arcname=out.name)
print(hashlib.sha256((out / 'manifest.json').read_bytes()).hexdigest())
