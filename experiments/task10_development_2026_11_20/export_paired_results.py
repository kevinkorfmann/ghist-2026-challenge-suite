"""Retain paired diagnostic rows and SHA mappings for independent rescoring."""
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = Path('/tmp/task10-evidence-20261120/paired_phase_diagnostic')
DEST = ROOT / 'paired_phase_evidence'
shutil.copytree(SOURCE, DEST, dirs_exist_ok=True)
labels = json.loads((ROOT / 'phase_evidence/phase_training/result.json').read_text())['labels']
for prediction in sorted((DEST / 'cohort_0').glob('*_predictions.json')):
    variant = prediction.name.removesuffix('_predictions.json')
    artifacts, cohorts = {}, []
    for index in range(4):
        truth_path = ROOT / f'phase_evidence/phase_development_cohort_{index}/truth.json'
        prediction_path = DEST / f'cohort_{index}/{prediction.name}'
        cohorts.append({'cohort_id': f'phase_development_{index}',
            'truth': json.loads(truth_path.read_text()), 'predictions': json.loads(prediction_path.read_text())})
        for path in [truth_path, prediction_path]:
            artifacts[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    scorer_path = ROOT / f'paired_phase_scorer_inputs/{variant}.json'
    scorer_path.parent.mkdir(exist_ok=True)
    scorer_path.write_text(json.dumps({'labels': labels, 'cohorts': cohorts}, indent=2) + '\n')
    artifacts[str(scorer_path.relative_to(ROOT))] = hashlib.sha256(scorer_path.read_bytes()).hexdigest()
    scorer_path.with_suffix('.manifest.json').write_text(json.dumps({'stage': 'PAIRED_DEVELOPMENT_DIAGNOSTIC',
        'scorer_input': str(scorer_path.relative_to(ROOT)), 'files_sha256': artifacts}, indent=2) + '\n')
print('Exported paired diagnostic scorer inputs')
