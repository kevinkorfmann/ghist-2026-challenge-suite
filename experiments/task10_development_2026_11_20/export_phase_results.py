"""Retain raw phase-development artifacts and verifiable scorer manifests."""
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = Path('/tmp/task10-evidence-20261120')
DEST = ROOT / 'phase_evidence'
for name in ['phase_training', 'phase_calibration', 'phase_predictions',
             *[f'phase_development_cohort_{index}' for index in range(4)]]:
    shutil.copytree(SOURCE / name, DEST / name, dirs_exist_ok=True)
for variant in ['phase', 'phase_free', 'baseline']:
    artifacts = {}
    scorer_path = DEST / f'phase_calibration/{variant}_scorer_input.json'
    scorer_input = json.loads(scorer_path.read_text())
    rows_by_cohort = {cohort['cohort_id']: cohort for cohort in scorer_input['cohorts']}
    for index in range(4):
        directory = DEST / f'phase_development_cohort_{index}'
        truth_path = directory / 'truth.json'
        prediction_path = (directory / 'representation_comparison/raw_haplotype.json' if variant == 'baseline' else
            DEST / f'phase_predictions/cohort_{index}/{variant}/predictions.json')
        for kind, path in [('truth', truth_path), ('predictions', prediction_path)]:
            if rows_by_cohort[f'phase_development_{index}'][kind] != json.loads(path.read_text()):
                raise ValueError(f'Scorer input differs from raw {path}')
            artifacts[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    artifacts[str(scorer_path.relative_to(ROOT))] = hashlib.sha256(scorer_path.read_bytes()).hexdigest()
    target = ROOT / 'phase_scorer_inputs' / f'{variant}.manifest.json'
    target.parent.mkdir(exist_ok=True)
    target.write_text(json.dumps({'stage': 'DEVELOPMENT', 'scorer_input': str(scorer_path.relative_to(ROOT)),
        'files_sha256': artifacts}, indent=2) + '\n')
    print(variant, scorer_path)
