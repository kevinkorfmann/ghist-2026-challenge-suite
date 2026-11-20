import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parent
SOURCE = Path('/tmp/task10-evidence-20261120')
DEST = ROOT / 'evidence'
for name in ['diagnostic', 'whole_family_diagnostic', 'development_cohort_0', 'development_cohort_1', 'tree_analysis_0', 'tree_analysis_1', 'tree_analysis_2']:
    shutil.copytree(SOURCE / name, DEST / name, dirs_exist_ok=True)
labels = ['siblings', 'avuncular', 'avuncular_once_removed', 'avuncular_twice_removed',
    'full_1st_cousin', 'full_1st_cousin_once_removed', 'full_1st_cousin_twice_removed',
    'full_1st_cousin_thrice_removed', 'full_2nd_cousin', 'full_2nd_cousin_once_removed',
    'full_2nd_cousin_twice_removed', 'full_3rd_cousin', 'full_3rd_cousin_once_removed', 'full_4th_cousin']
variants = ['raw_haplotype', 'diploid_state', 'exact_diploid_oracle', 'whole_family_bonsai', 'NJ_5', 'NJ_10']
for variant in variants:
    cohorts, artifacts = [], {}
    for index in range(2):
        directory = DEST / f'development_cohort_{index}'
        truth_path = directory / 'truth.json'
        if variant == 'whole_family_bonsai':
            predictions_path = directory / 'consensus/predictions.json'
        elif variant.startswith('NJ_'):
            predictions_path = DEST / f'tree_analysis_{index + 1}/predictions_{variant[3:]}.json'
        else:
            predictions_path = directory / f'representation_comparison/{variant}.json'
        truth = json.loads(truth_path.read_text())
        predictions = json.loads(predictions_path.read_text())
        cohorts.append({'cohort_id': f'task10_development_{index}', 'truth': truth, 'predictions': predictions})
        for path in [truth_path, predictions_path]:
            artifacts[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    target = ROOT / 'scorer_inputs' / f'{variant}.json'
    target.parent.mkdir(exist_ok=True)
    target.write_text(json.dumps({'labels': labels, 'cohorts': cohorts}, indent=2) + '\n')
    artifacts[str(target.relative_to(ROOT))] = hashlib.sha256(target.read_bytes()).hexdigest()
    (target.parent / f'{variant}.manifest.json').write_text(json.dumps({'stage': 'DEVELOPMENT', 'files_sha256': artifacts}, indent=2) + '\n')
    print(variant, 'retained input', target)
