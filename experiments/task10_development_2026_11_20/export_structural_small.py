"""Export exact small qualification evidence and raw-binding mapping for independent review."""
from pathlib import Path
import hashlib
import json
import shutil
import tarfile

root=Path('/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120')
assessment=root/'structural_small_assessment_v3'
if not (assessment/'result.json').exists():
    raise RuntimeError('Full frozen assessment has not completed')
out=root/'structural_small_export_v3_incomplete'
out.mkdir()
mapping={}
missing=[]
paths=list(assessment.glob('*.json'))
paths += [root/'structural_small_accounting_v3.json',root/'structural_calibration_v3/confidence_freeze.json',root/'structural-qualification-source-v3/source_manifest.json',root/'structural-launch-v3-small/launch_gate.json']
for index in range(8):
    paths += [root/'structural_qualification_v3'/f'cohort_{index:03d}'/name for name in ['result.json','identifiers.json','truth.json','structural_predictions.json','pair_only_predictions.json','baseline_predictions.json','structural_tree.json','missingness_contract/result.json','preregistered_design.json']]
paths.append(root/'structural-small-10858090_0.out')
for original in paths:
    if not original.exists():
        missing.append(str(original))
        continue
    relative=original.relative_to(root)
    target=out/relative
    target.parent.mkdir(parents=True,exist_ok=True)
    digest=hashlib.sha256(original.read_bytes()).hexdigest()
    shutil.copy2(original,target)
    if hashlib.sha256(target.read_bytes()).hexdigest()!=digest or hashlib.sha256(original.read_bytes()).hexdigest()!=digest:
        raise ValueError('Artifact changed during export')
    mapping[str(relative)]={'remote_path':str(original),'sha256':digest}
(out/'missing_files.json').write_text(json.dumps(missing,indent=2)+'\n')
(out/'manifest.json').write_text(json.dumps(mapping,indent=2)+'\n')
with tarfile.open(root/'structural_small_export_v3_incomplete.tar.gz','w:gz') as tar:
    tar.add(out,arcname=out.name)
print(hashlib.sha256((out/'manifest.json').read_bytes()).hexdigest())
