"""Parity and full-denominator development evidence for explicit fallback policy."""
from pathlib import Path
import json
import os
import joblib
from structural_fallback_policy import infer_with_fallback
from structural_pair_candidate import pair_predictions
from structural_confidence import attach_global_confidence
from scripts.develop_task10_detected_ibd import save_json,sha256,verify_inputs

if not os.environ.get('SLURM_JOB_ID'):
    raise RuntimeError('Scheduled parity required')
r=Path('/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120');source=Path(__file__).resolve().parent
verify_inputs(json.loads((source/'source_manifest.json').read_text()))
e=r/'detected_panel_v2/model_evaluation';bundle=joblib.load(e/'pair_only_model.joblib');freeze=json.loads((e/'model_calibration_freeze.json').read_text());confidence=json.loads((r/'structural_calibration_v3/confidence_freeze.json').read_text())
out=r/'structural_fallback_parity_v1';out.mkdir();bound={};calibration=[];development=[];small=[]
for i in range(32,48):
 f=r/'detected_panel_v2'/f'cohort_{i:03d}'/'features.json';bound[str(f)]=sha256(f)
 pair=pair_predictions(json.loads(f.read_text()),bundle,freeze['temperatures']['pair_only']['temperature'])
 rows,tree=infer_with_fallback(pair)
 if tree['pragmatic_policy']['used_fallback']:
  raise ValueError('Previously valid calibration/development changed inferability')
 if i<40:
  old=r/'structural_calibration_v3'/f'cohort_{i:03d}_inference.json';expected=json.loads(old.read_text())['predictions'];bound[str(old)]=sha256(old)
  if rows!=expected:raise ValueError('Calibration labels changed')
  calibration.append({'index':i,'pairs':len(rows),'unchanged':True})
 else:
  original=json.loads((r/'integer_edge_decoder_v1/pair_only_scorer_input.json').read_text())['cohorts'][i-40]['predictions']
  if [(x['indID1'],x['indID2'],x['relation']) for x in rows]!=[(x['indID1'],x['indID2'],x['relation']) for x in original]:raise ValueError('Development labels changed')
  development.append({'index':i,'pairs':len(rows),'unchanged':True})
for i in range(8):
 d=r/'structural_qualification_v3'/f'cohort_{i:03d}'
 for name in ['pair_only_predictions.json','truth.json','identifiers.json']:
  bound[str(d/name)]=sha256(d/name)
 pair=json.loads((d/'pair_only_predictions.json').read_text());rows,tree=infer_with_fallback(pair)
 rows=attach_global_confidence(rows,list(map(str,bundle['model'].classes_)),confidence)
 save_json(out/f'cohort_{i:03d}_inference.json',{'predictions':rows,'tree':tree})
 small.append({'cohort_id':f'nested_{i:03d}','truth':json.loads((d/'truth.json').read_text()),'predictions':rows})
verify_inputs(bound)
save_json(out/'scorer_input.json',{'labels':list(map(str,bundle['model'].classes_)),'cohorts':small})
save_json(out/'result.json',{'state':'PARITY_PASS_NEW_POLICY_SMALL_DEVELOPMENT_ONLY','calibration':calibration,'development':development,'small_pairs':sum(len(c['truth']) for c in small),'strict_v3_status':'FAILED','input_sha256':bound,'source_manifest_sha256':sha256(source/'source_manifest.json'),'output_sha256':{str(p):sha256(p) for p in out.iterdir() if p.is_file()}})
