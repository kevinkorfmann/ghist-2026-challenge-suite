"""Freeze inference handoff before actual final-data inspection/execution."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import shutil
import sys

r=Path('/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120')
s=r/'detected-panel-source'
q=r/'structural-qualification-source-v3'
h=r/'structural-inference-handoff-v1'
sys.path[:0]=[str(q),str(s/'experiments/task10_development_2026_11_20'),str(s),str(s/'src'),str(s/'scripts'),str(r)]
import structural_inference_handoff  # noqa: E402,F401
from scripts.develop_task10_detected_ibd import read_map,sha256  # noqa: E402
map_path=s/'assets/external/ped-sim-map/refined_average.plink.map'
protocol={'status':'FROZEN_INFERENCE_HANDOFF_NO_FINAL_INPUT_ACCESSED',
 'original_source':str(s),'qualification_source':str(q),'cpus':8,'wall_minutes':90,'memory_GiB':44,
 'model_path':str(r/'detected_panel_v2/model_evaluation/pair_only_model.joblib'),
 'model_sha256':'bdfecf8dc908d1222b43c2b16731f4cbe2f465ed7dc93a8a090ecd4b224d4159',
 'model_calibration_path':str(r/'detected_panel_v2/model_evaluation/model_calibration_freeze.json'),
 'confidence_path':str(r/'structural_calibration_v3/confidence_freeze.json'),
 'confidence_sha256':'3fd057e02e827cdd85d5aca2a5959700f8eed034b4dae451dd9225f969efbc8f',
 'global_correctness_confidence':.7744,
 'runtime':{'python':platform.python_version(),'packages':{name:importlib.metadata.version(name) for name in ['numpy','scipy','scikit-learn','msprime','joblib']}},
 'map_sha256':sha256(map_path),'map_bounds':{chrom:[int(pos[0]),int(pos[-1])] for chrom,(pos,_) in read_map(map_path).items()},
 'minimum_markers_per_chromosome':100,'default_people':229,'required_pairs':26106,
 'crop_rule':'Inclusive fixed training-map intervals; outside-domain records excluded and counted; no cM extrapolation. Training simulations contained no outside-domain markers. Preserve original IDs/alleles. Reject inadequate markers. 100 is operational caller minimum, not validation of marker-density equivalence.',
 'missing_policy_default':'error','pragmatic_missingness_protocol_sha256':'3648c748c8ce1663a1a11167671e1c88890c3a235fad7c1eb68caa48213b9abc',
 'optional_drop_site':'Disabled until frozen 25% thinning diagnostic passes exact declared gates; ≥75% eligible markers on every chromosome. MCAR proxy only; not structured-missingness qualification. No major-allele filling.',
 'confidence_interpretation':'One global calibration correctness estimate; zero individual resolution; not root posterior. Does not assert strict ECE gate has passed.',
 'csv':'Exactly indID1,indID2,relation; all26106 unordered pairs; official14labels; original229IDs; validate_relatedness challenge10 before atomic publication',
 'limitations':['No final labels or challenge truth used','Source domain has mu1.25e-8 versus official humanmu1e-8','Phased complete calls or separately validated pragmatic complete-site filtering','Tree topology and conditional root ambiguity remain; source reports are not individual posterior probabilities','Root owns one final selection/upload/submission; candidate generation itself makes no readiness claim'],
 'import_origins':{name:str(Path(module.__file__).resolve()) for name,module in sorted(sys.modules.items()) if getattr(module,'__file__',None) and (str(Path(module.__file__).resolve()).startswith(str(s)+'/') or str(Path(module.__file__).resolve()).startswith(str(q)+'/'))},
 'contract_tests':4,'tests_scope':'Inclusive endpoints, crop/preparation identity, invalid crop atomicity, all229ID/26106pair CSV completeness, exact deterministic25% thinning with unchanged alleles'}
h.mkdir()
files=['structural_inference_handoff.py','structural_inference_handoff.sbatch','handoff_crop_parity.sbatch','verify_handoff_crop_parity.py','test_handoff_contract.py','freeze_inference_handoff.py']
for name in files:shutil.copy2(r/name,h/name)
(h/'structural_inference_handoff_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
manifest={str(p):sha256(p) for p in h.iterdir()}
for path in [s/'source_manifest.json',q/'source_manifest.json',Path(protocol['model_path']),Path(protocol['confidence_path']),Path(protocol['model_calibration_path']),r/'site-thinning-source-v1/pragmatic_missingness_protocol.json']:
 manifest[str(path)]=sha256(path)
(h/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for path in h.iterdir():path.chmod(0o444)
h.chmod(0o555)
print(json.dumps({'protocol_sha256':sha256(h/'structural_inference_handoff_protocol.json'),'manifest_sha256':sha256(h/'source_manifest.json'),'origins':protocol['import_origins'],'runtime':protocol['runtime']}))
