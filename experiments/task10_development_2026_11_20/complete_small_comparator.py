"""Complete missing comparator only; preserve the structural inference failure."""
from pathlib import Path
import json
import os
from ghist_relatedness.labels import HUMAN_LABELS
from scripts.develop_task10_detected_ibd import complete_features,detected_segments,frozen_pair_predictions,genome_length,save_json,sha256,write_features

if not os.environ.get('SLURM_JOB_ID'):
    raise RuntimeError('Scheduled diagnostic required')
r=Path('/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120')
s=r/'detected-panel-source';out=r/'structural_small_comparator_v1';out.mkdir()
d=r/'structural_qualification_v3/cohort_000';assets=s/'assets';mp=assets/'external/ped-sim-map/refined_average.plink.map'
roster=sorted(json.loads((d/'identifiers.json').read_text()).values())
f=complete_features(detected_segments(d/'detected.ibd.gz',mp,'raw_haplotype'),roster,genome_length(mp))
write_features(f,out/'baseline_features.tsv')
frozen_pair_predictions(out/'baseline_features.tsv',assets/'artifacts/task10_large_family/clean_spectrum_v1/model.joblib',assets/'artifacts/task10_large_family/pair_model_min5/model.joblib',out/'baseline_predictions.json')
bound={};variants={name:[] for name in ['baseline','pair_only','structural']}
for i in range(8):
    folder=r/'structural_qualification_v3'/f'cohort_{i:03d}'
    truth_path=folder/'truth.json';truth=json.loads(truth_path.read_text());bound[str(truth_path)]=sha256(truth_path)
    roster_path=folder/'identifiers.json';bound[str(roster_path)]=sha256(roster_path)
    for name in variants:
        p=(out/'baseline_predictions.json') if i==0 and name=='baseline' else folder/f'{name}_predictions.json'
        rows=json.loads(p.read_text()) if p.exists() else []
        if p.exists():bound[str(p)]=sha256(p)
        variants[name].append({'cohort_id':f'nested_{i:03d}','truth':truth,'predictions':rows})
for name,cohorts in variants.items():save_json(out/f'{name}_scorer_input.json',{'labels':list(HUMAN_LABELS.values()),'cohorts':cohorts})
save_json(out/'result.json',{'state':'COMPARATOR_COMPLETION_ONLY','strict_failure_preserved':True,'structural_missing_predictions_count':120,'input_sha256':bound,'source_sha256':sha256(Path(__file__))})
