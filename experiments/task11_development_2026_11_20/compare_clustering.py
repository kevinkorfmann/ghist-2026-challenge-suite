"""Training-only comparison; never inspects calibration/confirmation labels."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from scripts.cluster_task11_families import affinity_matrix
from scripts.validate_task11_complete_pipeline import partition_metrics
root = Path(sys.argv[1])
rows = []
for i in range(16):
    out = root / 'results_v2' / f'train_{i:02d}'
    truth = json.loads((out / 'truth.json').read_text())
    stored = json.loads((out / 'result.json').read_text())
    if not stored['completed']:
        raise ValueError(f'Incomplete {out}')
    frame = pd.read_csv(out / 'summary/pair_features.csv')
    affinity = affinity_matrix(frame, truth['samples'], 'relatedness_ge_10')
    distance = affinity.max() - affinity
    np.fill_diagonal(distance, 0)
    options = {'spectral': stored['inference']['clusters']['families']}
    for linkage in ('single', 'average'):
        labels = AgglomerativeClustering(n_clusters=8, metric='precomputed', linkage=linkage).fit_predict(distance)
        options[linkage] = [[s for s, label in zip(truth['samples'], labels) if label == j] for j in range(8)]
    for method, families in options.items():
        rows.append({'cohort': i, 'method': method, **partition_metrics(truth['families'], families)})
summary = {}
for method in ('spectral', 'single', 'average'):
    selected = [r for r in rows if r['method'] == method]
    summary[method] = {name: sum(r[name] for r in selected) for name in (
        'exactly_recovered_families', 'missed_within_pairs', 'false_cross_family_pairs', 'true_within_pairs')}
result = {'scope': '16 training cohorts only; no calibration/confirmation labels accessed', 'summary': summary, 'cohorts': rows}
(root / 'clustering_training_comparison.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(summary, indent=2))
