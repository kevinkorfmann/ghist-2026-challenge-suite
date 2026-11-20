"""Read-only independent scoring and label/order metamorphic checks."""
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

from scripts.cluster_task11_families import affinity_matrix
from scripts.develop_task11_pedigree import dataset, digest

root = Path(sys.argv[1])
source = root / "source_v3"
protocol = json.loads((root / "protocol_v3.json").read_text())
if protocol["clustering_policy"] != "average":
    raise ValueError("Audit method mismatch")
subprocess.run([sys.executable, str(source / "scripts/score_relatedness_cohorts.py"),
               "--input", str(root / "independent_scoring_input.json"),
               "--out", str(root / "independent_score_v3.json")], check=True)
trials = []
for index in range(8):
    out = root / "results_v3" / f"confirmation_{index:02d}"
    result, truth, *_ = dataset(out / "result.json")
    frame = pd.read_csv(out / "summary/pair_features.csv")
    original = {tuple(sorted(family)) for family in json.loads((out / "clusters.json").read_text())["families"]}
    for permutation in range(20):
        rng = np.random.default_rng(20311100 + index * 100 + permutation)
        new_names = [f"synthetic11_{i}" for i in rng.permutation(len(truth["samples"]))]
        mapping = dict(zip(truth["samples"], new_names))
        reverse = {value: key for key, value in mapping.items()}
        samples = [mapping[s] for s in rng.permutation(truth["samples"])]
        renamed = frame.copy()
        renamed["indID1"] = renamed["indID1"].map(mapping)
        renamed["indID2"] = renamed["indID2"].map(mapping)
        affinity = affinity_matrix(renamed, samples, "relatedness_ge_10")
        distance = affinity.max() - affinity
        np.fill_diagonal(distance, 0)
        labels = AgglomerativeClustering(n_clusters=8, metric="precomputed", linkage="average").fit_predict(distance)
        groups = {tuple(sorted(reverse[s] for s, label in zip(samples, labels) if label == j)) for j in range(8)}
        trials.append({"cohort": index, "permutation": permutation, "partition_unchanged": groups == original})
report = {"scope": "Metamorphic robustness only; the same eight biological cohorts are reused, not 160 independent cohorts.",
          "label_and_order_permutations": len(trials), "unchanged": sum(t["partition_unchanged"] for t in trials),
          "all_passed": all(t["partition_unchanged"] for t in trials), "trials": trials,
          "protocol_sha256": digest(root / "protocol_v3.json"),
          "candidate_sha256": digest(root / "candidate_v3/candidate.json"),
          "independent_score_sha256": digest(root / "independent_score_v3.json"),
          "audit_source_sha256": digest(Path(__file__)),
          "clustering_source_sha256": digest(source / "scripts/develop_task11_pedigree.py"),
          "submission_count": 0}
(root / "confirmation_audit_v3.json").write_text(json.dumps(report, indent=2)+"\n")
print(json.dumps({"permutations": report["label_and_order_permutations"], "unchanged": report["unchanged"]}))
