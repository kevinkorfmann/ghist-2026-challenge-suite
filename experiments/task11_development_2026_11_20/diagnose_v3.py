"""Post-confirmation descriptive diagnostics; never a qualification score."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import sys

root = Path(sys.argv[1])
groups = defaultdict(lambda: {"truth": Counter(), "correct": Counter(), "confusion": defaultdict(Counter)})
rows = []
for index in range(8):
    cohort = root / "results_v3" / f"confirmation_{index:02d}"
    result = json.loads((cohort / "result.json").read_text())
    truth = json.loads((cohort / "truth.json").read_text())
    predictions = {tuple(sorted((r["indID1"], r["indID2"]))): r for r in json.loads((cohort / "predictions.json").read_text())}
    with (cohort / "summary/hbd_features.csv").open() as handle:
        hbd = {r["individual"]: float(r["inbreeding_fraction"]) for r in csv.DictReader(handle)}
    for row in truth["pairs"]:
        key = tuple(sorted((row["indID1"], row["indID2"])))
        prediction = predictions.get(key, {}).get("relation", "MISSING")
        endpoint_hbd = max(hbd.get(s, 0) for s in key)
        bin_name = "<0.05" if endpoint_hbd < .05 else "0.05-0.15" if endpoint_hbd < .15 else ">=0.15"
        topology = "complex" if result["cohort"]["complex_topology"] else "simple"
        kind = truth["pedigrees"][row["family"]]["kind"]
        for name in ("all", "topology:" + topology, "ancestral_mating:" + kind, "max_detected_hbd:" + bin_name):
            group = groups[name]
            group["truth"][row["relation"]] += 1
            group["correct"][row["relation"]] += prediction == row["relation"]
            group["confusion"][row["relation"]][prediction] += 1
        rows.append({"cohort": index, **row, "prediction": prediction, "topology": topology,
                     "ancestral_mating": kind, "max_detected_hbd": endpoint_hbd})
report = {"scope": "Post-failure descriptive development diagnostics; overlapping strata are not independent replicates",
          "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          "groups": {name: {"counts": group["truth"], "correct": group["correct"],
                           "recall": {label: group["correct"][label] / count for label, count in group["truth"].items()},
                           "confusion": group["confusion"]} for name, group in groups.items()}, "rows": rows}
(root / "diagnostics_v3.json").write_text(json.dumps(report, indent=2) + "\n")
label = "avuncular_once_removed"
for name, group in report["groups"].items():
    print(name, group["correct"].get(label, 0), group["counts"].get(label, 0), dict(group["confusion"].get(label, {})))
