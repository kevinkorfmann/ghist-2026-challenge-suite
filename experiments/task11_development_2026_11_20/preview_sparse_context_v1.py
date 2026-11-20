"""Preview registered pedigree sizes only; no genotype or detector outcomes."""
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1]).resolve()
source = root / "source_v5"
sys.path[:0] = [str(source), str(source / "src"), str(source / "scripts")]
from scripts import develop_task11_pedigree as dev

design_path = root / "sparse_context_stress_design_v1.json"
if dev.digest(design_path) != "2b4f32a1f8b1083e119ebde8351c16921028f5ce79f4ec4b13c5b4c969897b0c":
    raise ValueError("Stress design changed")
design = dev.load(design_path)
cohorts, total = [], Counter()
for config in design["cohorts"]:
    rng = dev.np.random.default_rng(config["seed"])
    families, counts = [], Counter()
    for index, label in enumerate(dev.LABELS):
        sizes = tuple(int(x) for x in rng.integers(config["min_branch_size"], config["max_branch_size"] + 1, 2))
        third = bool((config["index"] + index) % 2) or (config["complex_topology"] and label.startswith("avuncular"))
        pedigree = dev.focal_pedigree(label, sizes, third, config["complex_topology"])
        label_counts = Counter(relation for _, _, relation in pedigree["pairs"])
        counts.update(label_counts)
        families.append({"focal_label": label, "branch_sizes": sizes, "third_branch": third,
                         "observed_samples": len(pedigree["observed"]), "label_counts": dict(label_counts)})
    total.update(counts)
    cohorts.append({"config": config, "families": families, "label_counts": dict(counts)})
out = root / "sparse_context_topology_preview_v1.json"
if out.exists():
    raise FileExistsError(out)
dev.dump(out, {"scope": "Registered pedigree metadata only; no genotypes/detection/predictions generated",
    "design_sha256": dev.digest(design_path), "preview_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "cohorts": cohorts, "total_label_counts": dict(total), "total_pairs": sum(total.values())})
print(json.dumps({"label_counts": dict(total), "total_pairs": sum(total.values())}))
