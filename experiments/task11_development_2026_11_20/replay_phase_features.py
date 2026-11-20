"""Hash-bound development replay: reuse detector bytes, change phase features only."""
import ast
import copy
import os
import shutil
import sys
from pathlib import Path

from scripts import develop_task11_pedigree as dev

root, index = Path(sys.argv[1]).resolve(), int(sys.argv[2])
split, index = ("train", index) if index < 16 else ("calibration", index - 16)
if index not in range(dev.COUNTS[split]):
    raise ValueError("Only original development cohorts can be replayed")
protocol = dev.load(root / "protocol_v5.json")
old_protocol = dev.load(root / "protocol_v4.json")
dev.validate_protocol(protocol)
dev.require_runtime(protocol)
bounds = {**protocol["inputs"], str(root / "protocol_v5.json"): dev.digest(root / "protocol_v5.json"),
          str(root / "protocol_v4.json"): dev.digest(root / "protocol_v4.json")}
for source, record in ((root / "source_v4", old_protocol), (root / "source_v5", protocol)):
    bounds.update({str(source / path): h for path, h in record["source_sha256"].items()})
dev.require_hashes(bounds)
for name, value in old_protocol["source_sha256"].items():
    if name != "scripts/develop_task11_pedigree.py" and protocol["source_sha256"].get(name) != value:
        raise ValueError(f"A detector/base-feature dependency changed: {name}")


def functions(path):
    return {node.name: ast.dump(node, include_attributes=False) for node in ast.parse(path.read_text()).body
            if isinstance(node, ast.FunctionDef)}


old_functions = functions(root / "source_v4/scripts/develop_task11_pedigree.py")
new_functions = functions(root / "source_v5/scripts/develop_task11_pedigree.py")
preserved = ["detect_cluster_extract", "execute", "fit", "predict", "sibling_pool", "generate_genotypes"]
if any(old_functions[name] != new_functions[name] for name in preserved):
    raise ValueError("Detector/base-feature/generator/classifier code changed")
baseline = root / "results_v4" / f"{split}_{index:02d}"
result, truth, rows, old_features, *_ = dev.dataset(baseline / "result.json")
if dev.load(baseline / "protocol.json")["master_sha256"] != dev.digest(root / "protocol_v4.json"):
    raise ValueError("Baseline master protocol mismatch")
bounds.update({str(baseline / name): h for name, h in result["output_sha256"].items()})
bounds[str(baseline / "result.json")] = dev.digest(baseline / "result.json")
out = root / "results_v5" / f"{split}_{index:02d}"
shutil.copytree(baseline, out)
(out / "result.json").unlink()
phase = dev.haplotype_features(out / "segments.ibd.gz", out / "segments.hbd.gz",
          root.parents[1] / "external/dog-recombination/campbell2016_average.plink.map",
          [(r["indID1"], r["indID2"]) for r in rows])
features = dev.np.concatenate([old_features[:, :-18], phase], axis=1)
if not dev.np.isfinite(features).all():
    raise ValueError("Nonfinite recomputed phase features")
dev.np.savez_compressed(out / "features.npz", features=features)
record = dev.load(out / "protocol.json")
record.update(master_sha256=dev.digest(root / "protocol_v5.json"), source_sha256=protocol["source_sha256"],
              inputs=protocol["inputs"], job_id=os.environ.get("SLURM_JOB_ID"))
dev.dump(out / "protocol.json", record)
result = copy.deepcopy(result)
result["protocol_sha256"] = dev.digest(out / "protocol.json")
result["inference"]["phase_exposure"] = {str(t): {
    "zero_endpoint_denominators": int((phase[:, i*6+2:i*6+4] == 0).sum()),
    "retained_cm_quantiles": dev.np.quantile(phase[:, i*6+2:i*6+4], [0,.25,.5,.75,1]).tolist(),
    "retained_fraction_quantiles": dev.np.quantile(phase[:, i*6+4:i*6+6], [0,.25,.5,.75,1]).tolist()}
    for i, t in enumerate(dev.HBD)}
result["development_feature_replay"] = {"baseline": str(baseline), "result_sha256": dev.digest(baseline / "result.json"),
    "replay_script_sha256": dev.digest(Path(__file__)), "preserved_functions_checked": preserved,
    "detector_executed_this_run": False, "scope": "Original observed detector/base-feature bytes are preserved; only18phase features recomputed. Original stage commands/logs retained as predecessor evidence."}
dev.require_hashes(bounds)
result["output_sha256"] = {str(p.relative_to(out)): dev.digest(p) for p in out.rglob("*") if p.is_file() and p.name != "result.json"}
dev.dump(out / "result.json", result)
print(out)
