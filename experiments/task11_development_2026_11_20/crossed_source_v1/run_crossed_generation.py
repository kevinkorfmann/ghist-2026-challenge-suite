"""All registered crossed development cohorts execute the full synthetic pipeline."""
import datetime
import importlib.util
import os
import sys
import traceback
from pathlib import Path

root, number = Path(sys.argv[1]).resolve(), int(sys.argv[2])
source = root / "source_v5"
sys.path[:0] = [str(source), str(source / "src"), str(source / "scripts")]
from scripts import develop_task11_pedigree as dev

plan_path = root / "crossed_development_plan_v1.json"
plan = dev.load(plan_path)
if plan["schema"] != "task11-crossed-development-v1" or plan["gates"] != dev.GATES:
    raise ValueError("Unknown design or changed gates")
dev.require_hashes(plan["bound_inputs"])
if datetime.datetime.now(datetime.UTC) > datetime.datetime.fromisoformat(plan["budget"]["latest_generation_start_utc"]):
    raise RuntimeError("Hard generation start deadline passed")
if "SLURM_JOB_ID" not in os.environ:
    raise RuntimeError("Compute requires Slurm")
protocol = dev.load(root / "protocol_v5.json")
dev.require_runtime(protocol)
dev.require_hashes(protocol["inputs"])
dev.require_hashes({str(source / path): value for path, value in protocol["source_sha256"].items()})
generator_path = root / "crossed_source_v1/generate_crossed_genotypes.py"
spec = importlib.util.spec_from_file_location("crossed_generator", generator_path)
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)
config = plan["cohorts"][number]
if config["global_index"] != number:
    raise ValueError("Unexpected cohort index")
out = root / "results_crossed_v1" / f"{config['split']}_{config['index']:03d}"
out.mkdir(parents=True, exist_ok=False)
dev.dump(out / "protocol.json", {"master_sha256": dev.digest(plan_path), "cohort": config,
    "base_protocol_sha256": dev.digest(root / "protocol_v5.json"), "source_sha256": plan["bound_inputs"],
    "job_id": os.environ["SLURM_JOB_ID"], "final_challenge_access": False})
result = {"cohort": config, "protocol_sha256": dev.digest(out / "protocol.json"),
          "submission_ready": False, "submission_count": 0}
bounds = {**plan["bound_inputs"], str(plan_path): dev.digest(plan_path)}
try:
    project = root.parents[1]
    map_path = project / "external/dog-recombination/campbell2016_average.plink.map"
    truth = generator.generate_genotypes(map_path, out, config)
    bounds.update({str(out / name): dev.digest(out / name) for name in ("synthetic.vcf.gz", "truth.json")})
    inference = dev.detect_cluster_extract(out, map_path, project / "external/hap-ibd/hap-ibd.jar",
                                          truth["samples"], protocol["clustering_policy"])
    result["inference"] = inference
    result["partition"] = dev.partition_metrics(truth["families"], inference["clusters"]["families"])
    dev.require_hashes(bounds)
    result["completed"] = True
except Exception as error:  # noqa: BLE001 - retain failed-cohort audit and exit nonzero
    result.update(completed=False, failure={"type": type(error).__name__, "message": str(error), "traceback": traceback.format_exc()})
result["output_sha256"] = {str(path.relative_to(out)): dev.digest(path) for path in sorted(out.rglob("*")) if path.is_file()}
dev.dump(out / "result.json", result)
print(result["completed"])
if not result["completed"]:
    raise SystemExit(1)
