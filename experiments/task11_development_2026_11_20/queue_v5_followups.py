"""Queue dependent development fitting/ablation only, never confirmation."""
import json
from pathlib import Path
import subprocess


def submit(dependency, script):
    result = subprocess.run(["sbatch", "--parsable", f"--dependency=afterok:{dependency}", script],
                            capture_output=True, text=True, check=True)
    return int(result.stdout.strip().splitlines()[-1].split(";")[0])


fit = submit(10853295, "fit_v5.sbatch")
ablation = submit(fit, "ablation_v5.sbatch")
plan = {"development_feature_replay": 10853295, "fit": fit, "exposed_cohort_ablation": ablation,
        "fresh_confirmation_launched": False, "submission_count": 0}
with Path("job_plan_v5.json").open("x") as handle:
    json.dump(plan, handle, indent=2)
print(json.dumps(plan))
