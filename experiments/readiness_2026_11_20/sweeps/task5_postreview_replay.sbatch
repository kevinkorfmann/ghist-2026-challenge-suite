#!/bin/bash
#SBATCH --job-name=readiness-sweep-replay
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:30:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/scientific/readiness_sweeps/logs/replay_%j.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
CODE="$ROOT/reviews/readiness-sweeps-postreview-20261120"
OUT="$ROOT/scientific/readiness_sweeps/task5_postreview_replay"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$CODE/src:$CODE"
cd "$CODE"
mkdir "$OUT"
for case in 000 001 002; do
 "$ROOT/.venv/bin/python" -m scripts.predict_single_sweep \
  --vcf "$ROOT/scientific/readiness_sweeps/task5_v4/case_${case}/simulation.vcf" \
  --accessible "$ROOT/data/raw/shared/sweeps/testing.accessible.bed" \
  --dataset testing --bed "$OUT/case_${case}.bed" --out "$OUT/case_${case}.json"
done
"$ROOT/.venv/bin/python" -m scripts.summarize_sweep_readiness \
 --design experiments/readiness_2026_11_20/sweeps/task5_design_v4.json \
 --input-dir "$ROOT/scientific/readiness_sweeps/task5_v4" \
 --out "$OUT/initial_attempt_summary.json"
"$ROOT/.venv/bin/python" - <<'PY'
import hashlib
import json
from pathlib import Path
root = Path("/vast/projects/smathi/cohort/korfmann/ghist2026/scientific/readiness_sweeps")
out = root / "task5_postreview_replay"
checks = []
for index in range(3):
    original_path = root / "task5_v4" / f"case_{index:03d}" / "prediction.json"
    replay_path = out / f"case_{index:03d}.json"
    original, replay = (json.loads(p.read_text()) for p in (original_path, replay_path))
    keys = ("method", "center", "interval", "consensus_eligible", "localization", "contig_length", "haplotypes", "accepted_biallelic_sites")
    if any(original[key] != replay[key] for key in keys):
        raise ValueError("Postreview parser changed a canonical simulation prediction")
    checks.append({"case": index, "identical_prediction_fields": list(keys),
                   "original_sha256": hashlib.sha256(original_path.read_bytes()).hexdigest(),
                   "replay_sha256": hashlib.sha256(replay_path.read_bytes()).hexdigest()})
design = json.loads(Path("experiments/readiness_2026_11_20/sweeps/task5_design_v4.json").read_text())
result = {"status": "PASS", "scope": "Input guard compatibility on the same three completed simulations; not new independent biological replicates", "checks": checks,
          "source_sha256": {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in design["source_sha256"]}}
(out / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
PY
