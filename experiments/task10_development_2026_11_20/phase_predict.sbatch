#!/bin/bash
#SBATCH --job-name=task10-phase-predict
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:20:00
#SBATCH --array=0-3
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task10-development-20261120/phase-predict-%A_%a.out
source /etc/profile
set -euo pipefail
module load python/3.11.11
A=/vast/projects/smathi/cohort/korfmann/ghist2026
D="$A/experiments/task10-development-20261120"
S="$D/phase-source"
export PYTHONPATH="$S:$S/src:$S/scripts"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
for variant in phase phase_free; do
"$A/.venv/bin/python" "$S/scripts/develop_task10_detected_ibd.py" phase-predict --asset-root "$A" --model "$D/phase_training/${variant}_model.joblib" --cohort "$D/phase_development_cohort_$SLURM_ARRAY_TASK_ID" --out-dir "$D/phase_predictions/cohort_$SLURM_ARRAY_TASK_ID/$variant"
done
