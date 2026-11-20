#!/usr/bin/env bash
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:05:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/crossed_export_%j.out
#SBATCH --error=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120/logs/crossed_export_%j.err
source /etc/profile
set -euo pipefail
module load python/3.11.11
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/task11-development-20261120
PROJECT_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
"$PROJECT_ROOT/.venv/bin/python" -I "$TASK_ROOT/export_crossed_evidence_v1.py" "$TASK_ROOT"
