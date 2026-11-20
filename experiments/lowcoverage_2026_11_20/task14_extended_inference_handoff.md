# Conditional cap500 Task14 final inference handoff

This package is prepared, not qualified or executed. The root agent owns challenge input selection, matching current references/indices, final execution and submission. Both cap500 development cases and all six fresh confirmation cohorts must pass unchanged gates, and independent audits must finish before use. The wrapper checks the frozen fresh18-target report before input access. Earlier cap200 and original-six failures remain retained.

Use a 4-CPU Slurm job with 5632M per CPU, a strict 55-minute limit, and an actual start no later than 2026-11-20 22:45 UTC. The wrapper enforces that start deadline. The source/runtime/tools/formatter and explicit input hashes are verified before and after execution. It writes a draft YAML only if all final numerical gates pass, and contains no downloader or submission transport.

```sh
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
export PYTHONPATH="$TASK_ROOT/source_v9/src:$TASK_ROOT/source_v9"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python "$TASK_ROOT/lowcoverage_predict_extended_v1.py" \
  --frozen-root "$TASK_ROOT" \
  --directory <ROOT_OWNED_BAM_DIRECTORY> \
  --reference <MATCHING_CURRENT_REFERENCE_FASTA> \
  --output <NEW_ROOT_OWNED_OUTPUT_DIRECTORY>
```

The exact wrapper, source-manifest and freeze hashes, required qualification-report name and resource cap are in `task14_extended_inference_handoff_v1.json`. No actual challenge paths have been selected or read by this agent.
