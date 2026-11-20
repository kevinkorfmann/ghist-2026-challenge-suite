# Task 14 frozen inference handoff

Prepared for the coordinating root agent. No challenge input has been read by the low-coverage agent. Do not execute on challenge data until all six frozen confirmation cohorts and the independent audit pass. This wrapper does not fetch data or submit files.

The immutable estimator is remote `source_v7`, manifest SHA256 `43a608f4bde3334fc8239f48b4dd4834ee3a87ee547acb15b379bfc4838e21f4`. Confirmation freeze SHA256 is `1f61090830f5d33cd62eb117488d32f906e670d8819add5f30edb1b26b0b11bc`. The explicit-input wrapper `lowcoverage_predict_window_v1.py` has SHA256 `9f3e7d15aee30b18b81399e296825a9cd38ab068c094e29d78f92a3796dff247`. The wrapper is new orchestration around unchanged estimator functions; it has five focused gate tests, and requires independent review before challenge execution.

Run only as a Betty Slurm job with four CPUs, `5632M` per CPU, partition `genoa-std-mem`, and a 45-minute initial cap. Coordinate the two-job low-coverage limit with any remaining development jobs. Supply verified absolute paths for `TASK14_INPUT_DIR`, `TASK14_REFERENCE`, and a new `TASK14_OUTPUT`; the input directory must contain exactly the twenty sorted/indexed BAMs and the official two-column `popfile.txt`. The reference must be the matching current dataset-specific reference with an existing FAI. Do not silently mix October-regenerated inputs with older versions. Record official entity versions and downloaded-file checksums outside the estimator.

```bash
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
TASK_PYTHON=/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python
export PYTHONPATH="$TASK_ROOT/source_v7/src:$TASK_ROOT/source_v7"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$TASK_ROOT/source_v7"
"$TASK_PYTHON" "$TASK_ROOT/lowcoverage_predict_window_v1.py" \
  --frozen-root "$TASK_ROOT" \
  --directory "$TASK14_INPUT_DIR" \
  --reference "$TASK14_REFERENCE" \
  --output "$TASK14_OUTPUT"
```

The command performs both initial seeds and initialized continuation, using exactly the frozen GL2/SAF filters, 500 blocks, window 100, tolerance `1e-4`, maximum 50 epochs, global folding, and eight starts each for L-BFGS-B and Powell. The reported primary estimate is always seed 1. There is no alternate-seed selection or truth access.

The wrapper rechecks frozen source/runtime/binaries/template and raw inputs before and after execution, verifies that imported estimator modules actually come from `source_v7`, and preserves raw spectra, logs and all optimizer records. It always writes a numerical-gate execution record after a complete run; it writes `task14.submission.draft.yaml` only when every numerical gate passes. A failed numerical gate leaves the candidate unqualified. Do not relax a tolerance, choose a different seed, or omit a failed check under this frozen-policy label.

Independent root audit should verify raw numerical evidence, input/source bindings, exact official schema, and confirmation qualification before any submission. The target-specific scientific recovery gate cannot be measured on challenge data; it is established, or fails, on the independent synthetic confirmation set. Submission transport and the final choice remain owned by root.
