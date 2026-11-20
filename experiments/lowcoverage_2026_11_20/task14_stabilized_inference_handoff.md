# Stabilized Task 14 inference handoff

This prospective wrapper has not accessed challenge inputs. Root owns any final-input access and submission. The original v7 policy failed all six confirmation stability checks and must not be described as qualified. This handoff applies only if the new frozen stabilization policy passes both development cases, all six fresh 963xxx confirmation cohorts, and independent audit.

Remote wrapper: `lowcoverage_predict_stabilized_v1.py`, SHA256 `a73d7cc7419370eb1303c2de7c5a6f72a207eab651bc601d0c01234894427f98`.

Frozen estimator: `source_v8/source_manifest.json`, SHA256 `441c1a03441a84dfa942700d5a391a08da353109ddcf09e2c3b9bd28615855c0`. Freeze: `task14_stabilization_freeze_v1.json`, SHA256 `35c2cb3576befed1c8bc000ede4baa4c37fcae901a5c9e4c91aa0d8a9f6111f2`.

Before input audit, the wrapper requires `task14_stabilized_confirmation_gate_audit_v1.json` with matching freeze, confirmation mode, exactly 18 primary targets, and every registered gate passed. Its checksum is then bound through execution. Independent root review is additionally required outside that machine check; the report itself remains immutable.

Supply verified absolute paths for `TASK14_INPUT_DIR`, `TASK14_REFERENCE`, and a new `TASK14_OUTPUT`. The directory must contain exactly the twenty sorted/indexed BAMs and `popfile.txt`; the reference must match the current official dataset and have an FAI. Record official input entity versions and hashes separately. Root has reported a current download AccessDenied issue; this wrapper does not bypass or resolve that transport problem.

Use Betty Slurm, partition `genoa-std-mem`, four CPUs, `5632M` per CPU, and a 45-minute initial cap. Coordinate the overall low-coverage job limit before dispatch.

```bash
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
TASK_PYTHON=/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python
export PYTHONPATH="$TASK_ROOT/source_v8/src:$TASK_ROOT/source_v8"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$TASK_ROOT/source_v8"
"$TASK_PYTHON" "$TASK_ROOT/lowcoverage_predict_stabilized_v1.py" \
  --frozen-root "$TASK_ROOT" \
  --directory "$TASK14_INPUT_DIR" \
  --reference "$TASK14_REFERENCE" \
  --output "$TASK14_OUTPUT"
```

The wrapper uses unchanged GL2/SAF filters and the frozen source_v8 full-window estimator: 500 blocks, window 500, tolerance `1e-6`, cap 200 epochs, fixed seeds 761001/761002 plus continuation initialized from seed 1, global folding, and the unchanged eight-start L-BFGS-B/Powell fitter. The primary estimate is seed 1. It does not choose among estimator seeds or read truth values.

Before and after execution it verifies source/runtime/binaries/template, checks actual imported module paths, and hashes raw input BAMs, reference, popfile and indices. It retains raw spectra, logs and all optimizer starts. A complete run writes `execution_binding.json`; only passage of all unchanged numerical gates produces `task14.submission.draft.yaml`. It never uploads or submits. A numerical failure leaves the candidate unqualified and must not trigger a silent tolerance change or seed selection.
