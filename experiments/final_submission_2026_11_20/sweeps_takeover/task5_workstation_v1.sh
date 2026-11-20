#!/bin/bash
set -uo pipefail
SWEEP_ROOT=/home/kkor/ghist-final-sweeps-20261120-2120
cd "$SWEEP_ROOT"
ulimit -v 12582912
export PYTHONPATH="$SWEEP_ROOT/src:$SWEEP_ROOT/scripts"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export PYTHONDONTWRITEBYTECODE=1
date -u +%FT%TZ > logs/task5_started_utc.txt
/usr/bin/time -v timeout --signal=TERM --kill-after=30s 3600s /home/kkor/ghist-2026-review-20261120/.venv/bin/python "$SWEEP_ROOT/experiments/submission_package_2026_11_20/sweeps_best_effort/run_final_sweep.py" --source-root "$SWEEP_ROOT" --source-freeze "$SWEEP_ROOT/experiments/submission_package_2026_11_20/sweeps_best_effort/source_freeze.json" --raisd "$SWEEP_ROOT/bin/RAiSD" --task 5 --vcf "$SWEEP_ROOT/inputs/onesweep.final.vcf.gz" --vcf-sha256 1a11d96cc0a2742e6bb5428c2f72a21278bdb9b3eafc0e5612788eda9100eff0 --accessible "$SWEEP_ROOT/inputs/9.accessible.bed" --mask-sha256 9f5235a687683eec6a7f19b0cf1334f545bda23497ff8e715658b1086778ec33 --out-dir "$SWEEP_ROOT/task5_v1"
TASK5_EXIT=$?
echo "$TASK5_EXIT" > logs/task5_exit_code.txt
date -u +%FT%TZ > logs/task5_ended_utc.txt
exit "$TASK5_EXIT"
