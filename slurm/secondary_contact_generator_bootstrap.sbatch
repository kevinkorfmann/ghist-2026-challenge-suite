#!/usr/bin/env bash
#SBATCH --job-name=sc-gen-bootstrap
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
module load python/3.11.11
PROJECT_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
cd "$PROJECT_ROOT"
OUTPUT_ROOT="${OUTPUT_ROOT:-artifacts/task2_secondary_contact/generator_bootstrap}"
if [[ -n "${DESIGN:-}" ]]; then
  mapfile -t TRUTH_VALUES < <(jq -r ".rows[${SLURM_ARRAY_TASK_ID}].truth[]" "$DESIGN")
  if [[ "${#TRUTH_VALUES[@]}" -ne 6 ]]; then
    echo "could not read six truth values for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
  fi
else
  TRUTH_VALUES=(
    "${ANCESTRAL:-32865}" "${MAINLAND:-20100}" "${ISLAND:-79400}"
    "${SPLIT:-4350}" "${CONTACT:-1390}" "${MIGRATION:-0.000208}"
  )
fi
EXTRA_ARGS=()
if [[ "${FULL_PROJECTION:-0}" == "1" ]]; then
  EXTRA_ARGS+=(--full-projection)
fi
if [[ "${ADAPT_DT:-0}" == "1" ]]; then
  EXTRA_ARGS+=(--adaptive-integration)
fi
if [[ -n "${INITIAL_SCALED:-}" ]]; then
  read -r -a INITIAL_VALUES <<< "$INITIAL_SCALED"
  EXTRA_ARGS+=(--initial-scaled "${INITIAL_VALUES[@]}")
fi
if [[ -n "${RECOMBINATION_MULTIPLIERS:-}" ]]; then
  read -r -a RECOMBINATION_VALUES <<< "$RECOMBINATION_MULTIPLIERS"
  EXTRA_ARGS+=(--recombination-multipliers "${RECOMBINATION_VALUES[@]}")
fi
if [[ -n "${MUTATION_RATE_MULTIPLIERS:-}" ]]; then
  read -r -a MUTATION_VALUES <<< "$MUTATION_RATE_MULTIPLIERS"
  EXTRA_ARGS+=(--mutation-rate-multipliers "${MUTATION_VALUES[@]}")
fi
if [[ -n "${LOCAL_SIZE_MULTIPLIERS:-}" ]]; then
  read -r -a LOCAL_SIZE_VALUES <<< "$LOCAL_SIZE_MULTIPLIERS"
  EXTRA_ARGS+=(--local-size-multipliers "${LOCAL_SIZE_VALUES[@]}")
fi
.venv/bin/python scripts/bootstrap_secondary_contact_generator_fit.py \
  --out "$OUTPUT_ROOT/replicate_${SLURM_ARRAY_TASK_ID}.json" \
  --replicate "$SLURM_ARRAY_TASK_ID" \
  --truth "${TRUTH_VALUES[@]}" \
  --blocks "${BLOCKS:-10}" \
  --block-length "${BLOCK_LENGTH:-10000000}" \
  --projection ${PROJECTION:-40 30} \
  --restarts "${RESTARTS:-12}" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --maxiter "${MAXITER:-400}" \
  --seed "${SEED:-20260919}" \
  --integration-dt-fac "${DT_FAC:-0.02}" \
  --migration-asymmetry "${MIGRATION_ASYMMETRY:-1.0}" \
  --polarization-error "${POLARIZATION_ERROR:-0.0}" \
  "${EXTRA_ARGS[@]}"
