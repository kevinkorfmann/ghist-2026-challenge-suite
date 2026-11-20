#!/usr/bin/env bash
#SBATCH --job-name=task9-package-v2-preflight
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:05:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/task9-offline-package-v2-20261120/job-%j.log
source /etc/profile >/dev/null 2>&1
set -eo pipefail
module load python/3.11.11
set -u
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026
STAGE="$TASK_ROOT/reviews/task9-offline-package-v2-20261120"
PACKAGE="$STAGE/task9-package-v2"
cd "$STAGE"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export PYTHONPATH="$PACKAGE/src:$PACKAGE"
test ! -e preflight.json
# These independently pinned hashes came from the reviewed local checkout.
printf '%s  %s\n' \
  f1df68f7e5f2397fd5b7d2540d330c601034ea66f54d199e93b7a2c3c95df0de task9-package-v2.tar.gz \
  4c0ec8f9ae12b00b77506b803812c3378f1d43bf873ff206a1db2bffeca5da9e task9-package-v2/package.json \
  | sha256sum --check
"$TASK_ROOT/.venv/bin/python" -B \
  "$PACKAGE/scripts/prepare_task9_submission_package.py" preflight \
  --package "$PACKAGE" --artifact-root "$TASK_ROOT" \
  --sex-map "$TASK_ROOT/external/ped-sim-map/refined_mf.simmap" \
  --expected-package-sha256 4c0ec8f9ae12b00b77506b803812c3378f1d43bf873ff206a1db2bffeca5da9e \
  > preflight.json
sha256sum task9-package-v2.tar.gz task9-package-v2/package.json > transfer_checksums.txt
