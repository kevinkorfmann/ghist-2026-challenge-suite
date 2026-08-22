#!/usr/bin/env bash
#SBATCH --job-name=torchao-probe
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

module load python/3.11.11

PROBE_ROOT="$(mktemp -d /tmp/ghist-torchao-probe.XXXXXX)"
trap 'rm -rf "$PROBE_ROOT"' EXIT

for VERSION in 0.10.0 0.11.0 0.12.0 0.13.0 0.14.0 0.15.0 0.16.0 0.17.0; do
  rm -f "$PROBE_ROOT"/*.whl
  if ! python -m pip download \
    --disable-pip-version-check \
    --no-deps \
    --only-binary=:all: \
    --dest "$PROBE_ROOT" \
    "torchao==$VERSION" >/dev/null; then
    printf '%s unavailable\n' "$VERSION"
    continue
  fi
  WHEEL="$(find "$PROBE_ROOT" -maxdepth 1 -type f -name 'torchao-*.whl' -print -quit)"
  python - "$VERSION" "$WHEEL" <<'PY'
import sys
import zipfile

version, wheel = sys.argv[1:]
with zipfile.ZipFile(wheel) as archive:
    names = set(archive.namelist())
print(version, "nf4tensor=" + str("torchao/dtypes/nf4tensor.py" in names))
PY
done
