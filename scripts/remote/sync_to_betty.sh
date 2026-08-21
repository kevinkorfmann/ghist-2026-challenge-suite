#!/usr/bin/env bash
set -euo pipefail
REMOTE_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
ssh betty "mkdir -p '$REMOTE_ROOT'/logs '$REMOTE_ROOT'/data/processed '$REMOTE_ROOT'/data/simulations"
rsync -az \
  --exclude '.git/' --exclude '.venv/' --exclude 'artifacts/' --exclude 'checkpoints/' \
  --exclude 'data/' --exclude 'logs/' \
  ./ "betty:$REMOTE_ROOT/"
