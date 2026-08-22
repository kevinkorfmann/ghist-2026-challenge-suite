#!/usr/bin/env bash
set -euo pipefail
REMOTE_ROOT="/vast/projects/smathi/cohort/korfmann/ghist2026"
if [[ -n "$(git status --porcelain)" ]]; then
  echo "Refusing to sync an uncommitted worktree" >&2
  exit 1
fi
ssh betty "mkdir -p '$REMOTE_ROOT'/logs '$REMOTE_ROOT'/data/processed '$REMOTE_ROOT'/data/simulations"
rsync -az \
  --exclude '.git/' --exclude '.venv/' --exclude 'artifacts/' --exclude 'checkpoints/' \
  --exclude 'data/' --exclude 'logs/' \
  ./ "betty:$REMOTE_ROOT/"
REVISION="$(git rev-parse HEAD)"
ssh betty "printf '%s\\n' '$REVISION' > '$REMOTE_ROOT/.source-revision'"
