#!/bin/bash
#SBATCH --job-name=ghist-lowcov-winsfs-build
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/build-winsfs-%j.log
source /etc/profile
module load rust/1.85.0
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
export CARGO_HOME="$TASK_ROOT/external/cargo-cache"
cd "$TASK_ROOT/external/winsfs"
cargo build --locked --release -j 4 --package winsfs-cli
./target/release/winsfs --version
sha256sum ./target/release/winsfs
