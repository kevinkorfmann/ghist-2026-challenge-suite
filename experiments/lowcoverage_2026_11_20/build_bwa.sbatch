#!/bin/bash
#SBATCH --job-name=ghist-lowcov-bwa
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:10:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/build-bwa-%j.log
source /etc/profile
set -euo pipefail
cd /vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/external/bwa
make -j4
sha256sum bwa
git rev-parse HEAD
