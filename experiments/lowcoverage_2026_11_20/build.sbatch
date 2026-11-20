#!/bin/bash
#SBATCH --job-name=ghist-lowcov-build
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=00:15:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/build-%j.log
source /etc/profile
set -euo pipefail
cd /vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/external/htslib
git submodule update --init --recursive
autoreconf -i
./configure --disable-lzma --disable-bz2 --disable-libcurl
make -j4 libhts.a
cd ../angsd
make -j4 angsd HTSSRC=../htslib LIBS="/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/external/htslib/libhts.a -lz -lm -lpthread"
make -C misc -j4 realSFS HTSSRC=../../htslib LIBS="/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/external/htslib/libhts.a -lz -lm -lpthread"
./angsd --version
sha256sum angsd misc/realSFS
