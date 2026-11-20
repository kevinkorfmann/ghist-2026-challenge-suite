#!/bin/bash
#SBATCH --job-name=ghist-vcf-preflight-10
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=00:05:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/final_runs/task10-root-20261120-2028/preflight-%j.log
source /etc/profile
set -euo pipefail
python3 /vast/projects/smathi/cohort/korfmann/ghist2026/final_runs/preflight_final_vcf_20261120.py --input /vast/projects/smathi/cohort/korfmann/ghist2026/final_runs/task10-root-20261120-2028/inputs/relatedness2_final.vcf.gz --output /vast/projects/smathi/cohort/korfmann/ghist2026/final_runs/task10-root-20261120-2028/vcf_preflight.json --expected-sha256 9e1a20af4907befce07ead7b80aa61331c0b9a3ea2c4c451fc4acc9f6a131bc4
