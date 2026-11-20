#!/bin/bash
#SBATCH --job-name=ghist-vcf-preflight-11
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=00:05:00
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/final_runs/task11-root-20261120-2002/preflight-%j.log
source /etc/profile
set -euo pipefail
python3 /vast/projects/smathi/cohort/korfmann/ghist2026/final_runs/preflight_final_vcf_20261120.py --input /vast/projects/smathi/cohort/korfmann/ghist2026/final_runs/task11-root-20261120-2002/inputs/relatedness3.final.vcf.gz --output /vast/projects/smathi/cohort/korfmann/ghist2026/final_runs/task11-root-20261120-2002/vcf_preflight.json --expected-sha256 5bbdd1a850f65308ed2376a159cbe63c7540e76eab9e91bb77b75c127a2fc553
