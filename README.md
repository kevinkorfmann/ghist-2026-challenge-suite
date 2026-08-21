<p align="center">
  <img src="assets/ghist-2026-logo.png" alt="GHIST 2026 logo" width="280">
</p>

# GHIST 2026 challenge suite

This private workspace develops and benchmarks methods for every available GHIST
2026 demographic-history, sweep-detection, and relatedness challenge. It does
**not** submit anything to Synapse. Each challenge has an isolated directory,
registered datasets, reference methods, validation protocol, and frozen holdout.

The first end-to-end target is demographic-history Task 1 (Bottleneck).

Its estimator deliberately combines four signals:

1. a folded site-frequency spectrum (SFS) plus log-distance LD decay,
2. spatial/window summaries of diversity and the SFS, and
3. conditional coalescence hazards from available or inferred genealogies, and
4. an ensemble of convolutional-transformer posterior models trained on coalescent
   simulations with nuisance-parameter randomization.

Public testing answers are isolated under `benchmarks/holdout/`, excluded from
training and selection, and reserved for one-time audits after choices are frozen.

See [the challenge registry](configs/challenges.yaml), [validation protocol](docs/validation_protocol.md),
and [method comparison](docs/method_comparison.md).

## Current status

| Area | State | Evidence / next gate |
|---|---|---|
| Repository | ready | private remote, neutral branch naming, isolated challenge directories |
| Data intake | waiting for sign-in | authenticated downloader, checksums, and forensic VCF audit are ready |
| Scheduled compute | CPU verified | environment bootstrap and end-to-end bottleneck smoke completed through SLURM |
| GPU runtime | verifying | CUDA 12.8 wheel pinned; scheduled B200 gate must pass before full fitting |
| Bottleneck model | prototype verified | SFS, 19-bin LD, window context, genealogical hazards, mixture posterior, tree reference |
| Other challenges | registered | specifications and data entities recorded; implementations proceed after Task 1 freeze |

The recent coalescent time-stratification result from
[Deng, Pritchard, and Spence (2026)](https://doi.org/10.64898/2026.08.11.744210)
is implemented as conditional coalescence-hazard features. These are corrupted and
masked during training to represent inferred-genealogy uncertainty. The posterior
design follows the genome-scale SFS+LD findings of
[Smith et al. (2026)](https://doi.org/10.1093/genetics/iyag107), while retaining
raw spatial context and independent mechanistic references.

## Repository layout

```text
benchmarks/holdout/     public testing answers, isolated from selection
challenges/             one directory for each demographic, sweep, and relatedness task
configs/                dataset/entity registry and frozen task metadata
docs/                   method comparison, validation, input audit, and data access
scripts/                intake, audits, merging, smoke gates, and remote synchronization
slurm/                  scheduled CPU/GPU workflows for betty
src/                    reusable inference implementations
tests/                  deterministic unit and contract checks
```

## Setup

The checked-out workspace already uses `.venv`. To recreate it:

```bash
python3.12 -m venv .venv
.venv/bin/pip install numpy scipy scikit-learn msprime tskit pysam pyyaml pandas matplotlib tqdm torch
```

## Workflow

```bash
# 1. Put the authenticated Synapse download here:
mkdir -p data
# data/bottleneck.testing.vcf.gz

# 2. Inspect and extract the observed feature tensor.
.venv/bin/python -m ghist_bottleneck.cli extract \
  --vcf data/bottleneck.testing.vcf.gz \
  --out artifacts/observed.npz

# 3. Generate simulations matched to the observed sample count.
.venv/bin/python -m ghist_bottleneck.cli simulate \
  --observed artifacts/observed.npz \
  --out artifacts/simulations.npz \
  --n-sims 5000

# 4. Train five independently seeded models and evaluate held-out coverage.
.venv/bin/python -m ghist_bottleneck.cli train \
  --dataset artifacts/simulations.npz \
  --out-dir checkpoints

# 5. Produce a local prediction and YAML (no upload/submission).
.venv/bin/python -m ghist_bottleneck.cli infer \
  --observed artifacts/observed.npz \
  --checkpoint-dir checkpoints \
  --out artifacts/bottleneck.prediction.yaml
```

Use `--quick` on `simulate` and `train` for smoke tests. The full run uses
chromosome-held-out parameter blocks and reports RMSE, median absolute error,
90% interval coverage, and ensemble disagreement.
