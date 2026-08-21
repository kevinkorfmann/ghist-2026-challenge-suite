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
