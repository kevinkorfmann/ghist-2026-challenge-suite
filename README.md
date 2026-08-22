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

Nominal-rate and confounder-perturbed simulations are labeled and evaluated
separately. Perturbations cover mutation/recombination uncertainty, local
recombination-rate variation, missing calls, and smooth diversity suppression;
the nominal domain retains half of the training mass for precision.

Public testing answers are isolated under `benchmarks/holdout/`, excluded from
training and selection, and reserved for one-time audits after choices are frozen.

See [the challenge registry](configs/challenges.yaml), [validation protocol](docs/validation_protocol.md),
and [method comparison](docs/method_comparison.md).

## Current status

| Area | State | Evidence / next gate |
|---|---|---|
| Repository | ready | private remote, neutral branch naming, isolated challenge directories |
| Data intake | testing data profiled | Synapse checksums match; 1,194,350 fully phased SNPs and 114 haplotypes audited in job `7731234` |
| Scheduled compute | CPU and GPU verified | CPU end-to-end job `7731188`; B200 forward/backward gate `7731200` |
| GPU runtime | ready | PyTorch 2.11.0, CUDA 12.8, NVIDIA B200 MIG 2g.45gb |
| CXT reference | ready | pinned 19,484,400-parameter checkpoint loaded in job `7731231`; SHA-256 recorded |
| FastCXT reference | translated | grouped-cell training job `7731776`; validation-selected checkpoint translated all 100 Mb in job `7731790` |
| Bottleneck model | Task 1 frozen, not submitted | primary estimate 4,587 generations and 0.52242 post-decline fraction; exact chromosome recovery and three mechanistic engines agree |
| Secondary Contact | active | exact hidden-epoch SFS fit stable across polarization and projection; matched score-aligned SFS and full 100 Mb SINGER/POLEGON campaigns running |
| Remaining challenges | data secured | all registered inputs, templates, maps, accessibility masks, and sealed finals checksum-verified |

The recent coalescent time-stratification result from
[Deng, Pritchard, and Spence (2026)](https://doi.org/10.64898/2026.08.11.744210)
is implemented as conditional coalescence-hazard features. These are corrupted and
masked during training to represent inferred-genealogy uncertainty. The posterior
design follows the genome-scale SFS+LD findings of
[Smith et al. (2026)](https://doi.org/10.1093/genetics/iyag107), while retaining
raw spatial context and independent mechanistic references.

Two independent translation branches are supported. The autoregressive
[CXT](https://github.com/kr-colab/cxt) bridge converts posterior pairwise TMRCA
samples into conditional hazards. The bidirectional Mamba
[FastCXT](https://github.com/kevinkorfmann/fastcxt) bridge analytically integrates
its Gaussian log-TMRCA means and variances over each time interval before the
same hazard transform. Both retain blockwise uncertainty and can be evaluated
alone, masked, or ensembled without treating one inferred genealogy as truth.

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
.venv/bin/pip install -e '.[references]'
```

## Workflow

```bash
# 1. Put the authenticated Synapse download here:
mkdir -p data
# data/bottleneck.testing.vcf.gz

# 2. Inspect and extract the observed feature tensor.
.venv/bin/python -m ghist_bottleneck.cli sfs-reference \
  --vcf data/bottleneck.testing.vcf.gz \
  --out artifacts/bottleneck.sfs-reference.json

.venv/bin/python -m ghist_bottleneck.cli extract-blocks \
  --vcf data/bottleneck.testing.vcf.gz \
  --out artifacts/observed.blocks.npz

# 3. Generate simulations matched to the observed sample count.
.venv/bin/python -m ghist_bottleneck.cli simulate \
  --observed artifacts/observed.blocks.npz \
  --out artifacts/simulations.npz \
  --n-sims 5000

# 4. Train five independently seeded models and evaluate held-out coverage.
.venv/bin/python -m ghist_bottleneck.cli train \
  --dataset artifacts/simulations.npz \
  --out-dir checkpoints

# 5. Produce a local prediction and YAML (no upload/submission).
.venv/bin/python -m ghist_bottleneck.cli infer-blocks \
  --observed artifacts/observed.blocks.npz \
  --checkpoint-dir checkpoints \
  --out artifacts/bottleneck.prediction.yaml
```

Use `--quick` on `simulate` and `train` for smoke tests. The full run uses
chromosome-held-out parameter blocks and reports RMSE, median absolute error,
90% interval coverage, and ensemble disagreement.

The context model consumes all fifty contiguous 2 Mb blocks, matching the
physical scale of full simulations. A robust precision aggregation limits the
leverage of localized recombination, accessibility, or diversity anomalies,
while the whole 100 Mb VCF independently feeds the mechanistic SFS fit.

## Task 1 frozen result

The frozen local estimate is 4,587.06 generations and a post-decline fraction of
0.5224206. A 200-replicate bootstrap over the fifty physical 2 Mb blocks gives
95% intervals of 4,183.89–4,982.07 generations and 0.51306–0.53080. Nothing has
been uploaded or submitted.

A separate score-aligned contingency was trained on 900 full 100 Mb simulations:
100 space-filling points across the complete prior and 50 denser points near the
observed fit, with six genomic replicates per point and parameter-point-isolated
calibration, selection, and untouched evaluation. Its selected two-layer
calibrator predicts 4,338.37 and 0.49587. It cut untouched broad-prior RRMSE from
5.789 to 0.898 and focused RRMSE from 0.506 to 0.206, but worsened typical focused
natural errors and point-mean biases. It is therefore retained for extreme
held-out estimates rather than replacing the mechanistic primary for this VCF.

On 29 untouched nominal whole-chromosome cells, the raw moments estimator had a
13.28% median absolute relative generation error and a 0.01132 median absolute
fraction error. dadi recovered 13.70% and 0.01173; the learned context ensemble
recovered 23.44% and 0.03014 and remains secondary. A selection-chosen blend was
rejected because it worsened both natural-scale medians on untouched evaluation.

Observed cross-checks are tightly grouped: dadi folded gives 4,643.18 and
0.52144, joint SFS+LD gives 4,514.50 and 0.52070, and LD alone gives 4,314.16 and
0.54242. The exact SFS result remains primary because it won the frozen recovery
gate. Full hashes, job identifiers, and gate decisions are recorded in
[`experiments/bottleneck/task1_freeze.json`](experiments/bottleneck/task1_freeze.json).
The publication figure is in [`figures/task1_bottleneck/`](figures/task1_bottleneck/).
