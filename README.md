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
| Secondary Contact | Task 2 frozen, not submitted | corrected 40-by-30 hidden-epoch SFS estimate: 21,139 and 82,967 individuals, split 4,544 generations, contact 1,455 generations, migration 0.000199006 |
| Growth with background selection | Task 3 frozen, not submitted | 5,000 generations and two-fold exponential growth; projection, linked-block, broad recovery, and selection-strength stress gates pass |
| Split migration with background selection | final numerical gate running | leading candidate 16,000 / 24,000 individuals, 4,500-generation split, and 0.00015 symmetric migration; projection, recovery, arm-holdout, and linked-block gates pass |
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

## Task 2 frozen result

The Secondary Contact generator-matched spectrum is now the mechanistic
reference. On a known-parameter simulation it reduced competition RRMSE from
3.195 for the constant-size approximation to 0.309. On the observed unfolded
spectrum it improves multinomial log likelihood by about 3,851 units relative
to that approximation. Folded and unfolded fits agree closely within each
resolution.

The resolution gate is closed. On 32 matched finite-site chromosomes, the
unfolded 40-by-30 projection has mean and median competition RRMSE of 0.1050 and
0.0989, compared with 0.1091 and 0.1030 at 50-by-35. It wins 26 of 32 paired
replicates; its mean score advantage is 0.00408 with a paired bootstrap interval
of 0.00187–0.00620. Eight-replicate comparisons also place it ahead of 60-by-40
and the full sample spectrum. An independent 48-chromosome comparison across 24
varying demographic settings also selects 40-by-30: it wins 38 of 48 pairs,
with mean RRMSE 0.1106 versus 0.1143 and a paired mean-difference interval of
−0.00514 to −0.00214. The refined raw fit is 20,127.66 and 79,105.52 for
the current population sizes, 4,346.33 generations for the split, 1,395.75
generations for contact, and 0.000207236 for migration. Tighter adaptive
integration changes every raw parameter by far below 0.01%, and independent
refined finite-site cells remain within 0.1% of the standard integration.

A point-grouped transfer gate spans 24 nearby demographic settings with two
100 Mb chromosomes per setting: 14 settings for calibration, five for method
selection, and five untouched for evaluation. A score-selected log-ridge
correction, anchored by two independent 32-chromosome center campaigns, lowers
untouched mean RRMSE from 0.1172 to 0.0478 and median RRMSE from 0.1073 to
0.0397. All five parameter-wise median errors improve and fall between 0.81%
and 1.81%. Its current accepted estimate is 21,138.89 and 82,967.40 for the
population sizes, 4,544.20 and 1,455.04 generations for split and contact, and
0.000199006 for migration.

The full-length stress gate contains eight replicates each under nominal,
alternating 0.25×/4× recombination, 0.8×/1.2× mutation rate,
0.9×/1.0× local size, 2:1 directional migration, and 0.1% polarization
error. No optimizer boundaries occurred. The correction wins all 40 replicates
across nominal, recombination, mutation, local-size, and polarization conditions,
reducing mean RRMSE to 0.033–0.083. Under deliberately extreme directional
migration, both symmetric-model estimates are poor and correction changes mean
RRMSE from 0.321 to 0.330; this limitation is retained rather than hidden. The
challenge generator specifies symmetric migration, so the primary estimate is
unchanged.

Two learned spectrum representations improved untouched synthetic evaluation,
including a variance-stabilized ensemble with RRMSE 0.507 and median absolute
relative errors of 3.9–11.9% across parameters. Both produced observed
predictions inconsistent with the mechanistic fit and are therefore rejected
as primary estimators. The complete 100 Mb SINGER and POLEGON time-sliced
genealogy campaign is retained as an independent check. Its current held-out
results support population-size and split-time scale but show that marginal
coalescence hazards alone are too weak to identify the contact–migration pair.
The accepted prediction, source hashes, job identifiers, and all gates are in
[`experiments/secondary_contact/task2_freeze.json`](experiments/secondary_contact/task2_freeze.json),
with publication figures in
[`figures/task2_secondary_contact/`](figures/task2_secondary_contact/). No result
has been uploaded or submitted.

## Task 4 leading result

The leading Split Migration with background selection estimate is 16,000 AFR and
24,000 EUR individuals, a split 4,500 generations ago, and symmetric migration
of 0.00015 per generation. The raw median across three joint-spectrum
projections is 16,192 / 22,949 individuals, 4,388 generations, and 0.000146697.
The accepted values are generator-level rounding supported by a matched
four-arm selection simulation, not a fitted correction.

Across six held-out truth points, the uncorrected estimator beats every tested
calibration. A chromosome-stratified 64-replicate 1 Mb block bootstrap contains
the accepted candidate in all four 95% intervals, while leave-one-arm-out fits
remain compatible with the whole-genome estimate. The final dataset stays
sealed and nothing has been uploaded or submitted. One stricter numerical-
integration gate remains before the result is frozen. The validation figure is in
[`figures/task4_split_migration_bgs/`](figures/task4_split_migration_bgs/).
