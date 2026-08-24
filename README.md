<p align="center">
  <img src="assets/ghist-2026-logo.png" alt="GHIST 2026 logo" width="240">
</p>

# GHIST 2026 challenge suite

Private, reproducible work on all 12 GHIST 2026 challenges. Five methods are
frozen and the remaining inputs are secured.
Nothing has been submitted. Final holdouts remain sealed.

[Challenge registry](configs/challenges.yaml) ·
[Validation protocol](docs/validation_protocol.md) ·
[Method comparison](docs/method_comparison.md)

## Results

“Closeness” below always refers to held-out simulations or numerical checks,
never to a revealed testing answer.

| # | Challenge | State | Local result | Closeness on validation |
|---:|---|---|---|---|
| 1 | Bottleneck | frozen | 4,587 generations; post-decline fraction 0.52242 | 13.3% median time error and 0.0113 median fraction error on 29 untouched chromosomes; physical-block 95% intervals: 4,184–4,982 and 0.5131–0.5308 |
| 2 | Secondary contact | frozen | N₁ 21,139; N₂ 82,967; split 4,544 generations; contact 1,455 generations; migration 0.000199006 | 4.78% mean error score on untouched demographic settings; parameter-wise median errors 0.81–1.81% |
| 3 | Growth with background selection | frozen | 5,000 generations; growth factor 2.0 | matched recovery errors 2.38% and 0.03%; 128-block 95% intervals: 3,984–5,913 and 1.851–2.122 |
| 4 | Split migration with background selection | frozen | AFR 16,000; EUR 24,000; split 4,500 generations; migration 0.00015 | matched recovery errors 1.36–2.18% for sizes/time and 10.7% for migration; strict integration changed every parameter by <0.0006% |
| 5 | Single hard sweep | frozen | chr7:38,196,951–38,316,951; center 38,256,951 bp | untouched exact-scale recovery: 27.0 kb overall median error; all 12 observed-like cells hit the interval; mean challenge score 0.625 |
| 6–12 | Remaining sweep and relatedness challenges | queued | testing inputs, masks, maps, templates, and sealed finals secured | method development not yet frozen |

Frozen means that the estimator, prediction, and validation gates are locked
locally. It does not mean submitted.

## How each result was obtained

1. Bottleneck — folded frequency-spectrum inference was selected as primary.
   Linkage decay, independent spectrum software, genealogy summaries, and a
   chromosome-held-out context model were used as cross-checks. Physical-block
   bootstrap intervals quantify spatial uncertainty.

2. Secondary contact — a generator-matched two-population spectrum was fitted
   at resolution selected on independent chromosomes. A point-grouped held-out
   calibration reduced the error score from 11.7% to 4.78%, followed by
   recombination, mutation, polarization, and migration-direction stress tests.

3. Growth with background selection — distance-stratified unfolded spectra used
   the supplied coding annotation and recombination map. Projection stability,
   linked-block bootstrap, known-truth forward simulations, and selection-strength
   perturbations all support the rounded generator-level result.

4. Split migration with background selection — three joint-spectrum projections
   were combined over the well-powered noncoding distance band. Six held-out
   truth points favored no calibration. Arm holdouts, a 64-replicate linked-block
   bootstrap, matched linked-selection recovery, and stricter integration all pass.

5. Single hard sweep — nested diversity, frequency-spectrum, haplotype, linkage,
   and variant-gap scans converge on the same region; RAiSD provides an independent
   reference scan. The candidate center is the median of three fine localizers.
   Interval width is selected against the published scoring rule. The locked
   60 kb half-width passed the untouched exact-scale recovery gate.

## Evidence

| Challenge | Frozen record | Figure |
|---|---|---|
| Bottleneck | [`task1_freeze.json`](experiments/bottleneck/task1_freeze.json) | [`figures/task1_bottleneck/`](figures/task1_bottleneck/) |
| Secondary contact | [`task2_freeze.json`](experiments/secondary_contact/task2_freeze.json) | [`figures/task2_secondary_contact/`](figures/task2_secondary_contact/) |
| Growth with background selection | [`task3_freeze.json`](experiments/growth_bgs/task3_freeze.json) | [`figures/task3_growth_bgs/`](figures/task3_growth_bgs/) |
| Split migration with background selection | [`task4_freeze.json`](experiments/split_migration_bgs/task4_freeze.json) | [`figures/task4_split_migration_bgs/`](figures/task4_split_migration_bgs/) |
| Single hard sweep | [`task5_freeze.json`](experiments/single_sweep/task5_freeze.json) | [`figures/task5_single_sweep/`](figures/task5_single_sweep/) |

## Layout

```text
challenges/   task definitions and submission schemas
configs/      dataset registry and frozen metadata
experiments/  concise evidence and freeze records
figures/      publication-quality validation figures
scripts/      audits, fitting, simulation, validation, and freezing
slurm/        Betty compute workflows
src/          reusable inference code
tests/        deterministic unit and contract checks
```

Heavy computation runs on Betty through SLURM. Source revisions, dataset hashes,
job identifiers, validation decisions, and submission count are recorded in each
freeze file.
