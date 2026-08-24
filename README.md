<p align="center">
  <img src="assets/ghist-2026-logo.png" alt="GHIST 2026 logo" width="240">
</p>

# GHIST 2026 challenge suite

Private, reproducible work on all 12 GHIST 2026 challenges. Eleven methods are
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
| 7 | Single sweep with background selection | frozen | chr7:106,416,500–106,576,500; center 106,496,500 bp | 10 localizers span 130.4 kb; independent peak is 5.04× background; transferred untouched observed-like recovery has 18.7 kb median and 45.7 kb maximum error |
| 8 | Multiple sweeps with background selection | frozen | six 120 kb intervals centered at 76,361,767; 88,298,355; 90,666,410; 130,070,549; 138,440,689; and 155,737,529 bp | exactly six peaks exceed the matched threshold; weakest selected μ 23.41 versus next μ 19.34; exon snaps ≤16.1 kb; transferred recovery maximum error 45.7 kb |
| 9 | Relatedness from IBD segments | frozen | 989 testing pairs; official relation and meiosis-count schema | fresh untouched simulation: 82.34% exact-label, 77.81% balanced, and 85.35% meiosis-count accuracy across 12,032 pairs |
| 10 | Large human family | frozen | 26,106 pairs across 229 individuals; official relation-only schema | untouched pair accuracy 45.74%, balanced accuracy 45.58%, and meiosis accuracy 68.37%; two topology searches agree on 98.12% of the candidate |
| 11 | Multifamily cohort with close inbreeding | frozen | 275 within-family dog pairs across eight families; all eight accepted labels represented | family partition ARI 1.0 across 16 runs; untouched pair accuracy 63.58%; 98.18% HBD-threshold stability and zero sibling-component conflicts |
| 12 | Large family with background relatedness | frozen | 50,086 human pairs across 317 individuals; official relation-only schema | 40.94% exact and 60.95% meiosis accuracy without added background; 37.18% exact and 55.20% meiosis accuracy across 12 structured-background holdouts |
| 6 | Soft sweep | active | testing input, accessibility mask, and sealed final secured | locked 48-replicate recovery gate in progress |

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

7. Single sweep with background selection — five nested composite scans, three
   split-flank linkage scales, the largest variant gap, and RAiSD all converge.
   The median center and 80 kb half-width contain all ten localizers. Exact-scale
   hard-sweep recovery transfers cleanly, while the observed independent peak is
   fivefold stronger than the best distant background-selection peak.

8. Multiple sweeps with background selection — one-megabase-separated RAiSD
   peaks were required to exceed twice the matched Task 7 background maximum.
   Consensus centers were snapped to Havana release 104 exons, and the 60 kb
   half-width was transferred from untouched exact-scale sweep recovery.

9. Relatedness from IBD segments — exact Ped-sim calibration uses the supplied
   sex-specific map and crossover-interference model. The selected forest combines
   total sharing with ordered segment shape, binned IBD1/IBD2 spectra, chromosome
   sharing, and sex-map features. Model selection and the final evaluation use
    separate simulation seeds; the candidate follows the official four-column schema.

10. Large human family — eight human-map Ped-sim replicates separate training,
    model selection, and untouched evaluation. The selected 5 cM forest uses
    total sharing and ordered segment shape. Two constrained branch reconstructions
    may change a call only when they agree exactly and preserve its meiosis degree.

11. Multifamily cohort with close inbreeding — stable long-segment graph
    clustering recovers eight families. Dog-map simulations calibrate pairwise
    IBD1/IBD2 shape, union-HBD masking is varied across three thresholds, and
    accepted pedigree tuples plus sibling-component consistency resolve the
    remaining same-degree ambiguities.

12. Large family with background relatedness — human-map Ped-sim calibration
    supplies pairwise long-segment and ordered-shape probabilities. A lower-quartile
    background correction and shortest-background shape trim were selected on 12
    structured holdouts. Branch topology may refine relation type only when its
    meiosis degree agrees with the calibrated pair model.

## Evidence

| Challenge | Frozen record | Figure |
|---|---|---|
| Bottleneck | [`task1_freeze.json`](experiments/bottleneck/task1_freeze.json) | [`figures/task1_bottleneck/`](figures/task1_bottleneck/) |
| Secondary contact | [`task2_freeze.json`](experiments/secondary_contact/task2_freeze.json) | [`figures/task2_secondary_contact/`](figures/task2_secondary_contact/) |
| Growth with background selection | [`task3_freeze.json`](experiments/growth_bgs/task3_freeze.json) | [`figures/task3_growth_bgs/`](figures/task3_growth_bgs/) |
| Split migration with background selection | [`task4_freeze.json`](experiments/split_migration_bgs/task4_freeze.json) | [`figures/task4_split_migration_bgs/`](figures/task4_split_migration_bgs/) |
| Single hard sweep | [`task5_freeze.json`](experiments/single_sweep/task5_freeze.json) | [`figures/task5_single_sweep/`](figures/task5_single_sweep/) |
| Single sweep with background selection | [`task7_freeze.json`](experiments/single_sweep_bgs/task7_freeze.json) | [`figures/task7_single_sweep_bgs/`](figures/task7_single_sweep_bgs/) |
| Multiple sweeps with background selection | [`task8_freeze.json`](experiments/multiple_sweeps_bgs/task8_freeze.json) | [`figures/task8_multiple_sweeps_bgs/`](figures/task8_multiple_sweeps_bgs/) |
| Relatedness from IBD segments | [`task9_freeze.json`](experiments/ibd_segments/task9_freeze.json) | [`figures/task9_ibd_segments/`](figures/task9_ibd_segments/) |
| Large human family | [`task10_freeze.json`](experiments/large_family/task10_freeze.json) | [`figures/task10_large_family/`](figures/task10_large_family/) |
| Multifamily cohort with close inbreeding | [`task11_freeze.json`](experiments/multifamily_inbreeding/task11_freeze.json) | [`figures/task11_multifamily_inbreeding/`](figures/task11_multifamily_inbreeding/) |
| Large family with background relatedness | [`task12_freeze.json`](experiments/large_family_background/task12_freeze.json) | [`figures/task12_large_family_background/`](figures/task12_large_family_background/) |

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
