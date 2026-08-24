<p align="center">
  <img src="assets/ghist-2026-logo.png" alt="GHIST 2026 logo" width="240">
</p>

# GHIST 2026 challenge suite

Reproducible solutions for all 12 challenges.

**12/12 methods frozen · 0 submissions · final holdouts sealed**

[Challenge registry](configs/challenges.yaml) ·
[Validation protocol](docs/validation_protocol.md) ·
[Method comparison](docs/method_comparison.md)

## How close are we?

Closeness is measured only on untouched simulations or numerical checks—not on
revealed testing answers. ↓ is better; ↑ is better.

| # | Challenge | Frozen estimate or output | Held-out closeness |
|---:|---|---|---|
| 1 | [Bottleneck](experiments/bottleneck/task1_freeze.json) · [plots](figures/task1_bottleneck/) | 4,587 generations; fraction 0.52242 | ↓ 13.3% time error; ↓ 0.011 fraction error |
| 2 | [Secondary contact](experiments/secondary_contact/task2_freeze.json) · [plots](figures/task2_secondary_contact/) | N₁ 21,139; N₂ 82,967; split 4,544; contact 1,455; migration 1.99×10⁻⁴ | ↓ 4.78% mean error |
| 3 | [Growth + background selection](experiments/growth_bgs/task3_freeze.json) · [plots](figures/task3_growth_bgs/) | 5,000 generations; 2.0× growth | ↓ 2.38% time; ↓ 0.03% growth error |
| 4 | [Split migration + background selection](experiments/split_migration_bgs/task4_freeze.json) · [plots](figures/task4_split_migration_bgs/) | AFR 16k; EUR 24k; split 4,500; migration 1.5×10⁻⁴ | ↓ 1.36–2.18% sizes/time; ↓ 10.7% migration |
| 5 | [Single hard sweep](experiments/single_sweep/task5_freeze.json) · [plots](figures/task5_single_sweep/) | chr7:38.197–38.317 Mb | ↓ 27.0 kb median; ↑ score 0.625 |
| 6 | [Single soft sweep](experiments/soft_sweep/task6_freeze.json) · [plots](figures/task6_soft_sweep/) | chr7:16.532–16.602 Mb | ↓ 16.3 kb median; ↓ 95.5 kb p90; ↑ score 0.539¹ |
| 7 | [Sweep + background selection](experiments/single_sweep_bgs/task7_freeze.json) · [plots](figures/task7_single_sweep_bgs/) | chr7:106.417–106.577 Mb | ↓ 18.7 kb median; ↓ 45.7 kb max |
| 8 | [Multiple sweeps + background selection](experiments/multiple_sweeps_bgs/task8_freeze.json) · [plots](figures/task8_multiple_sweeps_bgs/) | Six 120 kb intervals | ↑ 6/6 peaks; ↓ 45.7 kb max transferred error |
| 9 | [IBD pairs](experiments/ibd_segments/task9_freeze.json) · [plots](figures/task9_ibd_segments/) | 989 pair labels | ↑ 82.34% exact; ↑ 85.35% degree |
| 10 | [Large human family](experiments/large_family/task10_freeze.json) · [plots](figures/task10_large_family/) | 26,106 pairs; 229 people | ↑ 45.74% exact; ↑ 68.37% degree |
| 11 | [Inbred multifamily cohort](experiments/multifamily_inbreeding/task11_freeze.json) · [plots](figures/task11_multifamily_inbreeding/) | 275 pairs; eight families | ↑ 63.58% exact; ↑ family ARI 1.00 |
| 12 | [Family + background relatedness](experiments/large_family_background/task12_freeze.json) · [plots](figures/task12_large_family_background/) | 50,086 pairs; 317 people | ↑ 37.18% exact; ↑ 55.20% degree² |

¹ Signal-matched untouched recovery. ² Structured-background holdouts.

## Approach

- **Demography (1–4):** frequency spectra, linked-block uncertainty, and held-out simulations.
- **Sweeps (5–8):** multi-scale diversity, haplotype, linkage, variant-gap, and independent RAiSD scans.
- **Relatedness (9–12):** map-matched Ped-sim calibration, segment-shape forests, and guarded pedigree consistency.

“Frozen” means the estimator, prediction, and validation gates are locked
locally. It does not mean submitted.
