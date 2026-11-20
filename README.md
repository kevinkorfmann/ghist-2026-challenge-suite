<p align="center">
  <img src="assets/ghist-2026-logo.png" alt="GHIST 2026 logo" width="240">
</p>

# GHIST 2026 challenge suite

Inference methods for the original 12 GHIST 2026 challenges. The two later
low-coverage challenges also have scored final submissions; their scientific
qualification limitations remain documented.

**November 20: all 14 final entries scored.** [Results and participant ranking](docs/ghist_participant_report_2026_11_20.md) · [Detailed Challenge 9 report](docs/task9_method_and_results_2026_11_20.md). Task 9 is second of four on both metrics in the 23:52 Berlin snapshot. Scientific qualification limitations remain documented in the reports below.

[Current readiness and follow-up validation](docs/readiness_2026_11_20.md) ·
[Complete pedigree pipelines and input checks](docs/parallel_followup_2026_11_20.md) ·
[Task 9 package and offline final-run handoff](docs/task9_submission_handoff_2026_11_20.md) ·
[Parallel scientific validation and measured limitations](docs/scientific_validation_2026_11_20.md) ·
[Latest reliability fixes and limitations](docs/reliability_2026_11_20.md) ·
[Initial review and submission blockers](docs/review_2026_11_20.md) ·
[Reproducible setup](docs/setup.md) ·
[Post-freeze testing-answer audit](docs/testing_answer_audit_2026_11_20.json)

The August freezes below are historical testing results. The scientific review
found scaling, simulation/extraction mismatches and validation leakage that
invalidate several historical promotion gates. Their archived metrics are not
current evidence of submission readiness. The old Challenge 9 candidate uses
two retired labels; its corrected candidate now passes the separately versioned
[Task 9 release checks](docs/task9_readiness_2026_11_20.md). This review covers
code and scientific qualification; final-data runs and submission activity are
managed separately and are not established by these validation reports.

[Challenge registry](configs/challenges.yaml) ·
[Validation protocol](docs/validation_protocol.md) ·
[Method comparison](docs/method_comparison.md)

## Historical simulation validation (August 2026)

These are the originally reported metrics, preserved for audit. Subsequent
review found truth-assisted sweep localization and an unavailable clean
reference in background-relatedness evaluation, so the old "held-out" claims
require the qualifications in the linked scientific reports. ↓ is better; ↑ is better.

| # | Challenge | Frozen estimate or output | Held-out closeness |
|---:|---|---|---|
| 1 | [Bottleneck](experiments/bottleneck/task1_freeze.json) · [plots](figures/task1_bottleneck/) | 4,587 generations; fraction 0.52242 | ↓ 13.3% time error; ↓ 0.011 fraction error |
| 2 | [Secondary contact](experiments/secondary_contact/task2_freeze.json) · [plots](figures/task2_secondary_contact/) | N₁ 21,139; N₂ 82,967; split 4,544; contact 1,455; migration 1.99×10⁻⁴ | ↓ 4.78% mean error |
| 3 | [Growth + background selection](experiments/growth_bgs/task3_freeze.json) · [plots](figures/task3_growth_bgs/) | 5,000 generations; 2.0× growth | ↓ 2.38% time; ↓ 0.03% growth error |
| 4 | [Split migration + background selection](experiments/split_migration_bgs/task4_freeze.json) · [plots](figures/task4_split_migration_bgs/) | AFR 16k; EUR 24k; split 4,500; migration 1.5×10⁻⁴ | ↓ 1.36–2.18% sizes/time; ↓ 10.7% migration |
| 5 | [Single hard sweep](experiments/single_sweep/task5_freeze.json) · [plots](figures/task5_single_sweep/) | chr7:38.197–38.317 Mb | ↓ 27.0 kb median; ↑ score 0.625 |
| 6 | [Single soft sweep](experiments/soft_sweep/task6_freeze.json) · [plots](figures/task6_soft_sweep/) | chr7:16.532–16.602 Mb | ↓ 16.3 kb median; ↓ 95.5 kb p90; ↑ score 0.539¹ |
| 7 | [Sweep + background selection](experiments/single_sweep_bgs/task7_freeze.json) · [plots](figures/task7_single_sweep_bgs/) | chr7:106.417–106.577 Mb | ↓ 18.7 kb median; ↓ 45.7 kb max |
| 8 | [Multiple sweeps + background selection](experiments/multiple_sweeps_bgs/task8_freeze.json) · [plots](figures/task8_multiple_sweeps_bgs/) | Six 120 kb intervals | Transferred single-sweep evidence; insufficient for Task 8 |
| 9 | [IBD pairs](experiments/ibd_segments/task9_freeze.json) · [plots](figures/task9_ibd_segments/) | 989 pair labels | ↑ 82.34% exact; ↑ 85.35% degree |
| 10 | [Large human family](experiments/large_family/task10_freeze.json) · [plots](figures/task10_large_family/) | 26,106 pairs; 229 people | ↑ 47.20% exact; ↑ 69.08% degree; ↑ 0.70 points vs previous³ |
| 11 | [Inbred multifamily cohort](experiments/multifamily_inbreeding/task11_freeze.json) · [plots](figures/task11_multifamily_inbreeding/) | 275 pairs; eight families | ↑ 63.58% exact; ↑ family ARI 1.00 |
| 12 | [Family + background relatedness](experiments/large_family_background/task12_freeze.json) · [plots](figures/task12_large_family_background/) | 50,086 pairs; 317 people | ↑ 37.53% exact; ↑ 53.24% degree; ↑ 2.68 points vs previous² |

¹ Signal-matched untouched recovery. ² Empirical short-segment background holdouts. ³ Eight untouched pedigree simulations; all eight improved.

## Frozen test values

<details>
<summary>Show exact locked outputs</summary>

| # | Test value |
|---:|---|
| 1 | generations = 4587.059040; fraction = 0.522420594 |
| 2 | N₁ = 21138.888967; N₂ = 82967.404032; split = 4544.198811; contact = 1455.041806; migration = 0.000199006115 |
| 3 | generations = 5000; growth = 2.0 |
| 4 | AFR = 16000; EUR = 24000; split = 4500; migration = 0.00015 |
| 5 | chr7:38196951–38316951 |
| 6 | chr7:16532298–16602298 |
| 7 | chr7:106416500–106576500 |
| 8 | chr7:76301767–76421767; 88238355–88358355; 90606410–90726410; 130010549–130130549; 138380689–138500689; 155677529–155797529 |
| 9 | 989 pair labels |
| 10 | 26,106 pair labels |
| 11 | 275 pair labels across eight families ([record](experiments/multifamily_inbreeding/task11_freeze.json)) |
| 12 | 50,086 pair labels ([record](experiments/large_family_background/task12_freeze.json)) |

</details>

## Approach

- **Demography (1–4):** frequency spectra, linked-block uncertainty, and held-out simulations.
- **Sweeps (5–8):** multi-scale diversity, haplotype, linkage, variant-gap, and independent RAiSD scans.
- **Relatedness (9–12):** map-matched Ped-sim calibration, cohort-matched segment spectra, and guarded pedigree consistency.

“Frozen” means the estimator, prediction, and validation gates are locked
locally. It does not mean submitted.
