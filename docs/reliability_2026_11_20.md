# Second reliability pass — November 20, 2026

The code review covers the original twelve implemented challenges. It does not
certify perfect inference or final-submission readiness. Final-data runs remain
deferred at the user's request. The two later low-coverage challenges remain
unimplemented. No predictions were uploaded and no competition submission was
made.

## Corrections

| Challenges | Change and reason |
|---|---|
| 1 — Bottleneck | Reject zero/out-of-range coordinates and incomplete first-site genotypes in fixed-sample block extraction. Preserve continuous predictions so a small positive decline fraction cannot round to zero. Reject nonfinite fits and continue past an individual failed optimizer restart. |
| 2 — Secondary contact | Derive population membership and sizes from the VCF header, including reordered samples and phased/unphased diploids. Reject missing headers, invalid ancestral annotations, duplicate positions, and mixed contigs. Fix the archive-fitting CLI's missing default convention. |
| 3 — Growth + BGS | Remove preset 5,000-generation/2× outputs. Add chromosome-group holdout model comparison, keep fitted estimates continuous, and require the observed model comparison to support the convention used by legacy calibration. Check the full independent recovery grid rather than one favored point. |
| 4 — Split migration + BGS | Export the actual projection median rather than four preset values. Derive sample counts from the AFR/EUR header groups. Validate coordinates and complete genotypes. |
| 1–4 — Spectrum fits | Reject empty, negative, nonfinite, wrong-dimensional spectra; invalid projections; and invalid rate, exposure, restart, or integration settings. Report optimizer-boundary and restart-agreement diagnostics. These diagnostics do not establish scientific accuracy. |
| 5–7 — Single sweeps | Existing scan, coordinate, accessibility, localization, and recovery tests remain in the full suite. This pass does not replace their estimators or claim new large-data validation. |
| 8 — Multiple sweeps + BGS | Remove the exactly-six-peaks freeze and transferred Task 5/7 certification. Add a variable-count candidate detector with a finite-sample threshold from independent matched null chromosome maxima. Require explicit matched multi-sweep recovery before freezing. Merge overlapping BED intervals and clip to chromosome boundaries. |
| 9 — IBD pairs | Retain the prior ten-label correction, canonical pair validation, and recorded model evaluation; covered again by the full suite. |
| 10–12 — Family relatedness | Check exact unordered pair coverage against actual VCF sample IDs; Task 11 checks the partition and all within-family pairs. Remove preset cohort/family/pair counts and the requirement that every permitted dog relationship must appear. |

Fixed-sample BGS extractors now reject missing genotypes rather than silently
mixing sample sizes or retaining full physical exposure after filtering calls.
A missing-data analysis requires a projected spectrum and corresponding callable
exposure. The existing generic projected extractors remain available.

The secondary-contact simulator documentation now correctly describes migration
as backward-time lineage rates, following the
[msprime API](https://tskit.dev/msprime/docs/stable/api.html#msprime.Demography.set_migration_rate).
The numerical model itself was not changed by this terminology correction.

Version 0.1.1 and the lockfile force editable installations to discover the new
shared inference package. Refresh existing environments with the setup command:

```bash
uv sync --locked --python 3.12 --extra dev --extra references --extra figures
```

## Verification

All **188 tests pass** in the clean locked macOS and Linux ARM64 environments,
and in Betty's original model environment (job **10848387**, exit 0). The
additional 35 tests cover the failures above, model recovery at an independent
parameter point, held-out chromosome isolation, sample reordering, family
membership, and variable sweep counts. Critical Ruff checks, 151 shell syntax
checks, and ten installed CLI help entry points pass. The committed
[verification records](../experiments/reliability_2026_11_20/verification.json)
include source hashes, logs, and simulation job IDs.

The new macOS environment is outside the synced Documents directory because
its hidden-file flags repeatedly broke Python's editable import paths. The old
environment is preserved, `.venv` points to the working replacement, and a fresh
process verifies the installed entry points. Betty still has dependency
warnings from its older environment; these did not fail the suite.

## Independent growth checks

The [chromosome holdout record](../experiments/reliability_2026_11_20/growth_chromosome_cv.json)
uses the archived, unlabeled middle-distance block spectra, leaving one complete
chromosome arm out at a time. Public testing answers were not used to fit or
choose models. At projection 40, summed held-out log likelihood is -1,076.211 for
exponential growth and -1,065.712 for sudden growth. Sudden growth wins three of
four arms. Its full-data candidate is approximately 2,840.94 generations and
1.81584× growth. **This candidate is not frozen**: disagreement across arms and
BGS recovery uncertainty remain.

The new SLiM smoke checks independently generate exponential and sudden growth
at an arbitrary new parameter point: ancestral size 8,000, duration 7,600
generations, growth ratio 3.2. They use exonic Gamma_H17 background selection,
the Comeron map, a 1 Mb segment, 20 diploid samples, independent seeds, scaling
factor 10, and burn-in 1. Results and execution records are under
[experiments/reliability_2026_11_20](../experiments/reliability_2026_11_20/).
The exponential-generated run had 12.25% timing error and 1.40% growth-ratio
error; the sudden-generated run had 2.56% timing error and 6.27% growth-ratio
error when each was fitted with its matching model. Even on the exponential
simulation, sudden growth had a slightly higher likelihood (0.57 log units).
This reinforces the convention ambiguity rather than establishing a winner.
Two small runs are diagnostic checks, not a recovery grid or evidence of final
accuracy. SLiM's scaling approximation also requires broader stress evaluation.

## Changed entry points

- `python -m scripts.compare_growth_models --blocks INPUT.npz --out comparison.json`
  performs chromosome holdouts. `freeze_growth_bgs.py` now requires
  `--model-comparison`; legacy exponential calibration cannot certify a sudden
  model or inconsistent chromosome evidence.
- `python -m scripts.simulate_growth_bgs_selection --convention sudden ...`
  generates sudden growth and records recovery for its matching fitted model.
  The default exponential convention remains available.
- `python -m scripts.detect_multiple_sweeps --help` lists the required matched
  null reports, context, reference exons, chromosome length, and localization
  settings. At alpha 0.05 at least 19 independent null maxima are required.
  The error-rate interpretation assumes exchangeable whole-chromosome null
  scans. Merely declaring matching context is not proof that simulations match.
- `python -m scripts.freeze_task8_multiple_sweeps --candidate CANDIDATE.json
  --recovery RECOVERY.json --dataset testing --bed OUTPUT.bed --out FREEZE.json`
  requires recovery identifying the exact candidate SHA-256, independent
  evaluation without public answers, matched bottleneck/growth+BGS+multiple
  sweeps, and precision/recall/coverage meeting prespecified minimums. There is
  no implicit testing path or fixed sweep count. Empty detections cannot be
  converted into a fabricated submission interval.
- Challenge 10 and 12 historical freeze entry points now require
  `--sample-vcf`. Challenge 11 uses its explicitly named testing VCF. These
  historical testing writers are not final inference workflows.

The old `select_background_sweep_candidates.py` remains a **provisional
historical heuristic**. Its Task 7 background threshold and fixed separation
must not be treated as matched Task 8 validation. Archived freeze JSON and
figures are preserved as historical records; their old `passed` fields do not
supersede this review.

## Submission readiness

The [first review](review_2026_11_20.md) still describes the unresolved scientific
and data requirements for all twelve challenges. In particular, Challenge 2's
hidden-history assumptions, Challenge 3's growth/BGS ambiguity, and Challenge
8's matched-null sensitivity and localization need stronger independent
validation. Family relationship accuracy also remains probabilistic and limited
by the recorded holdout results. New guards prevent unsupported freezes; they
are not substitutes for better estimators or validated current final runs.
