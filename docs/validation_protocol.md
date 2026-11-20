# Validation protocol

## Non-negotiable rules

1. Public testing answers are isolated under `benchmarks/holdout/` and excluded from
   training and model selection.
2. Simulation seeds, parameter blocks, chromosomes, pedigrees, or families define
   validation groups. Random windows or individual pairs from the same generating
   unit never cross splits.
3. Every learned result is compared with at least one mechanistically distinct
   reference method.
4. Calibration, interval coverage, failure strata, and distribution-shift checks are
   reported alongside point accuracy.
5. Final predictions are generated locally and never uploaded automatically.

## Gates

- Synthetic parameter recovery across the full prior and boundary strata.
- Recovery under misspecified mutation/recombination rates and missingness.
- Matched physical window size and blockwise aggregation across the 100 Mb target.
- Linkage-bin occupancy and recovery after streaming extraction.
- Held-out chromosome or family performance.
- Agreement analysis between mechanistic and learned estimators.
- One-time testing-answer audit after all choices are frozen.

## Release checks added November 20, 2026

- Passing unit tests and matching archived checksums are separate from current
  data provenance, scientific accuracy, and final-submission readiness.
- Before a current-data run, record the Synapse entity version and SHA-256. Do
  not silently replace archived inputs or apply legacy 0-based VCF assumptions
  to the corrected October files.
- Validate every candidate with `scripts/validate_submission.py`, specifying the
  challenge and dataset. Supply `--pairs` for a known pair list or `--samples`
  when every sample combination is required. Task 11 instead needs the union of
  within-family pair sets, not all cohort pairs.
- Testing sweep BEDs must use chromosome 7; final BEDs must use chromosome 9.
  Task 9 labels must match the current accepted list. Historical Task 9 freezes
  contain retired classes and are not accepted by the current validator.
- The November public-answer audit is isolated in `benchmarks/holdout/` and is
  used only to report the already frozen predictions. It does not alter fitting,
  thresholds, or method selection. The marked weaknesses require new independent
  validation before any revised estimator is promoted.
- Final outputs require a fresh run on the correct final inputs. Exporting a
  historical testing freeze is never a final inference run.
