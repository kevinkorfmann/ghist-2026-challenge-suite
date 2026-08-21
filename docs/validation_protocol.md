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
