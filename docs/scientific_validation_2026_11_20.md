# Scientific validation of the original 12 GHIST challenges

The review found and repaired material scientific implementation defects. New
independent simulations also exposed estimation failures that software tests
cannot resolve. **None of the twelve complete challenge pipelines is certified
ready for final submission by this review.** Final-data runs remain deferred;
no submission or replacement of a historical freeze was made. The two later
low-coverage challenges remain unimplemented.

Three parallel reviewers covered Tasks 2–4, 5–8 and 9–12. The coordinating review
covered Task 1, integration and evidence reconciliation. A second reviewer also
audited the Task 1 changes and found two further defects, both fixed: leftover
model files could bypass new normalization provenance, and larger diagnostic
panels could reuse random seeds. No public challenge answers were used to fit,
calibrate or select methods during this scientific review.

## Per-challenge findings and release gates

| Task | What the review established | What still prevents certification |
|---|---|---|
| 1 — Bottleneck | Corrected diploid population scaling; repaired the dadi optimizer on current SciPy; rejected incompatible learned artifacts. On 12 matched 10 Mb chromosomes, median moments timing error was 19.0%, with a severe weak-decline failure exceeding 3,900%. | Regenerate intended-prior learned artifacts and calibration; assess full-length recovery, identifiability and interval coverage. Solver agreement does not remove shared model bias. |
| 2 — Secondary contact | Retained between-population fixed differences and corrected hidden demographic epochs. Independent genealogies and diffusion agree for the specified history. Some asymmetric-history fits have large contact-time errors. | Regenerate training/calibration; assess plausible current histories and optimizer robustness on independent complete chromosomes. |
| 3 — Growth with BGS | Corrected physical-unit SLiM scaling and likelihood-domain handling. Unscaled growth simulation has +27.5% timing error and weakly favors the wrong growth convention. | Replicated model discrimination, selection/map shifts, longer burn-in, full exposure and uncertainty coverage. |
| 4 — Split/migration with BGS | Corrected physical scaling, removed truth initialization, and validated all observed likelihood bins. Corrected smoke runs still have substantial parameter errors. | Broad independent BGS recovery and calibration, projection stability, full-chromosome coverage. |
| 5 — Hard sweep | Removed truth-centered fine searches; repaired RAiSD task/coordinate matching. Recomputed archived recovery gives 25.4 kb median and 118.6 kb p90 error on its last 64 cases. | That reevaluation is retrospective, not a new untouched promotion set. Broad matched independent localization coverage remains required. |
| 6 — Soft sweep | Corrected localization exposes a p90 error around 2.19 Mb for the recalibrated retrospective composite. Its gate fails. | New independent soft-sweep recovery and a validated localizer; do not reuse the historical accuracy claim. |
| 7 — Sweep with BGS | Freeze requires matched BGS evidence. In eight unscaled single-sweep cases, the ensemble median error is 37.9 kb and only 5/8 are within 80 kb. | Broad matched BGS localization evidence and coverage; the small one-region panel is insufficient. |
| 8 — Multiple sweeps with BGS | New unscaled scan detects 8/16 sites across eight two-sweep segments, with 8/11 reported regions hitting a site. | Whole-chromosome matched-null calibration, count/completeness, false positives and localization coverage. Segment-scale nulls cannot certify chromosome-wide error. |
| 9 — IBD pairs | New data exposed undefined probabilities from retired labels. A ten-label refit on original calibration data has no such failures on two untouched confirmation sets; exact accuracy is 89.1% matched and 87.2% shifted. | Candidate remains unpromoted. Validate current input, zero-IBD pairs and biological shifts; similar accuracy is not an improvement claim. |
| 10 — Large human family | Independent exact-IBD template accuracy averages 47.5% matched and 44.8% under Poisson crossover shift. Split guards and grouped evaluation are strengthened. | Simulate complete families through genotypes, IBD detection and pedigree consensus; pair-template accuracy does not validate the complete estimator. |
| 11 — Inbred multifamily cohort | Independent dog-map exact-IBD template accuracy averages 63.2%; plotting handles absent classes and actual family counts. | Matched inbreeding/founders, HBD masking, family discovery and full pedigree inference with entire families held out. |
| 12 — Family with background relatedness | Removed the unavailable clean evaluation reference and corrected split/blend handling. Independent-reference contamination checks average 38.5% and 36.4% exact accuracy. | Retrain with the corrected reference protocol; validate correlated background genealogy and complete families through detected IBD. |

The demographic panels use linked chromosomes; relatedness intervals resample
whole pedigrees. Repeated parameter settings, map perturbations and contamination
regimes are not counted as independent biological samples. Small diagnostic
panels cannot certify full-prior accuracy, tail risk or nominal interval coverage.
Legacy artifacts remain preserved; repaired provenance guards deliberately
require regeneration where the old evidence is scientifically incompatible.

The neutral Tasks 2–4 panel completed 13 of 22 planned datasets within a
20-minute diagnostic budget. Seven Task 2 and two Task 4 cases remained
unfinished; the unscaled BGS Task 4 fit also reached its separate 15-minute
limit after simulation completed. These are recorded as time-limited outcomes,
not successful recoveries. Their absence can systematically bias completed-only
accuracy summaries, which therefore cannot describe the full planned panel.

## Detailed evidence

- [Task 1 report](scientific_bottleneck_2026_11_20.md) and [records](../experiments/scientific_2026_11_20/bottleneck/).
- [Tasks 2–4 report](scientific_demography_2026_11_20.md) and [records](../experiments/scientific_2026_11_20/demography/).
- [Tasks 5–8 report](scientific_sweeps_2026_11_20.md) and [records](../experiments/scientific_2026_11_20/sweeps/).
- [Tasks 9–12 report](scientific_relatedness_2026_11_20.md) and [records](../experiments/scientific_2026_11_20/relatedness/).
- [Software verification](../experiments/scientific_2026_11_20/verification.json), including source hashes and platform-specific test logs.

Scientific experiments ran on Betty through Slurm and on the separately
authorized Linux host. These were independent validation simulations, not
current final-data inference runs. No diagnostic failure was converted into a
success by discarding it, retuning on its truth or inserting a fallback answer.

The final regression suite passes **241 tests** on macOS/Python 3.12, Linux
ARM64/Python 3.12 and Betty/Python 3.11. Ruff, Python compilation and syntax
checks for all 18 changed/new Slurm scripts also pass. These checks establish
the exercised software behavior; they do not certify scientific accuracy.

## Duplicate Claude reviews

Three read-only Claude reviews were attempted against the baseline source.
All terminated at the account's five-hour usage limit before producing a final
report. Their prompts and failure metadata are retained in
[the attempt record](../experiments/scientific_2026_11_20/claude/attempts.json).
These attempts do **not** count as completed independent reviews. Access
restoration is pending; the parallel code and scientific checks above are
independent of that limitation.

## Before any submission

Regenerate invalidated training/calibration and pass each task's independent
release gates before promoting new estimators. Once final runs resume, record
current input versions and hashes, execute the complete pipelines, validate
every submission schema, and inspect diagnostics. Competition submission is a
separate action and has not been performed or scheduled.
