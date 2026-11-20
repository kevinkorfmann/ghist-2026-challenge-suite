# Readiness continuation ownership — November 20, 2026

The user requested continued work until readiness, while avoiding work already
owned by other chats. The active chat **Review 2026 ghist challenges** is
integrating existing shared repairs and owns Tasks 1–8, Task 9 packaging and
Task 12. The **Delegate three unassigned tasks** chat completed its independent
reviews and is monitoring readiness notifications. Their existing files and
experiment snapshots remain outside this continuation's editing scope.

The **Check submission status** chat owns the following new work areas:

| Agent | New work | Owned files/directories |
|---|---|---|
| `task10_readiness_development` | Diagnose Task 10 whole-family accuracy and develop a separately versioned candidate | `scripts/develop_task10_detected_ibd.py`, `tests/test_task10_detected_ibd.py`, `docs/task10_development_2026_11_20.md`, `experiments/task10_development_2026_11_20/` |
| `task11_readiness_development` | New dog founder/pedigree simulations, pair-first development and untouched confirmation | `scripts/develop_task11_pedigree.py`, `tests/test_task11_pedigree_development.py`, `docs/task11_development_2026_11_20.md`, `experiments/task11_development_2026_11_20/` |
| `lowcoverage_readiness` | Tasks 13–14 BAM likelihood pipeline and independent synthetic validation | new `scripts/lowcoverage_*.py`, `src/ghist_lowcoverage/`, `tests/test_lowcoverage_*.py`, `docs/lowcoverage_readiness_2026_11_20.md`, `experiments/lowcoverage_2026_11_20/` |
| Coordinator | Independent cohort scoring and final cross-review | `scripts/score_relatedness_cohorts.py`, `tests/test_relatedness_cohort_scoring.py`, `experiments/relatedness_qualification_2026_11_20/`, this document |

All remote computations use separate source snapshots and Slurm work areas.
There are no edits to previous result files, shared models, existing production
scripts or release defaults in this continuation. At 19:13 UTC the user authorized
starting submissions for the strongest working candidates while development
continues. The coordinator owns the new Task 9 final-input run and Synapse upload/
submission workflow, using the other chat's completed v2 package unchanged.
The separate run root is `final_runs/task9-root-20261120-1914/` on Betty;
local working files are under `/private/tmp/ghist-final-submissions-20261120/task9/`.
The existing submission-access chat is idle; the main review chat continues
shared integration and Tasks 1–8/12. Task 14 final work remains conditional on
confirmation and independent review. A changed candidate must be frozen
before independent confirmation outcomes are opened; unsuccessful confirmation
cannot be turned into success by lowering its previously recorded criteria.

The coordinator's independent scorer retains missing/abstained pairs in the
accuracy denominator and fails coverage on extra cross-family pairs. It reports
every accepted label, including absent classes, and resamples whole cohorts
for descriptive uncertainty. It does not compute an all-label interval by
silently dropping bootstrap samples lacking a class. Calibration uses the
probability of the actual emitted label and is unavailable if any expected
pair lacks a resolved prediction or a complete probability distribution.
These metrics do not establish biological independence or authorize a release.

The scorer accepts JSON containing `labels` and `cohorts`, each cohort containing
`cohort_id`, `truth` and `predictions`. Pair rows have `indID1`, `indID2` and
`relation`; prediction rows may also contain the complete `probabilities`
mapping. Run with `scripts/score_relatedness_cohorts.py --input INPUT --out NEW_OUTPUT`.
Its initial 18 focused tests and Ruff checks pass.


## Scored closeout, 21:44 UTC

Tasks 11 and 14 are SCORED (submissions 9785361 and 9785360). Exact scores, receipt and candidate hashes, evidence scope, and preserved scientific limitations are recorded in [the scored closeout](continuation_scored_closeout_2026_11_20.md). Task 10's first submission was INVALID and is being corrected; the coordinator continues remaining submissions.


## Final scored closeout, 22:01 UTC

All 14 tasks now have valid server scores. [The consolidated scored summary](final_scored_summary_2026_11_20.md) lists every valid submission and score; its linked machine audit verifies exact candidates and email evidence. All current receipts are SCORED, and six initial invalid attempts remain preserved separately. Scientific qualification records remain unchanged. No further submission is required for these 14 tasks.
