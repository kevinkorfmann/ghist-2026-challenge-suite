# Scored submission closeout — November 20, 2026

**Final status:** all 14 tasks are now SCORED. See the [consolidated final summary](final_scored_summary_2026_11_20.md) and its machine audit. The chronology below preserves intermediate states, including invalid attempts that were subsequently corrected.

Tasks 11 and 14 are **SCORED**. The coordinator verified the server score emails and reported the following receipts. This closeout records that evidence; the submission screenshots show the success popup and uploaded entity, rather than the scores themselves.

| Task | Submission | Queue | Uploaded entity | Server score email time (UTC) | Score |
|---|---|---|---|---|---|
| 11 | 9785361 | 9619676 | syn78780263, version 1 | 21:35:51 | relation 0.42540795348737687; meioses_count 0.5724422540498211 |
| 14 | 9785360 | 9620086 | syn78780264, version 1 | 21:34:50 | RRMSE 0.042643298573969445 |

Task 11 was submitted as “Task11 corrected family context 2026-11-20”; score email `1a1fae3204484254`, scoring logs `syn78780307`. Task 14 was submitted as “Task14 likelihood spectrum 2026-11-20”; score email `1a1fae233fe71b87`, scoring logs `syn78780293`. Relatedness scores use the challenge's truth-class-size-weighted Jaccard, not ordinary accuracy.

The exact submitted file hashes were independently rechecked against each readiness record and upload record. Task 11 CSV SHA-256 is `8ab80bd7d7c39a7ec95de49290b71ab1edcfa53c206ee431386625fe4fea951f`; Task 14 YAML SHA-256 is `350de6e7b9f48463c3bab3f4ea2b2aa266530d80ae3f17cc0cb7b9ea94f2ccbe`. Models and predictions were not changed after scoring. Historical `UPLOADED_NOT_SUBMITTED` upload records remain preserved and are superseded by the new `SCORED` receipts.

The repository copies of [Task 11 receipt](../experiments/final_submission_2026_11_20/task11/submission_receipt.json), [Task 11 score](../experiments/final_submission_2026_11_20/task11/score_record.json), [Task 14 receipt](../experiments/final_submission_2026_11_20/task14/submission_receipt.json), and [Task 14 score](../experiments/final_submission_2026_11_20/task14/score_record.json) bind the candidate, readiness, original upload record, and submission screenshot hashes. Identical copies and screenshots are staged in `/private/tmp/ghist-final-submissions-20261120/task11/` and `/private/tmp/ghist-final-submissions-20261120/task14/`.

Submission success does not change the recorded scientific qualification results. Task 11's corrected context candidate failed strict qualification; reserved stress and confirmation cohorts were not run. Task 14 passed its final numerical checks, while earlier strict development failures and absence of fresh cap-500 confirmation remain documented. Both were submitted under the user's deadline-first authorization.

At the coordinator's closeout report, Task 10's first submission was **INVALID**, and correction was underway. It is not marked complete here. The coordinator continues the remaining submissions and owns all submission UI actions. This note requires no additional computation, tuning, or submission action for Tasks 11 and 14.

This resident continuation note supplements the existing continuation documents. Reading `continuation_review_2026_11_20.md` encountered an iCloud dataless hydration error during closeout; its existing content was preserved.

## Subsequent email confirmations

Read-only Gmail retrieval independently confirmed the Task 11 and 14 scores above. It also confirmed Task 13 submission 9785362 as SCORED at 21:42:26 UTC, RRMSE 0.6172118592711575. The schema-corrected Task 10 retry 9785368 is SCORED at 21:46:36 UTC: relation 0.12225889450060641 and meioses_count 0.1283610601352386. Original invalid attempts remain preserved. No models or predictions were changed by this monitoring work.

Task 5 submission 9785369 scored recall 0, size 1.2, F1 0; Task 7 submission 9785371 scored recall 1, size 1.6, F1 0.5555555555555556; Task 8 submission 9785372 scored recall 0.06060606060606061, size 2.4, F1 0.10695187165775401. Task 6 was still processing at this check.

Tasks 1–4 first submissions 9785364–9785367 were INVALID because the server could not find a `parameters` section in their YAML. The coordinator assigned schema-only corrections, preserving parameter values. Task 12 first submission 9785363 was INVALID because its queue rejected two relationship labels; its correction is separately owned. The compact [email monitor](../experiments/final_submission_2026_11_20/score_email_monitor.json) records exact email IDs, submission IDs, queues, scores and invalid reasons, and is updated as later receipts arrive.

At 21:47:17 UTC, Task 6 submission 9785370 scored recall 0, size 1.5, F1 0. At 21:49:43 UTC, the schema-corrected Task 12 retry 9785373 scored relation 0.11235515689383353 and meioses_count 0.13564060673696537. These valid receipts supersede the pending states above; initial invalid receipts remain retained.
