# All 14 tasks scored — November 20, 2026

All 14 GHIST final tasks have valid server scores. The final outstanding score emails arrived by 21:56:24 UTC, before the 23:59 UTC deadline. The [consolidated machine audit](../experiments/final_submission_2026_11_20/all_tasks_scored_audit.json) verifies 14 unique scored submission IDs, exact score emails, all submitted candidate SHA-256/MD5 hashes, matching temporary/repository records, and 93 referenced artifact hashes. Six initial invalid attempts remain preserved.

| Task | Valid submission | Server score |
|---|---|---|
| 1 | 9785374 | RRMSE 0.733000615505 |
| 2 | 9785375 | RRMSE 2.13767297778 |
| 3 | 9785376 | RRMSE 0.796883521581 |
| 4 | 9785377 | RRMSE 0.516921092855 |
| 5 | 9785369 | recall 0; size 1.2; F1 0 |
| 6 | 9785370 | recall 0; size 1.5; F1 0 |
| 7 | 9785371 | recall 1; size 1.6; F1 0.555555555556 |
| 8 | 9785372 | recall 0.0606060606061; size 2.4; F1 0.106951871658 |
| 9 | 9785353 | relation 0.833706007793; meioses 0.842384548401 |
| 10 | 9785368 | relation 0.122258894501; meioses 0.128361060135 |
| 11 | 9785361 | relation 0.425407953487; meioses 0.57244225405 |
| 12 | 9785373 | relation 0.112355156894; meioses 0.135640606737 |
| 13 | 9785362 | RRMSE 0.617211859271 |
| 14 | 9785360 | RRMSE 0.042643298574 |

The JSON audit retains full floating-point values and each email ID, entity/version, queue, candidate hash, receipt path and score-record path. Relatedness scores are the challenge’s truth-class-size-weighted Jaccard, not ordinary accuracy. RRMSE is an error measure; these different metric families should not be compared as a single score.

Tasks 1–4 were retried after adding the official YAML `parameters` wrapper with parameter values unchanged. Tasks 10 and 12 were retried after the server rejected relationship labels; the recorded same-meiosis label projection preserves pair IDs and implied meioses but merges biologically different relationships. Initial invalid receipts are separate from valid scored receipts. Task 9’s original pending UI receipt and rounded dashboard score remain preserved; its current receipt is SCORED and `score_record_exact.json` records the exact email result.

Server acceptance does not erase scientific limitations. The original release failures, incomplete confirmation, optimizer convergence limits and synthetic-to-final transfer limitations remain in the bound readiness and method records. These entries were submitted under the user’s deadline-first authorization. This closeout made no changes to models, parameters or predictions and did not use final scores to choose a new candidate.

Current records are under `/Users/kevinkorfmann/Documents/ChatGPT/ghist/experiments/final_submission_2026_11_20/taskN/` and `/private/tmp/ghist-final-submissions-20261120/taskN/`. The coordinator owns public communication and any remaining administrative actions; no further submission is required for these 14 tasks.

## Current relative results and delivery

These comparisons use visible official dashboard rows checked between 21:45 and 22:05 UTC on November 20, 2026. They are current per-task comparisons, not a final overall ranking. Task 9 leads all three displayed entries on both metrics; the next best displayed scores are relation 0.828 and meioses 0.840. Most other ranked entries trail the leaders. No accessible public rank was found for Tasks 13 and 14. RRMSE is lower-is-better; F1 and relatedness weighted Jaccard are higher-is-better.

| Task | Score (rounded) | Current displayed rank |
|---|---|---|
| 1 | RRMSE 0.733 | [6/6](https://www.synapse.org/Synapse:syn74803798/wiki/640737) |
| 2 | RRMSE 2.1377 | [3/4](https://www.synapse.org/Synapse:syn74803798/wiki/640747) |
| 3 | RRMSE 0.79688 | [3/3](https://www.synapse.org/Synapse:syn74803798/wiki/640742) |
| 4 | RRMSE 0.51692 | [3/3](https://www.synapse.org/Synapse:syn74803798/wiki/640750) |
| 5 | F1 0.000 | [Tied last (five-way tie among six)](https://www.synapse.org/Synapse:syn74803798/wiki/640748) |
| 6 | F1 0.000 | [Tied last (two-way tie among six)](https://www.synapse.org/Synapse:syn74803798/wiki/640744) |
| 7 | F1 0.556 | [4/6](https://www.synapse.org/Synapse:syn74803798/wiki/640749) |
| 8 | F1 0.107 | [5/6](https://www.synapse.org/Synapse:syn74803798/wiki/640745) |
| 9 | Relation 0.834; meioses 0.842 | [1/3 on both](https://www.synapse.org/Synapse:syn74803798/wiki/640981) |
| 10 | Relation 0.122; meioses 0.128 | [2/2 on both](https://www.synapse.org/Synapse:syn74803798/wiki/640982) |
| 11 | Relation 0.425; meioses 0.572 | [2/2 on both](https://www.synapse.org/Synapse:syn74803798/wiki/640984) |
| 12 | Relation 0.112; meioses 0.136 | [2/2 on both](https://www.synapse.org/Synapse:syn74803798/wiki/640983) |
| 13 | RRMSE 0.61721 | Not available |
| 14 | RRMSE 0.042643 | Not available |

Task 11 and Task 12 dashboard headings appear swapped; comparisons were bound to our exact submission IDs and corresponding score emails, rather than page headings. Each task's `relative_score_record_current.json` preserves competitor rows, metric direction and evidence hashes.

The [all-task methods writeup](https://www.synapse.org/Synapse:syn78774887/wiki/644290) was published and submitted by the other owning chat. Its [receipt](../experiments/final_submission_2026_11_20/writeup_deadline/writeup_submission_receipt.json) confirms queue 9619664 and UI submission success at 22:00:42 UTC; no writeup score or email submission ID was yet available. Root's local alternative drafts were not published.

The full score/rank report was sent to Kevin's Telegram Saved Messages at 23:05 Europe/Berlin (22:05 UTC), with the sent message visibly present and composer cleared. The [delivery receipt](../experiments/final_submission_2026_11_20/telegram_all14_score_delivery.json) preserves the exact text and screenshot evidence. No further candidate or submission work remains for these 14 tasks.
