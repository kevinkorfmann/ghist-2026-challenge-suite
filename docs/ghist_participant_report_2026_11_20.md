# GHIST 2026 — Kevin Korfmann’s results and participant ranking

**Kevin Korfmann (`ancestrum`) has valid scored submissions to all 14 challenges.** The 12 public final leaderboards were refreshed on November 20, 2026, between **23:52:36 and 23:53:23 Europe/Berlin** (22:52:36–22:53:23 UTC). This is an unofficial snapshot before the deadline; it is not a closed final result.

**By number of challenges won, paiksca, wisezeal, and austin.t.daigle are tied with four wins each across the 12 visible boards.** Kevin is second of four in Challenge 9 on both published metrics. His entries in Challenges 13 and 14 are scored, but comparable standings are not available in this record.

[Detailed Challenge 9 method, validation and submitted result](task9_method_and_results_2026_11_20.md) · [Machine-readable evidence and calculations](ghist_report_evidence_2026_11_20.json)

## Ranking by challenges won

| Wins rank | Participant | Wins | Winning tasks | Entries on 12 visible boards | Total rank points | Mean rank points per entry |
|---:|---|---:|---|---:|---:|---:|
| 1 | paiksca | 4 | 4, 10, 11, 12 | 12 | 8.150 | 67.9% |
| 1 | wisezeal | 4 | 1, 2, 3, 5 | 8 | 5.333 | 66.7% |
| 1 | austin.t.daigle | 4 | 6, 7, 8, 9 | 5 | 4.417 | 88.3% |
| 4 | rgollnisch | 0 | — | 7 | 2.950 | 42.1% |
| 4 | srong | 0 | — | 4 | 2.083 | 52.1% |
| 4 | ancestrum | 0 | — | 12 | 2.000 | 16.7% |
| 4 | tstentella | 0 | — | 1 | 0.800 | 80.0% |
| 4 | lzong | 0 | — | 4 | 0.667 | 16.7% |
| 4 | candicewu | 0 | — | 1 | 0.600 | 60.0% |

Wins determine the primary rank and ties are preserved. For an optional overall ordering, rows break equal win counts by total rank points. On a task with n > 1 entrants, rank points are `(n − average rank)/(n − 1)`, from zero to one. Tied scores share the average occupied rank. Unentered tasks add zero to total points. The mean uses entered tasks only, so it measures performance conditional on participation, not breadth.

A sole valid entrant earns **one full win and one point**. No minimum competitor count is imposed. If Kevin is confirmed as the sole valid entrant on both Tasks 13 and 14, those add **two wins**, placing him fourth by wins relative to these visible-board results. This remains conditional because competitor counts are unknown. A complete 14-task ranking also requires verifying everyone else’s results there.

Lower RRMSE is better; higher sweep F1 and relatedness weighted Jaccard are better. Relation and meiosis ranks are averaged within each relatedness task, so a task cannot count twice. Both metrics give the same participant ordering here. Rankings use the displayed precision for everyone, rather than comparing Kevin’s private full precision with competitors’ rounded values. Raw scores from different task families are never averaged.

## Kevin’s fourteen scored entries

These exact scores come from retained scored-submission receipts. They establish accepted scoring, not scientific qualification of every method. In particular, the low-coverage methods retained unmet scientific qualification gates.

| Task | Submission | Exact server scores | Current visible standing |
|---:|---|---|---|
| 1 | `9785374` | RRMSE: 0.7330006155047631 | 6 / 6 |
| 2 | `9785375` | RRMSE: 2.1376729777837356 | 3 / 4 |
| 3 | `9785376` | RRMSE: 0.796883521581379 | 3 / 3 |
| 4 | `9785377` | RRMSE: 0.5169210928553797 | 3 / 3 |
| 5 | `9785369` | f1: 0.0; recall: 0.0; size: 1.2 | 4.5 / 7 (tied midrank) |
| 6 | `9785370` | f1: 0.0; recall: 0.0; size: 1.5 | 6.5 / 7 (tied midrank) |
| 7 | `9785371` | f1: 0.5555555555555556; recall: 1.0; size: 1.6 | 5 / 7 |
| 8 | `9785372` | f1: 0.10695187165775401; recall: 0.06060606060606061; size: 2.4 | 6 / 7 |
| 9 | `9785353` | meioses_count: 0.8423845484005662; relation: 0.8337060077932977 | 2 / 4 |
| 10 | `9785368` | meioses_count: 0.1283610601352386; relation: 0.12225889450060641 | 2 / 2 |
| 11 | `9785361` | meioses_count: 0.5724422540498211; relation: 0.42540795348737687 | 2 / 2 |
| 12 | `9785373` | meioses_count: 0.13564060673696537; relation: 0.11235515689383353 | 2 / 2 |
| 13 | `9785362` | RRMSE: 0.6172118592711575 | Unverified; scored entry confirmed |
| 14 | `9785360` | RRMSE: 0.042643298573969445 | Unverified; scored entry confirmed |

## Challenge 9 contribution and result

Kevin directed a simulation-trained random-forest entry that converted exact IBD segments into **218 fixed features** and predicted ten relationship labels plus their deterministic meiosis counts. It trained on **60,266 pairs**, was independently confirmed on **20,023 observed pairs**, and produced all **997 final pair predictions** from **39,916 input segments**. The detailed report explains the accepted-label support repair, the rejected alternative model, validation design, runtime checks, source hashes, and known weaknesses for distant cousins. AI agents assisted with analysis, implementation, execution, review and documentation under Kevin’s direction.

The exact final scores were **0.8337060077932977 for relationship** and **0.8423845484005662 for meioses**. Both are weighted Jaccard scores, not ordinary accuracy. The earlier three-entry lead was superseded by austin.t.daigle’s entry. The refreshed four-entry order is austin.t.daigle (0.838 / 0.849), ancestrum (0.834 / 0.842), paiksca (0.828 / 0.840), and rgollnisch (0.822 / 0.827), listed as relationship / meioses.

## Source tables and scope

| Task | Leader | Official table | Checked at, Berlin |
|---:|---|---|---|
| 1 | wisezeal | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640737) | 23:52:43 |
| 2 | wisezeal | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640747) | 23:52:58 |
| 3 | wisezeal | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640742) | 23:53:21 |
| 4 | paiksca | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640750) | 23:53:21 |
| 5 | wisezeal | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640748) | 23:53:21 |
| 6 | austin.t.daigle | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640744) | 23:53:22 |
| 7 | austin.t.daigle | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640749) | 23:53:22 |
| 8 | austin.t.daigle | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640745) | 23:53:22 |
| 9 | austin.t.daigle | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640981) | 23:52:36 |
| 10 | paiksca | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640982) | 23:53:22 |
| 11 | paiksca | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640984) | 23:53:23 |
| 12 | paiksca | [Leaderboard](https://www.synapse.org/Synapse:syn74803798/wiki/640983) | 23:53:23 |

The public headings for Tasks 11 and 12 are swapped; this report maps the rows using submission IDs and queue-specific receipts. The demographic dashboard lists no public standings for Tasks 13 and 14 in the retained check; earlier anonymous evaluation requests returned HTTP 403. Their absence from the wins table is an evidence gap, not a discount for uncontested wins.

The evidence file retains all refreshed rendered tables, parsed entries, ranks and point totals, exact score extracts for all fourteen submissions, and the earlier detailed Task 9 audit. Historical observations remain timestamped so that subsequent entries do not silently rewrite what was known earlier.
