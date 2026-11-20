# GHIST login and submission check — November 20, 2026

The user selected **individual** participation for the intended submission.
No prediction was uploaded, no final inference was started, and no competition
submission was made during this access check.

Browser sign-in succeeded. After the user's explicit confirmation of the
certification agreement and challenge rules, the certification quiz passed
**15/15** and the 2026 home page confirmed successful registration on November
20. A project was created with **Private** visibility selected:
[GHIST 2026 - Kevin Korfmann](https://www.synapse.org/Synapse:syn78774887/wiki/).

The authenticated Challenge 9 final-submission dialog opened successfully and
identified its destination as **GHIST 2026 Relatedness 1 Final**. It was closed
without selecting a file or submitting. This verifies access to the submission
flow; it does not establish a remaining-quota count. No final attempt was
consumed by this work. Individual participation is the user's choice, but no
submission has yet committed that choice on Synapse.

## Deadline and submission mechanics

The live [GHIST home page](https://www.synapse.org/Synapse:syn74803798/wiki/640571)
and [file tutorial](https://www.synapse.org/Synapse:syn74803798/wiki/640752)
give a deadline of **November 20, 2026, 23:59 UTC**, which is **November 21,
2026, 00:59 Europe/Berlin**. Testing attempts are unlimited; each final
challenge permits one submission. The first submission fixes individual/team
participation, including when submitting a testing attempt.

The browser workflow is to upload a prediction file to a personal Synapse
project, select that file in the relevant challenge's submission dialog, and
check its receipt and submission status. Uploading a file alone does not submit
it. The existing Python downloader requires `SYNAPSE_AUTH_TOKEN`; browser login
does not set that environment variable. No token or `~/.synapseConfig` was
available to this shell during the check. Manual browser submission does not
require configuring a Python token.

Synapse's [participation guide](https://www.synapse.org/Synapse:syn74803798/wiki/640743)
says upload/submission requires Certified User status, challenge registration,
and any applicable data-use conditions. Certification, registration and creation
of the personal project are now verified. Prior submissions, remaining quota
and file-specific conditions have not been separately audited.

## Registration link discrepancy

Use **Click Here to Register on the GHIST 2026 home page** if registration is
needed. Its group is `3588302`, named `GHIST 2026 Participants`. The project
challenge API confirms challenge `4848`, participant group `3588302`.
The [Participants & Teams page](https://www.synapse.org/Synapse:syn74803798/wiki/640746)
agrees with those IDs.

The current How to Participate page still embeds registration group `3537096`
(`GHIST 2025 Participants`) and submission-team challenge `4788`. Those are
stale links. Do not use that page's registration widgets for the 2026 event.

## Most promising submission candidate

Challenge 9 (IBD pairs) is the first candidate to assess. Its ten-label refit
has untouched confirmation accuracy of 89.10% under matched interference and
87.20% under the Poisson crossover shift, without undefined accepted-label
probabilities. These are simulation results conditional on observed IBD.
The candidate has not been promoted or used for challenge predictions.
See the [scientific relatedness report](scientific_relatedness_2026_11_20.md).

The live [relatedness specification](https://www.synapse.org/Synapse:syn74803798/wiki/640950)
describes true IBD segments and the same ten accepted labels as the repaired
model. Current final metadata for `syn75202631` reports version **3**, file
`relatedness1.final.csv`, file handle `177273793`, modified September 14, 2026.
This metadata check did not download or inspect the final CSV.

Before a final file exists: approve/promote the corrected model based on its
recorded evidence, pin and verify the current input, check its domain and pair
coverage, run inference, validate the output schema and all input pairs, and
review diagnostics and the exact file to be submitted.

The [latest scientific review](scientific_validation_2026_11_20.md) does not
certify any complete pipeline ready for final submission. Saved testing files
under `artifacts/review_2026_11_20/testing/` are not final predictions. Tasks 1
and 5 have encouraging historical evidence but material outstanding release
checks; they must not be promoted on the strength of their old freezes alone.

## Queue IDs advertised by the live challenge pages

These IDs were read from the current public wiki submission widgets. Direct
queue access initially returned HTTP 403 anonymously. After registration, the
authenticated Task 9 final dialog was verified as described above. The other
queues and per-account quota counts have not been verified.

| Task | Testing queue | Final queue |
|---|---|---|
| 1 — Bottleneck | 9619647 | 9619646 |
| 2 — Secondary contact | 9619649 | 9619648 |
| 3 — Growth with BGS | 9619651 | 9619650 |
| 4 — Split/migration with BGS | 9619653 | 9619652 |
| 5 — Hard sweep | 9619657 | 9619656 |
| 6 — Soft sweep | 9619679 | 9619680 |
| 7 — Hard sweep with BGS | 9619659 | 9619658 |
| 8 — Multiple sweeps with BGS | 9619663 | 9619662 |
| 9 — IBD pairs | 9619672 | 9619673 |
| 10 — Large human family | 9619667 | 9619674 |
| 11 — Inbred multifamily cohort | 9619675 | 9619676 |
| 12 — Family with background relatedness | 9619677 | 9619678 |

Sources: [demography](https://www.synapse.org/Synapse:syn74803798/wiki/640740),
[sweeps](https://www.synapse.org/Synapse:syn74803798/wiki/640755),
[relatedness](https://www.synapse.org/Synapse:syn74803798/wiki/640950).

The [writeup tutorial](https://www.synapse.org/Synapse:syn74803798/wiki/640753)
advertises queue `9619664`. Consortium-level authorship eligibility requires a
public methods writeup submitted as the same individual or team as the final
prediction. Its queue also needs authenticated verification. Prepare the
methods, limitations and source-code references before publishing the project.

## Optional command-line access

For automated current-data downloads, Synapse supports a personal access token
with `view` and `download` scopes. Upload/submission via a client additionally
requires the relevant write permissions and token scope. Configure secrets
privately; do not put them into this repository or the chat. See the official
[client configuration](https://docs.synapse.org/synapse-docs/client-configuration)
and [account management](https://docs.synapse.org/synapse-docs/managing-your-account)
documentation. No new token was created during this check.
