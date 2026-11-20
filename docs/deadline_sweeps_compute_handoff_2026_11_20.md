# Tasks 5–8 compute handoff — 2026-11-20

At approximately 21:23 UTC, the user assigned compute/output completion for Tasks 5–8 to chat `01a1fad1-38c4-7b93-9e1d-a0ff8fd414f2`. The central submission chat owns external uploads. This agent has stopped new launches. **Do not duplicate the existing Task 6 and Task 7 jobs.** This agent will monitor those two jobs, preserve failures, collect their terminal evidence and leave their outputs at the paths below. Tasks 5 and 8 have not been launched.

The target candidate handoff is **22:30 UTC**, before the live 23:59 UTC deadline. These are explicitly user-authorized best-effort predictions, not scientific qualification. Earlier validation failures remain valid; no challenge answers or holdout truth may be used.

## Existing jobs and output locations

SSH alias: `betty`.

Base scientific output: `/vast/projects/smathi/cohort/korfmann/ghist2026/scientific/deadline_sweeps_20261120`

| Task | Slurm job | Current state at handoff | Output directory |
| --- | --- | --- | --- |
| 5 | none | Not launched; root downloading final input | No output directory created |
| 6 | `10859516` | COMPLETED32s, exit0; candidate verified | Base + `/task6` |
| 7 | `10859524` | FAILED45s, exit1; RAiSD loader missing GSL path | Base + `/task7` |
| 8 | none | Not launched; root downloading final input/exons | No output directory created |

Each existing job has 2 CPUs, 11 GiB RAM (`5632M` per CPU), and a 60-minute cap in `genoa-std-mem`. No automatic retries. The method lock permits at most 8 concurrent sweep CPUs. Do not launch duplicate jobs. Logs are `task6-10859516.log` and `task7-10859524.log` in the base directory.

On success, the only submission file is `<output>/deadline_submission.bed`. `<output>/completion.json` records the final intervals, strict BED validation, source/input/runtime hashes and geometry adjustments. `raw/prediction.json` preserves the unchanged method prediction. Task 7 also retains `raw/RAiSD_Report.task7`, `RAiSD_Info.task7`, `masked.vcf`, and its raw `candidate.bed`. A `failure.json` means failure; do not mistake an intermediate BED for completed inference.

## Frozen methods and software

Safe local prefix:

`/Users/kevinkorfmann/.cache/ghist-2026/complete-integration-20261120/experiments/submission_package_2026_11_20/sweeps_best_effort`

The same files are copied to the shared repository under `experiments/submission_package_2026_11_20/sweeps_best_effort`.

Remote frozen repository root:

`/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/sweep-deadline-inference-20261120`

Inside that root, the runner and lock are under the same `experiments/submission_package_2026_11_20/sweeps_best_effort` prefix.

| Artifact | SHA256 |
| --- | --- |
| `method_protocol.json` | `acc0ce4de87e4e99c282b97947dbe87b5748a2d8034dd519eb33410237a4fd17` |
| `run_final_sweep.py` | `76feb2e2e27bba5bfc60cc7d7cdf690e49c7c93bb38d0d3b108673a09ee09209` |
| `source_freeze.json` | `1fdcaf559e3f6233001438f5680414ba4bbfef29f6b5da99484872fa9701aeaf` |
| `source_snapshot_v2.tar.gz` | `27322a45d93f33b859c76a595b392ef81337880b2457a3f70c30cf78982431d4` |

The method lock predates all final input inspection. `operational_amendment.json` separately advances the handoff target to 22:30; the original method protocol remains immutable. Preparation revision 1 is retained; revision 2 only allows an absent Task 7 reference to trigger the already-locked Task 8 one-peak fallback. A supplied invalid reference still fails closed.

- Task 5: unchanged `blind_callable_median_three_v1`, half-width 60 kb.
- Task 6: existing development-selected `h12_snps50_maf0`, 50 SNP windows, MAF cutoff 0, half-width 75 kb. No independent confirmation claim.
- Task 7: unchanged masked `blind_callable_median_ten_v1`, RAiSD 50-SNP windows, half-width 80 kb.
- Task 8: explicitly provisional legacy RAiSD transfer. Fully callable windows, threshold twice the maximum callable Task 7 score more than 2 Mb from its own highest peak; descending score/lower-position ties; peaks at least 1 Mb apart. Snap only to callable exon bases within 60 kb. This transferred background is **not a calibrated or matched null**. If no reference background or no eligible candidate, the locked fallback is exactly one highest callable Task 8 peak, optionally snapped within 60 kb. It is not true-count inference.

For every task, the raw result is retained. Final BED geometry moves an inaccessible center to the nearest callable base (ties lower), clips its unchanged half-width to the containing callable component and contig bounds, and merges Task 8 overlaps. All changes are audited. Expected final chromosome is 9, with zero-based half-open intervals. The official mask starts at 40,910; do not assume the whole chromosome is accessible.

Twelve focused software tests passed locally and on Betty Slurm job `10859412` (8-second allocation, tests 1.23 seconds). They include actual small Task 6 execution, input mutation rejection, hash mismatches, geometry edges/gaps/ties, exon containment, Task 8 overlap/fallback and rejecting failed references. Two prior Linux pytest collection attempts failed because the ancestor old repository's pytest configuration inserted old source paths. Isolated `-c /dev/null --rootdir="$PWD" --confcutdir="$PWD"` fixed test discovery; production imports already resolved frozen source and the runner checks that explicitly. Raw check logs and failures are retained. No biological validation or accuracy claim is attached to these tests.

## Verified final inputs already staged

Root authenticated downloads have a small JSON provenance file beside each local input. These include entity/version, metadata MD5, verified MD5/SHA256, download route, and timestamp. The runner checks the supplied SHA256 again before and after inference.

| Input | Current entity/version | SHA256 |
| --- | --- | --- |
| Task 6 `softsweep.final.vcf.gz` | `syn75346851`, v2 | `5fc827ab97550ab91412bc51f5c62ebf4cd23f521f2fceab394ca402880bfdf0` |
| Task 7 `onesweep_bgs.final.vcf.gz` | `syn75346861`, v2 | `ef4f6fb9a795df476f6553e5890db1b3094fc026f288f11378b6c05c1ef67baa` |
| Final mask `9.accessible.bed` | `syn75821834`, v1 | `9f5235a687683eec6a7f19b0cf1334f545bda23497ff8e715658b1086778ec33` |

Local Task 6: `/private/tmp/ghist-final-submissions-20261120/task6/softsweep.final.vcf.gz` (5,968,095 bytes). Local Task 7: `/private/tmp/ghist-final-submissions-20261120/task7/onesweep_bgs.final.vcf.gz` (6,162,716 bytes). Local mask: `/private/tmp/ghist-final-submissions-20261120/taskshared/9.accessible.bed` (18 bytes).

Remote locations are `<base>/inputs/task6/softsweep.final.vcf.gz`, `<base>/inputs/task7/onesweep_bgs.final.vcf.gz`, and `<base>/inputs/shared/9.accessible.bed`, with provenance JSON alongside. Task 6 header was independently checked: exactly chr9, length 138,394,717, 25 diploid samples. The mask is `9\t40910\t138394717`.

Task 5 final entity is `syn75346852`, Task 8 is `syn75346876`, final exons are `syn77027379`. **Root is still downloading/verifying these; obtain actual paths/version/hash from root.** Do not substitute archived testing candidates, cached public answers, or old final versions. BED inputs must be plain text; preserve any compressed download and record decompression hashes if needed. VCF gzip is supported.

## Invocation contract for the new owner

Use the existing frozen `run.sbatch`. It loads Python 3.11.11, sets the frozen source/script paths and bounds numerical threading to 2 CPUs. All paths below are absolute; replace only task/input/output fields and authenticated download hashes. This is a template for unlaunched Tasks 5/8, **not permission to repeat existing Tasks 6/7**.

```bash
source /etc/profile >/dev/null 2>&1
sbatch --job-name=sweep-finalTASK \
  --output=/vast/projects/smathi/cohort/korfmann/ghist2026/scientific/deadline_sweeps_20261120/taskTASK-%j.log \
  /vast/projects/smathi/cohort/korfmann/ghist2026/reviews/sweep-deadline-inference-20261120/experiments/submission_package_2026_11_20/sweeps_best_effort/run.sbatch \
  --task TASK --vcf /absolute/verified/final.vcf.gz --vcf-sha256 VERIFIED_SHA256 \
  --accessible /vast/projects/smathi/cohort/korfmann/ghist2026/scientific/deadline_sweeps_20261120/inputs/shared/9.accessible.bed \
  --mask-sha256 9f5235a687683eec6a7f19b0cf1334f545bda23497ff8e715658b1086778ec33 \
  --out-dir /vast/projects/smathi/cohort/korfmann/ghist2026/scientific/deadline_sweeps_20261120/taskTASK
```

For Task 8 append `--exons /absolute/verified/exons.bed --exons-sha256 VERIFIED_SHA256 --task7-reference /vast/projects/smathi/cohort/korfmann/ghist2026/scientific/deadline_sweeps_20261120/task7`. Prefer the completed Task 7 reference. If unavailable before the deadline, omitting `--task7-reference` invokes the predeclared unsupported one-peak fallback and records why. Do not supply an incomplete/stale reference; it fails closed. Task 8 verifies the reference's source/mask/exposure and every bound artifact hash.

The runner refuses existing output directories. Preserve any failed attempt and coordinate a distinct directory for any justified operational rerun. Do not tune biological settings or thresholds from final outputs.

Existing Task 6 invocation used the template with `TASK=6`, its staged VCF and SHA256 above, and `<base>/task6`. Existing Task 7 invocation used `TASK=7`, its staged VCF and SHA256 above, and `<base>/task7`; no additional task-specific options were passed.

## Handoff files and verification

The preparation owned-file manifest is `/tmp/ghist-deadline-sweeps-owned.json`; it will be refreshed with terminal evidence for existing jobs. Cache is authoritative; exact owned files are copied to shared. Do not commit final VCFs or the large canonical VCFs. Their authenticated input provenance and content hashes suffice for the handoff; raw data remain local/Betty.

Before uploading, require Slurm terminal success, `completion.json` state `COMPLETED`, matching `deadline_submission.bed` hash, correct single-row Tasks 5/6/7 or nonoverlapping Task 8 rows, chr9 and within callable/contig bounds. Copy candidates and completion/provenance records to the central submission chat. This agent performs no external submission.

Task 6 terminal update: Slurm `10859516` completed in 32 seconds (runner26.36s), peak batch RSS416052K. All bound artifacts rehashed locally, BED revalidated within `[40910,138394717)`. Candidate copied to `/private/tmp/ghist-final-submissions-20261120/task6/sweeps_agent_result/deadline_submission.bed` and the safe/shared `sweeps_best_effort/results/task6/` directory. BED SHA256 `44836daae988246d2c42af96dc5402c66358ce5f743193901be5b99d19f780d1`. This is a provisional prediction, not a recovery result.

Task 7 terminal update: Slurm `10859524` failed after45seconds (batch RSS598808K) before RAiSD could produce a report. `raw/raisd.log` states `libgsl.so.28: cannot open shared object file`; `ldd` also showed missing `libgslcblas.so.0`. This was an environment omission in the new batch wrapper, not a data/model failure. The prior successful frozen Task7 smoke uses this line after loading Python:

```bash
export LD_LIBRARY_PATH="/vast/projects/smathi/cohort/korfmann/ghist2026/external/raisd/gsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
```

A read-only `ldd` check with that path resolves both GSL libraries to the existing installation. **The new owner must version/freeze this environment-only batch amendment before Task7/8 execution**, retaining the original job and failed output directory. No runner or biological method change is needed. Use a distinct Task7 output directory for a retry and point Task8's reference at the eventual completed directory. Do not infer success from the82MB canonical VCF left by the failed attempt. Raw canonical VCF remains on Betty; failure/started/RAiSD/Slurm records are copied to safe/shared `sweeps_best_effort/results/task7_failed/`. No retry was launched by this agent.

Both jobs launched by this agent are now terminal. There are no remaining active jobs under this agent's ownership. New owner owns all further Tasks5–8 computation. Task6 is complete and reusable; Task7 requires only the documented runtime fix and a distinct preserved retry.

Task5 recovery update at21:27UTC: new owner recovered the authenticated actual currentv2 `/private/tmp/ghist-final-submissions-20261120/task5/onesweep.final.vcf.gz`, MD5 `f8a29f8bb6b05273e81b290f6d0763cd`, SHA256 `1a11d96cc0a2742e6bb5428c2f72a21278bdb9b3eafc0e5612788eda9100eff0`. Use that original current artifact. An unused fallback had reconstructed the official POS+1 release note from authenticatedv1; independent paired audit confirmed629,667records, all6headers and all non-POS fields preserved, but ordinary-gzip output is incompatible with pysam. The actual currentv2 was subsequently compared to the reconstruction: every decompressed byte is identical, while actualv2 opens correctly in pysam. Both reports and original provenance are preserved under `sweeps_best_effort/input_reconstruction_audit/`; no reconstruction was used for inference and this reviewer launched no jobs. Actual-current comparison SHA256 `be631f22a4a15d0ae2f7189ec2d75b6a2523882e7fd5cbcf237ee6a1bc673f26`.
