# Task 14 pragmatic inference handoff, frozen before outcomes

This is a separate practical submission path authorized by the coordinating root after the user requested more submissions if somewhat OK. It does not qualify the failed original six, cap200, or unlaunched strict cap500 policies. No final input has been accessed by the lowcoverage agent.

The frozen source is `/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/source_v12`; all 25 preceding source files are unchanged. New development array `10858309` runs both old regime0 Q30/Q20 cohorts with full window500, tolerance1e-6, cap500, and the same three estimator runs and eight starts for each of two demographic optimizers. Seed1 is always the primary estimate. The two workers actually passed their start guards at20:30:09/20:30:10 UTC. No fresh963xxx seeds are consumed.

Bindings:

- Freeze `task14_pragmatic_freeze_v1.json`: `7a98a0e446cf5b5221fef48bace8743b5a604f7cd807fdff0c26df3df001b1bb`.
- Protocol `task14_pragmatic_protocol_v1.json`: `df5b72848dda79fbf95b474abe83f6240f9689903c0730e666bdab951956dc57`.
- Source manifest: `15725e52355e9a7138e3f345d5838518a81ca3ecb6d335accf66b4fa1286be66`.
- Final wrapper `source_v12/scripts/lowcoverage_predict_pragmatic.py`: `ef415c585e151b457b929a39fb6037165bbde92e92f86e4ecce06be3a9ba2f91`.
- Separate prospective practical addendum `task14_pragmatic_usefulness_addendum_v1.json`: `6318852e95ddc59004470a5e4971e3625fcb68bf60975b111a1191b272204f13`.

Root must first independently audit both completed development cases and `task14_pragmatic_development_gate_audit_v1.json`. The practical rule requires complete valid estimates, successful interior optimization and optimizer agreement in all three runs, exposure>=99%, all development target errors<=30%, and every seed2/continued target within5% of fixed seed1. The original strict1% target gate, spectrum gates, and observed stopping remain unchanged and reported. A practical pass is not a scientific qualification.

Only root can create the explicit decision JSON after its audit. Required fields are `authorize_final_inference: true`, `strict_qualification: false`, the exact `freeze_sha256`, `wrapper_sha256`, `development_report_sha256`, and `primary_estimate: seed1`. Root should additionally bind its independent audit and practical-addendum hash and explain any strict failures accepted. The wrapper hashes the decision before and after execution. This document is not that decision.

Root-owned Slurm execution must use4CPUs,5632M/CPU,55minutes, source `/etc/profile` before strict shell options, `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, and `PYTHONPATH` selecting `source_v12/src:source_v12`. The wrapper enforces an actual final start no later than22:00UTC, leaving at least64minutes after the55minute wall cap before23:59UTC.

```text
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python \
  /vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/source_v12/scripts/lowcoverage_predict_pragmatic.py \
  --frozen-root /vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120 \
  --freeze-sha 7a98a0e446cf5b5221fef48bace8743b5a604f7cd807fdff0c26df3df001b1bb \
  --root-decision ROOT_OWNED_DECISION_JSON --decision-sha VERIFIED_DECISION_SHA256 \
  --directory ROOT_VERIFIED_BAM_AND_POPFILE_DIRECTORY \
  --reference ROOT_VERIFIED_CURRENT_REFERENCE_FASTA --output NEW_ROOT_OWNED_OUTPUT_DIRECTORY
```

The output includes all numerical gates, the fixed primary estimate, pre/post source/runtime/input/SAF bindings, and `task14.pragmatic.submission.draft.yaml`. A draft is written after valid serialization even when a strict numerical gate fails; it remains explicitly unsubmitted and scientifically unqualified. Root must inspect the actual final all-three-run sensitivities and provenance before any submission. The wrapper has no downloader, challenge truth access, or submission transport.

The historical window100 estimator remains the broader six-cohort evidence base: primary18-target median absolute error5.354%, worst all-run error11.331%, but every cohort fails strict1% target stability (1.08-3.93%). The cap200 full-window development reduced target sensitivity below0.04% in two regime0 cases but failed stopping in five of six estimator runs. No bias correction or seed selection is authorized by these observations.
