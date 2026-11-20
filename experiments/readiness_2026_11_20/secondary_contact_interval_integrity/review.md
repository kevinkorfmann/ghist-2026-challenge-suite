# Secondary-contact interval domain correction

The actual Task 2 freeze entrypoint accepted an empty interval object, reversed migration bounds `[2, 1]`, and negative population bounds `[-3, -1]` in synthetic reproductions against merged production source. These artifacts could be labeled method-frozen despite malformed sampling intervals.

The correction validates exactly the five parameter names emitted by the existing bootstrap producer, two finite ordered bounds per parameter, positive population sizes and split time, and nonnegative contact time and migration. It preserves the existing mainland/island and generations-since-migration output vocabulary. Zero contact-time and migration intervals remain valid. Validation precedes atomic writing, preserving an existing output on failure.

This is an interval shape/domain guard. It does not establish bootstrap denominator provenance, empirical coverage, generating-model adequacy, selected-fit endpoint provenance, or scientific readiness. The historical Task 2 recovery and misspecification failures remain unresolved. No model, optimizer, acceptance threshold, data, or simulation changed.

The fixed source passed 102 focused tests and Ruff. Fifty-four new actual-entrypoint tests cover malformed bounds, complete parameter sets, valid boundaries, and unchanged existing outputs. The previously valid complete interval reproduction remains valid with the same values; the three malformed examples now fail. Broader integration and independent review are recorded separately.
