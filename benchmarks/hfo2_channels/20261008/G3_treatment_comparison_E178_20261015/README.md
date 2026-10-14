# E178 — signed three-treatment comparison, including incomplete observations

This is a postprocessor for the two existing centers, both original Hessian
probe steps and five already registered Q coordinates per center. All20
model-point rows are retained. It reads the byte-bound E176 pre-label
expectations, E173 frozen/joint observations and an explicitly hash-bound
E176 native atomic-audit report. It never launches a calculator, contacts a
scheduler, generates a geometry or reads the physical+.005holdout.

## What is compared

F: frozen transverse coordinates; A: release only the31atomic fixed-Q
directions while freezing that Q seed's cell; J: joint atomic/open-cell
release. Q itself contains cell motion, so neither F nor A defines one
common fixed cell across all Q points. The signed A-F, J-A and J-F values
are candidate-treatment differences, not unique mode energies or proven
differences between global nested minima.

For each pair, retain its absolute difference and its change relative to
the same pair's observed origins. Keep the raw origin offsets separately.
The pointwise identity is checked both for treatment differences and for
model-minus-DFT error differences. Do not subtract maxima from different
points to claim a decomposition. Positive release is retained, not clipped.

Missing A does not erase known J-F. Missing only the A origin suppresses
A-related shape comparisons, not their available absolute differences or
J-F shape. A processed hole is distinct from an unprocessed point. Every
summary lists its actual IDs and denominator4, excluding the origin;
empty observations have null error, not a perfect zero. No ranking of
incomplete treatment cohorts is performed.

## Actual evidence available for this report

Atomic input is the **empty audited prefix at06:11:47CST**, SHA b7027bd8;
there is no measured A candidate in this version. The live job remained
pending when checked at06:30:19CST. The report therefore contains available
F/J comparisons and explicit A abstentions, not a new atomic material result.
E173's previously measured four axial observations per center are reused:

| Center | Probe step (A) | Max frozen shape discrepancy | Max joint shape discrepancy | Max missed joint-release shape change | Measured A axial points |
|---|---:|---:|---:|---:|---:|
| T | 0.01 | 0.416264 | 1.327177 | 1.168761 | 0/4 |
| T | 0.02 | 0.403160 | 1.346071 | 1.174551 | 0/4 |
| Preserving | 0.01 | 0.075880 | 4.671404 | 4.606959 | 0/4 |
| Preserving | 0.02 | 0.116535 | 4.720758 | 4.606393 | 0/4 |

All discrepancies above are meV per12atomHf4O8cell, **not per HfO2**;
divide by4for the latter. They are training discrepancies, not independent
prediction errors or measured physical uncertainty. These numbers reproduce
the earlier two-treatment evidence, not another four independent datasets.
This localizes the preserving discrepancy to a release-response term but
does not yet isolate atomic versus cell release, anharmonicity, branch
changes or finite residuals as its microscopic cause.

The report separates8old frozen statics/8.136814completed-call proxy
core-hours,0observed new atomic SCFs and296additional original curvature
SCFs. Cached preparation and original joint runs are explicitly additional
and not recounted here. It is not a full campaign cost or speedup estimate.
Atomic costs cover only completed nodes in the supplied audit; in-flight
or unaudited calculations are not silently counted as free. The source
accounting state and audit timestamp/mode remain explicit in the output.

## Reproduce or consume a later native audit

From the repository root:

```powershell
python benchmarks/hfo2_channels/20261008/G3_treatment_comparison_E178_20261015/analyze.py --atomic-audit benchmarks/hfo2_channels/20261008/G3_atomic_audit_E176_20261015/native_prefix.json --atomic-audit-sha256 b7027bd8db935fbc9e87f8b56a92999dbedeab50454e75538c2e5442c3eccd45 --output NEW_COMPARISON.json
```

For later observations, first run the unchanged E176 raw-native audit on
the original job, then supply its exact fresh report path/hash. The
postprocessor checks handle, accounting identity, ordered Q prefix,
budget/candidate counts, replay flag and atomic norm. It does **not**
reparse raw SCFs itself and does not replace that audit. A large joint
gradient does not invalidate an otherwise valid atomic-only candidate;
it prevents calling it a joint stationary point. A complete audit requires
normal job completion and all10processed nodes, not10successful candidates.

## Regression failure retained and corrected

First related run103pass/1fail, repeated with the previous startup helper;
the previous75test set alone also gave74pass/1fail. The failing old test
required full bitwise equality of recomputed floating model values with
the sealed E176 JSON. Diagnostic comparison finds124floating differences,
maximum2.842170943040401e-14; all nonfloating metadata and source pins agree.
The exact low-level cause is not proved by these observations.

Keep the original prediction bytes unchanged and assert their SHA256
explicitly. Only formal-model arithmetic replay now uses1e-12absolute/
zero-relative tolerance, with exact types, keys, IDs and metadata. Negative
tests reject a1e-9numerical change and a source-label change. This is a
unit-test reproducibility criterion, **not a DFT tolerance or physical
error bound**. Native exact cache/optimizer checks, six physical assets,
100Ry/full10auDZP and.03atomic/.10ordinary-NEB gates are not changed.
Both failed related-test XML files and the numeric diagnostic are retained.

Current verification scope is106related checks over5modules, including29
new comparison checks and2new negative formal-replay checks, not a full
repository regression. Independent archive replay and exact output-byte
comparison are recorded in validation.json:106local/106independent pass,
zero failures/errors/skips. Independent tested source tree5b710adf,
archive d1d6ac49; complete report bytes match exactly e467df51. This does
not upgrade the empty atomic prefix into a measured material result.

Next: consume actual B3native candidates/holes and compare the signed
A-F/J-A factors with both unchanged formal models and actual costs. Full
T/Cmma/B0–B5 controls, released-branch stability, final training-only
representatives/domain and independent predictions remain separate work.
