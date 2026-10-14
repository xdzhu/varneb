# E180 — complete-reference equality and shared-Q comparability limits

This zero-new-DFT experiment tests a necessary strong control before any
claim that mode reduction improves prediction. It uses the original two
centers, both raw Hessian amplitudes and all25registered Q coordinates per
center. The complete audited T/preserving sheets supply25/17seen candidates;
their other points stay unobserved, not last-trial labels. All800model-point
rows are retained; no coefficient, reference rotation, rank or best step is
selected by outcome. No physical+.005or independent material point is read.

## Two precisely different controls

1. Complete T optical rank33plus3open-cell directions spans the same measured
   36Dspace as the original two coupled directions plus their34Dcomplement.
   Impose exactly the same two Q constraints, retain the nonzero raw gradient,
   and compare direct KKT solution with the existing Schur solution.
2. For the smaller registered reference spans, allow their atomic directions
   plus3open-cell directions and **clamp all omitted directions**. Impose
   those same Q constraints, check constraint rank/free curvature, then check
   the unchanged.5A/minimum-distance geometry domain.

The second experiment is not E159's release of all outside-reference
directions; its previously negative outside-reference blocks are not made
positive or reclassified. Here the allowed null space is smaller because
omitted modes are clamped. It is a different, frozen-subspace control.

## Actual full-space result

All energies below are meV per12atomHf4O8cell, not per HfO2. Divide by4for
the latter. Shape discrepancies subtract each model's and DFT cohort's own
already seen origin; origin bias remains separately available in the data.

| Center | Step (A) | Seen non-origin DFT points | Full T max shape discrepancy | Adaptive full-space max shape discrepancy |
|---|---:|---:|---:|---:|
| T | .01 | 24 | 28.776664 | 28.776664 |
| T | .02 | 24 | 28.872042 | 28.872042 |
| Preserving | .01 | 16 | 20.521839 | 20.521839 |
| Preserving | .02 | 16 | 20.720640 | 20.720640 |

Across all25Qpoints for each center/amplitude, maximum full-T/KKT versus
adaptive/Schur energy difference is2.637e-16eV/cell, coordinate difference
1.559e-14A. This is numerical confirmation of the standard same-problem
identity, not independent predictive accuracy or physical resolution.
The finite-domain model discrepancies remain. The larger T maximum here
uses all25seen points, unlike E165's earlier six-point common snapshot;
the old table is not overwritten or counted as an independent dataset.

## Compatibility does not establish a winning small reference

Eligible non-origin **seen** point counts below are the same at both probe
steps; original denominator24remains in every statistic. Cmma frames are
kept individually even when equal. Nominal reference rank and allowed
dimension are different: final shared-observable count is2in every row.

| Atomic reference | Rank | Allowed dimension | T eligible seen points | Preserving eligible seen points |
|---|---:|---:|---:|---:|
| T geometric triplet | 3 | 6 | 8 | 0 |
| T complete doublets | 4 | 7 | 10 | 8 |
| All T optical | 33 | 36 | 24 | 16 |
| Cmma103 | 4 | 7 | 6 | 4 |
| Cmma106 | 4 | 7 | 6 | 4 |
| Cmma114 | 4 | 7 | 6 | 4 |
| Cmma127 | 4 | 7 | 6 | 4 |

There is **no non-origin coordinate jointly geometry-eligible for every
registered reference**, even though the linear constraints have formal
full row rank. Statistics on that common set are null, not zero error.
Outside-guard model values are retained only as diagnostics. Their enormous
harmonic extrapolations must not be counted as poor material predictions.
Nor do different own-eligible subsets support an accuracy ranking.

These particular Qaxes were chosen earlier from the coupled raw negative
space. A small reference may weakly overlap them and require a large
amplitude to meet them. This tests compatibility with **that operational
chart**, not all possible coordinates/models based on Cmma or T. A strong
reference in its own appropriate chart remains an uncompleted comparison;
the present limitation cannot replace it or establish adaptive superiority.

## Implementation and source boundaries

`vcneb.affine_quadratic` solves the standard equality-constrained quadratic,
with explicit raw-positive constraint-null check and no clipping/pseudoinverse.
It is backend-independent, not a new optimizer or default. See the
[API and mathematical boundary](../../../../docs/AFFINE_QUADRATIC_CONTROLS.md).
Allowed spans are transported with the unchanged registered mapping/bank;
both electronic settings and authors' geometry/eigenvectors are not retuned.
Source LDA frequencies do not supply our PBE curvature. Private author files
are hashed before and after reading, never written; no raw or transformed
basis/response arrays or model geometries are exported.

Original296curvature SCFs,118/172conditional-sheet SCFs and Gamma/canary/G2/
reference preparation remain additional. Every model sees the same measured
matrix here; this is not equal preparation cost, a speedup or an independent
accuracy result. Restricted and complete allowed spaces are explicitly unequal.

Verification:81local/81independent related checks across6modules,29new,
0failures/errors/skips; not full regression. The independent private-input
CLI reproduces the entire scalar report byte-for-byte39e31d79. Tested tree
40bf6241/archiveae30bf74. The same authorized private files are read-only
inputs in that replay, not included in its Git archive. Raw source/native
pins remain exact; this is code/artifact reproducibility, not the protocol's
independent material test. Earlier exploratory reports with development
metadata fields remain in Cscratch; no prospective forecast is overwritten.

Reproduce with authorized private sources and a fresh output path:

```text
python benchmarks/hfo2_channels/20261008/G3_shared_Q_controls_E180_20261015/analyze.py --source-root PRIVATE_SOURCE --output NEW_REPORT.json
```

Next actual evidence is atomic-release results and nonlinear, stable-branch
tests followed by the registered independent forecasts. The full-space
linear-algebra identity itself cannot be the JCTC novelty claim. This case
does not complete B0–B5/T/Cmma prospective comparisons or certify the surfaces.
