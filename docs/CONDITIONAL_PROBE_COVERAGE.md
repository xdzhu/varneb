# Finite independent checks of a conditional mode surface

`register_conditional_probe_plan` and `audit_conditional_probe_coverage` report
the finite test budget honestly. They supplement, not replace or weaken,
`screen_conditional_interpolation`. No calculator or scheduler is called.
The low-level predictor tested here is bilinear interpolation of the four
measured corners, not B2--B5 fitting or a new scientific prediction model.

## Registration before results

Use the same `ModePlane`, `ConditionalModeSurface` and per-node
`ConditionalPointEvidence` as the existing evidence gate. Register named,
unique, strictly interior test locations. Include one unique evaluation ID
per training node, explicit common energy zero, start/branch policy,
selection rule, protocol and all four screening thresholds.

```python
from vcneb import (
    ConditionalProbeSpec, ConditionalProbePlan, ConditionalProbeRecord,
    register_conditional_probe_plan, audit_conditional_probe_coverage,
)

plan = register_conditional_probe_plan(
    plane, surface, training_point_evidence, registered_locations,
    training_evaluation_ids=training_evaluation_ids,
    protocol_id=protocol_id, start_policy_id=start_policy_id,
    selection_rule_id=selection_rule_id, probe_labels_unread=True,
    energy_zero_id=energy_zero_id,
    gradient_tolerance=gradient_tolerance,
    maximum_neighbor_orthogonal_jump=neighbor_jump_limit,
    maximum_holdout_energy_error=energy_error_limit,
    maximum_holdout_coordinate_error=metric_coordinate_error_limit,
)
manifest = plan.to_dict()
# Save in a fresh namespace and publish the manifest, code/input/geometry
# hashes, branch starts and selection rule BEFORE measuring the test points.
# Later load that exact, published manifest:
plan = ConditionalProbePlan.from_dict(published_manifest)

report = audit_conditional_probe_coverage(
    plane, surface, training_point_evidence, plan, observations_by_probe_id,
)
```

The manifest binds reference, basis, metric, signed axes, pressure/boundary
labels, complete training energies/coordinates/diagnostics/evidence, probe
positions, selection policy and thresholds. A changed training snapshot is
rejected; an altered serialized manifest fails its checksum. A checksum is
**not proof of a publication time or blinding**. Verify it against the
previously published registration, record code and physical input hashes
upstream, and retain real source audits. IDs catch declared calculation
reuse but cannot authenticate independently computed raw results. A shared
batch-audit hash may legitimately describe several distinct calculations.

For newly sampled grids, the snapshot also binds the complete per-node
`point_records`, including unselected/unconverged starts. Changing or removing
those records after registration fails the same training-snapshot guard.
Internal four-corner patch screens retain their corresponding point records.
Previously registered grids lacking these optional records keep their exact
legacy snapshot semantics/hash; no old branch provenance is fabricated.
Neither diagnostic retention nor serialization grants a curvature/TS or
whole-grid interpolation certificate.

## Measurements, failures and missing points

Each observation is `ConditionalProbeRecord(evaluation_id, point=holdout)`.
The existing `ConditionalHoldout` holds Q, energy, ordered coordinates,
orthogonal gradient and the raw/branch/curvature/resolution evidence.
Unknown IDs, moved locations or a repeated training/test evaluation ID are
rejected, not silently dropped. Keep unsuccessful calculations as
`ConditionalProbeRecord(evaluation_id, failure="actual failure reason")`;
missing observations remain pending. Neither is replaced by zero energy.

Every measured point gets the unchanged evidence screen on its four-corner
patch: raw inputs/common energy zero, resolved positive orthogonal curvature,
orthogonal-gradient threshold, neighbor jump, Q drift, energy discrepancy
and mass/declared-metric coordinate discrepancy. All results, including
failed checks, remain in the report and in maximum observed errors.

Energy discrepancy is **measured minus predicted**, in `surface.energy_unit`.
Coordinate discrepancy is the norm in the declared metric, in
`plane.amplitude_unit`; it is not necessarily Cartesian Angstrom. Preserve
the cell/f.u. normalization and E versus H convention at the evaluator
boundary. No unit conversion, variance, Gaussian confidence interval or
continuous-domain error bound is inferred from these deterministic errors.

## What coverage means

- A 9x9 node grid has 64 cells. Eight passing tests in eight cells mean
  eight passing **point checks**, not 64 validated cells.
- Two tests inside one cell do not count as two covered cells.
- A passing point does not bound the error everywhere in its cell or prove
  branch uniqueness, a global lower envelope, a TS or a barrier.
- A failed point is not omitted from the maximum observed error; calculation
  failure and pending status are reported separately from numerical mismatch.
- The original full-grid screen is returned unchanged. It still requires a
  test in every cell and all its other gates. Even a full-grid exploratory
  pass is not a rigorous continuous-area certificate.

This reporting API does not authorize a global smooth contour from eight
tests, nor does it refit a failed model. Future localized plotting must state
its sampled/tested scope and independent evidence; absent a justified regional
error analysis, the code reports point checks only.

## Reproducible bounded verification

```text
python -m scripts.verify_conditional_probe_coverage --workdir NEW_DIRECTORY
python -m pytest -q tests/test_conditional_coverage.py tests/test_conditional_evidence.py tests/test_mode_surface.py tests/test_clamped_mode_evaluator.py
```

The synthetic verifier writes registration before constructing its eight
test labels; it tests 8-of-64 coverage, signed/metric discrepancies, one
failed calculation, one pending point and rejection of refits/reuse/unknown
observations. It makes **zero calculator calls**. These tests are software
validation, not an HfO2 conditional PES, material holdout or evidence of H2.

For the HfO2 study, finish all G2 paths before selecting G3 bottlenecks or
modes; each of at most two surfaces still has at most 81 modeling nodes and
eight independent probes. Conditional nodes may require multiple SCFs.
The reserved +0.005 path holdout remains unopened until the full registered
models/strong controls/error/cost package is frozen. No matrix, physical
parameter, NEB threshold, source deployment or allocation changes here.

For optional numerical curvature initialization during the underlying
relaxations, see [conditional optimization](CONDITIONAL_OPTIMIZATION.md).
It changes neither these independent-validation gates nor calculator budgets.
