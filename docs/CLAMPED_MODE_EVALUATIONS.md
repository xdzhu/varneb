# Mode evaluations under the same clamped-plane boundary as VARNEB

`ClampedModeCoordinates` bridges the existing exact active-cell chart to
`CalculatorModeEvaluator` and the ordinary fixed-Q orthogonal optimizer.
It is calculator independent: energy, forces and stress come from the same
ASE-compatible factory and raw-output validator, not a second backend API.
This implementation is not a material conditional-PES/TS/prediction certificate.

## Geometry and units

For row cells, with explicit scale L in Angstrom:

```text
H = H0 @ [I + sum_j (y_cell,j/L) D_j].T
s = s0 + u @ inv(H0)
```

The D matrices are exactly the one/three Frobenius-orthonormal cell directions
of the caller's `ClampedPlaneBoundary`, also used by endpoint and path
optimization. Two substrate row vectors remain fixed, not just their stresses.
Open tilts are generally nonsymmetric deformation, not six symmetric strains.
All atomic u and scaled cell y coordinates have units Angstrom; physical
gradient is eV/Angstrom. H=E+PV uses compressive-positive P in eV/Angstrom³,
raw ASE tensile-positive stress and actual volume. NEB/spring/CI forces are
never used. An externally prescribed substrate strain is not a relaxable
coordinate in this chart.

The historical `ReferenceCellCoordinates` convention (atomic Angstrom plus
six dimensionless symmetric strains) and its exact serialized evaluator
contract remain unchanged. Do not reuse its 42-coordinate Hf4O8 plane as a
39-coordinate clamped plane, or copy dimensionless-strain metric weights
without the required L conversion. For a mass-metric mode plane here, atomic
weights are masses in amu and scaled-cell weights must also be explicitly
declared in amu; Q then has sqrt(amu)*Angstrom units. There is no universal
relative atomic/cell metric or automatic phonon assignment.

## Use the declared common reference and periodic lift

```python
from vcneb import ClampedModeCoordinates, CalculatorModeEvaluator

# boundary is exactly the boundary already used by endpoints/VARNEB.
chart = ClampedModeCoordinates(reference_atoms, boundary=boundary,
                               cell_scale_A=registered_cell_scale)
coordinates = chart.from_atoms(ordered_unwrapped_path_image)
evaluator = CalculatorModeEvaluator(
    chart, existing_fixed_ASE_factory, workdir=fresh_directory,
    calculator_id=registered_physical_contract_id,
    result_validator=existing_raw_SCF_validator,
    validation_id=registered_validator_id,
    pressure_eV_per_A3=registered_pressure,
    minimum_distance_A=registered_geometry_limit,
)
potential_eV, physical_gradient_eV_A = evaluator(coordinates)
```

`from_atoms` requires an already registered atom order, frame, origin and
continuous periodic lift. It does no nearest-image wrapping, symmetry repair,
rotation fit or atom rematching. A full-cell displacement remains a full-cell
displacement. Incompatible substrate/species/constraints are rejected.
Reference accessors return copies; calculator data are not attached to chart
geometry copies. A common reference permits both escape candidates to be
represented without an image-specific gauge or reference fit.

Use a `ModePlane` of matching dimension, metric and reference. Source modes
must first be transported into that declared coordinate chart and normalized
explicitly; an arbitrary linear combination is not a new phonon eigenvector.
The existing `relax_orthogonal_at_q` accepts the evaluator, explicit rigid
translation gauge from `chart.rigid_translation_directions()`, and a trial
geometry validator. For safeguarded BFGS that validator is mandatory and
runs before any calculator call. Orthogonal-gradient convergence is not a
positive-curvature, competing-branch or global-minimum certificate; all
`branch_outcomes` must be retained.

## Recovery does not change the physical experiment

The evaluator records the substrate, tilt policy, deformation basis, L,
reference geometry, declared calculator and validator IDs, pressure and
geometry limit. Changed tilt/scale/reference/backend IDs reject resume.
Exact completed coordinates reuse cache; failed/partial directories remain
and retried evaluations use fresh directories. Upstream raw-input/log audits
must still prove the declared calculator ID: a string alone does not do that.
Only complete evaluations enter `n_new_evaluations`; failed attempts and
preparation costs must also be counted from their original provenance.

Native 6-component ASE stress is accepted; raw 3x3 stress must be symmetric.
It is checked before ASE's conversion can average a nonsymmetric tensor.
Geometry mutation or invalid E/F/stress never becomes a valid cached result.

## Validation and HfO2 study boundary

```text
python -m scripts.verify_clamped_mode_evaluator --workdir NEW_WORKDIR --output NEW_REPORT.json
python -m pytest -q tests/test_clamped_mode_evaluator.py tests/test_mode_evaluator.py tests/test_clamped_joint_curvature.py tests/test_mode_surface.py
```

The verifier performs 232 EMT energy/gradient work evaluations for eight
oblique/rotated, normal-only/open-tilt and P0/nonzero-P cases, then reparses
18 portable native terminal records from the existing T/M training chains.
It checks exact cache recovery, common-reference ordered roundtrips and
agreement with the unchanged active-gradient primitive. EMT finite differences
are independent implementation checks, not new native DFT finite differences
or uncertainty bounds. Native frames are replayed, not recomputed.

G2 still must finish before G3 bottleneck/mode selection. The original
100 Ry/full 10 au DZP/PBE/KPT/SCF/no-CI contract and held-out +0.005 boundary remain
unchanged. This bridge creates no new chain, surface, local Hessian or holdout.
The existing whole-grid interpolation screen still requires every cell's
independent validation; eight independent G3 samples cannot be misreported
as 64 validated cells. [Finite probe reporting](CONDITIONAL_PROBE_COVERAGE.md)
now records registered point errors, failed/pending tests and the limited
observed-cell coverage without weakening that gate. Continuous-region error
bounds and material G3 validation remain separate, uncompleted requirements.
