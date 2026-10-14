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

The safeguarded check also applies if the two Q coordinates and explicit
frozen directions span the full chart. There is then no orthogonal
optimization: an invalid geometry is rejected before any energy call, and
an accepted point is evaluated exactly once. An empty projected gradient
is not evidence that an unmeasured physical direction is stable.

The whole-grid `sample_conditional_mode_surface` now forwards the same
`optimizer`, `trial_validator`, `initial_trust_radius` and
`include_frozen_start` options. Thus an explicitly guarded clamped calculation
does not lose its pre-evaluation geometry screen when moved from one point
to a grid. The default remains L-BFGS-B for compatibility; use safeguarded
BFGS with the registered geometry validator for geometry-limited DFT trials.
This does not alter the calculator's electronic inputs or select material
modes, branches or tolerances.

```python
surface = sample_conditional_mode_surface(
    plane, registered_q1, registered_q2, evaluator,
    energy_unit="eV", potential_kind="E", boundary_condition=boundary_id,
    frozen_directions=chart.rigid_translation_directions(),
    starts=registered_starts, include_frozen_start=True,
    optimizer="safeguarded_bfgs", trial_validator=registered_geometry_validator,
    initial_trust_radius=registered_metric_radius,
)
```

Grid `n_evaluations` now includes callbacks from **every** supplied start,
including unconverged candidates, not just the selected candidate. The
single-point legacy `n_evaluations` still describes its selected start;
`total_n_evaluations` sums complete branch outcomes. Older multi-start records
without those outcomes have an unknown total, not a selected-only estimate.
Callback counts include exact-cache hits and final verification calls; they
are not unique new SCFs or core-hours. Audit the evaluator's source IDs,
failed attempts and actual preparation resources separately.

The grid also retains `surface.point_records[j_q2][i_q1]`: complete defensive,
read-only copies of each `ConditionalModePoint` and its `branch_outcomes`,
including unselected or gradient-unconverged starts. The legacy selected-start
count and total grid count remain distinct. No extra evaluator calls are made
to retain this evidence. Manually imported older grids default to an empty
record tuple, meaning unavailable provenance, not zero attempted starts.
Every recorded node must agree with its grid Q, selected coordinates/energy/
residual/start and all-start count; incomplete or inconsistent records fail.
A node where no start converges still raises the original error instead of
inventing a point. Raw backend failures and preparation costs must be kept in
the calculator's separate audit ledger. These records are not curvature,
global-minimum or branch-continuity certificates.

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

## Optional cumulative evaluation-attempt budget

`CalculatorModeEvaluator(..., max_evaluation_attempts=registered_limit)`
limits fresh reserved evaluation attempts across the entire workdir, including
failed/partial output and process restarts. The default `None` retains the
legacy uncapped flow and its exact physical cache contract. A capped namespace
stores a separate immutable operational policy; no electronic setting,
geometry chart or physical cache identity changes. Choose/publish the limit
with the experimental plan, not after viewing held-out labels.

Geometry validation and exact completed cache reuse precede admission. A
fresh attempt reserves an exclusive indexed marker before the factory is
invoked. Failure/interruption does not refund the slot, even if the factory
never starts. Exhaustion raises before any new calculator call. Resume must
keep the same policy and complete reservation history; missing, corrupted,
inconsistent or over-budget records are rejected, not regenerated. The cap
cannot be retroactively attached to an old uncapped workdir, removed or raised
through ordinary resume. Any reviewed continuation requires a separately
audited operational decision; this API does not make it automatically.

`n_evaluation_attempts` reports cumulative conservative reservations (legacy
uncapped workdirs count preserved eval directories). `n_new_evaluations` keeps
its old meaning of completed new evaluations in the current process. Cache
hits consume no attempt; an interrupted reservation can consume one without
spending DFT time. Neither count is a backend SCF count, electronic iteration
count, core-hours, full preparation cost or successful material prediction.
Continue to audit native outputs and actual resource use separately.

For measured performance, distinguish an attempt reservation, a completed
calculator callback, a fresh native SCF and an allocation's core-hours.
Charge the full same-chain history, not just its final converged segment;
count shared exact source calls once, but do not merge different new calls
merely because their structures or result hashes agree. Endpoint/reference
preparation, failed starts and local conditional relaxations are separate
cost categories, not free work.

The [seven-chain HfO2 cost example](../benchmarks/hfo2_channels/20261008/clamped_chain_costs_E147_20261014/README.md)
reuses the existing `ResponseEvaluationCost` API and byte-pinned native
receipts to separate 1918 actual production SCFs from cached frames and
sixteen allocation durations. Run `python -m scripts.audit_hfo2_clamped_chain_costs
--output NEW_REPORT.json` from repository root to reproduce that declared
completed-production subset. It explicitly excludes other study costs and
does not constitute an optimizer speed-up or prediction benchmark.

Use one writer per workdir. Exclusive markers refuse a stale writer rather
than silently sharing slots. They are not a cross-workdir/distributed budget
allocator or proof against external deletion/tampering or storage power loss.
The caller still owns job/resource limits, preservation of original records
and raw physical audits. Budget refusal never changes optimization/electronic
parameters or fills a failed grid point with a fabricated energy. For this
HfO2 study, no material limit or G3 mode/centre is chosen before complete G2.

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

The [2026-10-14 eight-path G2 audit](../benchmarks/hfo2_channels/20261008/clamped_plus1_reversing_terminal_E148_20261014/README.md)
now completes the ordinary training-path matrix. This is not G3 validation:
physical barrier-error bounds, stable conditional branches and independent
predictions are still absent. At that audit checkpoint, `BiaxialClampedCurvatureCoordinates`
provided directional sampling/physical gradients, not the mode evaluator's
`to_atoms`/`energy_and_gradient` interface. A controlled-strain conditional
surface must bridge that interface and keep the prescribed substrate parameter
out of the released complement; passing a same-size array is not sufficient.

## Explicit prescribed-strain conditional coordinates

`BiaxialClampedModeCoordinates` now supplies that interface while reusing
`BiaxialClampedCurvatureCoordinates` for physical geometry and gradients.
It does not replace the fixed-substrate chart or change its cache contract.
For an anchor strain epsilon0 and a declared scale L_strain in Angstrom,
the extra coordinate is `t = L_strain*(epsilon-epsilon0)`. The substrate
rows are `H0[:2]*(1+epsilon)/(1+epsilon0)`; open-row directions and ordered
atomic coordinates retain the original chart. There is no automatic metric,
strain range, mode selection or material surface construction.

```python
from vcneb import BiaxialClampedModeCoordinates, CalculatorModeEvaluator

chart = BiaxialClampedModeCoordinates(
    reference_atoms, boundary=boundary, cell_scale_A=registered_cell_scale,
    anchor_strain=registered_anchor_strain, strain_scale_A=registered_strain_scale,
)
evaluator = CalculatorModeEvaluator(
    chart, existing_fixed_ASE_factory, workdir=fresh_directory,
    calculator_id=registered_physical_contract_id,
    result_validator=existing_raw_SCF_validator, validation_id=registered_validator_id,
    max_evaluation_attempts=registered_limit,
)
```

The strain scale, anchor, controlled row map and active-cell basis are part
of the immutable physical cache identity. Changed definitions reject resume.
`from_atoms` requires the same ordered continuous lift and one common positive
biaxial row factor; it never fits a rotation or wraps atoms.
Raw calculation geometry must agree with the declared unwrapped point.
Stress retains substrate reactions: `L_strain*gradient[-1]` is the frozen
partial d(E+PV)/d epsilon, not an equilibrium response until internal stationarity
is independently established. The controlled gradient need not vanish.

The evaluator exposes `fixed_coordinate_covectors`. Before relaxation or
orthogonal-curvature sampling, the solver checks `C.T @ Z = 0` for every
released column Z and prescribed covector C. Merely mixing a strain direction
into a Q axis is insufficient and is rejected before any calculator call.
Retain t as a pure axis, include it in the full constrained covector span, or
explicitly freeze its direction with a compatible metric/plane. The three
rigid translations should still be treated as a declared gauge. A transparent
callable wrapper must forward this metadata; a bare lambda cannot carry an
undeclared mechanical contract. No control is silently frozen or a supplied
mode plane silently changed. Whole-grid calls use the same point-level guard.

This is implementation validation, not a new HfO2 G3 surface, stationary
branch/TS or prediction certificate. EMT numerical derivatives/conditional
relaxations check the full physical interface; retained native zero/+1%
frames test common-reference replay without new DFT or physical+.005access.
The original ABACUS100Ry/full10auDZP and finite study budgets stay unchanged.

The [observed cross-strain work example](../benchmarks/hfo2_channels/20261008/cross_strain_work_E151_20261014/README.md)
reuses this controlled chart with the existing0/+1% native structures. Its
dual endpoint contractions independently match Cartesian force/stress work,
without intermediate DFT or held-out structures. The reported connection-
dependent work terms and quadrature residuals are not unique mode energies,
stationary barrier derivatives or material uncertainty bounds. Complete
harmonic release is also checked for invariance under retained/reference
basis changes: coordinate renaming alone cannot establish prediction gain.
