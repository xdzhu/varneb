# Joint-coordinate directional curvature: sampling versus inference

This backend-independent interface stages calculator-free paired geometries
and assembles physical enthalpy-gradient differences. It changes no DFT input,
performs no SCF, selects no phonon automatically and certifies no saddle.

## Define one chart and a declared direction set

Use `ActiveJointCurvatureCoordinates.for_clamped_plane` with the exact boundary
used by the endpoint and path. Other explicit allowed deformation bases,
including fixed cell, are supported by that chart. All atomic and scaled cell
coordinates are in Angstrom; directions are orthonormal columns in this same
metric. A mass-weighted phonon vector is not automatically a chart direction:
convert its convention explicitly before constructing the basis.

```python
from vcneb import (ActiveJointCurvatureCoordinates, joint_curvature_probes,
                   assemble_joint_directional_curvature)

chart = ActiveJointCurvatureCoordinates.for_clamped_plane(
    reference, cell_scale_A=cell_scale, boundary=boundary,
)
B = declared_training_basis  # shape (chart.size, k), B.T @ B = I
probes = joint_curvature_probes(chart, B, step_A=h)
# Ordered: direction0+, direction0-, direction1+, direction1-, ...
# For each static result, form the *full chart* gradient:
# g = chart.enthalpy_gradient(atoms, forces, stress_matrix, pressure_eV_A3)
# After raw-input/geometry/SCF audits, arrange g_plus/g_minus as (k, chart.size).
curvature = assemble_joint_directional_curvature(B, g_plus, g_minus, step_A=h)
```

`center_delta_A` may specify a nonzero chart center. Probe displacement tags are
copied/read-only; ASE atoms remain editable so the caller can attach a static
calculator after staging. Neither basis normalization nor per-image wrapping
is performed. Extra ASE atomic constraints are rejected. The optional
`candidate_validator` can reject steric/ill-conditioned geometries but cannot
repair them, change atom order or attach a calculator during staging.
Invalid steps, incomplete gradients and nonrepresentable differences fail closed.

## What the measured matrices mean

For orthonormal B and paired full-coordinate gradients,

`A[:,j] = [g(y0+h*B[:,j]) - g(y0-h*B[:,j])] / (2*h)`.

This estimates `H*B`. `raw_projected = B.T @ A` estimates `B.T*H*B`.
Its raw reciprocity defect and the operator change due to symmetrization are
retained; symmetrization must not hide an inconsistent force/stress adapter.
`transverse_action = A - B @ raw_projected` preserves coupling outside B.
Nonzero transverse action can be physical, not a numerical defect.

A k-direction experiment does **not** measure the self-curvature of the
unsampled complement. It cannot justify eliminating that complement, count
all negative directions, claim a conditional minimum or certify a full-space
TS. For harmonic elimination, measure the required retained+eliminated blocks
and use `relax_orthogonal_curvature` only when every eliminated eigenvalue
exceeds a declared, measured resolution floor. A pseudoinverse is not a cure
for an unstable branch. A nonstationary image still yields only local curvature.

Use physical `d(E+PV)/dy`, never NEB-projected or spring forces. Raw evaluations
must prove one ordered geometry line, chart center, pressure and calculator
contract. Array dimensions alone do not prove this provenance.

## Numerical and cost checks

Repeat the same basis at two declared steps. The existing
`audit_directional_curvature_consistency` cross-checks signed energies and
directional gradients. Compare projected operators and reciprocity; retain
disagreement rather than changing cutoff or electronic convergence settings.
Two-step spread is an operational resolution diagnostic, not a rigorous total
DFT error bar. If the basis changes, explicitly transform the common subspace
before comparing matrices; do not subtract mismatched coordinate arrays.

One amplitude requires exactly2*k static probes. Two amplitudes require4*k
plus a center if no identical audited cache exists. Twelve-atom tilt-released
HfO2 has36 translation-free joint directions:72 probes at one amplitude,
144 at two. Selected-direction sampling reduces calls but weakens inference;
it is not evidence that omitted directions are stable. These counts are not
a new production budget or a submitted G3 calculation matrix.

`python -m scripts.verify_joint_curvature_stencil --report /a/new/receipt.json`
runs without pytest or DFT. It checks720 inert geometries from the ten existing
unrelaxed HfO2 starters, a coupled analytic matrix and180 cheap Cu/EMT evaluations
at two steps in normal/tilt and oblique/rotated cases. Geometric and EMT step
sizes are implementation checks, not automatically chosen HfO2 SCF amplitudes.
This validation does not supply a material joint Hessian, independent barrier
prediction, convergence acceleration or TS certificate.
