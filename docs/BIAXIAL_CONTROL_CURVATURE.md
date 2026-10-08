# Controlled substrate strain versus relaxed internal coordinates

This is a prerequisite for the registered HfO2 boundary-response experiment,
not a new elasticity theorem or a material prediction result. The existing
[fixed-plane interface](EPITAXIAL_BOUNDARY.md) remains the optimizer boundary.
No calculator setting, production image, or existing source archive is changed.

## Why the extra coordinate is needed

An endpoint or band at a prescribed epsilon may relax atoms and the allowed
third lattice vector, but it must not relax epsilon itself. The earlier
`ActiveJointCurvatureCoordinates.for_clamped_plane` correctly samples those
internal directions at **one fixed substrate**. That chart alone cannot sample
the mixed response to changing the substrate control parameter.

`BiaxialClampedCurvatureCoordinates` adds a separately tagged control direction
for derivative/probe preparation. It is not an ASE filter or an optimizer.
Continue using `ClampedPlaneFilter` and the corresponding VCNEB boundary at
each fixed epsilon. Never make the augmented space a freely movable NEB space.

## One explicit reference and metric

Let H0 be the reference row-cell matrix at anchor epsilon0, B_a the existing
exact open-cell deformation matrices, l the cell scale and L the explicitly
registered control scale. All coordinates y=(u_atom,u_cell,t) are in Angstrom:

```text
epsilon = epsilon0 + t/L
H[:2]   = H0[:2] * (1+epsilon)/(1+epsilon0)
H[2]    = (H0 + sum_a u_cell[a] * H0 @ B_a.T/l)[2]
s       = s0 + u_atom @ inv(H0)
r       = s @ H
```

The last coordinate is prescribed, not released. The first two rows are set
exactly; open-cell increments cannot produce an unintended substrate drift.
The anchor may be 0 or +1% and H0 may be oblique, tilted or rigidly rotated.
Atoms retain their order and unwrapped periodic lifts. Probe construction
drops the calculator; it never reuses free-cell energies at a changed substrate.
L changes the coordinate metric, not a cutoff, Hamiltonian or strain value;
record it before training/holdout comparisons rather than tuning it on labels.

For ASE tensile-positive stress sigma and compressive-positive P, the full
physical gradient uses the actual current cell:

```text
g_atom = -forces @ H.T @ inv(H0.T)
g_cell[a] = cell_work_derivative(sigma,H,dH/du_cell[a],pressure=P)
g_t       = cell_work_derivative(sigma,H,dH/dt,pressure=P)
d(E+PV)/d epsilon = L * g_t
```

These are partial derivatives. A stationary smooth internal branch is still
required before applying the envelope theorem to a relaxed energy or barrier.
The augmented Hessian's controlled direction is **not** an admissible physical
instability for a fixed-epsilon TS: saddle checks use the internal active
space, with the declared reaction direction and translations handled separately.

## Reuse the same probe and reduction machinery

```python
from vcneb import BiaxialClampedCurvatureCoordinates, joint_curvature_probes

chart = BiaxialClampedCurvatureCoordinates(
    audited_reference, boundary=registered_boundary,
    cell_scale_A=registered_cell_scale,
    anchor_strain=registered_epsilon,
    strain_scale_A=registered_control_scale,
)
# directions are predeclared orthonormal columns in the augmented metric.
probes = joint_curvature_probes(chart, directions, step_A=registered_step)
# After genuine same-contract evaluations, for each probe:
gradient = chart.enthalpy_gradient(
    probe.delta_A, evaluated_atoms, physical_forces, physical_stress, pressure,
)
```

Use the existing `assemble_joint_directional_curvature` with complete augmented
gradients. It keeps H*B, the raw projected matrix and couplings outside the
sampled directions; it does not certify unsampled self-curvatures. The API
refuses a geometry that does not match its declared point, including a changed
periodic lift, rather than silently projecting or remapping it.

`controlled_direction()` identifies the last column.
`internal_relaxation_basis()` has a zero last row and removes the three atomic
translations. It can supply the admissible internal space for the existing
quadratic reduction, after the actual eliminated block passes its stability
and resolution checks. Retained reaction/mode directions must be removed from
that release space explicitly. No pseudoinverse or automatic branch discovery
is introduced; B2--B5 prediction selection/freezing remains unfinished.

Substrate reaction stress is excluded from **internal convergence**, but must
be retained in **external strain work**. Dropping the reactions from g_t would
erase the measured boundary response even when all open tractions are zero.

## Validation and limits

`tests/test_biaxial_curvature.py` checks all coordinate derivatives against
independent Cu/EMT energy differences in eight combinations of tilt release,
rotation and deformed centre, including a nonzero anchor strain and pressure.
Two step sizes independently check mixed internal/control reciprocity without
discarding transverse gradient response. Other tests cover the control-scale
chain rule, the old fixed-strain limit, reaction work, stable synthetic release,
bad inputs, wrong geometries, and all ten actual uncomputed HfO2 training seeds.

The nonzero Cu/EMT pressure is a mathematical check, not a change to HfO2 P=0.
Ten real seed geometries are checked with **zero DFT calls**. No clamped HfO2
Hessian, boundary-response barrier, favorable window, certified TS or independent
material forecast follows from these tests. Source hashes, isolated regression
and current jobs are recorded in the [E027 receipt](../benchmarks/hfo2_channels/20261008/biaxial_control/validation_delivery.json).

The original G1 -> G2 -> G3 gates, +0.5% path holdout, two active chains and
finite DFT budgets stay unchanged. The sampler does not authorize evaluation
of a held-out condition or any additional production matrix.
