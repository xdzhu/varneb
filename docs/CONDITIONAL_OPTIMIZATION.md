# Conditional optimization: optional curvature initialization

`relax_orthogonal_at_q` and `sample_conditional_mode_surface` accept an
optional `initial_hessian` for `optimizer="safeguarded_bfgs"`. Default None
preserves identity initialization and the existing optimizer defaults.
This is standard numerical preconditioning, not a novel physical algorithm.

```python
point = relax_orthogonal_at_q(
    plane, q, fixed_physical_evaluator,
    frozen_directions=declared_gauge_and_control_directions,
    optimizer="safeguarded_bfgs",
    trial_validator=declared_geometry_and_domain_guard,
    initial_trust_radius=declared_step_radius,
    initial_hessian=approximate_full_coordinate_hessian,
    gradient_tolerance=declared_conditional_tolerance,
)
```

The supplied matrix uses the full physical coordinates of the plane and the
dual physical energy gradient. For a local length chart its units are
eV/Angstrom². It is not a mass-weighted phonon matrix, inverse Hessian, NEB
force Jacobian or matrix already reduced to an unspecified complement basis.
With metric-orthonormal open columns U, initialize the inverse of U.T H0 U.
No extra multiplication by the coordinate metric belongs in that Hessian
restriction. Nonunit-metric constrained quadratic tests verify this convention.

Reject wrong shape, complex/nonfinite/asymmetric inputs, a nonpositive open
restriction or a numerically unrepresentable inverse before energy callbacks.
Negative curvature in the fixed Q/gauge/control directions is allowed. This
numerical admissibility of the supplied seed is **not** independent evidence
that the physical released block is positive, a center is stationary, or a
conditional branch is a minimum. Do not treat unmeasured Hessian columns as
physical zeros. Measured-model coverage and physical stability/error audits
remain separate from choosing an approximate optimizer preconditioner.

The seed does not modify the callback energy/force/stress or final tolerance.
Existing geometry/domain guards, Armijo backtracking, step radius, fixed-Q
constraints and per-start outcomes still apply. It is rejected for L-BFGS-B
or an entirely prescribed point, where it would otherwise be silently unused.

Use the existing `CalculatorModeEvaluator.max_evaluation_attempts` to bound
cumulative calculator attempts, including failures; budgets are not refunded
or enlarged by this optimizer. Final independent callback checks and exact
cache hits need not be new SCFs. Report callbacks, true backend calls and
source preparation costs separately, including unsuccessful starts.

The E157 HfO2 replay uses the measured E155 H/g only to define an **exact
quadratic surrogate**. Identity and exact-curvature seeds reach the same
analytic conditional solution with different callback counts. This validates
implementation, not DFT acceleration or the actual conditional surface.
Real finite conditional branches and same-accuracy total-cost tests remain.

The [E158 HfO2 canary worker](../benchmarks/hfo2_channels/20261008/conditional_canary_worker_E158_20261014/README.md)
provides a bounded same-input native test route. It refuses incomplete
two-step data, reuses the audited origin exactly and keeps two independent
optimizer arms/cost caps. Implementation checks alone do not start DFT or
complete a conditional face; source/admission and physical audits still apply.
