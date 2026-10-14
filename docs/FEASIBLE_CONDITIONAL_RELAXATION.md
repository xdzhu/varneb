# Experimental ball-feasible conditional relaxation

`vcneb.feasible_relaxation.relax_orthogonal_in_ball` is an opt-in,
backend-independent, single-start experiment. Existing optimizers, defaults
and frozen production traces are unchanged. It consumes the same full
physical energy/dual-gradient interface and returns the existing
`ConditionalModePoint` representation. It neither launches DFT nor selects
a calculator, pressure, mechanical boundary or electronic setting.

The experiment follows the native HfO2 [E168 diagnosis](../benchmarks/hfo2_channels/20261008/G3_bound_diagnostic_E168_20261015/README.md):
many proposed steps hit the declared search-radius guard while a substantial
constant-radius tangent gradient remained. This motivates a feasible search
step, not a larger domain or looser gradient target. Projection-based steps
and BFGS are established optimization ingredients, not a new theorem here.

## Contract and example

```python
from vcneb.feasible_relaxation import relax_orthogonal_in_ball

candidate = relax_orthogonal_in_ball(
    plane, q, energy_and_gradient,
    start=audited_start,
    orthogonal_radius_bound=declared_released_radius,
    trial_validator=original_geometry_guard,
    initial_hessian=original_full_coordinate_hessian,
    initial_trust_radius=0.05,
    gradient_tolerance=0.03,
)
```

The example's geometry, Hessian, radius and target must be supplied from a
declared study. It is not a runnable material recipe with guessed values.
The original start must already preserve Q, frozen directions and the
geometry/domain guard. It is not silently projected into a valid structure.
Prescribed controls must be excluded from the release space before callbacks.

The ball bounds the Euclidean norm of coefficients in the
metric-orthonormal **released** basis about `plane.frozen_coordinates(q)`.
It is not automatically a total atomic/cell displacement limit. If the
original total metric ball is centered at `plane.reference`, orthogonality
gives the available released radius

`R_released = sqrt(R_total**2 - ||frozen(q) - reference||_metric**2)`.

The fixed contribution must be smaller than R_total. The same original
total-radius guard must still check every physical proposal. This formula
does not apply without its common center, metric and orthogonal-embedding
premises. The numerical search domain is not a support reaction or a reason
to remove a physically movable gradient component.

## Step and stopping behavior

A projected quasi-Newton endpoint stays within the declared released ball.
If it is not a descent displacement for the actual released gradient, the
solver falls back to projected steepest descent. Armijo backtracking uses
the actual displacement and retains the independent original geometry guard.
The Hessian seed uses the full physical-coordinate convention; its released
restriction must pass Cholesky without eigenvalue clipping.

Success still requires the **full released-gradient norm** below the target,
not a projected gradient or boundary KKT residual. A lower-energy boundary
point with nonzero physical gradient is not returned as a successful
conditional candidate. The original evaluator controls unique-call budgets;
budget/backend errors propagate without retry or input repair. Callback cost
includes final verification, while exact-cache hits may not cost new SCFs.
Positive released curvature, branch stability and uniqueness remain separate
checks, just as for the existing conditional solver.

## Current evidence and boundary

The [frozen four-case analytic comparison](../benchmarks/hfo2_channels/20261008/feasible_steps_E169_20261015/README.md)
uses matched starts, seeds, targets, domains and twelve unique-call caps.
One deliberately misaligned near-wall quadratic uses7 versus6 oracle calls;
interior and nonlinear controls show no call reduction. Both arms refuse an
exterior-well boundary point as convergence. An analytic ASE energy/force/
stress adapter also retains its original conservative factory budget.
These are implementation/control tests, not HfO2 DFT acceleration, a new
production default, or independent material prediction. Actual same-input,
bounded material experiments are still required before making such claims.

The [E170 matched cache/input contract](../benchmarks/hfo2_channels/20261008/G3_feasible_material_E170_20261015/README.md)
and [E171 bounded execution supplement](../benchmarks/hfo2_channels/20261008/G3_feasible_execution_E171_20261015/README.md)
now provide a concrete HfO2 training experiment. They reuse the exact old
start/cache and retain the original physics, target and total-domain guard.
The native zero-DFT probe passes; the prospective arm has a separate12new-
attempt cap and must wait for the original eight frozen-R controls. Local
and independent execution checks are not material acceleration evidence.
The existing production/default solver and live source remain unchanged.
