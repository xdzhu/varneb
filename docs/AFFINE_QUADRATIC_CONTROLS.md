# Shared-coordinate harmonic controls

`vcneb.affine_quadratic.constrain_quadratic` is a calculator-independent,
analysis-only utility for a standard equality-constrained quadratic:

```text
minimize g.x + x.H.x/2
subject to x in span(U), A.T x = q.
```

U contains orthonormal allowed directions in one declared physical metric.
A defines the observables; it need not be orthonormal and is not rescaled.
Directions outside U stay clamped. Energy, gradient, Hessian and coordinates
must share one reference, Hamiltonian, boundary and unit convention.

```python
import numpy as np
from vcneb.affine_quadratic import constrain_quadratic

H = np.diag([-1., 2., 3.])
g = np.array([.1, -.2, .3])
model = constrain_quadratic(
    H, g, np.eye(3), np.eye(3)[:, :1], stability_floor=0.)
energy_change = model.energy_change([.1])
x = model.coordinates([.1])
gradient_q = model.retained_gradient([.1])
residuals = model.residuals([.1])
```

This analytic example uses a raw zero floor, not a DFT stability certificate.
The constraint-null Hessian must exceed the caller's explicit floor. An
unstable or unresolved free direction is rejected, not pseudoinverted,
clipped or repaired. An empty free space is a constraint identity, not a
physical stability measurement. The affine gradient is retained. Array
fields are immutable copies; caller mutations do not change the model.

The KKT solution supplies allowed-space stationarity and dE/dq=-lambda.
Residuals are algebraic diagnostics, not physical energy-error bounds.
The utility does not evaluate a calculator, check a geometry domain,
discover a nonlinear branch or certify a TS. Callers must perform their
actual trust-domain/minimum-distance checks and keep invalid predictions
as abstentions. It does not alter NEB/optimizer/backend defaults.

## Why the full-space control matters

An orthogonal coordinate change cannot create a new prediction from the
same complete quadratic. In particular, the complete fixed-reference
KKT solution and an equivalent stable Schur elimination are the same
constrained problem. Their equality is standard mathematics, not an
adaptive-model accuracy improvement or a new theory.

For a small reference span, say explicitly whether missing directions
are clamped or released. These are different models. Imposing observables
selected for another reference can also require large displacements even
when the formal constraint matrix has full rank. That is a compatibility/
domain limitation, not evidence that the reference is generally inadequate.
Do not rank errors from different eligible subsets as a matched benchmark.

The [HfO2 training comparison](../benchmarks/hfo2_channels/20261008/G3_shared_Q_controls_E180_20261015/README.md)
keeps both original finite-difference steps, the complete T optical span,
both complete T doublets, the geometric triplet and all four Cmma
registrations. It exports scalar diagnostics only. It is not an independent
material-prediction benchmark or the entire B0–B5 control matrix.
