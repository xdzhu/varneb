"""Geometry-limited branch relaxation with an explicit directional stop test.

This numerical optimizer neither detects a stable branch nor certifies a TS.
The evaluator owns physical settings, provenance, caching and call budgets.
"""
from __future__ import annotations

import numpy as np


def relax_branch_candidate(evaluate, start, free_basis, branch_direction, *,
                           initial_hessian, valid, gradient_tolerance=0.03,
                           branch_gradient_tolerance=1e-4,
                           initial_trust_radius=0.05, max_iterations=200):
    """Relax in a linear free space, also resolving one shallow direction.

    All coordinates use an explicit Euclidean length chart. ``evaluate``
    returns energy and its full dual gradient. ``initial_hessian`` is a
    positive numerical preconditioner on the free space, not necessarily a
    measured physical Hessian. Constraints remain fixed by construction.
    Both stopping tests must pass. Native evaluator exceptions propagate;
    no retries or parameter changes are performed here.
    """
    if not callable(evaluate) or not callable(valid):
        raise TypeError('evaluator and geometry guard must be callable')
    if any(np.iscomplexobj(a) for a in (start, free_basis, branch_direction, initial_hessian)):
        raise ValueError('real coordinates, directions and preconditioner required')
    anchor = np.array(start, dtype=float, copy=True)
    basis = np.asarray(free_basis, dtype=float)
    direction = np.asarray(branch_direction, dtype=float)
    hessian = np.asarray(initial_hessian, dtype=float)
    n = anchor.size
    if (anchor.shape != (n,) or basis.ndim != 2 or basis.shape[0] != n
            or basis.shape[1] < 1 or direction.shape != (n,) or hessian.shape != (n, n)
            or not all(np.isfinite(a).all() for a in (anchor, basis, direction, hessian))
            or not np.allclose(basis.T @ basis, np.eye(basis.shape[1]), atol=2e-12, rtol=0)
            or not np.allclose(basis @ (basis.T @ direction), direction, atol=2e-12, rtol=0)
            or not np.isclose(direction @ direction, 1., atol=2e-12, rtol=0)
            or not np.allclose(hessian, hessian.T, atol=2e-12, rtol=0)):
        raise ValueError('finite orthonormal free space and unit open branch direction required')
    scalars = (gradient_tolerance, branch_gradient_tolerance, initial_trust_radius)
    if any(not np.isfinite(v) or v <= 0 for v in scalars):
        raise ValueError('positive finite tolerances and trust radius required')
    if isinstance(max_iterations, bool) or int(max_iterations) != max_iterations or max_iterations < 1:
        raise ValueError('positive integer iteration cap required')
    reduced = basis.T @ hessian @ basis
    factor = np.linalg.cholesky(0.5 * (reduced + reduced.T))
    seed = np.linalg.solve(factor.T, np.linalg.solve(factor, np.eye(factor.shape[0])))
    if not np.isfinite(seed).all():
        raise ValueError('finite positive preconditioner inverse required')
    inverse = seed.copy()
    z = np.zeros(basis.shape[1])
    calls = 0
    trace = []

    def at(value):
        nonlocal calls
        x = anchor + basis @ value
        if not valid(x):
            raise ValueError('evaluation geometry must pass before calculator invocation')
        energy, full_gradient = evaluate(x)
        full_gradient = np.asarray(full_gradient, dtype=float)
        if (not np.isfinite(energy) or full_gradient.shape != anchor.shape
                or not np.isfinite(full_gradient).all()):
            raise ValueError('finite energy and full physical gradient required')
        calls += 1
        return float(energy), full_gradient

    energy, full_gradient = at(z)
    radius = float(initial_trust_radius)
    status = 'iteration_cap'
    for iteration in range(int(max_iterations) + 1):
        gradient = basis.T @ full_gradient
        norm = float(np.linalg.norm(gradient))
        monitored = float(direction @ full_gradient)
        trace.append(dict(iteration=iteration, energy=energy,
                          released_gradient_norm=norm, branch_gradient=monitored))
        if norm <= gradient_tolerance and abs(monitored) <= branch_gradient_tolerance:
            status = 'gradient_candidate_stability_unchecked'
            break
        if iteration == max_iterations:
            break
        step_direction = -inverse @ gradient
        slope = float(step_direction @ gradient)
        if slope >= -1e-12 * np.linalg.norm(step_direction) * norm:
            inverse = seed.copy()
            step_direction = -gradient
        length = np.linalg.norm(step_direction)
        if length == 0:
            status = 'zero_step_without_convergence'
            break
        step_direction *= min(1., radius / length)
        accepted = None
        step = 1.
        for _ in range(20):
            proposal = z + step * step_direction
            if valid(anchor + basis @ proposal):
                next_energy, next_full = at(proposal)
                if next_energy <= energy + 1e-4 * step * float(gradient @ step_direction):
                    accepted = proposal, next_energy, next_full
                    break
            step *= .5
        if accepted is None:
            status = 'line_search_stopped_without_convergence'
            break
        next_z, next_energy, next_full = accepted
        displacement = next_z - z
        change = basis.T @ (next_full - full_gradient)
        curvature = float(displacement @ change)
        if curvature > 1e-10 * np.linalg.norm(displacement) * np.linalg.norm(change):
            rho = 1. / curvature
            transform = np.eye(z.size) - rho * np.outer(displacement, change)
            inverse = transform @ inverse @ transform.T + rho * np.outer(displacement, displacement)
        z, energy, full_gradient = next_z, next_energy, next_full
        if step >= .75:
            radius = min(2 * initial_trust_radius, 1.25 * radius)
    return dict(status=status, coordinates=(anchor + basis @ z).tolist(), energy=energy,
                released_gradient_norm=float(np.linalg.norm(basis.T @ full_gradient)),
                branch_gradient=float(direction @ full_gradient), evaluator_calls=calls,
                trace=trace, stability_or_TS_certified=False)
