"""Opt-in ball-feasible conditional relaxation; not a new default or TS test.

The ball is a declared metric search domain, not an additional mechanical
constraint. Boundary KKT convergence never replaces the full released-gradient
target. The existing mode_surface solvers and their frozen traces are unchanged.
"""
from __future__ import annotations
from typing import Callable
import numpy as np
from .mode_surface import (
    ModePlane,
    ConditionalModePoint,
    ConditionalBranchOutcome,
    _orthogonal_directions,
    _require_fixed_controls,
)


def _ball_projection(z, radius):
    norm = float(np.linalg.norm(z))
    return z.copy() if norm <= radius else z * np.nextafter(radius / norm, 0.0)


def relax_orthogonal_in_ball(
    plane: ModePlane,
    q,
    energy_and_gradient: Callable,
    *,
    start,
    orthogonal_radius_bound: float,
    trial_validator: Callable,
    frozen_directions=None,
    gradient_tolerance: float = 1e-6,
    max_iterations: int = 200,
    initial_trust_radius: float = 0.25,
    initial_hessian=None,
) -> ConditionalModePoint:
    """One explicitly supplied branch with projected, guarded BFGS steps.

    The radius bounds metric-normalized released coefficients about the
    frozen-Q point, not the full physical displacement. For a total-radius R
    domain around plane.reference, callers must derive and document the
    available released radius after subtracting fixed-coordinate contributions.
    A physical geometry guard remains mandatory and precedes every callback.
    Initial starts are checked, not repaired. Q and frozen controls are retained.

    Projection supplies a feasible endpoint. Armijo uses the actual projected
    displacement; a non-descent quasi-Newton proposal falls back to projected
    steepest descent. Full released gradients determine success, including at
    the boundary. A missing interior stationary candidate remains a failure,
    never a certified constrained minimum. No calculator or scheduler is used.
    """
    for name, value in (
        ("orthogonal_radius_bound", orthogonal_radius_bound),
        ("gradient_tolerance", gradient_tolerance),
        ("initial_trust_radius", initial_trust_radius),
    ):
        if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if (
        isinstance(max_iterations, bool)
        or int(max_iterations) != max_iterations
        or max_iterations < 1
    ):
        raise ValueError("max_iterations must be a positive integer")
    if not callable(energy_and_gradient) or not callable(trial_validator):
        raise ValueError("energy_and_gradient and trial_validator must be callable")
    if np.iscomplexobj(q) or np.iscomplexobj(start):
        raise ValueError("real Q and start required")
    q = np.asarray(q, dtype=float)
    if q.shape != (2,) or not np.isfinite(q).all():
        raise ValueError("two finite Q coordinates required")
    base = plane.frozen_coordinates(q)
    u = _orthogonal_directions(plane, frozen_directions)
    _require_fixed_controls(energy_and_gradient, u)
    if not u.shape[1]:
        raise ValueError("ball relaxation requires at least one released direction")
    x = np.asarray(start, dtype=float)
    if x.shape != base.shape or not np.isfinite(x).all():
        raise ValueError("finite full-coordinate start required")
    z = u.T @ (plane.metric_weights * (x - base))
    if not np.allclose(base + u @ z, x, atol=2e-12, rtol=0):
        raise ValueError("start must already preserve fixed Q and frozen directions")
    if np.linalg.norm(z) > orthogonal_radius_bound or not trial_validator(base + u @ z):
        raise ValueError("start violates original domain or geometry guard; no repair")
    inverse = np.eye(u.shape[1])
    if initial_hessian is not None:
        if np.iscomplexobj(initial_hessian):
            raise ValueError("initial_hessian must be real")
        h = np.asarray(initial_hessian, dtype=float)
        if (
            h.shape != (base.size, base.size)
            or not np.isfinite(h).all()
            or not np.allclose(h, h.T, atol=1e-12, rtol=1e-10)
        ):
            raise ValueError(
                "finite symmetric full-coordinate initial_hessian required"
            )
        reduced = u.T @ h @ u
        if not np.isfinite(reduced).all():
            raise ValueError("seed restriction not representable")
        try:
            factor = np.linalg.cholesky(0.5 * reduced + 0.5 * reduced.T)
        except np.linalg.LinAlgError as exc:
            raise ValueError(
                "released seed must be positive definite; no clipping"
            ) from exc
        inverse = np.linalg.solve(factor.T, np.linalg.solve(factor, np.eye(u.shape[1])))
        if not np.isfinite(inverse).all():
            raise ValueError("seed inverse not representable")
    seed = inverse.copy()
    calls = 0

    def evaluate(value):
        nonlocal calls
        point = base + u @ value
        if np.linalg.norm(value) > orthogonal_radius_bound or not trial_validator(
            point
        ):
            raise ValueError("refuse unsafe callback; no geometry repair")
        energy, gradient = energy_and_gradient(point)
        calls += 1
        if np.iscomplexobj(gradient):
            raise ValueError("real physical gradient required")
        gradient = np.asarray(gradient, dtype=float)
        if (
            gradient.shape != base.shape
            or not np.isfinite(energy)
            or not np.isfinite(gradient).all()
        ):
            raise ValueError("finite energy and full physical gradient required")
        return float(energy), u.T @ gradient

    energy, g = evaluate(z)
    radius = initial_trust_radius
    for _ in range(int(max_iterations)):
        if np.linalg.norm(g) <= gradient_tolerance:
            break
        direction = -inverse @ g
        norm = float(np.linalg.norm(direction))
        if not np.isfinite(norm):
            raise ValueError("nonfinite numerical search direction")
        if norm == 0:
            break
        direction *= min(1.0, radius / norm)
        direction = _ball_projection(z + direction, orthogonal_radius_bound) - z
        if g @ direction >= -1e-12 * np.linalg.norm(g) * np.linalg.norm(direction):
            inverse = seed.copy()
            direction = -g * min(1.0, radius / np.linalg.norm(g))
            direction = _ball_projection(z + direction, orthogonal_radius_bound) - z
        if not g @ direction < 0:
            break
        trial = None
        step = 1.0
        for _ in range(20):
            proposed = z + step * direction
            if np.linalg.norm(proposed) <= orthogonal_radius_bound and trial_validator(
                base + u @ proposed
            ):
                ne, ng = evaluate(proposed)
                if ne <= energy + 1e-4 * (g @ (proposed - z)):
                    trial = (proposed, ne, ng)
                    break
            step *= 0.5
        if trial is None:
            break
        next_z, ne, ng = trial
        s = next_z - z
        y = ng - g
        curvature = float(s @ y)
        if curvature > 1e-10 * np.linalg.norm(s) * np.linalg.norm(y):
            rho = 1.0 / curvature
            t = np.eye(z.size) - rho * np.outer(s, y)
            inverse = t @ inverse @ t.T + rho * np.outer(s, s)
        z, energy, g = next_z, ne, ng
        if step >= 0.75:
            radius = min(2 * initial_trust_radius, 1.25 * radius)
    energy, g = evaluate(z)
    residual = float(np.linalg.norm(g))
    if residual > gradient_tolerance:
        raise RuntimeError("no full released-gradient convergence within declared ball")
    outcome = ConditionalBranchOutcome(
        start_index=0,
        coordinates=base + u @ z,
        energy=energy,
        orthogonal_gradient_norm=residual,
        n_evaluations=calls,
        converged=True,
    )
    return ConditionalModePoint(
        q=q.copy(),
        coordinates=base + u @ z,
        energy=energy,
        orthogonal_gradient_norm=residual,
        n_evaluations=calls,
        n_starts=1,
        selected_start=0,
        branch_outcomes=(outcome,),
    )
