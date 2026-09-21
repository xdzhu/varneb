"""Calculator-independent VARNEB path-optimization strategies.

The calculator/backend boundary ends at ``energy``, ``forces`` and ``stress``.
This registry makes the other side of that boundary explicit: optimizer choice
and convergence controls are properties of the VARNEB path controller and can
be changed without changing the calculator profile or launcher.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class OptimizerSpec:
    """Public metadata for one calculator-independent path optimizer."""

    name: str
    family: str
    notes: str


_OPTIMIZERS: tuple[OptimizerSpec, ...] = (
    OptimizerSpec("FIRE", "fire", "robust default for ordinary VC-NEB"),
    OptimizerSpec("BlockFIRE", "fire", "independent image trust radii"),
    OptimizerSpec("SplitFIRE", "fire", "independent atomic and cell trust radii"),
    OptimizerSpec("ImageScaledFIRE", "fire", "image-local force scaling"),
    OptimizerSpec("StagedFIRE", "fire", "coarse-to-fine FIRE stages"),
    OptimizerSpec("BFGS", "quasi_newton", "quasi-Newton path updates"),
    OptimizerSpec("LBFGS", "quasi_newton", "limited-memory quasi-Newton updates"),
    OptimizerSpec("BFGSLineSearch", "quasi_newton", "line-search quasi-Newton updates"),
)


def optimizer_specs() -> tuple[OptimizerSpec, ...]:
    """Return the immutable list of available path strategies."""

    return _OPTIMIZERS


def get_optimizer_spec(name: str) -> OptimizerSpec:
    """Resolve a user-facing optimizer name, accepting legacy separators."""

    key = str(name).strip().upper().replace("-", "").replace("_", "")
    for spec in _OPTIMIZERS:
        if spec.name.upper().replace("-", "").replace("_", "") == key:
            return spec
    choices = ", ".join(spec.name for spec in _OPTIMIZERS)
    raise ValueError(f"unknown VARNEB optimizer {name!r}; choose one of: {choices}")


def optimizer_capability_matrix() -> list[dict[str, str]]:
    """Return JSON-friendly strategy metadata independent of backend choice."""

    return [{"name": spec.name, "family": spec.family, "notes": spec.notes, "backend_independent": "true"}
            for spec in _OPTIMIZERS]


__all__ = ["OptimizerSpec", "get_optimizer_spec", "optimizer_specs", "optimizer_capability_matrix"]
