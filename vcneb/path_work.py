"""Signed discrete physical-gradient work, not a mode energy partition.

Coordinates and covectors must be a conjugate pair in a caller-declared chart.
Endpoint trapezoids only approximate its segment work; disagreement with raw
energies is a quadrature/consistency diagnostic, never an uncertainty bound,
barrier prediction or certification of the unsampled path.
"""
from __future__ import annotations

from collections.abc import Mapping
import numpy as np


def sampled_path_work(coordinates, gradients, energies, subspaces: Mapping[str, np.ndarray]):
    """Integrate each signed subspace contribution along the sampled polygon.

Columns define orthonormal, mutually orthogonal subspaces in the supplied
coordinate representation. Their union must cover it completely. Do not mix
incommensurate units in a basis; separate atomic mass-coordinate and scaled
cell blocks can nevertheless be contracted with their respective covectors.
Negative work is retained. Basis rotations within each span leave the result
unchanged. No coordinate gauge, periodic lift or calculator data is repaired.
"""
    if any(np.iscomplexobj(value) for value in (coordinates, gradients, energies)):
        raise ValueError("explicit real coordinates, covectors and energies required")
    x, g, energy = [np.asarray(value, dtype=float) for value in (coordinates, gradients, energies)]
    if (x.ndim != 2 or min(x.shape) < 1 or len(x) < 2 or g.shape != x.shape
            or energy.shape != (len(x),)
            or any(not np.isfinite(a).all() for a in (x, g, energy))):
        raise ValueError("matching finite sampled coordinates/covectors and energies required")
    if not isinstance(subspaces, Mapping) or not subspaces:
        raise ValueError("explicit complete named subspaces required")
    bases = {}
    for name, value in subspaces.items():
        if not isinstance(name, str) or not name.strip() or name != name.strip():
            raise ValueError("nonempty trimmed subspace names required")
        if np.iscomplexobj(value):
            raise ValueError("explicit real subspace columns required")
        basis = np.asarray(value, dtype=float)
        if (basis.ndim != 2 or basis.shape[0] != x.shape[1] or not np.isfinite(basis).all()
                or not np.allclose(basis.T@basis, np.eye(basis.shape[1]), atol=1e-8, rtol=0)):
            raise ValueError("finite orthonormal subspace columns required")
        bases[name] = basis
    union = np.concatenate(list(bases.values()), axis=1)
    if (union.shape != (x.shape[1], x.shape[1])
            or not np.allclose(union.T@union, np.eye(x.shape[1]), atol=1e-8, rtol=0)):
        raise ValueError("subspaces must be a complete mutually orthogonal partition")
    delta, average = np.diff(x, axis=0), .5*g[:-1]+.5*g[1:]
    total = np.einsum("ni,ni->n", average, delta)
    parts = {name: np.einsum("nk,nk->n", average@basis, delta@basis)
             for name, basis in bases.items()}
    part_sum = np.sum(list(parts.values()), axis=0)
    defect = part_sum-total
    increments = np.diff(energy)
    residual = increments-total
    if any(not np.isfinite(v).all() for v in (total, defect, increments, residual, *parts.values())):
        raise ValueError("sampled-work arithmetic overflow")
    # This is an identity check, not a physical accuracy/energy-error tolerance.
    scale = max(1., float(np.max(np.abs(total))),
                float(np.max(np.sum(np.abs(list(parts.values())), axis=0))))
    if float(np.max(np.abs(defect))) > 1e-8*scale:
        raise ValueError("subspace work sum does not reproduce full sampled work")
    cumulative = np.r_[0., np.cumsum(total)]
    energy_change = energy-energy[0]
    cumulative_parts = {name:np.r_[0.,np.cumsum(value)] for name,value in parts.items()}
    if any(not np.isfinite(v).all() for v in (cumulative,energy_change,energy_change-cumulative,*cumulative_parts.values())):
        raise ValueError("cumulative sampled-work arithmetic overflow")
    return dict(
        method="endpoint-gradient trapezoids along the declared coordinate polygon",
        segment_energy_increments=increments.tolist(), segment_full_work=total.tolist(),
        segment_partition_work={name:value.tolist() for name,value in parts.items()},
        segment_energy_minus_work=residual.tolist(),
        cumulative_full_work=cumulative.tolist(),
        cumulative_partition_work={name:value.tolist() for name,value in cumulative_parts.items()},
        cumulative_energy_change=energy_change.tolist(),
        cumulative_energy_minus_work=(energy_change-cumulative).tolist(),
        maximum_partition_sum_defect=float(np.max(np.abs(defect))),
        maximum_segment_energy_minus_work=float(np.max(np.abs(residual))),
        maximum_cumulative_energy_minus_work=float(np.max(np.abs(energy_change-cumulative))),
        endpoint_energy_minus_work=float(energy_change[-1]-cumulative[-1]),
        physical_error_bound_measured=False,
        limitations=["signed path-dependent work, not state-function energies assigned to modes",
                     "not an error bound, NEB-force integral, conditional surface or stationary saddle",
                     "a coarse quadrature residual does not isolate stress/force error from unsampled nonlinearity"],
    )
