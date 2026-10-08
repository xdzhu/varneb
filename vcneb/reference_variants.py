"""Apply a declared parent operation without erasing atom-path identity.

This is a geometric variant constructor, not a polarization calculation or
an irreducible-representation assignment. Parent sites define the permutation;
nearest-site remapping of the distorted product is deliberately not performed.
"""

from __future__ import annotations

from dataclasses import dataclass

from ase import Atoms
import numpy as np
from ase.geometry import find_mic


@dataclass(frozen=True)
class ReferenceVariant:
    atoms: Atoms
    source_to_target: tuple[int, ...]
    rotation: np.ndarray
    translation: np.ndarray


def apply_parent_operation(
    parent: Atoms,
    structure: Atoms,
    rotation,
    translation,
    *,
    mapping_tolerance_A: float = 1e-6,
    metric_tolerance: float = 1e-8,
) -> ReferenceVariant:
    """Return a variant in the original cell with parent-defined site labels.

    Fractional column coordinates obey ``f' = rotation @ f + translation``.
    The operation must map every parent site bijectively to a same-species
    parent site, and preserve the product-cell metric. Different product-cell
    orientations need a separately declared cell transformation, not this API.
    """
    if (len(parent) == 0 or parent.get_chemical_symbols() != structure.get_chemical_symbols()
            or not parent.pbc.all() or not structure.pbc.all()):
        raise ValueError("parent/product need identical ordered species and 3D periodicity")
    if (not np.isfinite(mapping_tolerance_A) or mapping_tolerance_A <= 0
            or not np.isfinite(metric_tolerance) or metric_tolerance <= 0):
        raise ValueError("positive finite mapping/metric tolerances required")
    for atoms in (parent, structure):
        if (not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell.array).all()
                or np.linalg.det(atoms.cell.array) <= 0):
            raise ValueError("finite, right-handed nonsingular geometry required")
    r = np.asarray(rotation, dtype=float)
    t = np.asarray(translation, dtype=float)
    if (r.shape != (3, 3) or t.shape != (3,) or not np.isfinite(r).all()
            or not np.isfinite(t).all() or not np.allclose(r, np.rint(r), atol=1e-12, rtol=0)
            or not np.isclose(abs(np.linalg.det(r)), 1, atol=1e-12, rtol=0)):
        raise ValueError("integer unimodular rotation and finite fractional translation required")
    for atoms in (parent, structure):
        metric = atoms.cell.array @ atoms.cell.array.T
        if not np.allclose(r.T @ metric @ r, metric, atol=metric_tolerance, rtol=metric_tolerance):
            raise ValueError("operation does not preserve parent/product metric")
    sites = parent.get_scaled_positions(wrap=False)
    transformed_parent = sites @ r.T + t
    mapping = []
    for i, position in enumerate(transformed_parent):
        vectors = (sites - position) @ parent.cell.array
        _, distances = find_mic(vectors, parent.cell.array, pbc=True)
        candidates = np.flatnonzero((distances <= mapping_tolerance_A)
                                   & (parent.numbers == parent.numbers[i]))
        if len(candidates) != 1:
            raise ValueError("parent operation lacks an unambiguous species-preserving site match")
        mapping.append(int(candidates[0]))
    if len(set(mapping)) != len(parent):
        raise ValueError("parent operation is not a bijective site mapping")
    transformed = structure.get_scaled_positions(wrap=False) @ r.T + t
    ordered = np.empty_like(transformed)
    ordered[np.asarray(mapping)] = transformed
    variant = structure.copy()
    variant.calc = None
    variant.set_scaled_positions(ordered % 1.0)
    return ReferenceVariant(variant, tuple(mapping), r.copy(), t.copy())
