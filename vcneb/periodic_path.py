"""Explicit ordered periodic lifts, never an automatic path repair.

The chosen policy is component-wise fractional minimum image at each adjacent
segment. It cannot infer an undersampled winding or an atom permutation. A
half-cell step is ambiguous and must be resolved by the caller. Calculators
are detached from returned copies; audited periodic-equivalent caches can be
reattached separately. The input chain is not modified.
"""

from __future__ import annotations

import numpy as np

from .core import validate_path_geometry


def minimum_image_path_lift(images):
    """Return copied images and recorded integer lattice shifts.

    Use only when short adjacent steps are the declared input convention.
    This is not Cartesian nearest-image search in an arbitrary skew cell.
    Already continuous, sufficiently sampled winding paths retain their
    winding: endpoints are not independently remapped to the first image.
    """
    if len(images) < 2:
        raise ValueError("periodic path lift requires at least two images")
    before = validate_path_geometry(images)
    symbols, pbc = images[0].get_chemical_symbols(), np.asarray(images[0].pbc, dtype=bool)
    result, shifts = [], []
    previous = None
    for i, original in enumerate(images):
        if original.get_chemical_symbols() != symbols or not np.array_equal(original.pbc, pbc):
            raise ValueError("periodic path lift requires identical ordered atoms and PBC")
        q = original.get_scaled_positions(wrap=False).copy()
        shift = np.zeros_like(q, dtype=np.int64)
        if previous is not None:
            delta = q - previous
            reduced = delta[:, pbc] - np.rint(delta[:, pbc])
            if np.any(np.isclose(np.abs(reduced), .5, atol=1e-10, rtol=0)):
                raise ValueError(f"ambiguous half-cell step at image {i}; provide a denser explicit lift")
            shift[:, pbc] = -np.rint(delta[:, pbc]).astype(np.int64)
        lifted = original.copy()
        lifted.calc = None
        if np.any(shift):
            lifted.set_scaled_positions(q + shift)
        previous = lifted.get_scaled_positions(wrap=False)
        difference = previous - q - shift
        if not np.allclose(difference, 0, atol=1e-12, rtol=0):
            raise ValueError("periodic lift changed geometry beyond integer lattice translations")
        result.append(lifted)
        shifts.append(shift.tolist())
    after = validate_path_geometry(result)
    report = {
        "policy": "explicit_adjacent_fractional_minimum_image",
        "integer_lattice_shifts_by_image_atom": shifts,
        "original_segment_lengths_A": before["segment_lengths_A"],
        "lifted_segment_lengths_A": after["segment_lengths_A"],
        "atomic_permutation_applied": False,
        "physical_geometry_changed": False,
        "calculator_results_implicitly_reused": False,
        "limitations": "short adjacent fractional steps assumed; undersampled winding is not inferred",
    }
    return result, report


def validate_periodic_path_lift(images):
    """Reject a chain inconsistent with the explicit short-step convention."""
    _, report = minimum_image_path_lift(images)
    if np.any(np.asarray(report["integer_lattice_shifts_by_image_atom"])):
        raise ValueError("supplied path lacks a continuous periodic lift; explicitly prepare and audit its integer lattice shifts before DFT")
    return report
