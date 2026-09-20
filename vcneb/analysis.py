"""Small calculator-free diagnostics for completed VCNEB paths."""

from __future__ import annotations

import numpy as np

from ase import Atoms

from .core import cell_matrix, deformation_from_cell


def path_reaction_coordinate(
    images: list[Atoms] | tuple[Atoms, ...],
    *,
    cell_scale: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return cumulative VCNEB metric coordinates and segment lengths.

    Atomic coordinates are unwrapped in fractional space with the reference
    endpoint cell; the cell block uses ``cell_scale * (F-I)``.  The function is
    intentionally calculator-free and therefore works on an incomplete path,
    but it does not infer atom mappings between different chemical orders.
    """

    if len(images) < 2:
        raise ValueError("at least two images are required")
    reference = images[0]
    ref_cell = cell_matrix(reference)
    if np.linalg.det(ref_cell) <= 0.0:
        raise ValueError("the reference image must have a positive-volume cell")
    if cell_scale is None:
        cell_scale = abs(float(np.linalg.det(ref_cell))) ** (1.0 / 3.0)
    if not np.isfinite(cell_scale) or cell_scale <= 0.0:
        raise ValueError("cell_scale must be finite and positive")

    previous_q = reference.get_scaled_positions(wrap=False).copy()
    vectors = []
    for image_index, image in enumerate(images):
        if len(image) != len(reference):
            raise ValueError(f"image {image_index} atom count differs from reference")
        if image.get_chemical_symbols() != reference.get_chemical_symbols():
            raise ValueError(f"image {image_index} atom order differs from reference")
        q = image.get_scaled_positions(wrap=False).copy()
        if image_index:
            delta = q - previous_q
            for axis, periodic in enumerate(np.asarray(reference.pbc, dtype=bool)):
                if periodic:
                    delta[:, axis] -= np.rint(delta[:, axis])
            q = previous_q + delta
        previous_q = q
        deform = deformation_from_cell(cell_matrix(image), ref_cell)
        vectors.append(
            np.concatenate([(q @ ref_cell).reshape(-1), ((deform - np.eye(3)) * cell_scale).reshape(-1)])
        )

    values = np.asarray(vectors, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("path contains non-finite coordinates")
    segments = np.linalg.norm(np.diff(values, axis=0), axis=1)
    return np.concatenate([[0.0], np.cumsum(segments)]), segments


def mode_contribution_fractions(coordinates: np.ndarray) -> np.ndarray:
    """Normalize squared modal amplitudes independently at every image."""

    values = np.asarray(coordinates, dtype=float)
    if values.ndim != 2 or not values.size or not np.all(np.isfinite(values)):
        raise ValueError("coordinates must be a finite, non-empty 2D array")
    squared = values * values
    denominator = squared.sum(axis=1, keepdims=True)
    fractions = np.zeros_like(squared)
    np.divide(squared, denominator, out=fractions, where=denominator > 0.0)
    return fractions


def dominant_mode_indices(coordinates: np.ndarray, *, top: int = 3) -> list[list[int]]:
    """Return the mode indices ranked by absolute amplitude for each image."""

    if top < 1:
        raise ValueError("top must be positive")
    values = np.asarray(coordinates, dtype=float)
    if values.ndim != 2 or not np.all(np.isfinite(values)):
        raise ValueError("coordinates must be a finite 2D array")
    limit = min(int(top), values.shape[1])
    return [np.argsort(-np.abs(row), kind="stable")[:limit].astype(int).tolist() for row in values]


__all__ = [
    "dominant_mode_indices",
    "mode_contribution_fractions",
    "path_reaction_coordinate",
]
