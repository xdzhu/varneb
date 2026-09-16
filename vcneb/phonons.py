"""Calculator-independent Gamma-point normal-mode analysis.

This module deliberately analyzes already available force constants.  It does
not launch a calculator or claim that a non-stationary NEB image has phonons.
For a stationary reference structure it converts a Cartesian force-constant
matrix into mass-weighted Gamma modes and projects a separately aligned path
onto those normal coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


Array = np.ndarray
_FREQUENCY_FACTOR_CM1 = 521.4708986  # sqrt(eV / (Angstrom^2 amu)) -> cm^-1


def _as_force_constant_matrix(force_constants: Array, n_atoms: int) -> Array:
    matrix = np.asarray(force_constants, dtype=float)
    width = 3 * n_atoms
    if matrix.shape == (n_atoms, 3, n_atoms, 3):
        matrix = matrix.reshape(width, width)
    if matrix.shape != (width, width):
        raise ValueError(
            "force_constants must have shape (3N, 3N) or (N, 3, N, 3)"
        )
    if not np.all(np.isfinite(matrix)):
        raise ValueError("force_constants contains non-finite values")
    return 0.5 * (matrix + matrix.T)


def _validated_masses(masses: Array) -> Array:
    values = np.asarray(masses, dtype=float).reshape(-1)
    if not len(values) or not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("masses must contain one finite positive value per atom")
    return values


def translation_basis(masses: Array) -> Array:
    """Return orthonormal mass-weighted Cartesian translation columns."""

    values = _validated_masses(masses)
    basis = np.zeros((3 * len(values), 3), dtype=float)
    for axis in range(3):
        basis[axis::3, axis] = np.sqrt(values)
    basis /= np.linalg.norm(basis, axis=0, keepdims=True)
    return basis


def mass_weighted_dynamical_matrix(force_constants: Array, masses: Array) -> Array:
    """Build the symmetric Gamma dynamical matrix in eV/(A^2 amu)."""

    values = _validated_masses(masses)
    matrix = _as_force_constant_matrix(force_constants, len(values))
    weights = np.repeat(np.sqrt(values), 3)
    return matrix / np.outer(weights, weights)


@dataclass(frozen=True)
class GammaModes:
    """Mass-weighted Gamma modes for one stationary reference configuration.

    ``eigenvectors`` are unit vectors in mass-weighted Cartesian coordinates;
    the signed frequency is negative for a negative Hessian eigenvalue.  Thus
    a negative value represents an imaginary harmonic frequency but does not,
    by itself, prove a NEB path image is a first-order saddle.
    """

    masses_amu: Array
    eigenvalues_eV_per_A2_amu: Array
    frequencies_cm1: Array
    eigenvectors: Array
    translations_projected: bool

    def __post_init__(self) -> None:
        masses = _validated_masses(self.masses_amu)
        width = 3 * len(masses)
        values = np.asarray(self.eigenvalues_eV_per_A2_amu, dtype=float).reshape(-1)
        frequencies = np.asarray(self.frequencies_cm1, dtype=float).reshape(-1)
        vectors = np.asarray(self.eigenvectors, dtype=float)
        if values.size != width or frequencies.size != width or vectors.shape != (width, width):
            raise ValueError("GammaModes arrays must have dimensions (3N,) and (3N, 3N)")
        if not all(np.all(np.isfinite(item)) for item in (values, frequencies, vectors)):
            raise ValueError("GammaModes contains non-finite values")
        object.__setattr__(self, "masses_amu", masses.copy())
        object.__setattr__(self, "eigenvalues_eV_per_A2_amu", values.copy())
        object.__setattr__(self, "frequencies_cm1", frequencies.copy())
        object.__setattr__(self, "eigenvectors", vectors.copy())

    @property
    def n_atoms(self) -> int:
        return len(self.masses_amu)

    def cartesian_eigenvectors(self) -> Array:
        """Return displacement eigenvectors normalized in the mass metric."""

        weights = np.repeat(np.sqrt(self.masses_amu), 3)
        return self.eigenvectors / weights[:, None]


def diagonalize_gamma_modes(
    force_constants: Array,
    masses: Array,
    *,
    project_translations: bool = True,
) -> GammaModes:
    """Diagonalize a mass-weighted Gamma dynamical matrix.

    Translation projection removes numerical acoustic contamination before
    diagonalization.  It intentionally leaves three near-zero eigenvalues in
    the returned complete basis so downstream arrays always have size ``3N``.
    """

    values = _validated_masses(masses)
    matrix = mass_weighted_dynamical_matrix(force_constants, values)
    if project_translations:
        translations = translation_basis(values)
        projector = np.eye(3 * len(values)) - translations @ translations.T
        matrix = projector @ matrix @ projector
        matrix = 0.5 * (matrix + matrix.T)
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * _FREQUENCY_FACTOR_CM1
    return GammaModes(
        masses_amu=values,
        eigenvalues_eV_per_A2_amu=eigenvalues,
        frequencies_cm1=frequencies,
        eigenvectors=eigenvectors,
        translations_projected=project_translations,
    )


def project_displacements_onto_gamma_modes(displacements: Array, modes: GammaModes) -> Array:
    """Project aligned Cartesian displacements onto mass-weighted modes.

    ``displacements`` has shape ``(N, 3)`` or ``(n_frames, N, 3)`` and must
    already use a common reference cell, atom mapping, and periodic gauge.
    The result has shape ``(3N,)`` or ``(n_frames, 3N)`` in
    ``sqrt(amu) * Angstrom`` normal-coordinate units.
    """

    values = np.asarray(displacements, dtype=float)
    single = values.ndim == 2
    if single:
        values = values[None, :, :]
    if values.ndim != 3 or values.shape[1:] != (modes.n_atoms, 3):
        raise ValueError("displacements must have shape (N, 3) or (n_frames, N, 3)")
    if not np.all(np.isfinite(values)):
        raise ValueError("displacements contains non-finite values")
    weighted = (values * np.sqrt(modes.masses_amu)[None, :, None]).reshape(len(values), -1)
    coordinates = weighted @ modes.eigenvectors
    return coordinates[0] if single else coordinates


def tangent_mode_overlaps(coordinates: Array, *, tolerance: float = 1e-14) -> Array:
    """Return signed overlap of each path segment tangent with every mode."""

    values = np.asarray(coordinates, dtype=float)
    if values.ndim != 2 or values.shape[0] < 2:
        raise ValueError("coordinates must have shape (n_frames >= 2, n_modes)")
    if not np.all(np.isfinite(values)) or tolerance <= 0.0:
        raise ValueError("coordinates must be finite and tolerance positive")
    tangents = np.diff(values, axis=0)
    norms = np.linalg.norm(tangents, axis=1)
    result = np.zeros_like(tangents)
    active = norms > tolerance
    result[active] = tangents[active] / norms[active, None]
    return result


def load_gamma_force_constants(path: str | Path) -> tuple[Array, Array]:
    """Load standard VARNEB ``.npz`` force constants and atomic masses.

    The archive must contain ``force_constants`` and ``masses_amu`` arrays;
    calculator, supercell, displacement and structure provenance stay as
    additional metadata and are deliberately not guessed here.
    """

    with np.load(Path(path), allow_pickle=False) as data:
        if "force_constants" not in data or "masses_amu" not in data:
            raise ValueError("force-constant archive requires force_constants and masses_amu arrays")
        return np.asarray(data["force_constants"], dtype=float), np.asarray(data["masses_amu"], dtype=float)


__all__ = [
    "GammaModes",
    "diagonalize_gamma_modes",
    "load_gamma_force_constants",
    "mass_weighted_dynamical_matrix",
    "project_displacements_onto_gamma_modes",
    "tangent_mode_overlaps",
    "translation_basis",
]
