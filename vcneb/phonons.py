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
_BOHR_TO_ANGSTROM = 0.529177210903
_THZ_TO_CM1 = 33.35640951981521


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


def force_constants_to_eV_per_A2(force_constants: Array, *, unit: str = "eV/angstrom^2") -> Array:
    """Convert a documented force-constant unit to ``eV / Angstrom^2``.

    Phonopy's ABACUS interface stores force constants as
    ``eV / (Angstrom * au)`` because ABACUS structures use Bohr lattice
    units. Treating those values as ``eV / Angstrom^2`` underestimates every
    harmonic frequency by ``sqrt(Bohr)``. No unit is guessed: callers must
    retain the documented source convention.
    """

    values = np.asarray(force_constants, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("force_constants contains non-finite values")
    normalized = unit.strip().lower().replace("ångström", "angstrom").replace("å", "angstrom")
    if normalized in {"ev/angstrom^2", "ev/a^2"}:
        return values.copy()
    if normalized in {"ev/angstrom.au", "ev/(angstrom*au)", "ev/a.au"}:
        return values / _BOHR_TO_ANGSTROM
    raise ValueError(
        "unsupported force-constant unit; use eV/angstrom^2 or eV/angstrom.au"
    )


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


def anchored_gamma_axis_weights(
    modes: GammaModes,
    mode_indices: Array,
    cartesian_anchors: Array,
    *,
    tolerance: float = 1e-10,
) -> Array:
    """Build two signed axes from physical anchors inside a mode subspace.

    The two columns of ``cartesian_anchors`` are atomic displacement patterns
    in Angstrom. Each is projected onto the selected *mass-weighted* mode
    subspace, then orthogonalized in the given order. The resulting physical
    axes are invariant to arbitrary rotations/sign changes of an exactly
    degenerate eigenvector basis. The anchors fix orientation; individual
    mode numbers do not. A vanishing or dependent projected anchor is an
    error rather than an arbitrary axis choice.
    """

    indices = np.asarray(mode_indices)
    anchors = np.asarray(cartesian_anchors, dtype=float)
    count = modes.eigenvectors.shape[1]
    if indices.ndim != 1 or len(indices) < 2 or not np.issubdtype(indices.dtype, np.integer):
        raise ValueError("mode_indices must contain at least two integer indices")
    if np.any(indices < 0) or np.any(indices >= count) or len(np.unique(indices)) != len(indices):
        raise ValueError("mode_indices must be unique and within the Gamma spectrum")
    if anchors.shape != (count, 2) or not np.all(np.isfinite(anchors)):
        raise ValueError("cartesian_anchors must contain two finite atomic displacement columns")
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be finite and positive")
    sqrt_masses = np.repeat(np.sqrt(modes.masses_amu), 3)
    subspace = modes.eigenvectors[:, indices]
    if not np.allclose(subspace.T @ subspace, np.eye(len(indices)), rtol=1e-8, atol=1e-8):
        raise ValueError("selected Gamma eigenvectors must be orthonormal")
    coefficients = subspace.T @ (sqrt_masses[:, None] * anchors)
    if np.any(np.linalg.norm(coefficients, axis=0) <= tolerance):
        raise ValueError("an anchor has negligible projection onto the selected mode subspace")
    first = coefficients[:, 0] / np.linalg.norm(coefficients[:, 0])
    second = coefficients[:, 1] - first * np.dot(first, coefficients[:, 1])
    if np.linalg.norm(second) <= tolerance * np.linalg.norm(coefficients[:, 1]):
        raise ValueError("projected anchors are linearly dependent")
    second /= np.linalg.norm(second)
    weights = np.zeros((count, 2))
    weights[indices, 0] = first
    weights[indices, 1] = second
    return weights


@dataclass(frozen=True)
class PhonopyGammaEigenpairs:
    """Raw, auditable Gamma-point eigenpairs obtained from Phonopy.

    The eigenvectors are Phonopy's mass-weighted dynamical-matrix
    eigenvectors, with their arbitrary complex phase fixed deterministically
    at Gamma. They are therefore directly suitable for comparison with the
    VARNEB mass-weighted normal-coordinate convention. Degenerate-mode
    rotations remain physically arbitrary and must be analyzed as subspaces.
    """

    frequencies_thz: Array
    eigenvectors_mass_weighted: Array

    def __post_init__(self) -> None:
        frequencies = np.asarray(self.frequencies_thz, dtype=float).reshape(-1)
        vectors = np.asarray(self.eigenvectors_mass_weighted, dtype=float)
        if not len(frequencies) or vectors.shape != (len(frequencies), len(frequencies)):
            raise ValueError("Phonopy Gamma eigenpairs must contain a square eigenvector matrix")
        if not np.all(np.isfinite(frequencies)) or not np.all(np.isfinite(vectors)):
            raise ValueError("Phonopy Gamma eigenpairs contain non-finite values")
        if not np.allclose(vectors.T @ vectors, np.eye(len(frequencies)), rtol=1e-8, atol=1e-8):
            raise ValueError("Phonopy Gamma eigenvectors are not orthonormal")
        object.__setattr__(self, "frequencies_thz", frequencies.copy())
        object.__setattr__(self, "eigenvectors_mass_weighted", vectors.copy())


def phonopy_gamma_eigenpairs(phonon: object, *, imaginary_tolerance: float = 1e-10) -> PhonopyGammaEigenpairs:
    """Obtain phase-fixed real Gamma eigenvectors directly from a Phonopy object.

    ``phonon`` is deliberately duck-typed so VARNEB remains calculator and
    Phonopy-installation independent at import time. It must provide the
    standard ``run_qpoints`` and ``get_qpoints_dict`` methods. At exact Gamma
    without a directional NAC request the dynamical matrix is real; a material
    imaginary component after deterministic phase fixing is rejected rather
    than silently discarded.
    """

    if imaginary_tolerance <= 0.0 or not np.isfinite(imaginary_tolerance):
        raise ValueError("imaginary_tolerance must be finite and positive")
    run_qpoints = getattr(phonon, "run_qpoints", None)
    get_qpoints_dict = getattr(phonon, "get_qpoints_dict", None)
    if not callable(run_qpoints) or not callable(get_qpoints_dict):
        raise TypeError("phonon must provide Phonopy run_qpoints and get_qpoints_dict methods")
    run_qpoints([[0.0, 0.0, 0.0]], with_eigenvectors=True)
    result = get_qpoints_dict()
    try:
        frequencies = np.asarray(result["frequencies"], dtype=float)
        vectors = np.asarray(result["eigenvectors"], dtype=complex)
    except (KeyError, TypeError) as exc:
        raise ValueError("Phonopy Gamma result requires frequencies and eigenvectors") from exc
    if frequencies.shape[0] != 1 or vectors.shape[0] != 1:
        raise ValueError("Phonopy Gamma query returned an unexpected q-point shape")
    frequencies = frequencies[0]
    vectors = vectors[0]
    if vectors.shape != (len(frequencies), len(frequencies)):
        raise ValueError("Phonopy Gamma eigenvectors have an unexpected shape")

    # Eigenvectors have an arbitrary U(1) phase. Fix it using the largest
    # component of each vector; for a real Gamma dynamical matrix this makes
    # each vector real up to numerical noise without choosing a direction in a
    # degenerate subspace.
    phase_fixed = vectors.copy()
    for mode in range(phase_fixed.shape[1]):
        pivot = int(np.argmax(np.abs(phase_fixed[:, mode])))
        phase_fixed[:, mode] *= np.exp(-1j * np.angle(phase_fixed[pivot, mode]))
    scale = max(1.0, float(np.max(np.abs(phase_fixed.real))))
    if float(np.max(np.abs(phase_fixed.imag))) > imaginary_tolerance * scale:
        raise ValueError(
            "Phonopy Gamma eigenvectors retain a material imaginary component; "
            "use an exact non-directional Gamma calculation"
        )
    return PhonopyGammaEigenpairs(
        frequencies_thz=frequencies,
        eigenvectors_mass_weighted=phase_fixed.real,
    )


def save_phonopy_gamma_eigenpairs(path: str | Path, eigenpairs: PhonopyGammaEigenpairs) -> None:
    """Save directly generated Phonopy Gamma eigenpairs in a portable NPZ archive."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        target,
        frequencies_thz=eigenpairs.frequencies_thz,
        eigenvectors_mass_weighted=eigenpairs.eigenvectors_mass_weighted,
    )


def load_phonopy_gamma_eigenpairs(path: str | Path) -> PhonopyGammaEigenpairs:
    """Load an NPZ archive written by :func:`save_phonopy_gamma_eigenpairs`."""

    with np.load(Path(path), allow_pickle=False) as data:
        if "frequencies_thz" not in data or "eigenvectors_mass_weighted" not in data:
            raise ValueError("Phonopy Gamma archive requires frequencies_thz and eigenvectors_mass_weighted")
        return PhonopyGammaEigenpairs(
            frequencies_thz=np.asarray(data["frequencies_thz"], dtype=float),
            eigenvectors_mass_weighted=np.asarray(data["eigenvectors_mass_weighted"], dtype=float),
        )


def gamma_modes_from_phonopy_eigenpairs(
    eigenpairs: PhonopyGammaEigenpairs,
    masses: Array,
) -> GammaModes:
    """Convert directly exported Phonopy Gamma eigenpairs to ``GammaModes``.

    Phonopy reports frequencies in THz. Their conversion to signed harmonic
    eigenvalues is used only to preserve VARNEB's common report schema; the
    supplied Phonopy eigenvectors themselves are retained for projection.
    """

    values = _validated_masses(masses)
    width = 3 * len(values)
    if len(eigenpairs.frequencies_thz) != width:
        raise ValueError("Phonopy Gamma eigenpair count is incompatible with masses")
    frequencies_cm1 = eigenpairs.frequencies_thz * _THZ_TO_CM1
    eigenvalues = np.sign(frequencies_cm1) * (np.abs(frequencies_cm1) / _FREQUENCY_FACTOR_CM1) ** 2
    return GammaModes(
        masses_amu=values,
        eigenvalues_eV_per_A2_amu=eigenvalues,
        frequencies_cm1=frequencies_cm1,
        eigenvectors=eigenpairs.eigenvectors_mass_weighted,
        translations_projected=False,
    )


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
    "PhonopyGammaEigenpairs",
    "anchored_gamma_axis_weights",
    "diagonalize_gamma_modes",
    "force_constants_to_eV_per_A2",
    "gamma_modes_from_phonopy_eigenpairs",
    "load_gamma_force_constants",
    "load_phonopy_gamma_eigenpairs",
    "mass_weighted_dynamical_matrix",
    "phonopy_gamma_eigenpairs",
    "project_displacements_onto_gamma_modes",
    "save_phonopy_gamma_eigenpairs",
    "tangent_mode_overlaps",
    "translation_basis",
]
