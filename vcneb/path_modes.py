"""Calculator-free empirical atomic directions fitted to an aligned path.

These singular vectors summarize sampled displacements. They are *not* phonon
eigenvectors, harmonic frequencies, transition-state unstable modes, or a
predictive model of the off-path potential-energy surface.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


Array = np.ndarray


def _weighted_path(
    displacements_A: Array,
    masses_amu: Array,
    *,
    remove_translations: bool,
) -> tuple[Array, Array, Array]:
    displacements = np.asarray(displacements_A, dtype=float)
    masses = np.asarray(masses_amu, dtype=float).reshape(-1)
    if (displacements.ndim != 3 or displacements.shape[2] != 3
            or displacements.shape[0] < 1 or displacements.shape[1] != len(masses)
            or not np.all(np.isfinite(displacements))):
        raise ValueError("displacements must be a finite (n_images, n_atoms, 3) array")
    if not len(masses) or not np.all(np.isfinite(masses)) or np.any(masses <= 0.0):
        raise ValueError("masses must contain one positive finite value per atom")
    if remove_translations and len(masses) < 2:
        raise ValueError("translation-free atomic modes require at least two atoms")
    if remove_translations:
        translations = np.average(displacements, axis=1, weights=masses)
        centered = displacements - translations[:, None, :]
    else:
        translations = np.zeros((len(displacements), 3))
        centered = displacements.copy()
    weighted = (centered * np.sqrt(masses)[None, :, None]).reshape(len(displacements), -1)
    return weighted, translations, masses


@dataclass(frozen=True)
class PathAtomicSVD:
    """Mass-metric orthonormal path directions and fit diagnostics.

    ``basis_mass_weighted`` has unit-normalized columns. ``coefficients`` and
    residuals have units ``sqrt(amu) * Angstrom``. The reference atom mapping,
    periodic unwrapping, and cell gauge belong to the caller and must be
    audited separately.
    """

    masses_amu: Array
    basis_mass_weighted: Array
    coefficients_sqrt_amu_A: Array
    singular_values_sqrt_amu_A: Array
    residuals_sqrt_amu_A: Array
    removed_translations_A: Array
    captured_squared_norm_fraction: float
    remove_translations: bool

    @property
    def n_components(self) -> int:
        return self.basis_mass_weighted.shape[1]

    def project(self, displacements_A: Array) -> tuple[Array, Array]:
        """Project held-out aligned displacements; return Q and residual norm."""

        weighted, _, masses = _weighted_path(
            displacements_A, self.masses_amu,
            remove_translations=self.remove_translations,
        )
        if len(masses) != len(self.masses_amu):
            raise ValueError("atom count differs from fitted path")
        coefficients = weighted @ self.basis_mass_weighted
        residuals = np.linalg.norm(
            weighted - coefficients @ self.basis_mass_weighted.T, axis=1
        )
        return coefficients, residuals


def fit_path_atomic_svd(
    displacements_A: Array,
    masses_amu: Array,
    *,
    n_components: int,
    remove_translations: bool = True,
    rank_tolerance: float = 1e-10,
) -> PathAtomicSVD:
    """Fit a declared number of empirical path axes in the atomic mass metric.

    Fitting uses the same path for discovery, so captured norm is descriptive,
    not independent validation. Use :func:`leave_one_image_out_residuals` for a
    modest interpolation check; neither statistic identifies phonon modes.
    """

    weighted, translations, masses = _weighted_path(
        displacements_A, masses_amu, remove_translations=remove_translations
    )
    if weighted.shape[0] < 2:
        raise ValueError("fitting requires at least two path images")
    max_components = min(
        weighted.shape[0], weighted.shape[1] - (3 if remove_translations else 0)
    )
    if not isinstance(n_components, int) or not 1 <= n_components <= max_components:
        raise ValueError("n_components exceeds images or atomic internal dimensions")
    if not np.isfinite(rank_tolerance) or rank_tolerance <= 0.0:
        raise ValueError("rank_tolerance must be finite and positive")
    _, singular_values, vectors_t = np.linalg.svd(weighted, full_matrices=False)
    if singular_values[n_components - 1] <= rank_tolerance * singular_values[0]:
        raise ValueError("requested empirical path axis is numerically rank deficient")
    basis = vectors_t[:n_components].T.copy()
    for column in range(n_components):
        pivot = int(np.argmax(np.abs(basis[:, column])))
        if basis[pivot, column] < 0.0:
            basis[:, column] *= -1.0
    coefficients = weighted @ basis
    residuals = np.linalg.norm(weighted - coefficients @ basis.T, axis=1)
    total = float(np.dot(singular_values, singular_values))
    captured = float(np.dot(singular_values[:n_components], singular_values[:n_components])) / total
    return PathAtomicSVD(
        masses_amu=masses.copy(),
        basis_mass_weighted=basis,
        coefficients_sqrt_amu_A=coefficients,
        singular_values_sqrt_amu_A=singular_values.copy(),
        residuals_sqrt_amu_A=residuals,
        removed_translations_A=translations,
        captured_squared_norm_fraction=captured,
        remove_translations=remove_translations,
    )


def leave_one_image_out_residuals(
    displacements_A: Array,
    masses_amu: Array,
    *,
    n_components: int,
    remove_translations: bool = True,
) -> Array:
    """Refit without each image and measure its held-out atomic residual."""

    displacements = np.asarray(displacements_A, dtype=float)
    if displacements.ndim != 3 or len(displacements) < 3:
        raise ValueError("leave-one-image-out validation requires at least three images")
    residuals = np.empty(len(displacements), dtype=float)
    for index in range(len(displacements)):
        training = np.delete(displacements, index, axis=0)
        fit = fit_path_atomic_svd(
            training, masses_amu, n_components=n_components,
            remove_translations=remove_translations,
        )
        _, held_out = fit.project(displacements[index:index + 1])
        residuals[index] = held_out[0]
    return residuals


__all__ = ["PathAtomicSVD", "fit_path_atomic_svd", "leave_one_image_out_residuals"]
