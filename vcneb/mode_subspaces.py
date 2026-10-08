"""Explicit real subspaces of complex reference directions, without ASR.

The result is a geometric span, not a set of new phonon eigenpairs. Complex
phases and real mixing are allowed, but an unexpected real rank or excessive
projection error is rejected rather than silently discarding directions.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RealModeSubspace:
    basis_columns: np.ndarray
    real_singular_values: np.ndarray
    complex_column_norms: np.ndarray
    relative_column_projection_errors: np.ndarray
    orthogonality_max_defect: float
    expected_rank: int
    relative_rank_tolerance: float


def real_mode_subspace(columns, *, expected_rank, relative_rank_tolerance=5e-5,
                       max_relative_projection_error=5e-5):
    """Represent a declared complex geometric span by real orthonormal columns.

    Input columns already use the chosen flattened coordinate/metric convention.
    SVD of [Re(V), Im(V)] makes the phase policy explicit. Only singular values
    above the declared relative threshold count toward rank; exactly the
    requested rank must be present. Every original complex column must be
    reconstructed within the stated relative error. No mass conversion,
    translation removal, ASR, symmetry operation or eigensystem is inferred.
    Basis coordinates may rotate under phases or degenerate singular values;
    compare its projector and total subspace weight, not individual SVD axes.
    """
    values = np.array(columns, dtype=complex, copy=True)
    if (values.ndim != 2 or not 1 <= values.shape[0] <= 384
            or not 1 <= values.shape[1] <= values.shape[0]
            or not np.isfinite(values).all()):
        raise ValueError("finite bounded complex columns required")
    if (isinstance(expected_rank, (bool, np.bool_))
            or not isinstance(expected_rank, (int, np.integer))
            or not 1 <= expected_rank <= min(values.shape[0], 2 * values.shape[1])):
        raise ValueError("positive feasible integer expected rank required")
    for tolerance in (relative_rank_tolerance, max_relative_projection_error):
        if not np.isfinite(tolerance) or not 0 < tolerance < 1:
            raise ValueError("finite rank/projection tolerances within (0,1) required")
    # Independent column scales are immaterial to a span. Scaling first avoids
    # overflow in norms while retaining each phase and direction.
    scales = np.max(np.abs(values), axis=0)
    if np.any(scales <= 0) or not np.isfinite(scales).all():
        raise ValueError("nonzero numerically resolvable columns required")
    scaled = values / scales
    scaled_norms = np.linalg.norm(scaled, axis=0)
    norms = scales * scaled_norms
    if not np.isfinite(norms).all() or np.any(norms <= 0):
        raise ValueError("column norms exceed numerical bounds")
    normalized = scaled / scaled_norms
    basis, singular, _ = np.linalg.svd(np.c_[normalized.real, normalized.imag],
                                      full_matrices=False)
    rank = int(np.count_nonzero(singular > relative_rank_tolerance * singular[0]))
    if rank != expected_rank:
        raise ValueError(f"real rank {rank} differs from declared rank {expected_rank}")
    basis = basis[:, :rank]
    errors = np.linalg.norm(normalized - basis @ (basis.T @ normalized), axis=0)
    if np.max(errors) > max_relative_projection_error:
        raise ValueError("real-span projection error exceeds declared precision")
    defect = float(np.max(np.abs(basis.T @ basis - np.eye(rank))))
    for array in (basis, singular, norms, errors):
        array.flags.writeable = False
    return RealModeSubspace(basis, singular, norms, errors, defect,
                            int(expected_rank), float(relative_rank_tolerance))


def mass_translation_overlaps(basis_columns, masses_amu):
    """Measure, but do not remove, rigid-translation overlap in a mass metric."""
    basis = np.asarray(basis_columns)
    masses = np.asarray(masses_amu, dtype=float)
    if (masses.ndim != 1 or not 1 <= len(masses) <= 128
            or not np.isfinite(masses).all() or np.any(masses <= 0)
            or basis.ndim != 2 or basis.shape[0] != 3 * len(masses)
            or not basis.shape[1] or not np.isfinite(basis).all()):
        raise ValueError("finite compatible basis and positive ordered masses required")
    roots = np.sqrt(masses)
    roots /= roots.max()
    translations = np.zeros((3 * len(masses), 3))
    for direction in range(3):
        translations[direction::3, direction] = roots
    translations /= np.linalg.norm(translations, axis=0)
    return np.linalg.norm(translations.T @ basis, axis=0)
