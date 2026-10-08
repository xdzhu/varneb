"""Calculator-free path projections preserving a declared periodic lift.

Only the *initial* image chooses integer representatives relative to a
reference. The same shifts are applied to every subsequent image; no
image-wise minimum-image folding, atom mapping, or cell rotation is performed.
The resulting amplitudes describe geometry, not energy contributions.
"""

from __future__ import annotations

import numpy as np

from .core import deformation_from_cell
from .periodic_path import validate_periodic_path_lift


def continuous_reference_coordinates(images, reference):
    """Return non-affine atomic displacements and full finite cell strains.

    ``u=(q+n-q_ref)@H_ref`` is measured in Angstrom in the fixed reference
    metric. ``n`` is fixed at image zero. With ASE row cells,
    ``H=H_ref@F.T``; Green strain is ``(F.T@F-I)/2``. Atomic and cell
    coordinates are separate. The short-adjacent-step lift convention must
    already be valid; it cannot infer an undersampled winding path.
    """
    if (not len(reference) or not np.isfinite(reference.cell.array).all()
            or np.linalg.det(reference.cell.array) <= 1e-12):
        raise ValueError("nonempty reference with a finite right-handed cell required")
    lift = validate_periodic_path_lift(images)
    if (images[0].get_chemical_symbols() != reference.get_chemical_symbols()
            or not np.array_equal(images[0].pbc, reference.pbc)):
        raise ValueError("reference and path require identical ordered atoms and PBC")
    q_ref = reference.get_scaled_positions(wrap=False)
    q = np.asarray([image.get_scaled_positions(wrap=False) for image in images])
    if not np.isfinite(q).all() or not np.isfinite(q_ref).all():
        raise ValueError("reference/path coordinates must be finite")
    pbc = np.asarray(reference.pbc, dtype=bool)
    initial = q[0] - q_ref
    reduced = initial[:, pbc] - np.rint(initial[:, pbc])
    if np.any(np.isclose(np.abs(reduced), .5, rtol=0, atol=1e-10)):
        raise ValueError("initial half-cell reference image is ambiguous; register its gauge explicitly")
    shifts = np.zeros_like(q_ref, dtype=np.int64)
    shifts[:, pbc] = -np.rint(initial[:, pbc]).astype(np.int64)
    delta = q + shifts[None, :, :] - q_ref[None, :, :]
    deformation = np.asarray([deformation_from_cell(a.cell.array, reference.cell.array) for a in images])
    green = .5 * (np.einsum("nji,njk->nik", deformation, deformation) - np.eye(3))
    return {
        "displacements_A": delta @ reference.cell.array,
        "fractional_offsets": delta,
        "fixed_integer_lattice_shifts_by_atom": shifts,
        "deformation_gradients": deformation,
        "green_strains": green,
        "periodic_lift_policy": lift["policy"],
    }


def project_reference_basis(displacements_A, basis_columns, *, metric_weights=None):
    """Project translation-free displacements in a declared positive metric.

    The columns must be orthonormal in the *weighted* flattened Cartesian
    space. Unit weights give Angstrom amplitudes; masses in amu give
    sqrt(amu)-Angstrom amplitudes (e.g. Gamma dynamical-matrix eigenvectors).
    Degenerate eigenvector coordinates depend on their arbitrary rotation;
    squared norms summed over a complete degenerate subspace do not.
    """
    u = np.asarray(displacements_A, dtype=float)
    if u.ndim != 3 or u.shape[2] != 3 or not u.shape[0] or not u.shape[1] or not np.isfinite(u).all():
        raise ValueError("finite (n_images, n_atoms, 3) displacements required")
    n_atoms = u.shape[1]
    weights = np.ones(n_atoms) if metric_weights is None else np.asarray(metric_weights, dtype=float)
    if weights.shape != (n_atoms,) or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("metric_weights require one positive finite value per atom")
    basis = np.asarray(basis_columns, dtype=float)
    if (basis.ndim != 2 or basis.shape[0] != 3 * n_atoms or not basis.shape[1]
            or not np.isfinite(basis).all()
            or not np.allclose(basis.T @ basis, np.eye(basis.shape[1]), rtol=0, atol=1e-8)):
        raise ValueError("finite orthonormal weighted basis columns required")
    translation = np.average(u, axis=1, weights=weights)
    centered = u - translation[:, None, :]
    weighted = (centered * np.sqrt(weights)[None, :, None]).reshape(len(u), -1)
    amplitudes = weighted @ basis
    residual = weighted - amplitudes @ basis.T
    norm_squared = np.sum(weighted**2, axis=1)
    captured = np.zeros(len(u))
    np.divide(np.sum(amplitudes**2, axis=1), norm_squared, out=captured, where=norm_squared > 1e-24)
    return {"amplitudes": amplitudes, "residual_norm": np.linalg.norm(residual, axis=1),
            "total_norm": np.sqrt(norm_squared), "captured_squared_norm_fraction": np.clip(captured, 0, 1),
            "removed_translations_A": translation,
            "zero_displacement": norm_squared <= 1e-24}
