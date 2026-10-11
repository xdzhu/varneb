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


def project_reference_atomic_gradient(forces_eV_A, current_cell_A, reference_cell_A,
                                      basis_columns, *, masses_amu):
    """Pull back a fixed-cell energy covector, then resolve its internal span.

    For row cells, reference displacement ``u`` produces current displacement
    ``u @ inv(H_ref) @ H``. Thus ``g_q = -F @ (inv(H_ref) @ H).T / sqrt(m)``
    is conjugate to ``q=sqrt(m)*u``, in eV/(sqrt(amu)*Angstrom). Applying the
    displacement transform to the force instead would violate virtual work.

    Rigid translation is an explicitly removed coordinate gauge, not an ASR
    correction to source modes or forces. Restrict the original weighted basis
    to that gauge and require its rank to survive. Return its measured source
    translation overlap and the raw gradient's translation component as well.
    Scalar span norms are invariant to rotations within the basis. They are
    neither energy fractions nor ordinary NEB convergence norms; the cell
    gradient and mechanically fixed-coordinate reactions are not included.
    """
    forces = np.asarray(forces_eV_A, dtype=float)
    masses = np.asarray(masses_amu, dtype=float)
    if (forces.ndim != 2 or forces.shape[1] != 3 or not forces.shape[0]
            or not np.isfinite(forces).all() or masses.shape != (len(forces),)
            or not np.isfinite(masses).all() or np.any(masses <= 0)):
        raise ValueError("finite atomic forces in eV/A and positive ordered masses in amu required")
    current, reference = [np.asarray(h, dtype=float) for h in (current_cell_A, reference_cell_A)]
    for cell in (current, reference):
        if (cell.shape != (3, 3) or not np.isfinite(cell).all()
                or not np.isfinite(np.linalg.det(cell)) or np.linalg.det(cell) <= 1e-12
                or np.linalg.cond(cell) > 1e12):
            raise ValueError("finite well-conditioned right-handed cells in Angstrom required")
    basis = np.asarray(basis_columns, dtype=float)
    if (basis.ndim != 2 or basis.shape[0] != forces.size or not basis.shape[1]
            or not np.isfinite(basis).all()
            or not np.allclose(basis.T @ basis, np.eye(basis.shape[1]), rtol=0, atol=1e-8)):
        raise ValueError("finite orthonormal mass-weighted reference basis required")
    roots = np.sqrt(masses)
    normalized_roots = roots / roots.max()
    normalized_roots /= np.linalg.norm(normalized_roots)
    translations = np.zeros((forces.size, 3))
    for axis in range(3):
        translations[axis::3, axis] = normalized_roots
    source_overlap = translations.T @ basis
    restricted = basis - translations @ source_overlap
    allowed_basis, singular, _ = np.linalg.svd(restricted, full_matrices=False)
    if np.count_nonzero(singular > 1e-10) != basis.shape[1]:
        raise ValueError("reference rank lost on the translation-free chart; no silent mode removal")
    covector = (-forces @ np.linalg.solve(reference, current).T / roots[:, None]).ravel()
    rigid = translations @ (translations.T @ covector)
    internal = covector - rigid
    retained = allowed_basis @ (allowed_basis.T @ internal)
    omitted = internal - retained
    squared = float(np.dot(internal, internal))
    return {
        "gradient_unit": "eV/(sqrt(amu)*Angstrom)",
        "source_basis_rank": basis.shape[1], "translation_free_basis_rank": len(singular),
        "source_translation_overlap_Frobenius_norm": float(np.linalg.norm(source_overlap)),
        "restricted_basis_min_singular_value": float(singular.min()),
        "source_modes_or_forces_ASR_corrected": False,
        "translation_free_chart_restriction_applied": True,
        "raw_reference_covector_norm": float(np.linalg.norm(covector)),
        "rigid_translation_gradient_norm": float(np.linalg.norm(rigid)),
        "internal_gradient_norm": float(np.sqrt(squared)),
        "retained_gradient_norm": float(np.linalg.norm(retained)),
        "omitted_gradient_norm": float(np.linalg.norm(omitted)),
        "retained_squared_gradient_fraction": float(np.dot(retained, retained)/squared) if squared > 1e-24 else None,
        "zero_internal_gradient": squared <= 1e-24,
        "norm_identity_defect": float(np.dot(retained, retained)+np.dot(omitted, omitted)-squared),
    }
