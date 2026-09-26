"""Compare GaN's joint atomic--strain unstable direction with endpoint Gamma modes.

This is a *direction* decomposition, not an energy/barrier decomposition.  A
joint Hessian eigenvector uses Cartesian atomic displacements plus scaled cell
strain; only its atomic block can be compared with fixed-cell Gamma phonons.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def project_atomic_direction(
    atomic_cartesian: np.ndarray,
    masses_amu: np.ndarray,
    gamma_basis: np.ndarray,
    groups: list[list[int]],
) -> dict:
    """Remove *mass-weighted* common translation before optical projection.

    The joint Hessian removes equal-Cartesian translation.  With unlike Ga/N
    masses this is not the acoustic gauge of a mass-weighted phonon basis.
    """

    atomic = np.asarray(atomic_cartesian, dtype=float)
    masses = np.asarray(masses_amu, dtype=float)
    basis = np.asarray(gamma_basis, dtype=float)
    n_atoms = len(masses)
    if (atomic.shape != (n_atoms, 3) or basis.shape != (3 * n_atoms, 3 * n_atoms)
            or not np.all(np.isfinite(atomic)) or not np.all(np.isfinite(masses))
            or np.any(masses <= 0)):
        raise ValueError("invalid atomic direction, masses, or Gamma basis")
    if not np.allclose(basis.T @ basis, np.eye(3 * n_atoms), atol=1e-7):
        raise ValueError("Gamma eigenvectors are not orthonormal")
    optical_indices = [index for group in groups for index in group]
    if sorted(optical_indices) != list(range(3, 3 * n_atoms)):
        raise ValueError("Gamma groups must partition the optical modes")

    raw = (atomic * np.sqrt(masses)[:, None]).ravel()
    if np.linalg.norm(raw) <= 1e-12:
        raise ValueError("joint unstable direction has no atomic component")
    raw_q = basis.T @ (raw / np.linalg.norm(raw))
    center_of_mass = np.sum(masses[:, None] * atomic, axis=0) / np.sum(masses)
    optical = ((atomic - center_of_mass) * np.sqrt(masses)[:, None]).ravel()
    if np.linalg.norm(optical) <= 1e-12:
        raise ValueError("joint atomic direction is only a common translation")
    q = basis.T @ (optical / np.linalg.norm(optical))
    optical_norm = float(np.sum(q[3:] ** 2))
    fractions = [float(np.sum(q[group] ** 2) / optical_norm) for group in groups]
    return {
        "raw_acoustic_squared_fraction_before_com_removal": float(np.sum(raw_q[:3] ** 2)),
        "acoustic_squared_fraction_after_com_removal": float(np.sum(q[:3] ** 2)),
        "removed_center_of_mass_displacement_in_joint_direction_A": center_of_mass.tolist(),
        "group_fractions_of_mass_weighted_atomic_optical_direction": fractions,
    }


def _rotation_to_center(endpoint_cell: np.ndarray, center_cell: np.ndarray) -> np.ndarray:
    deformation = np.linalg.solve(endpoint_cell, center_cell)
    left, _, right = np.linalg.svd(deformation)
    rotation = left @ right
    if np.linalg.det(rotation) < 0:
        raise ValueError("endpoint-to-center cell mapping reverses handedness")
    return rotation


def block_spectrum(atomic: np.ndarray, strain: np.ndarray, coupling: np.ndarray) -> dict:
    """Audit whether a joint negative curvature survives either frozen block.

    The Schur complement is the local atomic curvature when the positive
    strain block is minimized at quadratic order. It is not a phonon spectrum.
    """

    atomic = np.asarray(atomic, dtype=float)
    strain = np.asarray(strain, dtype=float)
    coupling = np.asarray(coupling, dtype=float)
    if (atomic.ndim != 2 or strain.ndim != 2 or atomic.shape[0] != atomic.shape[1]
            or strain.shape[0] != strain.shape[1]
            or coupling.shape != (len(atomic), len(strain))
            or not all(np.isfinite(array).all() for array in (atomic, strain, coupling))
            or not np.allclose(atomic, atomic.T, atol=1e-8)
            or not np.allclose(strain, strain.T, atol=1e-8)):
        raise ValueError("invalid symmetric atomic/strain Hessian blocks")
    atomic_eigenvalues = np.linalg.eigvalsh(atomic)
    strain_eigenvalues = np.linalg.eigvalsh(strain)
    if min(atomic_eigenvalues[0], strain_eigenvalues[0]) <= 0:
        raise ValueError("Schur diagnostic requires positive frozen blocks")
    joint = np.block([[atomic, coupling], [coupling.T, strain]])
    joint_eigenvalues = np.linalg.eigvalsh(joint)
    relaxed_atomic = atomic - coupling @ np.linalg.solve(strain, coupling.T)
    return {
        "lowest_frozen_atomic_eV_per_A2": float(atomic_eigenvalues[0]),
        "lowest_frozen_strain_eV_per_A2": float(strain_eigenvalues[0]),
        "lowest_joint_two_eV_per_A2": joint_eigenvalues[:2].tolist(),
        "lowest_strain_relaxed_atomic_schur_eV_per_A2": float(
            np.linalg.eigvalsh(relaxed_atomic)[0]
        ),
        "atomic_strain_coupling_frobenius_eV_per_A2": float(np.linalg.norm(coupling)),
    }


def _translation_free_blocks(hessian: np.ndarray, n_atoms: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    translation = np.kron(np.ones((1, n_atoms)) / np.sqrt(n_atoms), np.eye(3))
    atomic_basis = np.linalg.svd(translation, full_matrices=True)[2][3:].T
    n_cartesian = 3 * n_atoms
    atomic = atomic_basis.T @ hessian[:n_cartesian, :n_cartesian] @ atomic_basis
    strain = hessian[n_cartesian:, n_cartesian:]
    coupling = atomic_basis.T @ hessian[:n_cartesian, n_cartesian:]
    return atomic, strain, coupling


def analyze(gamma_root: Path, joint_root: Path) -> dict:
    gamma_manifest_path = gamma_root / "manifest.json"
    gamma_audit_path = gamma_root / "audit.json"
    gamma_subspace_path = gamma_root / "subspace_audit_v2.json"
    gamma_manifest = json.loads(gamma_manifest_path.read_text(encoding="utf-8"))
    gamma_audit = json.loads(gamma_audit_path.read_text(encoding="utf-8"))
    gamma_subspace = json.loads(gamma_subspace_path.read_text(encoding="utf-8"))
    if (gamma_manifest.get("pressure_GPa") != 45.7
            or gamma_manifest.get("supercell_matrix") != np.eye(3, dtype=int).tolist()
            or len(gamma_manifest.get("cases", [])) != 32
            or gamma_audit.get("status")
            != "computed_GaN_endpoint_atomic_Gamma_1x1x1_requires_numerical_interpretation"
            or gamma_subspace.get("status")
            != "audited_GaN_endpoint_Gamma_atomic_subspaces_not_TS_modes"
            or gamma_subspace["source_sha256"].get("audit.json") != sha256(gamma_audit_path)):
        raise ValueError("endpoint Gamma provenance or audit status changed")

    hessian = {}
    source_paths = [gamma_manifest_path, gamma_audit_path, gamma_subspace_path]
    center_sha = None
    for label, distance in (("0p01", 0.01), ("0p02", 0.02)):
        root = joint_root / f"refined_hessian_{label}_v1"
        audit_root = joint_root / f"refined_hessian_{label}_audit_v1"
        manifest_path = root / "manifest.json"
        center_path = root / "center_POSCAR"
        audit_path = audit_root / "audit.json"
        npz_path = audit_root / "joint_hessian.npz"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if (manifest.get("phonon_supercell_matrix") != np.eye(3, dtype=int).tolist()
                or manifest.get("pressure_GPa") != 45.7
                or manifest.get("encut_eV") != 1000
                or manifest.get("step_A") != distance
                or audit.get("n_static_displacements") != 36
                or audit.get("step_A") != distance
                or not np.isfinite(audit.get("raw_hessian_reciprocity_relative_defect", np.nan))
                or audit["raw_hessian_reciprocity_relative_defect"] >= 0.01
                or not np.isfinite(audit.get("energy_gradient_max_abs_difference_eV_per_A", np.nan))
                or audit["source_sha256"].get("manifest") != sha256(manifest_path)
                or manifest["source_sha256"].get("center_POSCAR") != sha256(center_path)):
            raise ValueError(f"joint Hessian provenance changed at {label}")
        if center_sha is None:
            center_sha = sha256(center_path)
        elif sha256(center_path) != center_sha:
            raise ValueError("two joint Hessians use different center structures")
        with np.load(npz_path, allow_pickle=False) as data:
            hessian_matrix = np.asarray(data["hessian"], dtype=float)
            eigenvalues = np.asarray(data["eigenvalues"], dtype=float)
            eigenvectors = np.asarray(data["eigenvectors"], dtype=float)
        if (hessian_matrix.shape != (18, 18)
                or not np.allclose(hessian_matrix, hessian_matrix.T, atol=1e-8)
                or eigenvalues.shape != (15,) or eigenvectors.shape != (18, 15)
                or not np.allclose(eigenvalues, audit["translation_free_eigenvalues_eV_per_A2"], atol=1e-8)
                or not np.allclose(eigenvectors.T @ eigenvectors, np.eye(15), atol=1e-8)
                or not eigenvalues[0] < 0 < eigenvalues[1]):
            raise ValueError(f"joint index-one eigensystem changed at {label}")
        blocks = block_spectrum(*_translation_free_blocks(hessian_matrix, n_atoms=4))
        if (not np.allclose(blocks["lowest_joint_two_eV_per_A2"], eigenvalues[:2], atol=1e-7)
                or not blocks["lowest_strain_relaxed_atomic_schur_eV_per_A2"] < 0):
            raise ValueError(f"joint block decomposition differs from audited eigenvalues at {label}")
        hessian[label] = {
            "center": read(center_path, format="vasp"),
            "eigenvalue": float(eigenvalues[0]),
            "direction": eigenvectors[:, 0],
            "blocks": blocks,
            "numerical_integrity": {
                "raw_hessian_reciprocity_relative_defect": float(
                    audit["raw_hessian_reciprocity_relative_defect"]
                ),
                "energy_gradient_max_abs_difference_eV_per_A": float(
                    audit["energy_gradient_max_abs_difference_eV_per_A"]
                ),
                "center_translation_free_gradient_eV_per_A": float(
                    audit["center_gradient_translation_free_eV_per_A"]
                ),
            },
        }
        source_paths.extend((manifest_path, center_path, audit_path, npz_path))
    cross_step_overlap = float(abs(np.dot(
        hessian["0p01"]["direction"], hessian["0p02"]["direction"],
    )))
    if cross_step_overlap < 0.99:
        raise ValueError("joint unstable direction is not stable across displacement steps")

    phases = {}
    for phase in ("B4", "B1"):
        endpoint_path = gamma_root / f"{phase}_CONTCAR"
        endpoint = read(endpoint_path, format="vasp")
        if (endpoint.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]
                or gamma_manifest["source"][phase]["sha256"].get("CONTCAR")
                != sha256(endpoint_path)):
            raise ValueError(f"{phase} endpoint atom identity changed")
        selected_groups = [sorted(record["mode_indices"]) for record in
                           gamma_subspace["phases"][phase]["groups_ranked_by_path_amplitude"]]
        phase_results = {}
        for phonon_label in ("d0.01", "d0.02"):
            npz_path = gamma_root / f"{phase}_{phonon_label}.npz"
            projection_path = gamma_root / f"{phase}_{phonon_label}_path_projection.json"
            if (gamma_subspace["source_sha256"].get(npz_path.name) != sha256(npz_path)
                    or gamma_subspace["source_sha256"].get(projection_path.name)
                    != sha256(projection_path)):
                raise ValueError(f"{phase} Gamma force-constant provenance changed")
            projection = json.loads(projection_path.read_text(encoding="utf-8"))
            groups = [sorted(record["mode_indices"]) for record in
                      projection["degenerate_mode_subspaces"]
                      if min(record["mode_indices"]) >= 3]
            with np.load(npz_path, allow_pickle=False) as data:
                masses = np.asarray(data["masses_amu"], dtype=float)
                gamma_basis = np.asarray(data["eigenvectors_mass_weighted"], dtype=float)
                frequencies = np.asarray(data["frequencies_thz"], dtype=float)
            if (masses.shape != (4,) or frequencies.shape != (12,)
                    or projection["n_atoms"] != 4
                    or sorted(index for group in groups for index in group) != list(range(3, 12))):
                raise ValueError(f"{phase} Gamma eigensystem is incomplete")
            by_hessian_step = {}
            for hessian_label in ("0p01", "0p02"):
                item = hessian[hessian_label]
                rotation = _rotation_to_center(endpoint.cell.array, item["center"].cell.array)
                # Row Cartesian vectors: endpoint @ R = center; invert for projection.
                atomic_endpoint = item["direction"][:12].reshape(4, 3) @ rotation.T
                projected = project_atomic_direction(atomic_endpoint, masses, gamma_basis, groups)
                lookup = {tuple(group): fraction for group, fraction in zip(
                    groups, projected["group_fractions_of_mass_weighted_atomic_optical_direction"]
                )}
                projected.update({
                    "endpoint_to_center_rotation_frobenius_from_identity": float(
                        np.linalg.norm(rotation - np.eye(3))
                    ),
                    "path_selected_three_group_atomic_optical_fraction": float(sum(
                        lookup[tuple(group)] for group in selected_groups
                    )),
                    "optical_groups": [{
                        "mode_indices": group,
                        "mean_frequency_THz": float(np.mean(frequencies[group])),
                        "fraction": fraction,
                    } for group, fraction in zip(groups, projected[
                        "group_fractions_of_mass_weighted_atomic_optical_direction"
                    ])],
                })
                by_hessian_step[hessian_label] = projected
            phase_results[phonon_label] = by_hessian_step
            source_paths.extend((npz_path, projection_path))
        source_paths.append(endpoint_path)
        phase_results["path_selected_three_groups"] = selected_groups
        phases[phase] = phase_results

    return {
        "status": "audited_GaN_joint_unstable_direction_endpoint_Gamma_atomic_projection",
        "scope": "45.7 GPa GaN B4-to-B1, four-atom path cell; 1000 eV joint TS candidate and 600 eV endpoint Gamma bases",
        "joint_hessian_negative_eigenvalues_eV_per_A2": {
            label: item["eigenvalue"] for label, item in hessian.items()
        },
        "joint_hessian_atomic_strain_block_diagnostics": {
            label: item["blocks"] for label, item in hessian.items()
        },
        "joint_hessian_numerical_integrity": {
            label: item["numerical_integrity"] for label, item in hessian.items()
        },
        "joint_unstable_direction_cross_step_absolute_overlap": cross_step_overlap,
        "joint_unstable_direction_atomic_cartesian_squared_fraction": float(
            np.linalg.norm(hessian["0p01"]["direction"][:12]) ** 2
        ),
        "joint_unstable_direction_strain_coordinate_squared_fraction": float(
            np.linalg.norm(hessian["0p01"]["direction"][12:]) ** 2
        ),
        "phases": phases,
        "source_sha256": {(path.relative_to(gamma_root.parent) if path.is_relative_to(gamma_root)
                           else path.relative_to(joint_root.parent)).as_posix(): sha256(path)
                          for path in source_paths},
        "limitations": [
            "Only the atomic block of the joint unstable direction is projected; the strain block remains separate.",
            "Frozen atomic and frozen strain block positivity plus joint negative curvature identify coupling-driven instability within the declared local enthalpy coordinate chart; numerical eigenvalues depend on the cell-coordinate scale.",
            "The Hessians are symmetrized finite-difference estimates with nonzero raw reciprocity and energy-gradient residuals; two-step agreement supports a qualitative sign conclusion, not exact curvature values.",
            "Squared fractions describe mode direction, not energy, barrier contributions, or phonon softening at the TS.",
            "Endpoint Gamma phonons use the 600 eV production static contract; the joint Hessian is a separate 1000 eV diagnostic, so this is a geometric comparison only.",
            "Path-selected three groups were identified on the original 29-image chain and are not an independent path-training set.",
            "The four-atom cell Gamma spectrum does not address finite-q modes or polar non-analytic LO-TO corrections.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gamma-root", type=Path, required=True)
    parser.add_argument("--joint-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = analyze(args.gamma_root, args.joint_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "cross_step_overlap": result["joint_unstable_direction_cross_step_absolute_overlap"],
        "phases": {phase: result["phases"][phase]["d0.01"]["0p01"][
            "path_selected_three_group_atomic_optical_fraction"
        ] for phase in ("B4", "B1")},
    }))


if __name__ == "__main__":
    main()
