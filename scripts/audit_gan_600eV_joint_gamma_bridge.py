"""Audit a same-600-eV GaN joint-Hessian / endpoint-Gamma comparison.

Only the atomic component of the joint direction is projected onto the
endpoint optical subspaces. Fractions are geometric, not barrier shares.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.audit_gan_ts_endpoint_gamma_bridge import (
    _rotation_to_center,
    _translation_free_blocks,
    block_spectrum,
    project_atomic_direction,
)
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(gamma_root: Path, hessian_root: Path, center_outcar: Path) -> dict:
    manifest_path = gamma_root / "manifest.json"
    gamma_audit_path = gamma_root / "audit.json"
    subspace_path = gamma_root / "subspace_audit_v2.json"
    hessian_audit_path = hessian_root / "audit.json"
    hessian_npz_path = hessian_root / "joint_hessian.npz"
    center_audit_path = ROOT / "benchmarks/numerical_integrity/gan_600eV_ts_newton_canary_20260928.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    gamma_audit = json.loads(gamma_audit_path.read_text(encoding="utf-8"))
    subspace = json.loads(subspace_path.read_text(encoding="utf-8"))
    audit = json.loads(hessian_audit_path.read_text(encoding="utf-8"))
    center_audit = json.loads(center_audit_path.read_text(encoding="utf-8"))
    if (manifest.get("pressure_GPa") != 45.7
            or manifest.get("supercell_matrix") != np.eye(3, dtype=int).tolist()
            or gamma_audit.get("status")
            != "computed_GaN_endpoint_atomic_Gamma_1x1x1_requires_numerical_interpretation"
            or subspace.get("status")
            != "audited_GaN_endpoint_Gamma_atomic_subspaces_not_TS_modes"
            or subspace["source_sha256"].get("audit.json") != sha256(gamma_audit_path)
            or audit.get("status")
            != "GaN_600eV_near_stationary_joint_hessian_one_step_not_TS_certificate"
            or audit.get("n_static_displacements") != 36
            or audit.get("step_A") != 0.02
            or audit.get("negative_count_below_minus_0p1_eV_per_A2") != 1
            or audit["source_sha256"].get("center_OUTCAR") != sha256(center_outcar)
            or audit["source_sha256"].get("center_audit") != sha256(center_audit_path)
            or center_audit.get("input_contract", {}).get("ENCUT_eV") != 600
            or center_audit.get("OUTCAR_sha256") != sha256(center_outcar)
            or any(center_audit.get("input_sha256", {}).get(name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("Gamma or joint-Hessian provenance changed")
    for phase in ("B4", "B1"):
        if any(manifest["source"][phase]["sha256"].get(name) != digest
               for name, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"{phase} endpoint is not the 600-eV contract")
    if any(case["input_sha256"].get(name) != digest
           for case in audit["cases"] for name, digest in PRODUCTION_INPUT_SHA256.items()):
        raise ValueError("joint-Hessian displacement has a different electronic contract")

    center = read(center_outcar)
    with np.load(hessian_npz_path, allow_pickle=False) as source:
        hessian = np.asarray(source["hessian"], dtype=float)
        eigenvalues = np.asarray(source["eigenvalues"], dtype=float)
        eigenvectors = np.asarray(source["eigenvectors"], dtype=float)
    if (hessian.shape != (18, 18) or eigenvectors.shape != (18, 15)
            or not np.allclose(eigenvectors.T @ eigenvectors, np.eye(15), atol=1e-8)
            or not np.allclose(eigenvalues, audit["eigenvalues_eV_per_A2"], atol=1e-8)
            or not np.isclose(audit["path_tangent_squared_overlap_with_modes"][0],
                              0.9957555646436711, atol=1e-7)):
        raise ValueError("600-eV joint eigensystem changed")
    blocks = block_spectrum(*_translation_free_blocks(hessian, n_atoms=4))
    direction = eigenvectors[:, 0]
    phases: dict[str, dict] = {}
    sources = [manifest_path, gamma_audit_path, subspace_path,
               hessian_audit_path, hessian_npz_path, center_audit_path, center_outcar]
    for phase in ("B4", "B1"):
        endpoint_path = gamma_root / f"{phase}_CONTCAR"
        npz_path = gamma_root / f"{phase}_d0.01.npz"
        projection_path = gamma_root / f"{phase}_d0.01_path_projection.json"
        if (sha256(endpoint_path) != manifest["source"][phase]["sha256"]["CONTCAR"]
                or sha256(npz_path) != subspace["source_sha256"][npz_path.name]
                or sha256(projection_path) != subspace["source_sha256"][projection_path.name]):
            raise ValueError(f"{phase} endpoint-Gamma source changed")
        endpoint = read(endpoint_path, format="vasp")
        if endpoint.get_chemical_symbols() != center.get_chemical_symbols():
            raise ValueError("endpoint and candidate atom identity differ")
        projection = json.loads(projection_path.read_text(encoding="utf-8"))
        groups = [sorted(item["mode_indices"])
                  for item in projection["degenerate_mode_subspaces"]
                  if min(item["mode_indices"]) >= 3]
        with np.load(npz_path, allow_pickle=False) as source:
            masses = np.asarray(source["masses_amu"], dtype=float)
            gamma_basis = np.asarray(source["eigenvectors_mass_weighted"], dtype=float)
            frequencies = np.asarray(source["frequencies_thz"], dtype=float)
        rotation = _rotation_to_center(endpoint.cell.array, center.cell.array)
        atomic_endpoint = direction[:12].reshape(4, 3) @ rotation.T
        projected = project_atomic_direction(atomic_endpoint, masses, gamma_basis, groups)
        fractions = projected["group_fractions_of_mass_weighted_atomic_optical_direction"]
        phases[phase] = {
            "optical_groups": [{
                "mode_indices": group,
                "mean_frequency_THz": float(np.mean(frequencies[group])),
                "fraction": fraction,
            } for group, fraction in zip(groups, fractions)],
            "acoustic_leakage_after_center_of_mass_removal":
                projected["acoustic_squared_fraction_after_com_removal"],
            "endpoint_to_candidate_rotation_frobenius_from_identity":
                float(np.linalg.norm(rotation - np.eye(3))),
        }
        sources.extend((endpoint_path, npz_path, projection_path))
    return {
        "status": "GaN_same_600eV_joint_endpoint_Gamma_direction_audited_not_TS_certificate",
        "scope": "GaN B4-to-B1 at 45.7 GPa, VASP 6.3.2, all DFT at production ENCUT=600 eV",
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "joint_hessian_step_A": 0.02,
        "joint_block_lowest_curvatures_eV_per_A2": blocks,
        "joint_eigenvalues_first_two_eV_per_A2": eigenvalues[:2].tolist(),
        "path_tangent_squared_overlap_with_negative_mode":
            audit["path_tangent_squared_overlap_with_modes"][0],
        "joint_direction_atomic_cartesian_squared_fraction": float(np.sum(direction[:12] ** 2)),
        "joint_direction_strain_coordinate_squared_fraction": float(np.sum(direction[12:] ** 2)),
        "numerical_integrity": {
            "center_gradient_eV_per_A": audit["center_gradient_translation_free_eV_per_A"],
            "energy_gradient_max_abs_difference_eV_per_A":
                audit["energy_gradient_max_abs_difference_eV_per_A"],
            "energy_hessian_diagonal_max_abs_difference_eV_per_A2":
                audit["energy_hessian_diagonal_max_abs_difference_eV_per_A2"],
        },
        "phases": phases,
        "source_sha256": {str(path): sha256(path) for path in sources},
        "limitations": [
            "Single 0.02 A full joint-Hessian step; directional half-step energy tests do not replace a full second Hessian.",
            "Energy-gradient/stress discrepancy prevents strict transition-state certification.",
            "Only the atomic part of the joint direction is projected; fractions are not energy shares.",
            "Endpoint Gamma modes are not transition-state phonons or finite-q modes.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gamma-root", type=Path, required=True)
    parser.add_argument("--hessian-root", type=Path, required=True)
    parser.add_argument("--center-outcar", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true",
                        help="replace only the explicitly selected generated audit")
    args = parser.parse_args()
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(args.output)
    result = analyze(args.gamma_root, args.hessian_root, args.center_outcar)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"],
                      "lowest_joint": result["joint_eigenvalues_first_two_eV_per_A2"][0]}))


if __name__ == "__main__":
    main()
