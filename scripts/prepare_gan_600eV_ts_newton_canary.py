"""Stage one same-contract 600-eV GaN saddle-refinement static candidate.

The archived nonstationary 600-eV image-15 Hessian proposes a Newton step.
This writes geometry and provenance only; it never runs VASP or changes its
electronic settings. The candidate is not a transition-state certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read, write

from vcneb.joint_curvature import JointCurvatureCoordinates


PRODUCTION_INPUT_SHA256 = {
    "INCAR": "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772",
    "KPOINTS": "b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee",
    "POTCAR": "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def basal_projection_A(cell: np.ndarray) -> float:
    """Maximum projection of the third cell vector onto the first two."""

    matrix = np.asarray(cell, dtype=float)
    if matrix.shape != (3, 3) or not np.isfinite(matrix).all():
        raise ValueError("invalid cell")
    return float(max(abs(np.dot(matrix[2], matrix[i])) / np.linalg.norm(matrix[i])
                     for i in (0, 1)))


def prepare(source: Path, trajectory_path: Path, manifest_path: Path,
            audit_path: Path, hessian_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    summary_path = source / "source_summary.json"
    center_poscar_path = manifest_path.parent / "step_0p02_v1" / "center_POSCAR"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    params = summary.get("calculator_parameters", {})
    if (params.get("encut") != 600.0 or params.get("kpts") != [8, 8, 6]
            or params.get("isym") != -1 or params.get("symprec") != 1e-4
            or params.get("ediff") != 1e-7 or params.get("setups") != {"Ga": "_d", "N": ""}
            or manifest.get("pressure_GPa") != 45.7
            or manifest.get("source_sha256", {}).get("trajectory") != sha256(trajectory_path)
            or manifest.get("center_POSCAR_sha256") != sha256(center_poscar_path)
            or audit.get("source_sha256", {}).get("manifest") != sha256(manifest_path)
            or audit.get("status") != "GaN_image15_local_joint_curvature_at_nonstationary_NEB_candidate"
            or len(audit.get("cases", [])) != 36
            or any({key: case.get("input_sha256", {}).get(key) for key in PRODUCTION_INPUT_SHA256}
                   != PRODUCTION_INPUT_SHA256 for case in audit["cases"])):
        raise ValueError("600-eV GaN source or production input contract changed")
    center = read(trajectory_path, index=15)
    if center.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]:
        raise ValueError("source atom identity changed")
    archived_center = read(center_poscar_path, format="vasp")
    delta_scaled = center.get_scaled_positions(wrap=False) - archived_center.get_scaled_positions(wrap=False)
    delta_scaled -= np.rint(delta_scaled)
    if (not np.allclose(center.cell.array, archived_center.cell.array, atol=2e-5, rtol=0)
            or not np.allclose(delta_scaled, 0, atol=2e-5, rtol=0)):
        raise ValueError("compact final-chain image 15 differs from the Hessian center")
    chart = JointCurvatureCoordinates(center, float(manifest["cell_scale_A"]))
    with np.load(hessian_path, allow_pickle=False) as data:
        eigenvalues = np.asarray(data["eigenvalues"], dtype=float)
        eigenvectors = np.asarray(data["eigenvectors"], dtype=float)
        gradient = np.asarray(data["center_gradient"], dtype=float)
    if (eigenvalues.shape != (15,) or eigenvectors.shape != (18, 15)
            or gradient.shape != (18,) or not np.isfinite(gradient).all()
            or not np.allclose(eigenvalues, audit["translation_free_eigenvalues_eV_per_A2"],
                               atol=1e-8, rtol=0)
            or not np.allclose(eigenvectors.T @ eigenvectors, np.eye(15), atol=1e-8)
            or not eigenvalues[0] < 0 < eigenvalues[1]
            or not np.isclose(np.linalg.norm(chart.translation_free_basis().T @ gradient),
                              audit["center_gradient_translation_free_eV_per_A"], atol=1e-8)):
        raise ValueError("archived 600-eV Hessian/gradient does not match its audit")
    correction = -eigenvectors @ ((eigenvectors.T @ gradient) / eigenvalues)
    atom_max = float(np.max(np.linalg.norm(correction[:12].reshape(4, 3), axis=1)))
    cell_norm = float(np.linalg.norm(correction[12:]))
    if (not np.isfinite(correction).all() or np.linalg.norm(correction) > 0.025
            or atom_max > 0.01 or cell_norm > 0.02):
        raise ValueError("Newton proposal exceeds the prespecified local trust radius")
    candidate = chart.displaced(correction)
    distances = candidate.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    basal = basal_projection_A(candidate.cell.array)
    if (candidate.get_volume() <= 0 or float(np.min(distances)) < 1.4
            or basal > 1e-4):
        raise ValueError("Newton proposal failed the local geometry-risk gate")
    output.mkdir(parents=True)
    write(output / "POSCAR", candidate, format="vasp", direct=True, sort=False)
    reread = read(output / "POSCAR", format="vasp")
    if (not np.allclose(reread.cell.array, candidate.cell.array, atol=1e-10, rtol=0)
            or not np.allclose(reread.get_scaled_positions(wrap=False)
                               - candidate.get_scaled_positions(wrap=False), 0,
                               atol=1e-10, rtol=0)):
        raise ValueError("POSCAR roundtrip changed the candidate")
    result = {
        "status": "geometry_preflight_passed_no_DFT",
        "purpose": "GaN_45p7_600eV_image15_one_step_joint_Newton_static_canary",
        "claim_limit": "Nonstationary image-15 Hessian proposal only; no TS certification",
        "electronic_contract": {
            "ENCUT_eV": 600, "PBE_PAW": "Ga_d+N", "kpoints": "Gamma 8x8x6",
            "ISYM": -1, "SYMPREC": 1e-4, "EDIFF": 1e-7,
            "original_input_sha256": PRODUCTION_INPUT_SHA256,
        },
        "pressure_GPa": 45.7,
        "center_gradient_translation_free_eV_per_A": audit["center_gradient_translation_free_eV_per_A"],
        "newton_joint_norm_A": float(np.linalg.norm(correction)),
        "newton_max_atom_displacement_A": atom_max,
        "newton_cell_scaled_norm_A": cell_norm,
        "candidate_volume_A3": float(candidate.get_volume()),
        "candidate_minimum_distance_A": float(np.min(distances)),
        "candidate_basal_projection_A": basal,
        "geometry_risk_gate_basal_projection_A": 1e-4,
        "POSCAR_sha256": sha256(output / "POSCAR"),
        "source_sha256": {
            "summary": sha256(summary_path), "trajectory": sha256(trajectory_path),
            "center_POSCAR": sha256(center_poscar_path),
            "manifest": sha256(manifest_path), "audit": sha256(audit_path),
            "hessian": sha256(hessian_path), "preparer": sha256(Path(__file__)),
        },
        "limitations": [
            "Geometry-risk gate is empirical for this local cell family, not a VASP Bravais proof.",
            "No electronic input is staged or VASP job submitted by this script.",
        ],
    }
    (output / "manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "trajectory", "manifest", "audit", "hessian", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source, args.trajectory, args.manifest, args.audit,
                             args.hessian, args.output)))


if __name__ == "__main__":
    main()
