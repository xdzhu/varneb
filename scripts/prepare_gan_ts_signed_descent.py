"""Stage two same-contract VASP statics along GaN's audited unstable mode.

This checks local downhill behavior before any basin-following optimization.
It never treats a signed static probe as a B4/B1 basin connection.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256
from vcneb.joint_curvature import JointCurvatureCoordinates


def basal_projection_A(cell: np.ndarray) -> float:
    """This GaN input's near-orthogonal third-vector Bravais-risk diagnostic."""

    return float(max(abs(np.dot(cell[2], cell[i])) / np.linalg.norm(cell[i])
                     for i in (0, 1)))


def prepare(center_dir: Path, hessian_manifest: Path, hessian_audit: Path,
            hessian_npz: Path, step_comparison: Path, work_root: Path) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    manifest = json.loads(hessian_manifest.read_text(encoding="utf-8"))
    audit = json.loads(hessian_audit.read_text(encoding="utf-8"))
    comparison = json.loads(step_comparison.read_text(encoding="utf-8"))
    source = manifest.get("source_sha256", {})
    source_comparison = comparison.get("source_sha256", {})
    archive_comparison = comparison.get("archived_mode_comparison", {})
    if (manifest.get("purpose") != "GaN_near_stationary_candidate_1000eV_joint_atomic_strain_Hessian_not_TS_certificate"
            or manifest.get("step_A") != 0.01
            or manifest.get("phonon_supercell_matrix") != np.eye(3, dtype=int).tolist()
            or manifest.get("electronic_kpoints") != [8, 8, 6]
            or audit.get("status") != "GaN_near_stationary_1000eV_joint_curvature_one_step_audited_not_TS_certificate"
            or audit.get("source_sha256", {}).get("manifest") != sha256(hessian_manifest)
            or source.get("center_POSCAR") != sha256(center_dir / "POSCAR")
            or source.get("center_OUTCAR") != sha256(center_dir / "OUTCAR")
            or source_comparison.get("second_audit") != sha256(hessian_audit)
            or archive_comparison.get("second_hessian_npz_sha256") != sha256(hessian_npz)
            or archive_comparison.get("lowest_mode_absolute_overlap", 0) < 0.99
            or not all(comparison.get("negative_index_stable_by_cutoff", {}).values())
            or comparison.get("kind") != "gan"):
        raise ValueError("same-center, two-step audited GaN unstable mode is required")
    with np.load(hessian_npz) as archived:
        eigenvalues = np.asarray(archived["eigenvalues"], dtype=float)
        eigenvectors = np.asarray(archived["eigenvectors"], dtype=float)
        center_gradient = np.asarray(archived["center_gradient"], dtype=float)
    if (eigenvalues.shape != (15,) or eigenvectors.shape != (18, 15)
            or center_gradient.shape != (18,)
            or not np.allclose(eigenvalues, audit["translation_free_eigenvalues_eV_per_A2"], atol=1e-8, rtol=0)
            or not eigenvalues[0] < -3.5 or not eigenvalues[1] > 2.0):
        raise ValueError("near-TS archived eigensystem does not match raw audit")
    center = read(center_dir / "POSCAR", format="vasp")
    scale = float(manifest["cell_scale_A"])
    chart = JointCurvatureCoordinates(center, scale)
    amplitude = 0.10
    unstable = eigenvectors[:, 0]
    work_root.mkdir(parents=True)
    shutil.copy2(center_dir / "POSCAR", work_root / "center_POSCAR")
    records = []
    for sign, label in ((1, "plus"), (-1, "minus")):
        name = f"negative_mode_{label}"
        atoms = chart.displaced(sign * amplitude * unstable)
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        minimum = float(np.min(distances))
        risk = basal_projection_A(atoms.cell.array)
        if (len(atoms) != 4 or minimum < 1.4 or np.linalg.cond(atoms.cell.array) > 50
                or risk > 0.00018081):
            raise ValueError(f"unsafe or Bravais-risk signed trial: {name}")
        directory = work_root / "cases" / name
        directory.mkdir(parents=True)
        write(directory / "POSCAR", atoms, format="vasp", direct=True, sort=False)
        if not same_geometry(atoms, read(directory / "POSCAR", format="vasp")):
            raise ValueError(f"signed trial changed during POSCAR roundtrip: {name}")
        for filename in ("INCAR", "KPOINTS", "POTCAR"):
            shutil.copy2(center_dir / filename, directory / filename)
        hashes = {filename: sha256(directory / filename)
                  for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
        (directory / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8",
        )
        records.append({
            "name": name, "sign": sign, "input_sha256": hashes,
            "POSCAR_sha256": hashes["POSCAR"],
            "volume_A3": float(atoms.get_volume()),
            "minimum_distance_A": minimum,
            "basal_projection_A": risk,
            "quadratic_predicted_delta_enthalpy_eV": float(
                sign * amplitude * (unstable @ center_gradient)
                + 0.5 * amplitude * amplitude * eigenvalues[0]
            ),
        })
    result = {
        "purpose": "GaN_near_TS_1000eV_signed_unstable_mode_static_downhill_probe_not_basin_link",
        "status": "inputs_finalized_no_DFT", "n_atoms": 4,
        "pressure_GPa": 45.7, "encut_eV": 1000,
        "phonon_supercell_matrix": np.eye(3, dtype=int).tolist(),
        "electronic_kpoints": [8, 8, 6],
        "cell_scale_A": scale, "amplitude_A": amplitude,
        "lowest_eigenvalue_eV_per_A2": float(eigenvalues[0]),
        "center_mode_gradient_eV_per_A": float(unstable @ center_gradient),
        "center_enthalpy_eV_per_cell": float(audit["center_enthalpy_eV_per_cell"]),
        "bravais_risk_empirical_limit_A": 0.00018081,
        "mode_sign_is_arbitrary": True,
        "cases": records,
        "source_sha256": {
            "center_POSCAR": sha256(center_dir / "POSCAR"),
            "center_OUTCAR": sha256(center_dir / "OUTCAR"),
            "INCAR": sha256(center_dir / "INCAR"),
            "KPOINTS": sha256(center_dir / "KPOINTS"),
            "POTCAR": sha256(center_dir / "POTCAR"),
            "hessian_manifest": sha256(hessian_manifest),
            "hessian_audit": sha256(hessian_audit),
            "hessian_npz": sha256(hessian_npz),
            "step_comparison": sha256(step_comparison),
            "preparer": sha256(Path(__file__)),
        },
        "limitations": [
            "The two signed structures are statics, not relaxed basin endpoints.",
            "The Bravais-risk limit is empirical for this GaN input only.",
            "This does not alter the original 600 eV VCNEB path.",
        ],
    }
    (work_root / "manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("center-dir", "hessian-manifest", "hessian-audit", "hessian-npz",
                 "step-comparison", "work-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.center_dir, args.hessian_manifest, args.hessian_audit,
                     args.hessian_npz, args.step_comparison, args.work_root)
    print(json.dumps({"status": result["status"], "n_cases": len(result["cases"]),
                      "predicted_deltas_eV": [c["quadratic_predicted_delta_enthalpy_eV"]
                                              for c in result["cases"]]}))


if __name__ == "__main__":
    main()
