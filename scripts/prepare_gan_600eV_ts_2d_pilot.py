"""Prepare a guarded same-600-eV 2D joint-mode cut around the GaN candidate.

The first plane is deliberately local: eight 3x3-grid points plus four
half-step axial holdouts. All points share the original VASP electronic
contract. This is a frozen enthalpy cut, not a relaxed conditional PES.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sampling_coordinates() -> list[tuple[str, float, float]]:
    points = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if i == j == 0:
                continue
            name = f"grid_u{'m' if i < 0 else 'p' if i > 0 else 'z'}_v{'m' if j < 0 else 'p' if j > 0 else 'z'}"
            points.append((name, 0.02 * i, 0.0125 * j))
    points.extend((
        ("hold_um", -0.01, 0.0), ("hold_up", 0.01, 0.0),
        ("hold_vm", 0.0, -0.00625), ("hold_vp", 0.0, 0.00625),
    ))
    return points


def prepare(center_work: Path, hessian_root: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    center_outcar = center_work / "OUTCAR"
    hessian_audit_path = hessian_root / "audit.json"
    hessian_npz = hessian_root / "joint_hessian.npz"
    hessian_audit = json.loads(hessian_audit_path.read_text(encoding="utf-8"))
    if (hessian_audit.get("status") != "GaN_600eV_near_stationary_joint_hessian_one_step_not_TS_certificate"
            or hessian_audit.get("source_sha256", {}).get("center_OUTCAR") != sha256(center_outcar)
            or hessian_audit.get("n_static_displacements") != 36
            or hessian_audit.get("negative_count_below_minus_0p1_eV_per_A2") != 1
            or hessian_audit.get("center_gradient_translation_free_eV_per_A", 1) >= 0.01
            or any(sha256(center_work / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("same-600-eV center/Hessian or input contract is not audited")
    center = read(center_outcar)
    chart = JointCurvatureCoordinates(center, cell_scale_A=3.3982714330050063)
    with np.load(hessian_npz, allow_pickle=False) as archive:
        eigenvalues = np.asarray(archive["eigenvalues"], dtype=float)
        modes = np.asarray(archive["eigenvectors"], dtype=float)
    if (eigenvalues.shape != (15,) or modes.shape != (18, 15)
            or not np.allclose(eigenvalues, hessian_audit["eigenvalues_eV_per_A2"], atol=1e-8)
            or not np.allclose(modes.T @ modes, np.eye(15), atol=1e-8)
            or not eigenvalues[0] < 0 < eigenvalues[1]):
        raise ValueError("600-eV local mode basis differs from its audit")
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    write(output / "center_POSCAR", center, format="vasp", direct=True, sort=False)
    if not same_geometry(center, read(output / "center_POSCAR", format="vasp")):
        raise ValueError("2D center geometry roundtrip failed")
    records = []
    for name, q_u, q_v in sampling_coordinates():
        atoms = chart.displaced(q_u * modes[:, 0] + q_v * modes[:, 1])
        geometry = geometry_preflight(atoms)
        if geometry["empirical_near_symmetry_warning"]:
            raise ValueError(f"near-symmetry Bravais risk at {name}")
        case = output / "cases" / name
        case.mkdir()
        write(case / "POSCAR", atoms, format="vasp", direct=True, sort=False)
        if not same_geometry(atoms, read(case / "POSCAR", format="vasp")):
            raise ValueError(f"2D POSCAR roundtrip failed at {name}")
        for input_name in PRODUCTION_INPUT_SHA256:
            shutil.copy2(center_work / input_name, case / input_name)
        hashes = {input_name: sha256(case / input_name)
                  for input_name in ("POSCAR", *PRODUCTION_INPUT_SHA256)}
        if any(hashes[input_name] != digest
               for input_name, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError("2D VASP electronic setting changed while staging")
        (case / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
        )
        records.append({
            "name": name, "q_u_A": q_u, "q_v_A": q_v,
            **geometry, "input_sha256": hashes,
        })
    result = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_local_joint_mode_2D_frozen_pilot",
        "pressure_GPa": 45.7, "n_atoms": 4,
        "axes": {
            "u": "negative joint atom-strain eigenvector at the 600-eV candidate",
            "v": "lowest positive orthogonal joint atom-strain eigenvector at the same candidate",
            "eigenvalues_eV_per_A2": eigenvalues[:2].tolist(),
            "cell_scale_A": chart.cell_scale_A,
        },
        "grid_q_u_A": [-0.02, 0.0, 0.02],
        "grid_q_v_A": [-0.0125, 0.0, 0.0125],
        "center_reused_from_job": 27792422,
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "center_POSCAR_sha256": sha256(output / "center_POSCAR"),
        "source_sha256": {
            "center_OUTCAR": sha256(center_outcar),
            "hessian_audit": sha256(hessian_audit_path),
            "hessian_npz": sha256(hessian_npz),
            "preparer": sha256(Path(__file__)),
        },
        "cases": records,
        "limitations": [
            "Frozen local enthalpy cut, not a whole B4-to-B1 path surface.",
            "The 0.02-A Hessian has not yet been step-size validated.",
            "Symmetry-risk screening is empirical and cannot guarantee VASP parser success.",
        ],
    }
    (output / "manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("center-work", "hessian-root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.center_work, args.hessian_root, args.output)
    print(json.dumps({"status": result["status"], "n_cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
