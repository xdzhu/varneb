"""Audit a single GaN 600-eV Newton static against its original VCNEB image."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from vcneb.joint_curvature import JointCurvatureCoordinates


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(work: Path, trajectory: Path, source: Path,
          manifest_path: Path, audit_path: Path, hessian_path: Path,
          slurm_job: int) -> dict:
    staged = json.loads((work / "manifest.json").read_text(encoding="utf-8"))
    source_summary = source / "source_summary.json"
    center_poscar = manifest_path.parent / "step_0p02_v1" / "center_POSCAR"
    required_source = {
        "summary": sha256(source_summary), "trajectory": sha256(trajectory),
        "center_POSCAR": sha256(center_poscar), "manifest": sha256(manifest_path),
        "audit": sha256(audit_path), "hessian": sha256(hessian_path),
        "preparer": sha256(Path(__file__).with_name("prepare_gan_600eV_ts_newton_canary.py")),
    }
    if (staged.get("status") != "geometry_preflight_passed_no_DFT"
            or staged.get("purpose") != "GaN_45p7_600eV_image15_one_step_joint_Newton_static_canary"
            or staged.get("source_sha256") != required_source
            or staged.get("electronic_contract", {}).get("ENCUT_eV") != 600):
        raise ValueError("staged 600-eV source/contract changed")
    expected_inputs = dict(staged["electronic_contract"]["original_input_sha256"])
    expected_inputs["POSCAR"] = staged["POSCAR_sha256"]
    if any(sha256(work / name) != digest for name, digest in expected_inputs.items()):
        raise ValueError("VASP static input differs from the original 600-eV contract")
    outcar = work / "OUTCAR"
    log = outcar.read_text(encoding="utf-8", errors="replace")
    if ("General timing and accounting informations" not in log
            or "aborting loop because EDIFF is reached" not in log
            or "Inconsistent Bravais lattice types" in log):
        raise ValueError("VASP did not finish a converged electronic static")
    atoms = read(outcar)
    proposed = read(work / "POSCAR", format="vasp")
    source_center = read(trajectory, index=15)
    if (atoms.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]
            or proposed.get_chemical_symbols() != atoms.get_chemical_symbols()
            or not np.allclose(atoms.cell.array, proposed.cell.array, atol=2e-5, rtol=0)):
        raise ValueError("VASP output geometry or atom identity changed")
    dq = atoms.get_scaled_positions(wrap=False) - proposed.get_scaled_positions(wrap=False)
    dq -= np.rint(dq)
    if not np.allclose(dq, 0, atol=2e-5, rtol=0):
        raise ValueError("VASP output coordinates differ from the proposed static input")
    pressure = 45.7 * GPa
    chart = JointCurvatureCoordinates(source_center, float(
        json.loads(manifest_path.read_text(encoding="utf-8"))["cell_scale_A"]
    ))
    basis = chart.translation_free_basis()
    with np.load(hessian_path, allow_pickle=False) as archive:
        center_gradient = np.asarray(archive["center_gradient"], dtype=float)
    original_audit = json.loads(audit_path.read_text(encoding="utf-8"))
    center_norm = float(np.linalg.norm(basis.T @ center_gradient))
    if not np.isclose(center_norm,
                      original_audit["center_gradient_translation_free_eV_per_A"],
                      rtol=0, atol=1e-8):
        raise ValueError("original 600-eV center gradient no longer matches its audit")
    forces = np.asarray(atoms.get_forces(), dtype=float)
    stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    if (forces.shape != (4, 3) or stress.shape != (3, 3)
            or not np.isfinite(forces).all() or not np.isfinite(stress).all()):
        raise ValueError("static forces or stress incomplete")
    gradient = chart.enthalpy_gradient(atoms, forces, stress, pressure)
    current_norm = float(np.linalg.norm(basis.T @ gradient))
    original_h = float(source_center.get_potential_energy() + pressure * source_center.get_volume())
    current_h = float(atoms.get_potential_energy() + pressure * atoms.get_volume())
    return {
        "status": "GaN_600eV_one_step_static_audited_not_TS_certificate",
        "slurm_job": slurm_job,
        "input_contract": staged["electronic_contract"],
        "center_gradient_translation_free_eV_per_A": center_norm,
        "candidate_gradient_translation_free_eV_per_A": current_norm,
        "gradient_ratio": current_norm / center_norm,
        "center_enthalpy_eV_per_cell": original_h,
        "candidate_enthalpy_eV_per_cell": current_h,
        "enthalpy_change_eV_per_cell": current_h - original_h,
        "candidate_max_atomic_force_eV_per_A": float(np.max(np.linalg.norm(forces, axis=1))),
        "candidate_gradient_eV_per_A": gradient.tolist(),
        "input_sha256": expected_inputs,
        "source_sha256": required_source,
        "OUTCAR_sha256": sha256(outcar),
        "auditor_sha256": sha256(Path(__file__)),
        "interpretation": (
            "A single static gradient change is not saddle convergence. Recompute a "
            "same-600-eV Hessian at a stationary candidate and verify both basin links."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "trajectory", "source", "manifest", "audit", "hessian", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--slurm-job", type=int, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work, args.trajectory, args.source, args.manifest,
                   args.audit, args.hessian, args.slurm_job)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "gradient_ratio": result["gradient_ratio"]}))


if __name__ == "__main__":
    main()
