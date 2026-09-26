"""Audit the isolated GaN follow-up static; never call it a certified TS."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.audit_gan_ts_newton_probe import failed_bravais_case
from scripts.prepare_gan_ts_bracketed_trial import cross_axis_skew_angstrom
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256
from vcneb.joint_curvature import JointCurvatureCoordinates


def audit(work_root: Path, center_poscar: Path, newton_manifest: Path,
          bracketed_manifest: Path, bracketed_audit: Path,
          completed_source_outcar: Path, hessian_archive: Path) -> dict:
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    first = json.loads(newton_manifest.read_text(encoding="utf-8"))
    prior = json.loads(bracketed_audit.read_text(encoding="utf-8"))
    sources = {
        "center_POSCAR": center_poscar, "newton_manifest": newton_manifest,
        "bracketed_manifest": bracketed_manifest, "bracketed_audit": bracketed_audit,
        "completed_OUTCAR": completed_source_outcar, "hessian_archive": hessian_archive,
        "preparer": Path(__file__).with_name("prepare_gan_ts_followup_trial.py"),
    }
    if (manifest["purpose"] != "GaN_image15_followup_1000eV_local_stationarity_probe_not_TS_certificate"
            or manifest["status"] != "inputs_finalized_no_DFT"
            or manifest["evaluation_encut_eV"] != 1000
            or len(manifest["cases"]) != 1
            or any(sha256(path) != manifest["source_sha256"][key] for key, path in sources.items())):
        raise ValueError("follow-up input or source provenance changed")
    reference = read(center_poscar, format="vasp")
    coordinates = JointCurvatureCoordinates(reference, float(first["cell_scale_A"]))
    step = 0.625 * np.asarray(first["correction_A"], dtype=float)
    step += np.asarray(manifest["correction_A"], dtype=float)
    expected = coordinates.displaced(step)
    case = manifest["cases"][0]
    directory = work_root / "cases" / case["name"]
    supplied = read(directory / "POSCAR", format="vasp")
    if (case["name"] != "newton_followup_full" or not same_geometry(expected, supplied)
            or not np.isclose(cross_axis_skew_angstrom(supplied.cell.array),
                              manifest["empirical_bravais_risk"]["candidate_skew_A"],
                              rtol=0, atol=1e-10)):
        raise ValueError("follow-up candidate geometry changed")
    outcar_text = (directory / "OUTCAR").read_text(encoding="utf-8", errors="replace")
    if "Inconsistent Bravais lattice types found for crystalline and" in outcar_text:
        result = failed_bravais_case(directory, case)
    else:
        atoms, result = completed_case(
            directory, {**case, "POSCAR_sha256": case["input_sha256"]["POSCAR"]},
        )
        gradient = coordinates.enthalpy_gradient(
            atoms, atoms.get_forces(), atoms.get_stress(voigt=False), 45.7 * GPa,
        )
        norm = float(np.linalg.norm(coordinates.translation_free_basis().T @ gradient))
        result.update({
            "status": "complete_converged_static_not_TS_certificate",
            "enthalpy_eV_per_cell": float(atoms.get_potential_energy() + 45.7 * GPa * atoms.get_volume()),
            "translation_free_gradient_eV_per_A": norm,
            "gradient_ratio_to_prior": norm / prior["result"]["translation_free_gradient_eV_per_A"],
            "gradient_components_eV_per_A": gradient.tolist(),
        })
    return {
        "status": "GaN_followup_local_trial_audited_not_TS_certificate",
        "result": result,
        "source_sha256": {"manifest": sha256(manifest_path),
                          "completed_source_OUTCAR": sha256(completed_source_outcar),
                          "auditor": sha256(Path(__file__))},
        "interpretation_boundary": "A lower gradient supports local refinement, not a certified variable-cell TS.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work-root", "center-poscar", "newton-manifest", "bracketed-manifest",
                 "bracketed-audit", "completed-source-outcar", "hessian-archive", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work_root, args.center_poscar, args.newton_manifest,
                   args.bracketed_manifest, args.bracketed_audit,
                   args.completed_source_outcar, args.hessian_archive)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "result": result["result"]["status"],
                      "gradient": result["result"].get("translation_free_gradient_eV_per_A")}))


if __name__ == "__main__":
    main()
