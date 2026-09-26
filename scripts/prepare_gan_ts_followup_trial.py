"""Stage one audited GaN 1000 eV local saddle trial after the 0.625 step.

The 600 eV Hessian is only an approximate Jacobian. This script never changes
the VASP physics contract, source trajectory, or a completed calculation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_ts_bracketed_trial import cross_axis_skew_angstrom
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256
from vcneb.joint_curvature import JointCurvatureCoordinates


def prepare(
    *, center_poscar: Path, newton_manifest: Path, bracketed_manifest: Path,
    bracketed_audit: Path, completed_case_dir: Path, hessian_archive: Path,
    work_root: Path,
) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    original = json.loads(newton_manifest.read_text(encoding="utf-8"))
    bracket = json.loads(bracketed_manifest.read_text(encoding="utf-8"))
    audit = json.loads(bracketed_audit.read_text(encoding="utf-8"))
    if (original["status"] != "inputs_finalized_no_DFT"
            or bracket["status"] != "inputs_finalized_no_DFT"
            or bracket["evaluation_encut_eV"] != 1000
            or bracket["source_hessian_encut_eV"] != 600
            or bracket["source_sha256"]["previous_manifest"] != sha256(newton_manifest)
            or audit["source_sha256"]["manifest"] != sha256(bracketed_manifest)
            or original["source_sha256"]["hessian_archive"] != sha256(hessian_archive)
            or original["center_POSCAR_sha256"] != sha256(center_poscar)
            or len(bracket["cases"]) != 1
            or bracket["cases"][0]["name"] != "newton_0p625"
            or audit["result"]["status"] != "complete_converged_static_not_TS_certificate"):
        raise ValueError("source trajectory, Hessian, or completed trial provenance failed")
    case = bracket["cases"][0]
    source, parsed = completed_case(completed_case_dir, case)
    if (sha256(completed_case_dir / "OUTCAR") != audit["result"]["outcar_sha256"]
            or parsed["input_sha256"] != audit["result"]["input_sha256"]):
        raise ValueError("completed source is not the audited 0.625 calculation")

    reference = read(center_poscar, format="vasp")
    coordinates = JointCurvatureCoordinates(reference, float(original["cell_scale_A"]))
    first = np.asarray(original["correction_A"], dtype=float)
    expected = coordinates.displaced(0.625 * first)
    if not same_geometry(expected, source, tolerance=2e-5):
        raise ValueError("completed trial is not the declared 0.625 geometry")
    with np.load(hessian_archive, allow_pickle=False) as archived:
        eigenvalues = archived["eigenvalues"].copy()
        eigenvectors = archived["eigenvectors"].copy()
    if (eigenvalues.shape != (15,) or eigenvectors.shape != (18, 15)
            or np.count_nonzero(eigenvalues < -0.1) != 1
            or eigenvalues[0] > -3.5 or eigenvalues[1] < 2.0
            or not np.allclose(eigenvectors.T @ eigenvectors, np.eye(15), atol=1e-8)):
        raise ValueError("the archived single-negative Hessian is not trusted")
    gradient = coordinates.enthalpy_gradient(
        source, source.get_forces(), source.get_stress(voigt=False), 45.7 * GPa,
    )
    physical_gradient = float(np.linalg.norm(coordinates.translation_free_basis().T @ gradient))
    if (not np.allclose(gradient, audit["result"]["gradient_components_eV_per_A"], atol=1e-7, rtol=0)
            or not np.isclose(physical_gradient, audit["result"]["translation_free_gradient_eV_per_A"], atol=1e-7, rtol=0)):
        raise ValueError("raw VASP gradient does not reproduce the prior audit")
    correction = -eigenvectors @ ((eigenvectors.T @ gradient) / eigenvalues)
    if (not np.isfinite(correction).all() or np.linalg.norm(correction) > 0.012
            or np.max(np.linalg.norm(correction[:12].reshape(4, 3), axis=1)) > 0.005
            or np.linalg.norm(correction[12:]) > 0.008):
        raise ValueError("second correction exceeds the prespecified local trust radius")
    candidate = coordinates.displaced(0.625 * first + correction)
    distances = candidate.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    skew = cross_axis_skew_angstrom(candidate.cell.array)
    risk_limit = float(bracket["empirical_bravais_risk"]["allowed_skew_A"])
    if (np.min(distances) < 1.4 or not skew < risk_limit
            or len(candidate) != 4 or candidate.get_chemical_symbols() != reference.get_chemical_symbols()):
        raise ValueError("unsafe geometry or empirical VASP Bravais-risk limit exceeded")
    source_incar = (completed_case_dir / "INCAR").read_text(encoding="utf-8")
    source_kpoints = (completed_case_dir / "KPOINTS").read_text(encoding="utf-8")
    if ("ENCUT = 1000.000000" not in source_incar or "ISYM = -1" not in source_incar
            or "SYMPREC = 1.00e-04" not in source_incar or "8 8 6" not in source_kpoints):
        raise ValueError("1000 eV diagnostic VASP settings changed")

    case_dir = work_root / "cases" / "newton_followup_full"
    case_dir.mkdir(parents=True)
    write(case_dir / "POSCAR", candidate, format="vasp", direct=True, sort=False)
    if not same_geometry(candidate, read(case_dir / "POSCAR", format="vasp")):
        raise ValueError("staged POSCAR differs from proposed geometry")
    for name in ("INCAR", "KPOINTS", "POTCAR"):
        if sha256(completed_case_dir / name) != case["input_sha256"][name]:
            raise ValueError(f"completed source {name} changed")
        shutil.copy2(completed_case_dir / name, case_dir / name)
    hashes = {name: sha256(case_dir / name) for name in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
    (case_dir / "sha256.inputs.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf-8")
    result = {
        "purpose": "GaN_image15_followup_1000eV_local_stationarity_probe_not_TS_certificate",
        "status": "inputs_finalized_no_DFT", "pressure_GPa": 45.7,
        "evaluation_encut_eV": 1000, "source_hessian_encut_eV": 600,
        "current_gradient_eV_per_A": physical_gradient,
        "correction_joint_norm_A": float(np.linalg.norm(correction)),
        "correction_A": correction.tolist(),
        "predicted_full_gradient_eV_per_A": float(np.linalg.norm(
            eigenvectors.T @ gradient + eigenvalues * (eigenvectors.T @ correction)
        )),
        "empirical_bravais_risk": {
            "source_skew_A": cross_axis_skew_angstrom(source.cell.array),
            "candidate_skew_A": skew, "allowed_skew_A": risk_limit,
            "warning": "Case-local bracket is not a guarantee VASP will initialize.",
        },
        "cases": [{"name": "newton_followup_full", "input_sha256": hashes,
                   "minimum_distance_A": float(np.min(distances)),
                   "volume_A3": float(candidate.get_volume())}],
        "source_sha256": {
            "center_POSCAR": sha256(center_poscar),
            "newton_manifest": sha256(newton_manifest),
            "bracketed_manifest": sha256(bracketed_manifest),
            "bracketed_audit": sha256(bracketed_audit),
            "completed_OUTCAR": sha256(completed_case_dir / "OUTCAR"),
            "hessian_archive": sha256(hessian_archive),
            "preparer": sha256(Path(__file__)),
        },
        "limitations": "Approximate 600 eV Hessian and one 1000 eV static; not a TS certificate.",
    }
    (work_root / "manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("center-poscar", "newton-manifest", "bracketed-manifest",
                 "bracketed-audit", "completed-case-dir", "hessian-archive", "work-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(**vars(args))
    print(json.dumps({"status": result["status"],
                      "correction_A": result["correction_joint_norm_A"],
                      "candidate_skew_A": result["empirical_bravais_risk"]["candidate_skew_A"]}))


if __name__ == "__main__":
    main()
