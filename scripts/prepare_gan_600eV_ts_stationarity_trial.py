"""Prepare one bounded, same-contract GaN saddle-stationarity trial.

The predeclared question is whether a small energy-informed correction to the
fixed VCNEB center can meet both the VASP force/stress and independent
energy-secant gates. A failed gate must not be relabeled a certified TS.
No electronic setting, external pressure, or source path is changed here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from vcneb.joint_curvature import JointCurvatureCoordinates


CENTER_OUTCAR_SHA256 = "66e9ddf2bf5ee1f3f510bf71b8ec8dea2e05ae3d8be60e9593e741a87dabcfe5"
HESSIAN_SHA256 = "730aa496c937e9089c84d0aa43d15238d4d05ed1b28ac0084a7305c6f5303e8f"
INPUT_SHA256 = {
    "INCAR": "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772",
    "KPOINTS": "b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee",
    "POTCAR": "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168",
}
CELL_SCALE_A = 3.3982714330050063
PRESSURE_GPA = 45.7
NORMAL_STEP_A = 0.005
UNSTABLE_STEP_A = 0.01


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hybrid_newton_step(hessian_npz: Path, normal_report: Path) -> dict:
    """Use measured small-step normal energy slopes and center force/stress elsewhere."""
    if sha256(hessian_npz) != HESSIAN_SHA256:
        raise ValueError("local joint Hessian basis changed")
    report = json.loads(normal_report.read_text(encoding="utf-8"))
    if (report.get("status") != "GaN_600eV_normal_strain_step_dependence_raw_audited"
            or report.get("pressure_GPa") != PRESSURE_GPA
            or report.get("source_sha256", {}).get("center_OUTCAR") != CENTER_OUTCAR_SHA256):
        raise ValueError("normal-strain report lacks the original raw audit")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        eig = np.asarray(archive["eigenvalues"], dtype=float)
        modes = np.asarray(archive["eigenvectors"], dtype=float)
        gradient = np.asarray(archive["center_gradient"], dtype=float)
    if (eig.shape != (15,) or modes.shape != (18, 15)
            or gradient.shape != (18,) or not eig[0] < 0 < eig[1]
            or not np.allclose(modes.T @ modes, np.eye(15), atol=1e-8)):
        raise ValueError("invalid translation-free index-one local basis")
    hybrid = gradient.copy()
    normal = {}
    for axis in (12, 13, 14):
        candidates = [row for row in report["rows"]
                      if row["axis"] == axis and row["step_A"] == NORMAL_STEP_A]
        if len(candidates) != 1:
            raise ValueError(f"missing 0.005-A normal-strain secant for axis {axis}")
        hybrid[axis] = float(candidates[0]["energy_enthalpy_secant_eV_per_A"])
        normal[str(axis)] = hybrid[axis]
    modal_gradient = modes.T @ hybrid
    shift = -modes @ (modal_gradient / eig)
    if (not np.isfinite(shift).all() or np.linalg.norm(shift) > 0.004
            or np.max(np.abs(shift[:12])) > 0.0015
            or np.max(np.abs(shift[12:])) > 0.0025):
        raise ValueError("one-step correction leaves predeclared trust radius")
    return {
        "center_stress_force_gradient_eV_per_A": gradient.tolist(),
        "hybrid_energy_informed_gradient_eV_per_A": hybrid.tolist(),
        "normal_energy_secants_eV_per_A": normal,
        "eigenvalues_eV_per_A2": eig.tolist(),
        "modal_gradient_eV_per_A": modal_gradient.tolist(),
        "newton_shift_A": shift.tolist(),
        "newton_shift_norm_A": float(np.linalg.norm(shift)),
        "unstable_mode": modes[:, 0].tolist(),
    }


def _same_geometry(a, b, tolerance: float = 2e-5) -> bool:
    if a.get_chemical_symbols() != b.get_chemical_symbols():
        return False
    if np.max(np.abs(a.cell.array - b.cell.array)) > tolerance:
        return False
    d = a.get_scaled_positions(wrap=False) - b.get_scaled_positions(wrap=False)
    d -= np.rint(d)
    return bool(np.max(np.linalg.norm(d @ b.cell.array, axis=1)) < tolerance)


def prepare(source: Path, hessian_npz: Path, normal_report: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    if sha256(source / "OUTCAR") != CENTER_OUTCAR_SHA256:
        raise ValueError("source center OUTCAR changed")
    if any(sha256(source / name) != digest for name, digest in INPUT_SHA256.items()):
        raise ValueError("source VASP input contract changed")
    outcar = (source / "OUTCAR").read_text(encoding="utf-8", errors="replace")
    if ("aborting loop because EDIFF is reached" not in outcar
            or "General timing and accounting informations for this job" not in outcar):
        raise ValueError("source center electronic calculation incomplete")
    center = read(source / "OUTCAR")
    if center.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]:
        raise ValueError("unexpected center composition/order")
    correction = hybrid_newton_step(hessian_npz, normal_report)
    original_chart = JointCurvatureCoordinates(center, CELL_SCALE_A)
    shifted = original_chart.displaced(np.asarray(correction["newton_shift_A"]))
    local_chart = JointCurvatureCoordinates(shifted, CELL_SCALE_A)
    if shifted.get_volume() <= 0:
        raise ValueError("invalid shifted center")
    directions = [("center", np.zeros(local_chart.size))]
    for axis in (12, 13, 14):
        for sign, label in ((-1, "minus"), (1, "plus")):
            displacement = np.zeros(local_chart.size)
            displacement[axis] = sign * NORMAL_STEP_A
            directions.append((f"normal{axis:02d}_{label}", displacement))
    mode = np.asarray(correction["unstable_mode"], dtype=float)
    for sign, label in ((-1, "minus"), (1, "plus")):
        directions.append((f"unstable_{label}", sign * UNSTABLE_STEP_A * mode))
    output.mkdir(parents=True)
    cases_root = output / "cases"
    cases_root.mkdir()
    cases = []
    for name, displacement in directions:
        atoms = local_chart.displaced(displacement)
        if (atoms.get_volume() <= 0 or atoms.get_chemical_symbols() != center.get_chemical_symbols()
                or not np.isfinite(atoms.positions).all()):
            raise ValueError(f"invalid geometry for {name}")
        case = cases_root / name
        case.mkdir()
        write(case / "POSCAR", atoms, format="vasp", direct=True, sort=False)
        if not _same_geometry(atoms, read(case / "POSCAR", format="vasp")):
            raise ValueError(f"POSCAR roundtrip changed {name}")
        for file in INPUT_SHA256:
            shutil.copy2(source / file, case / file)
        hashes = {file: sha256(case / file) for file in ("POSCAR", *INPUT_SHA256)}
        if any(hashes[file] != INPUT_SHA256[file] for file in INPUT_SHA256):
            raise ValueError(f"input contract changed in {name}")
        cases.append({"name": name, "joint_displacement_A": displacement.tolist(),
                      "volume_A3": float(atoms.get_volume()), "input_sha256": hashes})
    manifest = {
        "status": "predeclared_stationarity_trial_inputs_ready_no_DFT",
        "purpose": "test whether a nearby energy-informed center satisfies independent energy and VASP stress stationarity gates",
        "pressure_GPa": PRESSURE_GPA,
        "electronic_contract": "VASP6.3.2/PBE/Ga_d+N/600eV/Gamma8x8x6/EDIFF1e-7/ISYM-1/SYMPREC1e-4",
        "input_sha256": INPUT_SHA256,
        "source_center_OUTCAR_sha256": CENTER_OUTCAR_SHA256,
        "hessian_npz_sha256": HESSIAN_SHA256,
        "normal_report_sha256": sha256(normal_report),
        "preparer_sha256": sha256(Path(__file__)),
        "cell_scale_A": CELL_SCALE_A,
        "normal_step_A": NORMAL_STEP_A,
        "unstable_step_A": UNSTABLE_STEP_A,
        "n_cases": len(cases),
        "correction": correction,
        "cases": cases,
        "acceptance": {
            "raw_SCF_complete": True,
            "max_atomic_force_eV_per_A": 0.01,
            "max_hydrostatic_and_normal_stress_residual_kbar": 2.0,
            "max_energy_derived_normal_pressure_residual_kbar": 2.0,
            "normal_energy_stress_consistency_kbar": 2.0,
            "unstable_direction_energy_curvature_negative": True,
            "full_joint_Hessian_index_one_required_for_final_TS_certificate": True,
        },
        "claim_limit": "This bounded trial alone cannot certify an index-one TS; passing stationarity gates would trigger a fresh full Hessian and basin audit, while failure must be reported as non-certification.",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("source", "hessian-npz", "normal-report", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.source, args.hessian_npz, args.normal_report, args.output)
    print(json.dumps({"status": result["status"], "n_cases": result["n_cases"],
                      "newton_shift_norm_A": result["correction"]["newton_shift_norm_A"]}))


if __name__ == "__main__":
    main()
