"""Raw-audit the bounded GaN VASP stationary-point trial and apply its fixed gates.

The result is a decision on this *trial*, not by itself a full Hessian or a
transition-state certificate. Licensed POTCAR bytes are never exported.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from vcneb.joint_curvature import JointCurvatureCoordinates


EV_A3_TO_KBAR = 10.0 / GPa


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _same_geometry(a, b, tolerance=3e-5) -> bool:
    if a.get_chemical_symbols() != b.get_chemical_symbols():
        return False
    if np.max(np.abs(a.cell.array - b.cell.array)) > tolerance:
        return False
    d = a.get_scaled_positions(wrap=False) - b.get_scaled_positions(wrap=False)
    d -= np.rint(d)
    return bool(np.max(np.linalg.norm(d @ b.cell.array, axis=1)) < tolerance)


def _read_case(root: Path, record: dict, pressure: float) -> dict:
    directory = root / "cases" / record["name"]
    for name, digest in record["input_sha256"].items():
        if sha256(directory / name) != digest:
            raise ValueError(f"{record['name']}: {name} input hash changed")
    outcar = directory / "OUTCAR"
    content = outcar.read_text(encoding="utf-8", errors="replace")
    if ("aborting loop because EDIFF is reached" not in content
            or "General timing and accounting informations for this job" not in content):
        raise ValueError(f"{record['name']}: incomplete electronic/output footer")
    raw = read(outcar)
    if not _same_geometry(raw, read(directory / "POSCAR")):
        raise ValueError(f"{record['name']}: VASP changed static geometry")
    energy = float(raw.get_potential_energy())
    free_energy = float(raw.get_potential_energy(force_consistent=True))
    forces = np.asarray(raw.get_forces(), dtype=float)
    stress = np.asarray(raw.get_stress(voigt=False), dtype=float)
    enthalpy = energy + pressure * raw.get_volume()
    if (not np.isfinite([energy, free_energy, enthalpy]).all()
            or not np.isfinite(forces).all() or not np.isfinite(stress).all()
            or forces.shape != (4, 3) or stress.shape != (3, 3)):
        raise ValueError(f"{record['name']}: missing energy/force/stress")
    return {
        "name": record["name"],
        "atoms": raw,
        "energy_eV_cell": energy,
        "free_energy_eV_cell": free_energy,
        "enthalpy_eV_cell": float(enthalpy),
        "volume_A3": float(raw.get_volume()),
        "forces_eV_A": forces,
        "stress_eV_A3": stress,
        "OUTCAR_sha256": sha256(outcar),
    }


def audit(root: Path) -> dict:
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "predeclared_stationarity_trial_inputs_ready_no_DFT"
            or manifest.get("n_cases") != 9 or manifest.get("pressure_GPa") != 45.7):
        raise ValueError("invalid predeclared trial manifest")
    pressure = manifest["pressure_GPa"] * GPa
    records = [_read_case(root, entry, pressure) for entry in manifest["cases"]]
    data = {entry["name"]: entry for entry in records}
    center = data["center"]
    chart = JointCurvatureCoordinates(center["atoms"], manifest["cell_scale_A"])

    def gradient(entry: dict) -> np.ndarray:
        return chart.enthalpy_gradient(entry["atoms"], entry["forces_eV_A"],
                                       entry["stress_eV_A3"], pressure)

    center_gradient = gradient(center)
    force_max = float(np.linalg.norm(center["forces_eV_A"], axis=1).max())
    stress_residual = center["stress_eV_A3"] + np.eye(3) * pressure
    stress_kbar = float(np.max(np.abs(stress_residual)) * EV_A3_TO_KBAR)
    normal_rows = []
    for axis in (12, 13, 14):
        minus = data[f"normal{axis:02d}_minus"]
        plus = data[f"normal{axis:02d}_plus"]
        step = manifest["normal_step_A"]
        secant = (plus["enthalpy_eV_cell"] - minus["enthalpy_eV_cell"]) / (2 * step)
        integrated = (gradient(minus)[axis] + 4 * center_gradient[axis]
                      + gradient(plus)[axis]) / 6
        dvolume = (plus["volume_A3"] - minus["volume_A3"]) / (2 * step)
        if abs(dvolume) < 1e-6:
            raise ValueError(f"axis {axis}: vanishing volume derivative")
        normal_rows.append({
            "axis": axis,
            "enthalpy_secant_eV_per_A": float(secant),
            "stress_gradient_Simpson_eV_per_A": float(integrated),
            "dV_dq_A2": float(dvolume),
            "energy_derived_pressure_residual_kbar": float(abs(secant / dvolume) * EV_A3_TO_KBAR),
            "energy_stress_mismatch_kbar": float(abs((secant - integrated) / dvolume) * EV_A3_TO_KBAR),
        })
    unstable_minus = data["unstable_minus"]
    unstable_plus = data["unstable_plus"]
    u_step = manifest["unstable_step_A"]
    unstable_secant = (unstable_plus["enthalpy_eV_cell"]
                        - unstable_minus["enthalpy_eV_cell"]) / (2 * u_step)
    unstable_curvature = (unstable_plus["enthalpy_eV_cell"]
                          - 2 * center["enthalpy_eV_cell"]
                          + unstable_minus["enthalpy_eV_cell"]) / u_step ** 2
    acceptance = manifest["acceptance"]
    max_energy_kbar = max(row["energy_derived_pressure_residual_kbar"]
                          for row in normal_rows)
    max_mismatch_kbar = max(row["energy_stress_mismatch_kbar"] for row in normal_rows)
    gates = {
        "raw_SCF_complete": True,
        "center_atomic_force": force_max <= acceptance["max_atomic_force_eV_per_A"],
        "center_stress": stress_kbar <= acceptance["max_hydrostatic_and_normal_stress_residual_kbar"],
        "normal_energy_stationarity": max_energy_kbar <= acceptance["max_energy_derived_normal_pressure_residual_kbar"],
        "normal_energy_stress_consistency": max_mismatch_kbar <= acceptance["normal_energy_stress_consistency_kbar"],
        "unstable_energy_curvature": unstable_curvature < 0,
    }
    case_reports = []
    for entry in records:
        case_reports.append({
            "case": entry["name"],
            "energy_eV_cell": entry["energy_eV_cell"],
            "free_energy_eV_cell": entry["free_energy_eV_cell"],
            "enthalpy_eV_cell": entry["enthalpy_eV_cell"],
            "volume_A3": entry["volume_A3"],
            "max_atomic_force_eV_per_A": float(np.linalg.norm(entry["forces_eV_A"], axis=1).max()),
            "OUTCAR_sha256": entry["OUTCAR_sha256"],
        })
    return {
        "status": "stationarity_trial_raw_audited",
        "all_predeclared_trial_gates_pass": bool(all(gates.values())),
        "claim_limit": ("Even if all current gates pass, a fresh complete 15D Hessian and basin links are required for a final index-one TS certificate."),
        "pressure_GPa": manifest["pressure_GPa"],
        "center_max_atomic_force_eV_per_A": force_max,
        "center_stress_residual_max_kbar": stress_kbar,
        "center_stress_residual_tensor_kbar": (stress_residual * EV_A3_TO_KBAR).tolist(),
        "center_stress_force_gradient_norm_eV_per_A": float(np.linalg.norm(center_gradient)),
        "max_normal_energy_derived_pressure_residual_kbar": max_energy_kbar,
        "max_normal_energy_stress_mismatch_kbar": max_mismatch_kbar,
        "unstable_secant_eV_per_A": float(unstable_secant),
        "unstable_secant_within_0p01_eV_per_A_diagnostic": bool(abs(unstable_secant) <= 0.01),
        "unstable_energy_curvature_eV_per_A2": float(unstable_curvature),
        "free_energy_vs_energy_max_abs_eV_cell": max(abs(x["free_energy_eV_cell"] - x["energy_eV_cell"])
                                                  for x in records),
        "normal_rows": normal_rows,
        "gates": gates,
        "cases": case_reports,
        "manifest_sha256": sha256(manifest_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.case)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "all_predeclared_trial_gates_pass", "center_stress_residual_max_kbar",
        "max_normal_energy_derived_pressure_residual_kbar",
        "max_normal_energy_stress_mismatch_kbar",
        "unstable_energy_curvature_eV_per_A2")}))


if __name__ == "__main__":
    main()
