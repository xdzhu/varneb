"""Raw-audit the 72 prospective GaN q nodes before extending the 2D chart."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_600eV_atomic_tube_q9 import IMAGE_INDICES, Q_FULL, Q_NEW
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(work: Path, dense_report: Path, dense_manifest: Path,
          trajectory: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    previous = json.loads(dense_report.read_text(encoding="utf-8"))
    old_manifest = json.loads(dense_manifest.read_text(encoding="utf-8"))
    if (manifest.get("status") != "72_interleaved_inputs_finalized_no_DFT"
            or manifest.get("purpose")
            != "GaN_45p7_600eV_central_18x9_frozen_atomic_transverse_cut"
            or manifest.get("n_old_measured_points") != 90
            or manifest.get("n_new_statics") != 72
            or manifest.get("central_image_indices") != list(IMAGE_INDICES)
            or manifest.get("q_atom_A_full") != list(Q_FULL)
            or manifest.get("q_atom_A_new") != list(Q_NEW)
            or manifest.get("prospective_error_gate_meV_per_GaN") != 1.0
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or manifest.get("source_sha256", {}).get("dense_report") != sha256(dense_report)
            or manifest.get("source_sha256", {}).get("dense_manifest") != sha256(dense_manifest)
            or manifest.get("source_sha256", {}).get("trajectory") != sha256(trajectory)
            or previous.get("status")
            != "GaN_600eV_central_dense_atomic_tube_raw_audited"
            or previous.get("source_sha256", {}).get("manifest") != sha256(dense_manifest)
            or previous.get("central_image_indices") != list(IMAGE_INDICES)
            or previous.get("q_atom_A") != old_manifest.get("q_atom_A")
            or len(manifest.get("cases", [])) != 72):
        raise ValueError("GaN 18x9 source, input, or pressure contract changed")
    frames = read(trajectory, index=":")
    if len(frames) != 29:
        raise ValueError("expected 29 archived path images")
    pressure = 45.7 * GPa
    old = np.asarray(previous["excess_enthalpy_meV_per_GaN"], dtype=float)
    centers = np.asarray(previous["path_enthalpy_eV_per_cell"], dtype=float)
    if (old.shape != (18, 5) or centers.shape != (18,)
            or not np.isfinite(old).all() or not np.isfinite(centers).all()
            or any(not np.isclose(
                frames[index].get_potential_energy() + pressure * frames[index].get_volume(),
                centers[row], atol=1e-8, rtol=0,
            ) for row, index in enumerate(IMAGE_INDICES))):
        raise ValueError("audited 600-eV path-center enthalpies changed")
    matrix = {index: {round(float(q), 6): float(old[row, column])
                      for column, q in enumerate(previous["q_atom_A"])}
              for row, index in enumerate(IMAGE_INDICES)}
    expected = {(index, round(q, 6)) for index in IMAGE_INDICES for q in Q_NEW}
    observed_keys = set()
    records = []
    for case in manifest["cases"]:
        index, q = int(case["image_index"]), round(float(case["q_atom_A"]), 6)
        key = (index, q)
        if (key not in expected or key in observed_keys
                or case.get("role") != "prospective_interleaved_q_holdout"
                or any(case.get("input_sha256", {}).get(name) != digest
                       for name, digest in PRODUCTION_INPUT_SHA256.items())):
            raise ValueError(f"duplicate, unexpected, or nonproduction case: {case['name']}")
        observed_keys.add(key)
        atoms, detail = completed_case(
            work / "cases" / case["name"],
            {**case, "POSCAR_sha256": case["input_sha256"]["POSCAR"]},
        )
        if not np.allclose(atoms.cell.array, frames[index].cell.array,
                           atol=2e-5, rtol=0):
            raise ValueError(f"off-path cell changed: {case['name']}")
        enthalpy = float(atoms.get_potential_energy() + pressure * atoms.get_volume())
        excess = (enthalpy - centers[index - IMAGE_INDICES[0]]) * 500.0
        prediction = float(case["prediction_excess_meV_per_GaN"])
        if not np.isfinite([enthalpy, excess, prediction]).all():
            raise ValueError(f"nonfinite measurement or frozen prediction: {case['name']}")
        matrix[index][q] = excess
        detail.update({
            "image_index": index, "q_atom_A": q,
            "enthalpy_eV_per_cell": enthalpy,
            "excess_enthalpy_meV_per_GaN": excess,
            "prospective_prediction_meV_per_GaN": prediction,
            "prediction_minus_DFT_meV_per_GaN": prediction - excess,
        })
        records.append(detail)
    if observed_keys != expected or any(set(matrix[index]) != set(Q_FULL)
                                         for index in IMAGE_INDICES):
        raise ValueError("18x9 GaN grid is incomplete")
    errors = np.array([row["prediction_minus_DFT_meV_per_GaN"] for row in records])
    max_error = float(np.max(np.abs(errors)))
    return {
        "status": "GaN_600eV_central_18x9_raw_audited",
        "claim_limit": "45.7-GPa frozen central (s,q_atom) cut; not a whole-path, relaxed, or two-phonon PES",
        "pressure_GPa": 45.7, "formula_units_per_cell": 2,
        "n_reused_measured_points": 90, "n_new_raw_audited_statics": 72,
        "central_image_indices": list(IMAGE_INDICES),
        "arc_fraction_s": previous["arc_fraction_s"],
        "path_enthalpy_eV_per_cell": centers.tolist(),
        "q_atom_A": list(Q_FULL),
        "excess_enthalpy_meV_per_GaN": [
            [matrix[index][round(q, 6)] for q in Q_FULL] for index in IMAGE_INDICES
        ],
        "prospective_max_abs_error_meV_per_GaN": max_error,
        "prospective_rms_error_meV_per_GaN": float(np.sqrt(np.mean(errors**2))),
        "predeclared_error_gate_meV_per_GaN": 1.0,
        "prospective_gate_pass": max_error <= 1.0,
        "new_cases": records,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "dense_report": sha256(dense_report),
            "dense_manifest": sha256(dense_manifest),
            "trajectory": sha256(trajectory),
            "auditor": sha256(Path(__file__)),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "dense-report", "dense-manifest", "trajectory", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.work, args.dense_report, args.dense_manifest,
                   args.trajectory)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "status", "n_new_raw_audited_statics",
        "prospective_max_abs_error_meV_per_GaN", "prospective_gate_pass",
    )}))


if __name__ == "__main__":
    main()
