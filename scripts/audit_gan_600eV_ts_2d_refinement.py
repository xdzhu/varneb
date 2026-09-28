"""Audit twelve new VASP statics and join the original GaN 13-point cut.

This is a frozen local enthalpy surface at the production 600-eV/45.7-GPa
contract; it does not certify a whole-path PES or a stationary TS.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from scripts.prepare_gan_600eV_ts_2d_refinement import coordinate_key


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(work: Path, regenerated_manifest_path: Path, pilot_audit_path: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    regenerated = json.loads(regenerated_manifest_path.read_text(encoding="utf-8"))
    pilot = json.loads(pilot_audit_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_45p7_600eV_local_joint_mode_5x5_frozen_refinement"
            or manifest.get("pressure_GPa") != 45.7
            or manifest.get("formula_units_per_cell") != 2
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or manifest.get("source_sha256", {}).get("pilot_audit") != sha256(pilot_audit_path)
            or regenerated.get("purpose") != manifest["purpose"]
            or len(manifest.get("cases", [])) != 12
            or len(regenerated.get("cases", [])) != 12
            or pilot.get("status") != "GaN_600eV_local_2D_frozen_enthalpy_pilot_raw_audited"
            or len(pilot.get("cases", [])) != 12):
        raise ValueError("joint 5x5 refinement does not join its original 600-eV pilot")
    prior = {coordinate_key(float(case["q_u_A"]), float(case["q_v_A"]))
             for case in pilot["cases"]}
    prior.add((0.0, 0.0))
    if len(prior) != 13:
        raise ValueError("original center and twelve pilot points are not unique")
    pressure = 45.7 * GPa
    h0 = float(pilot["center_enthalpy_eV_per_cell"])
    gradient = np.asarray(pilot["center_gradient_uv_eV_per_A"], dtype=float)
    eig = np.asarray(pilot["hessian_eigenvalues_uv_eV_per_A2"], dtype=float)
    results = []
    fresh = set()
    for item, repeat in zip(manifest["cases"], regenerated["cases"]):
        if (item["name"] != repeat["name"]
                or item["grid_index"] != repeat["grid_index"]
                or item["q_u_A"] != repeat["q_u_A"]
                or item["q_v_A"] != repeat["q_v_A"]
                or item["input_sha256"] != repeat["input_sha256"]):
            raise ValueError("current preparer does not regenerate identical VASP inputs")
        if any(item["input_sha256"].get(name) != digest
               for name, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError("VASP electronic contract differs from original 600 eV")
        q_u, q_v = float(item["q_u_A"]), float(item["q_v_A"])
        key = coordinate_key(q_u, q_v)
        if key in prior or key in fresh:
            raise ValueError("duplicate or old GaN joint-grid coordinate")
        fresh.add(key)
        atoms, detail = completed_case(
            work / "cases" / item["name"],
            {**item, "POSCAR_sha256": item["input_sha256"]["POSCAR"]},
        )
        h = float(atoms.get_potential_energy() + pressure * atoms.get_volume())
        model = h0 + gradient[0] * q_u + gradient[1] * q_v + 0.5 * (
            eig[0] * q_u**2 + eig[1] * q_v**2
        )
        detail.update({
            "grid_index": item["grid_index"], "q_u_A": q_u, "q_v_A": q_v,
            "enthalpy_eV_per_cell": h,
            "delta_enthalpy_meV_per_GaN": float((h - h0) * 500),
            "quadratic_model_residual_meV_per_GaN": float((h - model) * 500),
        })
        results.append(detail)
    full = {coordinate_key(float(q_u), float(q_v))
            for q_u in manifest["q_u_grid_A"] for q_v in manifest["q_v_grid_A"]}
    if len(results) != 12 or len(full) != 25 or prior | fresh != full:
        raise ValueError("original and new GaN points do not complete one 5x5 grid")
    errors = np.asarray([item["quadratic_model_residual_meV_per_GaN"]
                         for item in results], dtype=float)
    if not np.isfinite(errors).all():
        raise ValueError("nonfinite new-point residual")
    return {
        "status": "GaN_600eV_local_joint_5x5_refinement_raw_audited",
        "pressure_GPa": 45.7, "formula_units_per_cell": 2,
        "n_old_measured_coordinates": 13, "n_new_raw_audited_statics": 12,
        "n_total_measured_coordinates": 25,
        "new_point_quadratic_model_max_abs_error_meV_per_GaN": float(np.max(np.abs(errors))),
        "new_point_quadratic_model_rms_error_meV_per_GaN": float(np.sqrt(np.mean(errors**2))),
        "cases": results,
        "source_sha256": {
            "submitted_manifest": sha256(manifest_path),
            "current_source_regenerated_manifest": sha256(regenerated_manifest_path),
            "original_pilot_audit": sha256(pilot_audit_path),
            "auditor": sha256(Path(__file__)),
        },
        "claim_limit": (
            "A frozen local joint atom-strain enthalpy cut; these extra statics do not "
            "prove orthogonal relaxation, whole-path coverage, or strict TS stationarity."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--regenerated-manifest", type=Path, required=True)
    parser.add_argument("--pilot-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work, args.regenerated_manifest, args.pilot_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "n_total_measured_coordinates",
        "new_point_quadratic_model_max_abs_error_meV_per_GaN",
    )}))


if __name__ == "__main__":
    main()
