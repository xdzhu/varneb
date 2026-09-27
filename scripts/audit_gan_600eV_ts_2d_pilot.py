"""Audit the GaN 600-eV frozen two-joint-mode enthalpy pilot and holdouts."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def centered_curvature(plus: float, zero: float, minus: float, step: float) -> float:
    if not np.isfinite([plus, zero, minus, step]).all() or step <= 0:
        raise ValueError("nonfinite energies or nonpositive step")
    return float((plus + minus - 2 * zero) / step**2)


def audit(work: Path, center_work: Path, hessian_root: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    hessian_audit_path = hessian_root / "audit.json"
    hessian_npz = hessian_root / "joint_hessian.npz"
    hessian_audit = json.loads(hessian_audit_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_45p7_600eV_local_joint_mode_2D_frozen_pilot"
            or len(manifest.get("cases", [])) != 12
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or manifest.get("source_sha256", {}).get("center_OUTCAR") != sha256(center_work / "OUTCAR")
            or manifest.get("source_sha256", {}).get("hessian_audit") != sha256(hessian_audit_path)
            or manifest.get("source_sha256", {}).get("hessian_npz") != sha256(hessian_npz)
            or manifest.get("source_sha256", {}).get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_ts_2d_pilot.py"))
            or hessian_audit.get("negative_count_below_minus_0p1_eV_per_A2") != 1):
        raise ValueError("2D pilot provenance or 600-eV input contract failed")
    center = read(center_work / "OUTCAR")
    chart = JointCurvatureCoordinates(center, float(manifest["axes"]["cell_scale_A"]))
    pressure = 45.7 * GPa
    h0 = float(center.get_potential_energy() + pressure * center.get_volume())
    with np.load(hessian_npz, allow_pickle=False) as archive:
        eig = np.asarray(archive["eigenvalues"], dtype=float)
        modes = np.asarray(archive["eigenvectors"], dtype=float)
    if (not np.allclose(eig[:2], manifest["axes"]["eigenvalues_eV_per_A2"], atol=1e-8)
            or modes.shape != (18, 15)
            or not np.allclose(modes.T @ modes, np.eye(15), atol=1e-8)):
        raise ValueError("2D axes changed since the full Hessian audit")
    center_gradient = chart.enthalpy_gradient(
        center, center.get_forces(), center.get_stress(voigt=False), pressure
    )
    g_uv = modes[:, :2].T @ center_gradient
    cases = []
    energy_by_name = {"center": h0}
    for record in manifest["cases"]:
        if any(record.get("input_sha256", {}).get(name) != digest
               for name, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"case {record['name']} changed electronic contract")
        parsed, detail = completed_case(
            work / "cases" / record["name"],
            {**record, "POSCAR_sha256": record["input_sha256"]["POSCAR"]},
        )
        h = float(parsed.get_potential_energy() + pressure * parsed.get_volume())
        qu = float(record["q_u_A"])
        qv = float(record["q_v_A"])
        model_h = h0 + g_uv[0] * qu + g_uv[1] * qv + 0.5 * (eig[0] * qu**2 + eig[1] * qv**2)
        detail.update({
            "q_u_A": qu, "q_v_A": qv,
            "enthalpy_eV_per_cell": h,
            "delta_enthalpy_meV_per_GaN": float((h - h0) * 500),
            "quadratic_model_residual_meV_per_GaN": float((h - model_h) * 500),
        })
        cases.append(detail)
        energy_by_name[record["name"]] = h
    emap = energy_by_name
    du = float(manifest["grid_q_u_A"][2])
    dv = float(manifest["grid_q_v_A"][2])
    curvatures = {
        "u_grid_eV_per_A2": centered_curvature(emap["grid_up_vz"], h0, emap["grid_um_vz"], du),
        "u_half_eV_per_A2": centered_curvature(emap["hold_up"], h0, emap["hold_um"], du / 2),
        "v_grid_eV_per_A2": centered_curvature(emap["grid_uz_vp"], h0, emap["grid_uz_vm"], dv),
        "v_half_eV_per_A2": centered_curvature(emap["hold_vp"], h0, emap["hold_vm"], dv / 2),
        "uv_mixed_grid_eV_per_A2": float((emap["grid_up_vp"] - emap["grid_up_vm"]
                                          - emap["grid_um_vp"] + emap["grid_um_vm"])
                                         / (4 * du * dv)),
    }
    grid = [case for case in cases if case["case"].startswith("grid_")]
    holdouts = [case for case in cases if case["case"].startswith("hold_")]
    max_grid = float(max(abs(case["quadratic_model_residual_meV_per_GaN"]) for case in grid))
    max_holdout = float(max(abs(case["quadratic_model_residual_meV_per_GaN"]) for case in holdouts))
    gate = {
        "negative_u_positive_v_at_both_scales": bool(
            curvatures["u_grid_eV_per_A2"] < 0
            and curvatures["u_half_eV_per_A2"] < 0
            and curvatures["v_grid_eV_per_A2"] > 0
            and curvatures["v_half_eV_per_A2"] > 0
        ),
        "holdout_model_error_below_0p2_meV_per_GaN": max_holdout <= 0.2,
    }
    return {
        "status": "GaN_600eV_local_2D_frozen_enthalpy_pilot_raw_audited",
        "claim_limit": "Local frozen enthalpy cut; not a full-path PES or a certified index-one TS",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "center_enthalpy_eV_per_cell": h0,
        "center_gradient_uv_eV_per_A": g_uv.tolist(),
        "hessian_eigenvalues_uv_eV_per_A2": eig[:2].tolist(),
        "energy_curvatures": curvatures,
        "quadratic_model_max_grid_error_meV_per_GaN": max_grid,
        "quadratic_model_max_axial_holdout_error_meV_per_GaN": max_holdout,
        "predeclared_gate": gate,
        "cases": cases,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "center_OUTCAR": sha256(center_work / "OUTCAR"),
            "hessian_audit": sha256(hessian_audit_path),
            "hessian_npz": sha256(hessian_npz),
            "auditor": sha256(Path(__file__)),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "center-work", "hessian-root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work, args.center_work, args.hessian_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "gate": result["predeclared_gate"]}))


if __name__ == "__main__":
    main()
