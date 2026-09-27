"""Audit GaN central atomic-tube refinement and independent q half-steps."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_600eV_path_tube import leave_one_anchor_out_errors
from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_600eV_atomic_tube_refinement import (
    HALFSTEP_HOLDOUTS,
    NEW_ANCHORS,
)
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(work: Path, trajectory: Path, atomic_feasibility_path: Path,
          joint_feasibility_path: Path, transport_path: Path,
          first_audit_path: Path, hessian_npz: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    atomic_feasibility = json.loads(atomic_feasibility_path.read_text(encoding="utf-8"))
    joint_feasibility = json.loads(joint_feasibility_path.read_text(encoding="utf-8"))
    transport = json.loads(transport_path.read_text(encoding="utf-8"))
    first = json.loads(first_audit_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose")
            != "GaN_45p7_600eV_atomic_only_central_tube_targeted_refinement"
            or manifest.get("new_anchors") != list(NEW_ANCHORS)
            or manifest.get("transverse_holdout_anchors") != list(HALFSTEP_HOLDOUTS)
            or manifest.get("anchor_q_A") != [-0.05, 0.05]
            or manifest.get("holdout_q_A") != [-0.025, 0.025]
            or manifest.get("predeclared_LOO_gate_meV_per_GaN") != 1.0
            or manifest.get("predeclared_halfstep_curvature_relative_gate") != 0.05
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or len(manifest.get("cases", [])) != 12
            or manifest["source_sha256"].get("trajectory") != sha256(trajectory)
            or manifest["source_sha256"].get("atomic_feasibility")
            != sha256(atomic_feasibility_path)
            or manifest["source_sha256"].get("joint_feasibility")
            != sha256(joint_feasibility_path)
            or manifest["source_sha256"].get("transport_comparison")
            != sha256(transport_path)
            or manifest["source_sha256"].get("first_audit") != sha256(first_audit_path)
            or manifest["source_sha256"].get("hessian_npz") != sha256(hessian_npz)
            or manifest["source_sha256"].get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_atomic_tube_refinement.py"))
            or first.get("status")
            != "GaN_600eV_central_atomic_transverse_tube_raw_audited"
            or first["source_sha256"].get("trajectory") != sha256(trajectory)
            or first["source_sha256"].get("atomic_feasibility")
            != sha256(atomic_feasibility_path)
            or atomic_feasibility.get("status")
            != "GaN_600eV_atomic_transverse_tube_geometry_only_no_DFT"
            or joint_feasibility.get("status")
            != "GaN_600eV_path_adapted_joint_normal_tube_geometry_only"
            or transport.get("status")
            != "GaN_600eV_atomic_transverse_transport_offline_comparison"):
        raise ValueError("refinement provenance or original 600-eV contract failed")
    frames = read(trajectory, index=":")
    pressure = 45.7 * GPa
    arc = first["arc_fraction_s"]
    h0 = first["path_enthalpy_eV_per_cell"]
    if (len(frames) != 29
            or not np.allclose(arc, joint_feasibility["arc_fraction_s"], atol=1e-12)
            or any(not np.isclose(
                float(frame.get_potential_energy() + pressure * frame.get_volume()),
                h0[index], atol=1e-8, rtol=0
            ) for index, frame in enumerate(frames))):
        raise ValueError("audited q=0 path changed")
    by_anchor: dict[int, dict[str, dict]] = {index: {} for index in first["anchors"]}
    for record in first["cases"]:
        index = int(record["image_index"])
        side = "plus" if record["q_atom_A"] > 0 else "minus"
        if side in by_anchor[index]:
            raise ValueError("duplicate original signed point")
        by_anchor[index][side] = record
    holdout_by_anchor: dict[int, dict[str, dict]] = {index: {} for index in HALFSTEP_HOLDOUTS}
    new_cases = []
    for record in manifest["cases"]:
        if any(record.get("input_sha256", {}).get(name) != digest
               for name, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"{record['name']} changed electronic inputs")
        atoms, detail = completed_case(
            work / "cases" / record["name"],
            {**record, "POSCAR_sha256": record["input_sha256"]["POSCAR"]},
        )
        index = int(record["image_index"])
        if not np.allclose(atoms.cell.array, frames[index].cell.array, atol=2e-5, rtol=0):
            raise ValueError(f"atomic-only cell changed at {record['name']}")
        q = float(record["q_atom_A"])
        side = "plus" if q > 0 else "minus"
        h = float(atoms.get_potential_energy() + pressure * atoms.get_volume())
        detail.update({
            "role": record["role"],
            "image_index": index,
            "arc_fraction_s": float(arc[index]),
            "q_atom_A": q,
            "enthalpy_eV_per_cell": h,
            "delta_enthalpy_meV_per_GaN": float((h - h0[index]) * 500),
        })
        bucket = by_anchor.setdefault(index, {}) if record["role"] == "interpolation_anchor" \
            else holdout_by_anchor[index]
        if side in bucket:
            raise ValueError(f"duplicate refinement side at {record['name']}")
        bucket[side] = detail
        new_cases.append(detail)
    anchors = sorted(by_anchor)
    if any(set(by_anchor[index]) != {"minus", "plus"} for index in anchors):
        raise ValueError("combined q=0.05 grid lacks a signed pair")
    holdouts = leave_one_anchor_out_errors(anchors, arc, by_anchor)
    max_loo = max(item["absolute_error_meV_per_GaN"] for item in holdouts)
    first_curvature = {item["image_index"]: item["frozen_atomic_transverse_curvature_eV_per_A2"]
                       for item in first["anchor_diagnostics"]}
    q_checks = []
    for index in HALFSTEP_HOLDOUTS:
        signed = holdout_by_anchor[index]
        if set(signed) != {"minus", "plus"}:
            raise ValueError(f"half-step q pair missing at image {index}")
        half_curvature = (
            signed["minus"]["enthalpy_eV_per_cell"]
            + signed["plus"]["enthalpy_eV_per_cell"]
            - 2 * h0[index]
        ) / 0.025**2
        full_curvature = float(first_curvature[index])
        relative = abs(half_curvature - full_curvature) / abs(full_curvature)
        q_checks.append({
            "image_index": index,
            "halfstep_curvature_eV_per_A2": float(half_curvature),
            "fullstep_curvature_eV_per_A2": full_curvature,
            "relative_difference": float(relative),
        })
    max_relative = max(item["relative_difference"] for item in q_checks)
    return {
        "status": "GaN_600eV_central_atomic_tube_targeted_refinement_raw_audited",
        "claim_limit": "Central path segment only; frozen atomic-only q, not globally continuous endpoint mode",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "n_path_images": 29,
        "n_original_offpath_cases": len(first["cases"]),
        "n_new_offpath_cases": len(new_cases),
        "anchors": anchors,
        "arc_fraction_s": arc,
        "path_enthalpy_eV_per_cell": h0,
        "anchor_q_A": [-0.05, 0.05],
        "holdout_q_A": [-0.025, 0.025],
        "leave_one_anchor_out_excess_enthalpy_errors": holdouts,
        "maximum_leave_one_out_error_meV_per_GaN": float(max_loo),
        "predeclared_LOO_gate_meV_per_GaN": 1.0,
        "path_interpolation_gate_pass": bool(max_loo <= 1.0),
        "halfstep_curvature_checks": q_checks,
        "maximum_halfstep_curvature_relative_difference": float(max_relative),
        "predeclared_halfstep_curvature_relative_gate": 0.05,
        "transverse_curvature_gate_pass": bool(max_relative <= 0.05),
        "central_contour_gate_pass": bool(max_loo <= 1.0 and max_relative <= 0.05),
        "new_cases": new_cases,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "trajectory": sha256(trajectory),
            "atomic_feasibility": sha256(atomic_feasibility_path),
            "joint_feasibility": sha256(joint_feasibility_path),
            "transport_comparison": sha256(transport_path),
            "first_audit": sha256(first_audit_path),
            "hessian_npz": sha256(hessian_npz),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "Path holdout tests linear interpolation of off-path excess enthalpy, not physical branch uniqueness.",
            "Two half-step anchors check even transverse curvature, not all q-dependent anharmonicity.",
            "Any figure must show actual DFT sample locations and the central-segment/domain boundary.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "trajectory", "atomic-feasibility", "joint-feasibility",
                 "transport", "first-audit", "hessian-npz", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.work, args.trajectory, args.atomic_feasibility,
                   args.joint_feasibility, args.transport, args.first_audit,
                   args.hessian_npz)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "max_LOO_meV_per_GaN": report["maximum_leave_one_out_error_meV_per_GaN"],
        "max_halfstep_curvature_relative_difference": report[
            "maximum_halfstep_curvature_relative_difference"
        ],
        "contour_gate_pass": report["central_contour_gate_pass"],
    }))


if __name__ == "__main__":
    main()
