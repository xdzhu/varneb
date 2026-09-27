"""Audit the targeted 600-eV GaN path-tube refinement and combined LOO error."""

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
from scripts.prepare_gan_600eV_path_tube_refinement import REFINEMENT_ANCHORS
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(work: Path, trajectory: Path, feasibility_path: Path,
          first_audit_path: Path, transport_path: Path, screen_path: Path,
          hessian_npz: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    first = json.loads(first_audit_path.read_text(encoding="utf-8"))
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_45p7_600eV_path_tube_error_targeted_refinement"
            or manifest.get("anchors") != list(REFINEMENT_ANCHORS)
            or manifest.get("q_perp_A") != [-0.015, 0.015]
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or len(manifest.get("cases", [])) != 18
            or manifest["source_sha256"].get("trajectory") != sha256(trajectory)
            or manifest["source_sha256"].get("feasibility") != sha256(feasibility_path)
            or manifest["source_sha256"].get("first_tube_audit") != sha256(first_audit_path)
            or manifest["source_sha256"].get("transport_comparison") != sha256(transport_path)
            or manifest["source_sha256"].get("geometry_screen") != sha256(screen_path)
            or manifest["source_sha256"].get("hessian_npz") != sha256(hessian_npz)
            or manifest["source_sha256"].get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_path_tube_refinement.py"))
            or first.get("status")
            != "GaN_600eV_path_adapted_frozen_transverse_tube_raw_audited"
            or first["source_sha256"].get("trajectory") != sha256(trajectory)
            or first["source_sha256"].get("feasibility") != sha256(feasibility_path)
            or first["source_sha256"].get("hessian_npz") != sha256(hessian_npz)
            or screen.get("status") != "GaN_600eV_offpath_refinement_geometry_only_no_DFT"):
        raise ValueError("same-contract refinement provenance failed")
    frames = read(trajectory, index=":")
    if (len(frames) != 29
            or not np.allclose(first["arc_fraction_s"], feasibility["arc_fraction_s"], atol=1e-12)
            or not np.allclose(first["path_enthalpy_eV_per_cell"],
                               feasibility["path_enthalpy_eV_per_cell"], atol=1e-8)):
        raise ValueError("audited centerline changed")
    pressure = 45.7 * GPa
    arc = first["arc_fraction_s"]
    h0 = first["path_enthalpy_eV_per_cell"]
    if any(not np.isclose(float(frame.get_potential_energy() + pressure * frame.get_volume()),
                          h0[i], atol=1e-8, rtol=0)
           for i, frame in enumerate(frames)):
        raise ValueError("q=0 enthalpy is not the archived 600-eV path")
    by_anchor: dict[int, dict] = {index: {} for index in first["anchors"]}
    for record in first["cases"]:
        index = int(record["image_index"])
        side = "plus" if record["q_perp_A"] > 0 else "minus"
        if side in by_anchor[index]:
            raise ValueError("duplicate first-batch anchor")
        by_anchor[index][side] = record
    new_cases = []
    for record in manifest["cases"]:
        if any(record.get("input_sha256", {}).get(name) != digest
               for name, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"{record['name']} changed electronic parameters")
        atoms, detail = completed_case(
            work / "cases" / record["name"],
            {**record, "POSCAR_sha256": record["input_sha256"]["POSCAR"]},
        )
        index = int(record["image_index"])
        q = float(record["q_perp_A"])
        side = "plus" if q > 0 else "minus"
        if index in by_anchor and side in by_anchor[index]:
            raise ValueError("refinement repeats an audited point")
        h = float(atoms.get_potential_energy() + pressure * atoms.get_volume())
        detail.update({
            "image_index": index,
            "arc_fraction_s": float(arc[index]),
            "q_perp_A": q,
            "enthalpy_eV_per_cell": h,
            "delta_enthalpy_meV_per_GaN": float((h - h0[index]) * 500),
        })
        by_anchor.setdefault(index, {})[side] = detail
        new_cases.append(detail)
    anchors = sorted(by_anchor)
    if any(set(by_anchor[index]) != {"minus", "plus"} for index in anchors):
        raise ValueError("combined grid lacks a signed pair")
    holdouts = leave_one_anchor_out_errors(anchors, arc, by_anchor)
    maximum = max(item["absolute_error_meV_per_GaN"] for item in holdouts)
    return {
        "status": "GaN_600eV_path_tube_refined_raw_audited",
        "claim_limit": "Measured frozen tube; smooth contour requires holdout gate below 0.10 meV/GaN",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "n_path_images": 29,
        "n_original_offpath_cases": len(first["cases"]),
        "n_new_offpath_cases": len(new_cases),
        "anchors": anchors,
        "arc_fraction_s": arc,
        "path_enthalpy_eV_per_cell": h0,
        "transverse_amplitude_A": 0.015,
        "leave_one_anchor_out_excess_enthalpy_errors": holdouts,
        "maximum_leave_one_out_error_meV_per_GaN": float(maximum),
        "predeclared_interpolation_gate_meV_per_GaN": 0.10,
        "smooth_contour_gate_pass": bool(maximum <= 0.10),
        "new_cases": new_cases,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "trajectory": sha256(trajectory),
            "feasibility": sha256(feasibility_path),
            "first_audit": sha256(first_audit_path),
            "transport_comparison": sha256(transport_path),
            "geometry_screen": sha256(screen_path),
            "hessian_npz": sha256(hessian_npz),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "Leave-one-anchor-out linear interpolation is a deterministic interpolation check, not a DFT error estimate.",
            "The transverse normal is path-dependent and is not one global Gamma phonon.",
            "The tube is frozen, not orthogonally relaxed or a global two-phonon surface.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "trajectory", "feasibility", "first-audit", "transport",
                 "screen", "hessian-npz", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work, args.trajectory, args.feasibility,
                   args.first_audit, args.transport, args.screen, args.hessian_npz)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "n_new_cases": result["n_new_offpath_cases"],
        "max_loo_meV_per_GaN": result["maximum_leave_one_out_error_meV_per_GaN"],
        "smooth_contour_gate_pass": result["smooth_contour_gate_pass"],
    }))


if __name__ == "__main__":
    main()
