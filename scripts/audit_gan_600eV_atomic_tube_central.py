"""Audit the bounded 600-eV GaN central atomic-transverse tube."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.analyze_gan_600eV_atomic_tube_feasibility import atomic_normals
from scripts.audit_gan_600eV_path_tube import leave_one_anchor_out_errors
from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_600eV_atomic_tube_central import AMPLITUDE_A, ANCHORS
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(work: Path, trajectory: Path, atomic_feasibility_path: Path,
          joint_feasibility_path: Path, transport_path: Path,
          hessian_npz: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    atomic_feasibility = json.loads(atomic_feasibility_path.read_text(encoding="utf-8"))
    joint_feasibility = json.loads(joint_feasibility_path.read_text(encoding="utf-8"))
    transport = json.loads(transport_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_45p7_600eV_atomic_only_central_path_tube_pilot"
            or manifest.get("anchors") != list(ANCHORS)
            or manifest.get("q_atom_A") != [-AMPLITUDE_A, 0.0, AMPLITUDE_A]
            or manifest.get("predeclared_LOO_gate_meV_per_GaN") != 1.0
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or len(manifest.get("cases", [])) != 2 * len(ANCHORS)
            or manifest["source_sha256"].get("trajectory") != sha256(trajectory)
            or manifest["source_sha256"].get("atomic_feasibility")
            != sha256(atomic_feasibility_path)
            or manifest["source_sha256"].get("joint_feasibility")
            != sha256(joint_feasibility_path)
            or manifest["source_sha256"].get("transport_comparison")
            != sha256(transport_path)
            or manifest["source_sha256"].get("hessian_npz") != sha256(hessian_npz)
            or manifest["source_sha256"].get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_atomic_tube_central.py"))
            or atomic_feasibility.get("status")
            != "GaN_600eV_atomic_transverse_tube_geometry_only_no_DFT"
            or joint_feasibility.get("status")
            != "GaN_600eV_path_adapted_joint_normal_tube_geometry_only"
            or transport.get("status")
            != "GaN_600eV_atomic_transverse_transport_offline_comparison"
            or transport["fixed_seed_projection"]["minimum_central_5_to_22_overlap"] <= 0.95):
        raise ValueError("central path, chart, or 600-eV electronic provenance failed")
    frames = read(trajectory, index=":")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        seed = np.asarray(archive["eigenvectors"][:, 1], dtype=float)
    normals, metrics = atomic_normals(frames, seed)
    if len(frames) != 29 or any(
            not np.isclose(value, atomic_feasibility["transport"][key], atol=1e-10, rtol=0)
            for key, value in metrics.items()):
        raise ValueError("atomic normal changed since input generation")
    pressure = 45.7 * GPa
    arc = joint_feasibility["arc_fraction_s"]
    h0 = joint_feasibility["path_enthalpy_eV_per_cell"]
    if any(not np.isclose(float(frame.get_potential_energy() + pressure * frame.get_volume()),
                          h0[index], atol=1e-8, rtol=0)
           for index, frame in enumerate(frames)):
        raise ValueError("archived q=0 enthalpy changed")
    by_anchor: dict[int, dict[str, dict]] = {index: {} for index in ANCHORS}
    cases = []
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
            "image_index": index,
            "arc_fraction_s": float(arc[index]),
            "q_atom_A": q,
            "enthalpy_eV_per_cell": h,
            "delta_enthalpy_meV_per_GaN": float((h - h0[index]) * 500),
        })
        if side in by_anchor[index]:
            raise ValueError("duplicate signed atomic point")
        by_anchor[index][side] = detail
        cases.append(detail)
    diagnostics = []
    for index in ANCHORS:
        if set(by_anchor[index]) != {"minus", "plus"}:
            raise ValueError("atomic central tube signed pair incomplete")
        low = by_anchor[index]["minus"]["enthalpy_eV_per_cell"]
        high = by_anchor[index]["plus"]["enthalpy_eV_per_cell"]
        analytic_slope = -float(np.dot(frames[index].get_forces().ravel(), normals[index]))
        finite_slope = float((high - low) / (2 * AMPLITUDE_A))
        diagnostics.append({
            "image_index": index,
            "arc_fraction_s": float(arc[index]),
            "analytic_atomic_transverse_slope_eV_per_A": analytic_slope,
            "energy_finite_difference_slope_eV_per_A": finite_slope,
            "slope_mismatch_eV_per_A": finite_slope - analytic_slope,
            "frozen_atomic_transverse_curvature_eV_per_A2":
                float((high + low - 2 * h0[index]) / AMPLITUDE_A**2),
        })
    holdouts = leave_one_anchor_out_errors(list(ANCHORS), arc, by_anchor)
    maximum = max(item["absolute_error_meV_per_GaN"] for item in holdouts)
    return {
        "status": "GaN_600eV_central_atomic_transverse_tube_raw_audited",
        "claim_limit": "Central segment only, frozen atomic-only q; not a global two-phonon PES",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "n_path_images": 29,
        "n_offpath_static_cases": len(cases),
        "anchors": list(ANCHORS),
        "arc_fraction_s": arc,
        "path_enthalpy_eV_per_cell": h0,
        "atomic_transverse_amplitude_A": AMPLITUDE_A,
        "predeclared_LOO_gate_meV_per_GaN": 1.0,
        "leave_one_anchor_out_excess_enthalpy_errors": holdouts,
        "maximum_leave_one_out_error_meV_per_GaN": float(maximum),
        "central_contour_gate_pass": bool(maximum <= 1.0),
        "maximum_energy_force_slope_mismatch_eV_per_A": float(
            max(abs(item["slope_mismatch_eV_per_A"]) for item in diagnostics)
        ),
        "anchor_diagnostics": diagnostics,
        "cases": cases,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "trajectory": sha256(trajectory),
            "atomic_feasibility": sha256(atomic_feasibility_path),
            "joint_feasibility": sha256(joint_feasibility_path),
            "transport_comparison": sha256(transport_path),
            "hessian_npz": sha256(hessian_npz),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "Leave-one-anchor-out error measures linear interpolation of signed off-path excess enthalpy, not DFT convergence or global branch uniqueness.",
            "The atomic normal has poor identity continuity outside the central segment.",
            "The transverse q coordinate holds each variable-cell path image cell fixed; orthogonal degrees of freedom are not relaxed.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "trajectory", "atomic-feasibility", "joint-feasibility",
                 "transport", "hessian-npz", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.work, args.trajectory, args.atomic_feasibility,
                   args.joint_feasibility, args.transport, args.hessian_npz)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "max_LOO_meV_per_GaN": report["maximum_leave_one_out_error_meV_per_GaN"],
        "contour_gate_pass": report["central_contour_gate_pass"],
    }))


if __name__ == "__main__":
    main()
