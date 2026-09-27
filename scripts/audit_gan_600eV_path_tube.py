"""Audit the 600-eV GaN path-adapted frozen transverse enthalpy tube."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.analyze_gan_600eV_path_tube_feasibility import ANCHORS, transported_normals
from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def leave_one_anchor_out_errors(anchors: list[int], arc: list[float],
                                by_anchor: dict[int, dict]) -> list[dict]:
    """Interpolate off-path *excess* enthalpy, not the curved 1D barrier."""

    errors = []
    for k in range(1, len(anchors) - 1):
        left, mid, right = anchors[k - 1:k + 2]
        weight = (arc[mid] - arc[left]) / (arc[right] - arc[left])
        for side in ("minus", "plus"):
            expected = (1 - weight) * by_anchor[left][side]["delta_enthalpy_meV_per_GaN"] \
                + weight * by_anchor[right][side]["delta_enthalpy_meV_per_GaN"]
            actual = by_anchor[mid][side]["delta_enthalpy_meV_per_GaN"]
            errors.append({
                "image_index": mid, "side": side,
                "interpolated_excess_meV_per_GaN": float(expected),
                "actual_excess_meV_per_GaN": float(actual),
                "absolute_error_meV_per_GaN": float(abs(actual - expected)),
            })
    return errors


def audit(work: Path, trajectory: Path, feasibility_path: Path,
          hessian_root: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    hessian_npz = hessian_root / "joint_hessian.npz"
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_45p7_600eV_path_adapted_transverse_enthalpy_tube_pilot"
            or manifest.get("anchors") != list(ANCHORS)
            or manifest.get("q_perp_A") != [-0.015, 0.0, 0.015]
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or len(manifest.get("cases", [])) != 18
            or manifest.get("source_sha256", {}).get("trajectory") != sha256(trajectory)
            or manifest.get("source_sha256", {}).get("feasibility") != sha256(feasibility_path)
            or manifest.get("source_sha256", {}).get("hessian_npz") != sha256(hessian_npz)
            or manifest.get("source_sha256", {}).get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_path_tube.py"))):
        raise ValueError("path-tube geometry or production input provenance failed")
    frames = read(trajectory, index=":")
    if len(frames) != 29 or sha256(trajectory) != feasibility["source_sha256"]["trajectory"]:
        raise ValueError("audited q=0 path changed")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        modes = np.asarray(archive["eigenvectors"], dtype=float)
    scale = 3.3982714330050063
    normals, transport = transported_normals(frames, modes[:, 1], scale)
    if any(not np.isclose(transport[key], feasibility["transport"][key], atol=1e-10, rtol=0)
           for key in transport):
        raise ValueError("transported normal changed from feasibility analysis")
    pressure = 45.7 * GPa
    arc = feasibility["arc_fraction_s"]
    h0 = feasibility["path_enthalpy_eV_per_cell"]
    if any(not np.isclose(float(frame.get_potential_energy() + pressure * frame.get_volume()),
                          h0[i], rtol=0, atol=1e-8)
           for i, frame in enumerate(frames)):
        raise ValueError("q=0 path enthalpy is not the audited chain")
    by_anchor: dict[int, dict] = {i: {} for i in ANCHORS}
    cases = []
    for record in manifest["cases"]:
        if any(record.get("input_sha256", {}).get(name) != digest
               for name, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"case {record['name']} changed 600-eV input contract")
        atoms, detail = completed_case(
            work / "cases" / record["name"],
            {**record, "POSCAR_sha256": record["input_sha256"]["POSCAR"]},
        )
        index = int(record["image_index"])
        q = float(record["q_perp_A"])
        side = "plus" if q > 0 else "minus"
        h = float(atoms.get_potential_energy() + pressure * atoms.get_volume())
        detail.update({
            "image_index": index, "arc_fraction_s": float(arc[index]),
            "q_perp_A": q,
            "enthalpy_eV_per_cell": h,
            "delta_enthalpy_meV_per_GaN": float((h - h0[index]) * 500),
        })
        if side in by_anchor[index]:
            raise ValueError(f"duplicate off-path side at anchor {index}")
        by_anchor[index][side] = detail
        cases.append(detail)
    diagnostics = []
    amplitude = 0.015
    for index in ANCHORS:
        if set(by_anchor[index]) != {"minus", "plus"}:
            raise ValueError(f"missing signed transverse pair at image {index}")
        low = by_anchor[index]["minus"]["enthalpy_eV_per_cell"]
        high = by_anchor[index]["plus"]["enthalpy_eV_per_cell"]
        chart = JointCurvatureCoordinates(frames[index], scale)
        gradient = chart.enthalpy_gradient(
            frames[index], frames[index].get_forces(),
            frames[index].get_stress(voigt=False), pressure,
        )
        analytic_slope = float(np.dot(gradient, normals[index]))
        finite_slope = float((high - low) / (2 * amplitude))
        curvature = float((high + low - 2 * h0[index]) / amplitude**2)
        diagnostics.append({
            "image_index": index,
            "arc_fraction_s": float(arc[index]),
            "analytic_transverse_slope_eV_per_A": analytic_slope,
            "energy_finite_difference_slope_eV_per_A": finite_slope,
            "slope_mismatch_eV_per_A": finite_slope - analytic_slope,
            "frozen_transverse_curvature_eV_per_A2": curvature,
        })
    holdouts = leave_one_anchor_out_errors(list(ANCHORS), arc, by_anchor)
    return {
        "status": "GaN_600eV_path_adapted_frozen_transverse_tube_raw_audited",
        "claim_limit": "Local path tube, not a global two-phonon PES or orthogonally relaxed conditional minimum",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "n_path_images": 29,
        "n_offpath_static_cases": 18,
        "anchors": list(ANCHORS),
        "arc_fraction_s": arc,
        "path_enthalpy_eV_per_cell": h0,
        "transverse_amplitude_A": amplitude,
        "anchor_diagnostics": diagnostics,
        "leave_one_anchor_out_excess_enthalpy_errors": holdouts,
        "maximum_leave_one_out_error_meV_per_GaN": float(
            max(item["absolute_error_meV_per_GaN"] for item in holdouts)
        ),
        "maximum_transverse_slope_mismatch_eV_per_A": float(
            max(abs(item["slope_mismatch_eV_per_A"]) for item in diagnostics)
        ),
        "cases": cases,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "trajectory": sha256(trajectory),
            "feasibility": sha256(feasibility_path),
            "hessian_npz": sha256(hessian_npz),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "The off-path grid has only nine path anchors and two transverse levels.",
            "Leave-one-anchor-out interpolation errors should guide adaptive densification.",
            "Finite-difference energy slope and analytic stress/force slope need not agree below their numerical error budget.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "trajectory", "feasibility", "hessian-root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work, args.trajectory, args.feasibility, args.hessian_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "max_loo_meV_per_GaN": result["maximum_leave_one_out_error_meV_per_GaN"],
        "max_slope_mismatch": result["maximum_transverse_slope_mismatch_eV_per_A"],
    }))


if __name__ == "__main__":
    main()
