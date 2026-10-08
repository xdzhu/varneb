"""Audit complete fixed-contract probes and compare energy/gradient curvature."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.stress import voigt_6_to_full_3x3_stress

from scripts.audit_hfo2_static_replica import INPUT_FILES, audited_results, sha256
from vcneb.joint_curvature import JointCurvatureCoordinates


def audit_point(root: Path, index: int) -> dict:
    manifest = json.loads((root / "probe_manifest.json").read_text(encoding="utf-8"))
    point = manifest["points"][index]
    directory = root / "calculations" / f"{index:02d}"
    if any(sha256(directory / name) != point["input_sha256"][name] for name in INPUT_FILES):
        raise ValueError("planned input hash changed before/after SCF")
    raw = audited_results(directory)
    center = read(Path(manifest["source_directory"]) / "STRU", format="abacus")
    if sha256(Path(manifest["source_directory"]) / "STRU") != manifest["source_STRU_sha256"]:
        raise ValueError("source geometry changed")
    chart = JointCurvatureCoordinates(center, manifest["cell_scale_A"])
    directions = np.array(manifest["directions"])
    atoms = chart.displaced(point["sign"] * point["step_A"] * directions[point["direction"]])
    gradient = chart.enthalpy_gradient(atoms, raw["forces"], voigt_6_to_full_3x3_stress(raw["stress"]), 0)
    return {"status": "SCF_and_input_contract_passed", "index": index,
            "energy_eV_cell": float(raw["energy"]), "gradient_eV_A": gradient.tolist(),
            "directional_gradient_eV_A": float(gradient @ directions[point["direction"]]),
            "log_sha256": sha256(directory / "OUT.ABACUS/running_scf.log")}


def score_pairs(manifest: dict, results: dict[int, dict]) -> list[dict]:
    scored = []
    for direction in (0, 1):
        for step in (0.01, 0.02):
            pair = {point["sign"]: results.get(point["index"]) for point in manifest["points"]
                    if point["direction"] == direction and point["step_A"] == step}
            if not pair.get(-1) or not pair.get(1):
                continue
            minus, plus = pair[-1], pair[1]
            e_minus, e_plus = minus["energy_eV_cell"], plus["energy_eV_cell"]
            energy_gradient = (e_plus - e_minus) / (2 * step)
            gradient_curvature = (plus["directional_gradient_eV_A"] - minus["directional_gradient_eV_A"]) / (2 * step)
            energy_curvature = (e_plus + e_minus - 2 * manifest["center_energy_eV_cell"]) / step ** 2
            scored.append({"direction": direction, "step_A": step,
                           "central_energy_gradient_eV_A": energy_gradient,
                           "center_force_stress_gradient_eV_A": manifest["predicted_center_directional_gradients_eV_A"][direction],
                           "gradient_difference_eV_A": energy_gradient - manifest["predicted_center_directional_gradients_eV_A"][direction],
                           "energy_curvature_eV_A2": energy_curvature,
                           "gradient_curvature_eV_A2": gradient_curvature,
                           "curvature_difference_eV_A2": energy_curvature - gradient_curvature,
                           "interpretation": "local direction only; nonstationary center, not full Hessian or TS"})
    return scored


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    if args.index is not None:
        result = audit_point(args.root, args.index)
        output = args.root / "calculations" / f"{args.index:02d}" / "point_audit.json"
        output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"index": args.index, "status": result["status"]}))
        return
    manifest = json.loads((args.root / "probe_manifest.json").read_text(encoding="utf-8"))
    results = {}
    for point in manifest["points"]:
        path = args.root / "calculations" / f"{point['index']:02d}" / "point_audit.json"
        if path.exists():
            results[point["index"]] = audit_point(args.root, point["index"])
    report = {"status": "complete" if len(results) == 8 else "partial",
              "n_audited_points": len(results), "pairs": score_pairs(manifest, results),
              "TS_certified": False, "parameters_changed": False}
    (args.root / "work_probe_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
