"""Reuse eight G0 static probes; no DFT or transition-state certification."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.stress import voigt_6_to_full_3x3_stress

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.audit_hfo2_static_replica import INPUT_FILES, audited_results, sha256
from vcneb import assemble_joint_directional_curvature
from vcneb.joint_curvature import JointCurvatureCoordinates


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "benchmarks/hfo2_channels/20261008"


def analyze(manifest: dict, results: dict[int, dict]) -> dict:
    """Assemble selected columns without treating missing directions as stable."""
    basis = np.array(manifest["directions"], dtype=float).T
    center_gradient = np.array(manifest["center_gradient_eV_A"], dtype=float)
    if (basis.shape != (42, 2) or center_gradient.shape != (42,)
            or not np.isfinite(center_gradient).all() or manifest["pressure_eV_A3"] != 0
            or not np.isfinite(manifest["cell_scale_A"]) or manifest["cell_scale_A"] <= 0):
        raise ValueError("expected the registered zero-pressure 12-atom free-cell G0 chart")
    expected = {(direction, step, sign) for direction in (0, 1)
                for step in (.01, .02) for sign in (-1, 1)}
    tags = [(p["direction"], p["step_A"], p["sign"]) for p in manifest["points"]]
    indices = [p["index"] for p in manifest["points"]]
    if len(tags) != 8 or set(tags) != expected or set(indices) != set(range(8)) or set(results) != set(range(8)):
        raise ValueError("exactly eight uniquely tagged complete G0 results required")
    for index, point in results.items():
        if (point["index"] != index or point["status"] != "SCF_and_input_contract_passed"
                or not np.isfinite(point["energy_eV_cell"])
                or len(point["log_sha256"]) != 64):
            raise ValueError("incomplete or incorrectly tagged audited point")
    matrices, records = [], []
    for step in (.01, .02):
        pairs = {(p["direction"], p["sign"]): results[p["index"]]
                 for p in manifest["points"] if p["step_A"] == step}
        plus = np.array([pairs[(i, 1)]["gradient_eV_A"] for i in (0, 1)])
        minus = np.array([pairs[(i, -1)]["gradient_eV_A"] for i in (0, 1)])
        curvature = assemble_joint_directional_curvature(basis, plus, minus, step_A=step)
        matrices.append(curvature)
        records.append({
            "step_A": step, "raw_projected_eV_A2": curvature.raw_projected.tolist(),
            "symmetric_projected_eV_A2": curvature.symmetric_projected.tolist(),
            "projected_eigenvalues_eV_A2": np.linalg.eigvalsh(curvature.symmetric_projected).tolist(),
            "reciprocity_relative_defect": curvature.reciprocity_relative_defect,
            "symmetrization_operator_change_eV_A2": curvature.symmetrization_operator_change_eV_A2,
            "full_hessian_action_eV_A2": curvature.full_hessian_action.tolist(),
            "transverse_action_eV_A2": curvature.transverse_action.tolist(),
            "transverse_column_norms_eV_A2": np.linalg.norm(curvature.transverse_action, axis=0).tolist(),
            "cell_response_to_atomic_direction_norm_eV_A2": float(np.linalg.norm(curvature.full_hessian_action[36:, 0])),
            "atomic_response_to_cell_direction_norm_eV_A2": float(np.linalg.norm(curvature.full_hessian_action[:36, 1])),
        })
    return {
        "format_version": 1, "status": "selected_G0_Hessian_columns_assembled",
        "scope": "exploratory reuse of nonstationary free-cell G0 probes, not G2 clamped training or prospective prediction",
        "pressure_eV_A3": 0, "cell_scale_A": manifest["cell_scale_A"],
        "center_full_gradient_norm_eV_A": float(np.linalg.norm(center_gradient)),
        "direction_semantics": manifest["direction_semantics"], "directions": basis.tolist(),
        "step_records": records,
        "two_step_full_action_operator_spread_eV_A2": float(np.linalg.norm(matrices[0].full_hessian_action - matrices[1].full_hessian_action, ord=2)),
        "two_step_projected_operator_spread_eV_A2": float(np.linalg.norm(matrices[0].symmetric_projected - matrices[1].symmetric_projected, ord=2)),
        "source_log_sha256": manifest["source_log_sha256"],
        "probe_log_sha256": {str(i): results[i]["log_sha256"] for i in range(8)},
        "new_DFT_calls": 0, "parameters_changed": False, "raw_source_audit_performed": False,
        "full_joint_Hessian_measured": False, "unsampled_stability_known": False,
        "saddle_certified": False, "conditional_surface_established": False,
        "independent_barrier_prediction_validated": False,
        "interpretation": "Nonzero transverse response requires extra measurements before elimination; positive two-direction eigenvalues do not prove full stability. Two-step spread is not a rigorous total DFT error bar.",
    }


def audit_raw(probe_root: Path) -> tuple[dict, dict[int, dict]]:
    """Read historical sources only; check contract, geometry and old gradients."""
    manifest = json.loads((probe_root / "probe_manifest.json").read_text())
    archived = json.loads((EVIDENCE / "work_probe_manifest.json").read_text())
    if manifest != archived:
        raise ValueError("historical G0 manifest content changed")
    source = Path(manifest["source_directory"])
    if (sha256(source / "STRU") != manifest["source_STRU_sha256"]
            or sha256(source / "OUT.ABACUS/running_scf.log") != manifest["source_log_sha256"]):
        raise ValueError("historical center geometry/log changed")
    if any(sha256(source / name) != digest for name, digest in CONTRACT.items()):
        raise ValueError("historical center calculator contract changed")
    center = read(source / "STRU", format="abacus")
    chart = JointCurvatureCoordinates(center, manifest["cell_scale_A"])
    raw_center = audited_results(source)
    gradient0 = chart.enthalpy_gradient(center, raw_center["forces"], voigt_6_to_full_3x3_stress(raw_center["stress"]), 0)
    if (abs(raw_center["energy"] - manifest["center_energy_eV_cell"]) > 1e-10
            or not np.allclose(gradient0, manifest["center_gradient_eV_A"], atol=1e-12, rtol=0)):
        raise ValueError("historical center results differ")
    results = {}
    for item in manifest["points"]:
        index = item["index"]
        directory = probe_root / "calculations" / f"{index:02d}"
        if (any(item["input_sha256"].get(name) != digest for name, digest in CONTRACT.items())
                or any(sha256(directory / name) != item["input_sha256"][name] for name in INPUT_FILES)):
            raise ValueError("probe calculator/geometry bytes changed")
        expected = chart.displaced(item["sign"] * item["step_A"] * np.array(manifest["directions"][item["direction"]]))
        atoms = read(directory / "STRU", format="abacus")
        if not same_ordered_geometry(atoms, expected):
            raise ValueError("probe geometry does not match its signed chart coordinate")
        raw = audited_results(directory)
        gradient = chart.enthalpy_gradient(expected, raw["forces"], voigt_6_to_full_3x3_stress(raw["stress"]), 0)
        archived_point = json.loads((EVIDENCE / "point_audits/calculations" / f"{index:02d}" / "point_audit.json").read_text())
        log_hash = sha256(directory / "OUT.ABACUS/running_scf.log")
        if (log_hash != archived_point["log_sha256"]
                or abs(raw["energy"] - archived_point["energy_eV_cell"]) > 1e-10
                or not np.allclose(gradient, archived_point["gradient_eV_A"], atol=1e-12, rtol=0)):
            raise ValueError("raw probe results differ from archived evidence")
        results[index] = {"index": index, "status": "SCF_and_input_contract_passed",
                          "energy_eV_cell": float(raw["energy"]), "gradient_eV_A": gradient.tolist(), "log_sha256": log_hash}
    return manifest, results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing an existing analysis receipt")
    if args.raw_root is None:
        manifest = json.loads((EVIDENCE / "work_probe_manifest.json").read_text())
        results = {i: json.loads((EVIDENCE / "point_audits/calculations" / f"{i:02d}" / "point_audit.json").read_text()) for i in range(8)}
    else:
        manifest, results = audit_raw(args.raw_root)
    report = analyze(manifest, results)
    report.update({"raw_source_audit_performed": args.raw_root is not None,
                   "raw_root": str(args.raw_root) if args.raw_root is not None else None,
                   "raw_probe_manifest_sha256": sha256(args.raw_root / "probe_manifest.json") if args.raw_root is not None else None,
                   "archived_manifest_sha256": sha256(EVIDENCE / "work_probe_manifest.json"),
                   "analysis_source_sha256": sha256(Path(__file__)),
                   "joint_stencil_source_sha256": sha256(ROOT / "vcneb/joint_stencil.py")})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "new_DFT_calls": 0, "raw_audit": report["raw_source_audit_performed"]}))


if __name__ == "__main__":
    main()
