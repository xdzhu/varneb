"""Raw-log audit and prospective 59-to-81 point BTO frozen-cut check."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.units import GPa
from scipy.interpolate import CloughTocher2DInterpolator


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _number_rows(lines: list[str], marker: str, count: int, offset: int) -> np.ndarray:
    indices = [index for index, line in enumerate(lines) if marker in line]
    if len(indices) != 1:
        raise ValueError(f"expected one {marker} table")
    values = []
    for line in lines[indices[0] + offset:indices[0] + offset + count]:
        fields = line.split()
        try:
            values.append([float(value) for value in fields[-3:]])
        except ValueError as exc:
            raise ValueError(f"invalid {marker} numeric row") from exc
    result = np.asarray(values, dtype=float)
    if result.shape != (count, 3) or not np.isfinite(result).all():
        raise ValueError(f"incomplete {marker} table")
    return result


def audit(dense_path: Path, old_path: Path, points_dir: Path,
          result_path: Path, old_result_path: Path, work: Path) -> dict:
    dense, old, staged, result, old_result = map(_load, (
        dense_path, old_path, points_dir / "manifest.json", result_path, old_result_path,
    ))
    if (dense.get("n_points") != 81 or len(dense.get("points", [])) != 81
            or old.get("n_real_DFT_points") != 59 or len(old.get("samples", [])) != 59
            or old.get("source_sha256", {}).get("dense_manifest") != sha256(dense_path)
            or staged.get("kind") != "inert_frozen_cubic_cell_transverse_soft_missing22_not_T_to_C_barrier"
            or staged.get("n_points") != 22 or len(staged.get("points", [])) != 22
            or staged.get("source_dense_manifest_sha256") != sha256(dense_path)
            or staged.get("source_59_analysis_sha256") != sha256(old_path)
            or result.get("status") != "converged"
            or result.get("source_manifest_sha256") != sha256(points_dir / "manifest.json")
            or result.get("source_kind") != staged["kind"]
            or result.get("dft_requested") is not True
            or len(result.get("points", [])) != 22
            or result.get("mpi_ranks") != 32
            or old_result.get("status") != "converged"
            or result.get("calculator_parameters") != old_result.get("calculator_parameters")
            or result.get("asset_sha256") != old_result.get("asset_sha256")
            or result.get("abacus_binary_sha256") != old_result.get("abacus_binary_sha256")):
        raise ValueError("BTO 59-to-81 result/protocol provenance is incomplete")
    parameters = result["calculator_parameters"]
    if (parameters.get("calculation") != "scf" or parameters.get("basis_type") != "lcao"
            or parameters.get("dft_functional") != "pbe"
            or float(parameters.get("ecutwfc", -1)) != 100.0
            or list(parameters.get("kpts", [])) != [4, 4, 4]
            or parameters.get("cal_force") != 1 or parameters.get("cal_stress") != 1):
        raise ValueError("BTO production single-point contract changed")
    old_input = {point["input_sha256"]["INPUT"] for point in old_result["points"]}
    old_kpt = {point["input_sha256"]["KPT"] for point in old_result["points"]}
    if len(old_input) != 1 or len(old_kpt) != 1:
        raise ValueError("old BTO INPUT/KPT hashes are not uniform")
    dense_by_index = {tuple(point["grid_index_q1_q2"]): point for point in dense["points"]}
    old_by_index = {(point["dense_i_q1"], point["dense_j_q2"]): point
                    for point in old["samples"]}
    result_by_name = {point["name"]: point for point in result["points"]}
    if len(dense_by_index) != 81 or len(old_by_index) != 59 or len(result_by_name) != 22:
        raise ValueError("duplicate BTO dense, old, or new coordinate")
    for index, item in old_by_index.items():
        if dense_by_index[index]["structure_sha256"] != item["source_structure_sha256"]:
            raise ValueError(f"old frozen structure changed at {index}")
    old_xy = np.asarray([[item["q_parallel_sqrt_amu_A"],
                          item["q_transverse_sqrt_amu_A"]] for item in old["samples"]])
    old_energy = np.asarray([item["energy_minus_C_eV_per_BTO"] for item in old["samples"]])
    predictor = CloughTocher2DInterpolator(old_xy, old_energy)
    new_samples = []
    for point in staged["points"]:
        index = tuple(point["grid_index_q1_q2"])
        name = point["name"]
        recorded = result_by_name.get(name)
        directory = work / name
        if (index in old_by_index or dense_by_index[index]["structure_sha256"] != point["structure_sha256"]
                or recorded is None or recorded.get("status") != "converged"
                or recorded.get("source_structure_sha256") != point["structure_sha256"]
                or sha256(directory / "source_POSCAR") != point["structure_sha256"]
                or recorded["input_sha256"]["INPUT"] not in old_input
                or recorded["input_sha256"]["KPT"] not in old_kpt):
            raise ValueError(f"new BTO static record failed source audit: {name}")
        if any(sha256(directory / key) != digest
               for key, digest in recorded["input_sha256"].items()):
            raise ValueError(f"BTO INPUT/KPT/STRU changed after run: {name}")
        log_path = directory / "OUT.ABACUS" / "running_scf.log"
        if sha256(log_path) != recorded["log_sha256"]:
            raise ValueError(f"BTO log hash differs from result: {name}")
        log = log_path.read_text(encoding="utf-8", errors="replace")
        if log.count("charge density convergence is achieved") != 1:
            raise ValueError(f"BTO static SCF was not uniquely converged: {name}")
        mpi = re.findall(r"\bDSIZE\s*=\s*(\d+)", log)
        final = re.findall(r"!FINAL_ETOT_IS\s+([-+\d.eE]+)\s+eV", log)
        if mpi != ["32"] or len(final) != 1:
            raise ValueError(f"BTO MPI or final-energy marker invalid: {name}")
        lines = log.splitlines()
        forces = _number_rows(lines, "TOTAL-FORCE (eV/Angstrom)", 5, 2)
        raw_stress_kbar = _number_rows(lines, "TOTAL-STRESS (KBAR)", 3, 2)
        energy = float(final[0])
        stress = -raw_stress_kbar * 0.1 * GPa
        if (abs(energy - float(recorded["energy_eV"])) > 1e-6
                or abs(np.max(np.linalg.norm(forces, axis=1))
                       - float(recorded["max_atomic_force_eV_per_A"])) > 1e-6
                or not np.allclose(stress, recorded["stress_eV_per_A3"], atol=1e-7, rtol=0)):
            raise ValueError(f"BTO raw energy/force/stress differs from ASE result: {name}")
        x, y = float(point["q1"]), float(point["q2"])
        if (x != float(dense["grid"]["q1"][index[0]])
                or y != float(dense["grid"]["q2"][index[1]])):
            raise ValueError(f"BTO grid coordinate differs from manifest: {name}")
        measured = energy - float(old["reference_energy_eV_per_BTO"])
        prediction = float(predictor(x, y))
        if not np.isfinite([measured, prediction]).all():
            raise ValueError(f"BTO predictor or DFT energy nonfinite: {name}")
        new_samples.append({
            "name": name, "dense_i_q1": index[0], "dense_j_q2": index[1],
            "q_parallel_sqrt_amu_A": x, "q_transverse_sqrt_amu_A": y,
            "energy_eV_per_BTO": energy, "energy_minus_C_eV_per_BTO": measured,
            "prior_59_point_prediction_eV_per_BTO": prediction,
            "prospective_error_meV_per_BTO": float(1000 * (prediction - measured)),
            "source_structure_sha256": point["structure_sha256"],
            "log_sha256": recorded["log_sha256"],
        })
    combined = list(old["samples"]) + new_samples
    if len({(item["dense_i_q1"], item["dense_j_q2"]) for item in combined}) != 81:
        raise ValueError("BTO 81-point grid is incomplete")
    errors = np.asarray([item["prospective_error_meV_per_BTO"] for item in new_samples])
    return {
        "kind": "BTO_transverse_soft_frozen_C_cell_81_real_DFT_points_not_conditional_PES_or_MEP",
        "status": "22_new_nodes_raw_energy_force_stress_audited",
        "n_prior_DFT_points": 59, "n_new_DFT_points": 22, "n_total_DFT_points": 81,
        "slurm_job_id": result["slurm_job_id"],
        "reference_energy_eV_per_BTO": old["reference_energy_eV_per_BTO"],
        "boundary_condition": old["boundary_condition"],
        "axis_labels": old["axis_labels"], "axis_mode_weights": old["axis_mode_weights"],
        "prospective_22_point_max_abs_error_meV_per_BTO": float(np.max(np.abs(errors))),
        "prospective_22_point_rms_error_meV_per_BTO": float(np.sqrt(np.mean(errors**2))),
        "new_holdouts": new_samples, "samples": combined,
        "source_sha256": {
            "dense_manifest": sha256(dense_path), "prior_59_audit": sha256(old_path),
            "new_points_manifest": sha256(points_dir / "manifest.json"),
            "new_result_manifest": sha256(result_path),
            "old_result_manifest": sha256(old_result_path),
            "auditor": sha256(Path(__file__)),
        },
        "claim_limit": "Frozen C-cell atomic-mode slice only; contour, minimum and path shadow are not a relaxed T-to-C barrier.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("dense-manifest", "analysis59", "points-dir", "result", "old-result", "work", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.dense_manifest, args.analysis59, args.points_dir,
                   args.result, args.old_result, args.work)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "status", "n_total_DFT_points", "prospective_22_point_max_abs_error_meV_per_BTO",
    )}))


if __name__ == "__main__":
    main()
