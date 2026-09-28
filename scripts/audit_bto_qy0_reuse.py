"""Reaudit archived Qy=0 branches for reuse on the restricted BTO sheet.

This script performs no DFT. It checks the four existing open-Qy runs against
their original logs and evaluates stationarity in the *smaller*, fixed-Qy=0
subspace. A stationary Qy=0 branch is not claimed to be the unrestricted
conditional minimum or a curvature-certified restricted minimum.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.bto_q1q2_reference import (  # noqa: E402
    bto_transverse_soft_plane, load_bto_q1q2_reference, sha256,
    validate_bto_gamma_source,
)
from examples.preflight_bto_transverse_soft_conditional import (  # noqa: E402
    remaining_soft_y_direction,
)
from scripts.audit_bto_q1q2_conditional_pilot import (  # noqa: E402
    FINAL_ENERGY, _gradient_from_raw, _raw_force_stress,
)
from vcneb.mode_surface import _orthogonal_directions  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite audit: {args.output}")
    sources = json.loads(args.sources.read_text(encoding="utf-8"))
    if (sources.get("kind") != "bto_qy0_restricted_sheet_existing_branch_source_paths_not_new_DFT"
            or len(sources.get("cases", [])) not in (4, 5)):
        raise ValueError("expected four grid nodes and at most one archived holdout")
    common = Path(sources["common_inputs"])
    preflight = json.loads((common / "bto_qy0_q000_q000_preflight_2026-09-28.json").read_text())
    if preflight["kind"] != "bto_symmetry_restricted_soft_qy_zero_preflight_no_dft":
        raise ValueError("common BTO mode contract is not Qy-restricted")
    names = {
        "report": "bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json",
        "reference": "cubic_CONTCAR",
        "force_constants": "bto_cubic_gamma_force_constants.npz",
        "phonopy_eigenpairs": "bto_cubic_phonopy_gamma_eigenpairs.npz",
        "gamma_provenance": "bto_cubic_gamma_phonon_provenance.json",
        "force_sets": "bto_cubic_gamma_FORCE_SETS",
        "eigenpairs_provenance": "bto_cubic_phonopy_gamma_eigenpairs_provenance.json",
        "grid_result_manifest": "result_manifest_hf_27777454.json",
    }
    paths = {key: common / name for key, name in names.items()}
    if any(sha256(path) != preflight["source_sha256"][key]
           for key, path in paths.items()):
        raise ValueError("common BTO mode/DFT source files changed")
    validate_bto_gamma_source(
        paths["gamma_provenance"], paths["force_sets"],
        paths["eigenpairs_provenance"],
    )
    loaded = load_bto_q1q2_reference(
        paths["report"], paths["reference"], paths["force_constants"],
        paths["phonopy_eigenpairs"],
        strain_metric_weights_amu_A2=np.asarray(preflight["strain_metric_weights_amu_A2"]),
    )
    chart = loaded.chart
    plane = bto_transverse_soft_plane(
        loaded, strain_metric_weights_amu_A2=np.asarray(preflight["strain_metric_weights_amu_A2"]),
    )
    third = np.asarray(preflight["remaining_soft_y_metric_unit_direction"], dtype=float)
    if not np.allclose(third, remaining_soft_y_direction(loaded.modes, plane), atol=1e-10, rtol=0):
        raise ValueError("the third soft direction changed")
    open_directions = _orthogonal_directions(
        plane, np.column_stack([chart.rigid_translation_directions(), third]),
    )
    if open_directions.shape[1] != 15:
        raise ValueError("restricted open subspace must have 15 dimensions")
    grid = json.loads(paths["grid_result_manifest"].read_text())
    cubic = [row for row in grid["points"] if abs(row["q1"]) < 1e-12
             and abs(row["q2"]) < 1e-12]
    if grid["status"] != "converged" or grid["mpi_ranks"] != 32 or len(cubic) != 1:
        raise ValueError("the common cubic static reference is invalid")
    parameters = grid["calculator_parameters"]
    if (parameters.get("calculation") != "scf"
            or parameters.get("basis_type") != "lcao"
            or parameters.get("dft_functional") != "pbe"
            or float(parameters.get("ecutwfc", -1)) != 100.0
            or list(parameters.get("kpts", [])) != [4, 4, 4]
            or parameters.get("cal_force") != 1
            or parameters.get("cal_stress") != 1
            or any("10au" not in name for name in parameters.get("basis", {}).values())):
        raise ValueError("the common ABACUS 100-Ry/10-au-DZP static contract changed")
    expected_inputs = grid["points"][0]["input_sha256"]
    output_rows = []
    for case in sources["cases"]:
        workdir = Path(case["workdir"])
        replay_path, audit_path = Path(case["replay"]), Path(case["raw_audit"])
        if (sha256(replay_path) != case["replay_sha256"]
                or sha256(audit_path) != case["raw_audit_sha256"]):
            raise ValueError(f"archived replay/raw audit changed for {case['q']}")
        replay = json.loads(replay_path.read_text())
        audit = json.loads(audit_path.read_text())
        q = np.asarray(case["q"], dtype=float)
        branches = [row for row in replay.get("branch_outcomes", replay.get("branches", []))
                    if row["label"] == "frozen"]
        if (replay.get("status") not in {
                    "all_three_branch_outcomes_reproduced_from_raw_audited_cache",
                    "all_three_terminal_mappings_reconstructed_without_DFT",
                }
                or audit.get("status") != "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"
                or len(branches) != 1
                or not np.allclose(audit["q_parallel_q_transverse_sqrt_amu_A"], q, atol=1e-12, rtol=0)
                or not branches[0]["converged"]):
            raise ValueError(f"archived branch/audit status or Q differs for {case['q']}")
        branch = branches[0]
        terminal = workdir / branch["final_evaluation_directory"] / "result.json"
        if sha256(terminal) != case["terminal_result_sha256"]:
            raise ValueError(f"archived terminal result changed for {case['q']}")
        result = json.loads(terminal.read_text())
        raw_rows = [row for row in audit["evaluations"]
                    if row["directory"] == branch["final_evaluation_directory"]]
        if len(raw_rows) != 1:
            raise ValueError(f"terminal branch absent from raw audit for {case['q']}")
        raw = raw_rows[0]
        directory = terminal.parent
        log_path = directory / "OUT.ABACUS" / "running_scf.log"
        log = log_path.read_text(encoding="utf-8", errors="replace")
        if (result["status"] != "complete"
                or result["validation"]["mpi_dsize"] != 32
                or sha256(log_path) != raw["raw_log_sha256"]
                or raw["raw_log_sha256"] != result["validation"]["log_sha256"]
                or any(sha256(directory / name) != raw["raw_input_sha256"][name]
                       or raw["raw_input_sha256"][name] != result["validation"]["input_sha256"][name]
                       for name in ("INPUT", "KPT", "STRU"))
                or any(raw["raw_input_sha256"][name] != expected_inputs[name]
                       for name in ("INPUT", "KPT"))):
            raise ValueError(f"raw log, MPI or common DFT input changed for {case['q']}")
        coordinates = np.asarray(result["coordinates"], dtype=float)
        gradient = np.asarray(result["gradient"], dtype=float)
        atoms = chart.to_atoms(coordinates)
        if (coordinates.shape != (chart.coordinate_count,)
                or gradient.shape != coordinates.shape
                or not np.allclose(plane.project(coordinates), q, atol=1e-8, rtol=0)
                or abs(float(third @ (plane.metric_weights * coordinates))) > 1e-8):
            raise ValueError(f"terminal coordinates leave the Qy=0 sheet for {case['q']}")
        energies = FINAL_ENERGY.findall(log)
        forces, stress_kbar = _raw_force_stress(log, symbols=atoms.get_chemical_symbols())
        stress = -stress_kbar / 1602.176634
        relative = float(result["enthalpy_eV"] - cubic[0]["energy_eV"])
        projected = float(np.linalg.norm(open_directions.T @ gradient))
        max_stress = float(np.max(np.abs(stress_kbar)))
        if (len(energies) != 1
                or abs(float(energies[0]) - result["enthalpy_eV"]) > 1e-6
                or not np.allclose(_gradient_from_raw(chart, coordinates, forces, stress),
                                   gradient, atol=1e-6, rtol=0)
                or abs(relative - branch["energy_minus_c_eV_per_BTO"]) > 1e-8
                or abs(relative - raw["energy_minus_c_eV_per_BTO"]) > 1e-8
                or abs(float(np.linalg.norm(_orthogonal_directions(
                    plane, chart.rigid_translation_directions()).T @ gradient))
                       - branch["orthogonal_gradient_norm_eV_per_sqrt_amu_A"]) > 1e-8
                or projected > 0.003 or max_stress > 2.0):
            raise ValueError(f"Qy=0 restricted stationarity or raw output failed for {case['q']}")
        output_rows.append({
            "role": case.get("role", "grid_node"),
            "q_parallel_q_transverse_sqrt_amu_A": q.tolist(),
            "energy_minus_c_eV_per_BTO": relative,
            "restricted_15D_gradient_norm_eV_per_sqrt_amu_A": projected,
            "open_16D_gradient_norm_eV_per_sqrt_amu_A": branch["orthogonal_gradient_norm_eV_per_sqrt_amu_A"],
            "maximum_absolute_stress_kbar": max_stress,
            "third_soft_y_amplitude_sqrt_amu_A": float(third @ (plane.metric_weights * coordinates)),
            "coordinates_u_A_eta_voigt": coordinates.tolist(),
            "raw_log_sha256": raw["raw_log_sha256"],
            "terminal_result_sha256": case["terminal_result_sha256"],
            "replay_sha256": case["replay_sha256"],
            "raw_audit_sha256": case["raw_audit_sha256"],
        })
    output = {
        "kind": "bto_qy0_restricted_sheet_existing_raw_audited_stationary_branches_not_PES",
        "status": "existing_Qy0_points_eligible_for_restricted_sheet_reuse",
        "n_reused_DFT_points": len(output_rows),
        "n_new_DFT_evaluations": 0,
        "source_manifest_sha256": sha256(args.sources),
        "auditor_sha256": sha256(Path(__file__)),
        "reference_id": plane.reference_id,
        "electronic_kpoints": [4, 4, 4],
        "ecutwfc_Ry": 100,
        "phonon_supercell": [1, 1, 1],
        "fixed_third_soft_mode": "Q_y=0",
        "points": output_rows,
        "limitations": "stationarity/stress and raw-output reuse only; curvature, branch continuity and 2D interpolation remain open",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": output["status"], "n_reused_DFT_points": len(output_rows),
                      "points": [{"q": row["q_parallel_q_transverse_sqrt_amu_A"],
                                  "energy": row["energy_minus_c_eV_per_BTO"],
                                  "restricted_gradient": row["restricted_15D_gradient_norm_eV_per_sqrt_amu_A"]}
                                 for row in output_rows]}, indent=2))


if __name__ == "__main__":
    main()
