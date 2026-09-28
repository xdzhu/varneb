"""Read-only scientific audit of one BTO conditional-Q1/Q2 ABACUS pilot.

Requires original calculator logs, input files, result cache, and the saved
mode/reference contract. It does not run DFT or infer a 2D PES from one point.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.bto_q1q2_reference import (
    bto_transverse_soft_plane, load_bto_q1q2_reference, sha256,
    validate_bto_gamma_source,
)
from examples.preflight_bto_transverse_soft_conditional import remaining_soft_y_direction
from vcneb.mode_surface import _orthogonal_directions


FINAL_ENERGY = re.compile(r"!FINAL_ETOT_IS\s*=?\s*([-+\d.eE]+)\s*eV")
RAW_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?"
FORCE_ROW = re.compile(rf"\s*([A-Z][a-z]?\d+)\s+({RAW_FLOAT})\s+({RAW_FLOAT})\s+({RAW_FLOAT})\s*")
STRESS_ROW = re.compile(rf"\s*({RAW_FLOAT})\s+({RAW_FLOAT})\s+({RAW_FLOAT})\s*")


def _raw_force_stress(log: str, *, symbols: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Read the unique final ABACUS force and stress tables, without cache data."""

    lines = log.splitlines()
    force_heads = [i for i, line in enumerate(lines) if "TOTAL-FORCE (eV/Angstrom)" in line]
    stress_heads = [i for i, line in enumerate(lines) if "TOTAL-STRESS (KBAR)" in line]
    if len(force_heads) != 1 or len(stress_heads) != 1 or force_heads[0] >= stress_heads[0]:
        raise ValueError("expected one ordered final force and stress table")
    force_rows = []
    labels = []
    for line in lines[force_heads[0] + 1:stress_heads[0]]:
        match = FORCE_ROW.fullmatch(line)
        if match:
            labels.append(match.group(1))
            force_rows.append([float(v.replace("D", "E").replace("d", "e"))
                               for v in match.groups()[1:]])
    stress_rows = []
    for line in lines[stress_heads[0] + 1:]:
        match = STRESS_ROW.fullmatch(line)
        if match:
            stress_rows.append([float(v.replace("D", "E").replace("d", "e"))
                                for v in match.groups()])
        elif stress_rows:
            break
    forces = np.asarray(force_rows, dtype=float)
    stress_kbar = np.asarray(stress_rows, dtype=float)
    if (forces.shape != (len(symbols), 3) or stress_kbar.shape != (3, 3)
            or not np.all(np.isfinite(forces)) or not np.all(np.isfinite(stress_kbar))
            or [re.sub(r"\d+$", "", label) for label in labels] != symbols
            or len(set(labels)) != len(symbols)
            or not np.allclose(stress_kbar, stress_kbar.T, atol=1e-8, rtol=0.0)):
        raise ValueError("raw ABACUS force/stress tables are incomplete or nonphysical")
    return forces, stress_kbar


def _gradient_from_raw(chart, coordinates: np.ndarray, forces: np.ndarray,
                       stress_eV_per_A3: np.ndarray) -> np.ndarray:
    """Rebuild the atomic-plus-six-strain dual gradient from raw DFT tables."""

    atoms = chart.to_atoms(coordinates)
    deformation = np.linalg.solve(chart.reference_cell, np.asarray(atoms.cell.array))
    atomic = -forces @ deformation.T
    cell = atoms.get_volume() * np.linalg.solve(
        deformation.T, 0.5 * (stress_eV_per_A3 + stress_eV_per_A3.T),
    )
    cell = 0.5 * (cell + cell.T)
    return np.concatenate([atomic.reshape(-1), np.array([
        cell[0, 0], cell[1, 1], cell[2, 2],
        2 * cell[1, 2], 2 * cell[0, 2], 2 * cell[0, 1],
    ])])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--force-constants", type=Path, required=True)
    parser.add_argument("--phonopy-eigenpairs", type=Path, required=True)
    parser.add_argument("--gamma-provenance", type=Path)
    parser.add_argument("--force-sets", type=Path)
    parser.add_argument("--eigenpairs-provenance", type=Path)
    parser.add_argument("--canary-audit", type=Path)
    parser.add_argument("--grid-result-manifest", type=Path, required=True)
    parser.add_argument("--summary-name", default="conditional_result.json")
    parser.add_argument("--curvature-audit", type=Path,
                        help="optional audited negative-curvature directions for branch projection")
    parser.add_argument("--refined-result", type=Path,
                        help="the stationary parent result paired with --curvature-audit")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite existing audit: {args.output}")
    if Path(args.summary_name).name != args.summary_name or not args.summary_name.endswith(".json"):
        raise ValueError("summary-name must be a .json basename")
    if (args.curvature_audit is None) != (args.refined_result is None):
        raise ValueError("curvature-audit and refined-result must be supplied together")
    summary_path = args.workdir / args.summary_name
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    grid = json.loads(args.grid_result_manifest.read_text(encoding="utf-8"))
    lock_soft_y = preflight.get("kind") == "bto_symmetry_restricted_soft_qy_zero_preflight_no_dft"
    physical_soft_plane = (preflight.get("kind") == "bto_transverse_soft_conditional_preflight_no_dft"
                           or lock_soft_y)
    expected_kind = ("bto_symmetry_restricted_qy_zero_variable_cell_local_candidate_not_PES_or_barrier"
                     if lock_soft_y else
                     "bto_transverse_soft_variable_cell_conditional_local_candidate_not_PES_or_barrier"
                     if physical_soft_plane else
                     "bto_fixed_q1q2_variable_cell_conditional_local_candidate_not_T_to_C_barrier")
    if (summary.get("status") not in {
                "orthogonal_gradient_converged_curvature_unchecked",
                "orthogonal_gradient_and_stress_converged_curvature_unchecked",
                "orthogonal_gradient_converged_stress_target_failed_curvature_unchecked",
            }
            or summary.get("kind") != expected_kind
            or summary.get("preflight_sha256") != sha256(args.preflight)):
        raise ValueError("conditional result or preflight hash is invalid")
    if physical_soft_plane:
        if any(path is None for path in (args.gamma_provenance, args.force_sets,
                                         args.eigenpairs_provenance)):
            raise ValueError("soft/soft audit requires the Gamma source")
        gamma_hashes = validate_bto_gamma_source(
            args.gamma_provenance, args.force_sets, args.eigenpairs_provenance,
        )
        if any(preflight["source_sha256"].get(key) != value for key, value in gamma_hashes.items()):
            raise ValueError("soft/soft Gamma source changed after preflight")
        if not lock_soft_y:
            if args.canary_audit is None:
                raise ValueError("open-Q_y pilot requires the independent three-start canary audit")
            canary = json.loads(args.canary_audit.read_text(encoding="utf-8"))
            if (canary.get("status") != "all_three_raw_DFT_points_validated"
                    or canary.get("source_sha256", {}).get("preflight") != sha256(args.preflight)):
                raise ValueError("independent three-start canary audit does not match this preflight")
    loaded = load_bto_q1q2_reference(
        args.report, args.reference, args.force_constants, args.phonopy_eigenpairs,
        strain_metric_weights_amu_A2=np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float),
    )
    chart = loaded.chart
    plane = (bto_transverse_soft_plane(
        loaded, strain_metric_weights_amu_A2=np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float),
    ) if physical_soft_plane else loaded.plane)
    third = (np.asarray(preflight["remaining_soft_y_metric_unit_direction"], dtype=float)
             if physical_soft_plane else None)
    frozen_directions = chart.rigid_translation_directions()
    if lock_soft_y:
        rebuilt = remaining_soft_y_direction(loaded.modes, plane)
        if (preflight.get("n_fixed_third_soft_axes") != 1
                or preflight.get("n_relaxed_orthogonal_coordinates") != chart.coordinate_count - 6
                or third.shape != rebuilt.shape
                or not np.allclose(third, rebuilt, rtol=0.0, atol=1e-10)
                or summary.get("third_soft_mode_constraint")
                   != "Q_y=0 at every evaluation; transverse stability is not implied"):
            raise ValueError("restricted third-soft-mode source or constraint is invalid")
        frozen_directions = np.column_stack([frozen_directions, third])
    open_directions = _orthogonal_directions(plane, frozen_directions)
    branch_basis = None
    branch_center = None
    branch_audit_hash = None
    if physical_soft_plane and args.curvature_audit is not None:
        raise ValueError("old soft/stable curvature directions cannot audit the transverse-soft plane")
    if args.curvature_audit is not None:
        branch_audit = json.loads(args.curvature_audit.read_text(encoding="utf-8"))
        refined = json.loads(args.refined_result.read_text(encoding="utf-8"))
        branch_seed = summary.get("branch_seed")
        if (branch_audit.get("status") != "negative_orthogonal_curvature_rejects_conditional_local_minimum"
                or branch_audit.get("source_sha256", {}).get("refined_result") != sha256(args.refined_result)
                or not isinstance(branch_seed, dict)
                or branch_seed.get("curvature_audit_sha256") != sha256(args.curvature_audit)
                or summary.get("start_from_result_sha256") != sha256(args.refined_result)
                or branch_seed.get("negative_eigenvalue_eV_per_amu_A2") >= 0):
            raise ValueError("branch summary does not match the audited negative-curvature parent")
        branch_center = np.asarray(refined["coordinates_u_A_eta_voigt"], dtype=float)
        branch_basis = np.column_stack([
            np.asarray(item["metric_normalized_chart_direction"], dtype=float)
            for item in branch_audit["negative_directions"]
        ])
        if (branch_center.shape != (chart.coordinate_count,)
                or branch_basis.shape[0] != chart.coordinate_count
                or not np.allclose(branch_basis.T @ (plane.metric_weights[:, None] * branch_basis),
                                   np.eye(branch_basis.shape[1]), rtol=0.0, atol=1e-6)
                or not np.allclose(branch_audit["q1_q2_sqrt_amu_A"], summary["q1_q2_sqrt_amu_A"],
                                   rtol=0.0, atol=1e-12)):
            raise ValueError("branch directions, center, metric or fixed Q changed")
        branch_audit_hash = sha256(args.curvature_audit)
    q_key = "q_parallel_q_transverse_sqrt_amu_A" if physical_soft_plane else "q1_q2_sqrt_amu_A"
    q = np.asarray(preflight[q_key], dtype=float)
    if not np.allclose(summary["q1_q2_sqrt_amu_A"], q, atol=1e-12, rtol=0.0):
        raise ValueError("summary changed the fixed Q coordinates")
    if physical_soft_plane and not np.allclose(summary[q_key], q, atol=1e-12, rtol=0.0):
        raise ValueError("soft/soft summary lacks matching physical-axis coordinates")
    final_coordinates = np.asarray(summary["coordinates_u_A_eta_voigt"], dtype=float)
    if not np.allclose(plane.project(final_coordinates), q, atol=1e-8, rtol=0.0):
        raise ValueError("final geometry drifted from fixed Q1/Q2")
    final_atoms = chart.to_atoms(final_coordinates)
    distances = final_atoms.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    final_minimum_distance = float(np.min(distances))
    if final_minimum_distance < preflight["minimum_allowed_atomic_distance_A"]:
        raise ValueError("final geometry violates the atomic-distance guard")
    cubic = [point for point in grid["points"] if abs(point["q1"]) < 1e-12 and abs(point["q2"]) < 1e-12]
    if len(cubic) != 1:
        raise ValueError("no unique cubic energy reference")
    reference_energy = float(cubic[0]["energy_eV"])
    expected_input = grid["points"][0]["input_sha256"]
    records = []
    frozen_start = plane.frozen_coordinates(q)
    for directory in sorted(args.workdir.glob("eval-*-*")):
        if not directory.is_dir():
            continue
        result_path = directory / "result.json"
        if not result_path.is_file():
            raise ValueError(f"incomplete or failed evaluation directory: {directory}")
        result = json.loads(result_path.read_text(encoding="utf-8"))
        evidence = result["validation"]
        log_path = directory / "OUT.ABACUS" / "running_scf.log"
        log = log_path.read_text(encoding="utf-8", errors="replace")
        if ("charge density convergence is achieved" not in log
                or "!FINAL_ETOT_IS" not in log or "PMI server not found" in log
                or sha256(log_path) != evidence["log_sha256"]):
            raise ValueError(f"SCF/log hash invalid: {directory}")
        mpi = re.findall(r"\bDSIZE\s*=\s*(\d+)", log)
        if mpi != ["32"] or evidence["mpi_dsize"] != 32:
            raise ValueError(f"MPI size invalid: {directory}")
        for name in ("INPUT", "KPT", "STRU"):
            if sha256(directory / name) != evidence["input_sha256"][name]:
                raise ValueError(f"input file hash changed: {directory / name}")
        for name in ("INPUT", "KPT"):
            if evidence["input_sha256"][name] != expected_input[name]:
                raise ValueError(f"conditional point changed the reviewed BTO {name}")
        coordinates = np.asarray(result["coordinates"], dtype=float)
        if coordinates.shape != final_coordinates.shape or not np.all(np.isfinite(coordinates)):
            raise ValueError(f"coordinate vector invalid: {directory}")
        digest = hashlib.sha256(coordinates.tobytes()).hexdigest()
        if (result.get("status") != "complete"
                or result.get("coordinate_sha256") != digest
                or result.get("contract_sha256") != summary["evaluator_contract_sha256"]):
            raise ValueError(f"result cache identity/contract invalid: {directory}")
        if not np.allclose(plane.project(coordinates), q, atol=1e-8, rtol=0.0):
            raise ValueError(f"evaluation drifted from fixed Q1/Q2: {directory}")
        if lock_soft_y and abs(float(third @ (plane.metric_weights * (coordinates - frozen_start)))) > 1e-8:
            raise ValueError(f"evaluation drifted from Q_y=0: {directory}")
        atoms = chart.to_atoms(coordinates)
        point_distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(point_distances, np.inf)
        point_minimum = float(np.min(point_distances))
        gradient = np.asarray(result["gradient"], dtype=float)
        raw_energy = FINAL_ENERGY.findall(log)
        raw_forces, raw_stress_kbar = _raw_force_stress(
            log, symbols=atoms.get_chemical_symbols(),
        )
        raw_stress = -raw_stress_kbar / 1602.176634
        raw_chart_gradient = _gradient_from_raw(chart, coordinates, raw_forces, raw_stress)
        if (gradient.shape != coordinates.shape or not np.all(np.isfinite(gradient))
                or point_minimum < preflight["minimum_allowed_atomic_distance_A"]
                or abs(point_minimum - result["minimum_distance_A"]) > 1e-8
                or abs(atoms.get_volume() - result["volume_A3"]) > 1e-8
                or not np.isfinite(result["maximum_atomic_force_eV_per_A"])
                or result["maximum_atomic_force_eV_per_A"] < 0.0
                or not np.isfinite(result["enthalpy_eV"])
                or np.asarray(result["stress_eV_per_A3"]).shape != (3, 3)
                or not np.all(np.isfinite(result["stress_eV_per_A3"]))
                or len(raw_energy) != 1
                or abs(float(raw_energy[0]) - result["enthalpy_eV"]) > 1e-6
                or abs(np.max(np.linalg.norm(raw_forces, axis=1))
                       - result["maximum_atomic_force_eV_per_A"]) > 1e-6
                or not np.allclose(raw_stress, result["stress_eV_per_A3"], atol=1e-8, rtol=0.0)
                or not np.allclose(raw_chart_gradient, gradient, atol=1e-6, rtol=0.0)):
            raise ValueError(f"energy/gradient/stress/geometry invalid: {directory}")
        item = {
            "directory": directory.name,
            "raw_log_sha256": evidence["log_sha256"],
            "raw_input_sha256": evidence["input_sha256"],
            "energy_minus_c_eV_per_BTO": result["enthalpy_eV"] - reference_energy,
            "orthogonal_gradient_norm_eV_per_sqrt_amu_A": float(
                np.linalg.norm(open_directions.T @ gradient)
            ),
            "minimum_distance_A": result["minimum_distance_A"],
            "maximum_atomic_force_eV_per_A": result["maximum_atomic_force_eV_per_A"],
            "maximum_absolute_stress_kbar": float(np.max(np.abs(result["stress_eV_per_A3"])) * 1602.176634),
            "mpi_dsize": 32,
            "slurm_job_id": evidence["slurm_job_id"],
        }
        if third is not None:
            item["third_soft_y_amplitude_sqrt_amu_A"] = float(
                third @ (plane.metric_weights * (coordinates - frozen_start))
            )
        if branch_basis is not None:
            displacement = coordinates - branch_center
            item["negative_curvature_amplitudes_sqrt_amu_A"] = (
                branch_basis.T @ (plane.metric_weights * displacement)
            ).tolist()
        records.append(item)
    if not records:
        raise ValueError("conditional pilot has no completed DFT evaluations")
    if physical_soft_plane:
        for seed in preflight["branch_starts"]:
            expected = np.asarray(seed["coordinates_u_A_eta_voigt"], dtype=float)
            found = [directory for directory in args.workdir.glob("eval-*-*")
                     if (directory / "result.json").is_file() and (
                         np.allclose(
                             np.asarray(json.loads((directory / "result.json").read_text(encoding="utf-8"))["coordinates"],
                                        dtype=float), expected, rtol=0.0, atol=1e-10)
                         if lock_soft_y else np.array_equal(
                             np.asarray(json.loads((directory / "result.json").read_text(encoding="utf-8"))["coordinates"],
                                        dtype=float), expected))]
            if len(found) != 1:
                raise ValueError(f"a reviewed soft/soft branch start is absent: {seed['label']}")
    matches = [record for record in records if (
        (args.workdir / record["directory"] / "result.json").is_file()
        and np.array_equal(np.asarray(json.loads((args.workdir / record["directory"] / "result.json").read_text())
                                          ["coordinates"], dtype=float), final_coordinates)
    )]
    if len(matches) != 1 or abs(matches[0]["energy_minus_c_eV_per_BTO"] - summary["energy_minus_c_eV_per_BTO"]) > 1e-8:
        raise ValueError("final optimizer result is not an audited DFT point")
    if summary.get("final_evaluation_directory") is not None and matches[0]["directory"] != summary["final_evaluation_directory"]:
        raise ValueError("reported final DFT evaluation directory differs from audited point")
    if summary.get("maximum_absolute_stress_kbar") is not None and not np.isclose(
        matches[0]["maximum_absolute_stress_kbar"], summary["maximum_absolute_stress_kbar"], atol=1e-8, rtol=0.0
    ):
        raise ValueError("reported final raw stress differs from audited DFT result")
    stress_target = summary.get("stress_target_kbar")
    if stress_target is not None and bool(summary.get("stress_target_passed")) != (
        matches[0]["maximum_absolute_stress_kbar"] <= stress_target
    ):
        raise ValueError("reported stress-target status differs from audited DFT result")
    if summary["orthogonal_gradient_norm_eV_per_sqrt_amu_A"] > summary["gradient_tolerance_eV_per_sqrt_amu_A"]:
        raise ValueError("reported orthogonal gradient did not converge")
    audited_branch_outcomes = None
    if summary.get("branch_outcomes") is not None:
        outcomes = summary["branch_outcomes"]
        labels = summary.get("branch_start_labels")
        if (not isinstance(outcomes, list) or len(outcomes) != summary["n_starts"]
                or not isinstance(labels, list) or len(labels) != len(outcomes)
                or not 0 <= summary["selected_start"] < len(outcomes)):
            raise ValueError("branch-outcome count or labels differ from the optimizer starts")
        audited_branch_outcomes = []
        for index, (label, outcome) in enumerate(zip(labels, outcomes)):
            basename = outcome["final_evaluation_directory"]
            if (Path(basename).name != basename
                    or not any(row["directory"] == basename for row in records)):
                raise ValueError(f"branch terminal point is not in the audited DFT cache: {label}")
            directory = args.workdir / basename
            cached = json.loads((directory / "result.json").read_text(encoding="utf-8"))
            coords = np.asarray(outcome["coordinates_u_A_eta_voigt"], dtype=float)
            projected = float(np.linalg.norm(open_directions.T @ np.asarray(cached["gradient"], dtype=float)))
            converged = projected <= summary["gradient_tolerance_eV_per_sqrt_amu_A"]
            if (outcome["label"] != label or outcome["start_index"] != index
                    or coords.shape != final_coordinates.shape
                    or not np.array_equal(coords, np.asarray(cached["coordinates"], dtype=float))
                    or abs(outcome["energy_minus_c_eV_per_BTO"]
                           - (cached["enthalpy_eV"] - reference_energy)) > 1e-8
                    or abs(outcome["orthogonal_gradient_norm_eV_per_sqrt_amu_A"] - projected) > 1e-8
                    or bool(outcome["converged"]) != converged):
                raise ValueError(f"branch outcome differs from an audited DFT point: {label}")
            audited_branch_outcomes.append({
                "label": label,
                "converged": converged,
                "energy_minus_c_eV_per_BTO": cached["enthalpy_eV"] - reference_energy,
                "orthogonal_gradient_norm_eV_per_sqrt_amu_A": projected,
                "final_evaluation_directory": directory.name,
            })
        chosen = audited_branch_outcomes[summary["selected_start"]]
        if (not chosen["converged"] or chosen["final_evaluation_directory"] != matches[0]["directory"]
                or any(row["converged"] and row["energy_minus_c_eV_per_BTO"]
                       < chosen["energy_minus_c_eV_per_BTO"] - 1e-8
                       for row in audited_branch_outcomes)):
            raise ValueError("selected branch is not the lowest audited converged branch")
    final_branch_projection = None
    branch_job_points = None
    signed_basin_candidates = None
    if branch_basis is not None:
        final_branch_projection = matches[0]["negative_curvature_amplitudes_sqrt_amu_A"]
        branch_job_points = [item for item in records if item["slurm_job_id"] == summary["slurm_job_id"]]
        if not branch_job_points or not np.all(np.isfinite(final_branch_projection)):
            raise ValueError("branch result lacks completed DFT points from its own Slurm job")
        index = summary["branch_seed"]["negative_direction_index"]
        amplitude = float(summary["branch_seed"]["seed_amplitude_sqrt_amu_A"])
        if index < 0 or index >= branch_basis.shape[1] or amplitude <= 0.0:
            raise ValueError("branch seed direction or amplitude is invalid")
        signed = [item["negative_curvature_amplitudes_sqrt_amu_A"][index]
                  for item in branch_job_points]
        if not (any(value > 0.0 for value in signed) and any(value < 0.0 for value in signed)):
            raise ValueError("the reported two-sign branch job did not evaluate both signed basins")
        signed_basin_candidates = {}
        for label, selector in (("positive", lambda value: value > 0.0),
                                ("negative", lambda value: value < 0.0)):
            eligible = [item for item in branch_job_points
                        if selector(item["negative_curvature_amplitudes_sqrt_amu_A"][index])
                        and item["orthogonal_gradient_norm_eV_per_sqrt_amu_A"]
                        <= summary["gradient_tolerance_eV_per_sqrt_amu_A"]
                        and (stress_target is None or item["maximum_absolute_stress_kbar"] <= stress_target)]
            signed_basin_candidates[label] = (
                None if not eligible else min(eligible, key=lambda item: item["energy_minus_c_eV_per_BTO"])
            )
        all_eligible = [item for item in signed_basin_candidates.values() if item is not None]
        if all_eligible and min(item["energy_minus_c_eV_per_BTO"] for item in all_eligible) < (
            summary["energy_minus_c_eV_per_BTO"] - 1e-8
        ):
            raise ValueError("reported selected branch is higher than an audited converged signed candidate")
    output = {
        "kind": ("bto_symmetry_restricted_qy_zero_single_point_audit_not_PES_or_barrier"
                 if lock_soft_y else
                 "bto_transverse_soft_conditional_single_point_audit_not_PES_or_barrier"
                 if physical_soft_plane else "bto_q1q2_conditional_single_point_audit_not_PES_or_barrier"),
        "status": ("verified_gradient_stationary_stress_target_failed_curvature_unchecked"
                   if summary.get("stress_target_passed") is False else
                   "verified_gradient_stationary_branch_candidate_local_curvature_unchecked"
                   if branch_basis is not None else
                   "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"),
        "q1_q2_sqrt_amu_A": q.tolist(),
        "axis_kind": preflight.get("axis_kind", "archived_soft_stable"),
        "q_parallel_q_transverse_sqrt_amu_A": q.tolist() if physical_soft_plane else None,
        "final_energy_minus_c_eV_per_BTO": summary["energy_minus_c_eV_per_BTO"],
        "orthogonal_gradient_norm_eV_per_sqrt_amu_A": summary["orthogonal_gradient_norm_eV_per_sqrt_amu_A"],
        "final_cell_strain_voigt": final_coordinates[-6:].tolist(),
        "final_minimum_distance_A": final_minimum_distance,
        "final_maximum_absolute_stress_kbar": matches[0]["maximum_absolute_stress_kbar"],
        "stress_target_kbar": stress_target,
        "stress_target_passed": summary.get("stress_target_passed"),
        "n_individually_audited_DFT_points": len(records),
        "evaluations": records,
        "final_third_soft_y_amplitude_sqrt_amu_A": (
            None if third is None else float(third @ (plane.metric_weights * (final_coordinates - frozen_start)))
        ),
        "third_soft_mode_restricted_at_zero": lock_soft_y,
        "final_negative_curvature_amplitudes_sqrt_amu_A": final_branch_projection,
        "n_current_branch_job_DFT_points": None if branch_job_points is None else len(branch_job_points),
        "signed_basin_converged_candidates": signed_basin_candidates,
        "branch_outcomes": audited_branch_outcomes,
        "branch_outcomes_preserved_in_source": audited_branch_outcomes is not None,
        "source_sha256": {
            "auditor": sha256(Path(__file__)),
            "summary": sha256(summary_path),
            "preflight": sha256(args.preflight),
            "grid_result_manifest": sha256(args.grid_result_manifest),
            "curvature_audit": branch_audit_hash,
            "canary_audit": sha256(args.canary_audit) if args.canary_audit is not None else None,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: output[key] for key in (
        "status", "final_energy_minus_c_eV_per_BTO",
        "orthogonal_gradient_norm_eV_per_sqrt_amu_A", "n_individually_audited_DFT_points",
    )}, indent=2))


if __name__ == "__main__":
    main()
