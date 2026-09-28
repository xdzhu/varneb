"""Independently audit a nonconverged BTO Qy=0 cache before warm continuation.

This is deliberately separate from the stationary-point audit: a complete
electronic calculation is not proof that the conditional geometry converged.
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
    for name in ("workdir", "preflight", "report", "reference",
                 "force_constants", "phonopy_eigenpairs", "gamma_provenance",
                 "force_sets", "eigenpairs_provenance", "grid_result_manifest",
                 "output"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    parser.add_argument("--expected-job-id", required=True)
    parser.add_argument("--gradient-tolerance", type=float, required=True)
    parser.add_argument("--orthogonal-amplitude-bound", type=float, required=True)
    parser.add_argument("--stress-target-kbar", type=float, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite audit: {args.output}")
    if (args.gradient_tolerance <= 0 or args.orthogonal_amplitude_bound <= 0
            or args.stress_target_kbar <= 0
            or not np.isfinite(args.gradient_tolerance)
            or not np.isfinite(args.orthogonal_amplitude_bound)
            or not np.isfinite(args.stress_target_kbar)):
        raise ValueError("invalid optimizer tolerance or bound")
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    if preflight.get("kind") != "bto_symmetry_restricted_soft_qy_zero_preflight_no_dft":
        raise ValueError("only the fixed-Qy=0 pilot can use this failed-cache audit")
    source_paths = {
        "report": args.report,
        "reference": args.reference,
        "force_constants": args.force_constants,
        "phonopy_eigenpairs": args.phonopy_eigenpairs,
        "gamma_provenance": args.gamma_provenance,
        "force_sets": args.force_sets,
        "eigenpairs_provenance": args.eigenpairs_provenance,
        "grid_result_manifest": args.grid_result_manifest,
    }
    if any(sha256(path) != preflight["source_sha256"][name]
           for name, path in source_paths.items()):
        raise ValueError("the preflight source hashes changed")
    validate_bto_gamma_source(
        args.gamma_provenance, args.force_sets, args.eigenpairs_provenance,
    )
    weights = np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float)
    loaded = load_bto_q1q2_reference(
        args.report, args.reference, args.force_constants,
        args.phonopy_eigenpairs, strain_metric_weights_amu_A2=weights,
    )
    chart = loaded.chart
    plane = bto_transverse_soft_plane(
        loaded, strain_metric_weights_amu_A2=weights,
    )
    third = np.asarray(preflight["remaining_soft_y_metric_unit_direction"], dtype=float)
    if not np.allclose(third, remaining_soft_y_direction(loaded.modes, plane), atol=1e-10, rtol=0):
        raise ValueError("the fixed third soft-mode direction changed")
    frozen = np.column_stack([chart.rigid_translation_directions(), third])
    open_directions = _orthogonal_directions(plane, frozen)
    base = plane.frozen_coordinates(
        np.asarray(preflight["q_parallel_q_transverse_sqrt_amu_A"], dtype=float),
    )
    if (open_directions.shape[1] != 15
            or preflight["n_relaxed_orthogonal_coordinates"] != 15):
        raise ValueError("the restricted subspace has the wrong dimension")
    contract_path = args.workdir / "contract.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if (contract["calculator_id"] != preflight["calculator_id"]
            or contract["minimum_distance_A"] != preflight["minimum_allowed_atomic_distance_A"]
            or not contract["validation_id"].startswith("bto-abacus-converged-static-v1:")
            or not isinstance(contract["contract_sha256"], str)):
        raise ValueError("evaluator contract differs from the preflight")
    grid = json.loads(args.grid_result_manifest.read_text(encoding="utf-8"))
    if grid["status"] != "converged" or grid["mpi_ranks"] != 32:
        raise ValueError("the common BTO static reference is invalid")
    cubic = [row for row in grid["points"] if abs(row["q1"]) < 1e-12
             and abs(row["q2"]) < 1e-12]
    if len(cubic) != 1:
        raise ValueError("the cubic energy reference is not unique")
    input_hash = grid["points"][0]["input_sha256"]
    directories = sorted(args.workdir.glob("eval-*-*"))
    if not directories or any(not item.is_dir() for item in directories):
        raise ValueError("there are no complete evaluation directories")
    rows = []
    for index, directory in enumerate(directories):
        if not directory.name.startswith(f"eval-{index:06d}-"):
            raise ValueError("the failed run has missing or reordered evaluations")
        result_path = directory / "result.json"
        result = json.loads(result_path.read_text(encoding="utf-8"))
        evidence = result["validation"]
        log_path = directory / "OUT.ABACUS" / "running_scf.log"
        log = log_path.read_text(encoding="utf-8", errors="replace")
        if (result.get("status") != "complete"
                or result.get("contract_sha256") != contract["contract_sha256"]
                or evidence.get("slurm_job_id") != args.expected_job_id
                or evidence.get("mpi_dsize") != 32
                or re.findall(r"\bDSIZE\s*=\s*(\d+)", log) != ["32"]
                or "charge density convergence is achieved" not in log
                or "PMI server not found" in log
                or sha256(log_path) != evidence["log_sha256"]):
            raise ValueError(f"incomplete SCF, MPI or cache contract: {directory}")
        for name in ("INPUT", "KPT", "STRU"):
            if sha256(directory / name) != evidence["input_sha256"][name]:
                raise ValueError(f"input hash changed: {directory / name}")
        if any(evidence["input_sha256"][name] != input_hash[name]
               for name in ("INPUT", "KPT")):
            raise ValueError(f"common ABACUS INPUT/KPT changed: {directory}")
        coordinates = np.asarray(result["coordinates"], dtype=float)
        gradient = np.asarray(result["gradient"], dtype=float)
        if (coordinates.shape != (chart.coordinate_count,)
                or gradient.shape != coordinates.shape
                or not np.all(np.isfinite(coordinates))
                or not np.all(np.isfinite(gradient))
                or hashlib.sha256(coordinates.tobytes()).hexdigest()
                   != result["coordinate_sha256"]
                or not np.allclose(plane.project(coordinates),
                                   preflight["q_parallel_q_transverse_sqrt_amu_A"],
                                   atol=1e-8, rtol=0)
                or abs(float(third @ (plane.metric_weights * (coordinates - base)))) > 1e-8):
            raise ValueError(f"coordinates or fixed modes changed: {directory}")
        atoms = chart.to_atoms(coordinates)
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        raw_energy = FINAL_ENERGY.findall(log)
        forces, raw_stress_kbar = _raw_force_stress(
            log, symbols=atoms.get_chemical_symbols(),
        )
        raw_stress = -raw_stress_kbar / 1602.176634
        if (len(raw_energy) != 1
                or abs(float(raw_energy[0]) - result["enthalpy_eV"]) > 1e-6
                or abs(float(np.min(distances)) - result["minimum_distance_A"]) > 1e-8
                or np.min(distances) < preflight["minimum_allowed_atomic_distance_A"]
                or abs(atoms.get_volume() - result["volume_A3"]) > 1e-8
                or abs(np.max(np.linalg.norm(forces, axis=1))
                       - result["maximum_atomic_force_eV_per_A"]) > 1e-6
                or not np.allclose(raw_stress, result["stress_eV_per_A3"], atol=1e-8, rtol=0)
                or not np.allclose(_gradient_from_raw(chart, coordinates, forces, raw_stress),
                                   gradient, atol=1e-6, rtol=0)):
            raise ValueError(f"raw energy/force/stress/geometry mismatch: {directory}")
        z = open_directions.T @ (plane.metric_weights * (coordinates - base))
        rows.append({
            "directory": directory.name,
            "result_sha256": sha256(result_path),
            "energy_minus_c_eV_per_BTO": float(result["enthalpy_eV"] - cubic[0]["energy_eV"]),
            "orthogonal_gradient_norm_eV_per_sqrt_amu_A": float(np.linalg.norm(open_directions.T @ gradient)),
            "maximum_absolute_stress_kbar": float(np.max(np.abs(raw_stress_kbar))),
            "maximum_absolute_orthogonal_amplitude_sqrt_amu_A": float(np.max(np.abs(z))),
            "minimum_distance_A": float(np.min(distances)),
        })
    candidate = rows[-1]
    if (candidate["orthogonal_gradient_norm_eV_per_sqrt_amu_A"] <= args.gradient_tolerance
            or candidate["maximum_absolute_orthogonal_amplitude_sqrt_amu_A"]
               < args.orthogonal_amplitude_bound - 1e-5
            or candidate["maximum_absolute_orthogonal_amplitude_sqrt_amu_A"]
               > args.orthogonal_amplitude_bound + 1e-5
            or min(row["energy_minus_c_eV_per_BTO"] for row in rows)
               < candidate["energy_minus_c_eV_per_BTO"] - 1e-8):
        raise ValueError("the final cache is not a nonstationary bound-hit candidate")
    output = {
        "kind": "bto_qy0_failed_optimizer_cache_raw_audit_not_stationary_not_PES",
        "status": "all_raw_DFT_points_validated_nonstationary_bound_hit",
        "q_parallel_q_transverse_sqrt_amu_A": preflight["q_parallel_q_transverse_sqrt_amu_A"],
        "n_individually_audited_DFT_points": len(rows),
        "expected_job_id": args.expected_job_id,
        "gradient_tolerance_eV_per_sqrt_amu_A": args.gradient_tolerance,
        "orthogonal_amplitude_bound_sqrt_amu_A": args.orthogonal_amplitude_bound,
        "stress_target_kbar": args.stress_target_kbar,
        "candidate_evaluation_directory": candidate["directory"],
        "candidate_result_sha256": candidate["result_sha256"],
        "candidate_energy_minus_c_eV_per_BTO": candidate["energy_minus_c_eV_per_BTO"],
        "candidate_orthogonal_gradient_norm_eV_per_sqrt_amu_A": candidate["orthogonal_gradient_norm_eV_per_sqrt_amu_A"],
        "candidate_maximum_absolute_stress_kbar": candidate["maximum_absolute_stress_kbar"],
        "candidate_maximum_absolute_orthogonal_amplitude_sqrt_amu_A": candidate["maximum_absolute_orthogonal_amplitude_sqrt_amu_A"],
        "contract_sha256": contract["contract_sha256"],
        "source_sha256": {
            "auditor": sha256(Path(__file__)),
            "preflight": sha256(args.preflight),
            "grid_result_manifest": sha256(args.grid_result_manifest),
            "evaluator_contract": sha256(contract_path),
        },
        "evaluations": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: output[key] for key in (
        "status", "n_individually_audited_DFT_points",
        "candidate_orthogonal_gradient_norm_eV_per_sqrt_amu_A",
        "candidate_maximum_absolute_stress_kbar",
        "candidate_maximum_absolute_orthogonal_amplitude_sqrt_amu_A",
    )}, indent=2))


if __name__ == "__main__":
    main()
