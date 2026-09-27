"""Independently audit four BTO soft-direction static DFT probe outputs.

Reads original INPUT/KPT/STRU and ABACUS SCF logs plus the evaluator cache.
It does not run a calculator and does not certify a 2D conditional PES.
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
    bto_transverse_soft_plane,
    load_bto_q1q2_reference,
    sha256,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("workdir", "preflight", "soft_preflight", "branch_result", "report",
                 "reference", "force_constants", "phonopy_eigenpairs",
                 "grid_result_manifest", "summary", "output"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    parser.add_argument("--all-points-audit", type=Path,
                        help="required independent full raw DFT audit for transverse-soft probes")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite audit: {args.output}")
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    soft = json.loads(args.soft_preflight.read_text(encoding="utf-8"))
    branch = json.loads(args.branch_result.read_text(encoding="utf-8"))
    grid = json.loads(args.grid_result_manifest.read_text(encoding="utf-8"))
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    physical_soft_plane = preflight.get("kind") == "bto_transverse_soft_conditional_preflight_no_dft"
    if (preflight.get("kind") not in {
            "bto_q1q2_conditional_preflight_no_dft",
            "bto_transverse_soft_conditional_preflight_no_dft",
        }
            or soft.get("kind") != (
                "bto_transverse_soft_mixed_eigenvector_probe_preflight_no_dft"
                if physical_soft_plane else
                "bto_fixed_q1q2_soft_mixed_eigenvector_probe_preflight_no_dft")
            or summary.get("kind") != (
                "bto_transverse_soft_mixed_direction_four_DFT_probe_crosscheck_not_PES_or_barrier"
                if physical_soft_plane else
                "bto_fixed_q1q2_soft_mixed_direction_four_DFT_probe_crosscheck_not_PES_or_barrier")
            or summary.get("status") != "four_static_probes_complete_pending_independent_SCF_audit"
            or summary.get("preflight_sha256") != sha256(args.soft_preflight)
            or summary.get("evaluator_contract_sha256") != branch["evaluator_contract_sha256"]
            or soft.get("source_sha256", {}).get("branch_result") != sha256(args.branch_result)
            or soft.get("source_sha256", {}).get("preflight") != sha256(args.preflight)
            or len(soft.get("probes", [])) != 4 or len(summary.get("probes", [])) != 4
            or summary.get("n_new_evaluations") != 4):
        raise ValueError("four-probe result does not match the audited BTO branch/preflight")
    raw_audited = None
    if physical_soft_plane:
        if args.all_points_audit is None:
            raise ValueError("transverse-soft probes require the full independent raw DFT audit")
        raw_audited = json.loads(args.all_points_audit.read_text(encoding="utf-8"))
        names = {item["directory"] for item in raw_audited["evaluations"]}
        if (raw_audited.get("status") != "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"
                or raw_audited.get("source_sha256", {}).get("summary") != sha256(args.branch_result)
                or raw_audited.get("n_individually_audited_DFT_points") != len(names)
                or len(names) < 144):
            raise ValueError("full original ABACUS force/stress audit is incomplete")
    elif args.all_points_audit is not None:
        raise ValueError("full raw audit argument is reserved for the transverse-soft case")
    loaded = load_bto_q1q2_reference(
        args.report, args.reference, args.force_constants, args.phonopy_eigenpairs,
        strain_metric_weights_amu_A2=np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float),
    )
    if any(soft.get("source_sha256", {}).get(key) != digest
           for key, digest in loaded.source_hashes.items()):
        raise ValueError("mode/reference input hashes differ from soft-probe preflight")
    chart = loaded.chart
    plane = (bto_transverse_soft_plane(
        loaded,
        strain_metric_weights_amu_A2=np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float),
    ) if physical_soft_plane else loaded.plane)
    q = np.asarray(soft["q1_q2_sqrt_amu_A"], dtype=float)
    center = np.asarray(soft["center_coordinates_u_A_eta_voigt"], dtype=float)
    direction = np.asarray(soft["soft_direction_metric_unit_chart"], dtype=float)
    expected_input = grid["points"][0]["input_sha256"]
    original = json.loads((args.workdir / branch["final_evaluation_directory"] / "result.json").read_text(encoding="utf-8"))
    center_energy = float(original["enthalpy_eV"])
    center_gradient = np.asarray(original["gradient"], dtype=float)
    mode_amplitudes = plane.modal_amplitudes(direction)
    if (center_gradient.shape != center.shape or not np.all(np.isfinite(center_gradient))
            or not np.allclose(mode_amplitudes, soft["soft_direction_mode_amplitudes"],
                               rtol=0.0, atol=1e-9)
            or abs(float(np.dot(mode_amplitudes, mode_amplitudes)) - 1.0) > 1e-8):
        raise ValueError("center gradient or soft-mode decomposition differs from preflight")
    center_directional_gradient = float(direction @ center_gradient)
    center_log_path = args.workdir / branch["final_evaluation_directory"] / "OUT.ABACUS" / "running_scf.log"
    center_log = center_log_path.read_text(encoding="utf-8", errors="replace")
    center_markers = re.findall(r"!FINAL_ETOT_IS\s+([-+\d.eE]+)\s+eV", center_log)
    if len(center_markers) != 1:
        raise ValueError("branch center has no unique full-precision ABACUS final energy marker")
    center_energy_marker = float(center_markers[0])
    if (not np.array_equal(np.asarray(original["coordinates"], dtype=float), center)
            or summary["center_result_sha256"] != sha256(args.workdir / branch["final_evaluation_directory"] / "result.json")
            or abs(center_energy - summary["center_energy_eV"]) > 1e-10):
        raise ValueError("central branch energy/coordinates changed")
    candidates = []
    for directory in sorted(args.workdir.glob("eval-*-*")):
        path = directory / "result.json"
        if not path.is_file():
            continue
        record = json.loads(path.read_text(encoding="utf-8"))
        if record.get("validation", {}).get("slurm_job_id") == str(summary["slurm_job_id"]):
            candidates.append((directory, record))
    if len(candidates) != 4:
        raise ValueError("the Slurm job must have exactly four completed cache points")
    if raw_audited is not None and any(directory.name not in names for directory, _ in candidates):
        raise ValueError("a new soft-direction DFT point is absent from the independent raw audit")
    audited = []
    for item in soft["probes"]:
        step, sign = float(item["step_sqrt_amu_A"]), int(item["sign"])
        expected = center + sign * step * direction
        matches = [(directory, record) for directory, record in candidates
                   if np.array_equal(np.asarray(record["coordinates"], dtype=float), expected)]
        if len(matches) != 1:
            raise ValueError("each signed preflight probe must match exactly one new DFT point")
        directory, record = matches[0]
        coords = np.asarray(record["coordinates"], dtype=float)
        gradient = np.asarray(record["gradient"], dtype=float)
        evidence = record["validation"]
        log_path = directory / "OUT.ABACUS" / "running_scf.log"
        log = log_path.read_text(encoding="utf-8", errors="replace")
        sizes = re.findall(r"\bDSIZE\s*=\s*(\d+)", log)
        energy_markers = re.findall(r"!FINAL_ETOT_IS\s+([-+\d.eE]+)\s+eV", log)
        if (record.get("status") != "complete"
                or record.get("contract_sha256") != branch["evaluator_contract_sha256"]
                or record.get("coordinate_sha256") != hashlib.sha256(coords.tobytes()).hexdigest()
                or gradient.shape != coords.shape or not np.all(np.isfinite(gradient))
                or not np.isfinite(float(record["enthalpy_eV"]))
                or np.asarray(record["stress_eV_per_A3"]).shape != (3, 3)
                or not np.all(np.isfinite(record["stress_eV_per_A3"]))
                or not np.isfinite(float(record["maximum_atomic_force_eV_per_A"]))
                or "charge density convergence is achieved" not in log
                or "!FINAL_ETOT_IS" not in log or "PMI server not found" in log
                or len(energy_markers) != 1
                or sizes != ["32"] or evidence.get("mpi_dsize") != 32
                or evidence.get("log_sha256") != sha256(log_path)
                or not np.allclose(plane.project(coords), q, rtol=0.0, atol=1e-8)):
            raise ValueError(f"energy/force/stress/SCF/MPI/cache audit failed: {directory}")
        for name in ("INPUT", "KPT", "STRU"):
            if sha256(directory / name) != evidence["input_sha256"][name]:
                raise ValueError(f"original DFT input differs from cached input hash: {directory / name}")
        for name in ("INPUT", "KPT"):
            if evidence["input_sha256"][name] != expected_input[name]:
                raise ValueError(f"DFT {name} differs from reviewed BTO grid: {directory}")
        atoms = chart.to_atoms(coords)
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        minimum = float(np.min(distances))
        if (minimum < preflight["minimum_allowed_atomic_distance_A"]
                or abs(minimum - record["minimum_distance_A"]) > 1e-8
                or abs(float(atoms.get_volume()) - record["volume_A3"]) > 1e-8):
            raise ValueError(f"geometry and original DFT record differ: {directory}")
        directional_gradient = float(direction @ gradient)
        atomic_count = 3 * chart.n_atoms
        atomic_contribution = float(direction[:atomic_count] @ gradient[:atomic_count])
        strain_contribution = float(direction[atomic_count:] @ gradient[atomic_count:])
        marker_energy = float(energy_markers[0])
        if (abs(marker_energy - float(record["enthalpy_eV"])) > 1e-6
                or abs(atomic_contribution + strain_contribution - directional_gradient) > 1e-12):
            raise ValueError("full-precision energy or directional gradient decomposition differs")
        reported = [value for value in summary["probes"]
                    if value["step_sqrt_amu_A"] == step and value["sign"] == sign]
        if (len(reported) != 1
                or abs(float(record["enthalpy_eV"]) - reported[0]["energy_eV"]) > 1e-10
                or abs(directional_gradient - reported[0]["directional_gradient_eV_per_sqrt_amu_A"]) > 1e-10):
            raise ValueError("runner summary differs from original DFT point")
        audited.append({
            "step_sqrt_amu_A": step, "sign": sign, "directory": directory.name,
            "energy_eV": float(record["enthalpy_eV"]),
            "full_precision_log_energy_eV": marker_energy,
            "directional_gradient_eV_per_sqrt_amu_A": directional_gradient,
            "atomic_directional_gradient_eV_per_sqrt_amu_A": atomic_contribution,
            "strain_directional_gradient_eV_per_sqrt_amu_A": strain_contribution,
            "maximum_atomic_force_eV_per_A": float(record["maximum_atomic_force_eV_per_A"]),
            "maximum_absolute_stress_kbar": float(np.max(np.abs(record["stress_eV_per_A3"])) * 1602.176634),
            "minimum_atomic_distance_A": minimum, "result_sha256": sha256(directory / "result.json"),
            "log_sha256": sha256(log_path),
            "input_sha256": evidence["input_sha256"],
        })
    step_values = sorted({entry["step_sqrt_amu_A"] for entry in audited})
    curvatures = []
    for step in step_values:
        plus = next(entry for entry in audited if entry["step_sqrt_amu_A"] == step and entry["sign"] == 1)
        minus = next(entry for entry in audited if entry["step_sqrt_amu_A"] == step and entry["sign"] == -1)
        energy_curvature = (plus["energy_eV"] + minus["energy_eV"] - 2 * center_energy) / step**2
        full_precision_energy_curvature = (plus["full_precision_log_energy_eV"]
                                           + minus["full_precision_log_energy_eV"]
                                           - 2 * center_energy_marker) / step**2
        gradient_curvature = (plus["directional_gradient_eV_per_sqrt_amu_A"]
                              - minus["directional_gradient_eV_per_sqrt_amu_A"]) / (2 * step)
        atomic_curvature = (plus["atomic_directional_gradient_eV_per_sqrt_amu_A"]
                            - minus["atomic_directional_gradient_eV_per_sqrt_amu_A"]) / (2 * step)
        strain_curvature = (plus["strain_directional_gradient_eV_per_sqrt_amu_A"]
                            - minus["strain_directional_gradient_eV_per_sqrt_amu_A"]) / (2 * step)
        reported = next(value for value in summary["curvatures"] if value["step_sqrt_amu_A"] == step)
        if (abs(energy_curvature - reported["energy_curvature_eV_per_amu_A2"]) > 1e-9
                or abs(gradient_curvature - reported["gradient_curvature_eV_per_amu_A2"]) > 1e-9):
            raise ValueError("independent energy/gradient curvature differs from runner summary")
        curvatures.append({"step_sqrt_amu_A": step,
                           "energy_curvature_eV_per_amu_A2": energy_curvature,
                           "full_precision_log_energy_curvature_eV_per_amu_A2": full_precision_energy_curvature,
                           "gradient_curvature_eV_per_amu_A2": gradient_curvature,
                           "atomic_gradient_curvature_eV_per_amu_A2": atomic_curvature,
                           "strain_gradient_curvature_eV_per_amu_A2": strain_curvature,
                           "energy_gradient_disagreement_eV_per_amu_A2": abs(energy_curvature - gradient_curvature)})
    output = {
        "kind": ("bto_transverse_soft_mixed_direction_four_DFT_probe_independent_audit_not_PES_or_barrier"
                 if physical_soft_plane else
                 "bto_fixed_q1q2_soft_mixed_direction_four_DFT_probe_independent_audit_not_PES_or_barrier"),
        "status": "four_original_static_DFT_points_verified_soft_direction_curvature_screen_only",
        "q1_q2_sqrt_amu_A": q.tolist(), "slurm_job_id": str(summary["slurm_job_id"]),
        "center_full_precision_log_energy_eV": center_energy_marker,
        "center_soft_directional_gradient_eV_per_sqrt_amu_A": center_directional_gradient,
        "soft_direction_metric_fraction": {
            "first_Gamma_triplet": float(np.sum(mode_amplitudes[:3] ** 2)),
            "second_Gamma_triplet": float(np.sum(mode_amplitudes[3:6] ** 2)),
            "third_Gamma_triplet": float(np.sum(mode_amplitudes[6:9] ** 2)),
            "strain": float(np.sum(mode_amplitudes[-6:] ** 2)),
        },
        "center_log_sha256": sha256(center_log_path),
        "n_individually_audited_DFT_points": len(audited),
        "probes": audited, "curvatures": curvatures,
        "energy_curvature_step_disagreement_eV_per_amu_A2": abs(curvatures[0]["energy_curvature_eV_per_amu_A2"] - curvatures[1]["energy_curvature_eV_per_amu_A2"]),
        "gradient_curvature_step_disagreement_eV_per_amu_A2": abs(curvatures[0]["gradient_curvature_eV_per_amu_A2"] - curvatures[1]["gradient_curvature_eV_per_amu_A2"]),
        "source_sha256": {key: sha256(path) for key, path in {
            "summary": args.summary, "soft_preflight": args.soft_preflight,
            "branch_result": args.branch_result, "preflight": args.preflight,
            "grid_result_manifest": args.grid_result_manifest,
            **({"all_points_audit": args.all_points_audit} if raw_audited is not None else {}),
        }.items()},
        "limitations": "one fixed Q and one mixed direction only; positive results do not establish a global or continuous conditional PES",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": output["status"], "curvatures": curvatures}, indent=2))


if __name__ == "__main__":
    main()
