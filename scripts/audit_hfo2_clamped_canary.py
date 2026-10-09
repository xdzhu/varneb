"""Audit the single registered PO+ four-step clamped BFGS canary.

A healthy step-limit canary is not an endpoint or G2 barrier convergence pass.
No source repair, symmetry restoration, DFT call or restart is performed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru
from scripts.analyze_hfo2_network_update import structure_audit
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.hfo2_switching_path_polarization import save
from scripts.relax_clamped_ase_endpoint import load_seed, stress_report


def validate_summary(summary, manifest, count):
    if (manifest["phase_label"] != "PO_plus" or manifest["strain"] != 0.
            or manifest["pressure_gpa"] != 0. or not manifest["allow_tilt"]
            or manifest["physical_contract_sha256"] != CONTRACT):
        raise ValueError("only registered zero-strain PO_plus tilt-released canary allowed")
    if (summary["steps_requested"] != 4 or type(summary["optimizer_steps"]) is not int
            or summary["optimizer_steps"] not in range(5)
            or summary["fmax_target_eV_per_A"] != .03 or summary["open_stress_target_kbar"] != 2.
            or summary["maxstep"] != .02 or summary["external_pressure_gpa"] != 0.
            or summary["strain"] != 0. or not summary["allow_tilt"]
            or summary["physical_contract_sha256"] != CONTRACT):
        raise ValueError("fixed four-step canary recipe changed")
    if count != summary["optimizer_steps"]+1 or not 1 <= count <= 5:
        raise ValueError("canary must have one fresh SCF per geometry, max5")
    converged = summary["converged"]
    if not isinstance(converged, bool) or summary["status"] != ("completed" if converged else "step_limit"):
        raise ValueError("step-cap versus convergence status inconsistent")
    if (not np.isfinite([summary["max_atomic_force_eV_per_A"], summary["open_traction_norm_kbar"]]).all()
            or min(summary["max_atomic_force_eV_per_A"], summary["open_traction_norm_kbar"]) < 0):
        raise ValueError("finite nonnegative physical convergence diagnostics required")
    physical = (summary["max_atomic_force_eV_per_A"] < .03
                and summary["open_traction_norm_kbar"] < 2.)
    if converged != physical:
        raise ValueError("physical atomic/open-traction convergence assertion inconsistent")


def replay_directories(calls, boundary, *, full_physical_bytes=True, first_index=0):
    """Replay only observed contiguous raw SCFs; never evaluate a calculator."""
    rows = []
    for i, directory in enumerate(calls, start=first_index):
        if directory.name != f"scf_{i:06d}":
            raise ValueError("fresh contiguous call sequence required")
        call = json.loads((directory/"call_audit.json").read_text())
        names = tuple(CONTRACT) if full_physical_bytes else ("INPUT", "KPT")
        if any(sha256(directory/name) != CONTRACT[name] for name in names):
            raise ValueError("canary physical bytes changed")
        if (any(call["input_sha256"][name] != value for name, value in CONTRACT.items())
                or sha256(directory/"STRU") != call["input_sha256"]["STRU"]
                or sha256(directory/"OUT.ABACUS/running_scf.log") != call["raw_log_sha256"]):
            raise ValueError("raw input/log record changed")
        atoms = read_fixed_hfo2_stru(directory/"STRU")
        boundary.validate_images([atoms])
        raw = audited_results(directory)
        if any(not np.allclose(raw[k], call["results"][k], atol=1e-12, rtol=0) for k in raw):
            raise ValueError("actual E/F/stress and call record differ")
        atoms.calc = SinglePointCalculator(atoms, **raw)
        elapsed = call["elapsed_seconds"]
        if not np.isfinite(elapsed) or elapsed <= 0:
            raise ValueError("finite positive actual call time required")
        rows.append({"call": i, "energy_eV_cell": raw["energy"],
                     "max_atomic_force_eV_A": float(np.linalg.norm(raw["forces"], axis=1).max()),
                     "substrate_plane_drift_A": float(np.abs(atoms.cell.array[:2]-boundary.reference_cell[:2]).max()),
                     "elapsed_seconds": elapsed, "raw_log_sha256": call["raw_log_sha256"],
                     "input_sha256": call["input_sha256"], "structure_audit": structure_audit(atoms),
                     **stress_report(atoms, boundary)})
    return rows


def audit(root, seed_manifest, *, full_physical_bytes=True):
    _, boundary, manifest = load_seed(seed_manifest)
    summary_path = root/"endpoint_relax_summary.json"
    summary = json.loads(summary_path.read_text())
    calls = sorted((root/"calculator/image_0000").glob("scf_*"))
    validate_summary(summary, manifest, len(calls))
    if summary["seed_manifest_sha256"] != sha256(seed_manifest):
        raise ValueError("summary seed changed")
    rows = replay_directories(calls, boundary, full_physical_bytes=full_physical_bytes)
    final = rows[-1]
    for key, value in (("potential_energy_eV", final["energy_eV_cell"]),
                       ("max_atomic_force_eV_per_A", final["max_atomic_force_eV_A"]),
                       ("open_traction_norm_kbar", final["open_traction_norm_kbar"])):
        if not np.isclose(summary[key], value, atol=2e-8, rtol=0):
            raise ValueError("terminal summary not reproduced by last raw geometry")
    return {"status": "bounded_canary_raw_audit_passed_not_G2_path_or_stability",
            "new_DFT_calls_for_analysis": 0, "new_SCF_calls": len(rows), "real_MPI_ranks": 32,
            "endpoint_converged": summary["converged"], "endpoint_status": summary["status"],
            "SCF_seconds": sum(r["elapsed_seconds"] for r in rows),
            "SCF_core_hours": sum(r["elapsed_seconds"] for r in rows)*32/3600,
            "raw_full_physical_bytes_checked_here": full_physical_bytes,
            "summary_sha256": sha256(summary_path), "seed_manifest_sha256": sha256(seed_manifest),
            "analysis_driver_sha256": sha256(Path(__file__)), "rows": rows,
            "automatically_restart_or_submit_matrix": False,
            "phase_stability_or_G2_barriers_certified": False,
            "limitations": "Four-step transport/boundary health test; phase labels at three tolerances are not a Hessian-stability certificate. Clamped reactions are informational, only open traction enters endpoint convergence."}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--seed-manifest", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--offline-omitted-basis", action="store_true")
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError("refusing prior canary audit")
    result = audit(a.root, a.seed_manifest, full_physical_bytes=not a.offline_omitted_basis)
    save(a.output, result)
    print(json.dumps({k: result[k] for k in ("status", "new_SCF_calls", "endpoint_converged", "SCF_core_hours")}))
