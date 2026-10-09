"""Replay the registered one-endpoint twenty-step segment, zero new DFT."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.audit_hfo2_clamped_canary import replay_directories
from scripts.audit_hfo2_static_replica import sha256
from scripts.hfo2_switching_path_polarization import save
from scripts.prepare_hfo2_clamped_continuation import CANARY_HF, ORIGINAL_CANARY_AUDIT_SHA256
from scripts.relax_clamped_ase_endpoint import load_seed


def validate_summary(summary, manifest, count):
    parent = manifest.get("continuation", {})
    if (manifest["phase_label"] != "PO_plus" or manifest["strain"] != 0.
            or manifest["pressure_gpa"] != 0. or not manifest["allow_tilt"]
            or manifest["physical_contract_sha256"] != CONTRACT
            or parent.get("parent_job") != "28446324" or parent.get("parent_BFGS_steps") != 4
            or parent.get("parent_audit_sha256") != ORIGINAL_CANARY_AUDIT_SHA256
            or parent.get("restore_BFGS_Hessian") is not False):
        raise ValueError("only the registered audited partial PO+ canary continuation is allowed")
    if (summary["steps_requested"] != 20 or type(summary["optimizer_steps"]) is not int
            or summary["optimizer_steps"] not in range(1, 21) or count != summary["optimizer_steps"]
            or not 1 <= count <= 20 or summary["fmax_target_eV_per_A"] != .03
            or summary["open_stress_target_kbar"] != 2. or summary["maxstep"] != .02
            or summary["external_pressure_gpa"] != 0. or summary["strain"] != 0.
            or not summary["allow_tilt"] or summary["physical_contract_sha256"] != CONTRACT):
        raise ValueError("fixed20steps/at most20new SCFs with one reused seed required")
    f, t = summary["max_atomic_force_eV_per_A"], summary["open_traction_norm_kbar"]
    if not np.isfinite([f, t]).all() or min(f, t) < 0:
        raise ValueError("finite physical diagnostics required")
    converged = summary["converged"]
    if (type(converged) is not bool or converged != (f < .03 and t < 2.)
            or summary["status"] != ("completed" if converged else "step_limit")):
        raise ValueError("physical endpoint convergence versus segment cap assertion inconsistent")


def audit(root, seed_manifest, parent_endpoint, *, full_physical_bytes=True):
    atoms, boundary, manifest = load_seed(seed_manifest)
    summary_path = root/"endpoint_relax_summary.json"
    summary = json.loads(summary_path.read_text())
    calls = sorted((root/"calculator/image_0000").glob("scf_*"))
    validate_summary(summary, manifest, len(calls))
    if summary["seed_manifest_sha256"] != sha256(seed_manifest):
        raise ValueError("continuation seed manifest changed")
    cache = json.loads((root/"calculator/image_0000/seed_cache_audit.json").read_text())
    original = parent_endpoint/"calculator/image_0000/scf_000004"
    if (sha256(parent_endpoint/"canary_audit.json") != ORIGINAL_CANARY_AUDIT_SHA256
            or cache["policy"] != "identical_ordered_clamped_seed_only_no_optimizer_history"
            or cache["new_DFT_calls"] != 0 or cache["raw_source"] != CANARY_HF+"/endpoint/calculator/image_0000/scf_000004"
            or cache["raw_log_sha256"] != manifest["continuation"]["parent_log_sha256"]
            or cache["input_sha256"]["STRU"] != manifest["source_sha256"]
            or not same_ordered_geometry(atoms, read_fixed_hfo2_stru(original/"STRU"))):
        raise ValueError("exact terminal clamped seed reuse provenance failed")
    seed_row = replay_directories([original], boundary, full_physical_bytes=full_physical_bytes, first_index=4)[0]
    if cache["input_sha256"] != seed_row["input_sha256"] or cache["raw_log_sha256"] != seed_row["raw_log_sha256"]:
        raise ValueError("cache and original raw seed disagree")
    rows = replay_directories(calls, boundary, full_physical_bytes=full_physical_bytes)
    final = rows[-1]
    for k, v in (("potential_energy_eV", final["energy_eV_cell"]),
                 ("max_atomic_force_eV_per_A", final["max_atomic_force_eV_A"]),
                 ("open_traction_norm_kbar", final["open_traction_norm_kbar"])):
        if not np.isclose(summary[k], v, atol=2e-8, rtol=0):
            raise ValueError("terminal continuation summary not reproduced by raw SCF")
    elapsed = sum(r["elapsed_seconds"] for r in rows)
    return {"status": "audited_bounded_PO_continuation_not_stability_or_G2_network",
        "new_DFT_calls_for_analysis": 0, "new_SCF_calls": len(rows), "reused_seed_SCFs": 1,
        "real_MPI_ranks": 32, "endpoint_converged": summary["converged"],
        "endpoint_status": summary["status"], "fresh_BFGS_steps": summary["optimizer_steps"],
        "parent_BFGS_steps": 4, "SCF_seconds": elapsed, "SCF_core_hours": elapsed*32/3600,
        "raw_full_physical_bytes_checked_here": full_physical_bytes,
        "summary_sha256": sha256(summary_path), "seed_manifest_sha256": sha256(seed_manifest),
        "analysis_driver_sha256": sha256(Path(__file__)), "reused_seed": seed_row, "rows": rows,
        "terminal_phase_symbols": [s["symbol"] for s in final["structure_audit"]["symmetry_sweep"]],
        "phase_stability_or_variants_or_G2_barriers_certified": False,
        "automatically_restart_or_submit_matrix": False,
        "limitations": "Physical endpoint gate and three-tolerance phase labels do not certify a Hessian minimum, electronic polarization or G2 channel barriers. Fresh BFGS Hessian; initial SCF reuse cost excluded from new segment cost."}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--seed-manifest", type=Path, required=True)
    p.add_argument("--parent-endpoint", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--offline-omitted-basis", action="store_true")
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError("refusing prior continuation audit")
    result = audit(a.root, a.seed_manifest, a.parent_endpoint,
                   full_physical_bytes=not a.offline_omitted_basis)
    save(a.output, result)
    print(json.dumps({k: result[k] for k in ("status", "new_SCF_calls", "endpoint_converged", "SCF_core_hours")}))
