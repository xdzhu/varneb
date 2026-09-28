"""Correct the atomic-mass convention in an immutable BTO holdout report.

The original three-point report used ASE tabulated O mass (15.999 amu)
instead of the force-constant source mass (15.9994 amu) for its *diagnostic*
full-coordinate interpolation error. It did not affect ABACUS or the frozen
prediction. This addendum leaves both earlier artifacts byte-for-byte intact.
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

from scripts.audit_bto_qy0_preregistered_holdouts import CASES, digest, fetch, gate_one


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--earlier-report", type=Path, required=True)
    parser.add_argument("--frozen-plan", type=Path, required=True)
    parser.add_argument("--force-constants", type=Path, required=True)
    parser.add_argument("--host", default="hf")
    parser.add_argument("--remote-root", required=True)
    args = parser.parse_args()
    old_bytes, plan_bytes = args.earlier_report.read_bytes(), args.frozen_plan.read_bytes()
    old, plan = json.loads(old_bytes), json.loads(plan_bytes)
    if (old.get("status") != "interpolation_gate_failed"
            or old.get("frozen_plan_sha256") != digest(plan_bytes)
            or len(old.get("points", [])) != len(CASES)
            or len(old.get("source_files", [])) != len(CASES)):
        raise ValueError("earlier report or frozen plan has changed")
    force_hash = digest(args.force_constants.read_bytes())
    with np.load(args.force_constants) as source:
        masses = np.asarray(source["masses_amu"], dtype=float)
    if masses.shape != (5,) or np.any(masses <= 0):
        raise ValueError("force-constant mass source must contain five positive masses")
    changes = []
    for i, (_, summary_path, audit_path, preflight_path) in enumerate(CASES):
        paths = (summary_path, audit_path, preflight_path)
        previous = old["source_files"][i]
        if tuple(previous[key] for key in ("summary", "audit", "preflight")) != paths:
            raise ValueError("source-file identity differs from the immutable report")
        blobs = {key: fetch(args.host, args.remote_root, path) for key, path in zip(
            ("summary", "audit", "preflight"), paths)}
        if any(digest(blobs[key]) != previous["sha256"][key] for key in blobs):
            raise ValueError("an already-audited remote source changed")
        docs = {key: json.loads(value) for key, value in blobs.items()}
        if docs["preflight"].get("source_sha256", {}).get("force_constants") != force_hash:
            raise ValueError("wrong force-constant mass source")
        corrected = gate_one(plan["holdouts"][i], docs["summary"], docs["audit"],
                             docs["preflight"], digest(blobs["summary"]),
                             digest(blobs["preflight"]), masses, plan)
        before = old["points"][i]
        if (abs(corrected["audited_minus_predicted_energy_meV_per_BTO"]
                - before["audited_minus_predicted_energy_meV_per_BTO"]) > 1e-12
                or corrected["checks"] != before["checks"]
                or corrected["all_gates_pass"] != before["all_gates_pass"]):
            raise ValueError("mass correction unexpectedly changes an original gate")
        changes.append({
            "q_sqrt_amu_A": corrected["q_sqrt_amu_A"],
            "original_ASE_mass_coordinate_error_sqrt_amu_A": before[
                "full_atom_plus_strain_metric_coordinate_error_sqrt_amu_A"],
            "correct_phonopy_mass_coordinate_error_sqrt_amu_A": corrected[
                "full_atom_plus_strain_metric_coordinate_error_sqrt_amu_A"],
            "coordinate_gate_unchanged_pass": corrected["checks"]["full_coordinate"],
        })
    print(json.dumps({
        "kind": "BTO_Qy0_holdout_metric_atomic_mass_source_addendum",
        "status": "diagnostic_coordinate_metric_corrected_original_gates_unchanged",
        "original_report_sha256": digest(old_bytes),
        "frozen_plan_sha256": digest(plan_bytes),
        "force_constants_sha256": force_hash,
        "atomic_masses_amu_in_phonopy_order": masses.tolist(),
        "script_sha256": digest(Path(__file__).read_bytes()),
        "points": changes,
        "limitations": "This corrects only offline structure-error bookkeeping; no ABACUS calculation, energy or frozen prediction changed.",
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
