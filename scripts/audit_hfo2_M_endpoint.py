"""Accept the new M endpoint only after auditing every retained real SCF."""

import argparse
import json
from pathlib import Path

from ase.io import read
import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.prepare_hfo2_reference_variants import symmetry_report
from vcneb import endpoint_structure_record


def audit(root, output):
    summary = json.loads((root / "endpoint_relax_summary.json").read_text(encoding="utf-8"))
    if summary["status"] != "completed" or not summary["converged"] or summary["external_pressure_gpa"] != 0:
        raise ValueError("endpoint is not a completed zero-pressure optimization")
    directories = sorted((root / "calculator/image_0000").glob("scf_*"))
    if len(directories) != summary["optimizer_steps"] + 1:
        raise ValueError("one retained raw SCF per endpoint BFGS geometry required")
    calls = []
    for source in directories:
        if any(sha256(source / n) != h for n, h in CONTRACT.items()):
            raise ValueError("an endpoint SCF changed the fixed electronic contract")
        raw = audited_results(source)
        calls.append({"source": str(source), "energy_eV_cell": float(raw["energy"]),
                      "max_atomic_force_eV_A": float(np.linalg.norm(raw["forces"], axis=1).max()),
                      "max_abs_stress_kbar": float(np.abs(raw["stress"]).max() * 1602.176634),
                      "input_sha256": {n: sha256(source / n) for n in (*CONTRACT, "STRU")},
                      "raw_log_sha256": sha256(source / "OUT.ABACUS/running_scf.log")})
    final = calls[-1]
    geometry = read(root / "CONTCAR", format="vasp")
    raw_geometry = read(directories[-1] / "STRU", format="abacus")
    if not same_ordered_geometry(geometry, raw_geometry):
        raise ValueError("final CONTCAR/raw SCF ordered geometry mismatch")
    if (abs(final["energy_eV_cell"] - summary["potential_energy_eV"]) > 1e-8
            or final["max_atomic_force_eV_A"] >= .03 or final["max_abs_stress_kbar"] >= 2):
        raise ValueError("final raw endpoint does not meet force/stress/energy gates")
    record = {"status": "raw_endpoint_gates_passed", "job_id": "28251302",
              "n_real_audited_SCFs": len(calls), "optimizer_steps": summary["optimizer_steps"],
              "energy_eV_cell": final["energy_eV_cell"],
              "max_atomic_force_eV_A": final["max_atomic_force_eV_A"],
              "max_abs_stress_kbar": final["max_abs_stress_kbar"],
              "endpoint": endpoint_structure_record(geometry), "symmetry_audit": symmetry_report(geometry),
              "final_static_source": str(directories[-1]), "calls": calls,
              "limitations": "force/stress acceptance and phase-identity tolerance sweep, not a full phonon stability certificate"}
    output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.root, args.output)
    print(json.dumps({k: report[k] for k in ("status", "n_real_audited_SCFs", "energy_eV_cell",
                     "max_atomic_force_eV_A", "max_abs_stress_kbar", "symmetry_audit")}))


if __name__ == "__main__":
    main()
