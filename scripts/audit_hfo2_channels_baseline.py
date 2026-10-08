"""Normalize historical HfO2 summaries without merging different protocols."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


CASES = (
    ("ordinary_n7", "vcneb_n7_fire_distributed_job27678218", "27678218", "ordinary"),
    ("ordinary_n9", "vcneb_n9_fire_distributed_cmp_job27678407", "27678407", "ordinary"),
    ("historical_ci_n7", "vcneb_n7_ci_refine_job27678507", "27678507", "historical_CI"),
    ("guided_n7", "vcneb_n7_mode_guided_loose_job27687189", "27687189", "ordinary_guided_seed"),
)


def normalize_summary(summary: dict) -> dict:
    n_atoms = summary["path_diagnostics"]["n_atoms"]
    if n_atoms != 12:
        raise ValueError("this historical registry is restricted to 12-atom Hf4O8")
    energies = np.asarray(summary["image_enthalpies_eV"], dtype=float)
    fmax = float(summary["final_max_generalized_force_eV_per_A"])
    target = float(summary["fmax_target_eV_per_A"])
    if (energies.shape != (summary["n_images"],) or len(energies) < 3
            or not np.isfinite(energies).all() or not np.isfinite(fmax)
            or not np.isfinite(target) or fmax < 0 or target <= 0):
        raise ValueError("incomplete/nonfinite chain or invalid convergence metadata")
    parameters = summary["calculator_parameters"]
    if (parameters.get("ecutwfc") != 100 or parameters.get("kpts") != [2, 2, 2]
            or parameters.get("scf_thr") != 1e-8):
        raise ValueError("historical HfO2 physical contract differs")
    n_fu = 4
    peak = float(np.max(energies))
    forward, reverse = peak - energies[0], peak - energies[-1]
    delta = energies[-1] - energies[0]
    if (abs(forward - summary["barrier_enthalpy_eV"]) > 1e-7
            or abs(delta - summary["reaction_enthalpy_eV"]) > 1e-7):
        raise ValueError("stored barrier/reaction energy disagrees with full chain")
    return {
        "n_atoms": n_atoms, "formula_units": n_fu,
        "n_images_total": len(energies), "n_images_interior": len(energies) - 2,
        "forward_meV_fu": float(forward * 1000 / n_fu),
        "reverse_meV_fu": float(reverse * 1000 / n_fu),
        "reaction_meV_fu": float(delta * 1000 / n_fu),
        "barrier_identity_error_eV_cell": float(forward - reverse - delta),
        "endpoint_enthalpies_eV_cell": [float(energies[0]), float(energies[-1])],
        "relative_enthalpies_meV_fu": ((energies - energies[0]) * 1000 / n_fu).tolist(),
        "summary_status": summary["status"],
        "final_generalized_force_eV_A": fmax, "original_target_eV_A": target,
        "meets_original_force_target": fmax <= target,
        "meets_current_ordinary_0p10_force_target": fmax <= 0.10,
        "pressure_eV_A3": summary["path_diagnostics"]["pressure_eV_per_A3"],
        "peak_is_stationary_certified": False,
        "mode_variant_identity": "not_yet_audited",
        "small_barrier_resolution": "not_yet_measured",
        "physical_parameter_record": parameters,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("outputs/hfo2_t_to_po_pbe100_dzp10au"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    results = []
    for label, directory, job_id, protocol in CASES:
        path = args.root / directory / "vcneb_summary.json"
        report = normalize_summary(json.loads(path.read_text(encoding="utf-8")))
        report.update(label=label, historical_job_id=job_id, protocol=protocol,
                      source_summary=str(path), source_summary_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        results.append(report)
    reference = results[0]["endpoint_enthalpies_eV_cell"]
    endpoints_match = all(np.allclose(row["endpoint_enthalpies_eV_cell"], reference, atol=1e-7, rtol=0)
                          for row in results)
    report = {
        "schema_version": 1, "status": "historical_summary_audit_only",
        "common_endpoint_energies_match": endpoints_match,
        "endpoint_geometry_and_variant_match": "requires_raw_structure_audit",
        "interpretation": "protocols are separate; matching endpoint energies does not certify common variants or saddles",
        "records": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"common_endpoint_energies_match": endpoints_match,
                      "records": [{key: row[key] for key in ("label", "forward_meV_fu", "reverse_meV_fu")}
                                  for row in results]}, indent=2))


if __name__ == "__main__":
    main()
