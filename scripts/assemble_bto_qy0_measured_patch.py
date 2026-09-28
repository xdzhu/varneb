"""Assemble only raw-audited BTO Qy=0 sheet nodes, with no interpolation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.bto_q1q2_reference import sha256  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output_json.exists() or args.output_csv.exists():
        raise FileExistsError("measured-patch outputs must both be new paths")
    sources = json.loads(args.sources.read_text(encoding="utf-8"))
    if (sources.get("kind") != "bto_qy0_restricted_sheet_measured_patch_sources_no_interpolation"
            or len(sources.get("new_cases", [])) != 5):
        raise ValueError("expected the five independently audited new Qy=0 nodes")
    reused_path = Path(sources["reused_report"])
    if sha256(reused_path) != sources["reused_report_sha256"]:
        raise ValueError("the five archived-branch reuse report changed")
    reused = json.loads(reused_path.read_text(encoding="utf-8"))
    if (reused.get("status") != "existing_Qy0_points_eligible_for_restricted_sheet_reuse"
            or reused.get("n_reused_DFT_points") != 5
            or [row.get("role") for row in reused["points"]].count(
                "independent_cell_center_holdout") != 1):
        raise ValueError("four grid nodes and one independent center holdout are required")
    rows = []
    for old in reused["points"]:
        rows.append({
            "q": old["q_parallel_q_transverse_sqrt_amu_A"],
            "role": old["role"],
            "energy_minus_c_eV_per_BTO": old["energy_minus_c_eV_per_BTO"],
            "restricted_gradient_eV_per_sqrt_amu_A": old[
                "restricted_15D_gradient_norm_eV_per_sqrt_amu_A"],
            "maximum_absolute_stress_kbar": old["maximum_absolute_stress_kbar"],
            "third_soft_y_amplitude_sqrt_amu_A": old["third_soft_y_amplitude_sqrt_amu_A"],
            "coordinates_u_A_eta_voigt": old["coordinates_u_A_eta_voigt"],
            "provenance_sha256": old["terminal_result_sha256"],
            "source_class": "archived_raw_audited_Qy0_branch",
        })
    new_contract = None
    for case in sources["new_cases"]:
        audit_path, summary_path = Path(case["audit"]), Path(case["summary"])
        if sha256(audit_path) != case["audit_sha256"]:
            raise ValueError(f"new audit hash changed for {case['q']}")
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        q = np.asarray(case["q"], dtype=float)
        if (audit.get("status") != "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"
                or summary.get("status") != "orthogonal_gradient_and_stress_converged_curvature_unchecked"
                or summary.get("kind") != "bto_symmetry_restricted_qy_zero_variable_cell_local_candidate_not_PES_or_barrier"
                or audit.get("source_sha256", {}).get("summary") != sha256(summary_path)
                or not np.allclose(audit.get("q_parallel_q_transverse_sqrt_amu_A"), q,
                                   rtol=0, atol=1e-12)
                or not np.allclose(summary.get("q_parallel_q_transverse_sqrt_amu_A"), q,
                                   rtol=0, atol=1e-12)
                or audit.get("third_soft_mode_restricted_at_zero") is not True
                or summary.get("stress_target_passed") is not True
                or summary.get("stress_target_kbar") != 2.0
                or abs(audit.get("final_third_soft_y_amplitude_sqrt_amu_A", float("inf"))) > 1e-8
                or abs(summary.get("q_y_sqrt_amu_A", float("inf"))) > 1e-8
                or summary.get("electronic_kpoints") != [4, 4, 4]
                or summary.get("phonon_supercell") != [1, 1, 1]
                or summary.get("orthogonal_gradient_norm_eV_per_sqrt_amu_A", float("inf")) > 0.003
                or audit.get("final_maximum_absolute_stress_kbar", float("inf")) > 2.0
                or abs(audit["final_energy_minus_c_eV_per_BTO"]
                       - summary["energy_minus_c_eV_per_BTO"]) > 1e-8):
            raise ValueError(f"new conditional point is not audited and eligible: {case['q']}")
        contract = summary["evaluator_contract_sha256"]
        if new_contract is None:
            new_contract = contract
        elif new_contract != contract:
            raise ValueError("new nodes used different calculator contracts")
        rows.append({
            "q": q.tolist(),
            "role": "grid_node",
            "energy_minus_c_eV_per_BTO": summary["energy_minus_c_eV_per_BTO"],
            "restricted_gradient_eV_per_sqrt_amu_A": summary[
                "orthogonal_gradient_norm_eV_per_sqrt_amu_A"],
            "maximum_absolute_stress_kbar": audit["final_maximum_absolute_stress_kbar"],
            "third_soft_y_amplitude_sqrt_amu_A": audit["final_third_soft_y_amplitude_sqrt_amu_A"],
            "coordinates_u_A_eta_voigt": summary["coordinates_u_A_eta_voigt"],
            "provenance_sha256": case["audit_sha256"],
            "source_class": "new_raw_audited_Qy0_pilot",
        })
    if (len(rows) != 10 or len({tuple(row["q"]) for row in rows}) != 10
            or any(len(row["coordinates_u_A_eta_voigt"]) != 21
                   or not np.all(np.isfinite(row["coordinates_u_A_eta_voigt"]))
                   or row["restricted_gradient_eV_per_sqrt_amu_A"] > 0.003
                   or row["maximum_absolute_stress_kbar"] > 2.0
                   for row in rows)):
        raise ValueError("measured-patch nodes are incomplete or duplicated")
    rows.sort(key=lambda row: (row["q"][0], row["q"][1]))
    output = {
        "kind": "bto_qy0_restricted_sheet_ten_raw_audited_measured_nodes_not_interpolated_PES",
        "status": "ten_measured_nodes_one_independent_center_holdout_no_contour_certificate",
        "n_measured_nodes": len(rows),
        "n_reused_archived_nodes": 5,
        "n_new_pilot_nodes": 5,
        "n_independent_center_holdouts": 1,
        "energy_unit": "eV/BTO relative to common cubic C",
        "coordinate_unit": "sqrt(amu) Angstrom",
        "fixed_mode": "Q_y=0",
        "open_coordinates": "15 orthogonal atomic-plus-symmetric-strain directions",
        "calculator_contract_sha256_new_nodes": new_contract,
        "source_manifest_sha256": sha256(args.sources),
        "reused_report_sha256": sha256(reused_path),
        "auditor_sha256": sha256(Path(__file__)),
        "points": rows,
        "limitations": "measured nodes only; no branch/curvature or whole-domain interpolation certificate",
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    fields = ["Qz_sqrt_amu_A", "Qx_sqrt_amu_A", "E_minus_C_eV_per_BTO",
              "restricted_gradient_eV_per_sqrt_amu_A", "maximum_absolute_stress_kbar",
              "Qy_sqrt_amu_A", "role", "source_class", "provenance_sha256"]
    with args.output_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        for row in rows:
            writer.writerow([*row["q"], row["energy_minus_c_eV_per_BTO"],
                             row["restricted_gradient_eV_per_sqrt_amu_A"],
                             row["maximum_absolute_stress_kbar"],
                             row["third_soft_y_amplitude_sqrt_amu_A"], row["role"],
                             row["source_class"], row["provenance_sha256"]])
    print(json.dumps({"status": output["status"], "n_measured_nodes": len(rows),
                      "n_independent_center_holdouts": 1}, indent=2))


if __name__ == "__main__":
    main()
