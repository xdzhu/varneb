"""Read-only line-integral check of audited BTO fixed-Q soft-direction probes.

This quantifies, but cannot identify the origin of, an energy-versus-stress
gradient mismatch. It neither runs DFT nor certifies a conditional minimum.
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

from examples.bto_q1q2_reference import sha256
from vcneb.mode_surface import audit_directional_work_consistency


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--soft-preflight", type=Path, required=True)
    parser.add_argument("--soft-summary", type=Path, required=True)
    parser.add_argument("--soft-audit", type=Path, required=True)
    parser.add_argument("--branch-result", type=Path, required=True)
    parser.add_argument("--center-result", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite audit: {args.output}")
    preflight = json.loads(args.soft_preflight.read_text(encoding="utf-8"))
    summary = json.loads(args.soft_summary.read_text(encoding="utf-8"))
    audit = json.loads(args.soft_audit.read_text(encoding="utf-8"))
    branch = json.loads(args.branch_result.read_text(encoding="utf-8"))
    center = json.loads(args.center_result.read_text(encoding="utf-8"))
    direction = np.asarray(preflight["soft_direction_metric_unit_chart"], dtype=float)
    center_coordinates = np.asarray(preflight["center_coordinates_u_A_eta_voigt"], dtype=float)
    center_gradient = np.asarray(center["gradient"], dtype=float)
    if (preflight.get("status") != "four_soft_direction_geometries_passed_no_dft"
            or summary.get("status") != "four_static_probes_complete_pending_independent_SCF_audit"
            or audit.get("status") != "four_original_static_DFT_points_verified_soft_direction_curvature_screen_only"
            or audit.get("source_sha256", {}).get("summary") != sha256(args.soft_summary)
            or audit.get("source_sha256", {}).get("soft_preflight") != sha256(args.soft_preflight)
            or audit.get("source_sha256", {}).get("branch_result") != sha256(args.branch_result)
            or summary.get("center_result_sha256") != sha256(args.center_result)
            or not branch.get("final_evaluation_directory", "").endswith(
                center.get("coordinate_sha256", "")[:12])
            or branch.get("evaluator_contract_sha256") != center.get("contract_sha256")
            or audit.get("center_log_sha256") != center.get("validation", {}).get("log_sha256")
            or branch.get("q1_q2_sqrt_amu_A") != preflight.get("q1_q2_sqrt_amu_A")
            or direction.shape != center_coordinates.shape or center_gradient.shape != direction.shape
            or not all(np.all(np.isfinite(array)) for array in
                       (direction, center_coordinates, center_gradient))
            or not np.array_equal(np.asarray(center["coordinates"], dtype=float), center_coordinates)
            or len(audit.get("probes", [])) != 4):
        raise ValueError("soft-direction work inputs differ from the independently audited five-point line")
    points = {}
    for record in audit["probes"]:
        signed = int(record["sign"]) * float(record["step_sqrt_amu_A"])
        if signed in points or not np.isfinite(signed):
            raise ValueError("duplicate or nonfinite signed probe coordinate")
        points[signed] = record
    if len(points) != 4:
        raise ValueError("four distinct signed probes are required")
    abscissas = np.asarray([*sorted(points), 0.0], dtype=float)
    abscissas.sort()
    center_energy = float(audit["center_full_precision_log_energy_eV"])
    energies = np.asarray([
        center_energy if value == 0 else points[value]["full_precision_log_energy_eV"]
        for value in abscissas
    ], dtype=float)
    gradients = np.asarray([
        float(direction @ center_gradient) if value == 0 else
        points[value]["directional_gradient_eV_per_sqrt_amu_A"]
        for value in abscissas
    ], dtype=float)
    result = audit_directional_work_consistency(abscissas, energies, gradients)
    report = {
        "kind": "BTO_fixed_Q_soft_direction_energy_gradient_line_integrability_audit_not_PES",
        "status": "work_residual_quantified_cause_unresolved",
        "q1_q2_sqrt_amu_A": preflight["q1_q2_sqrt_amu_A"],
        "abscissas_sqrt_amu_A": result.abscissas.tolist(),
        "full_precision_energies_eV": energies.tolist(),
        "directional_gradients_eV_per_sqrt_amu_A": gradients.tolist(),
        "intervals": [
            {
                "from_sqrt_amu_A": float(abscissas[index]),
                "to_sqrt_amu_A": float(abscissas[index + 2]),
                "observed_delta_E_eV": float(result.energy_increments[part]),
                "simpson_integrated_gradient_eV": float(result.simpson_gradient_work[part]),
                "E_minus_integrated_gradient_meV": float(
                    1000 * result.energy_minus_gradient_work[part]),
            }
            for part, index in enumerate((0, 2))
        ],
        "maximum_absolute_work_residual_meV": float(
            1000 * result.maximum_absolute_work_residual),
        "source_sha256": {key: sha256(path) for key, path in {
            "soft_preflight": args.soft_preflight, "soft_summary": args.soft_summary,
            "soft_audit": args.soft_audit,
            "branch_result": args.branch_result, "center_result": args.center_result,
        }.items()},
        "limitations": "five DFT points and one mixed direction only; Simpson is exact through a cubic gradient but higher-order anharmonicity, electronic noise, or stress-adapter errors remain possible",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "intervals": report["intervals"]}, indent=2))


if __name__ == "__main__":
    main()
