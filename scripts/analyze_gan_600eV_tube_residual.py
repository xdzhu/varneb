"""Diagnose why the measured GaN path tube fails the contour interpolation gate.

This report identifies correlations only. A force-residual sign change is not
by itself proof of a numerical branch switch or an incorrect VCNEB barrier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(first_path: Path, refined_path: Path, transport_path: Path) -> dict:
    first = json.loads(first_path.read_text(encoding="utf-8"))
    refined = json.loads(refined_path.read_text(encoding="utf-8"))
    transport = json.loads(transport_path.read_text(encoding="utf-8"))
    if (first.get("status")
            != "GaN_600eV_path_adapted_frozen_transverse_tube_raw_audited"
            or refined.get("status") != "GaN_600eV_path_tube_refined_raw_audited"
            or refined.get("smooth_contour_gate_pass") is not False
            or refined["source_sha256"].get("first_audit") != sha256(first_path)
            or refined["source_sha256"].get("transport_comparison") != sha256(transport_path)
            or transport.get("status")
            != "GaN_600eV_path_normal_transport_offline_comparison"):
        raise ValueError("600-eV tube audit sources changed")
    amplitude = float(refined["transverse_amplitude_A"])
    if not np.isclose(amplitude, 0.015) or first["anchors"] == refined["anchors"]:
        raise ValueError("unexpected path-tube refinement")
    by_anchor: dict[int, dict[str, float]] = {index: {} for index in refined["anchors"]}
    for record in [*first["cases"], *refined["new_cases"]]:
        index = int(record["image_index"])
        side = "plus" if record["q_perp_A"] > 0 else "minus"
        if side in by_anchor[index]:
            raise ValueError("duplicate measured signed point")
        by_anchor[index][side] = float(record["delta_enthalpy_meV_per_GaN"])
    analytic_slopes = transport["fixed_seed_projection"]["transverse_slopes_eV_per_A"]
    if len(analytic_slopes) != 29:
        raise ValueError("all-image analytic transverse slopes missing")
    rows = []
    for index in refined["anchors"]:
        values = by_anchor[index]
        if set(values) != {"minus", "plus"}:
            raise ValueError("signed pair missing")
        # The four-atom cell contains two GaN formula units: 1 meV/GaN =
        # 0.002 eV/cell. q=0 is the audited original path by construction.
        slope = (values["plus"] - values["minus"]) * 0.002 / (2 * amplitude)
        curvature = (values["plus"] + values["minus"]) * 0.002 / amplitude**2
        q_min = -slope / curvature if curvature > 0 else None
        rows.append({
            "image_index": index,
            "arc_fraction_s": refined["arc_fraction_s"][index],
            "delta_H_minus_meV_per_GaN": values["minus"],
            "delta_H_plus_meV_per_GaN": values["plus"],
            "energy_finite_difference_transverse_slope_eV_per_A": float(slope),
            "analytic_force_stress_transverse_slope_eV_per_A": float(analytic_slopes[index]),
            "frozen_transverse_curvature_eV_per_A2": float(curvature),
            "quadratic_local_minimum_q_perp_A": None if q_min is None else float(q_min),
            "quadratic_minimum_outside_measured_width": bool(
                q_min is not None and abs(q_min) > amplitude
            ),
        })
    errors = refined["leave_one_anchor_out_excess_enthalpy_errors"]
    worst = sorted(errors, key=lambda row: row["absolute_error_meV_per_GaN"], reverse=True)
    return {
        "status": "GaN_600eV_narrow_tube_contour_gate_failed_with_localized_transverse_slope_turns",
        "first_max_LOO_meV_per_GaN": first["maximum_leave_one_out_error_meV_per_GaN"],
        "refined_max_LOO_meV_per_GaN": refined["maximum_leave_one_out_error_meV_per_GaN"],
        "predeclared_gate_meV_per_GaN": refined["predeclared_interpolation_gate_meV_per_GaN"],
        "worst_refined_holdouts": worst[:10],
        "anchor_diagnostics": rows,
        "analytic_sign_turns": [
            {"left_image": i, "right_image": i + 1,
             "left_slope_eV_per_A": analytic_slopes[i],
             "right_slope_eV_per_A": analytic_slopes[i + 1]}
            for i in range(28) if analytic_slopes[i] * analytic_slopes[i + 1] < 0
        ],
        "source_sha256": {
            "first_audit": sha256(first_path),
            "refined_audit": sha256(refined_path),
            "transport_comparison": sha256(transport_path),
            "analyzer": sha256(Path(__file__)),
        },
        "interpretation_limit": (
            "The sharp off-path response correlates with local transverse slope sign changes; "
            "the existing 0.10 eV/A VCNEB acceptance does not certify a sub-meV smooth tube. "
            "Do not tune ENCUT or relax the predeclared contour gate. A distinct targeted "
            "MEP/conditional-coordinate study is needed before a smooth GaN full-path contour."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("first-audit", "refined-audit", "transport", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = analyze(args.first_audit, args.refined_audit, args.transport)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "before": report["first_max_LOO_meV_per_GaN"],
        "after": report["refined_max_LOO_meV_per_GaN"],
    }))


if __name__ == "__main__":
    main()
