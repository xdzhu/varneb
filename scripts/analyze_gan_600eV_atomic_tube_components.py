"""Resolve GaN tube interpolation error into even and odd q responses.

Uses only audited 600-eV off-path statics and the original path forces.
Force-informed estimates here are post-hoc diagnostics, not an independent
contour-validation test or a new electronic-parameter choice.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(first_path: Path, refinement_path: Path, transport_path: Path) -> dict:
    first = json.loads(first_path.read_text(encoding="utf-8"))
    refined = json.loads(refinement_path.read_text(encoding="utf-8"))
    transport = json.loads(transport_path.read_text(encoding="utf-8"))
    if (first.get("status") != "GaN_600eV_central_atomic_transverse_tube_raw_audited"
            or refined.get("status")
            != "GaN_600eV_central_atomic_tube_targeted_refinement_raw_audited"
            or transport.get("status")
            != "GaN_600eV_atomic_transverse_transport_offline_comparison"
            or refined["source_sha256"].get("first_audit") != sha256(first_path)
            or first["source_sha256"].get("transport_comparison") != sha256(transport_path)
            or refined["source_sha256"].get("transport_comparison") != sha256(transport_path)
            or refined.get("anchor_q_A") != [-0.05, 0.05]
            or refined.get("predeclared_LOO_gate_meV_per_GaN") != 1.0):
        raise ValueError("600-eV raw audit and force provenance do not match")
    anchors = [int(index) for index in refined["anchors"]]
    arc = np.asarray(refined["arc_fraction_s"], dtype=float)
    slopes = np.asarray(transport["fixed_seed_projection"]["transverse_slopes_eV_per_A"],
                        dtype=float)
    overlaps = np.asarray(transport["fixed_seed_projection"]["neighbor_normal_overlaps"],
                          dtype=float)
    if (len(anchors) != 12 or anchors != sorted(set(anchors)) or len(arc) != 29
            or len(slopes) != 29 or len(overlaps) != 28
            or not np.all(np.diff(arc) > 0)
            or not np.array_equal(arc, np.asarray(first["arc_fraction_s"], dtype=float))):
        raise ValueError("central path data changed")

    signed: dict[int, dict[str, float]] = {index: {} for index in anchors}
    for case in [*first["cases"], *refined["new_cases"]]:
        if case.get("role", "interpolation_anchor") != "interpolation_anchor":
            continue
        index = int(case["image_index"])
        q = float(case["q_atom_A"])
        side = "plus" if q > 0 else "minus"
        if (index not in signed or not np.isclose(abs(q), 0.05, atol=1e-12)
                or side in signed[index]
                or not np.isclose(case["arc_fraction_s"], arc[index], atol=1e-12)):
            raise ValueError("signed point does not match central chart")
        signed[index][side] = float(case["delta_enthalpy_meV_per_GaN"])
    if any(set(pair) != {"minus", "plus"} for pair in signed.values()):
        raise ValueError("signed pair missing")

    q = 0.05
    rows = []
    for index in anchors:
        minus, plus = signed[index]["minus"], signed[index]["plus"]
        even = (plus + minus) / 2
        odd = (plus - minus) / 2
        # meV/GaN -> eV/two-GaN-cell: multiply by 2/1000.
        curvature = 4 * even / (1000 * q**2)
        finite_slope = odd / (500 * q)
        force_odd = 500 * q * slopes[index]
        rows.append({
            "image_index": index,
            "arc_fraction_s": float(arc[index]),
            "minus_excess_meV_per_GaN": minus,
            "plus_excess_meV_per_GaN": plus,
            "even_excess_meV_per_GaN": even,
            "odd_excess_meV_per_GaN": odd,
            "frozen_atomic_curvature_eV_per_A2": curvature,
            "finite_energy_transverse_slope_eV_per_A": finite_slope,
            "path_force_transverse_slope_eV_per_A": float(slopes[index]),
            "odd_from_path_force_meV_per_GaN": float(force_odd),
            "energy_force_odd_residual_meV_per_GaN": float(odd - force_odd),
            "preceding_normal_overlap": None if index == 0 else float(overlaps[index - 1]),
            "following_normal_overlap": None if index == 28 else float(overlaps[index]),
        })

    first_diagnostics = {int(item["image_index"]): item
                         for item in first["anchor_diagnostics"]}
    for row in rows:
        original = first_diagnostics.get(row["image_index"])
        if original is None:
            continue
        if (not np.isclose(row["frozen_atomic_curvature_eV_per_A2"],
                           original["frozen_atomic_transverse_curvature_eV_per_A2"],
                           atol=1e-10, rtol=0)
                or not np.isclose(row["finite_energy_transverse_slope_eV_per_A"],
                                  original["energy_finite_difference_slope_eV_per_A"],
                                  atol=1e-10, rtol=0)
                or not np.isclose(row["path_force_transverse_slope_eV_per_A"],
                                  original["analytic_atomic_transverse_slope_eV_per_A"],
                                  atol=1e-10, rtol=0)):
            raise ValueError("even/odd energy and force units disagree with raw audit")

    components = {row["image_index"]: row for row in rows}
    loo = []
    for position in range(1, len(anchors) - 1):
        left, middle, right = anchors[position - 1:position + 2]
        weight = (arc[middle] - arc[left]) / (arc[right] - arc[left])
        predicted_even = ((1 - weight) * components[left]["even_excess_meV_per_GaN"]
                          + weight * components[right]["even_excess_meV_per_GaN"])
        predicted_odd = ((1 - weight) * components[left]["odd_excess_meV_per_GaN"]
                         + weight * components[right]["odd_excess_meV_per_GaN"])
        actual_even = components[middle]["even_excess_meV_per_GaN"]
        actual_odd = components[middle]["odd_excess_meV_per_GaN"]
        force_odd = components[middle]["odd_from_path_force_meV_per_GaN"]
        loo.append({
            "image_index": middle,
            "even_interpolation_error_meV_per_GaN": float(actual_even - predicted_even),
            "odd_interpolation_error_meV_per_GaN": float(actual_odd - predicted_odd),
            "force_informed_minus_absolute_error_meV_per_GaN": float(abs(
                components[middle]["minus_excess_meV_per_GaN"]
                - (predicted_even - force_odd))),
            "force_informed_plus_absolute_error_meV_per_GaN": float(abs(
                components[middle]["plus_excess_meV_per_GaN"]
                - (predicted_even + force_odd))),
        })
    maximum_even = max(abs(row["even_interpolation_error_meV_per_GaN"]) for row in loo)
    maximum_odd = max(abs(row["odd_interpolation_error_meV_per_GaN"]) for row in loo)
    maximum_force_informed = max(
        max(row["force_informed_minus_absolute_error_meV_per_GaN"],
            row["force_informed_plus_absolute_error_meV_per_GaN"])
        for row in loo
    )
    return {
        "status": "GaN_600eV_atomic_tube_even_odd_force_diagnostic",
        "electronic_contract": "GaN VASP 600 eV original path and same-input signed statics",
        "q_atom_A": q,
        "anchor_components": rows,
        "leave_one_anchor_out_components": loo,
        "maximum_even_LOO_absolute_error_meV_per_GaN": float(maximum_even),
        "maximum_odd_LOO_absolute_error_meV_per_GaN": float(maximum_odd),
        "maximum_force_informed_signed_LOO_absolute_error_meV_per_GaN": float(
            maximum_force_informed),
        "maximum_energy_force_odd_residual_meV_per_GaN": float(max(
            abs(row["energy_force_odd_residual_meV_per_GaN"]) for row in rows)),
        "predeclared_signed_linear_LOO_gate_meV_per_GaN": 1.0,
        "contour_certified": False,
        "limitations": [
            "Even/odd and force-informed decompositions were chosen after the first LOO failure.",
            "An original-path projected force is not an independently sampled off-path enthalpy.",
            "q is an atomic-only frozen displacement about each variable-cell NEB image.",
            "The chart is not continuous through the B4/B1 endpoints.",
        ],
        "source_sha256": {
            "first_raw_audit": sha256(first_path),
            "refinement_raw_audit": sha256(refinement_path),
            "transport": sha256(transport_path),
            "analysis_script": sha256(Path(__file__)),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-audit", type=Path, required=True)
    parser.add_argument("--refinement-audit", type=Path, required=True)
    parser.add_argument("--transport", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = analyze(args.first_audit, args.refinement_audit, args.transport)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "max_even_LOO_meV_per_GaN": result["maximum_even_LOO_absolute_error_meV_per_GaN"],
        "max_odd_LOO_meV_per_GaN": result["maximum_odd_LOO_absolute_error_meV_per_GaN"],
        "max_force_informed_LOO_meV_per_GaN": result[
            "maximum_force_informed_signed_LOO_absolute_error_meV_per_GaN"],
        "contour_certified": result["contour_certified"],
    }))


if __name__ == "__main__":
    main()
