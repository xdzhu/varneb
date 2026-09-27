"""Compare offline GaN atomic-tube interpolants without changing DFT inputs.

This is a post-hoc diagnostic of already audited 600-eV points. Selecting an
interpolant using these same leave-one-out errors cannot certify a contour;
an independent holdout would be needed for that claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import Akima1DInterpolator, CubicSpline, PchipInterpolator


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(first_path: Path, refinement_path: Path) -> dict:
    first = json.loads(first_path.read_text(encoding="utf-8"))
    refinement = json.loads(refinement_path.read_text(encoding="utf-8"))
    if (first.get("status") != "GaN_600eV_central_atomic_transverse_tube_raw_audited"
            or refinement.get("status")
            != "GaN_600eV_central_atomic_tube_targeted_refinement_raw_audited"
            or refinement["source_sha256"].get("first_audit") != sha256(first_path)
            or refinement.get("anchor_q_A") != [-0.05, 0.05]
            or refinement.get("predeclared_LOO_gate_meV_per_GaN") != 1.0
            or refinement.get("path_interpolation_gate_pass") is not False
            or len(first.get("cases", [])) != 16
            or len(refinement.get("new_cases", [])) != 12):
        raise ValueError("audited, same-600-eV comparison inputs are inconsistent")

    anchors = [int(index) for index in refinement["anchors"]]
    arc = np.asarray(refinement["arc_fraction_s"], dtype=float)
    if (anchors != sorted(set(anchors)) or len(anchors) != 12
            or len(arc) != 29 or not np.all(np.diff(arc) > 0)
            or not np.array_equal(arc, np.asarray(first["arc_fraction_s"], dtype=float))):
        raise ValueError("path abscissa or anchor map changed")
    values: dict[int, dict[str, float]] = {index: {} for index in anchors}
    for case in [*first["cases"], *refinement["new_cases"]]:
        if case.get("role", "interpolation_anchor") != "interpolation_anchor":
            continue
        index = int(case["image_index"])
        q = float(case["q_atom_A"])
        side = "plus" if q > 0 else "minus"
        value = float(case["delta_enthalpy_meV_per_GaN"])
        if (index not in values or not np.isclose(abs(q), 0.05, atol=1e-12)
                or side in values[index] or not np.isfinite(value)
                or not np.isclose(float(case["arc_fraction_s"]), arc[index], atol=1e-12)):
            raise ValueError("duplicate or nonmatching signed off-path point")
        values[index][side] = value
    if any(set(sides) != {"minus", "plus"} for sides in values.values()):
        raise ValueError("signed anchor pairs are incomplete")

    x = arc[anchors]
    methods = ("linear", "pchip", "akima", "natural_cubic")
    errors: dict[str, list[dict]] = {method: [] for method in methods}
    for holdout_position in range(1, len(anchors) - 1):
        image_index = anchors[holdout_position]
        keep = np.arange(len(anchors)) != holdout_position
        x_train = x[keep]
        x_test = float(x[holdout_position])
        for side in ("minus", "plus"):
            y = np.array([values[index][side] for index in anchors], dtype=float)
            y_train = y[keep]
            predictions = {
                "linear": float(np.interp(x_test, x_train, y_train)),
                "pchip": float(PchipInterpolator(x_train, y_train)(x_test)),
                "akima": float(Akima1DInterpolator(x_train, y_train)(x_test)),
                "natural_cubic": float(CubicSpline(x_train, y_train, bc_type="natural")(
                    x_test)),
            }
            for method, prediction in predictions.items():
                errors[method].append({
                    "image_index": image_index,
                    "side": side,
                    "predicted_excess_meV_per_GaN": prediction,
                    "actual_excess_meV_per_GaN": float(y[holdout_position]),
                    "absolute_error_meV_per_GaN": abs(prediction - y[holdout_position]),
                })
    summary = {}
    for method, records in errors.items():
        magnitudes = np.asarray([record["absolute_error_meV_per_GaN"]
                                 for record in records])
        worst = max(records, key=lambda record: record["absolute_error_meV_per_GaN"])
        summary[method] = {
            "n_holdouts": len(records),
            "maximum_absolute_error_meV_per_GaN": float(np.max(magnitudes)),
            "median_absolute_error_meV_per_GaN": float(np.median(magnitudes)),
            "rms_error_meV_per_GaN": float(np.sqrt(np.mean(magnitudes**2))),
            "worst_holdout": worst,
        }
    if not np.isclose(
            summary["linear"]["maximum_absolute_error_meV_per_GaN"],
            refinement["maximum_leave_one_out_error_meV_per_GaN"],
            atol=1e-9, rtol=0):
        raise ValueError("offline linear baseline does not reproduce the raw audit")
    return {
        "status": "GaN_600eV_atomic_tube_posthoc_interpolant_comparison",
        "pressure_GPa": 45.7,
        "electronic_contract": "Original 600-eV GaN VASP production inputs; no new DFT",
        "central_segment_images": [5, 22],
        "anchors": anchors,
        "q_atom_A": [-0.05, 0.05],
        "predeclared_linear_LOO_gate_meV_per_GaN": 1.0,
        "predeclared_linear_gate_pass": False,
        "method_summaries": summary,
        "individual_holdout_errors": errors,
        "contour_certified_by_this_comparison": False,
        "reason_not_certified": (
            "Interpolation families are compared after the existing holdout results were "
            "seen. The same points cannot both select a model and independently validate "
            "that selected model; the original predeclared linear gate remains failed."
        ),
        "source_sha256": {
            "first_raw_audit": sha256(first_path),
            "refinement_raw_audit": sha256(refinement_path),
            "comparison_script": sha256(Path(__file__)),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-audit", type=Path, required=True)
    parser.add_argument("--refinement-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = compare(args.first_audit, args.refinement_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["method_summaries"], indent=2))


if __name__ == "__main__":
    main()
