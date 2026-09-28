"""Pre-register three unsampled BTO Qy=0 interpolation holdouts.

Uses only independently audited measured nodes. No calculator is called.
The upper T-side holdout is inside a measured triangle; lower two are cell
centers. Thresholds are declared here before any of these three DFT runs.
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

from examples.bto_q1q2_reference import sha256  # noqa: E402


HOLDOUTS = (
    ((0.15, 0.15), ((0.0, 0.0), (0.0, 0.3), (0.3, 0.0), (0.3, 0.3)),
     (0.25, 0.25, 0.25, 0.25), "bilinear_center"),
    ((0.45, 0.15), ((0.3, 0.0), (0.3, 0.3), (0.6, 0.0), (0.6, 0.3)),
     (0.25, 0.25, 0.25, 0.25), "bilinear_center"),
    ((1.05, 0.20), ((0.9, 0.0), (0.9, 0.3), (1.2, 0.3)),
     (1.0 / 3.0, 1.0 / 6.0, 0.5), "triangle_barycentric"),
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--measured", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite holdout plan: {args.output}")
    measured = json.loads(args.measured.read_text(encoding="utf-8"))
    if (measured.get("status") != "ten_measured_nodes_one_independent_center_holdout_no_contour_certificate"
            or measured.get("n_measured_nodes") != 10):
        raise ValueError("only the archived ten-node measured patch can seed this plan")
    by_q = {tuple(row["q"]): row for row in measured["points"]}
    if len(by_q) != 10:
        raise ValueError("measured coordinates are not unique")
    predictions = []
    for q, corners, weights, scheme in HOLDOUTS:
        if q in by_q or any(key not in by_q for key in corners):
            raise ValueError(f"holdout is already measured or lacks source corners: {q}")
        w = np.asarray(weights, dtype=float)
        corner_q = np.asarray(corners, dtype=float)
        if (not np.isclose(np.sum(w), 1.0, atol=1e-12, rtol=0)
                or not np.allclose(w @ corner_q, q, atol=1e-12, rtol=0)):
            raise ValueError(f"interpolation weights do not reach holdout Q: {q}")
        energy = sum(float(weight) * by_q[corner]["energy_minus_c_eV_per_BTO"]
                     for corner, weight in zip(corners, weights))
        coordinates = np.tensordot(w, np.asarray([
            by_q[corner]["coordinates_u_A_eta_voigt"] for corner in corners
        ], dtype=float), axes=(0, 0))
        predictions.append({
            "q_parallel_q_transverse_sqrt_amu_A": list(q),
            "interpolation": scheme,
            "source_corners_q": [list(corner) for corner in corners],
            "source_corner_weights": list(weights),
            "predicted_energy_minus_c_eV_per_BTO": float(energy),
            "predicted_coordinates_u_A_eta_voigt": coordinates.tolist(),
        })
    output = {
        "kind": "predeclared_three_BTO_Qy0_restricted_sheet_interpolation_holdouts_no_DFT",
        "status": "three_unmeasured_holdouts_predicted_before_DFT",
        "source_measured_patch_sha256": sha256(args.measured),
        "energy_absolute_error_gate_meV_per_BTO": 2.0,
        "full_atom_plus_strain_metric_coordinate_error_gate_sqrt_amu_A": 0.10,
        "gradient_gate_eV_per_sqrt_amu_A": 0.003,
        "maximum_absolute_stress_gate_kbar": 2.0,
        "fixed_third_soft_mode": "Q_y=0",
        "n_holdouts": len(predictions),
        "holdouts": predictions,
        "interpretation": "local error screens only; passing does not certify global smoothness, positive curvature or branch continuity",
        "no_DFT_launched": True,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": output["status"],
                      "holdouts": [{"q": row["q_parallel_q_transverse_sqrt_amu_A"],
                                    "predicted_energy": row["predicted_energy_minus_c_eV_per_BTO"]}
                                   for row in predictions]}, indent=2))


if __name__ == "__main__":
    main()
