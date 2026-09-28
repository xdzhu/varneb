"""Freeze an even cubic-mode polynomial and two genuinely new BTO holdouts.

The model trains on nine audited grid nodes only. The four previously measured
interior points are reported retrospectively, never included in the fit. The
two new predictions must be saved before their DFT inputs are submitted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


NEW_HOLDOUTS = (
    ((0.45, 0.25), ((0.3, 0.0), (0.3, 0.3), (0.6, 0.0), (0.6, 0.3)),
     (1/12, 5/12, 1/12, 5/12), "bilinear_structure_seed"),
    ((1.05, 0.25), ((0.9, 0.0), (0.9, 0.3), (1.2, 0.3)),
     (1/6, 1/3, 1/2), "triangle_structure_seed"),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def even_basis(q: np.ndarray) -> np.ndarray:
    q = np.atleast_2d(np.asarray(q, dtype=float))
    z, x = q[:, 0], q[:, 1]
    return np.column_stack((
        z*z + x*x,
        z**4 + x**4,
        z*z*x*x,
        z**6 + x**6,
        z**4*x*x + z*z*x**4,
    ))


def plan(measured: dict, retrospective: dict) -> dict:
    if (measured.get("status") != "ten_measured_nodes_one_independent_center_holdout_no_contour_certificate"
            or measured.get("n_measured_nodes") != 10
            or retrospective.get("status") != "interpolation_gate_failed"
            or len(retrospective.get("points", [])) != 3):
        raise ValueError("requires the audited nine-node grid and three failed linear holdouts")
    rows = measured["points"]
    by_q = {tuple(row["q"]): row for row in rows}
    train = [row for row in rows if row["role"] == "grid_node"]
    if len(by_q) != 10 or len(train) != 9:
        raise ValueError("training set must be nine unique measured grid nodes")
    q = np.array([row["q"] for row in train], dtype=float)
    y = np.array([row["energy_minus_c_eV_per_BTO"] for row in train], dtype=float)
    matrix = even_basis(q)
    if np.linalg.matrix_rank(matrix) != 5:
        raise ValueError("even sixth-order basis is rank deficient")
    coefficients = np.linalg.lstsq(matrix, y, rcond=None)[0]
    train_rms = float(np.sqrt(np.mean((matrix @ coefficients - y) ** 2)) * 1000)
    old_holdouts = [
        {"q": row["q"], "observed_energy_minus_C_eV_per_BTO": row["energy_minus_c_eV_per_BTO"],
         "role": "preexisting_center_not_used_for_fit"}
        for row in rows if row["role"] == "independent_cell_center_holdout"
    ]
    old_holdouts += [
        {"q": row["q_sqrt_amu_A"],
         "observed_energy_minus_C_eV_per_BTO": row["audited_energy_minus_C_eV_per_BTO"],
         "role": "prior_linear_interpolation_holdout_seen_before_model_choice"}
        for row in retrospective["points"]
    ]
    for row in old_holdouts:
        prediction = float(even_basis(np.asarray(row["q"]))[0] @ coefficients)
        row["model_prediction_eV_per_BTO"] = prediction
        row["observed_minus_predicted_meV_per_BTO"] = 1000 * (
            row["observed_energy_minus_C_eV_per_BTO"] - prediction)
    new = []
    for q_new, corners, weights, scheme in NEW_HOLDOUTS:
        if q_new in by_q or any(tuple(item["q"]) == q_new for item in old_holdouts):
            raise ValueError("prospective model holdout is already measured")
        w = np.array(weights, dtype=float)
        if (not np.isclose(w.sum(), 1, rtol=0, atol=1e-12)
                or not np.allclose(w @ np.array(corners), q_new, rtol=0, atol=1e-12)
                or any(corner not in by_q for corner in corners)):
            raise ValueError("prospective holdout is outside its measured interpolation cell")
        seed = np.tensordot(w, [by_q[corner]["coordinates_u_A_eta_voigt"]
                                for corner in corners], axes=(0, 0))
        new.append({
            "q_sqrt_amu_A": list(q_new),
            "model_prediction_energy_minus_C_eV_per_BTO": float(
                even_basis(np.asarray(q_new))[0] @ coefficients),
            "predicted_structure_u_A_eta_voigt": seed.tolist(),
            "structure_seed_interpolation": scheme,
            "structure_seed_corners_q": [list(corner) for corner in corners],
            "structure_seed_weights": list(weights),
        })
    return {
        "kind": "BTO_Qy0_even_mode_sixth_order_model_two_new_pre_DFT_holdouts",
        "status": "two_model_predictions_frozen_before_new_DFT",
        "model": "E-E_C = a2*(Qz^2+Qx^2) + a4*(Qz^4+Qx^4) + b4*Qz^2*Qx^2 + a6*(Qz^6+Qx^6) + b6*(Qz^4*Qx^2+Qz^2*Qx^4)",
        "coefficient_order": ["a2", "a4", "b4", "a6", "b6"],
        "coefficients_eV_per_BTO": coefficients.tolist(),
        "zero_intercept_at_cubic": True,
        "n_training_grid_nodes": len(train),
        "training_rms_meV_per_BTO": train_rms,
        "retroactively_evaluated_not_independent_for_model_selection": old_holdouts,
        "energy_absolute_error_gate_meV_per_BTO": 2.0,
        "full_structure_seed_error_gate_sqrt_amu_A": 0.10,
        "gradient_gate_eV_per_sqrt_amu_A": 0.003,
        "maximum_absolute_stress_gate_kbar": 2.0,
        "fixed_third_soft_mode": "Q_y=0",
        "new_holdouts": new,
        "limitations": "exploratory model chosen after seeing the first three linear-interpolation holdouts; only the two NEW points are prospective validation. No positive-curvature or global-PES claim follows from this fit.",
        "no_DFT_launched": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--measured", type=Path, required=True)
    parser.add_argument("--linear-holdout-audit", type=Path, required=True)
    args = parser.parse_args()
    result = plan(json.loads(args.measured.read_text(encoding="utf-8")),
                  json.loads(args.linear_holdout_audit.read_text(encoding="utf-8")))
    result["source_sha256"] = {"measured": sha256(args.measured),
                               "linear_holdout_audit": sha256(args.linear_holdout_audit),
                               "planner": sha256(Path(__file__))}
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
