"""Bounded analytic fixed-control benchmark; no DFT or material forecast."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from vcneb.stationary_branch import stationary_quadratic_branch


def run_benchmark():
    rows = []
    errors = {"displacement": 0., "energy": 0., "internal_gradient": 0., "control_derivative": 0.}
    for index in (0, 1):
        for release in (False, True):
            for rotate in (False, True):
                h = np.array([[2. if index == 0 else -2., .4, .7, 0.],
                              [.4, 4., .3, 0.], [.7, .3, -3., 0.], [0., 0., 0., -6.]])
                g = np.array([.1, -.15, .04, .5])
                u = np.linalg.qr(np.random.default_rng(812).normal(size=(4, 4)))[0] if rotate else np.eye(4)
                h, g = u.T @ h @ u, u.T @ g
                q, c = u.T[:, :1], u.T[:, 2]
                r = u.T[:, 1:2] if release else np.empty((4, 0))
                allowed = u.T[:, [0, 1, 3]]
                internal = np.column_stack((q, r))
                branch = stationary_quadratic_branch(h, g, -10000., q, r, c,
                    admissible_internal_basis=allowed, expected_index=index,
                    stability_floor=.01, internal_curvature_floor=.01)
                for t in (-.12, 0., .15):
                    point = branch.evaluate(t)
                    solution = -np.linalg.solve(internal.T @ h @ internal, internal.T @ (g + h @ c * t))
                    direct = internal @ solution + c * t
                    expected_energy = float(g @ direct + .5 * direct @ h @ direct)
                    dt = 1e-5
                    derivative = (branch.evaluate(t+dt).energy_change - branch.evaluate(t-dt).energy_change)/(2*dt)
                    errors["displacement"] = max(errors["displacement"], float(np.max(np.abs(point.full_displacement-direct))))
                    errors["energy"] = max(errors["energy"], abs(point.energy_change-expected_energy))
                    errors["internal_gradient"] = max(errors["internal_gradient"], point.declared_internal_gradient_norm)
                    errors["control_derivative"] = max(errors["control_derivative"], abs(derivative-point.control_gradient))
                    rows.append({"index": index, "release": release, "rotated": rotate, "control": t,
                                 "energy_change": point.energy_change, "control_gradient": point.control_gradient,
                                 "control_curvature": branch.control_curvature,
                                 "unrepresented_internal_gradient_norm": point.unrepresented_internal_gradient_norm})
    if (max(errors[k] for k in ("displacement", "energy", "internal_gradient")) > 1e-12
            or errors["control_derivative"] > 1e-9):
        raise ValueError("analytic fixed-control benchmark failed; retain failure, do not adjust on material labels")
    root = Path(__file__).resolve().parents[1]
    sources = ["vcneb/stationary_branch.py", "vcneb/quadratic_reduction.py", "vcneb/relaxed_curvature.py",
               "scripts/check_stationary_branch.py"]
    return {"format_version": 1, "benchmark": "analytic_fixed_control_not_material_prediction",
            "groups": 8, "points": len(rows), "maximum_errors": errors, "rows": rows,
            "source_sha256": {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in sources},
            "numpy_version": np.__version__, "DFT_calls": 0, "calculator_calls": 0,
            "HfO2_curvature_or_independent_predictions_measured": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing an existing stationary-branch benchmark report")
    report = run_benchmark()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"groups": report["groups"], "points": report["points"],
                      "maximum_errors": report["maximum_errors"], "DFT_calls": 0}))


if __name__ == "__main__":
    main()
