"""Four frozen analytic controls, not HfO2 DFT or a material speed benchmark."""
from pathlib import Path
import argparse
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
sys.dont_write_bytecode = True
from vcneb import ModePlane, relax_orthogonal_at_q
from vcneb.feasible_relaxation import relax_orthogonal_in_ball

POLICY = dict(
    radius=0.5,
    gradient_target=0.03,
    initial_trust_radius=0.05,
    unique_call_cap=12,
    max_iterations=200,
)
CASES = (
    dict(
        name="interior_quadratic",
        start=[0.0, 0.0],
        target=[0.15, -0.1],
        seed="identity",
    ),
    dict(
        name="near_wall_misaligned_positive_seed",
        start=[0.4999, 0.0],
        target=[0.25, 0.15],
        seed="misaligned",
    ),
    dict(name="nonlinear_valley", start=[0.49, -0.08], seed="identity"),
    dict(
        name="exterior_quadratic_minimum",
        start=[0.0, 0.0],
        target=[0.8, 0.1],
        seed="identity",
    ),
)


def physical_value(case, x):
    z = x[2:]
    g = np.zeros(4)
    if case["name"] == "nonlinear_valley":
        u, v = z
        t = v - u * u
        energy = (u - 0.2) ** 2 + 10 * t * t
        g[2:] = [2 * (u - 0.2) - 40 * u * t, 20 * t]
    else:
        d = z - case["target"]
        energy = 0.5 * d @ d
        g[2:] = d
    return float(energy), g


def benchmark():
    p = ModePlane(
        np.zeros(4),
        np.eye(4),
        np.ones(4),
        np.eye(4)[:, :2],
        ("Q1", "Q2"),
        "analytic amplitude",
        "analytic/feasible-steps-not-material",
    )
    rows = []
    for case in CASES:
        for method in ("original_safeguarded_bfgs", "ball_feasible_bfgs"):
            cache = {}
            requests = 0
            guards = dict(checks=0, rejected=0)

            def guard(x):
                guards["checks"] += 1
                safe = bool(np.linalg.norm(x[2:]) <= POLICY["radius"])
                guards["rejected"] += int(not safe)
                return safe

            def evaluate(x):
                nonlocal requests
                requests += 1
                key = x.tobytes()
                if key not in cache:
                    if len(cache) >= POLICY["unique_call_cap"]:
                        raise RuntimeError(
                            "explicit analytic unique-evaluation cap exhausted"
                        )
                    if not np.linalg.norm(x[2:]) <= POLICY["radius"]:
                        raise ValueError(
                            "analytic oracle called outside unchanged domain"
                        )
                    energy, g = physical_value(case, x)
                    cache[key] = dict(x=x.copy(), energy=energy, gradient=g.copy())
                r = cache[key]
                return r["energy"], r["gradient"].copy()

            h = np.eye(4)
            if case["seed"] == "misaligned":
                h[2:, 2:] = np.linalg.inv(np.array([[1.0, 5.0], [5.0, 30.0]]))
            common = dict(
                start=np.r_[0.0, 0.0, case["start"]],
                trial_validator=guard,
                initial_hessian=h,
                initial_trust_radius=POLICY["initial_trust_radius"],
                gradient_tolerance=POLICY["gradient_target"],
                max_iterations=POLICY["max_iterations"],
            )
            result = None
            reason = None
            try:
                if method == "original_safeguarded_bfgs":
                    start = common.pop("start")
                    result = relax_orthogonal_at_q(
                        p,
                        [0, 0],
                        evaluate,
                        starts=[start],
                        include_frozen_start=False,
                        optimizer="safeguarded_bfgs",
                        **common
                    )
                else:
                    result = relax_orthogonal_in_ball(
                        p,
                        [0, 0],
                        evaluate,
                        orthogonal_radius_bound=POLICY["radius"],
                        **common
                    )
            except RuntimeError as exc:
                reason = str(exc)
                if reason not in (
                    "explicit analytic unique-evaluation cap exhausted",
                    "no starting point reached the orthogonal-gradient tolerance",
                    "no full released-gradient convergence within declared ball",
                ):
                    raise
            rows.append(
                dict(
                    case=case["name"],
                    method=method,
                    status="gradient_candidate"
                    if result is not None
                    else "no_candidate",
                    reason=reason,
                    callback_requests=requests,
                    unique_analytic_evaluations=len(cache),
                    guard_counts=guards,
                    candidate_coordinates=None
                    if result is None
                    else result.coordinates.tolist(),
                    candidate_released_gradient_norm=None
                    if result is None
                    else result.orthogonal_gradient_norm,
                    all_trial_gradient_norms=[
                        float(np.linalg.norm(r["gradient"][2:])) for r in cache.values()
                    ],
                    all_evaluated_radius_valid=all(
                        np.linalg.norm(r["x"][2:]) <= POLICY["radius"]
                        for r in cache.values()
                    ),
                )
            )
    return dict(
        status="explicit_analytic_control_comparison_not_material_validation",
        policy=POLICY,
        frozen_cases=list(CASES),
        rows=rows,
        new_DFT_or_scheduler_writes=0,
        HfO2_acceleration_independent_prediction_or_new_optimizer_theorem_claimed=False,
        original_production_source_domain_target_or_inputs_modified=False,
        physical_reserved_geometry_generated_or_read=False,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError("fresh analytic output required")
    result = benchmark()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=2)
        f.write("\n")
    print(
        json.dumps(
            [
                dict(
                    case=r["case"],
                    method=r["method"],
                    status=r["status"],
                    unique_calls=r["unique_analytic_evaluations"],
                    guard_rejections=r["guard_counts"]["rejected"],
                )
                for r in result["rows"]
            ]
        )
    )


if __name__ == "__main__":
    main()
