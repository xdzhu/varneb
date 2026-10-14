"""The frozen analytic cohort is not a material speedup experiment."""
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "benchmarks/hfo2_channels/20261008/feasible_steps_E169_20261015"
spec = importlib.util.spec_from_file_location(
    "feasible_steps_case", CASE / "benchmark.py"
)
b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b)


def test_all_four_frozen_cases_and_both_arms_remain_in_the_saved_result():
    result = b.benchmark()
    saved = json.loads((CASE / "analysis.json").read_bytes())
    assert result == saved
    assert result["policy"] == dict(
        radius=0.5,
        gradient_target=0.03,
        initial_trust_radius=0.05,
        unique_call_cap=12,
        max_iterations=200,
    )
    assert len(result["frozen_cases"]) == 4 and len(result["rows"]) == 8
    assert not result[
        "HfO2_acceleration_independent_prediction_or_new_optimizer_theorem_claimed"
    ]
    assert result["new_DFT_or_scheduler_writes"] == 0
    assert not result["physical_reserved_geometry_generated_or_read"]
    assert not result["original_production_source_domain_target_or_inputs_modified"]
    for row in result["rows"]:
        assert (
            row["all_evaluated_radius_valid"]
            and row["unique_analytic_evaluations"] <= 12
        )
        if row["status"] == "gradient_candidate":
            assert row["candidate_released_gradient_norm"] <= 0.03


def test_interior_and_nonlinear_controls_are_kept_when_cost_does_not_improve():
    rows = b.benchmark()["rows"]
    for name, count in [("interior_quadratic", 4), ("nonlinear_valley", 9)]:
        arms = [r for r in rows if r["case"] == name]
        assert len(arms) == 2 and all(
            r["unique_analytic_evaluations"] == count for r in arms
        )


def test_exterior_minimum_refuses_a_boundary_candidate_in_both_arms():
    rows = [
        r for r in b.benchmark()["rows"] if r["case"] == "exterior_quadratic_minimum"
    ]
    assert len(rows) == 2
    assert all(
        r["status"] == "no_candidate" and r["candidate_coordinates"] is None
        for r in rows
    )


def test_near_wall_observation_is_a_small_analytic_count_difference_not_ASE_or_DFT_speedup():
    old, new = [
        r
        for r in b.benchmark()["rows"]
        if r["case"] == "near_wall_misaligned_positive_seed"
    ]
    assert (
        old["unique_analytic_evaluations"] == 7
        and new["unique_analytic_evaluations"] == 6
    )
    assert old["guard_counts"]["rejected"] == 7 and new["guard_counts"]["rejected"] == 0
