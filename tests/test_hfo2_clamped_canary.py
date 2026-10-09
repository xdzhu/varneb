"""Synthetic protocol fixtures, not material convergence or DFT evidence."""
import copy
from pathlib import Path

import pytest

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_clamped_canary import validate_summary


def recipe():
    manifest = {"phase_label": "PO_plus", "strain": 0., "pressure_gpa": 0.,
                "allow_tilt": True, "physical_contract_sha256": CONTRACT}
    summary = {"steps_requested": 4, "optimizer_steps": 4, "fmax_target_eV_per_A": .03,
               "open_stress_target_kbar": 2., "maxstep": .02, "external_pressure_gpa": 0.,
               "strain": 0., "allow_tilt": True, "physical_contract_sha256": CONTRACT,
               "converged": False, "status": "step_limit", "max_atomic_force_eV_per_A": .11,
               "open_traction_norm_kbar": 5.}
    return summary, manifest


def test_healthy_step_limit_is_not_an_electronic_failure_or_endpoint_pass():
    s, m = recipe()
    validate_summary(s, m, 5)
    s.update(converged=True, status="completed", max_atomic_force_eV_per_A=.01, open_traction_norm_kbar=1.)
    validate_summary(s, m, 5)


@pytest.mark.parametrize("key,value", [("steps_requested", 5), ("maxstep", .03),
    ("open_stress_target_kbar", .1), ("strain", .005), ("converged", True),
    ("optimizer_steps", False), ("max_atomic_force_eV_per_A", float("nan")),
    ("open_traction_norm_kbar", float("inf"))])
def test_canary_recipe_or_false_physical_claims_refused(key, value):
    s, m = recipe()
    s[key] = value
    with pytest.raises(ValueError):
        validate_summary(s, m, 5)


@pytest.mark.parametrize("count", [0, 4, 6])
def test_missing_repeated_or_over_budget_calls_refused(count):
    s, m = recipe()
    with pytest.raises(ValueError, match="max5"):
        validate_summary(s, m, count)


def test_exact_generic_attachment_directory_and_no_alternative_namespace():
    root = Path(__file__).resolve().parents[1]
    body = (root/"scripts/audit_hfo2_clamped_canary.py").read_text()
    assert 'calculator/image_0000' in body
    assert 'calculator/00' not in body
    script = (root/"cluster/hf_hfo2_clamped_canary_E045_20261009.slurm").read_text()
    assert '--steps 4' in script and '#SBATCH --time=00:30:00' in script
    assert '-m scripts.audit_hfo2_clamped_canary' in script
    assert 'srun' not in script and 'sbatch' not in script
