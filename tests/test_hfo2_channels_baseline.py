"""Tests of normalization/identity, never a substitute for SCF audits."""

import pytest

from scripts.audit_hfo2_channels_baseline import normalize_summary


def summary():
    return {"path_diagnostics": {"n_atoms": 12, "pressure_eV_per_A3": 0},
            "image_enthalpies_eV": [-100, -99.96, -100.3], "n_images": 3,
            "barrier_enthalpy_eV": 0.04, "reaction_enthalpy_eV": -0.3,
            "final_max_generalized_force_eV_per_A": 0.09,
            "fmax_target_eV_per_A": 0.05, "status": "step_limit_reached",
            "calculator_parameters": {"ecutwfc": 100, "kpts": [2, 2, 2], "scf_thr": 1e-8}}


def test_units_reverse_identity_and_original_threshold_qualifier():
    result = normalize_summary(summary())
    assert result["forward_meV_fu"] == pytest.approx(10)
    assert result["reverse_meV_fu"] == pytest.approx(85)
    assert result["reaction_meV_fu"] == pytest.approx(-75)
    assert result["barrier_identity_error_eV_cell"] == pytest.approx(0)
    assert not result["meets_original_force_target"]
    assert result["meets_current_ordinary_0p10_force_target"]
    assert not result["peak_is_stationary_certified"]
    assert result["summary_status"] == "step_limit_reached"


@pytest.mark.parametrize("error", ["cutoff", "atoms", "barrier", "partial", "nonfinite"])
def test_mixed_or_incomplete_source_is_rejected(error):
    data = summary()
    if error == "cutoff":
        data["calculator_parameters"]["ecutwfc"] = 120
    elif error == "atoms":
        data["path_diagnostics"]["n_atoms"] = 24
    elif error == "barrier":
        data["barrier_enthalpy_eV"] = 0.05
    elif error == "partial":
        data["image_enthalpies_eV"].pop()
    else:
        data["image_enthalpies_eV"][1] = float("nan")
    with pytest.raises(ValueError):
        normalize_summary(data)
