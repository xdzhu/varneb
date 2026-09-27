"""The off-axis canary must not require a preferred branch energy sign."""

import pytest

from scripts.audit_bto_transverse_soft_canary import (
    audit,
    branch_energy_diagnostics,
)


def test_valid_stiff_off_axis_canary_is_not_rejected_for_missing_lowering():
    rows = [
        {"label": "frozen", "relative_energy_eV_per_BTO": -0.03},
        {"label": "+Q_y", "relative_energy_eV_per_BTO": -0.02},
        {"label": "-Q_y", "relative_energy_eV_per_BTO": -0.01},
    ]
    diagnostic = branch_energy_diagnostics(rows)
    assert diagnostic["plus_branch_lowers_frozen"] is False
    assert diagnostic["frozen_minus_plus_branch_meV_per_BTO"] == pytest.approx(-10)
    assert diagnostic["signed_branch_energy_difference_meV_per_BTO"] == pytest.approx(-10)


def test_preflight_cannot_escape_audited_input_directory(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    with pytest.raises(ValueError, match="direct child"):
        audit(inputs, tmp_path / "work", expected_job_id="1",
              preflight_path=tmp_path / "other.json")


def test_missing_or_nonfinite_three_branch_energy_is_rejected():
    with pytest.raises(ValueError, match="three reviewed branch labels"):
        branch_energy_diagnostics([
            {"label": "frozen", "relative_energy_eV_per_BTO": 0},
        ])
    with pytest.raises(ValueError, match="finite"):
        branch_energy_diagnostics([
            {"label": "frozen", "relative_energy_eV_per_BTO": 0},
            {"label": "+Q_y", "relative_energy_eV_per_BTO": float("nan")},
            {"label": "-Q_y", "relative_energy_eV_per_BTO": 0},
        ])
