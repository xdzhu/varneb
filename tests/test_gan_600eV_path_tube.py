import pytest

from scripts.audit_gan_600eV_path_tube import leave_one_anchor_out_errors


def test_linear_offpath_excess_has_zero_anchor_holdout_error():
    values = {
        0: {"minus": {"delta_enthalpy_meV_per_GaN": 0.0},
            "plus": {"delta_enthalpy_meV_per_GaN": 1.0}},
        1: {"minus": {"delta_enthalpy_meV_per_GaN": 1.0},
            "plus": {"delta_enthalpy_meV_per_GaN": 2.0}},
        2: {"minus": {"delta_enthalpy_meV_per_GaN": 2.0},
            "plus": {"delta_enthalpy_meV_per_GaN": 3.0}},
    }
    result = leave_one_anchor_out_errors([0, 1, 2], [0.0, 0.5, 1.0], values)
    assert len(result) == 2
    assert all(item["absolute_error_meV_per_GaN"] == pytest.approx(0.0)
               for item in result)


def test_curved_excess_is_detected_without_interpolating_barrier():
    values = {
        i: {side: {"delta_enthalpy_meV_per_GaN": (10.0 if i == 1 else 0.0)}
            for side in ("minus", "plus")}
        for i in (0, 1, 2)
    }
    result = leave_one_anchor_out_errors([0, 1, 2], [0.0, 0.5, 1.0], values)
    assert [item["absolute_error_meV_per_GaN"] for item in result] == [10.0, 10.0]
