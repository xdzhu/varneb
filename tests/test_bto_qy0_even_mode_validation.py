"""Reproduce the frozen BTO even-mode fit without DFT or network access."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from scripts.plan_bto_qy0_even_mode_validation import even_basis, plan


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "benchmarks" / "numerical_integrity"


def test_frozen_plan_recomputes_from_nine_grid_nodes_only() -> None:
    measured = json.loads((DATA / "bto_qy0_ten_measured_nodes_20260928.json").read_text())
    linear = json.loads((DATA / "bto_qy0_three_holdout_gate_20260928.json").read_text())
    frozen = json.loads((DATA / "bto_qy0_even_mode_two_holdout_plan_20260928.json").read_text())
    recomputed = plan(measured, linear)
    for key in recomputed:
        assert recomputed[key] == frozen[key]
    assert recomputed["n_training_grid_nodes"] == 9
    assert len(recomputed["retroactively_evaluated_not_independent_for_model_selection"]) == 4
    assert all(abs(point["observed_minus_predicted_meV_per_BTO"]) < 2
               for point in recomputed["retroactively_evaluated_not_independent_for_model_selection"])
    assert recomputed["new_holdouts"][0]["model_prediction_energy_minus_C_eV_per_BTO"] == pytest.approx(-0.0362966618383853)


def test_even_basis_respects_axis_sign_and_exchange() -> None:
    vectors = even_basis(np.array([[0.4, 0.2], [-0.4, 0.2], [0.2, 0.4]]))
    assert np.allclose(vectors[0], vectors[1])
    assert np.allclose(vectors[0], vectors[2])
