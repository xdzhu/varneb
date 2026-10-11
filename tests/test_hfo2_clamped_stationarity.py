"""Real native-gradient replay plus mask/subspace diagnostic semantics."""
from pathlib import Path

import numpy as np
import pytest

from scripts.analyze_hfo2_clamped_stationarity import analyze, FIELDS
from scripts.analyze_hfo2_clamped_reference import terminal_images
from tests.test_hfo2_clamped_residual import synthetic_chain

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/"benchmarks/hfo2_channels/20261008/clamped_M_terminal_E068_20261011"
AUDIT_SHA = "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c"


def test_actual_M_peak_reaction_is_not_allowed_gradient():
    result = analyze(CASE, AUDIT_SHA)
    peak = result["sampled_peak"]
    assert result["native_portable_frames_reparsed"] == 9 and result["new_DFT_calls"] == 0
    assert peak["image_index"] == 3
    assert peak["true_generalized_force_max_vector_eV_per_A"] == pytest.approx(2.065093534498197)
    assert peak["constraint_reaction_force_max_vector_eV_per_A"] == pytest.approx(2.063130393399425)
    assert peak["allowed_true_generalized_force_max_vector_eV_per_A"] == pytest.approx(.09002380828307668)
    assert peak["allowed_true_tangential_force_eV_per_A"] == pytest.approx(-.12426388097783847)
    assert result["ordinary_fmax_eV_A"] == pytest.approx(.09921408814764762)
    assert result["diagnostics_leave_ordinary_forces_byte_identical"]
    assert not result["G3_selection_or_holdout_access"]
    assert peak["true_perpendicular_force_eV_per_A"] > 2.6  # legacy active-raw field preserved
    assert peak["allowed_true_perpendicular_force_eV_per_A"] < .3


def test_full_projector_and_public_diagnostic_agree_without_force_change():
    chain = synthetic_chain()
    before = chain.get_forces().copy()
    saddle, path = chain.saddle_diagnostics(), chain.path_diagnostics()
    i = saddle["image_index"]
    mask = chain._active_x_mask()
    raw = chain._last_true_forces_x[i]*mask
    allowed = chain._project_constraint_force(raw)*mask
    reaction = raw-allowed
    assert abs(np.dot(allowed, reaction)) < 1e-10
    assert saddle["allowed_true_generalized_force_norm_eV_per_A"] == pytest.approx(np.linalg.norm(allowed))
    assert saddle["constraint_reaction_force_norm_eV_per_A"] == pytest.approx(np.linalg.norm(reaction))
    for field in FIELDS:
        assert saddle[field] == pytest.approx(path["images"][i][field])
    assert np.array_equal(before, chain.get_forces())


def test_unconstrained_active_space_has_no_projector_reaction():
    chain = synthetic_chain()
    # Cells are unchanged; removing the restriction only for this synthetic
    # diagnostic test does not alter any production input or material chain.
    chain.mode_basis = None
    saddle = chain.saddle_diagnostics()
    assert saddle["constraint_reaction_force_norm_eV_per_A"] == 0
    assert saddle["constraint_reaction_force_max_vector_eV_per_A"] == 0
    assert saddle["true_generalized_force_max_vector_eV_per_A"] == saddle["allowed_true_generalized_force_max_vector_eV_per_A"]
    assert saddle["true_tangential_force_eV_per_A"] == saddle["allowed_true_tangential_force_eV_per_A"]


def test_pinned_native_audit_rejects_changed_receipt():
    with pytest.raises(ValueError):
        terminal_images(CASE, "0"*64)
