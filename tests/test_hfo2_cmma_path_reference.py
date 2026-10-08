from pathlib import Path
from types import SimpleNamespace
import json

from ase.io import read
import numpy as np
import pytest

from scripts.analyze_hfo2_cmma_path_reference import (
    endpoint_assignment, lowest_complete_doublet_indices, prepare_frame,
    proper_signed_frames, register_frames,
)


ROOT = Path(__file__).resolve().parents[1]


def tetragonal():
    return read(ROOT / "benchmarks/hfo2_channels/20261008/reference_variants/T.vasp")


def test_all_proper_signed_frames_are_bounded_and_unique():
    rotations = list(proper_signed_frames())
    assert len(rotations) == 24
    assert len({tuple(r.ravel()) for r in rotations}) == 24
    for r in rotations:
        assert np.array_equal(r @ r.T, np.eye(3))
        assert round(np.linalg.det(r)) == 1


def test_endpoint_assignment_reports_global_alternative_and_preserves_inputs():
    target = tetragonal()
    source = target.copy()
    source.set_scaled_positions(source.get_scaled_positions(wrap=False) - .25)
    before = source.positions.copy()
    fit = endpoint_assignment(source, target, np.eye(3), [.25]*3, target.cell.array, second_best=True)
    assert fit["source_to_target"] == list(range(12))
    assert fit["maximum_site_offset_A"] < 1e-14
    assert fit["second_best_assignment_cost_gap_A"] > 1
    assert np.array_equal(source.positions, before)


def test_registration_scores_geometry_before_any_modes_and_keeps_ties():
    target = tetragonal()
    author_t = target.copy()
    author_t.set_scaled_positions(author_t.get_scaled_positions(wrap=False) - .25)
    cmma = author_t.copy()
    score, selected, gate = register_frames(author_t, cmma, target, target)
    assert len(score) == 192
    assert gate["T_minimum_cost_sum_A"] < 1e-13
    assert len(selected) > 1
    assert "mode coverage" in gate["not_used_for_selection"]
    assert all(s["registered_PO_assignment"]["second_best_assignment_cost_gap_A"] > 1 for s in selected)


@pytest.mark.parametrize("fault", ["improper", "sheared_metric", "species", "nan", "empty"])
def test_invalid_registration_never_fits_or_repairs(fault):
    target = tetragonal()
    source, r, cell = target.copy(), np.eye(3), target.cell.array.copy()
    if fault == "improper": r[0, 0] = -1
    elif fault == "sheared_metric": cell[0, 1] = .3
    elif fault == "species": source.numbers[0] = 8
    elif fault == "nan": source.positions[0, 0] = np.nan
    elif fault == "empty": source = source[:0]
    with pytest.raises(ValueError):
        endpoint_assignment(source, target, r, [.25]*3, cell)


def test_frame_and_common_chart_transport_do_not_modify_material_references():
    t = tetragonal()
    c = t.copy()
    c.set_cell(np.diag([4., 6., 7.]), scale_atoms=True)
    before = c.positions.copy(), c.cell.array.copy(), t.positions.copy()
    rng = np.random.default_rng(16)
    e = np.linalg.qr(rng.normal(size=(36, 4)))[0].astype(complex)
    modes = SimpleNamespace(mass_weighted_eigenvectors=e)
    p = [1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10]
    frame = {"rotation_and_integer_basis": [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
             "origin_fractional": [.25, .25, .25], "registered_PO_assignment": {"source_to_target": p}}
    author = {"source_to_target_indices": list(range(12)), "rotation_Cartesian_proper": np.eye(3)}
    reference, basis, audit = prepare_frame(modes, c, t, author, frame)
    assert basis.shape == (36, 4)
    assert np.allclose(basis.T @ basis, np.eye(4), atol=1e-14)
    assert np.array_equal(reference.cell.array, t.cell.array)
    assert np.array_equal(c.positions, before[0]) and np.array_equal(c.cell.array, before[1])
    assert np.array_equal(t.positions, before[2])
    assert audit["ASR_or_translation_removal_from_basis"] is False
    assert audit["maximum_relative_complex_projection_error"] < 1e-14


def test_lowest_four_modes_preserve_two_complete_doublets():
    assert np.array_equal(lowest_complete_doublet_indices([0, 1, 1, 2, 2, 3], np.arange(1, 6)), [1, 2, 3, 4])
    for frequencies in ([0, 1, 1, 1, 2, 2], [0, 1, 1, 2, 2, 2], [0, 1, 1, 1, 1, 3]):
        with pytest.raises(ValueError, match="doublets"):
            lowest_complete_doublet_indices(frequencies, np.arange(1, 6))


def test_archived_material_registration_keeps_ambiguity_metric_and_claim_limits():
    report = json.loads((ROOT / "benchmarks/hfo2_channels/20261008/cmma_path_mapping/analysis.json").read_text())
    assert len(report["frame_scores_geometry_only"]) == 192
    assert report["registration"]["T_compatible_frame_count"] == 16
    assert report["registration"]["selected_frame_count"] == 4
    assert report["new_DFT_calls"] == 0 and report["physical_parameters_changed"] is False
    observation = report["observations"][0]
    assert observation["audited_images"] == 10
    assert observation["replayed_fmax_eV_A"] == pytest.approx(.059881615690252875)
    assert observation["common_T_origin_zero_displacement"][0]
    assert observation["T_lowest_two_optical_doublets_indices_zero_based"] == [3, 4, 5, 6]
    for frame in observation["Cmma_reference_frames"]:
        assert frame["basis_audit"]["real_rank"] == 4
        assert frame["basis_audit"]["ASR_or_translation_removal_from_basis"] is False
        assert frame["basis_audit"]["rigid_translation_overlap_subspace_Frobenius_norm"] > .003
        assert .271 < frame["common_T_origin_fraction"][3] < .273
        assert .522 < frame["common_T_origin_fraction"][-1] < .523
        assert frame["Cmma_affine_origin_residual_sqrt_amu_A"][0] > 6
